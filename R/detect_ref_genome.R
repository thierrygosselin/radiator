#' @title Detect whether a dataset is reference-guided or de novo assembled
#'
#' @description
#' Determine if genomic data were produced using a **reference genome**
#' (reference-assisted assembly) or using a **de novo assembly**.
#'
#' The function integrates multiple sources of evidence depending on the input
#' provided:
#'
#' \strong{1. Radiator GDS (preferred)}
#' If a radiator-generated GDS is supplied, the function checks, in order:
#' \itemize{
#'   \item the \code{/radiator/reference.genome} node (logical value);
#'   \item the \code{/radiator/markers.meta/CHROM} field;
#'   \item the \code{SeqArray} \code{chromosome} node.
#' }
#'
#' \strong{2. VCF imported with SeqArray or read by radiator}
#' When a VCF has been converted to GDS via \pkg{SeqArray}, or when a VCF file
#' is provided directly, the function may detect reference-guided assembly from:
#' \itemize{
#'   \item the VCF header’s \code{##reference=} field;
#'   \item the presence of \code{##contig=<ID=...>} definitions;
#'   \item the structure of chromosome labels.
#' }
#'
#' \strong{3. Data-source-specific heuristics}
#' Additional rules are applied for specific pipelines:
#' \itemize{
#'   \item \strong{Stacks}: optionally, detection of \code{.} or \code{+} in
#'     \code{annotation/id} (see internal comments in the function body);
#'   \item \strong{GATK}: chromosome names containing \code{"contig"};
#'   \item \strong{FreeBayes}: inspection of reference and contig metadata
#'     in the VCF header.
#' }
#'
#' @param data
#' A radiator GDS object, a \pkg{SeqArray} GDS object, or a path to a file.
#' If provided, the function extracts information from radiator metadata,
#' SeqArray header fields, \code{markers.meta}, and other pipeline-specific
#' indicators.
#' Default: \code{data = NULL}.
#'
#' @param verbose (logical, optional)
#' When \code{TRUE}, the function prints messages describing the decision
#' process and final classification.
#' Default: \code{verbose = TRUE}.
#'
#' @return
#' A single logical value:
#' \itemize{
#'   \item \code{TRUE} – reference-assisted assembly;
#'   \item \code{FALSE} – de novo assembly.
#' }
#'
#' @examples
#' \dontrun{
#' # Using a radiator GDS file
#' gds <- radiator::read_rad("my_data.gds")
#' ref.genome <- radiator::detect_ref_genome(data = gds)
#'
#' # Using a VCF file directly
#' ref.genome <- radiator::detect_ref_genome(data = "variants.vcf.gz")
#' }
#'
#' @name detect_ref_genome
#' @rdname detect_ref_genome
#' @export
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#'
detect_ref_genome <- function(data = NULL, verbose = TRUE) {

  ref.genome   <- NULL
  sync.gds     <- FALSE
  radiator.gds <- NULL

  # No info at all -------------------------------------------------------------
  if (is.null(data)) {
    if (verbose) {
      message("\nUndetermined data will be tagged: de novo\n")
      message("Reads assembly: de novo")
    }
    return(FALSE)
  }

  # If we have data, figure out what it is -------------------------------------
  data.type <- radiator::detect_genomic_format(data)

  # VCF file -------------------------------------------------------------------
  if (identical(data.type, "vcf.file")) {
    ref.genome <- FALSE

    detect.source   <- check_header_source_vcf(vcf = data)
    ref.genome.info <- detect.source$check.header$reference
    contig.info     <- detect.source$check.header$contig
    n.contig.info   <- if (is.null(contig.info)) 0L else nrow(contig.info)

    # VCF has a reference header line -> reference-guided
    if (!is.null(ref.genome.info) && length(ref.genome.info) >= 1L) {
      ref.genome <- TRUE
    }

    # VCF has at least one contig definition -> reference-guided
    if (n.contig.info >= 1L) {
      ref.genome <- TRUE
    }

    if (isTRUE(ref.genome)) {
      if (verbose) message("Reads assembly: reference-assisted")
    } else {
      if (verbose) message("Reads assembly: de novo")
      ref.genome <- FALSE
    }
    return(ref.genome)
  }

  # DArT special case ----------------------------------------------------------
  if (identical(data.type, "dart")) {
    ref.genome <- radiator::detect_dart_format(data = data, verbose = FALSE)$ref.genome

    if (isTRUE(ref.genome)) {
      if (verbose) message("Reads assembly: reference-assisted")
    } else {
      if (verbose) message("Reads assembly: de novo")
      ref.genome <- FALSE
    }
    return(ref.genome)
  }

  # SeqArray / GDS -------------------------------------------------------------
  if (!data.type %in% c("SeqVarGDSClass", "gds.file")) {
    rlang::abort("Input not supported for this function: read function documentation")
  }

  radiator_packages_dep(package = "SeqArray", cran = FALSE, bioc = TRUE)

  # if it's a GDS file, import
  if (identical(data.type, "gds.file")) {
    data      <- radiator::read_rad(data, verbose = FALSE)
    data.type <- "SeqVarGDSClass"
  }

  # Check GDS for radiator node ------------------------------------------------
  try.seq.summary <- FALSE
  try.more        <- FALSE

  radiator.gds <- gdsfmt::index.gdsn(node = data, path = "radiator", silent = TRUE)
  if (is.null(radiator.gds)) {
    # no radiator node at all -> do not create one, just use SeqArray heuristics
    try.seq.summary <- TRUE
  } else {
    # try to read existing reference.genome flag -------------------------------
    ref.node <- gdsfmt::index.gdsn(
      node   = radiator.gds,
      path   = "reference.genome",
      silent = TRUE
    )

    if (!is.null(ref.node)) {
      ref.val <- tryCatch(
        gdsfmt::read.gdsn(ref.node),
        error = function(e) NULL  # empty node or no data -> treat as not set
      )

      if (is.logical(ref.val) && length(ref.val) == 1L) {
        # good, scalar logical already stored by radiator
        ref.genome <- ref.val
      } else {
        # FUTURE ME:
        # reference.genome node exists but is empty or weird -> recompute
        # and overwrite, but only for files that already have radiator structure.
        try.seq.summary <- TRUE
        ref.genome      <- NULL
        sync.gds        <- TRUE
      }
    } else {
      # FUTURE ME:
      # reference.genome node is missing inside /radiator.
      # We *do not* create it for external / non-radiator GDS:
      # we will compute ref.genome for this session only (sync.gds stays FALSE).
      try.seq.summary <- TRUE
      sync.gds        <- FALSE
    }
  }

  # try to get the info from the GDS using SeqArray::seqSummary-----------------
  if (try.seq.summary) {
    # 1. get the ref genome and contig info...
    ref.genome.info <- tryCatch(
      SeqArray::seqSummary(gdsfile = data, varname = "$reference", verbose = FALSE),
      error = function(e) NULL
    )

    n.contig.info <- tryCatch(
      nrow(SeqArray::seqSummary(gdsfile = data, varname = "$contig", verbose = FALSE)),
      error = function(e) NULL
    )

    # If the VCF header has a reference path and at least one contig,
    # it's definitely a reference-guided assembly.
    if (!is.null(ref.genome.info) && length(ref.genome.info) >= 1L) {
      ref.genome <- TRUE
    }

    if (!is.null(n.contig.info) && n.contig.info >= 1L) {
      ref.genome <- TRUE
    }

    # if still undecided, we will fall back on chromosome heuristics
    if (is.null(ref.genome)) try.more <- TRUE
  }

  # try a little more ....------------------------------------------------------
  if (try.more) {
    # FUTURE ME:
    # This block uses chromosome names + data.source to make a best-effort call.
    # It's deliberately heuristic and should not be used as the *only* source
    # of truth when better metadata is available.

    # we try to get the chromosomes ...
    if (is.null(radiator.gds)) {
      chromosome <- SeqArray::seqGetData(data, "chromosome")
    } else {
      # if still unknown, try markers.meta/CHROM
      chromosome <- NULL
      chrom.node <- gdsfmt::index.gdsn(
        node   = radiator.gds,
        path   = "markers.meta/CHROM",
        silent = TRUE
      )
      if (!is.null(chrom.node)) {
        chromosome <- gdsfmt::read.gdsn(chrom.node)
      } else {
        chromosome <- tryCatch(
          SeqArray::seqGetData(data, "chromosome"),
          error = function(e) NULL
        )
      }
    }

    if (is.null(chromosome)) {
      ref.genome <- FALSE
    } else {
      data.source <- radiator::extract_data_source(gds = data)

      # GATK -------------------------------------------------------------------
      if (TRUE %in% stringi::stri_detect_fixed(str = data.source, pattern = "GATK")) {

        ref.genome <- any(stringi::stri_detect_fixed(
          str     = unique(chromosome),
          pattern = "contig"
        ))

      } else if (TRUE %in% stringi::stri_detect_fixed(str = data.source, pattern = "freeBayes")) {
        # FreeBayes ------------------------------------------------------------

        ref.header <- tryCatch(
          SeqArray::seqSummary(gdsfile = data, varname = "$reference", verbose = FALSE),
          error = function(e) NULL
        )

        contig.info <- tryCatch(
          SeqArray::seqSummary(gdsfile = data, varname = "$contig", verbose = FALSE),
          error = function(e) NULL
        )

        # If the VCF header has a reference path and at least one contig,
        # it's definitely a reference-guided assembly.
        if (!is.null(ref.header) && length(ref.header) >= 1L) {
          ref.genome <- TRUE
        }

        if (!is.null(contig.info) && NROW(contig.info) >= 1L) {
          ref.genome <- TRUE
        }

      } else {

        # Generic chromosome heuristics ----------------------------------------
        chrom.sample <- sample(
          x       = unique(chromosome),
          size    = min(length(unique(chromosome)), 100L),
          replace = FALSE
        )

        chrom.unique          <- length(unique(chrom.sample)) == 1L
        chrom.unique.radiator <- any(unique(chrom.sample) == "CHROM_1")
        chrom.unique.stacks   <- any(unique(chrom.sample) == "un")

        chrom.sep <- TRUE %in%
          stringi::stri_detect_regex(str = chrom.sample, pattern = "[^[:alnum:]]+") %>%
          unique()

        chrom.alpha <- TRUE %in%
          stringi::stri_detect_regex(str = chrom.sample, pattern = "[[:alpha:]]+") %>%
          unique()

        # base decision
        if (!chrom.unique && chrom.alpha && chrom.sep) {
          ref.genome <- TRUE
        } else {
          ref.genome <- FALSE
        }

        if (chrom.unique.radiator || chrom.unique.stacks) {
          ref.genome <- FALSE
        }

        # FUTURE ME (Stacks heuristics):
        # If you want to resurrect the Stacks-specific logic, it can go here,
        # after data.source detection, using annotation/id to refine ref.genome.
        #
        # if (TRUE %in% stringi::stri_detect_fixed(str = data.source, pattern = "Stacks")) {
        #   locus.type    <- SeqArray::seqGetData(data, "annotation/id")
        #   locus.missing <- TRUE %in% stringi::stri_detect_fixed(str = locus.type, pattern = ".") %>% unique()
        #   locus.strands <- TRUE %in% stringi::stri_detect_fixed(str = locus.type, pattern = "+") %>% unique()
        #
        #   if (locus.missing) ref.genome <- FALSE
        #   if (locus.strands) ref.genome <- TRUE
        #   if (!locus.missing && locus.strands) ref.genome <- TRUE
        # }
      }
    }
  }

  # If still unknown------------------------------------------------------------
  if (is.null(ref.genome)) {
    if (verbose) message("\nUndetermined chromosome field, will be tagged: de novo\n")
    ref.genome <- FALSE
  }

  # sync gds -------------------------------------------------------------------
  # FUTURE ME:
  # Only sync back to GDS when:
  # - there is a /radiator node, and
  # - we decided earlier that it's safe to overwrite (sync.gds == TRUE).
  if (!is.null(data) && !is.null(radiator.gds) && sync.gds) {
    update_radiator_gds(
      gds       = data,
      node.name = "reference.genome",
      value     = ref.genome
    )
  }

  # Result ---------------------------------------------------------------------
  if (isTRUE(ref.genome)) {
    if (verbose) message("Reads assembly: reference-assisted")
  } else {
    if (verbose) message("Reads assembly: de novo")
    ref.genome <- FALSE
  }

  ref.genome
} # End detect_ref_genome

#' @title Extract the reference genome filename/path used for a dataset
#'
#' @description
#' Retrieve the reference genome filename or path associated with a dataset,
#' when available.
#'
#' This helper is intended for internal or downstream use. It does *not*
#' attempt to determine whether the dataset is reference-guided (see
#' \code{\link{detect_ref_genome}} for that). Instead, it simply attempts to
#' recover the reference information from:
#'
#' \strong{1. VCF files}
#' Parsed via \code{check_header_source_vcf()}, reading the \code{##reference=}
#' header line when present.
#'
#' \strong{2. SeqArray / radiator GDS files}
#' Extracted using \code{SeqArray::seqSummary(gds, "$reference")}. If this
#' header field exists, it is returned directly.
#'
#' If no reference path is detected, the function returns \code{NULL} quietly.
#'
#' @param data
#' A VCF file path, a radiator GDS object, a \pkg{SeqArray} GDS object,
#' or a path to a GDS file.
#' Default: \code{data = NULL}.
#'
#' @param verbose (logical)
#' When \code{TRUE}, the function prints optional diagnostics explaining
#' why no reference path could be found.
#' Default: \code{verbose = FALSE}.
#'
#' @return
#' A character scalar containing the reference genome filename or path,
#' or \code{NULL} if no reference could be detected.
#'
#' @examples
#' \dontrun{
#' # From a VCF:
#' ref.path <- radiator::extract_ref_genome("variants.vcf.gz")
#'
#' # From a radiator GDS:
#' gds <- radiator::read_rad("my_data.gds")
#' ref.path <- radiator::extract_ref_genome(gds)
#' }
#'
#' @name extract_ref_genome
#' @rdname extract_ref_genome
#' @keywords internal
#' @export
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#'
extract_ref_genome <- function(data = NULL, verbose = FALSE) {

  # No input -------------------------------------------------------------------
  if (is.null(data)) {
    if (verbose) message("No data provided: cannot extract reference genome filename.")
    return(NULL)
  }

  # Identify data type ---------------------------------------------------------
  data.type <- radiator::detect_genomic_format(data)

  # VCF ------------------------------------------------------------------------
  if (identical(data.type, "vcf.file")) {

    detect.header <- check_header_source_vcf(vcf = data)
    ref.path <- detect.header$check.header$reference

    if (!is.null(ref.path) && length(ref.path) >= 1L) {
      return(ref.path[1L])
    }

    if (verbose) message("VCF: No ##reference= header field found.")
    return(NULL)
  }

  # DArT (placeholder, may expand later) ---------------------------------------
  if (identical(data.type, "dart")) {

    dart.info <- radiator::detect_dart_format(data = data, verbose = FALSE)

    # FUTURE ME:
    # Add something like dart.info$ref.genome.file if DArT exposes it.
    if ("ref.genome.file" %in% names(dart.info)) {
      ref.path <- dart.info$ref.genome.file
      if (!is.null(ref.path)) return(ref.path)
    }

    if (verbose) message("DArT: no reference filename available.")
    return(NULL)
  }

  # GDS ------------------------------------------------------------------------
  if (!data.type %in% c("SeqVarGDSClass", "gds.file")) {
    rlang::abort("Unsupported input: see read function documentation.")
  }

  radiator_packages_dep(package = "SeqArray", cran = FALSE, bioc = TRUE)

  # If it's a file path → read as radiator/SeqArray
  if (identical(data.type, "gds.file")) {
    data      <- radiator::read_rad(data, verbose = FALSE)
    data.type <- "SeqVarGDSClass"
  }

  # Try SeqArray header $reference --------------------------------------------
  ref.path <- tryCatch(
    SeqArray::seqSummary(gdsfile = data, varname = "$reference", verbose = FALSE),
    error = function(e) NULL
  )

  if (!is.null(ref.path) && length(ref.path) >= 1L) {
    return(ref.path[1L])
  }

  # FUTURE ME:
  # You can store ref file inside /radiator/reference.genome.file for radiators,
  # then read it like this:
  #
  # rad <- gdsfmt::index.gdsn(data, "radiator", silent = TRUE)
  # fn  <- gdsfmt::index.gdsn(rad, "reference.genome.file", silent = TRUE)
  # value <- gdsfmt::read.gdsn(fn)
  #
  # But only after you implement saving it.

  if (verbose) message("No reference genome filename detected in GDS.")
  return(NULL)
}#END extract_ref_genome



