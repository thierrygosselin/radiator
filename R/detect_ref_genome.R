#' @title Detect reference genome
#' @description Detect if the dataset was de novo or reference genome assembled

#' @param chromosome (string) String with chromosome ID, unique or not.
#' Default: \code{chromosome = NULL}.
#' @param data (optional, path or object) radiator GDS object or file. When this
#' argument is used, the function look inside the GDS radiator node for:
#' \itemize{
#' \item \code{reference.genome} node and \code{TRUE/FALSE} value
#' \item \code{markers.meta} node and detect the information using the \code{CHROM}
#' column
#' \item \code{chromosome} node directly.
#' }
#' Default: \code{data = NULL}.
#'
#'
#' @param verbose (optional, logical) When \code{verbose = TRUE}
#' the function is a little more chatty during execution.
#' Default: \code{verbose = TRUE}.

#' @return TRUE/FALSE if dataset is assembled with a reference genome or not.

#' @examples
#' \dontrun{
#' ref.genome <- radiator::detect_ref_genome(chromosome = tidy.data$CHROM)
#' }


#' @name detect_ref_genome
#' @rdname detect_ref_genome
# @keywords internal
#' @export
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

detect_ref_genome <- function(chromosome = NULL, data = NULL, verbose = TRUE) {
  if (is.null(chromosome) && is.null(data)) {
    ref.genome <- NULL
  } else {
    ref.genome <- NULL

    if (!is.null(data)) {
      # Detect format and import if necessary
      data.type <- radiator::detect_genomic_format(data)
      if (!data.type %in% c("SeqVarGDSClass", "gds.file")) {
        rlang::abort("Input not supported for this function: read function documentation")
      }
      if (!"SeqVarTools" %in% utils::installed.packages()[,"Package"]) {
        rlang::abort('Please install SeqVarTools for this option:\n
           install.packages("BiocManager")
           BiocManager::install("SeqVarTools")')
      }
      if (data.type == "gds.file") {
        data <- radiator::read_rad(data, verbose = verbose)
        data.type <- "SeqVarGDSClass"
      }

      # check radiator node
      radiator.gds <- gdsfmt::index.gdsn(node = data, path = "radiator", silent = TRUE)

      if (is.null(radiator.gds)) {
        chromosome <- SeqArray::seqGetData(data, "chromosome")
        sync.gds <- FALSE
      } else {
        ref.genome <- gdsfmt::get.attr.gdsn(gdsfmt::index.gdsn(
          node = radiator.gds, path = "reference.genome", silent = TRUE))$R.class[1]
        sync.gds <- FALSE
        if (!is.null(ref.genome) && !is.logical(ref.genome)) {
          ref.genome <- NULL
          sync.gds <- TRUE
        }
        # Check the markers.meta field
        if (is.null(ref.genome)) {
          sync.gds <- TRUE
          chromosome <- gdsfmt::index.gdsn(
            node = radiator.gds, path = "markers.meta/CHROM", silent = TRUE)
          if (!is.null(chromosome)) {
            chromosome <- gdsfmt::read.gdsn(chromosome)
          } else {
            chromosome <- SeqArray::seqGetData(data, "chromosome")
          }
        }
      }
    }#End of data

    if (is.null(ref.genome)) {
      if (is.null(chromosome)) {
        if (verbose) message("\nUndetermined chromosome field, will be tagged: de novo\n")
        ref.genome <- FALSE
      } else {
        ref.genome <- sample(x = unique(chromosome),
                             size = min(length(unique(chromosome)), 100),
                             replace = FALSE)
        # if the chrom.unique > 1 more likely not to be de novo assembly (e.g. with old stacks version)
        chrom.unique <- length(unique(ref.genome)) == 1
        chrom.unique.radiator <- any(unique(ref.genome) == "CHROM_1")
        chrom.unique.stacks <- any(unique(ref.genome) == "un")

        # presence of underscore or other separator: more likely ref genome
        chrom.sep <- TRUE %in%
          stringi::stri_detect_regex(str = ref.genome, pattern = "[^[:alnum:]]+") %>%
          unique

        # presence of letters = more likely ref genome
        chrom.alpha <- TRUE %in%
          stringi::stri_detect_regex(str = ref.genome, pattern = "[[:alpha:]]+") %>%
          unique

        if (chrom.unique && chrom.unique.radiator) ref.genome <- FALSE
        if (chrom.unique.radiator) ref.genome <- FALSE
        # if (!chrom.unique && chrom.alpha || chrom.sep) {
        #   ref.genome <- TRUE
        # } else {
        #   ref.genome <- FALSE
        # }
        if (!chrom.unique && chrom.alpha && chrom.sep) {
          ref.genome <- TRUE
        } else {
          ref.genome <- FALSE
        }
        if (chrom.unique.radiator) ref.genome <- FALSE
        if (chrom.unique.stacks) ref.genome <- FALSE


        # stacks related
        if (!is.null(data)) {
          data.source <- radiator::extract_data_source(gds = data)
          if (TRUE %in% stringi::stri_detect_fixed(str = data.source, pattern = "Stacks")) {
            locus.type <- SeqArray::seqGetData(data, "annotation/id")
            locus.missing <- unique(stringi::stri_detect_fixed(
              str = locus.type,
              pattern = "."))
            locus.strands <- TRUE %in% (unique(stringi::stri_detect_fixed(
              str = locus.type,
              pattern = "+")))

            if (locus.missing) ref.genome <- FALSE
            if (locus.strands) ref.genome <- TRUE
            if (!locus.missing && locus.strands) ref.genome <- TRUE
          }
        }
      }
    }
    chrom.unique <- chrom.alpha <- chrom.sep <- chrom.unique.radiator <- NULL

    # sync gds -------------------------------------------------------------------
    if (!is.null(data) && !is.null(radiator.gds) && sync.gds) {
      update_radiator_gds(gds = data, node.name = "reference.genome", value = ref.genome)
    }

    # Result ---------------------------------------------------------------------
    if (ref.genome) {
      if (verbose) message("Reads assembly: reference-assisted")
    } else {
      if (verbose) message("Reads assembly: de novo")
    }
  }
  return(ref.genome)
}#End detect_ref_genome
