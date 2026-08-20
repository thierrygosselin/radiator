# filter_monomorphic_vcf -------------------------------------------------------
#' @title Filter monomorphic SNPs in a VCF using bcftools (AC/AN-based)
#'
#' @description
#' Wrapper around \code{bcftools view} and \code{bcftools query} to:
#' \itemize{
#'   \item keep only polymorphic markers in a VCF, and
#'   \item generate a blacklist of monomorphic sites.
#' }
#'
#' A site is considered polymorphic in the sample if at least one ALT allele has
#' an allele count strictly between 0 and the total allele count
#' (\code{INFO/AN}). In bcftools expression syntax:
#' \preformatted{
#' INFO/AC>0 && INFO/AC<INFO/AN
#' }
#'
#' This works for multi-ALT sites: as soon as one ALT allele is segregating
#' (some REF and some ALT\eqn{_i}), the variant is kept. Sites fixed for REF or
#' fixed for a single ALT allele are dropped as monomorphic.
#'
#' The function is intended for use on freshly created VCFs where \code{INFO/AC}
#' and \code{INFO/AN} are trustworthy (e.g. direct output from FreeBayes or
#' bcftools). It does not inspect genotypes or GDS, only the VCF INFO fields.
#'
#' @param vcf (character)
#' Path to the input VCF (\code{.vcf} or \code{.vcf.gz}).
#'
#' @param out.vcf (optional, character)
#' Output VCF filename for polymorphic markers.
#' If \code{NULL}, the function generates a filename from the input VCF root by
#' appending \code{".polymorphic.vcf"} and \code{".gz"} when
#' \code{compress = TRUE}.
#' Default: \code{out.vcf = NULL}.
#'
#' @param blacklist.file (optional, character)
#' Path to a TSV file storing monomorphic markers (\code{CHROM}, \code{POS},
#' \code{ID}).
#' If \code{NULL}, the function generates a filename from the VCF root by
#' appending \code{".monomorphic_blacklist.tsv"}.
#' Default: \code{blacklist.file = NULL}.
#'
#' @param bcftools.path (character)
#' Path or name of the \code{bcftools} executable (e.g. \code{"bcftools"} or a
#' full path inside a conda environment).
#' Default: \code{bcftools.path = "bcftools"}.
#'
#' @param compress (logical)
#' Compress the output VCF using bgzip (\code{-Oz}).
#' Default: \code{compress = TRUE}.
#'
#' @param index (logical)
#' Index the output VCF with tabix when \code{compress = TRUE}.
#' Default: \code{index = TRUE}.
#'
#' @param verbose (logical)
#' Show progress messages and bcftools command lines.
#' Default: \code{verbose = TRUE}.
#'
#' @details
#' \strong{Important distinction — allele-level vs genotype-level polymorphism}
#'
#' \code{\link{filter_monomorphic_vcf}} identifies polymorphic sites solely from VCF
#' INFO fields (\code{AC} and \code{AN}), using the bcftools expression:
#'
#' \preformatted{
#' INFO/AC>0 && INFO/AC<INFO/AN
#' }
#'
#' A variant is kept if \strong{at least one ALT allele} has a non-zero allele count
#' and is not fixed in the sample. This is an \strong{allele-level definition of
#' polymorphism}:
#'
#' If some chromosomes carry REF and some carry ALT\eqn{_i}, the variant is
#' considered polymorphic—even if all individuals share the same genotype
#' (e.g., all REF/ALT heterozygotes).
#'
#' This pre-filtering step is intentionally lightweight and designed for
#' \strong{VCF-level slimming before merging for example scaffolds}.
#' It does not inspect genotypes and does not attempt to detect genotype-level
#' monomorphism.
#'
#' As a result, the set of monomorphic markers detected here may differ from
#' those detected later inside a GDS by \code{\link{filter_monomorphic}}.
#'
#' @note
#' \strong{Why results differ between \code{\link{filter_monomorphic_vcf}} and
#' \code{\link{filter_monomorphic}}}
#'
#' These functions intentionally implement two different biological definitions:
#' \itemize{
#'   \item \code{\link{filter_monomorphic_vcf}} removes sites that are
#'   \strong{allele-monomorphic}, based on INFO/AC and INFO/AN.
#'   \item \code{\link{filter_monomorphic}} removes sites that are
#'   \strong{genotype-monomorphic}, based on the distribution of genotype
#'   phenotypes in the GDS.
#' }
#'
#' A variant can contain both REF and ALT alleles (allele-level polymorphism)
#' but still have only one genotype phenotype across all individuals.
#' Therefore, it is expected and correct that
#' \code{\link{filter_monomorphic}} may remove additional markers.
#'
#' @return
#' Invisibly returns a list with:
#' \itemize{
#'   \item \code{polymorphic.vcf} – path to the filtered VCF;
#'   \item \code{blacklist} – path to the monomorphic marker list;
#'   \item \code{log} – path to the bcftools log file.
#' }
#'
#' @examples
#' \dontrun{
#' filter_monomorphic_vcf(
#'   vcf            = "chr1_freebayes.vcf.gz",
#'   out.vcf        = "chr1_freebayes_poly.vcf.gz",
#'   blacklist.file = "chr1_freebayes_mono.tsv",
#'   bcftools.path  = "bcftools"
#' )
#' }
#'
#' @name filter_monomorphic_vcf
#' @rdname filter_monomorphic_vcf
#' @export
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
filter_monomorphic_vcf <- function(
    vcf,
    out.vcf        = NULL,
    blacklist.file = NULL,
    bcftools.path  = "bcftools",
    compress       = TRUE,
    index          = TRUE,
    verbose        = TRUE
) {

  # Basic checks ---------------------------------------------------------------
  if (!file.exists(vcf)) {
    rlang::abort(stringi::stri_join("VCF not found: ", vcf))
  }

  # Ensure bcftools is available -----------------------------------------------
  bcftools_require(bcftools.path = bcftools.path)

  # Definition of "polymorphic" for future-me:
  # at least one ALT with 0 < AC < AN.
  polymorphic_expr_raw <- "INFO/AC>0 && INFO/AC<INFO/AN"

  # Paths and names (everything lives next to the VCF) -------------------------
  vcf.base <- basename(vcf)
  vcf.dir  <- dirname(vcf)
  vcf.root <- sub("\\.(vcf|vcf\\.gz)$", "", vcf.base, ignore.case = TRUE)

  # Polymorphic VCF path
  if (is.null(out.vcf)) {
    out.vcf <- paste0(
      vcf.root,
      ".polymorphic.vcf",
      if (compress) ".gz" else ""
    )
  } else {
    has_gz <- stringi::stri_endswith_fixed(out.vcf, ".gz")
    if (compress && !has_gz) {
      out.vcf <- paste0(out.vcf, ".gz")
    }
  }
  out.vcf <- file.path(vcf.dir, out.vcf)

  # Blacklist path
  if (is.null(blacklist.file)) {
    blacklist.file <- paste0(vcf.root, ".monomorphic_blacklist.tsv")
  }
  blacklist.file <- file.path(vcf.dir, blacklist.file)

  # Log file (next to the VCF)
  log.file <- file.path(
    vcf.dir,
    paste0(
      vcf.root,
      ".bcftools_filter_",
      format(Sys.time(), "%Y%m%d@%H%M%S"),
      ".log"
    )
  )

  # 1) Generate polymorphic VCF (bcftools view -i 'expr') ----------------------
  args.poly <- c(
    "view",
    "-i", polymorphic_expr_raw,
    if (compress) "-Oz" else "-Ov",
    "-o", out.vcf,
    vcf
  )

  if (verbose) {
    message("\nGenerating polymorphic VCF \u2192 ", out.vcf)
  }

  bcftools_exec(
    bcftools.path  = bcftools.path,
    args           = args.poly,
    log.file       = log.file,
    label          = "bcftools view (polymorphic)",
    verbose        = verbose,
    fail_on_status = TRUE
  )

  # 2) Optional index (bcftools index -t out.vcf) ------------------------------
  if (compress && index) {
    args.idx <- c(
      "index",
      "-t",
      out.vcf
    )

    if (verbose) {
      message("Indexing polymorphic VCF \u2192 ", out.vcf, ".tbi")
    }

    # For indexing, we usually don't want to abort if it fails completely,
    # but still want the error logged. Let the warning be handled outside if needed.
    res.idx <- bcftools_exec(
      bcftools.path  = bcftools.path,
      args           = args.idx,
      log.file       = log.file,
      label          = "bcftools index",
      verbose        = verbose,
      fail_on_status = FALSE
    )

    if (!identical(res.idx$status, 0L)) {
      warning("bcftools index failed. See log: ", log.file, call. = FALSE)
    }
  }

  # 3) Query monomorphic markers from the original VCF -------------------------
  #    bcftools query -f '%CHROM\t%POS\t%ID\n' -e 'expr' in.vcf
  args.query <- c(
    "query",
    "-f", "%CHROM\\t%POS\\t%ID\\n",
    "-e", polymorphic_expr_raw,
    vcf
  )

  if (verbose) {
    message("Writing monomorphic blacklist \u2192 ", blacklist.file)
  }

  res.query <- bcftools_exec(
    bcftools.path  = bcftools.path,
    args           = args.query,
    log.file       = log.file,
    label          = "bcftools query (monomorphic blacklist)",
    verbose        = verbose,
    fail_on_status = TRUE
  )

  # stdout is a single blob of text with all lines
  out_txt <- res.query$stdout

  if (!nzchar(out_txt)) {
    # No monomorphic sites found
    writeLines(character(), blacklist.file)
  } else {
    bl_lines <- strsplit(out_txt, "\n", fixed = TRUE)[[1]]
    # Drop trailing empty line if present
    if (length(bl_lines) > 0L && identical(bl_lines[[length(bl_lines)]], "")) {
      bl_lines <- bl_lines[-length(bl_lines)]
    }
    writeLines(bl_lines, blacklist.file)
  }

  invisible(list(
    polymorphic.vcf = out.vcf,
    blacklist       = blacklist.file,
    log             = log.file
  ))
} # END filter_monomorphic_vcf
