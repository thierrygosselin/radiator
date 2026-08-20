# filter_genotyping_vcf --------------------------------------------------------
#' @title Filter SNPs in a VCF based on genotyping / missing rate (bcftools)
#'
#' @description
#' Wrapper around \code{bcftools +fill-tags} and \code{bcftools view/query} to:
#' \itemize{
#'   \item keep only markers with missing proportion below or equal to a
#'   user-defined threshold, and
#'   \item generate a blacklist of markers above that threshold.
#' }
#'
#' The missing proportion per site is taken from the \code{F_MISSING} INFO tag,
#' computed by \code{bcftools +fill-tags}. A variant is kept if:
#' \preformatted{
#' F_MISSING <= genotyping.threshold
#' }
#'
#' This makes the function well suited as a lightweight, VCF-level pruning step
#' (e.g. before merging scaffold-level VCFs), without opening a GDS.
#'
#' @param vcf (character)
#' Path to the input VCF (\code{.vcf} or \code{.vcf.gz}).
#'
#' @param genotyping.threshold (double)
#' Maximum allowed missing proportion per site (i.e. \code{F_MISSING}).
#' For example, \code{genotyping.threshold = 0.5} keeps markers with at most
#' 50\% missing genotypes. Must be between 0 and 1.
#'
#' @param out.vcf (optional, character)
#' Output VCF filename for filtered markers.
#' If \code{NULL}, the function generates a filename from the input VCF root by
#' appending \code{".genotyping.vcf"} and \code{".gz"} when
#' \code{compress = TRUE}.
#' Default: \code{out.vcf = NULL}.
#'
#' @param blacklist.file (optional, character)
#' Path to a TSV file storing blacklisted markers (\code{CHROM}, \code{POS},
#' \code{ID}).
#' If \code{NULL}, the function generates a filename from the VCF root by
#' appending \code{".genotyping_blacklist.tsv"}.
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
#' @param keep.filled.vcf (logical)
#' Keep the intermediate VCF with \code{F_MISSING} filled-in
#' (the output of \code{bcftools +fill-tags})?
#' If \code{FALSE}, the function removes this file and its index at the end.
#' Default: \code{keep.filled.vcf = FALSE}.
#'
#' @param verbose (logical)
#' Show progress messages and bcftools command lines.
#' Default: \code{verbose = TRUE}.
#'
#' @details
#' The function operates entirely at the VCF level:
#' \enumerate{
#'   \item \code{bcftools +fill-tags} is used to compute \code{F_MISSING} for
#'   each site (if not already present);
#'   \item \code{bcftools view -i 'F_MISSING<=X'} keeps well-genotyped markers;
#'   \item \code{bcftools query -i 'F_MISSING>X'} writes a blacklist of
#'   removed markers.
#' }
#'
#' This is conceptually similar to \code{\link{filter_genotyping}}, but:
#' \itemize{
#'   \item \code{\link{filter_genotyping_vcf}} works on the raw VCF using
#'   \code{F_MISSING} (per-site missing proportion),
#'   \item \code{\link{filter_genotyping}} works on a GDS and uses the full
#'   radiator statistics and plotting machinery.
#' }
#'
#' @return
#' Invisibly returns a list with:
#' \itemize{
#'   \item \code{genotyping.vcf} – path to the filtered VCF;
#'   \item \code{blacklist} – path to the blacklist of high-missing markers;
#'   \item \code{log} – path to the bcftools log file;
#'   \item \code{filled.vcf} – path to the intermediate VCF with \code{F_MISSING}
#'   (may have been deleted if \code{keep.filled.vcf = FALSE}).
#' }
#'
#' @examples
#' \dontrun{
#' # Keep markers with at most 50% missing genotypes
#' filter_genotyping_vcf(
#'   vcf               = "chr1_freebayes.vcf.gz",
#'   genotyping.threshold = 0.5,
#'   out.vcf           = "chr1_freebayes_geno.vcf.gz",
#'   blacklist.file    = "chr1_freebayes_geno_blacklist.tsv",
#'   bcftools.path     = "bcftools"
#' )
#' }
#'
#' @name filter_genotyping_vcf
#' @rdname filter_genotyping_vcf
#' @export
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
filter_genotyping_vcf <- function(
    vcf,
    genotyping.threshold,
    out.vcf         = NULL,
    blacklist.file  = NULL,
    bcftools.path   = "bcftools",
    compress        = TRUE,
    index           = TRUE,
    keep.filled.vcf = FALSE,
    verbose         = TRUE
) {

  # Basic checks ---------------------------------------------------------------
  if (!file.exists(vcf)) {
    rlang::abort(stringi::stri_join("VCF not found: ", vcf))
  }

  if (missing(genotyping.threshold) || is.null(genotyping.threshold)) {
    rlang::abort("`genotyping.threshold` is missing. Provide a value between 0 and 1.")
  }

  if (!is.numeric(genotyping.threshold) ||
      length(genotyping.threshold) != 1L ||
      is.na(genotyping.threshold) ||
      genotyping.threshold < 0 || genotyping.threshold > 1) {
    rlang::abort("`genotyping.threshold` must be a single numeric between 0 and 1.")
  }

  # Ensure bcftools is available -----------------------------------------------
  bcftools_require(bcftools.path = bcftools.path)

  thr <- genotyping.threshold
  thr.label <- sprintf("%.3f", thr)

  # Paths and names (everything lives next to the VCF) -------------------------
  vcf.base <- basename(vcf)
  vcf.dir  <- dirname(vcf)
  vcf.root <- sub("\\.(vcf|vcf\\.gz)$", "", vcf.base, ignore.case = TRUE)

  # Intermediate VCF with F_MISSING tag
  filled.vcf <- file.path(
    vcf.dir,
    paste0(vcf.root, ".filltags.F_MISSING.vcf.gz")
  )

  # Filtered VCF path
  if (is.null(out.vcf)) {
    out.vcf <- paste0(
      vcf.root,
      ".genotyping.vcf",
      if (compress) ".gz" else ""
    )
  } else {
    has_gz <- stringi::stri_endswith_fixed(out.vcf, ".gz")
    if (compress && !has_gz) out.vcf <- paste0(out.vcf, ".gz")
  }
  out.vcf <- file.path(vcf.dir, out.vcf)

  # Blacklist path
  if (is.null(blacklist.file)) {
    blacklist.file <- paste0(vcf.root, ".genotyping_blacklist.tsv")
  }
  blacklist.file <- file.path(vcf.dir, blacklist.file)

  # Safety: ensure blacklist.file is a single valid path
  if (!is.character(blacklist.file) || length(blacklist.file) != 1L || is.na(blacklist.file) || !nzchar(blacklist.file)) {
    rlang::abort("blacklist.file must be a length-1, non-NA character path.")
  }

  # Log file (next to the VCF)
  log.file <- file.path(
    vcf.dir,
    paste0(
      vcf.root,
      ".bcftools_genotyping_",
      format(Sys.time(), "%Y%m%d@%H%M%S"),
      ".log"
    )
  )

  # 1) Fill F_MISSING with +fill-tags ------------------------------------------
  args.fill <- c(
    "+fill-tags",
    vcf,
    "-Oz",
    "-o", filled.vcf,
    "--",
    "-t", "F_MISSING"
  )

  if (verbose) {
    message(
      "\nComputing F_MISSING with bcftools +fill-tags ",
      "(genotyping.threshold = ", thr.label, ")"
    )
  }

  bcftools_exec(
    bcftools.path  = bcftools.path,
    args           = args.fill,
    log.file       = log.file,
    label          = "bcftools +fill-tags (F_MISSING)",
    verbose        = verbose,
    fail_on_status = TRUE
  )

  # 2) Generate filtered VCF (bcftools view -i 'F_MISSING<=thr') ---------------
  expr_keep <- sprintf("F_MISSING<=%.6f", thr)

  args.view <- c(
    "view",
    "-i", expr_keep,
    if (compress) "-Oz" else "-Ov",
    "-o", out.vcf,
    filled.vcf
  )

  if (verbose) message("Generating genotyping-filtered VCF → ", out.vcf)

  bcftools_exec(
    bcftools.path  = bcftools.path,
    args           = args.view,
    log.file       = log.file,
    label          = "bcftools view (genotyping filter)",
    verbose        = verbose,
    fail_on_status = TRUE
  )

  # 3) Optional index (bcftools index -t out.vcf) ------------------------------
  if (compress && index) {
    if (verbose) message("Indexing genotyping-filtered VCF → ", out.vcf, ".tbi")

    res.idx <- bcftools_exec(
      bcftools.path  = bcftools.path,
      args           = c("index", "-t", out.vcf),
      log.file       = log.file,
      label          = "bcftools index (genotyping)",
      verbose        = verbose,
      fail_on_status = FALSE
    )

    if (!identical(res.idx$status, 0L)) {
      warning("bcftools index failed. See log: ", log.file, call. = FALSE)
    }
  }

  # 4) Query blacklist (markers with F_MISSING > thr) --------------------------
  expr_black <- sprintf("F_MISSING>%.6f", thr)

  # KEY FIX: let bcftools write the file directly (-o), no stdout parsing
  args.query <- c(
    "query",
    "-f", "%CHROM\\t%POS\\t%ID\\n",
    "-i", expr_black,
    "-o", blacklist.file,
    filled.vcf
  )

  if (verbose) message("Writing genotyping blacklist → ", blacklist.file)

  bcftools_exec(
    bcftools.path  = bcftools.path,
    args           = args.query,
    log.file       = log.file,
    label          = "bcftools query (genotyping blacklist)",
    verbose        = verbose,
    fail_on_status = TRUE
  )

  # 5) Optional cleanup of filled.vcf ------------------------------------------
  if (!keep.filled.vcf) {
    if (file.exists(filled.vcf)) unlink(filled.vcf)
    filled.tbi <- paste0(filled.vcf, ".tbi")
    if (file.exists(filled.tbi)) unlink(filled.tbi)
    filled.vcf <- NULL
  }


  invisible(list(
    genotyping.vcf = out.vcf,
    blacklist    = blacklist.file,
    log          = log.file,
    filled.vcf   = filled.vcf
  ))
} # END filter_genotyping_vcf
