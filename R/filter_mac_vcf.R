# filter_mac_vcf --------------------------------------------------------------
#' @title Filter low-MAC variants in a VCF using bcftools (AC-based)
#'
#' @description
#' Wrapper around \code{bcftools view} and \code{bcftools query} to:
#' \itemize{
#'   \item keep only variants with at least one ALT allele having a
#'         minor-allele count (MAC) above a user-defined threshold, and
#'   \item generate a blacklist of low-MAC variants.
#' }
#'
#' The filter is implemented using the \code{INFO/AC} field written by the
#' variant caller. For multi-ALT sites, the bcftools expression
#' \code{INFO/AC >= mac.threshold} evaluates to \code{TRUE} if any ALT allele
#' has \code{AC >= mac.threshold}. Variants that do not satisfy this condition
#' are considered low-MAC and are blacklisted.
#'
#' The function assumes that \code{INFO/AC} (and \code{INFO/AN} if present)
#' come from a recent calling step and are still valid (i.e. no sample
#' subsetting has occurred since calling). It does not recalculate AC/AN;
#' it only uses the INFO fields already in the VCF.
#'
#' This is a lightweight, VCF-level pre-filter that is particularly useful
#' to slim down large, noisy VCFs (e.g. per-scaffold outputs) before merging
#' or importing into a GDS. For GDS-level MAC filtering based on genotypes,
#' see \code{\link{filter_ma}}.
#'
#' @param vcf (character)
#' Path to the input VCF (\code{.vcf} or \code{.vcf.gz}).
#'
#' @param mac.threshold (integer)
#' Minimum ALT allele count required to keep a variant.
#' Variants where all ALT alleles have \code{AC < mac.threshold} are removed.
#' Default: \code{mac.threshold = 4L}.
#'
#' @param out.vcf (optional, character)
#' Output VCF filename for MAC-filtered variants.
#' If \code{NULL}, the function generates a filename from the input VCF root by
#' appending \code{".mac.vcf"} and \code{".gz"} when \code{compress = TRUE}.
#' Default: \code{out.vcf = NULL}.
#'
#' @param blacklist.file (optional, character)
#' Path to a TSV file storing low-MAC variants (\code{CHROM}, \code{POS},
#' \code{ID}).
#' If \code{NULL}, the function generates a filename from the VCF root by
#' appending \code{".mac_blacklist.tsv"}.
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
#' @return
#' Invisibly returns a list with:
#' \itemize{
#'   \item \code{mac.vcf} – path to the MAC-filtered VCF;
#'   \item \code{blacklist} – path to the low-MAC variant list;
#'   \item \code{log} – path to the bcftools log file.
#' }
#'
#' @examples
#' \dontrun{
#' filter_mac_vcf(
#'   vcf            = "chr1_freebayes_poly.vcf.gz",
#'   mac.threshold  = 4L,
#'   out.vcf        = NULL,   # will append ".mac.vcf.gz" by default
#'   blacklist.file = NULL,
#'   bcftools.path  = "bcftools"
#' )
#' }
#'
#' @name filter_mac_vcf
#' @rdname filter_mac_vcf
#' @export
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
filter_mac_vcf <- function(
    vcf,
    mac.threshold  = 4L,
    out.vcf        = NULL,
    blacklist.file = NULL,
    bcftools.path  = "bcftools",
    compress       = TRUE,
    index          = TRUE,
    verbose        = TRUE
) {

  # ---- Basic checks -------------------------------------------------------
  if (!file.exists(vcf)) {
    rlang::abort(stringi::stri_join("VCF not found: ", vcf))
  }

  if (!is.numeric(mac.threshold) || length(mac.threshold) != 1L) {
    rlang::abort("`mac.threshold` must be a single numeric/integer value.")
  }
  mac.threshold <- as.integer(mac.threshold)

  # Ensure bcftools is available
  bcftools_require(bcftools.path = bcftools.path)

  # ---- Build bcftools expression ------------------------------------------
  # Keep variants where any ALT allele has AC >= mac.threshold
  mac_expr <- stringi::stri_join("INFO/AC>=", mac.threshold)

  # ---- Paths and filenames ------------------------------------------------
  vcf.base <- basename(vcf)
  vcf.dir  <- dirname(vcf)
  vcf.root <- sub("\\.(vcf|vcf\\.gz)$", "", vcf.base, ignore.case = TRUE)

  # MAC-filtered VCF path
  if (is.null(out.vcf)) {
    out.vcf <- paste0(
      vcf.root,
      ".mac.vcf",
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
    blacklist.file <- paste0(vcf.root, ".mac_blacklist.tsv")
  }
  blacklist.file <- file.path(vcf.dir, blacklist.file)

  # Log file (next to the VCF)
  log.file <- file.path(
    vcf.dir,
    paste0(
      vcf.root,
      ".bcftools_mac_",
      format(Sys.time(), "%Y%m%d@%H%M%S"),
      ".log"
    )
  )

  # ---- Run bcftools view (MAC filter) -------------------------------------
  args.view <- c(
    "view",
    "-i", mac_expr,
    if (compress) "-Oz" else "-Ov",
    "-o", out.vcf,
    vcf
  )

  if (verbose) {
    message("\nGenerating MAC-filtered VCF \u2192 ", out.vcf)
    message("MAC threshold (AC): ", mac.threshold)
  }

  res.view <- bcftools_exec(
    bcftools.path  = bcftools.path,
    args           = args.view,
    log.file       = log.file,
    label          = "MAC filter (view)",
    verbose        = verbose,
    fail_on_status = TRUE
  )

  # ---- Optional index -----------------------------------------------------
  if (compress && index) {
    args.idx <- c(
      "index",
      "-t",
      out.vcf
    )

    res.idx <- bcftools_exec(
      bcftools.path  = bcftools.path,
      args           = args.idx,
      log.file       = log.file,
      label          = "index",
      verbose        = verbose,
      fail_on_status = FALSE
    )

    if (!identical(res.idx$status, 0L)) {
      warning(
        "bcftools index failed for MAC-filtered VCF. See log: ",
        log.file,
        call. = FALSE
      )
    }
  }

  # ---- Query low-MAC variants (blacklist) ---------------------------------
  # bcftools query -f '%CHROM\t%POS\t%ID\n' -e 'INFO/AC>=mac' in.vcf
  fmt_arg <- "%CHROM\\t%POS\\t%ID\\n"

  args.query <- c(
    "query",
    "-f", fmt_arg,
    "-e", mac_expr,
    vcf
  )

  res.query <- bcftools_exec(
    bcftools.path  = bcftools.path,
    args           = args.query,
    log.file       = log.file,
    label          = "MAC blacklist (query)",
    verbose        = verbose,
    fail_on_status = TRUE
  )

  bl_txt <- res.query$stdout
  if (is.null(bl_txt) || !nzchar(bl_txt)) {
    writeLines(character(), blacklist.file)
  } else {
    # split on newlines in case stdout is one long scalar
    bl_vec <- unlist(strsplit(bl_txt, "\n", fixed = TRUE), use.names = FALSE)
    bl_vec <- bl_vec[nzchar(bl_vec)]
    writeLines(bl_vec, blacklist.file)
  }

  # ---- Return -------------------------------------------------------------
  invisible(list(
    mac.vcf   = out.vcf,
    blacklist = blacklist.file,
    log       = log.file
  ))
} # END filter_mac_vcf
