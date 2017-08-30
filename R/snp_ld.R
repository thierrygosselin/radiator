#' @name snp_ld
#' @title GBS/RADseq short distance linkage disequilibrium pruning
#' @description SNP short distance linkage disequilibrium pruning. With anonymous markers from
#' RADseq/GBS de novo discovery, you can minimize linkage disequilibrium (LD) by
#' choosing among these 3 options (see argument below).
#'
#' For long distance linkage disequilibrium pruning, see details below.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data A tidy data set with ID (LOCUS) and POS (SNP) information.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.
#' Usually the LOCUS and POS info is taken from a VCF file.

#' @param snp.ld (character) 3 options:
#' \code{snp.ld = "random"} selection, \code{snp.ld = "first"} or
#' \code{snp.ld = "last"} for SNP on the same short read/haplotype.
#' Default: \code{snp.ld = "first"}.

#' @export
#' @rdname snp_ld


#' @importFrom stringi stri_replace_all_fixed stri_join
#' @importFrom tibble has_name
#' @importFrom dplyr select distinct group_by sample_n summarise semi_join n_distinct

snp_ld <- function(data, snp.ld = "first") {
  snp.ld <- match.arg(snp.ld, c("first", "random", "last"))


  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file missing")

  # Import data ---------------------------------------------------------------

  if (is.vector(data)) {
    input <- radiator::tidy_wide(data = data, import.metadata = FALSE)
  } else {
    input <- data
  }

  # check genotype column naming ---------------------------------------------
  if (tibble::has_name(input, "GENOTYPE")) {
    colnames(input) <- stringi::stri_replace_all_fixed(
      str = colnames(input),
      pattern = "GENOTYPE",
      replacement = "GT",
      vectorize_all = FALSE)
  }

  # Check that fiel format as ID and POS -------------------------------------
  if (!tibble::has_name(input, "LOCUS") && !tibble::has_name(input, "POS")) {
    stop("snp.ld is only available for VCF file and/or files with ID and POS info")
  }

  message("Minimizing short distance LD...")
  snp.locus <- dplyr::distinct(input, LOCUS, POS)

  locus.stats <- dplyr::group_by(.data = snp.locus, LOCUS) %>%
    dplyr::tally(.) %>% dplyr::rename(SNP_N = n) %>%
    dplyr::group_by(SNP_N) %>%
    dplyr::tally(.)

  if (nrow(locus.stats) > 1) {
    range.number.snp.locus <- range(locus.stats$SNP_N, na.rm = TRUE)
    message("    The range in the number of SNP/reads or locus is: ", stringi::stri_join(range.number.snp.locus, collapse = "-"))
  } else {
    message("    There is no variation in the number of SNP/reads or locus across the data")
  }

  # Random selection ---------------------------------------------------------
  if (snp.ld == "random") {
    snp.select <- snp.locus %>%
      dplyr::group_by(LOCUS) %>%
      dplyr::sample_n(tbl = ., size = 1, replace = FALSE)
    snp.before <- nrow(snp.locus)
    snp.after <- nrow(snp.select)
    message("    Number of SNP before = ", snp.before)
    message("    Number of SNP removed = ", snp.before - snp.after)
    message("    Number of SNP after = ", snp.after)
  }

  # Fist SNP on the read -----------------------------------------------------
  if (snp.ld == "first") {
    snp.select <- snp.locus %>%
      dplyr::group_by(LOCUS) %>%
      dplyr::summarise(POS = min(POS))
    snp.before <- nrow(snp.locus)
    snp.after <- nrow(snp.select)
    message("    Number of SNP before = ", snp.before)
    message("    Number of SNP removed = ", snp.before - snp.after)
    message("    Number of SNP after = ", snp.after)
  }


  # Last SNP on the read -----------------------------------------------------
  if (snp.ld == "last") {
    snp.select <- snp.locus %>%
      dplyr::group_by(LOCUS) %>%
      dplyr::summarise(POS = max(POS))
    snp.before <- nrow(snp.locus)
    snp.after <- nrow(snp.select)
    message("    Number of SNP before = ", snp.before)
    message("    Number of SNP removed = ", snp.before - snp.after)
    message("    Number of SNP after = ", snp.after)
  }

  # filtering the VCF to minimize LD -----------------------------------------
  input <- dplyr::semi_join(input, snp.select, by = c("LOCUS", "POS"))
  message("    Filtering the dataset to minimize LD by keeping only 1 SNP per short read/haplotype")
  return(input)
}#End snp_ld
