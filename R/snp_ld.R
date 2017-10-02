#' @name snp_ld
#' @title GBS/RADseq short distance linkage disequilibrium pruning
#' @description SNP short distance linkage disequilibrium pruning.
#' With anonymous markers from
#' RADseq/GBS de novo discovery, you can minimize linkage disequilibrium (LD) by
#' choosing among these 4 options (see argument below).
#'
#' For long distance linkage disequilibrium pruning, see details below.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data A tidy data set with ID (LOCUS) and POS (SNP) information.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.
#' Usually the LOCUS and POS info is taken from a VCF file.

#' @param snp.ld (character) 4 options:
#' \code{snp.ld = "random"} for a random selection of 1 SNP on the read,
#' \code{snp.ld = "first"} for the first one on the read...,
#' \code{snp.ld = "last"} for the last SNP on the read and
#' \code{snp.ld = "middle"} for locus with > 2 SNPs/read the option to select at random
#' one SNP between the first and the last SNP on the read. If the locus as <= 2
#' SNPs on the read, the first one is selected. Note that for that last option,
#' the numbers are reported.
#' Default: \code{snp.ld = "first"}.

#' @export
#' @rdname snp_ld


#' @importFrom stringi stri_replace_all_fixed stri_join
#' @importFrom tibble has_name
#' @importFrom dplyr select distinct group_by sample_n summarise semi_join n_distinct

snp_ld <- function(data, snp.ld = "first") {
  snp.ld <- match.arg(snp.ld, c("first", "random", "last", "middle"))

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
  message("    snp.ld = ", snp.ld)

  snp.locus <- dplyr::distinct(input, LOCUS, POS) %>% dplyr::arrange(LOCUS, POS)

  locus.stats <- dplyr::group_by(.data = snp.locus, LOCUS) %>%
    dplyr::tally(.) %>%
    dplyr::rename(SNP_N = n) %>%
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

  # Middle SNP on the read -----------------------------------------------------

  if (snp.ld == "middle") {
    snp.locus.prep <- dplyr::group_by(.data = snp.locus, LOCUS) %>%
      dplyr::tally(.)

    pick.middle <- snp.locus.prep %>%
      dplyr::filter(n > 2) %>%
      dplyr::select(LOCUS)

    if (nrow(pick.middle) == 0) {
      message("IMPORTANT: the data doesn't have more than 3 SNPs per read")
      message("    First SNP will be selected instead...")
      snp.select <- snp.locus %>%
        dplyr::group_by(LOCUS) %>%
        dplyr::summarise(POS = min(POS))
    } else {

      # For locus with <= 2 SNPs/read just keep the first one.
      keep.first <- snp.locus.prep %>%
        dplyr::filter(n <= 2) %>%
        dplyr::select(LOCUS)
      message("    Number of locus with first SNP selected: ", nrow(keep.first))
      keep.first.select <- snp.locus %>%
        dplyr::filter(LOCUS %in% keep.first$LOCUS) %>%
        dplyr::group_by(LOCUS) %>%
        dplyr::summarise(POS = min(POS))

      pick.middle.select <- snp.locus %>%
        dplyr::filter(LOCUS %in% pick.middle$LOCUS) %>%
        dplyr::group_by(LOCUS) %>%
        dplyr::filter(POS != min(POS)) %>% # remove the first SNP
        dplyr::filter(POS != max(POS)) %>% # remove the last SNP
        dplyr::sample_n(tbl = ., size = 1, replace = FALSE) # pick one at random

      message("    Number of locus with random middle SNP selected: ", nrow(pick.middle))
      snp.select <- dplyr::bind_rows(keep.first.select, pick.middle.select) %>%
        dplyr::arrange(LOCUS, POS)
    }
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
