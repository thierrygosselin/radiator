# Fis
#' @title Fis filter
#' @description
#' TODO
#' @param approach Character. By \code{"SNP"} or by \code{"haplotype"}.
#' The function will consider the SNP or haplotype statistics to filter the marker.
#' Default by \code{"haplotype"}.
#' @param fis.min.threshold Number.
#' @param fis.max.threshold Number.
#' @param fis.diff.threshold Number (0 - 1)
#' @param pop.threshold Fixed number of pop required to keep the locus.
#' @param percent Is the threshold a percentage ? TRUE or FALSE.
#' @param filename (optional) The function uses \code{\link[fst]{write.fst}},
#' to write the tidy data frame in
#' the folder created in the working directory. The file extension appended to
#' the \code{filename} provided is \code{.rad}.
#' Default: \code{filename = NULL}.
#' @inheritParams tidy_genomic_data
#' @rdname filter_fis
#' @export

filter_fis <- function(data, approach = "haplotype", fis.min.threshold, fis.max.threshold, fis.diff.threshold, pop.threshold, percent, filename) {

  if (is.vector(data)) {
    data <- readr::read_tsv(data, col_names = T)
  }
  pop.number <- dplyr::n_distinct(data$POP_ID)

  if(stringi::stri_detect_fixed(pop.threshold, ".") & pop.threshold < 1) {
    multiplication.number <- 1/pop.number
    message("Using a proportion threshold...")
    threshold.id <- "of proportion"
  } else if (stringi::stri_detect_fixed(percent, "T")) {
    multiplication.number <- 100/pop.number
    message("Using a percentage threshold...")
    threshold.id <- "percent"
  } else {
    multiplication.number <- 1
    message("Using a fixed threshold...")
    threshold.id <- "population as a fixed"
  }

  if (missing(approach) | approach == "haplotype"){
    message("Approach selected: haplotype")
    fis.filter <- data %>%
      dplyr::select(LOCUS, POS, POP_ID, FIS) %>%
      dplyr::group_by (LOCUS, POP_ID) %>%
      dplyr::summarise(
        FIS_MIN = min(FIS),
        FIS_MAX = max(FIS),
        FIS_DIFF = FIS_MAX-FIS_MIN
      ) %>%
      dplyr::filter(FIS_MIN >= fis.min.threshold) %>%
      dplyr::filter(FIS_MAX <= fis.max.threshold) %>%
      dplyr::filter(FIS_DIFF <= fis.diff.threshold) %>%
      dplyr::group_by(LOCUS) %>%
      dplyr::tally(.) %>%
      dplyr::filter((n * multiplication.number) >= pop.threshold) %>%
      dplyr::select(LOCUS) %>%
      dplyr::left_join(data, by="LOCUS") %>%
      dplyr::arrange(LOCUS, POP_ID)
  } else {
    message("Approach selected: SNP")
    fis.filter <- data %>%
      dplyr::select(LOCUS, POS, POP_ID, FIS) %>%
      dplyr::group_by(LOCUS, POS, POP_ID) %>%
      dplyr::summarise(
        FIS_MIN = min(FIS),
        FIS_MAX = max(FIS),
        FIS_DIFF = FIS_MAX-FIS_MIN
      ) %>%
      dplyr::filter(FIS_MIN >= fis.min.threshold) %>%
      dplyr::filter(FIS_MAX <= fis.max.threshold) %>%
      dplyr::filter(FIS_DIFF <= fis.diff.threshold) %>%
      dplyr::group_by(LOCUS, POS) %>%
      dplyr::tally(.) %>%
      dplyr::filter((n * multiplication.number) >= pop.threshold) %>%
      dplyr::select(LOCUS, POS) %>%
      dplyr::left_join(data, by = c("LOCUS", "POS")) %>%
      dplyr::arrange(LOCUS, POS, POP_ID)
}


  if (missing(filename) == "FALSE") {
    message("Saving the file in your working directory...")
    readr::write_tsv(fis.filter, filename, append = FALSE, col_names = TRUE)
    saving <- paste("Saving was selected, the filename:", filename, sep = " ")
  } else {
    saving <- "Saving was not selected..."
  }


  invisible(cat(sprintf(
    "Fis filter:
Fis min >= %s or Fis max <= %s or
difference along the read/haplotype between the max and min Fis > %s,
all in %s percent of the sampling sites/pop were removed\n
The number of SNP removed by the Fis filter = %s SNP
The number of LOCI removed by the Fis filter = %s LOCI
The number of SNP before -> after the Fis filter: %s -> %s SNP
The number of LOCI before -> after the Fis filter: %s -> %s LOCI\n
%s\n
Working directory:
%s",
    fis.min.threshold, fis.max.threshold, fis.diff.threshold, pop.threshold,
    dplyr::n_distinct(data$POS) - dplyr::n_distinct(fis.filter$POS),
    dplyr::n_distinct(data$LOCUS) - dplyr::n_distinct(fis.filter$LOCUS),
    dplyr::n_distinct(data$POS), dplyr::n_distinct(fis.filter$POS),
    dplyr::n_distinct(data$LOCUS), dplyr::n_distinct(fis.filter$LOCUS),
    saving, getwd()
  )))
  return(fis.filter)
}
