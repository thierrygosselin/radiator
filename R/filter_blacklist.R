#' @title filter_loci_blacklist
#' @description This function use a blacklist of loci to filter
#' a dataset (sumstats, tidy vcf, haplotype file).
#' @param data A data frame object or file (using ".tsv").
#' @param blacklist The blacklist is a data frame or file (using ".tsv")
#' with one column named `LOCUS`.
#' @param type A description of the type of blacklist filtering used
#' (e.g. "paralog filter").
#' @import dplyr
#' @import readr
#' @export
#' @rdname filter_loci_blacklist

filter_loci_blacklist <- function (data, blacklist, type) {
  
  LOCUS <- NULL
  POP_ID <- NULL
  
  if (is.vector(blacklist) == "TRUE") {
      blacklist.loci <- read_tsv(blacklist, col_names = T)
      message("Using the blacklist in your directory")

    } else {
      blacklist.loci <- blacklist
      message("Using the blacklist from your global environment")

      }
  
  if (is.vector(data) == "TRUE") {
    data <- read_tsv(data, col_names = T)
    message("Using the file in your directory")

    } else {
      data <- data
      message("Using the file from your global environment")

      }
  
  blacklist.filter <- data %>%
  anti_join(blacklist.loci, by = "LOCUS") %>%
  arrange(LOCUS, POP_ID)
  
invisible(cat(sprintf(
"Blacklist %s filter:
The number of markers removed by the blacklist %s filter: SNP = %s, LOCI = %s
The number of markers before -> after the blacklist %s filter: %s -> %s SNP, %s -> %s LOI",
type, type, (n_distinct(data$POS)-n_distinct(blacklist.filter$POS)), (n_distinct(data$LOCUS)-n_distinct(blacklist.filter$LOCUS)), type, n_distinct(data$POS), n_distinct(blacklist.filter$POS), n_distinct(data$LOCUS), n_distinct(blacklist.filter$LOCUS)
)))
blacklist.filter
}
