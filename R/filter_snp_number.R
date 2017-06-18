# SNP number per haplotype
#' @title SNP number filter
#' @description Filter the dataset (sumstats or tidy vcf) to remove loci with more than x SNP.
#' @param data A data frame object or file (using ".tsv")
#' of class sumstats or a tidy VCF summarised.
#' @param max.snp.number Number
#' @param pop.threshold Fixed number of pop required to keep the locus.
#' @param percent Is the threshold a percentage ? TRUE or FALSE.
#' @param filename Name of the file written to the working directory (optional).
#' @rdname filter_snp_number
#' @export
#' @import stringi
#' @import dplyr
#' @import readr



filter_snp_number <- function(data, max.snp.number, pop.threshold, percent, filename) {
  
  
  
  POP_ID <- NULL
  LOCUS <- NULL
  SNP_N <- NULL
  
  if (is.vector(data) == "TRUE") {
    data <- read_tsv(data, col_names = T)
    message("Using the file in your directory")
  } else {
    data <- data
    message("Using the file from your global environment")
  }
  
  pop.number <- n_distinct(data$POP_ID)
  
  if(stri_detect_fixed(pop.threshold, ".") & pop.threshold < 1) {
    multiplication.number <- 1/pop.number
    message("Using a proportion threshold...")
    threshold.id <- "of proportion"
  } else if (stri_detect_fixed(percent, "T")) {
    multiplication.number <- 100/pop.number
    message("Using a percentage threshold...")
    threshold.id <- "percent"
  } else {
    multiplication.number <- 1
    message("Using a fixed threshold...")
    threshold.id <- "population as a fixed"
  }
  
  snp.filter <- data %>%
    group_by (LOCUS, POP_ID) %>%
    summarise(SNP_N = n_distinct(POS)) %>%
    filter(SNP_N <= max.snp.number) %>%
    group_by(LOCUS) %>%
    tally() %>%
    filter((n * multiplication.number) >= pop.threshold) %>%
    select(LOCUS) %>%
    left_join(data, by="LOCUS") %>%
    arrange(LOCUS, POP_ID)
  
  
  
  if (missing(filename) == "FALSE") {
    message("Saving the file in your working directory...")
    write_tsv(snp.filter, filename, append = FALSE, col_names = TRUE)
    saving <- paste("Saving was selected, the filename:", filename, sep = " ")
  } else {
    saving <- "Saving was not selected..."
  }
  
  
  invisible(cat(sprintf(
    "SNP number per haplotype filter:
max number of SNP allowed = %s in %s percent 
of the sampling sites/pop\n
The number of SNP removed by the SNP number filter = %s SNP
The number of LOCI removed by the SNP number filter = %s LOCI
The number of SNP before -> after the SNP number filter: %s -> %s SNP
The number of LOCI before -> after the SNP number filter: %s -> %s LOCI\n
%s\n
Working directory:
%s",
    max.snp.number, pop.threshold,
    n_distinct(data$POS)-n_distinct(snp.filter$POS),
    n_distinct(data$LOCUS)-n_distinct(snp.filter$LOCUS),
    n_distinct(data$POS), n_distinct(snp.filter$POS),
    n_distinct(data$LOCUS), n_distinct(snp.filter$LOCUS),
    saving, getwd()
  )))
  
  return(snp.filter)
  
}
