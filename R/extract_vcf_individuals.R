#' @title extract_vcf_individuals
#' @description Function that returns the individuals presente in a vcf file.
#' then estimate allele frequencies, and the genotyping error rate.
#' Useful to create a strata file to work with function using vcf and that require
#' a strata file. This will ensure that you have the same individuals in your vcf file and strata...
#' @param data The vcf file.
#' @rdname extract_vcf_individuals
#' @export
#' @return A tibble with a column: \code{INDIVIDUALS}.
#' @author Thierry Gosselin \email{thierrygosselin@icloud.com}
extract_vcf_individuals <- function(data) {
  temp.file <- suppressWarnings(suppressMessages(readr::read_table(file = data, n_max = 200, col_names = "HEADER")))
  skip.number <- which(stringi::stri_detect_fixed(str = temp.file$HEADER,
                                                  pattern = "#CHROM")) - 1
  temp.file <- NULL
  remove <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
  id <- tibble::data_frame(INDIVIDUALS = colnames(readr::read_tsv(
    file = data,
    n_max = 1,
    skip = skip.number,
    col_types = readr::cols(.default = readr::col_character())) %>%
      dplyr::select(-dplyr::one_of(remove))))
  return(id)
}#End extract_vcf_individuals
