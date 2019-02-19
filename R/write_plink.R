# write a plinkfile from a tidy data frame

#' @name write_plink
#' @title Write a plink tped/tfam file from a tidy data frame

#' @description Write a plink file from a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.

#' @param filename (optional) The file name prefix for tped/tfam files
#' written to the working directory. With default: \code{filename = NULL},
#' the date and time is appended to \code{radiator_plink_}.

#' @export
#' @rdname write_plink

#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join
#' @importFrom stringi stri_join stri_replace_all_fixed stri_sub stri_pad_left
#' @importFrom tidyr spread gather unite
#' @importFrom readr write_delim

#' @references Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR,
#' Bender D, et al.
#' PLINK: a tool set for whole-genome association and population-based linkage
#' analyses.
#' American Journal of Human Genetics. 2007: 81: 559–575. doi:10.1086/519795


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

write_plink <- function(data, filename = NULL) {

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    data <- radiator::tidy_wide(data = data, import.metadata = TRUE)
  }

  tped <- data %>%
    dplyr::arrange(INDIVIDUALS) %>%
    dplyr::mutate(
      COL1 = rep("0", n()),
      COL3 = rep("0", n()),
      COL4 = rep("0", n())
    ) %>%
    dplyr::select(COL1, MARKERS, COL3, COL4, INDIVIDUALS, GT) %>%
    dplyr::mutate(
      A1 = stringi::stri_sub(str = GT, from = 1, to = 3),
      A2 = stringi::stri_sub(str = GT, from = 4, to = 6)
    ) %>%
    dplyr::select(-GT) %>%
    tidyr::gather(ALLELES, GENOTYPE, -c(COL1, MARKERS, COL3, COL4, INDIVIDUALS)) %>%
    dplyr::mutate(
      GENOTYPE = as.character(as.numeric(GENOTYPE)),
      GENOTYPE = stringi::stri_pad_left(GENOTYPE, width = 2, pad = "0")
    ) %>%
    dplyr::arrange(INDIVIDUALS, ALLELES) %>%
    tidyr::unite(INDIVIDUALS_ALLELES, INDIVIDUALS, ALLELES, sep = "_") %>%
    dplyr::group_by(COL1, MARKERS, COL3, COL4) %>%
    tidyr::spread(data = ., key = INDIVIDUALS_ALLELES, value = GENOTYPE) %>%
    dplyr::arrange(MARKERS)

  tfam <- dplyr::distinct(.data = data, POP_ID, INDIVIDUALS) %>%
    dplyr::arrange(INDIVIDUALS) %>%
    dplyr::mutate(
      COL3 = rep("0",n()),
      COL4 = rep("0",n()),
      COL5 = rep("0",n()),
      COL6 = rep("-9",n())
    )

  # Create a filename to save the output files ********************************
  if (is.null(filename)) {
    # Get date and time to have unique filenaming
    file.date <- stringi::stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "", vectorize_all = FALSE)
    file.date <- stringi::stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
    file.date <- stringi::stri_sub(file.date, from = 1, to = 13)
    filename.tped <- stringi::stri_join("radiator_plink_", file.date, ".tped")
    filename.tfam <- stringi::stri_join("radiator_plink_", file.date, ".tfam")
  } else {
    filename.tped <- stringi::stri_join(filename, ".tped")
    filename.tfam <- stringi::stri_join(filename, ".tfam")
  }
  readr::write_delim(x = tped, path = filename.tped, col_names = FALSE, delim = " ")
  readr::write_delim(x = tfam, path = filename.tfam, col_names = FALSE, delim = " ")
} # end write_plink
