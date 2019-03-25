# write a fineradstructure file from a tidy data frame

#' @name write_fineradstructure
#' @title Write a fineRADstructure file from a tidy data frame

#' @description Write a fineRADstructure file from a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.

#' @param filename (optional) The file name prefix for the fineRADstructure file
#' written to the working directory. With default: \code{filename = NULL},
#' the date and time is appended to \code{radiator_fineradstructure_}.

#' @return A fineRADstructure file is written in the working directory.
#' An object is also returned in the global environment.

#' @export
#' @rdname write_fineradstructure

#' @references Malinsky M, Trucchi E, Lawson D, Falush D (2018)
#' RADpainter and fineRADstructure: population inference from RADseq data.
#' bioRxiv, 057711.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_fineradstructure <- function(data, filename = NULL) {

  # Checking for missing and/or default arguments ******************************
  if (missing(data)) rlang::abort("Input file missing")

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    data <- radiator::tidy_wide(data = data, import.metadata = FALSE)
  }

  #Check the required columns --------------------------------------------------
  if (!tibble::has_name(data, "GT_VCF_NUC")) {
    rlang::abort("Wrong genotype format: nucleotide information is required\ntry running the function radiator::tidy_genomic_data")
  }

  # Filename -------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (is.null(filename)) {
    filename <- stringi::stri_join("radiator_fineradstructure_", file.date, ".tsv")
    dic.filename <- stringi::stri_join("radiator_fineradstructure_dictionary_", file.date, ".tsv")
  } else {
    filename.problem <- file.exists(filename)
    if (filename.problem) {
      dic.filename <- stringi::stri_join(filename, "_fineradstructure_dictionary_", file.date, ".tsv")
      filename <- stringi::stri_join(filename, "_fineradstructure_", file.date, ".tsv")
    } else {
      dic.filename <- stringi::stri_join(filename, "_fineradstructure_dictionary.tsv")
      filename <- stringi::stri_join(filename, "_fineradstructure.tsv")
    }
  }

  # Conversion  ----------------------------------------------------------------
  # dictionnary pop and id
  dictionary <- dplyr::distinct(data, POP_ID, INDIVIDUALS) %>%
    dplyr::arrange(POP_ID, INDIVIDUALS) %>%
    dplyr::mutate(
      POP_ID = factor(POP_ID),
      NEW_POP = as.integer(POP_ID),
      NEW_POP = LETTERS[NEW_POP],
      INDIVIDUALS = factor(INDIVIDUALS),
      NEW_INDIVIDUALS = as.integer(INDIVIDUALS),
      ID = stringi::stri_join(NEW_POP, NEW_INDIVIDUALS)
    ) %>%
    readr::write_tsv(x = ., path = dic.filename)

  message("Dictionary filename for id and pop: ", dic.filename)


  want <- c("LOCUS", "MARKERS", "INDIVIDUALS", "GT_VCF_NUC")
  data <- suppressWarnings(
    dplyr::select(data, dplyr::one_of(want)) %>%
    dplyr::left_join(dplyr::select(dictionary, ID, INDIVIDUALS), by = "INDIVIDUALS")) %>%
    dplyr::mutate(
      INDIVIDUALS = NULL,
      GT_VCF_NUC = stringi::stri_replace_all_fixed(
        str = GT_VCF_NUC, pattern = "./.", replacement = NA, vectorize_all = FALSE)
    ) %>%
    separate_gt(x = ., sep = "/", gt = "GT_VCF_NUC",
                gather = TRUE, exclude = c("LOCUS", "MARKERS", "ID")) %>%
    dplyr::group_by(LOCUS, ID, ALLELE_GROUP) %>%
    dplyr::summarise(HAPLOTYPES = stringi::stri_join(HAPLOTYPES, collapse = "")) %>%
    dplyr::mutate(HAPLOTYPES = stringi::stri_replace_all_fixed(
      str = HAPLOTYPES, pattern = ".", replacement = "N", vectorize_all = FALSE)) %>%
    dplyr::arrange(LOCUS, ID, ALLELE_GROUP) %>%
    dplyr::group_by(LOCUS, ID) %>%
    dplyr::summarise(HAPLOTYPES = stringi::stri_join(HAPLOTYPES, collapse = "/")) %>%
    dplyr::ungroup(.) %>%
    data.table::as.data.table(.) %>%
    data.table::dcast.data.table(
      data = .,
      formula = LOCUS ~ ID, value.var = "HAPLOTYPES"
    ) %>%
    dplyr::arrange(LOCUS) %>%
    dplyr::select(-LOCUS) %>%
    readr::write_tsv(x = ., path = filename, na = "")
  return(data)
} # End write_fineradstructure
