# write a fineradstructure file from a tidy data frame

#' @name write_fineradstructure
#' @title Write a fineRADstructure file from a tidy data frame

#' @description Write a \href{https://github.com/millanek/fineRADstructure}{fineRADstructure}
#' file from a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.

#' @inheritParams radiator_common_arguments
#' @inheritParams read_strata

#' @param filename (optional) The file name prefix for the fineRADstructure file
#' written to the working directory. With default: \code{filename = NULL},
#' the date and time is appended to \code{radiator_fineradstructure_}.

#' @return A fineRADstructure file is written in the working directory.
#' An object is also returned in the global environment.

#' @details fineRADstructure requires special formating for populations,
#' it needs to be ONLY LETTERS, not numbers. This information is then merged
#' with the sample id (converted to integers). radiator will generate a
#' dictionary file.


#' @section Life cycle:
#'
#' It become increasingly difficult for me to follow all the different naming
#' schemes researcher uses, if they're is any strategy... Consequently, I have
#' abandoned the idea of formating with letters the populations
#' to generate the fineRADstrucure file. If you get an error see the
#' the details section.

#' @seealso \href{https://github.com/millanek/fineRADstructure}{fineRADstructure}
#' @export
#' @rdname write_fineradstructure

#' @references Malinsky, M., Trucchi, E., Lawson, D., Falush, D. (2018).
#' RADpainter and fineRADstructure: population inference from RADseq data.
#' Mol. Biol. Evol.  35(5), 1284-1290.
#' https://dx.doi.org/10.1093/molbev/msy023

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_fineradstructure <- function(data, strata = NULL, filename = NULL) {

  # test
  # strata = NULL
  # filename = NULL


  # Checking for missing and/or default arguments ******************************
  if (missing(data)) rlang::abort("Input file missing")

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) data <- radiator::tidy_wide(data = data, import.metadata = FALSE)

  # Strata checks --------------------------------------------------------------
  if (is.null(strata)) {
    strata <- radiator::generate_strata(data = data, pop.id = TRUE)
  } else {
    strata <- radiator::read_strata(strata = strata, pop.id = TRUE) %$%
      strata
  }

  # unique(strata$POP_ID)
  message("Inspecting strata...")

  pop.string.1 <- stringi::stri_sub(str = strata$POP_ID, from = 1, to = 1) %>%
    unique

  check1 <- FALSE %in%
    (
      pop.string.1 %>%
        stringi::stri_detect_regex(str = ., pattern = "[:alpha:]") %>%
        unique
    )


  if (check1) {
    problem.in.pop <- pop.string.1 %>%
      stringi::stri_extract_all_regex(
        str = .,
        pattern = "[^a-zA-Z]",
        simplify = FALSE
      ) %>%
      unique %>%
      unlist %>%
      stringi::stri_omit_empty_na(x = .) %>%
      stringi::stri_join(., collapse = ", ")
    message("\n\nProblematic strata starting with: ", problem.in.pop)
    rlang::abort("Read the function documentation...twice! ?radiator::write_fineradstructure")
  }

  #Check the required columns --------------------------------------------------
  if (!tibble::has_name(data, "GT_VCF_NUC")) {
    message("\n\nfineRADstructure is best used with HAPLOTYPES information")
    message("In radiator the format is identified as column name: GT_VCF_NUC")
    rlang::abort("Wrong genotype format: nucleotide information in genotypes required")
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
  message("Generating a dictionary to go back to original ids...")
  dictionary <- strata %>%
    dplyr::arrange(POP_ID, INDIVIDUALS) %>%
    dplyr::mutate(
      # NEW_POP = as.integer(POP_ID),
      # NEW_POP = LETTERS[NEW_POP],
      ID = stringi::stri_join(POP_ID, as.integer(factor(INDIVIDUALS)))
    ) %>%
    readr::write_tsv(x = ., file = dic.filename)

  message("File written: ", dic.filename)

  want <- c("LOCUS", "MARKERS", "INDIVIDUALS", "GT_VCF_NUC")
  data %<>%
    dplyr::select(tidyselect::any_of(want)) %>%
    dplyr::left_join(dplyr::select(dictionary, ID, INDIVIDUALS), by = "INDIVIDUALS") %>%
    dplyr::mutate(
      INDIVIDUALS = NULL,
      GT_VCF_NUC = replace(x = GT_VCF_NUC, which(GT_VCF_NUC == "./."), NA)
    ) %>%
    separate_gt(
      x = .,
      gt = "GT_VCF_NUC",
      gather = TRUE,
      exclude = c("LOCUS", "MARKERS", "ID")
    ) %>%
    dplyr::group_by(LOCUS, ID, ALLELES_GROUP) %>%
    dplyr::summarise(ALLELES = stringi::stri_join(ALLELES, collapse = "")) %>%
    dplyr::mutate(
      ALLELES = stringi::stri_replace_all_fixed(
        str = ALLELES,
        pattern = ".",
        replacement = "N",
        vectorize_all = FALSE
      )
    ) %>%
    dplyr::arrange(LOCUS, ID, ALLELE_GROUP) %>%
    dplyr::group_by(LOCUS, ID) %>%
    dplyr::summarise(ALLELES = stringi::stri_join(ALLELES, collapse = "/")) %>%
    dplyr::ungroup(.) %>%
    rad_wide(x = ., formula = "LOCUS ~ ID", values_from = "ALLELES") %>%
    dplyr::arrange(LOCUS) %>%
    dplyr::select(-LOCUS)
  readr::write_tsv(x = data, file = filename, na = "")
  message("File written: ", filename)
  return(data)
} # End write_fineradstructure
