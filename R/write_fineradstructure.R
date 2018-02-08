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

#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join funs
#' @importFrom data.table as.data.table dcast.data.table
#' @importFrom stringi stri_join stri_replace_all_fixed
#' @importFrom tibble has_name as_data_frame
#' @importFrom adegenet genind

#' @references Malinsky M, Trucchi E, Lawson D, Falush D (2018)
#' RADpainter and fineRADstructure: population inference from RADseq data.
#' bioRxiv, 057711.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_fineradstructure <- function(data, filename = NULL) {

  # Checking for missing and/or default arguments ******************************
  if (missing(data)) stop("Input file missing")

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    data <- radiator::tidy_wide(data = data, import.metadata = FALSE)
  }

  #Check the required columns --------------------------------------------------
  if (!tibble::has_name(data, "GT_VCF_NUC")) {
    stop("Wrong genotype format: nucleotide information is required\ntry running the function radiator::tidy_genomic_data")
  }

  # Filename -------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (is.null(filename)) {
    filename <- stringi::stri_join("radiator_fineradstructure_", file.date, ".tsv")
  } else {
    filename.problem <- file.exists(filename)
    if (filename.problem) {
      filename <- stringi::stri_join(filename, "_fineradstructure_", file.date, ".tsv")
    } else {
      filename <- stringi::stri_join(filename, "_fineradstructure.tsv")
    }
  }

  # Conversion  ----------------------------------------------------------------
  want <- c("MARKERS", "INDIVIDUALS", "GT_VCF_NUC")
  data <- suppressWarnings(dplyr::select(data, dplyr::one_of(want))) %>%
    dplyr::mutate(
      GT_VCF_NUC = stringi::stri_replace_all_fixed(
        str = GT_VCF_NUC, pattern = "./.", replacement = NA, vectorize_all = FALSE)
      ) %>%
    data.table::as.data.table(.) %>%
    data.table::dcast.data.table(
      data = .,
      formula = MARKERS ~ INDIVIDUALS, value.var = "GT_VCF_NUC"
      ) %>%
    dplyr::select(-MARKERS) %>%
    readr::write_tsv(x = ., path = filename, na = "")
  return(data)
} # End write_fineradstructure
