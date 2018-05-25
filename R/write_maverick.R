# write a maverick file from a tidy data frame

#' @name write_maverick
#' @title Write a maverick file from a tidy data frame
#' @description Write a maverick file from a tidy data frame
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.

#' @param filename (optional) The file name prefix for the maverick file
#' written to the working directory. With default: \code{filename = NULL},
#' the date and time is appended to \code{radiator_maverick_}.

#' @param ... other parameters passed to the function.

#' @return A maverick file is saved to the working directory (with \code{.txt}).

#' @export
#' @rdname write_maverick
#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join
#' @importFrom stringi stri_join stri_replace_all_fixed stri_extract_all_fixed stri_sub
#' @importFrom purrr flatten_chr
#' @importFrom tidyr spread gather
#' @importFrom readr write_tsv

#' @references Verity R, Nichols RA (2016) Estimating the Number of
#' Subpopulations (K) in Structured Populations.
#' Genetics, 203, genetics.115.180992-1839.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

write_maverick <- function(
  data,
  filename = NULL,
  ...
) {

  # Checking for missing and/or default arguments ******************************
  if (missing(data)) stop("Input file missing")

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    data <- radiator::tidy_wide(data = data, import.metadata = FALSE)
  }

  # necessary steps to make sure we work with unique markers and not duplicated LOCUS
  if (tibble::has_name(data, "LOCUS") && !tibble::has_name(data, "MARKERS")) {
    data <- dplyr::rename(.data = data, MARKERS = LOCUS)
  }


  # Create a marker vector  ------------------------------------------------
  markers <- dplyr::distinct(.data = data, MARKERS) %>%
    dplyr::arrange(MARKERS) %>%
    purrr::flatten_chr(.)
  markers <- c("INDIVIDUALS", "POP_ID", markers)

  # maverick format ------------------------------------------------------------
  data <- dplyr::select(.data = data, POP_ID, INDIVIDUALS, MARKERS, GT) %>%
    tidyr::separate(col = GT, into = c("A1", "A2"), sep = 3, extra = "drop", remove = TRUE) %>%
    # tidyr::gather(data = ., key = ALLELES, value = GT, -c(POP_ID, INDIVIDUALS, MARKERS)) %>%
    data.table::as.data.table(.) %>%
      data.table::melt.data.table(
        data = .,
        id.vars = c("INDIVIDUALS", "POP_ID", "MARKERS"),
        variable.name = "ALLELES",
        value.name = "GT"
      ) %>%
      tibble::as_data_frame(.) %>%
      dplyr::mutate(
        GT = stringi::stri_replace_all_fixed(
          str = GT, pattern = "000", replacement = "-9", vectorize_all = FALSE),
        GT = as.integer(GT)
      ) %>%
      dplyr::select(INDIVIDUALS, POP_ID, MARKERS, ALLELES, GT) %>%
      dplyr::arrange(POP_ID, INDIVIDUALS, MARKERS, ALLELES) %>%
      # tidyr::spread(data = ., key = MARKERS, value = GT) %>%
      data.table::as.data.table(.) %>%
      data.table::dcast.data.table(
        data = .,
        formula = INDIVIDUALS + POP_ID + ALLELES ~ MARKERS,
        value.var = "GT"
      ) %>%
      tibble::as_data_frame(.) %>%
      dplyr::mutate(
        POP_ID = as.integer(POP_ID),
        ALLELES = NULL
        ) %>%
      dplyr::arrange(POP_ID, INDIVIDUALS)

  # Write the file in maverick format -----------------------------------------
  message("Generating MavericK output file")
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")

  # Filename
  if (is.null(filename)) {
    filename <- stringi::stri_join("radiator_maverick_", file.date, ".txt")
    parameter.filename <- stringi::stri_join("radiator_maverick_parameters", file.date, ".txt")
  } else {
    filename <- stringi::stri_join(filename, ".txt")
    parameter.filename <- stringi::stri_join(filename, "_parameters", ".txt")
  }
  p.v <- suppressWarnings(stringi::stri_join(packageVersion("radiator"), collapse = ""))
  maverick.header <- stringi::stri_join(
    "#MavericK input file generated with radiator v.", p.v,
    ", date and time: ", file.date, sep = "")

  filename.connection <- file(filename, "w") # open the connection to the file

  writeLines(text = maverick.header, con = filename.connection, sep = "\n")
  writeLines(text = stringi::stri_join(markers, sep = "\t", collapse = "\t"), con = filename.connection, sep = "\n")
  close(filename.connection) # close the connection
  readr::write_tsv(x = data, path = filename, append = TRUE, col_names = FALSE)

  # Parameter file -------------------------------------------------------------
  message("Generating MavericK parameter file")
  readr::write_lines(x = "\n#### Data properties", path = parameter.filename)
  tibble::tibble(KEY = c("headerRow_on", "popCol_on", "ploidyCol_on", "ploidy", "missingData"),
                              VALUE = c("t", "t", "f", "2", "-9")) %>%
    readr::write_tsv(x = ., path = parameter.filename, append = TRUE, col_names = FALSE)

  readr::write_lines(x = "\n\n#### Model parameters", path = parameter.filename, append = TRUE)
  tibble::tibble(KEY = c("Kmin", "Kmax", "admix_on", "fixAlpha_on", "alpha", "alphaPropSD"),
                              VALUE = c("3", "3", "t", "t", "1.0", "0.10")) %>%
    readr::write_tsv(x = ., path = parameter.filename, append = TRUE, col_names = FALSE)

  readr::write_lines(x = "\n\n#### Simulation parameters", path = parameter.filename, append = TRUE)
  tibble::tibble(KEY = c("exhaustive_on", "mainRepeats", "mainBurnin", "mainSamples",
                         "thermodynamic_on", "thermodynamicRungs",
                         "thermodynamicBurnin", "thermodynamicSamples",
                         "EMalgorithm_on", "EMrepeats", "EMiterations"),
                 VALUE = c("f", "1", "500", "5000", "f", "20", "500", "1000", "f", "100", "100")) %>%
    readr::write_tsv(x = ., path = parameter.filename, append = TRUE, col_names = FALSE)

  readr::write_lines(x = "\n\n#### Basic output properties", path = parameter.filename, append = TRUE)
  tibble::tibble(KEY = c("outputLog_on", "outputLikelihood_on", "outputQmatrix_ind_on",
                         "outputQmatrix_pop_on",
                         "outputQmatrixError_ind_on", "outputQmatrixError_pop_on",
                         "outputEvidence_on", "outputEvidenceNormalised_on",
                         "outputEvidenceDetails_on"),
                 VALUE = c("t", "t", "t", "t", "f", "f", "f", "f", "f")) %>%
    readr::write_tsv(x = ., path = parameter.filename, append = TRUE, col_names = FALSE)
return(NULL)
} # end write_maverick
