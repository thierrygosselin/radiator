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
#' @references Verity R, Nichols RA (2016) Estimating the Number of
#' Subpopulations (K) in Structured Populations.
#' Genetics, 203, genetics.115.180992-1839.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

write_maverick <- function(
  data,
  filename = NULL,
  ...
) {

  # for testing
  # filename = NULL

  # dotslist -------------------------------------------------------------------
  dotslist <- list(...)
  want <- c("subsample.markers", "whitelist.markers")
  unknowned_param <- setdiff(names(dotslist), want)

  if (length(unknowned_param) > 0) {
    rlang::abort("Unknowned \"...\" parameters ",
         stringi::stri_join(unknowned_param, collapse = " "))
  }

  radiator.dots <- dotslist[names(dotslist) %in% want]
  subsample.markers <- radiator.dots[["subsample.markers"]]
  whitelist.markers <- radiator.dots[["whitelist.markers"]]


  # Checking for missing and/or default arguments ******************************
  if (missing(data)) rlang::abort("Input file missing")

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    data <- radiator::tidy_wide(data = data, import.metadata = FALSE)
  }

  # necessary steps to make sure we work with unique markers and not duplicated LOCUS
  if (tibble::has_name(data, "LOCUS") && !tibble::has_name(data, "MARKERS")) {
    data <- dplyr::rename(.data = data, MARKERS = LOCUS)
  }



  # Create a marker vector  ----------------------------------------------------
  markers <- dplyr::distinct(.data = data, MARKERS) %>%
    dplyr::arrange(MARKERS) %>%
    purrr::flatten_chr(.)

  # subsample markers
  if (!is.null(subsample.markers)) {
    # subsample.markers <- NULL
    # subsample.markers <- 1000
    markers <- sample(x = markers, size = subsample.markers, replace = FALSE)
    data <- dplyr::filter(data, MARKERS %in% markers)
  }
  if (!is.null(whitelist.markers)) {
    if (is.vector(whitelist.markers)) {
      whitelist.markers <- readr::read_tsv(
        file = whitelist.markers,
        col_types = readr::cols(.default = readr::col_character()))
    }
    data <- dplyr::filter(data, MARKERS %in% whitelist.markers$MARKERS)
    markers <- dplyr::distinct(.data = data, MARKERS) %>%
      dplyr::arrange(MARKERS) %>%
      purrr::flatten_chr(.)
  }

  markers <- c("INDIVIDUALS", "POP_ID", markers)

  # Get info for parameter file
  n.pop <- dplyr::n_distinct(data$POP_ID)

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
      tibble::as_tibble(.) %>%
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
        formula = INDIVIDUALS + POP_ID ~ MARKERS + ALLELES,
        value.var = "GT"
      ) %>%
      tibble::as_tibble(.) %>%
      dplyr::mutate(POP_ID = as.integer(POP_ID)) %>%
      dplyr::arrange(POP_ID, INDIVIDUALS)

  # markers <- colnames
  # Write the file in maverick format -----------------------------------------
  message("Generating MavericK output file")
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")

  # Filename
  path <- getwd() # get working directory

  if (is.null(filename)) {
    folder.name <- stringi::stri_join("radiator_maverick_", file.date)
    path.folder.maverick <- file.path(path, folder.name)
  } else {
    folder.name <- stringi::stri_join(filename, file.date, sep = "_")
    path.folder.maverick <- file.path(path, folder.name)
  }
  dir.create(path.folder.maverick)
  message("Folder created: ", folder.name)
  output.folder <- file.path(path.folder.maverick, "ouput")
  input.folder <- file.path(path.folder.maverick, "input")
  dir.create(output.folder)
  dir.create(input.folder)

  filename <- stringi::stri_join(input.folder, "/data.txt")
  parameter.filename <- stringi::stri_join(input.folder, "/parameters.txt")



  p.v <- suppressWarnings(stringi::stri_join(utils::packageVersion("radiator"), collapse = ""))
  maverick.header <- stringi::stri_join(
    "#MavericK input file generated with radiator v.", p.v,
    ", date and time: ", file.date, sep = "")

  filename.connection <- file(filename, "w") # open the connection to the file

  writeLines(text = maverick.header, con = filename.connection, sep = "\n")
  # writeLines(text = stringi::stri_join(markers, sep = "\t", collapse = "\t"), con = filename.connection, sep = "\n")
  close(filename.connection) # close the connection
  readr::write_tsv(x = data, file = filename, append = TRUE, col_names = TRUE)

  # Parameter file -------------------------------------------------------------
  message("Generating MavericK parameter file")
  readr::write_lines(x = "\n#### Data properties", file = parameter.filename)
  tibble::tibble(KEY = c("headerRow_on", "popCol_on", "ploidyCol_on", "ploidy", "missingData", "dataFormat"),
                              VALUE = c("t", "t", "f", "2", "-9", "2")) %>%
    readr::write_tsv(x = ., file = parameter.filename, append = TRUE, col_names = FALSE)

  readr::write_lines(x = "\n\n#### Model parameters", file = parameter.filename, append = TRUE)
  tibble::tibble(KEY = c("Kmin", "Kmax", "admix_on", "fixAlpha_on", "alpha", "alphaPropSD"),
                              VALUE = c("1", n.pop + 2, "t", "t", "1.0", "0.10")) %>%
    readr::write_tsv(x = ., file = parameter.filename, append = TRUE, col_names = FALSE)

  readr::write_lines(x = "\n\n#### Simulation parameters", file = parameter.filename, append = TRUE)
  tibble::tibble(KEY = c("exhaustive_on", "mainRepeats", "mainBurnin", "mainSamples",
                         "thermodynamic_on", "thermodynamicRungs",
                         "thermodynamicBurnin", "thermodynamicSamples",
                         "EMalgorithm_on", "EMrepeats", "EMiterations"),
                 VALUE = c("f", "3", "500", "5000", "t", "20", "500", "1000", "t", "100", "100")) %>%
    readr::write_tsv(x = ., file = parameter.filename, append = TRUE, col_names = FALSE)

  readr::write_lines(x = "\n\n#### Basic output properties", file = parameter.filename, append = TRUE)
  tibble::tibble(KEY = c("outputLog_on", "outputLikelihood_on", "outputQmatrix_ind_on",
                         "outputQmatrix_pop_on",
                         "outputQmatrixError_ind_on", "outputQmatrixError_pop_on",
                         "outputEvidence_on", "outputEvidenceNormalised_on",
                         "outputEvidenceDetails_on"),
                 VALUE = c("t", "t", "t", "t", "t", "t", "t", "t", "t")) %>%
    readr::write_tsv(x = ., file = parameter.filename, append = TRUE, col_names = FALSE)

  readr::write_lines(x = "\n\n#### Additional arguments", file = parameter.filename, append = TRUE)
  # generate output folder
  tibble::tibble(KEY = c("outputRoot", "inputRoot"), VALUE = c(output.folder, input.folder)) %>%
    readr::write_tsv(x = ., file = parameter.filename, append = TRUE, col_names = FALSE)

return(NULL)
} # end write_maverick
