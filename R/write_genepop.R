# write a genepop file from a tidy data frame

#' @name write_genepop

#' @title Write a genepop file from a tidy data frame

#' @description Write a genepop file from a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.

#' @param pop.levels (optional, string) A character string with your populations ordered.
#' Default: \code{pop.levels = NULL}. Described in \code{\link{read_strata}}.

#' @param genepop.header The first line of the Genepop file.
#' Default: \code{genepop.header = NULL} will use "radiator genepop with date".

#' @param markers.line (optional, logical) In the genepop and structure
#' file, you can write the markers on a single line separated by
#' commas \code{markers.line = TRUE},
#' or have markers on a separate line, i.e. in one column, for the genepop file
#' (not very useful with thousands of markers) and not printed at all for the
#' structure file.
#' Default: \code{markers.line = TRUE}.

#' @param filename (optional) The file name prefix for the genepop file
#' written to the working directory. With default: \code{filename = NULL},
#' the date and time is appended to \code{radiator_genepop_}.

#' @param ... other parameters passed to the function.

#' @return A genepop file is saved to the working directory.

#' @export
#' @rdname write_genepop

#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join
#' @importFrom stringi stri_join stri_replace_all_fixed stri_extract_all_fixed stri_sub stri_pad_left
#' @importFrom purrr flatten_chr
#' @importFrom tidyr spread gather
#' @importFrom readr write_delim

#' @references Raymond M. & Rousset F, (1995).
#' GENEPOP (version 1.2): population genetics software for exact tests
#' and ecumenicism.
#' J. Heredity, 86:248-249
#' @references Rousset F.
#' genepop'007: a complete re-implementation of the genepop software
#' for Windows and Linux.
#' Molecular Ecology Resources.
#' 2008, 8: 103-106.
#' doi:10.1111/j.1471-8286.2007.01931.x

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_genepop <- function(
  data,
  pop.levels = NULL,
  genepop.header = NULL,
  markers.line = TRUE,
  filename = NULL,
  ...
) {


  # # For test
  # data
  # pop.levels = NULL
  # genepop.header = NULL
  # markers.line = TRUE
  # filename = NULL
  # Checking for missing and/or default arguments ******************************
  if (missing(data)) rlang::abort("Input file missing")

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    data <- radiator::tidy_wide(data = data, import.metadata = TRUE)
  }

  # necessary steps to make sure we work with unique markers and not duplicated LOCUS
  if (tibble::has_name(data, "LOCUS") && !tibble::has_name(data, "MARKERS")) {
    data <- dplyr::rename(.data = data, MARKERS = LOCUS)
  }

  data <- dplyr::select(.data = data, POP_ID, INDIVIDUALS, MARKERS, GT)

  data <- data %>%
    dplyr::mutate(
      GT = stringi::stri_replace_all_fixed(
        str = as.character(GT),
        pattern = c("/", ":", "_", "-", "."),
        replacement = "",
        vectorize_all = FALSE),
      GT = stringi::stri_pad_left(str = as.character(GT), pad = "0", width = 6)
    )

  # pop.levels -----------------------------------------------------------------
  if (!is.null(pop.levels)) {
    data <- dplyr::mutate(
      .data = data,
      POP_ID = factor(POP_ID, levels = pop.levels, ordered = TRUE),
      POP_ID = droplevels(POP_ID)
    ) %>%
      dplyr::arrange(POP_ID, INDIVIDUALS, MARKERS)
  } else {
    data <- dplyr::mutate(.data = data, POP_ID = factor(POP_ID)) %>%
      dplyr::arrange(POP_ID, INDIVIDUALS, MARKERS)
  }

  # Create a marker vector  ------------------------------------------------
  markers <- dplyr::distinct(.data = data, MARKERS) %>%
    dplyr::arrange(MARKERS) %>%
    purrr::flatten_chr(.)

  # Wide format ----------------------------------------------------------------
  data <- data %>%
    dplyr::arrange(MARKERS) %>%
    dplyr::group_by(POP_ID, INDIVIDUALS) %>%
    tidyr::spread(data = ., key = MARKERS, value = GT) %>%
    dplyr::arrange(POP_ID, INDIVIDUALS) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(INDIVIDUALS = paste(INDIVIDUALS, ",", sep = ""))

  # Write the file in genepop format -------------------------------------------
  # Date and time
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")

  # Filename -------------------------------------------------------------------
  if (is.null(filename)) {
    filename <- stringi::stri_join("radiator_genepop_", file.date, ".gen")
  } else {
    filename.problem <- file.exists(filename)
    if (filename.problem) {
      filename <- stringi::stri_join(filename, "_genepop_", file.date, ".gen")
    } else {
      filename <- stringi::stri_join(filename, "_genepop", ".gen")
    }
  }

  # genepop header  ------------------------------------------------------------
  if (is.null(genepop.header)) {
    genepop.header <- stringi::stri_join("radiator genepop ", file.date)
  }

  # genepop construction
  pop <- data$POP_ID # Create a population vector
  data <- split(select(.data = data, -POP_ID), pop) # split genepop by populations
  filename.connection <- file(filename, "w") # open the connection to the file
  writeLines(text = genepop.header, con = filename.connection, sep = "\n") # write the genepop header
  if (markers.line) { # write the markers on a single line
    writeLines(text = stringi::stri_join(markers, sep = ",", collapse = ", "), con = filename.connection, sep = "\n")
  } else {# write the markers on a single column (separate lines)
    writeLines(text = stringi::stri_join(markers, sep = "\n"), con = filename.connection, sep = "\n")
  }
  close(filename.connection) # close the connection
  for (i in 1:length(data)) {
    readr::write_delim(x = as.data.frame("pop"), path = filename, delim = "\n", append = TRUE, col_names = FALSE)
    readr::write_delim(x = data[[i]], path = filename, delim = " ", append = TRUE, col_names = FALSE)
  }
}# End write_genepop
