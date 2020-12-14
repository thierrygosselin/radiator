# write a structure file from a tidy data frame

#' @name write_structure
#' @title Write a structure file from a tidy data frame
#' @description Write a structure file from a tidy data frame
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.

#' @inheritParams read_strata

#' @param filename (optional) The file name prefix for the structure file
#' written to the working directory. With default: \code{filename = NULL},
#' the date and time is appended to \code{radiator_structure_}.

#' @param ... other parameters passed to the function.

#' @return A structure file is saved to the working directory.

#' @export
#' @rdname write_structure
#' @references Pritchard JK, Stephens M, Donnelly P. (2000)
#' Inference of population structure using multilocus genotype data.
#' Genetics. Genetics Society of America. 155: 945–959.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

write_structure <- function(
  data,
  pop.levels = NULL,
  filename = NULL,
  ...
) {

  # Checking for missing and/or default arguments ******************************
  if (missing(data)) rlang::abort("Input file missing")

  # Import data ---------------------------------------------------------------
  data %<>% radiator::tidy_wide(data = .) %>%
    dplyr::select(STRATA, INDIVIDUALS, MARKERS, GT)

  # pop.levels -----------------------------------------------------------------
  if (!is.null(pop.levels)) {
    data %<>%
      dplyr::mutate(
        STRATA = factor(STRATA, levels = pop.levels, ordered = TRUE),
        STRATA = droplevels(STRATA)
      ) %>%
      dplyr::arrange(STRATA, INDIVIDUALS, MARKERS)
  } else {
    data %<>%
      dplyr::mutate(STRATA = factor(STRATA)) %>%
      dplyr::arrange(STRATA, INDIVIDUALS, MARKERS)
  }

  # Create a marker vector  ------------------------------------------------
  markers <- dplyr::distinct(.data = data, MARKERS) %>%
    dplyr::arrange(MARKERS) %>%
    purrr::flatten_chr(.)

  # Structure format ----------------------------------------------------------------
  data %<>%
    radiator::separate_gt(x = ., gt = "GT", gather = TRUE, exclude = c("STRATA", "INDIVIDUALS", "MARKERS"), split.chunks = 1L) %>%
    dplyr::mutate(
      ALLELES = dplyr::recode(.x = ALLELES, "000" = "-9"),
      ALLELES = as.integer(ALLELES)
    ) %>%
    radiator::rad_wide(x = ., formula = "INDIVIDUALS + STRATA ~ MARKERS + ALLELES_GROUP", values_from = "ALLELES") %>%
    dplyr::mutate(STRATA = as.integer(STRATA)) %>%
    dplyr::arrange(STRATA, INDIVIDUALS)

  # Write the file in structure format -----------------------------------------

  # Filename
  if (is.null(filename)) {
    file.date <- format(Sys.time(), "%Y%m%d@%H%M")
    filename <- stringi::stri_join("radiator_structure_", file.date, ".str")
  } else {
    filename <- stringi::stri_join(filename, ".str")
  }

  filename.connection <- file(filename, "w") # open the connection to the file
  writeLines(text = stringi::stri_join(markers, sep = "\t", collapse = "\t"),
             con = filename.connection, sep = "\n")
  close(filename.connection) # close the connection
  readr::write_tsv(x = data, file = filename, append = TRUE, col_names = FALSE)
} # end write_structure
