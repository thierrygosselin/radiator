# write a related file from a tidy data frame

#' @name write_related

#' @title Write a related file from a tidy data frame

#' @description Write a related file from a tidy data frame.
#' This output file format enables to run the data in the
#' \href{https://github.com/timothyfrasier/related}{related} R package
#' (Pew et al. 2015), which is essantially the R version of
#' \href{https://www.zsl.org/science/software/coancestry}{COANCESTRY} fortran
#' program developed by Jinliang Wang.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.

#' @param filename (optional) The file name prefix for the related file
#' written to the working directory. With default: \code{filename = NULL},
#' the date and time is appended to \code{radiator_related_}.
#' @inheritParams tidy_genomic_data

#' @param ... other parameters passed to the function.

#' @return A related file is saved to the working directory.

#' @export
#' @rdname write_related
#' @references Pew J, Muir PH, Wang J, Frasier TR (2015)
#' related: an R package for analysing pairwise relatedness from codominant
#' molecular markers.
#' Molecular Ecology Resources, 15, 557-561.
#' @references Wang, J. 2011.
#' COANCESTRY: A program for simulating, estimating and analysing relatedness
#' and inbreeding coefficients.
#' Molecular Ecology Resources 11(1): 141-145.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_related <- function(
  data,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  ...
) {

  # Checking for missing and/or default arguments ******************************
  if (missing(data)) rlang::abort("Input file missing")

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    data <- radiator::tidy_wide(data = data, import.metadata = TRUE)
  }

  # Format for related package -------------------------------------------------
  split.chunks <- 1L
  if (parallel.core > 1) split.chunks <- 3L

  data %<>%
    dplyr::select(INDIVIDUALS, MARKERS, GT) %>%
    separate_gt(
      x = .,
      gt = "GT",
      gather = TRUE,
      exclude = c("MARKERS", "INDIVIDUALS"),
      split.chunks = split.chunks
      ) %>%
    dplyr::mutate(
      MARKERS_ALLELES = stringi::stri_join(MARKERS, ALLELES_GROUP, sep = "."),
      MARKERS = NULL, ALLELES_GROUP = NULL
    ) %>%
    # dplyr::arrange(MARKERS_ALLELES, INDIVIDUALS) %>%
    rad_wide(x = ., formula = "INDIVIDUALS ~ MARKERS_ALLELES", values_from = "ALLELES")

  # Write the file in related format -------------------------------------------
  # Date and time
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")

  # Filename -------------------------------------------------------------------
  if (is.null(filename)) {
    filename <- stringi::stri_join("radiator_related_", file.date, ".txt")
  } else {
    filename.problem <- file.exists(filename)
    if (filename.problem) {
      filename <- stringi::stri_join(filename, "_related_", file.date, ".txt")
    } else {
      filename <- stringi::stri_join(filename, "_related", ".txt")
    }
  }
  readr::write_delim(x = data, file = filename, delim = " ", col_names = FALSE)
}# End write_related
