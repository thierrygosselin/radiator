# write_gsi_sim ----------------------------------------------------------------

#' @name write_gsi_sim

#' @title Write a gsi_sim file from a data frame (wide or long/tidy).

#' @description Write a gsi_sim file from a data frame (wide or long/tidy).
#' Used internally in \href{https://github.com/thierrygosselin/assigner}{assigner} and
#' \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data A tidy genomic data set in the working directory tidy formats.
#' \emph{How to get a tidy data frame ?}
#' Look for \pkg{radiator} \code{\link{tidy_genomic_data}}.

#' @param pop.levels (option, string) This refers to the levels in a factor. In this
#' case, the id of the pop.
#' Use this argument to have the pop ordered your way instead of the default
#' alphabetical or numerical order. e.g. \code{pop.levels = c("QUE", "ONT", "ALB")}
#' instead of the default \code{pop.levels = c("ALB", "ONT", "QUE")}.
#' Default: \code{pop.levels = NULL}. If you find this too complicated, there is also the
#' \code{strata} argument that can do the same thing, see below.

#' @param pop.labels (optional, string) Use this argument to rename/relabel
#' your pop or combine your pop. e.g. To combine \code{"QUE"} and \code{"ONT"}
#' into a new pop called \code{"NEW"}:
#' (1) First, define the levels for your pop with \code{pop.levels} argument:
#' \code{pop.levels = c("QUE", "ONT", "ALB")}.
#' (2) then, use \code{pop.labels} argument:
#' \code{pop.levels = c("NEW", "NEW", "ALB")}.#'
#' To rename \code{"QUE"} to \code{"TAS"}:
#' \code{pop.labels = c("TAS", "ONT", "ALB")}.
#' Default: \code{pop.labels = NULL}. If you find this too complicated, there is also the
#' \code{strata} argument that can do the same thing, see below.
#'
#' @param strata (optional) A tab delimited file with 2 columns with header:
#' \code{INDIVIDUALS} and \code{STRATA}.
#' Default: \code{strata = NULL}. Use this argument to rename or change
#' the populations id with the new \code{STRATA} column.
#' The \code{STRATA} column can be any hierarchical grouping.

#' @param filename The name of the file written to the working directory.
# @param ... other parameters passed to the function in the future

#' @seealso \href{https://github.com/eriqande/gsi_sim}{gsi_sim} and
#' \href{https://github.com/eriqande/rubias}{rubias}: genetic stock
#' identification (GSI) in the tidyverse.

#' @return A gsi_sim input file is saved to the working directory.
#' @export
#' @rdname write_gsi_sim
#' @references Anderson, Eric C., Robin S. Waples, and Steven T. Kalinowski. (2008)
#' An improved method for predicting the accuracy of genetic stock identification.
#' Canadian Journal of Fisheries and Aquatic Sciences 65, 7:1475-1486.
#' @references Anderson, E. C. (2010) Assessing the power of informative subsets of
#' loci for population assignment: standard methods are upwardly biased.
#' Molecular ecology resources 10, 4:701-710.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_gsi_sim <- function(
  data,
  pop.levels = NULL,
  pop.labels = NULL,
  strata = NULL,
  filename = "gsi_sim.unname.txt"
) {

  ## TEST
  # pop.levels = NULL
  # pop.labels = NULL
  # strata = NULL
  # filename = "gsi_sim.unname.txt"

  # Checking for missing and/or default arguments
  if (missing(data)) rlang::abort("Input file necessary to write the gsi_sim file is missing")

  # POP_ID in gsi_sim does not like spaces, we need to remove space in everything touching POP_ID...
  # pop.levels, pop.labels, pop.select, strata, etc
  check <- radiator::check_pop_levels(pop.levels = pop.levels,
                            pop.labels = pop.labels,
                            pop.select = NULL)
  # list2env(x = ., globalenv())
  pop.levels <- check$pop.levels
  pop.labels <- check$pop.labels
  # check <- NULL

  # Import data
  data %<>% radiator::tidy_wide(data = .)
  if (rlang::has_name(data, "STRATA") & !rlang::has_name(data, "POP_ID")) {
    data %<>% dplyr::rename(POP_ID = STRATA)
  }

  # Info for gsi_sim input
  n.individuals <- dplyr::n_distinct(data$INDIVIDUALS)  # number of individuals
  n.markers <- dplyr::n_distinct(data$MARKERS)          # number of markers
  list.markers <- sort(unique(data$MARKERS))           # list of markers

  if (!rlang::has_name(data, "GT")) {
    data %<>% radiator::calibrate_alleles(data = ., gt = TRUE) %$% input
  }

  # Spread/dcast in wide format
  data %<>%
    dplyr::select(MARKERS, POP_ID, INDIVIDUALS, GT) %>%
    dplyr::mutate(
      A1 = stringi::stri_sub(str = GT, from = 1, to = 3),
      A2 = stringi::stri_sub(str = GT, from = 4, to = 6),
      GT = NULL
    ) %>%
    tidyr::pivot_longer(
      cols = -c("POP_ID", "INDIVIDUALS", "MARKERS"),
      names_to = "ALLELES",
      values_to = "GT"
    ) %>%
    dplyr::arrange(MARKERS) %>%
    dplyr::mutate(
      MARKERS_ALLELES = stringi::stri_join(MARKERS , ALLELES, sep = "_"),
      MARKERS = NULL, ALLELES = NULL
    ) %>%
    dplyr::arrange(POP_ID, INDIVIDUALS, MARKERS_ALLELES) %>%
    tidyr::pivot_wider(
      data = .,
      id_cols = c("POP_ID", "INDIVIDUALS"),
      names_from = MARKERS_ALLELES,
      values_from = "GT"
    ) %>%
    dplyr::ungroup(.)

  # population levels and strata
  if (is.null(strata)) {
    strata <- radiator::generate_strata(data = data, pop.id = TRUE)
  } else {
    strata <- radiator::read_strata(
      strata = strata,
      pop.id = TRUE,
      pop.levels = pop.levels,
      pop.labels = pop.labels,
      verbose = FALSE) %$%
      strata

    data <- radiator::join_strata(data = data, strata = strata, pop.id = TRUE, verbose = FALSE)
  }

  # write gsi_sim file
  filename.connection <- file(filename, "w") # open the connection

  # Line 1: number of individuals and the number of markers
  writeLines(text = stringi::stri_join(n.individuals, n.markers, sep = " "),
             con = filename.connection, sep = "\n")

  # Line 2 and + : List of markers
  writeLines(text = stringi::stri_join(list.markers, sep = "\n"),
             con = filename.connection, sep = "\n")
  close(filename.connection) # close the connection

  # remaining lines, individuals and genotypes
  pop <- data$POP_ID # Create a vector with the population
  suppressWarnings(data %<>% dplyr::select(-POP_ID))  # remove pop id
  gsi_sim.split <- split(data, pop)  # split gsi_sim by populations
  pop.string <- as.character(unique(pop))
  for (k in pop.string) {
    readr::write_delim(x = tibble::as_tibble(stringi::stri_join("pop", k, sep = " ")),
                       file = filename, delim = "\n", append = TRUE, col_names = FALSE)
    readr::write_delim(x = gsi_sim.split[[k]],
                       file = filename, delim = " ", append = TRUE, col_names = FALSE)
  }

  gsi_sim.split <- data <- pop <- pop.string <- NULL
  return(filename)
} # End write_gsi function

