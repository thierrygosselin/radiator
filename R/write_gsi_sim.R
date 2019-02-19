# write a gsi_sim file

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

#' @seealso /href{https://github.com/eriqande/gsi_sim}{gsi_sim} and
#' /href{https://github.com/eriqande/rubias}{rubias}: genetic stock
#' identification (GSI) in the tidyverse.

#' @return A gsi_sim input file is saved to the working directory.
#' @export
#' @rdname write_gsi_sim

#' @importFrom tibble as_data_frame
#' @importFrom tidyr separate gather unite
#' @importFrom dplyr n_distinct rename mutate select left_join arrange
#' @importFrom stringi stri_replace_all_regex stri_join stri_replace_all_fixed
#' @importFrom readr write_delim

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

  # Checking for missing and/or default arguments ******************************
  if (missing(data)) rlang::abort("Input file necessary to write the gsi_sim file is missing")

  # POP_ID in gsi_sim does not like spaces, we need to remove space in everything touching POP_ID...
  # pop.levels, pop.labels, pop.select, strata, etc
  check <- check_pop_levels(pop.levels = pop.levels,
                            pop.labels = pop.labels,
                            pop.select = NULL)
  # list2env(x = ., globalenv())
  pop.levels <- check$pop.levels
  pop.labels <- check$pop.labels
  # check <- NULL

  # Import data
  if (is.vector(data)) {
    data <- radiator::tidy_wide(data = data, import.metadata = FALSE)
  } else {
    colnames(data) <- stringi::stri_replace_all_fixed(
      str = colnames(data),
      pattern = "GENOTYPE",
      replacement = "GT",
      vectorize_all = FALSE)

    # remove space in POP_ID
    data$POP_ID <- clean_pop_names(data$POP_ID)
    # necessary steps to make sure we work with unique markers and not duplicated LOCUS
    if (tibble::has_name(data, "LOCUS") && !tibble::has_name(data, "MARKERS")) {
      data <- dplyr::rename(.data = data, MARKERS = LOCUS)
    }
    # change sep in individual name
    data$INDIVIDUALS <- clean_ind_names(data$INDIVIDUALS)
  }

  # Info for gsi_sim input -----------------------------------------------------
  n.individuals <- dplyr::n_distinct(data$INDIVIDUALS)  # number of individuals
  n.markers <- dplyr::n_distinct(data$MARKERS)          # number of markers
  list.markers <- order(unique(data$MARKERS))           # list of markers

  # Spread/dcast in wide format ------------------------------------------------------
  data <- dplyr::select(data, MARKERS, POP_ID, INDIVIDUALS, GT) %>%
    tidyr::separate(data = ., col = GT, into = c("A1", "A2"), sep = 3, remove = TRUE) %>%
    # tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS, POP_ID)) %>%
    data.table::as.data.table(.) %>%
    data.table::melt.data.table(
      data = .,
      id.vars = c("MARKERS", "INDIVIDUALS", "POP_ID"),
      variable.name = "ALLELES",
      value.name = "GT"
    ) %>%
    tibble::as_data_frame(.) %>%
    dplyr::arrange(MARKERS) %>%
    tidyr::unite(col = MARKERS_ALLELES, MARKERS , ALLELES, sep = "_") %>%
    dplyr::arrange(POP_ID, INDIVIDUALS, MARKERS_ALLELES) %>%
    # dplyr::group_by(POP_ID, INDIVIDUALS) %>%
    # tidyr::spread(data = ., key = MARKERS_ALLELES, value = GT) %>%
    data.table::as.data.table(.) %>%
    data.table::dcast.data.table(
      data = .,
      formula = POP_ID + INDIVIDUALS ~ MARKERS_ALLELES,
      value.var = "GT"
    ) %>%
    tibble::as_data_frame(.) %>%
    dplyr::ungroup(.)

  # population levels and strata------------------------------------------------
  if (is.null(strata)) {
    strata <- dplyr::distinct(data, INDIVIDUALS, POP_ID)
  }

  strata.df <- read_strata(
    strata = strata,
    pop.id = TRUE,
    pop.levels = pop.levels,
    pop.labels = pop.labels,
    verbose = FALSE) %$%
    strata

  strata <- NULL #no longer needed

  if (tibble::has_name(data, "POP_ID")) data <- dplyr::select(.data = data, -POP_ID)
  data <- dplyr::left_join(data, strata.df, by = "INDIVIDUALS")

  # using pop.levels and pop.labels info if present
  data <- change_pop_names(
    data = data,
    pop.levels = pop.levels,
    pop.labels = pop.labels)

  # write gsi_sim file ---------------------------------------------------------
  # open the connection to the file
  filename.connection <- file(filename, "w")

  # Line 1: number of individuals and the number of markers
  writeLines(text = stringi::stri_join(n.individuals, n.markers, sep = " "),
             con = filename.connection, sep = "\n")

  # Line 2 and + : List of markers
  writeLines(text = stringi::stri_join(list.markers, sep = "\n"),
             con = filename.connection, sep = "\n")

  # close the connection to the file
  close(filename.connection) # close the connection

  # remaining lines, individuals and genotypes
  pop <- data$POP_ID # Create a vector with the population
  data <- suppressWarnings(dplyr::select(.data = data, -POP_ID))  # remove pop id
  gsi_sim.split <- split(data, pop)  # split gsi_sim by populations
  pop.string <- as.character(unique(pop))
  for (k in pop.string) {
    # utils::write.table(x = as.data.frame(stringi::stri_join("pop", k, sep = " ")), file = filename, append = TRUE, quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
    # utils::write.table(x = gsi_sim.split[[k]], file = filename, append = TRUE, quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
    readr::write_delim(x = as.data.frame(stringi::stri_join("pop", k, sep = " ")),
                       path = filename, delim = "\n", append = TRUE, col_names = FALSE)
    readr::write_delim(x = gsi_sim.split[[k]],
                       path = filename, delim = " ", append = TRUE, col_names = FALSE)
  }

  gsi_sim.split <- data <-pop <- pop.string <- NULL
  return(filename)
} # End write_gsi function

