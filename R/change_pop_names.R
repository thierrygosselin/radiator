# change POP_ID names

#' @name change_pop_names
#' @title Transform into a factor the POP_ID column, change names and reorder the levels
#' @description Transform into a factor the POP_ID column, change names and reorder the levels

#' @inheritParams tidy_genomic_data

#' @rdname change_pop_names
#' @export
#'
#' @importFrom stringi stri_replace_all_regex stri_join stri_replace_all_fixed
#' @importFrom dplyr n_distinct

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


change_pop_names <- function(data, pop.levels = NULL, pop.labels = NULL) {

  # checks ---------------------------------------------------------------------
  if (missing(data)) stop("Input file missing")

  # POP_ID in gsi_sim does not like spaces, we need to remove space in everything touching POP_ID...

  # removing spaces in data$POP_ID, pop.levels and pop.labels
  if (!is.null(pop.levels)) {
    if (is.null(pop.labels)) {
      pop.labels <-pop.levels <- clean_pop_names(pop.levels)
    }
    if (dplyr::n_distinct(data$POP_ID) != length(pop.levels)) {
      stop("The number of strata/POP_ID in the data is different than the number of pop.levels: check argument and data")
    }
  }

  if (!is.null(pop.labels)) {
    if (is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")
    if (length(pop.labels) != length(pop.levels)) stop("pop.levels and pop.labels with different length: check arguments")
    pop.labels <- clean_pop_names(pop.labels)
  }

  # in the data
  data$POP_ID <- clean_pop_names(data$POP_ID)

  # convert POP_ID to factor and change names-----------------------------------

  if (is.null(pop.levels)) { # no pop.levels
    data$POP_ID <- factor(data$POP_ID)
  } else {# with pop.levels
    data$POP_ID <- factor(x = data$POP_ID, levels = pop.levels, ordered = FALSE)
    levels(data$POP_ID) <- pop.labels
  }
  data <- dplyr::arrange(data, POP_ID, INDIVIDUALS)
  return(data)
}# end function change_pop_names
