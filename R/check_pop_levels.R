# Check pop.levels and pop.labels

#' @name check_pop_levels
#' @title Check the use of pop.levels, pop.labels and pop.select arguments.
#' @description Check that pop.levels and pop.labels and pop.select arguments
#' are used correctly and that the values are cleaned for spaces.

#' @inheritParams tidy_genomic_data

#' @rdname check_pop_levels
#' @export
#'
#' @importFrom stringi stri_replace_all_regex stri_join stri_replace_all_fixed
#' @importFrom dplyr n_distinct

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


check_pop_levels <- function(pop.levels = NULL, pop.labels = NULL, pop.select = NULL) {

  # checks ---------------------------------------------------------------------
  # removing spaces in data$POP_ID, pop.levels and pop.labels
  if (!is.null(pop.levels) && is.null(pop.labels)) {
    pop.labels <-pop.levels <- clean_pop_names(pop.levels)
  }

  if (!is.null(pop.labels)) {
    if (is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")
    if (length(pop.labels) != length(pop.levels)) {
      stop("pop.levels and pop.labels with different length: check arguments")
    }
    pop.labels <- clean_pop_names(pop.labels)
  }
  if (!is.null(pop.select)) pop.select <- clean_pop_names(pop.select)
  return(res = list(pop.levels = pop.levels, pop.labels = pop.labels, pop.select = pop.select))
}# end function change_pop_names
