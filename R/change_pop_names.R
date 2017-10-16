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
  if (!is.null(pop.levels) & is.null(pop.labels)) {
    pop.levels <- stringi::stri_replace_all_fixed(pop.levels, pattern = " ", replacement = "_", vectorize_all = FALSE)
    pop.labels <- pop.levels
  }

  if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")

  if (!is.null(pop.labels)) {
    if (length(pop.labels) != length(pop.levels)) stop("pop.labels and pop.levels must have the same length (number of groups)")
    pop.labels <- stringi::stri_replace_all_fixed(pop.labels, pattern = " ", replacement = "_", vectorize_all = FALSE)
  }

  # in the data
  data$POP_ID <- stringi::stri_replace_all_fixed(data$POP_ID, pattern = " ", replacement = "_", vectorize_all = FALSE)


  # the number of pop in data == pop.levels
  if (!is.null(pop.levels)) {
    if (dplyr::n_distinct(data$POP_ID) != length(pop.levels)) {
      stop("The number of groups in your input file doesn't match the number of groups in pop.levels")
    }
  }
  # convert POP_ID to factor and change names-----------------------------------

  if (is.null(pop.levels)) { # no pop.levels
    data$POP_ID <- factor(data$POP_ID)
  } else {# with pop.levels
    data$POP_ID <- factor(x = data$POP_ID, levels = pop.levels, ordered = TRUE)
    levels(data$POP_ID) <- pop.labels
  }
  return(data)
}# end function change_pop_names
