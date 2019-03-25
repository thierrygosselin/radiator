#' @title Clean marker's names for radiator and other packages
#' @description Function to clean marker's name
#' of weird separators
#' that interfere with some packages
#' or codes. \code{/}, \code{:}, \code{-} and \code{.} are changed to an underscore
#' \code{_}.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.
#' @param x (character string) Markers character string.
#' @rdname clean_markers_names
#' @export
clean_markers_names <- function(x) {
  x <- stringi::stri_replace_all_fixed(
    str = as.character(x),
    pattern = c("/", ":", "-", "."),
    replacement = "_",
    vectorize_all = FALSE)
}#End clean_markers_names

#' @title Clean individual's names for radiator and other packages
#' @description function to clean individual's name
#' that interfere with some packages
#' or codes. \code{_} and \code{:} are changed to a dash \code{-}.
#' Whitespaces are removed.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.
#' @param x (character string) Individuals character string.
#' @rdname clean_ind_names
#' @export
clean_ind_names <- function(x) {
  x <- stringi::stri_replace_all_fixed(
    str = as.character(x),
    pattern = c("_", ":", " "),
    replacement = c("-", "-", ""),
    vectorize_all = FALSE)
}#End clean_ind_names

#' @title Clean population's names for radiator and other packages
#' @description Function to clean pop's name
#' that interfere with some packages
#' or codes. Space is changed to an underscore \code{_}.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.
#' @param x (character string) Population character string.
#' @rdname clean_pop_names
#' @export
clean_pop_names <- function(x) {
clean_pop <- function(x) {
  stringi::stri_replace_all_fixed(
    str = as.character(x),
    pattern = " ",
    replacement = "_",
    vectorize_all = FALSE)
}

  if (is.factor(x)) {
    pop.levels <- clean_pop(as.character(levels(x)))
  } else {
    pop.levels <- clean_pop(as.character(unique(x)))
  }
  x <- clean_pop(as.character(x))
  x <- factor(x, levels = pop.levels)
}#End clean_pop_names
