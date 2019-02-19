#' @name read_blacklist_id
#' @title read_blacklist_id
#' @description Read a file or object with blacklisted individuals.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param blacklist.id (optional, path or object) A blacklist file in the working directory
#' or object in the global environment. The data frame
#' as 1 column (named \code{INDIVIDUALS}) and is filled with the individual IDs
#' The ids are cleaned with \code{\link{clean_ind_names}} for separators,
#' only \code{-} are tolerated. Duplicates are removed automatically.
#' Default: \code{blacklist.id = NULL}.

#' @inheritParams radiator_common_arguments

#' @rdname read_blacklist_id
#' @export
#' @return A tibble with column \code{INDIVIDUALS}.
#' @examples
#' \dontrun{
#' bl <- radiator::read_blacklist_id("blacklist.tsv")
#' }
read_blacklist_id <- function(blacklist.id = NULL, verbose = TRUE) {
  if (!is.null(blacklist.id)) {# With blacklist of ID
    if (is.vector(blacklist.id)) {
      suppressMessages(blacklist.id <- readr::read_tsv(
        blacklist.id,
        col_names = TRUE,
        col_types = readr::cols(.default = readr::col_character())))
    } else {
      if (!rlang::has_name(blacklist.id, "INDIVIDUALS")) {
        rlang::abort("Blacklist of individuals should have 1 column named: INDIVIDUALS")
      }
      blacklist.id <- dplyr::mutate_all(.tbl = blacklist.id, .funs = as.character)
    }
    blacklist.id$INDIVIDUALS <- radiator::clean_ind_names(blacklist.id$INDIVIDUALS)

    # remove potential duplicate id
    dup <- dplyr::distinct(.data = blacklist.id, INDIVIDUALS)
    blacklist.id.dup <- nrow(blacklist.id) - nrow(dup)
    if (blacklist.id.dup >1) {
      if (verbose) message("Duplicate id's in blacklist: ", blacklist.id.dup)
      blacklist.id <- dup
    }
    dup <- blacklist.id.dup <- NULL
    if (verbose) message("Number of individuals in blacklist: ", nrow(blacklist.id))
  }
  return(blacklist.id)
}#End read_blacklist_id
