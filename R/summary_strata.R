#' @title summary of strata
#' @description Summarise the information of a strata file or object.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.
#' @param strata (path or object) The strata file or object.
#' @rdname summary_strata
#' @export
#' @return
#' \enumerate{
#' \item Number of strata/populations
#' \item Number of individuals
#' \item Number of individuals per populations
#' \item Number of duplicate ids.
#' @examples
#' \dontrun{
#' radiator::summary_strata(strata)
#' }
summary_strata <- function(strata) {
  if (is.vector(strata)) {
    strata <- readr::read_tsv(file = strata)
  }
  colnames(strata) <- stringi::stri_replace_all_fixed(
    colnames(strata), "POP_ID", "STRATA",
    vectorize_all = FALSE)
  strata <- dplyr::select(strata, INDIVIDUALS, STRATA)
  strata.stats <- strata %>% dplyr::group_by(STRATA) %>% dplyr::tally(.) %>%
    dplyr::mutate(POP_IND = stringi::stri_join(STRATA, n, sep = " = "))

   duplicate.id <- nrow(strata) - length(unique(strata$INDIVIDUALS))

  message("Number of populations: ", dplyr::n_distinct(strata$STRATA))
  message("Number of individuals: ", dplyr::n_distinct(strata$INDIVIDUALS))
  message("Number of ind/pop:\n", stringi::stri_join(strata.stats$POP_IND, collapse ="\n"))
  message("\nNumber of duplicate id: ", duplicate.id)
}#End summary_strata
