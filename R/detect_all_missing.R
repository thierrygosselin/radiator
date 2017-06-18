# Detect markers with all missing genotypes

#' @name detect_all_missing

#' @title Detect markers with all missing genotypes

#' @description Detect if markers in tidy dataset have no genotypes at all.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator} 
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment.

#' @return The filtered dataset if problematic markers were found. Otherwise,
#' the untouch dataset.

#' @export
# @keywords internal
#' @rdname detect_all_missing
#' @importFrom dplyr select mutate group_by filter distinct

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


detect_all_missing <- function(data) {
  # markers with all missing... yes I've seen it... breaks code...
  marker.problem <- data %>% 
    dplyr::group_by(MARKERS) %>%
    dplyr::distinct(GT) %>%
    dplyr::summarise(GENOTYPED = length(GT[GT != "000000"])) %>% 
    dplyr::filter(GENOTYPED == 0) %>% 
    dplyr::ungroup(.)
  
  problem <- nrow(marker.problem)
  if (problem > 0) {
    message("Data set contains ", problem," marker(s) with no genotypes (all missing)...")
    message("    removing problematic markers...")
    data <- dplyr::filter(data, !MARKERS %in% marker.problem$MARKERS)
  }
  marker.problem <- problem <- NULL
  return(data)
}#End detect_all_missing
