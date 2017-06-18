# Minor Allele Frequency Diagnostic
#' @title MAF diagnostic
#' @description Minor Allele Frequency diagnostic, help choose a filter threshold.
#' @param data A file in the working directory or object in the global environment 
#' in wide or long (tidy) formats. To import, the function uses 
#' \href{https://github.com/thierrygosselin/radiator}{radiator} 
#' \code{\link[radiator]{tidy_wide}}.
#' \emph{See details of this function for more info}.
#' @param group.rank (Number) The number of group to class the MAF.
#' @param filename (optional) Name of the file written to the working directory.
#' @rdname diagnostic_maf
#' @export
#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join funs ntile
#' @importFrom readr read_tsv write_tsv
#' @details Highly recommended to look at the distribution of MAF
#' \link{plot_density_distribution_maf}

#' @examples
#' \dontrun{
#' problem <- radiator::diagnostic_maf(
#' data = tidy.salmon.data, group.rank = 10)
#' }


#' @seealso \link{filter_maf}


diagnostic_maf <- function(data, group.rank, filename = NULL){
  
  if (missing(data)) stop("Input file missing")
  if (missing(group.rank)) stop("group.rank argument missing")
  
  data <- radiator::tidy_wide(data = data, import.metadata = TRUE)
  
  # if (is.vector(data)) {
  #   data <- radiator::tidy_wide(data = data, import.metadata = TRUE)
  #   message("Using the file in your directory")
  #   
  # } else {
  #   data <- data
  #   message("Using the file from your global environment")
  # }
  
  # Check if allele frequency column is found
  if (!tibble::has_name(data, "FREQ_ALT")) {
    data <- radiator::allele_frequencies(data, verbose = TRUE)$freq.long %>% 
      dplyr::rename(FREQ_ALT = MAF_LOCAL, GLOBAL_MAF = MAF_GLOBAL)
  }
  
  # Local 
  test.local <- data %>%
    dplyr::select(LOCUS, POP_ID, FREQ_ALT) %>%
    dplyr::group_by(LOCUS, POP_ID) %>%
    dplyr::summarise(
      MAF_P = min(FREQ_ALT, na.rm = TRUE)
    ) %>% 
    dplyr::group_by(LOCUS) %>%
    dplyr::summarise(MAF_L = mean(MAF_P, na.rm = TRUE)) %>%
    dplyr::group_by(RANK = dplyr::ntile(MAF_L, group.rank)) %>%
    dplyr::summarise(
      LOCAL_MAF = mean(MAF_L, na.rm = T),
      n = length(LOCUS)
    ) %>%
    dplyr::select(-n)
  
  # Global
  
  test.global <- data %>%
    dplyr::select(LOCUS, POP_ID, GLOBAL_MAF) %>%
    dplyr::group_by(LOCUS) %>%
    dplyr::summarise(GLOBAL_MAF = mean(GLOBAL_MAF, na.rm = TRUE)) %>%
    dplyr::group_by(RANK = dplyr::ntile(GLOBAL_MAF, group.rank)) %>%
    dplyr::summarise(
      GLOBAL_MAF = mean(GLOBAL_MAF, na.rm = T),
      n = length(LOCUS)
    ) %>%
    dplyr::select(GLOBAL_MAF, n)
  
  maf.diagnostic <- dplyr::bind_cols(test.local, test.global)
  
  if (is.null(filename)) {
    message("Saving file: not selected")
  } else {
    message("Saving file: selected")
    readr::write_tsv(maf.diagnostic, filename, append = FALSE, col_names = TRUE)
    message("filename in the working directory: ", filename)
  }
  return(maf.diagnostic)
}
