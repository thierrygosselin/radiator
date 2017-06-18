# Filter with a whitelist of individuals
#' @title Filter individuals with whitelist
#' @description This function use a whitelist of individuals to filter
#' a dataset.
#' @param data A data frame object or file (using ".tsv") with an 'INDIVIDUALS'
#' column.
#' @param whitelist.id The whitelist of individuals with no header column,
#' the first column is assumed to be the indivuals ID.
#' @param pop.id.start The start of your population id 
#' in the name of your individual sample.
#' @param pop.id.end The end of your population id 
#' in the name of your individual sample.
#' @param pop.levels A character string with your populations ordered.
#' @param filename The name of the file written to the directory.
#' @import dplyr
#' @import readr
#' @export
#' @rdname filter_whitelist_id

filter_id_whitelist <- function (data, whitelist.id, pop.id.start, pop.id.end, pop.levels, filename) {
  
  X1 <- NULL
  INDIVIDUALS <- NULL
  
  
  
  if (is.vector(whitelist.id) == "TRUE") {
    message("Using the whitelist in your directory")
    
    whitelist.id <- read_tsv(whitelist.id, col_names = F) %>%
      select(INDIVIDUALS=X1) %>%
      mutate(INDIVIDUALS = factor(INDIVIDUALS))
    
  } else {
    message("Using the whitelist from your global environment")
    
    whitelist.id <- whitelist.id %>%
      select(INDIVIDUALS) %>%
      mutate(INDIVIDUALS = factor(INDIVIDUALS)) %>%
      ungroup()
    
  }
  
  if (is.vector(data) == "TRUE") {
    message("Using the file in your directory")
    
    data <- read_tsv(data, col_names = T)
    #     data <- read_tsv(data, col_names = T) %>%
    #       mutate(INDIVIDUALS = as.character(INDIVIDUALS))
    
  } else {
    message("Using the file from your global environment")
    
    data <- data %>%
      #         mutate(INDIVIDUALS = as.character(INDIVIDUALS)) %>%
      ungroup()
  }
  
  whitelist.id.filter <- whitelist.id %>%
    left_join(data, by = "INDIVIDUALS") %>%
    mutate(
      POP_ID = factor(substr(INDIVIDUALS, pop.id.start, pop.id.end),
                      levels = pop.levels, ordered =T)
    ) %>%
    arrange(LOCUS)
  
  write.table(whitelist.id.filter, filename, sep = "\t", row.names = F,
              col.names = T, quote = F)
  
  invisible(cat(sprintf(
    "Whitelist individuals filter 
Filename:
%s\n
Written in the directory:
%s",
    filename, getwd()
  )))
  return(whitelist.id.filter)
}

