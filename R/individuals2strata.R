# Make strata file from individuals

#' @name individuals2strata
#' @title Create a strata file from a list of individuals
#' @description If your individuals have a consistent naming scheme 
#' (e.g. SPECIES-POPULATION-MATURITY-YEAR-ID = CHI-QUE-ADU-2014-020), 
#' use this function to rapidly create a strata file. 
#' Several functions in \pkg{radiator} and \pkg{assigner} requires 
#' a \code{strata} argument, i.e. a data frame with the individuals and 
#' associated groupings. If you have already run 
#' \href{http://catchenlab.life.illinois.edu/stacks/}{stacks} on your data, 
#' the strata file is similar to a stacks `population map file`, make sure you 
#' have the required column names  (\code{INDIVIDUALS} and \code{STRATA}).

#' @param data A file or data frame object with individuals in a column. The 
#' column name is \code{INDIVIDUALS}.

#' @param strata.start (integer) The start of your strata id. See details for more info.

#' @param strata.end (integer) The end of your strata id. See details for more info.

#' @param filename (optional) The file name for the strata object if you
#' want to save it in the working directory.
#' Default: \code{filename = NULL}, the starta objct is in the global 
#' environment only (i.e. not written in the working directory).

#' @return a strata object and file, if requested. The file is tab delimited
#' with 2 columns named:
#' \code{INDIVIDUALS} and \code{STRATA}. 
#' The \code{STRATA} column can be any hierarchical grouping.


#' @details 
#' \code{strata.start} and \code{strata.end}
#' The info must be found within the name of your individual sample. If not, 
#' you'll have to create a strata file by hand, the old fashion way.
#' e.g. if your individuals are identified 
#' in this form : SPECIES-POPULATION-MATURITY-YEAR-ID = CHI-QUE-ADU-2014-020,
#' then, to have the population id in the \code{STRATA} column, 
#' \code{strata.start = 5} and \code{strata.end = 7}. 
#' The \code{STRATA} column can be any hierarchical grouping.

#' @export
#' @rdname individuals2strata
#' @importFrom stringi stri_sub
#' @importFrom dplyr mutate
#' @importFrom data.table fread
#' @importFrom readr write_tsv
#' @importFrom tibble as_data_frame


#' @examples
#' \dontrun{
#' strata.abalone <- individuals2strata(
#' data = "individuals.abalone.tsv",
#' strata.start = 5, 
#' strata.end = 7,
#' filename = "strata.abalone.tsv"
#' )
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

# required to pass the R CMD check and have 'no visible binding for global variable'
# if (getRversion() >= "2.15.1") {
#   utils::globalVariables(
#     c("DP", "AD", "vcf.headers", "GT_VCF", "INDIVIDUALS2", ""
#     )
#   )
# }

individuals2strata <- function(
  data, 
  strata.start,
  strata.end,
  filename = NULL
  ) {
  
  
  # Checking for missing and/or default arguments ******************************
  if (missing(data)) stop("Input file missing")
  if (missing(strata.start)) stop("strata.start argument missing")
  if (missing(strata.end)) stop("strata.end argument missing")
  
  if (is.vector(data)) { # for file in the working directory
    input <- data.table::fread(
      input = data, 
      sep = "\t",
      stringsAsFactors = FALSE, 
      header = TRUE, 
      showProgress = TRUE,
      verbose = FALSE, 
      data.table = FALSE
    )
  } else {# object in global environment
    input <- data
  }
  
  input <- tibble::as_data_frame(input) %>% 
    dplyr::mutate(
      INDIVIDUALS =  as.character(INDIVIDUALS),
      STRATA = stringi::stri_sub(str = INDIVIDUALS, from = strata.start, to = strata.end)
    )
  
  
  # Write to working directory
  if (!is.null(filename)) {
    message("Writing the strata object to the working directory: \n", filename)
    readr::write_tsv(x = input, path = filename, col_names = TRUE)
  }
  
  return(input)
} # end individuals2strata
