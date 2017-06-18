# write a genepop file from a tidy data frame

#' @name write_genepop

#' @title Write a genepop file from a tidy data frame

#' @description Write a genepop file from a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator} 
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.

#' @param pop.levels (optional, string) A character string with your populations ordered.
#' Default: \code{pop.levels = NULL}.

#' @param genepop.header The first line of the Genepop file.
#' Default: \code{genepop.header = NULL} will use "radiator genepop with date".

#' @param markers.line (optional, logical) In the genepop and structure
#' file, you can write the markers on a single line separated by 
#' commas \code{markers.line = TRUE}, 
#' or have markers on a separate line, i.e. in one column, for the genepop file
#' (not very useful with thousands of markers) and not printed at all for the
#' structure file.
#' Default: \code{markers.line = TRUE}.

#' @param filename (optional) The file name prefix for the genepop file 
#' written to the working directory. With default: \code{filename = NULL}, 
#' the date and time is appended to \code{radiator_genepop_}.

#' @param ... other parameters passed to the function.

#' @return A genepop file is saved to the working directory. 

#' @export
#' @rdname write_genepop

#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join
#' @importFrom stringi stri_join stri_replace_all_fixed stri_extract_all_fixed stri_sub stri_pad_left
#' @importFrom purrr flatten_chr
#' @importFrom tidyr spread gather
#' @importFrom readr write_delim

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_genepop <- function(
  data,
  pop.levels = NULL, 
  genepop.header = NULL, 
  markers.line = TRUE, 
  filename = NULL,
  ...
) {
  
  # Checking for missing and/or default arguments ******************************
  if (missing(data)) stop("Input file missing")
  
  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    input <- radiator::tidy_wide(data = data, import.metadata = TRUE)
  } else {
    input <- data
  }
  
  # check genotype column naming
  colnames(input) <- stringi::stri_replace_all_fixed(
    str = colnames(input), 
    pattern = "GENOTYPE", 
    replacement = "GT", 
    vectorize_all = FALSE
  )
  
  # necessary steps to make sure we work with unique markers and not duplicated LOCUS
  if (tibble::has_name(input, "LOCUS") && !tibble::has_name(input, "MARKERS")) {
    input <- dplyr::rename(.data = input, MARKERS = LOCUS)
  }
  
  input <- dplyr::select(.data = input, POP_ID, INDIVIDUALS, MARKERS, GT)
  
  input <- input %>% 
    dplyr::mutate(
      GT = stringi::stri_replace_all_fixed(
        str = as.character(GT), 
        pattern = c("/", ":", "_", "-", "."), 
        replacement = "", 
        vectorize_all = FALSE),
      GT = stringi::stri_pad_left(str = as.character(GT), pad = "0", width = 6)
    )
  
  # pop.levels -----------------------------------------------------------------
  if (!is.null(pop.levels)) {
    input <- dplyr::mutate(
      .data = input,
      POP_ID = factor(POP_ID, levels = pop.levels, ordered = TRUE),
      POP_ID = droplevels(POP_ID)
    ) %>% 
      dplyr::arrange(POP_ID, INDIVIDUALS, MARKERS)
  } else {
    input <- dplyr::mutate(.data = input, POP_ID = factor(POP_ID)) %>% 
      dplyr::arrange(POP_ID, INDIVIDUALS, MARKERS)
  }
  
  # Create a marker vector  ------------------------------------------------
  markers <- dplyr::distinct(.data = input, MARKERS) %>%
    dplyr::arrange(MARKERS) %>%
    purrr::flatten_chr(.)
  
  # Wide format ----------------------------------------------------------------
  input <- input %>%
    dplyr::arrange(MARKERS) %>% 
    dplyr::group_by(POP_ID, INDIVIDUALS) %>%
    tidyr::spread(data = ., key = MARKERS, value = GT) %>% 
    dplyr::arrange(POP_ID, INDIVIDUALS) %>%
    dplyr::ungroup(.) %>% 
    dplyr::mutate(INDIVIDUALS = paste(INDIVIDUALS, ",", sep = ""))
  
  # Write the file in genepop format -------------------------------------------
  
  # Filename ------------------------------------------------------------------
  if (is.null(filename)) {
    # Get date and time to have unique filenaming
    file.date <- stringi::stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "", vectorize_all = FALSE)
    file.date <- stringi::stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
    file.date <- stringi::stri_sub(file.date, from = 1, to = 13)
    filename <- stringi::stri_join("radiator_genepop_", file.date, ".gen")
  } else {
    filename <- stringi::stri_join(filename, ".gen")
  }
  
  
  # genepop header  ------------------------------------------------------------
  if (is.null(genepop.header)) {
    # Get date and time to have unique filenaming
    file.date <- stringi::stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "", vectorize_all = FALSE)
    file.date <- stringi::stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
    file.date <- stringi::stri_sub(file.date, from = 1, to = 13)
    genepop.header <- stringi::stri_join("radiator genepop ", file.date)
  }
  
  # genepop construction
  pop <- input$POP_ID # Create a population vector
  input <- split(select(.data = input, -POP_ID), pop) # split genepop by populations
  filename.connection <- file(filename, "w") # open the connection to the file
  writeLines(text = genepop.header, con = filename.connection, sep = "\n") # write the genepop header
  if (markers.line) { # write the markers on a single line
    writeLines(text = stringi::stri_join(markers, sep = ",", collapse = ", "), con = filename.connection, sep = "\n") 
  } else {# write the markers on a single column (separate lines)
    writeLines(text = stringi::stri_join(markers, sep = "\n"), con = filename.connection, sep = "\n")
  }
  close(filename.connection) # close the connection
  for (i in 1:length(input)) {
    readr::write_delim(x = as.data.frame("pop"), path = filename, delim = "\n", append = TRUE, col_names = FALSE)
    readr::write_delim(x = input[[i]], path = filename, delim = " ", append = TRUE, col_names = FALSE)
  }
}# End write_genepop
