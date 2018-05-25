# write a structure file from a tidy data frame

#' @name write_structure
#' @title Write a structure file from a tidy data frame
#' @description Write a structure file from a tidy data frame
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.


#' @param pop.levels (optional, string) A character string with your populations ordered.
#' Default: \code{pop.levels = NULL}.

#' @param markers.line (optional, logical) In the structure
#' file, you can write the markers on a single line separated by
#' commas \code{markers.line = TRUE},
#' or have markers on a separate line, i.e. in one column, for the structure file
#' (not very useful with thousands of markers) and not printed at all for the
#' structure file.
#' Default: \code{markers.line = TRUE}.


#' @param filename (optional) The file name prefix for the structure file
#' written to the working directory. With default: \code{filename = NULL},
#' the date and time is appended to \code{radiator_structure_}.

#' @param ... other parameters passed to the function.

#' @return A structure file is saved to the working directory.

#' @export
#' @rdname write_structure
#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join
#' @importFrom stringi stri_join stri_replace_all_fixed stri_extract_all_fixed stri_sub
#' @importFrom purrr flatten_chr
#' @importFrom tidyr spread gather
#' @importFrom readr write_tsv


#' @references Pritchard JK, Stephens M, Donnelly P. (2000)
#' Inference of population structure using multilocus genotype data.
#' Genetics. Genetics Society of America. 155: 945–959.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

write_structure <- function(
  data,
  pop.levels = NULL,
  markers.line = TRUE,
  filename = NULL,
  ...
) {

  # Checking for missing and/or default arguments ******************************
  if (missing(data)) stop("Input file missing")

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    input <- radiator::tidy_wide(data = data, import.metadata = FALSE)
  } else {
    input <- data
  }

  # necessary steps to make sure we work with unique markers and not duplicated LOCUS
  if (tibble::has_name(input, "LOCUS") && !tibble::has_name(input, "MARKERS")) {
    input <- dplyr::rename(.data = input, MARKERS = LOCUS)
  }


  input <- dplyr::select(.data = input, POP_ID, INDIVIDUALS, MARKERS, GT)

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

  # Structure format ----------------------------------------------------------------
  input <- input %>%
    tidyr::separate(col = GT, into = c("A1", "A2"), sep = 3, extra = "drop", remove = TRUE) %>%
    tidyr::gather(data = ., key = ALLELES, value = GT, -c(POP_ID, INDIVIDUALS, MARKERS)) %>%
    dplyr::mutate(
      GT = stringi::stri_replace_all_fixed(str = GT, pattern = "000", replacement = "-9", vectorize_all = FALSE),
      GT = as.integer(GT)
    ) %>%
    dplyr::select(INDIVIDUALS, POP_ID, MARKERS, ALLELES, GT) %>%
    tidyr::spread(data = ., key = MARKERS, value = GT) %>%
    dplyr::mutate(POP_ID = as.integer(POP_ID)) %>%
    dplyr::select(-ALLELES) %>%
    dplyr::arrange(POP_ID, INDIVIDUALS)

  # Write the file in structure format -----------------------------------------

  # Filename
  if (is.null(filename)) {
    file.date <- format(Sys.time(), "%Y%m%d@%H%M")
    filename <- stringi::stri_join("radiator_structure_", file.date, ".str")
  } else {
    filename <- stringi::stri_join(filename, ".str")
  }

  filename.connection <- file(filename, "w") # open the connection to the file
  writeLines(text = stringi::stri_join(markers, sep = "\t", collapse = "\t"),
             con = filename.connection, sep = "\n")
  close(filename.connection) # close the connection
  readr::write_tsv(x = input, path = filename, append = TRUE, col_names = FALSE)
} # end write_structure
