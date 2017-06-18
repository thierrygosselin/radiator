# write a hierfstat file from a tidy data frame

#' @name write_hierfstat
#' @title Write a hierfstat file from a tidy data frame

#' @description Write a hierfstat file from a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.

#' @param filename (optional) The file name prefix for the hierfstat file
#' written to the working directory. With default: \code{filename = NULL},
#' the date and time is appended to \code{radiator_hierfstat_}.

#' @return A hierfstat file is saved to the working directory.



#' @export
#' @rdname write_hierfstat
#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join
#' @importFrom stringi stri_join stri_replace_all_fixed stri_sub
#' @importFrom purrr flatten_chr
#' @importFrom tidyr spread unite
#' @importFrom readr write_delim

#' @references Goudet, J. (1995) FSTAT (Version 1.2): A computer program to
#' calculate F- statistics. Journal of Heredity, 86, 485-486.
#' @references Goudet, J. (2005) hierfstat, a package for r to compute and test hierarchical F-statistics. Molecular Ecology Notes, 5, 184-186.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_hierfstat <- function(data, filename = NULL) {

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file necessary to write the hierfstat file is missing")

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

  input <- dplyr::select(.data = input, POP_ID, INDIVIDUALS, MARKERS, GT) %>%
    dplyr::arrange(MARKERS, POP_ID, INDIVIDUALS)

  # Create a marker vector  ------------------------------------------------
  markers <- dplyr::distinct(.data = input, MARKERS) %>%
    dplyr::arrange(MARKERS) %>%
    purrr::flatten_chr(.)

  # Get the number of sample (pop) for hierfstat -------------------------------
  np <- nlevels(droplevels(input$POP_ID))
  np.message <- stringi::stri_join("    * Number of sample pop, np = ", np, sep = "")
  message(np.message)

  # Get the number of loci -----------------------------------------------------
  nl <- length(markers)
  nl.message <- stringi::stri_join("    * Number of markers, nl = ", nl, sep = "")
  message(nl.message)

  input <- input %>%
    dplyr::select(MARKERS, POP_ID, INDIVIDUALS, GT) %>%
    dplyr::mutate(
      GT = replace(GT, which(GT == "000000"), NA),
      A1 = as.numeric(stringi::stri_sub(str = GT, from = 1, to = 3)),
      A2 = as.numeric(stringi::stri_sub(str = GT, from = 4, to = 6))
    ) %>%
    dplyr::select(-GT)

  # Get the highest number used to label an allele -----------------------------
  nu <- max(c(unique(input$A1), unique(input$A2)), na.rm = TRUE)
  nu.message <- stringi::stri_join("    * The highest number used to label an allele, nu = ",
                           nu, sep = "")
  message(nu.message)

  # prep the data  -------------------------------------------------------------
  input <- suppressWarnings(
    tidyr::unite(data = input, GT, A1, A2, sep = "") %>%
      dplyr::mutate(GT = as.numeric(GT)) %>%
      dplyr::group_by(POP_ID, INDIVIDUALS) %>%
      tidyr::spread(data = ., MARKERS, GT) %>%
      dplyr::ungroup(.) %>%
      dplyr::arrange(POP_ID, INDIVIDUALS) %>%
      dplyr::mutate(POP_ID = as.integer(POP_ID)) %>%
      dplyr::select(-INDIVIDUALS)
  )

  # allele coding --------------------------------------------------------------
  allele.coding <- 1
  message("    * The alleles are encoded with one digit number")

  # Filename -------------------------------------------------------------------
  if (is.null(filename)) {
    # Get date and time to have unique filenaming
    file.date <- stringi::stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "", vectorize_all = FALSE)
    file.date <- stringi::stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
    file.date <- stringi::stri_sub(file.date, from = 1, to = 13)
    filename <- stringi::stri_join("radiator_hierfstat_", file.date, ".dat")
  } else {
    filename <- stringi::stri_join(filename, ".dat")
  }


  # FSTAT: write the first line ------------------------------------------------
  fstat.first.line <- stringi::stri_join(np, nl, nu, allele.coding, sep = " ")
  fstat.first.line <- as.data.frame(fstat.first.line)
  readr::write_delim(x = fstat.first.line, path = filename, delim = "\n", append = FALSE,
              col_names = FALSE)

  # FSTAT: write the locus name to the file
  loci.table <- as.data.frame(markers)
  readr::write_delim(x = loci.table, path = filename, delim = "\n", append = TRUE,
              col_names = FALSE)

  # FSTAT: write the pop and genotypes
  readr::write_delim(x = input, na = "00", path = filename, delim = "\t", append = TRUE,
              col_names = FALSE)
  input <- as.data.frame(input) # required by hierfstat...
  return(input)
}# End write_hierfstat
