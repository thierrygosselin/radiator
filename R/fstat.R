# tidy_fstat -------------------------------------------------------------------
#' @name tidy_fstat
#' @title fstat file to tidy dataframe

#' @description Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.
#' The function \code{tidy_fstat} reads a file in the
#' \href{http://www2.unil.ch/popgen/softwares/fstat.htm}{fstat} file (Goudet, 1995)
#' into a wide or long/tidy data frame
#'
#' To manipulate and prune the dataset prior to tidying, use the functions
#' \code{\link[radiator]{tidy_genomic_data}} and
#' \code{\link[radiator]{genomic_converter}}, that uses blacklist and whitelist along
#' several other filtering options.

#' @param data A \href{http://www2.unil.ch/popgen/softwares/fstat.htm}{fstat}
#' filename with extension \code{.dat}.

#' @param strata (optional) A tab delimited file with 2 columns. Header:
#' \code{INDIVIDUALS} and \code{STRATA}.
#' The \code{STRATA} column can be any hierarchical grouping.
#' To create a strata file see \code{\link[radiator]{individuals2strata}}.
#' Default: \code{strata = NULL}.

#' @param tidy (optional, logical) With \code{tidy = FALSE},
#' the markers are the variables and the genotypes the observations (wide format).
#' With the default: \code{tidy = TRUE}, markers and genotypes are variables
#' with their own columns (long format).

#' @param filename (optional) The file name for the tidy data frame
#' written to the working directory.
#' Default: \code{filename = NULL}, the tidy data is
#' in the global environment only (i.e. not written in the working directory).

#' @return The output in your global environment is a wide or long/tidy data frame.
#' If \code{filename} is provided, the wide or long/tidy data frame is also
#' written to the working directory.

#' @export
#' @rdname tidy_fstat

#' @examples
#' \dontrun{
#' # We will use the fstat dataset provided with adegenet package
#' require("hierfstat")
#'
#' # The simplest form of the function:
#' fstat.file <- radiator::tidy_fstat(
#'     data = system.file(
#'     "extdata/diploid.dat",
#'     package = "hierfstat"
#'     )
#'  )
#'
#' # To output a data frame in wide format, with markers in separate columns:
#' nancycats.wide <- radiator::tidy_fstat(
#'     data = system.file(
#'         "extdata/diploid.dat",
#'         package = "hierfstat"
#'     ),
#' tidy = FALSE
#' )
#' }


#' @references Goudet J. (1995).
#' FSTAT (Version 1.2): A computer program to calculate F-statistics.
#' Journal of Heredity 86:485-486

#' @seealso \href{https://github.com/jgx65/hierfstat}{hierfstat}

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


tidy_fstat <- function(data, strata = NULL, tidy = TRUE, filename = NULL) {

  # Checking for missing and/or default arguments-------------------------------
  if (missing(data)) rlang::abort("fstat file missing")


  # Import data ------------------------------------------------------------------
  if (is.vector(data)) {
    data <- readr::read_delim(file = data, delim = "?", col_names = "data", col_types = "c")
  } else {
    data %<>% dplyr::rename(data = V1)
  }

  # replace potential white space character: [\t\n\f\r\p{Z}] -----------------------------
  data$data %<>%
    stringi::stri_replace_all_regex(
      str = .,
      pattern = "\\s+",
      replacement = "\t",
      vectorize_all = FALSE
    )

  # metadata -------------------------------------------------------------------
  fstat.first.line <- data %>%
    dplyr::slice(1) %>%
    tidyr::separate(col = data, into = c("np", "nl", "nu", "allele.coding"), sep = "\t")

  # create a data frame --------------------------------------------------------

  # markers
  markers <- data %>%
    dplyr::slice(2:(as.numeric(fstat.first.line$nl) + 1)) %>%
    purrr::flatten_chr(.x = .)

  # Isolate the genotypes and pop column
  data %<>%
    dplyr::slice(-(1:(as.numeric(fstat.first.line$nl) + 1)))

  # separate the dataset by tab
  data <- tibble::as_tibble(
    do.call(rbind, stringi::stri_split_fixed(str = data$data, pattern = "\t")),
    .name_repair = "minimal"
  ) %>%
    magrittr::set_colnames(x = ., c("POP_ID", markers))

  # Create a string of id
  id <- tibble::tibble(INDIVIDUALS = paste0("IND-", seq_along(1:length(data$POP_ID))))

  # bind with data
  data <- dplyr::bind_cols(id, data)

  # Population info ------------------------------------------------------------
  # Strata
  if (!is.null(strata)) {
    if (is.vector(strata)) {
      number.columns.strata <- max(utils::count.fields(strata, sep = "\t"))
      col.types <- stringi::stri_join(rep("c", number.columns.strata), collapse = "")
      strata.df <- readr::read_tsv(file = strata, col_names = TRUE, col_types = col.types) %>%
        dplyr::rename(POP_ID = STRATA)
    } else {
      # message("strata object: yes")
      colnames(strata) <- stringi::stri_replace_all_fixed(
        str = colnames(strata),
        pattern = "STRATA",
        replacement = "POP_ID",
        vectorize_all = FALSE
      )
      strata.df <- strata
    }

    # Remove potential whitespace in pop_id
    strata.df$POP_ID <- stringi::stri_replace_all_fixed(
      strata.df$POP_ID,
      pattern = " ",
      replacement = "_",
      vectorize_all = FALSE
    )

    # Remove unwanted character in individual names
    strata.df <- strata.df %>%
      dplyr::mutate(
        INDIVIDUALS =  stringi::stri_replace_all_fixed(
          str = INDIVIDUALS,
          pattern = c("_", ",", ":"),
          replacement = c("-", "", "-"),
          vectorize_all = FALSE
        ),
        INDIVIDUALS = stringi::stri_replace_all_regex(
          str = INDIVIDUALS,
          pattern = "\\s+",
          replacement = "",
          vectorize_all = FALSE
        )
      )

    #join strata and data
    individuals <- dplyr::left_join(x = individuals, y = strata.df, by = "INDIVIDUALS")
  }

  # Tidy -------------------------------------------------------------------------

  # work on the genotype field
  data %<>%
    radiator::rad_long(
      x = .,
      cols = c("POP_ID", "INDIVIDUALS"),
      names_to = "MARKERS",
      values_to = "GENOTYPE",
      variable_factor = FALSE
      ) %>%
    tidyr::separate(
      col = GENOTYPE, into = c("A1", "A2"), sep = as.numeric(fstat.first.line$allele.coding)
    ) %>%
    dplyr::mutate(
      A1 = stringi::stri_pad_left(str = A1, pad = "0", width = 3),
      A2 = stringi::stri_pad_left(str = A2, pad = "0", width = 3)
    ) %>%
    tidyr::unite(GENOTYPE, A1, A2, sep = "")

  # wide format
  if (!tidy) {
    data %<>%
      radiator::rad_wide(
        x = .,
        formula = "POP_ID + INDIVIDUALS ~ MARKERS",
        values_from = "GENOTYPE"
      ) %>%
      dplyr::arrange(POP_ID, INDIVIDUALS)
  }

  # writing to a file  ---------------------------------------------------------
  if (!is.null(filename)) {
    readr::write_tsv(x = data, file = filename, col_names = TRUE)
  }

  return(data)
} # end tidy_fstat


# write_hierfstat --------------------------------------------------------------

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
#' @references Goudet, J. (1995) FSTAT (Version 1.2): A computer program to
#' calculate F- statistics. Journal of Heredity, 86, 485-486.
#' @references Goudet, J. (2005) hierfstat, a package for r to compute and test hierarchical F-statistics. Molecular Ecology Notes, 5, 184-186.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_hierfstat <- function(data, filename = NULL) {
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("Input file necessary to write the hierfstat file is missing")

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) data %<>% radiator::tidy_wide(data = ., import.metadata = TRUE)

  if (!rlang::has_name(data, "GT")) {
    data %<>% calibrate_alleles(data = ., verbose = FALSE) %$% input
  }

  data <- dplyr::select(.data = data, POP_ID, INDIVIDUALS, MARKERS, GT) %>%
    dplyr::arrange(MARKERS, POP_ID, INDIVIDUALS)

  # Create a marker vector  ------------------------------------------------
  markers <- dplyr::distinct(.data = data, MARKERS) %>%
    dplyr::arrange(MARKERS) %>%
    purrr::flatten_chr(.)

  # Get the number of sample (pop) for hierfstat -------------------------------
  if (is.factor(data$POP_ID)) data$POP_ID <- droplevels(data$POP_ID)

  np <- nlevels(droplevels(data$POP_ID))
  np.message <- stringi::stri_join("    * Number of sample pop, np = ", np, sep = "")
  message(np.message)

  # Get the number of loci -----------------------------------------------------
  nl <- length(markers)
  nl.message <- stringi::stri_join("    * Number of markers, nl = ", nl, sep = "")
  message(nl.message)

  data <- suppressWarnings(
    data %>%
      dplyr::select(MARKERS, POP_ID, INDIVIDUALS, GT) %>%
      dplyr::mutate(
        GT = replace(GT, which(GT == "000000"), NA),
        A1 = as.numeric(stringi::stri_sub(str = GT, from = 1, to = 3)),
        A2 = as.numeric(stringi::stri_sub(str = GT, from = 4, to = 6)),
        GT = NULL
      )
  )
  # Get the highest number used to label an allele -----------------------------
  nu <- max(c(unique(data$A1), unique(data$A2)), na.rm = TRUE)
  nu.message <- stringi::stri_join("    * The highest number used to label an allele, nu = ",
                                   nu, sep = "")
  message(nu.message)

  # prep the data  -------------------------------------------------------------
  data <- suppressWarnings(
    tidyr::unite(data = data, GT, A1, A2, sep = "") %>%
      dplyr::mutate(GT = as.numeric(GT)) %>%
      radiator::rad_wide(
        x = .,
        formula = "POP_ID + INDIVIDUALS ~ MARKERS",
        values_from = "GT"
        ) %<>%
      dplyr::arrange(POP_ID, INDIVIDUALS) %>%
      dplyr::mutate(POP_ID = as.integer(POP_ID), INDIVIDUALS = NULL)
  )

  # allele coding --------------------------------------------------------------
  allele.coding <- 1
  message("    * The alleles are encoded with one digit number")

  # Filename -------------------------------------------------------------------
  if (is.null(filename)) {
    filename <- stringi::stri_join("radiator_hierfstat_", file.date, ".dat")
  } else {
    filename <- stringi::stri_join(filename, "_hierfstat.dat")
  }


  # FSTAT: write the first line ------------------------------------------------
  fstat.first.line <- stringi::stri_join(np, nl, nu, allele.coding, sep = " ")
  fstat.first.line <- as.data.frame(fstat.first.line)
  readr::write_delim(x = fstat.first.line, file = filename, delim = "\n", append = FALSE,
                     col_names = FALSE)

  # FSTAT: write the locus name to the file
  loci.table <- as.data.frame(markers)
  readr::write_delim(x = loci.table, file = filename, delim = "\n", append = TRUE,
                     col_names = FALSE)

  # FSTAT: write the pop and genotypes
  readr::write_delim(x = data, na = "00", file = filename, delim = "\t", append = TRUE,
                     col_names = FALSE)
  data <- as.data.frame(data) # required by hierfstat...
  return(data)
}# End write_hierfstat
