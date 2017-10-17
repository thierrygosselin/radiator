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


#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join slice bind_cols bind_rows
#' @importFrom purrr invoke flatten_chr
#' @importFrom stringi stri_split_fixed stri_replace_all_fixed stri_join stri_pad_left stri_replace_all_regex
#' @importFrom readr write_tsv read_tsv read_delim
#' @importFrom utils count.fields
#' @importFrom tidyr separate gather
#' @importFrom tibble as_data_frame

#' @examples
#' \dontrun{
#' We will use the fstat dataset provided with adegenet package
#' if (!require("hierfstat")) install.packages("hierfstat")
#'
#' The simplest form of the function:
#' fstat.file <- radiator::tidy_fstat(
#' data = system.file(
#' "extdata/diploid.dat",
#' package = "hierfstat"
#' )
#' )
#'
#' # To output a data frame in wide format, with markers in separate columns:
#' nancycats.wide <- radiator::tidy_fstat(
#' data = system.file(
#' "extdata/diploid.dat",
#' package = "hierfstat"
#' ), tidy = FALSE
#' )
#' }


#' @references Goudet J. (1995).
#' FSTAT (Version 1.2): A computer program to calculate F-statistics.
#' Journal of Heredity 86:485-486

#' @seealso \href{https://github.com/jgx65/hierfstat}{hierfstat}

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


tidy_fstat <- function(data, strata = NULL, tidy = TRUE, filename = NULL) {

  # Checking for missing and/or default arguments-------------------------------
  if (missing(data)) stop("fstat file missing")


  # Import data ------------------------------------------------------------------
  if (is.vector(data)) {
    data <- readr::read_delim(file = data, delim = "?", col_names = "data", col_types = "c")
  } else {
    data <- dplyr::rename(data, data = V1)
  }

  # replace potential white space character: [\t\n\f\r\p{Z}] -----------------------------
  data$data = stringi::stri_replace_all_regex(
    str = data$data,
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
  data <- dplyr::slice(.data = data, -(1:(as.numeric(fstat.first.line$nl) + 1)))

  # separate the dataset by tab
  data <- tibble::as_data_frame(
    purrr::invoke(#similar to do.call
      rbind, stringi::stri_split_fixed(str = data$data, pattern = "\t")
    )
  )

  # new colnames
  colnames(data) <- c("POP_ID", markers)

  # Create a string of id
  id <- dplyr::data_frame(INDIVIDUALS = paste0("IND-", seq_along(1:length(data$POP_ID))))

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
  data <- tidyr::gather(
    data = data, key = MARKERS, value = GENOTYPE, -c(POP_ID, INDIVIDUALS)
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
    data <- data %>%
      dplyr::group_by(POP_ID, INDIVIDUALS) %>%
      tidyr::spread(data = ., key = MARKERS, value = GENOTYPE) %>%
      dplyr::arrange(POP_ID, INDIVIDUALS)
  }

  # writing to a file  ---------------------------------------------------------
  if (!is.null(filename)) {
    readr::write_tsv(x = data, path = filename, col_names = TRUE)
  }

  return(data)
} # end tidy_fstat
