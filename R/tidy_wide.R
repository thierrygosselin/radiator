# read data frames in long/tidy and wide format

#' @name tidy_wide

#' @title Read/Import and tidy genomic data frames long or wide format

#' @description Read/Import and tidy genomic data frames long or wide format.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @param data A file in the working directory or object in the global environment
#' in wide or long (tidy) formats. See details for more info.
#'
#' \emph{How to get a tidy data frame ?}
#' \href{https://github.com/thierrygosselin/radiator}{radiator}
#' \code{\link{tidy_genomic_data}} can transform 6 genomic data formats
#' in a tidy data frame.

#' @param import.metadata (optional, logical) With \code{import.metadata = TRUE}
#' the metadata (anything else than the genotype) will be imported for the long
#' format exclusively. Default: \code{import.metadata = FALSE}, no metadata.

#' @param ... other parameters passed to the function.

#' @return A tidy data frame in the global environment.
#' @export
#' @rdname tidy_wide
#' @importFrom stringi stri_replace_all_fixed stri_pad_left
#' @importFrom dplyr mutate select
#' @importFrom tibble as_data_frame has_name
#' @importFrom tidyr gather spread

#' @details \strong{Input data:}
#'
#' To discriminate the long from the wide format,
#' the function \pkg{radiator} \code{\link[radiator]{tidy_wide}} searches
#' for \code{MARKERS or LOCUS} in column names (TRUE = long format).
#' The data frame is tab delimitted.

#' \strong{Wide format:}
#' The wide format cannot store metadata info.
#' The wide format starts with these 2 id columns:
#' \code{INDIVIDUALS}, \code{POP_ID} (that refers to any grouping of individuals),
#' the remaining columns are the markers in separate columns storing genotypes.
#'
#' \strong{Long/Tidy format:}
#' The long format is considered to be a tidy data frame and can store metadata info.
#' (e.g. from a VCF see \pkg{radiator} \code{\link{tidy_genomic_data}}). A minimum of 4 columns
#' are required in the long format: \code{INDIVIDUALS}, \code{POP_ID},
#' \code{MARKERS or LOCUS} and \code{GENOTYPE or GT}. The rest are considered metata info.
#'
#' \strong{2 genotypes formats are available:}
#' 6 characters no separator: e.g. \code{001002 or 111333} (for heterozygote individual).
#' 6 characters WITH separator: e.g. \code{001/002 or 111/333} (for heterozygote individual).
#' The separator can be any of these: \code{"/", ":", "_", "-", "."}.
#'
#' \emph{How to get a tidy data frame ?}
#' \pkg{radiator} \code{\link{tidy_genomic_data}} can transform 6 genomic data formats
#' in a tidy data frame.



#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

tidy_wide <- function(data, import.metadata = FALSE, ...) {

  # Checking for missing and/or default arguments ******************************
  if (missing(data)) stop("Input file argument is missing")

  if (is.vector(data)) {# for file in the working directory
    # Scan column names
    scan.colnames <- readr::read_tsv(
      file = data,
      col_types = readr::cols(.default = readr::col_character()),
      n_max = 1)

    # Determine long (tidy) or wide dataset
    if ("MARKERS" %in% colnames(scan.colnames) | "LOCUS" %in% colnames(scan.colnames) ) {
      long.format <- TRUE
    } else {
      long.format <- FALSE
    }


    if (long.format) { # long (tidy) format
      # to select columns while importing the file
      input <- readr::read_tsv(file = data, col_types = readr::cols(.default = readr::col_character()))


      # switch GENOTYPE for GT in colnames if found
      if ("GENOTYPE" %in% colnames(input)) {
        colnames(input) <- stringi::stri_replace_all_fixed(
          str = colnames(input),
          pattern = "GENOTYPE",
          replacement = "GT",
          vectorize_all = FALSE)
      }

      if (!import.metadata) {
        if ("MARKERS" %in% colnames(input) & "LOCUS" %in% colnames(input)) {
          input <- dplyr::select(.data = input, POP_ID, INDIVIDUALS, LOCUS = MARKERS, GT)
        } else if ("MARKERS" %in% colnames(input) & !"LOCUS" %in% colnames(input)) {
          input <- dplyr::select(.data = input, POP_ID, INDIVIDUALS, LOCUS = MARKERS, GT)
        } else {
          input <- dplyr::select(.data = input, POP_ID, INDIVIDUALS, LOCUS, GT)
        }
      }

    } else {# wide format
      input <- readr::read_tsv(
        file = data,
        col_types = readr::cols(.default = readr::col_character())) %>%
        tidyr::gather(data = ., key = LOCUS, value = GT, -c(POP_ID, INDIVIDUALS))
    }
  } else {# object in global environment
    input <- data

    # Determine long (tidy) or wide dataset
    if ("MARKERS" %in% colnames(input) | "LOCUS" %in% colnames(input) ) {
      long.format <- TRUE
    } else {
      long.format <- FALSE
    }


    if (long.format) { # long (tidy) format
      # switch GENOTYPE for GT in colnames if found
      if ("GENOTYPE" %in% colnames(input)) {
        colnames(input) <- stringi::stri_replace_all_fixed(
          str = colnames(input),
          pattern = "GENOTYPE",
          replacement = "GT",
          vectorize_all = FALSE
        )
      }

      if (!import.metadata) {
        if ("MARKERS" %in% colnames(input) & "LOCUS" %in% colnames(input)) {
          input <- dplyr::select(.data = input, POP_ID, INDIVIDUALS, LOCUS = MARKERS, GT)
        } else if ("MARKERS" %in% colnames(input) & !"LOCUS" %in% colnames(input)) {
          input <- dplyr::select(.data = input, POP_ID, INDIVIDUALS, LOCUS = MARKERS, GT)
        } else {
          input <- dplyr::select(.data = input, POP_ID, INDIVIDUALS, LOCUS, GT)
        }
      }

    } else {# wide format
      input <- tidyr::gather(data = input, key = LOCUS, value = GT, -c(POP_ID, INDIVIDUALS))
    }
  }

  # unused objects
  scan.colnames <- NULL

  # Remove unwanted sep in the genotype filed.
  input <- input %>%
    dplyr::mutate(
      GT = stringi::stri_replace_all_fixed(
        str = as.character(GT),
        pattern = c("/", ":", "_", "-", "."),
        replacement = "",
        vectorize_all = FALSE),
      GT = stringi::stri_pad_left(str = as.character(GT), pad = "0", width = 6)
    )
  return(input)
}
