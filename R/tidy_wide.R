# read data frames in long/tidy and wide format

#' @name tidy_wide

#' @title Read/Import a tidy genomic data frames.

#' @description Read/Import and tidy genomic data frames. If data is in
#' wide format, the functions will gather the data.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @param data A file in the working directory or object in the global environment
#' in wide or long (tidy) formats. See details for more info.
#'
#' \emph{How to get a tidy data frame ?}
#' \href{https://github.com/thierrygosselin/radiator}{radiator}
#' \code{\link{tidy_genomic_data}}.

#' @param import.metadata (optional, logical) With \code{import.metadata = TRUE}
#' the metadata (anything else than the genotype) will be imported for the long
#' format exclusively. Default: \code{import.metadata = FALSE}, no metadata.

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
#' for \code{MARKERS} in column names (TRUE = long format).
#' The data frame is tab delimitted.

#' \strong{Wide format:}
#' The wide format cannot store metadata info.
#' The wide format starts with these 2 id columns:
#' \code{INDIVIDUALS}, \code{POP_ID} (that refers to any grouping of individuals),
#' the remaining columns are the markers in separate columns storing genotypes.
#'
#' \strong{Long/Tidy format:}
#' The long format is considered to be a tidy data frame and can store metadata info.
#' (e.g. from a VCF see \pkg{radiator} \code{\link{tidy_genomic_data}}).
#' A minimum of 4 columns
#' are required in the long format: \code{INDIVIDUALS}, \code{POP_ID},
#' \code{MARKERS} and \code{GT} for the genotypes.
#' The remaining columns are considered metadata info.
#'
#' \strong{Genotypes with separators:}
#' ALL separators will be removed.
#' Genotypes should be coded with 3 integers for each alleles.
#' 6 integers in total for the genotypes.
#' e.g. \code{001002 or 111333} (for heterozygote individual).
#' 6 integers WITH separator: e.g. \code{001/002 or 111/333} (for heterozygote individual).
#' The separator can be any of these: \code{"/", ":", "_", "-", "."}, and will
#' be removed.
#'
#'
#' \strong{separators in POP_ID, INDIVIDUALS and MARKERS:}
#' Some separators can interfere with packages or codes and are cleaned by radiator.
#' \itemize{
#' \item MARKERS: \code{/}, \code{:}, \code{-} and \code{.} are changed to an
#' underscore
#' \code{_}.
#' \item POP_ID: white spaces in population names are replaced by underscore.
#' \item INDIVIDUALS: \code{_} and \code{:} are changed to a dash \code{-}
#' }
#'
#' \emph{How to get a tidy data frame ?}
#' \pkg{radiator} \code{\link{tidy_genomic_data}} can transform 6 genomic data formats
#' in a tidy data frame.


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

tidy_wide <- function(data, import.metadata = FALSE) {

  # Checking for missing and/or default arguments-------------------------------
  if (missing(data)) stop("Input file argument is missing")

  if (is.vector(data)) {# for file in the working directory
    if (stringi::stri_detect_fixed(
      str = stringi::stri_sub(str = data, from = -4, to = -1),
      pattern = ".tsv")) {
      data <- readr::read_tsv(file = data, col_types = readr::cols(.default = readr::col_character()))
    } else if (radiator::detect_genomic_format(data) == "fst.file") {
      data <- radiator::read_rad(data = data)
    }
  }

  # Determine long (tidy) or wide dataset
  if (!"MARKERS" %in% colnames(data) && !"LOCUS" %in% colnames(data)) {
    data <- tidyr::gather(data = data, key = MARKERS, value = GT, -c(POP_ID, INDIVIDUALS))
  }

  # necessary steps to make sure we work with unique markers and not duplicated LOCUS
  if (tibble::has_name(data, "LOCUS") && !tibble::has_name(data, "MARKERS")) {
    data <- dplyr::rename(.data = data, MARKERS = LOCUS)
  }

  # reproducibility for old format
  if (tibble::has_name(data, "GENOTYPE")) {
    colnames(data) <- stringi::stri_replace_all_fixed(
      str = colnames(data),
      pattern = "GENOTYPE",
      replacement = "GT",
      vectorize_all = FALSE
    )
  }
  if (!import.metadata) {
    want <- c("POP_ID", "INDIVIDUALS", "MARKERS", "CHROM", "LOCUS", "POS", "GT",
              "GT_VCF_NUC", "GT_VCF", "GT_BIN")
    data <- suppressWarnings(dplyr::select(data, dplyr::one_of(want)))
  }

  # Remove unwanted sep in the genotypes (if found)
  if (tibble::has_name(data, "GT")) {
    gt.sep <- unique(
      stringi::stri_detect_fixed(
        str = sample(x = data$GT, size = 5, replace = FALSE),
        pattern = c("/", ":", "_", "-", ".")))
    if (length(gt.sep) > 1) gt.sep <- TRUE
    if (gt.sep) {
      data <- data %>%
        dplyr::mutate(
          GT = stringi::stri_replace_all_fixed(
            str = as.character(GT),
            pattern = c("/", ":", "_", "-", "."),
            replacement = "",
            vectorize_all = FALSE),
          GT = stringi::stri_pad_left(str = as.character(GT), pad = "0", width = 6))
    }
  }

  # clean markers names
  if (tibble::has_name(data, "MARKERS")) {
    data$MARKERS <- clean_markers_names(data$MARKERS)
  }

  # clean id names
  data$INDIVIDUALS <- clean_ind_names(data$INDIVIDUALS)

  # clean pop id
  data$POP_ID <- clean_pop_names(data$POP_ID)

  # Make sure no data groupings exists
  if (!is.null(dplyr::groups(data))) data <- dplyr::ungroup(data)
  return(data)
}#End tidy_wide
