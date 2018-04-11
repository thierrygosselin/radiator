# Discard monomorphic markers

#' @name discard_monomorphic_markers

#' @title Discard monomorphic markers

#' @description Discard monomorphic markers.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.


#' @param verbose (optional, logical) \code{verbose = TRUE} to be chatty
#' during execution.
#' Default: \code{verbose = FALSE}.

#' @return A list with the filtered input file and the blacklist of markers removed.

#' @export
#' @rdname discard_monomorphic_markers
#' @importFrom dplyr select mutate group_by ungroup rename tally filter semi_join n_distinct
#' @importFrom stringi stri_replace_all_fixed stri_join
#' @importFrom tibble has_name

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

discard_monomorphic_markers <- function(data, verbose = FALSE) {

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file missing")

  # Import data ---------------------------------------------------------------
  data <- radiator::tidy_wide(data = data, import.metadata = TRUE)


  if (tibble::has_name(data, "CHROM")) {
    markers.df <- dplyr::distinct(.data = data, MARKERS, CHROM, LOCUS, POS)
  }
  if (verbose) message("Scanning for monomorphic markers...")
  if (verbose) message("    Number of markers before = ", dplyr::n_distinct(data$MARKERS))

  if (tibble::has_name(data, "GT_BIN")) {
    mono.markers <- dplyr::select(.data = data, MARKERS, GT_BIN) %>%
      dplyr::filter(!is.na(GT_BIN)) %>%
      dplyr::distinct(MARKERS, GT_BIN) %>%
      dplyr::count(x = ., MARKERS) %>%
      dplyr::filter(n == 1) %>%
      dplyr::distinct(MARKERS)
  } else {
    mono.markers <- dplyr::select(.data = data, MARKERS, GT) %>%
      dplyr::filter(GT != "000000") %>%
      dplyr::distinct(MARKERS, GT) %>%
      dplyr::mutate(
        A1 = stringi::stri_sub(GT, 1, 3),
        A2 = stringi::stri_sub(GT, 4,6)
      ) %>%
      dplyr::select(-GT) %>%
      tidyr::gather(data = ., key = ALLELES_GROUP, value = ALLELES, -MARKERS) %>%
      dplyr::distinct(MARKERS, ALLELES) %>%
      dplyr::count(x = ., MARKERS) %>%
      dplyr::filter(n == 1) %>%
      dplyr::distinct(MARKERS)
  }
  # Remove the markers from the dataset
  if (verbose) message("    Number of monomorphic markers removed = ", nrow(mono.markers))

  if (length(mono.markers$MARKERS) > 0) {
    data <- dplyr::anti_join(data, mono.markers, by = "MARKERS")
    if (verbose) message("    Number of markers after = ", dplyr::n_distinct(data$MARKERS))
    if (tibble::has_name(data, "CHROM")) {
      mono.markers <- dplyr::left_join(mono.markers, markers.df, by = "MARKERS")
    }
  } else {
    mono.markers <- tibble::data_frame(MARKERS = character(0))
  }

  want <- c("MARKERS", "CHROM", "LOCUS", "POS")
  whitelist.polymorphic.markers <- suppressWarnings(
    dplyr::select(data, dplyr::one_of(want)) %>%
      dplyr::distinct(MARKERS, .keep_all = TRUE))
  res <- list(input = data,
              blacklist.monomorphic.markers = mono.markers,
              whitelist.polymorphic.markers = whitelist.polymorphic.markers
  )
  return(res)
} # end discard mono markers

