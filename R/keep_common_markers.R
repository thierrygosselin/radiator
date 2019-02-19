# Keep markers in common between all populations

#' @name keep_common_markers

#' @title Keep markers in common between all populations

#' @description Keep markers in common between all populations.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.


#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.

#' @param plot (optional, logical) \code{plot = TRUE} will produce an
#' \href{https://github.com/hms-dbmi/UpSetR}{UpSet plot} to visualize the number
#' of markers between populations. The package is required for this to work...
#' Default: \code{plot = FALSE}.

#' @param verbose (optional, logical) \code{verbose = TRUE} to be chatty
#' during execution.
#' Default: \code{verbose = FALSE}.

#' @return A list with the filtered input data set,
#' the whitelist of markers in common and the blacklist of markers discarded.


#' @export
#' @rdname keep_common_markers
#' @importFrom dplyr select mutate group_by ungroup rename tally filter semi_join n_distinct
#' @importFrom stringi stri_replace_all_fixed stri_join
#' @importFrom tibble has_name
# @importFrom UpSetR upset

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @section Life cycle:
#'
#' As of radiator v.0.0.22, keep.common.markers function is deprecated and was replaced by
#' \code{\link{filter_common_markers}} in an effort to have more meaningful functions names.

keep_common_markers <- function(data, plot = FALSE, verbose = FALSE) {

  message("\n\nDeprecated function, update your code to use: filter_common_markers\n\n")


  if (plot) {
    if (!requireNamespace("UpSetR", quietly = TRUE)) {
      rlang::abort("UpSetR needed for this function to work
         Install with install.packages('UpSetR')")
    }
  }

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("Input file missing")

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    data <- radiator::tidy_wide(data = data, import.metadata = TRUE)
  }

  # check genotype column naming
  if (tibble::has_name(data, "GENOTYPE")) {
    colnames(data) <- stringi::stri_replace_all_fixed(
      str = colnames(data),
      pattern = "GENOTYPE",
      replacement = "GT",
      vectorize_all = FALSE)
  }

  # necessary steps to make sure we work with unique markers and not duplicated LOCUS
  if (tibble::has_name(data, "LOCUS") && !tibble::has_name(data, "MARKERS")) {
    data <- dplyr::rename(.data = data, MARKERS = LOCUS)
  }

  # markers.meta
  want <- c("MARKERS", "CHROM", "LOCUS", "POS")
  markers.meta <- suppressWarnings(dplyr::select(data, dplyr::one_of(want)) %>%
                                     dplyr::distinct(MARKERS, .keep_all = TRUE))

  if (verbose) message("Using markers common in all populations:")

  if (tibble::has_name(data, "GT_BIN")) {
    blacklist <- dplyr::select(.data = data, MARKERS, POP_ID, GT_BIN) %>%
      dplyr::filter(!is.na(GT_BIN))
  } else {
    blacklist <- dplyr::select(.data = data, MARKERS, POP_ID, GT) %>%
      dplyr::filter(GT != "000000")
  }

  blacklist <- dplyr::distinct(blacklist, MARKERS, POP_ID) %>%
    dplyr::count(x = ., MARKERS) %>%
    dplyr::filter(n != dplyr::n_distinct(data$POP_ID)) %>%
    dplyr::distinct(MARKERS) %>%
    dplyr::arrange(MARKERS)

  markers.data <- dplyr::n_distinct(data$MARKERS)
  blacklist.markers <- nrow(blacklist)
  markers.in.common <- markers.data - blacklist.markers

  if (verbose) {
    n.markers <- stringi::stri_join(markers.data, blacklist.markers, markers.in.common, sep = "/")
    message("    Number of markers before/blacklisted/after:", n.markers)
  }

  if (plot) {
    pops <- unique(data$POP_ID)

    if (length(pops) > 1) {
      if (tibble::has_name(data, "GT_BIN")) {
        plot.data <- dplyr::filter(data, !is.na(GT_BIN))
      } else {
        plot.data <- dplyr::filter(data, GT != "000000")
      }
      plot.data <- dplyr::distinct(plot.data, MARKERS, POP_ID) %>%
        dplyr::mutate(
          n = rep(1, n()),
          POP_ID = stringi::stri_join("POP_", POP_ID)
        ) %>%
        tidyr::spread(data = ., key = POP_ID, value = n, fill = 0) %>%
        data.frame(.)

      UpSetR::upset(plot.data, nsets = length(pops), order.by = "freq", empty.intersections = "on")
      message("    UpSet plot generated to visualize common markers")
    }
  }

  if (blacklist.markers > 0) {
    # system.time(data2 <- dplyr::anti_join(data, blacklist, by = "MARKERS"))
    data <- dplyr::filter(data, !MARKERS %in% blacklist$MARKERS)

    if (ncol(markers.meta) > 1) {
      blacklist <- dplyr::left_join(blacklist, markers.meta, by = "MARKERS") %>%
        readr::write_tsv(x = ., path = "blacklist.not.in.common.markers.tsv")
    }

  } else {
    blacklist <- tibble::data_frame(INDIVIDUALS = character(0))
  }
  want <- c("MARKERS", "CHROM", "LOCUS", "POS")
  whitelist <- suppressWarnings(dplyr::select(data, dplyr::one_of(want)) %>%
                                dplyr::distinct(MARKERS, .keep_all = TRUE))


  if (blacklist.markers > 0) {
    readr::write_tsv(x = whitelist, path = "whitelist.common.markers.tsv")
  }

  return(res = list(input = data,
                    whitelist.common.markers = whitelist,
                    blacklist.not.in.common.markers = blacklist))
}
