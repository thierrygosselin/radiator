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
#' of markers between populations.
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
#' @importFrom UpSetR upset

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

keep_common_markers <- function(data, plot = FALSE, verbose = FALSE) {

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file missing")

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    input <- radiator::tidy_wide(data = data, import.metadata = TRUE)
  } else {
    input <- data
  }

  # check genotype column naming
  if (tibble::has_name(input, "GENOTYPE")) {
    colnames(input) <- stringi::stri_replace_all_fixed(
      str = colnames(input),
      pattern = "GENOTYPE",
      replacement = "GT",
      vectorize_all = FALSE)
  }

  # necessary steps to make sure we work with unique markers and not duplicated LOCUS
  if (tibble::has_name(input, "LOCUS") && !tibble::has_name(input, "MARKERS")) {
    input <- dplyr::rename(.data = input, MARKERS = LOCUS)
  }

  # markers.meta
  want <- c("MARKERS", "CHROM", "LOCUS", "POS")
  markers.meta <- suppressWarnings(dplyr::select(input, dplyr::one_of(want)) %>%
    dplyr::distinct(MARKERS, .keep_all = TRUE))

  if (verbose) message("Using markers common in all populations:")
  blacklist <- dplyr::select(.data = input, MARKERS, POP_ID, GT) %>%
    dplyr::filter(GT != "000000") %>%
    dplyr::distinct(MARKERS, POP_ID) %>%
    dplyr::count(x = ., MARKERS) %>%
    dplyr::filter(n != dplyr::n_distinct(input$POP_ID)) %>%
    dplyr::distinct(MARKERS) %>%
    dplyr::arrange(MARKERS)

  markers.input <- dplyr::n_distinct(input$MARKERS)
  blacklist.markers <- nrow(blacklist)
  markers.in.common <- markers.input - blacklist.markers

  if (verbose) message("    Number of markers before = ", markers.input)
  if (verbose) message("    Number of markers removed = ", blacklist.markers)
  if (verbose) message("    Number of markers after (common between populations) = ", markers.in.common)

  if (plot) {
    pops <- unique(input$POP_ID)

    if (length(pops) > 1) {

      plot.data <- dplyr::filter(input, GT != "000000") %>%
        dplyr::distinct(MARKERS, POP_ID) %>%
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
    input <- dplyr::anti_join(input, blacklist, by = "MARKERS")

    if (ncol(markers.meta) > 1) {
      blacklist <- dplyr::left_join(blacklist, markers.meta, by = "MARKERS")
    }

  } else {
    blacklist <- tibble::data_frame(INDIVIDUALS = character(0))
  }
  want <- c("MARKERS", "CHROM", "LOCUS", "POS")
  whitelist <- suppressWarnings(dplyr::select(input, dplyr::one_of(want)) %>%
    dplyr::distinct(MARKERS, .keep_all = TRUE))

  return(res = list(input = input,
                    whitelist.common.markers = whitelist,
                    blacklist.not.in.common.markers = blacklist))
}
