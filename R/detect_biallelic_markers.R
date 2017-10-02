# Detect if markers are biallelic

#' @name detect_biallelic_markers

#' @title Detect biallelic data

#' @description Detect if markers in tidy dataset are biallelic.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.

#' @param verbose (optional, logical) \code{verbose = TRUE} to be chatty
#' during execution.
#' Default: \code{verbose = FALSE}.

#' @return A logical character string (TRUE/FALSE). That answer the question if
#' the data set is biallelic or not.

#' @export
#' @rdname detect_biallelic_markers
#' @importFrom dplyr select mutate group_by ungroup rename tally filter
#' @importFrom stringi stri_replace_all_fixed stri_sub
#' @importFrom tibble has_name
#' @importFrom tidyr gather
#' @importFrom purrr flatten_chr

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

detect_biallelic_markers <- function(data, verbose = FALSE) {

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

  # markers with all missing... yes I've seen it... breaks code...
  # input <- detect_all_missing(data = input)
  marker.problem <- radiator::detect_all_missing(data = input)
  if (marker.problem$marker.problem) {
    input <- marker.problem$data
  }
  marker.problem <- NULL

  # Detecting biallelic markers-------------------------------------------------
  if (verbose) message("Scanning for number of alleles per marker...")
  if (tibble::has_name(input, "ALT")) {
    alt.num <- max(unique(
      stringi::stri_count_fixed(str = unique(input$ALT), pattern = ","))) + 1

    if (alt.num > 1) {
      biallelic <- FALSE
      if (verbose) message("    Data is multi-allelic")
    } else {
      biallelic <- TRUE
      if (verbose) message("    Data is bi-allelic")
    }
    alt.num <- NULL
  } else {
    # 10% of the markers are randomly sampled
    sampled.markers <- unique(input$MARKERS)
    sampled.markers <- sample(x = sampled.markers,
                              size = length(sampled.markers) * 0.10,
                              replace = FALSE)

    # system.time(input2 <- dplyr::select(.data = input, MARKERS, GT) %>%
    #   dplyr::filter(GT != "000000") %>%
    #   dplyr::distinct(MARKERS, GT) %>%
    #   dplyr::mutate(A1 = stringi::stri_sub(GT, 1, 3), A2 = stringi::stri_sub(GT, 4,6)) %>%
    #   dplyr::select(-GT) %>%
    #   tidyr::gather(data = ., key = ALLELES_GROUP, value = ALLELES, -MARKERS) %>%
    #   dplyr::distinct(MARKERS, ALLELES) %>%
    #   dplyr::count(x = ., MARKERS) %>%
    #   dplyr::summarise(BIALLELIC = max(n, na.rm = TRUE)) %>%
    #   purrr::flatten_chr(.x = .))


    biallelic <- dplyr::select(.data = input, MARKERS, GT) %>%
      dplyr::filter(GT != "000000") %>%
      dplyr::filter(MARKERS %in% sampled.markers) %>%
      dplyr::distinct(MARKERS, GT) %>%
      dplyr::mutate(A1 = stringi::stri_sub(GT, 1, 3), A2 = stringi::stri_sub(GT, 4,6)) %>%
      dplyr::select(-GT) %>%
      tidyr::gather(data = ., key = ALLELES_GROUP, value = ALLELES, -MARKERS) %>%
      dplyr::distinct(MARKERS, ALLELES) %>%
      dplyr::count(x = ., MARKERS) %>%
      dplyr::select(n) #%>%
      # dplyr::summarise(BIALLELIC = max(n, na.rm = TRUE)) %>%
      # purrr::flatten_chr(.x = .) %>%
      # unique(.)

    biallelic <- max(biallelic$n)
    # if (length(biallelic) != 1) stop("Mix of bi- and multi-allelic markers is not supported")

    if (biallelic > 4) {
      biallelic <- FALSE
      if (verbose) message("    Data is multi-allelic")
    } else {
      biallelic <- TRUE
      if (verbose) message("    Data is bi-allelic")
    }
  }
  return(biallelic)
} # End detect_biallelic_markers
