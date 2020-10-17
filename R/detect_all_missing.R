# Detect markers with all missing genotypes

#' @name detect_all_missing

#' @title Detect markers with all missing genotypes

#' @description Detect if markers in tidy dataset have no genotypes at all.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment.

#' @return The filtered dataset if problematic markers were found. Otherwise,
#' the untouch dataset.

#' @export
# @keywords internal
#' @rdname detect_all_missing
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


detect_all_missing <- function(data) {

  # For GDS data...
  # test <- extract_genotypes_metadata(
  #   gds = gds,
  #   genotypes.meta.select = c("INDIVIDUALS", "VARIANT_ID", "GT_BIN")
  # )
  #
  # all.missing <- dplyr::group_by(test, VARIANT_ID) %>%
  #   dplyr::summarise(MISSING = all(is.na(GT_BIN))) %>%
  #   dplyr::filter(!MISSING)



  res <- list()

  # Allows to have either GT, GT_VCF_NUC, GT_VCF or GT_BIN
  # If more than 1 is discovered in data, keep 1 randomly.
  detect.gt <- purrr::keep(.x = colnames(data), .p = colnames(data) %in% c("GT", "GT_VCF_NUC", "GT_VCF","GT_BIN"))
  if (length(detect.gt) > 1) detect.gt <- sample(x = detect.gt, size = 1)
  want <- c("MARKERS", detect.gt)
  blacklist.markers <- suppressWarnings(dplyr::select(data, dplyr::one_of(want)))

  # markers with all missing... yes I've seen it... breaks code...
  if (tibble::has_name(blacklist.markers, "GT")) {
    blacklist.markers <- blacklist.markers %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::distinct(GT) %>%
      dplyr::summarise(GENOTYPED = length(GT[GT != "000000"])) %>%
      dplyr::filter(GENOTYPED == 0) %>%
      dplyr::ungroup(.)
  }

  if (tibble::has_name(blacklist.markers, "GT_BIN")) {
    blacklist.markers <- blacklist.markers %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::distinct(GT_BIN) %>%
      dplyr::summarise(GENOTYPED = length(GT_BIN[!is.na(GT_BIN)])) %>%
      dplyr::filter(GENOTYPED == 0) %>%
      dplyr::ungroup(.)
  }

  if (tibble::has_name(blacklist.markers, "GT_VCF_NUC")) {
    blacklist.markers <- blacklist.markers %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::distinct(GT_VCF_NUC) %>%
      dplyr::summarise(GENOTYPED = length(GT_VCF_NUC[GT_VCF_NUC != "./."])) %>%
      dplyr::filter(GENOTYPED == 0) %>%
      dplyr::ungroup(.)
  }

  if (tibble::has_name(blacklist.markers, "GT_VCF")) {
    blacklist.markers <- blacklist.markers %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::distinct(GT_VCF) %>%
      dplyr::summarise(GENOTYPED = length(GT_VCF[GT_VCF != "./."])) %>%
      dplyr::filter(GENOTYPED == 0) %>%
      dplyr::ungroup(.)
  }

  problem <- nrow(blacklist.markers)
  if (problem > 0) {
    message("Dataset contains ", problem," marker(s) with no genotypes (all missing)...")
    message("    removing problematic markers...")
    res$data <- dplyr::filter(data, !MARKERS %in% blacklist.markers$MARKERS)
    readr::write_tsv(
      x = blacklist.markers,
      file = "blacklist.markers.all.missing.tsv"
    )
    res$blacklist.markers.all.missing <- blacklist.markers
    res$marker.problem <- TRUE
  } else {
    res$marker.problem <- FALSE
  }


  return(res)
}#End detect_all_missing
