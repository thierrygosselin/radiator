# Detect allele problems

#' @name detect_allele_problems

#' @title Detect alternate allele problems

#' @description Detect alternate allele problems.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users. The function computes alternate allele
#' counts and looks for alternate alleles bellow a certain threshold.
#' The summary statistics for the
#' markers with problematic allele is computed based on coverage and
#' genotype likelihood. This function is very fast to highlight: i) bias in
#' representation of allelic copies and unequal coverage and ii) type I error
#' during genotyping of heterozygote with low coverage data.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.

#' @inheritParams tidy_genomic_data

#' @param allele.threshold (integer) Threshold of alternate allele copies.
#' Below this threshold markers are blacklisted. Choose this threshold based
#' on your tolerance for uninformative or unreliable information and
#' minor allele frequency. Summary statistics for this blacklist will be calculated.
#' See details for more info.

#' @param ... (optional) Advance mode that allows to pass further arguments
#' for fine-tuning the function (see details).
#'
#'
#' @return A list with summary, plot and blacklist.
#' Summary statistics for the blacklisted markers are generated and include
#' coverage and genotype likelihood information.
#'
#' \code{$alt.allele.count.plot}: Distribution of marker's alternate allele number.
#' The shape is usually skewed right.
#'
#' \code{$alt.depth.count.distribution.plot}: Distribution of markers alternate
#' allele depth (in number of read). The range showed on the plot is
#' between 1 and 10, and you decide your tolerance for low coverage
#' and good genotype calls.
#'
#' This plot is associated with \code{$number.markers}, a table that shows the
#' number of markers <= \code{allele.threshold }. Also in the table, the
#' range of alternate allele read depth (1 to 10), same as plot but showing here
#' the the cumulative number of markers. Again, would you trust a marker with only
#' 1 alternate allele (an heterozygote individual) with a read depth of 1 or 2 ?


#' @details
#' \strong{Under developement, use with caution}
#'
#'
#' \strong{Input files:} see \pkg{radiator} \code{\link[radiator]{tidy_genomic_data}}
#' for detailed information about supported file format.
#'
#'
#' \strong{allele.threshold}:
#'
#' An \code{allele.threshold = 1}, says that you want to check the statistics
#' for markers with only 1 copy of the alternate allele. Said differently,
#' \code{allele.threshold = 1} = only 1 heterozygote individual is making this
#' marker polymorphic, very thin line if this allele is backed by less than 3 reads.
#'
#' Similarly, a \code{allele.threshold = 2} can represent markers with
#' 2 heterozygote individuals calling the polymorphic markers or 1 individual
#' homozygote for the alternate allele.
#'
#' Look for problem with :
#' \itemize{
#' \item read depth/coverage problem
#' \item high imbalance between the alt and ref depth coverage
#' \item genotype likelihood abnormally lower than normal
#' \item alternate allele with NA for read depth information
#' }
#'
#' \strong{Alternate allele with no depth information:}
#' You can highlight those problematic genotypes with
#' alt.no.depth <- dplyr::filter(your.object$allele.summary, is.na(ALLELE_ALT_DEPTH))
#'
#'
#' Use the summary statistics of the blacklisted markers to refine update
#' the blacklist of markers. You can decide to discard the markers or
#' blacklist the problematic genotypes.
#' \pkg{radiator} \code{\link[radiator]{tidy_genomic_data}} and
#' \code{\link[radiator]{genomic_converter}} allows to erase problematic genotypes
#' using a blacklist of genotypes...



#' @examples
#' \dontrun{
#' problem <- detect_allele_problems(
#'         data = salamander,
#'         strata = "strata.salmon.tsv",
#'         allele.threshold = 3)
#'
#' # The default with this function:
#' # filter.monomorphic = TRUE # = discarded
#' # filter.common.markers = TRUE # markers not in common between pop are discarded
#' }

#' @export
#' @rdname detect_allele_problems
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

detect_allele_problems <- function(
  data,
  allele.threshold = 3,
  verbose = TRUE,
  ...
) {

  # test
  data <- "/Users/thierry/Dropbox/r_packages/package_testing/carol/p1.mac4.r0.8/filtered_carol/corals_20181120@1828_filtered.rad"
  allele.threshold = 3
  verbose = FALSE

  cat("#######################################################################\n")
  cat("################## radiator::detect_allele_problems ###################\n")
  cat("#######################################################################\n")
  timing <- proc.time()

  # Checking for missing and/or default arguments-------------------------------
  if (missing(data)) rlang::abort("Input file missing")


  # Import data ---------------------------------------------------------------
  if (is.vector(data)) data <- radiator::tidy_wide(data, import.metadata = TRUE)

  # Check that coverage info is there
  if (!tibble::has_name(data, "READ_DEPTH")) rlang::abort("This function requires read depth metadata")
  if (!tibble::has_name(data, "ALLELE_ALT_DEPTH")) rlang::abort("This function requires allele depth metadata")
  if (!(tibble::has_name(data, "GL") | tibble::has_name(data, "PL"))) rlang::abort("This function requires genotype likelihood metadata")


  message("\nComputing allele count...")

  # Allele count ---------------------------------------------------------------
  # For each marker, count the number of alternate allele
  allele.count <- data %>%
    dplyr::filter(GT_BIN %in% c(1,2)) %>%
    dplyr::group_by(MARKERS) %>%
    dplyr::summarise(n = sum(GT_BIN, na.rm = TRUE))

  # test
  # test <- dplyr::filter(.data = allele.count, n == 1) %>%
  # dplyr::select(MARKERS) %>%
  # dplyr::inner_join(data, by = "MARKERS")

  # Histogram of allele counts--------------------------------------------------
  alt.allele.count.plot <- ggplot2::ggplot(allele.count, ggplot2::aes(n)) +
    ggplot2::geom_bar() +
    ggplot2::labs(x = "Alternate allele (copy)") +
    ggplot2::labs(y = "Distribution (marker number)") +
    # scale_x_continuous(breaks = c(0, 0.05, 0.1, 0.2, 0.5, 1),
    # labels = c("0", "0.05", "0.1", "0.2", "0.5", "1.0")) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 12, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 12, face = "bold"),
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.text = ggplot2::element_text(size = 12, face = "bold"),
      strip.text.x = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 8)#, angle = 90, hjust = 1, vjust = 0.5)
    )
  # alt.allele.count.plot

  # Summary stats and Blacklist ------------------------------------------------
  blacklist <- dplyr::filter(.data = allele.count, n <= allele.threshold) %>%
    dplyr::select(ALLELE_COPIES = n, MARKERS)

  if (tibble::has_name(data, "GL")) {
    data.info.needed <- data %>%
      dplyr::select(MARKERS, CHROM, LOCUS, POS, REF, ALT, READ_DEPTH, ALLELE_REF_DEPTH, ALLELE_ALT_DEPTH, GL) %>%
      dplyr::group_by(MARKERS, CHROM, LOCUS, POS, REF, ALT) %>%
      dplyr::summarise(
        READ_DEPTH_MEAN = mean(READ_DEPTH, na.rm = TRUE),
        READ_DEPTH_MIN = min(READ_DEPTH, na.rm = TRUE),
        READ_DEPTH_MAX = max(READ_DEPTH, na.rm = TRUE),
        ALLELE_REF_DEPTH_MEAN = mean(ALLELE_REF_DEPTH, na.rm = TRUE),
        ALLELE_REF_DEPTH_MIN = min(ALLELE_REF_DEPTH, na.rm = TRUE),
        ALLELE_REF_DEPTH_MAX = max(ALLELE_REF_DEPTH, na.rm = TRUE),
        ALLELE_ALT_DEPTH_MEAN = mean(ALLELE_ALT_DEPTH, na.rm = TRUE),
        ALLELE_ALT_DEPTH_MIN = min(ALLELE_ALT_DEPTH, na.rm = TRUE),
        ALLELE_ALT_DEPTH_MAX = max(ALLELE_ALT_DEPTH, na.rm = TRUE),
        GL_MEAN = mean(GL, na.rm = TRUE),
        GL_MIN = min(GL, na.rm = TRUE),
        GL_MAX = max(GL, na.rm = TRUE)
      )
  }
  if (tibble::has_name(data, "PL")) {
    data.info.needed <- data %>%
      dplyr::select(MARKERS, CHROM, LOCUS, POS, REF, ALT, READ_DEPTH, ALLELE_REF_DEPTH, ALLELE_ALT_DEPTH, PROB_HOM_REF, PROB_HET, PROB_HOM_ALT) %>%
      dplyr::group_by(MARKERS, CHROM, LOCUS, POS, REF, ALT) %>%
      dplyr::summarise(
        READ_DEPTH_MEAN = mean(READ_DEPTH, na.rm = TRUE),
        READ_DEPTH_MIN = min(READ_DEPTH, na.rm = TRUE),
        READ_DEPTH_MAX = max(READ_DEPTH, na.rm = TRUE),
        ALLELE_REF_DEPTH_MEAN = mean(ALLELE_REF_DEPTH, na.rm = TRUE),
        ALLELE_REF_DEPTH_MIN = min(ALLELE_REF_DEPTH, na.rm = TRUE),
        ALLELE_REF_DEPTH_MAX = max(ALLELE_REF_DEPTH, na.rm = TRUE),
        ALLELE_ALT_DEPTH_MEAN = mean(ALLELE_ALT_DEPTH, na.rm = TRUE),
        ALLELE_ALT_DEPTH_MIN = min(ALLELE_ALT_DEPTH, na.rm = TRUE),
        ALLELE_ALT_DEPTH_MAX = max(ALLELE_ALT_DEPTH, na.rm = TRUE),
        PROB_HOM_REF_MEAN = mean(PROB_HOM_REF, na.rm = TRUE),
        PROB_HET_MEAN = mean(PROB_HET, na.rm = TRUE),
        PROB_HOM_ALT_MEAN = mean(PROB_HOM_ALT, na.rm = TRUE)
      )
  }

  blacklist.summary <- dplyr::left_join(blacklist, data.info.needed, by = "MARKERS") %>%
    dplyr::arrange(ALLELE_COPIES, CHROM, LOCUS, POS)
  data.info.needed <- NULL

  allele.summary <- dplyr::left_join(blacklist, data, by = "MARKERS") %>%
    dplyr::filter(GT_BIN %in% c(1, 2)) %>%
    dplyr::select(-GT) %>%
    dplyr::left_join(
      dplyr::select(.data = blacklist.summary, -c(ALLELE_COPIES, CHROM, LOCUS, POS, REF, ALT))
      , by = "MARKERS"
    ) %>%
    dplyr::arrange(ALLELE_COPIES)

  problem <- dplyr::full_join(
    allele.count %>%
      dplyr::rename(ALT_ALLELE_NUMBER = n) %>%
      dplyr::group_by(ALT_ALLELE_NUMBER) %>%
      dplyr::tally(.) %>% # count alleles, longest part of the block
      dplyr::ungroup(.) %>%
      dplyr::filter(ALT_ALLELE_NUMBER <= allele.threshold),
    dplyr::select(allele.summary, ALLELE_COPIES, MARKERS, INDIVIDUALS, ALLELE_ALT_DEPTH) %>%
      dplyr::filter(!is.na(ALLELE_ALT_DEPTH)) %>%
      dplyr::group_by(ALLELE_COPIES) %>%
      dplyr::summarise(
        "ALT_DEPTH = 1" = length(ALLELE_ALT_DEPTH[ALLELE_ALT_DEPTH == 1]),
        "ALT_DEPTH <= 2" = length(ALLELE_ALT_DEPTH[ALLELE_ALT_DEPTH <= 2]),
        "ALT_DEPTH <= 3" = length(ALLELE_ALT_DEPTH[ALLELE_ALT_DEPTH <= 3]),
        "ALT_DEPTH <= 4" = length(ALLELE_ALT_DEPTH[ALLELE_ALT_DEPTH <= 4]),
        "ALT_DEPTH <= 5" = length(ALLELE_ALT_DEPTH[ALLELE_ALT_DEPTH <= 5]),
        "ALT_DEPTH <= 6" = length(ALLELE_ALT_DEPTH[ALLELE_ALT_DEPTH <= 6]),
        "ALT_DEPTH <= 7" = length(ALLELE_ALT_DEPTH[ALLELE_ALT_DEPTH <= 7]),
        "ALT_DEPTH <= 8" = length(ALLELE_ALT_DEPTH[ALLELE_ALT_DEPTH <= 8]),
        "ALT_DEPTH <= 9" = length(ALLELE_ALT_DEPTH[ALLELE_ALT_DEPTH <= 9]),
        "ALT_DEPTH <= 10" = length(ALLELE_ALT_DEPTH[ALLELE_ALT_DEPTH <= 10])
      ) %>%
      dplyr::rename(ALT_ALLELE_NUMBER = ALLELE_COPIES)
    , by = "ALT_ALLELE_NUMBER"
  )

  alt.depth.count.distribution.plot <- dplyr::filter(allele.summary, ALLELE_ALT_DEPTH <= 10) %>%
    ggplot2::ggplot(data = ., ggplot2::aes(ALLELE_ALT_DEPTH)) +
    ggplot2::geom_bar(stat = "count", na.rm = TRUE) + #breaks = c(1:10), binwidth = 0.5, bins = 10) +
    ggplot2::labs(x = "Alternate allele depth (read number)") +
    ggplot2::labs(y = "Distribution (marker number)") +
    ggplot2::scale_x_continuous(name = ggplot2::waiver(), breaks = 1:10, labels = 1:10) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 12, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 12, face = "bold"),
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.text = ggplot2::element_text(size = 12, face = "bold"),
      strip.text.x = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 8)#, angle = 90, hjust = 1, vjust = 0.5)
    ) +
    ggplot2::facet_grid(~ ALLELE_COPIES)
  # alt.depth.count.distribution.plot

  readr::write_tsv(x = allele.summary, path = "problematic.allele.summary.stats.tsv")
  message("Written to the working directory:\nproblematic.allele.summary.stats.tsv")

  blacklist.markers <- dplyr::ungroup(allele.summary) %>%
    dplyr::distinct(MARKERS, CHROM, LOCUS, POS)

  readr::write_tsv(x = blacklist.markers, path = "blacklist.markers.allele.problem.tsv")
  message("Written to the working directory:\nblacklist.markers.allele.problem.tsv")

  # if the problem is systematically targetting some individuals
  # n_distinct(allele.summary$INDIVIDUALS)
  # would be very small proportion of individuals.
  # But if the proportion on the total indi is large, its affecting all ind.
  # so this might be more a problem of not enough filtering...
  # need to figure out a way to translate this into something meaningful


  # Highlight problem with coverage info for the Alt allele...
  alt.no.depth <- dplyr::filter(allele.summary, is.na(ALLELE_ALT_DEPTH))

  if (nrow(alt.no.depth) > 1) {
    message("\nWarning:")
    message(paste(nrow(alt.no.depth), "genotypes from the blacklisted markers are missing ALT allele depth information"))
    message("This information is found in the table: $allele.summary")
  }

  cat("############################### RESULTS ###############################\n")
  res = list(
    allele.summary = allele.summary,
    alt.allele.count.distribution.plot = alt.allele.count.plot,
    alt.depth.count.distribution.plot = alt.depth.count.distribution.plot,
    number.markers = problem,
    blacklist.markers = blacklist.markers
  )
  message(paste0("Number of markers with allele copies below threshold: ", nrow(blacklist.markers)))
  timing <- proc.time() - timing
  message(stringi::stri_join("Computation time: ", round(timing[[3]]), " sec"))
  cat("############################## completed ##############################\n")
  return(res)
}
