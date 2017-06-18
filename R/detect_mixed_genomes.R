# detect mixed genomes
#' @title Detect mixed genomes
#' @description Highlight outliers individual's observed heterozygosity for a quick
#' diagnostic of mixed samples or poor polymorphism discovery due to DNA quality,
#' sequencing effort, etc.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.



#' @param ind.heterozygosity.threshold (string, double, optional)
#' Blacklist individuals based on observed heterozygosity (averaged across markers).
#' 
#' 
#' The string contains 2 thresholds values (min and max). 
#' The values are proportions (0 to 1), where 0 turns off the min threshold and
#' 1 turns off the max threshold.
#' Individuals with mean observed heterozygosity higher (>) or lower (<) 
#' than the thresholds will be blacklisted. 
#' 
#' Default: \code{ind.heterozygosity.threshold = NULL} will turn off completely 
#' the filter and the function will only output the plots and table of heterozygosity.

#' @details To help discard an individual based on his observed heterozygosity
#' (averaged across markers),
#' use the manhanttan plot to:
#' \enumerate{
#' \item contrast the individual with population and overall samples.
#' \item visualize the impact of missigness information (based on population or 
#' overall number of markers) and the individual observed heterozygosity. The
#' larger the point, the more missing genotypes.
#' }
#' \strong{Outlier above average:}
#' \itemize{
#' \item potentially represent two samples mixed together (action: blacklist), or...
#' \item a sample with more sequecing effort (point size small): did you merge your replicates fq files ? (action : keep and monitor)
#' \item a sample with poor sequencing effort (point size large) where the genotyped markers are
#' all heterozygotes, verify this with missingness (action: discard)
#' }
#' In all cases, if there is no bias in the population sequencing effort,
#' the size of the point will usually be "average" based on the population or
#' overall number of markers.
#' 
#' 
#' You can visualize individual observed heterozygosity, choose thresholds and
#' then visualize, choose thresholds and filter markers based on observed 
#' heterozygosity in one run with: \pkg{radiator} \code{\link{filter_het}}.

#' 
#' \strong{Outlier below average:}
#' \itemize{
#' \item A point with a size larger than the population or overall average (= lots of missing): 
#' the poor polymorphism discovery of the sample is probably the result of bad
#' DNA quality, a bias in sequencing effort, etc. (action: blacklist)
#' \item A point with a size that looks average (not much missing): 
#' this sample requires more attention (action: blacklist) and conduct more tests.
#' e.g. for biallelic data, look for coverage imbalance between ALT/REF allele.
#' At this point you need to distinguish between an artifact of poor polymorphism discovery 
#' or a biological reason (highly inbred individual, etc.).
#' }


#' @return The function returns inside the global environment a list with
#' 5 objects:
#' 
#' \enumerate{
#' \item the individual's heterozigosity (\code{$individual.heterozygosity})
#' a dataframe containing for each individual, the population id, the number of 
#' genotyped markers, the number of missing genotypes (based on the number of 
#' markers of the population and overall), the number of markers genotyped as heterozygote
#' and it's proportion based on the number of genotyped markers.
#' \item the heterozygosity statistics per populations and overall:\code{$heterozygosity.statistics}
#' \item the blacklisted individuals if \code{ind.heterozygosity.threshold} was selected: \code{$blacklist.ind.het}
#' \item the boxplot of individual heterozygosity:\code{$individual.heterozygosity.boxplot}
#' \item the manhattan plot of individual heterozygosity (\code{$individual.heterozygosity.manhattan.plot})
#' contrasted with missingness proportion based on the number of markers (population or overall).
#' The 2 facets will be identical when the dataset as common markers between
#' the populations. The dotted lines are the mean hetegozygosities.
#' }

#' @rdname detect_mixed_genomes
#' @export

#' @import ggplot2
#' @importFrom stringi stri_join stri_replace_all_fixed stri_sub
#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join funs
#' @importFrom readr write_tsv
#' @importFrom tibble data_frame
#' @importFrom tidyr complete gather unite spread nesting
#' @importFrom stats median sd

#' @examples
#' \dontrun{
#' #Step1: highlight outlier individuals, the simplest way to run:
#' 
#' outlier.ind.het <- radiator::detect_mixed_genomes(data = "wombat_tidy.tsv")
#' 
#' # Or if you don't have a tidy df:
#' 
#' outlier.ind.het <- radiator::tidy_genomic_data(
#' data = "wombat.vcf",
#' strata = "strata.wombat.tsv",
#' common.markers = FALSE,
#' vcf.metadata = FALSE,
#' verbose = TRUE
#' ) %>% 
#' radiator::detect_mixed_genomes(.)
#' 
#' 
#' #This example, without threshold, will not produce a blacklist of individuals.
#' 
#' #To look at the table with individual's heterozygosity:
#' 
#' outlier.ind.het$individual.heterozygosity
#' 
#' # To view the manhattan plot:
#' 
#' outlier.ind.het$individual.heterozygosity.manhattan.plot
#' 
#' # To view the box plot
#' outlier.ind.het$individual.heterozygosity.boxplot
#' 
#' # To save the boxplot:
#' 
#' ggsave(
#' "individual.heterozygosity.boxplot.pdf",
#' width = 15, height = 10, dpi = 600, units = "cm",
#' useDingbats = FALSE
#' )
#' 
#' # prefer a PNG:
#' 
#' ggsave(
#' "individual.heterozygosity.boxplot.png",
#' width = 15, height = 10, dpi = 300, units = "cm"
#' )
#' 
#' #Step2: blacklist ind with min and max Het obs thresholds
#' # Based on the look of the distribution using both jitter and boxplot,
#' # choose a threshold to blacklist the outliers and re-run the function.
#' 
#' outlier.ind.het <- radiator::detect_mixed_genomes(
#' data = "wombat_tidy.tsv", ind.heterozygosity.threshold = c(0.02, 0.2)
#' )
#' 
#' blacklist.ind.het <- outlier.ind.het$blacklist.ind.het
#' 
#' # To keep the blacklist:
#' readr::write_tsv(
#' x = blacklist.ind.het,
#' path = "blacklist.individuals.heterozygosity.tsv", col_names = TRUE
#' )
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com} and Peter Grewe \email{peter.grewe@@csiro.au}

detect_mixed_genomes <- function(
  data,
  ind.heterozygosity.threshold = NULL
) {
  cat("#######################################################################\n")
  cat("#################### radiator::detect_mixed_genomes #####################\n")
  cat("#######################################################################\n")
  timing <- proc.time()
  message("Analyzing data...")
  # manage missing arguments ---------------------------------------------------
  if (missing(data)) stop("missing data argument")
  
  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    input <- radiator::tidy_wide(data = data, import.metadata = FALSE)
  } else {
    input <- data
  }
  
  # check genotype column naming
  colnames(input) <- stringi::stri_replace_all_fixed(
    str = colnames(input), 
    pattern = "GENOTYPE", 
    replacement = "GT", 
    vectorize_all = FALSE
  )
  
  # necessary steps to make sure we work with unique markers and not duplicated LOCUS
  if (tibble::has_name(input, "LOCUS") && !tibble::has_name(input, "MARKERS")) {
    input <- dplyr::rename(.data = input, MARKERS = LOCUS)
  }
  
  # highlight heterozygote and missing (optimized for speed depending on input)
  # you see the difference with > 30K SNP
  
  n.markers.pop <- dplyr::filter(input, GT != "000000") %>% 
    dplyr::distinct(MARKERS, POP_ID) %>%
    dplyr::count(x = ., POP_ID)
  
  n.markers.overall <- dplyr::n_distinct(input$MARKERS[input$GT != "000000"])

  if (tibble::has_name(input, "GT_BIN")) {
    het.summary <- dplyr::mutate(
      .data = input,
      HET = dplyr::if_else(GT_BIN == 1, 1, 0, missing = 0)
      ) %>%
      dplyr::group_by(INDIVIDUALS) %>% 
      dplyr::mutate(GENOTYPED = length(GT_BIN[!is.na(GT_BIN)])) %>% 
      dplyr::ungroup(.)
    
  } else if (tibble::has_name(input, "GT_VCF")) {
    het.summary <- dplyr::mutate(
      .data = input,
      HET = dplyr::if_else(GT_VCF %in% c("1/0", "0/1"), 1, 0)
    ) %>%
      dplyr::group_by(INDIVIDUALS) %>%
      dplyr::mutate(GENOTYPED = length(INDIVIDUALS[GT_VCF != "./."])) %>% 
      dplyr::ungroup(.)
  } else {
    het.summary <- dplyr::mutate(
      .data = input,
      HET = dplyr::if_else(
        stringi::stri_sub(GT, 1, 3) != stringi::stri_sub(GT, 4, 6), 1, 0
      )
    ) %>% 
      dplyr::group_by(INDIVIDUALS) %>%
      dplyr::mutate(GENOTYPED = length(INDIVIDUALS[GT != "000000"])) %>% 
      dplyr::ungroup(.)
  }
  
  
  # Step 1. Highlight individual's heterozygosity  -----------------------------
  # Heterozygosity at the individual level before looking at the markers level per population
  # It's a good way to do outlier diagnostic ... mixed individuals
  
  # Create a new df with heterozygote info
  
  het.ind <- dplyr::select(.data = het.summary, POP_ID, INDIVIDUALS, HET, GENOTYPED) %>% 
    dplyr::full_join(n.markers.pop, by = "POP_ID") %>% 
    dplyr::group_by(INDIVIDUALS) %>% 
    dplyr::mutate(
      MISSING_PROP_POP = (n - GENOTYPED) / n,
      MISSING_PROP_OVERALL = (n.markers.overall - GENOTYPED) / n.markers.overall
    ) %>% 
    dplyr::group_by(INDIVIDUALS, POP_ID) %>% 
    dplyr::summarise(
      GENOTYPED = unique(GENOTYPED),
      MISSING_PROP_POP = unique(MISSING_PROP_POP),
      MISSING_PROP_OVERALL = unique(MISSING_PROP_OVERALL),
      HET_NUMBER = length(HET[HET == 1]),
      HET_PROP = HET_NUMBER / GENOTYPED
    ) %>%
    dplyr::arrange(POP_ID, HET_PROP) %>%
    dplyr::ungroup(.)
  
  het.ind.overall <- dplyr::mutate(.data = het.ind, POP_ID = as.character(POP_ID)) %>%
    dplyr::bind_rows(dplyr::mutate(.data = het.ind, POP_ID = rep("OVERALL", n()))) %>%
    dplyr::mutate(POP_ID = factor(POP_ID, levels = c(levels(het.ind$POP_ID), "OVERALL"))) %>% 
    tidyr::gather(data = ., key = MISSING_GROUP, value = MISSING_PROP, -c(POP_ID, INDIVIDUALS, GENOTYPED, HET_NUMBER, HET_PROP)) %>% 
    dplyr::mutate(MISSING_GROUP = factor(MISSING_GROUP, levels = c("MISSING_PROP_POP", "MISSING_PROP_OVERALL")))
  
  
  #Stats------------------------------------------------------------------------
  message("Calculating statistics")
  het.ind.stats <- het.ind.overall %>%
    dplyr::group_by(POP_ID) %>%
    dplyr::summarise(
      HET_MEAN = mean(HET_PROP, na.rm = TRUE),
      HET_MEDIAN = stats::median(HET_PROP, na.rm = TRUE),
      HET_SD = stats::sd(HET_PROP, na.rm = TRUE),
      HET_MIN = min(HET_PROP, na.rm = TRUE),
      HET_MAX = max(HET_PROP, na.rm = TRUE)
    ) %>% 
    dplyr::mutate_if(.tbl = ., .predicate = is.numeric, .funs = round, digits = 4) %>% 
    tidyr::unite(data = ., HET_RANGE, HET_MIN, HET_MAX, sep = " - ") %>% 
    dplyr::arrange(POP_ID, HET_MEAN)
  
  # Plots ----------------------------------------------------------------------
  message("Generating plots")
  
  rounder <- function(x, accuracy, f = round) {
    f(x / accuracy) * accuracy
  }
  y.breaks.by <- rounder(max(het.ind$HET_PROP, na.rm = TRUE)/10, 0.001, ceiling)
  y.breaks.max <- rounder(max(het.ind$HET_PROP, na.rm = TRUE), 0.001, ceiling)
  y.breaks <- seq(0, y.breaks.max + y.breaks.by , by = y.breaks.by)
  
  # labeller to rename in the facet_grid or facet_wrap call:
  facet_names <- ggplot2::as_labeller(c(`MISSING_PROP_OVERALL` = "Missing (overall)", `MISSING_PROP_POP` = "Missing (populations)"))
  
  individual.heterozygosity.manhattan.plot <- ggplot2::ggplot(data = het.ind.overall, ggplot2::aes(x = POP_ID, y = HET_PROP, size = MISSING_PROP, colour = POP_ID)) + 
    ggplot2::geom_jitter(alpha = 0.6) + 
    ggplot2::labs(y = "Individual's Mean Observed Heterozygosity (proportion)") +
    # labs(x = "Populations") +
    # labs(colour = "Populations") +
    ggplot2::scale_y_continuous(name = ggplot2::waiver(), breaks = y.breaks) +#, limits = c(0, y.breaks.max), expand = c(0.1, 0)) +
    ggplot2::scale_color_discrete(guide = "none") +
    ggplot2::scale_size_continuous(name = "Missing proportion") +
    # theme_minimal() +
    ggplot2::theme(
      # legend.position = "none",
      # panel.grid.major.y = element_line(linetype = "solid"),
      # panel.grid.minor.y = element_line(linetype = "longdash", size = 1),
      # panel.background = element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      # axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"), 
      # axis.text.x = element_text(size = 10, family = "Helvetica"),
      axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"), 
      axis.text.y = ggplot2::element_text(size = 8, family = "Helvetica")
    ) +
    ggplot2::geom_hline(mapping = ggplot2::aes(yintercept = HET_MEAN), het.ind.stats, linetype = "dotted", size = 0.6) + #mean
    # geom_hline(mapping = aes(yintercept = HET_sig_minus), het.ind.stats.pop, linetype = "dashed") + #3 sigma -
    # geom_hline(mapping = aes(yintercept = HET_sig_plus), het.ind.stats.pop, linetype = "dashed") + #3 sigma +
    ggplot2::facet_grid(MISSING_GROUP ~ POP_ID, switch = "x", scales = "free", labeller = ggplot2::labeller(MISSING_GROUP = facet_names))
  # individual.heterozygosity.manhattan.plot
  
  individual.heterozygosity.boxplot <- ggplot2::ggplot(data = het.ind.overall, ggplot2::aes(x = POP_ID, y = HET_PROP, colour = POP_ID)) + 
    ggplot2::geom_boxplot() + 
    ggplot2::labs(y = "Individual's Mean Observed Heterozygosity (proportion)") +
    ggplot2::labs(x = "Populations") +
    ggplot2::labs(colour = "Populations") +
    ggplot2::scale_y_continuous(name = ggplot2::waiver(), breaks = y.breaks, limits = c(0, y.breaks.max), expand = c(0.06, 0)) +
    ggplot2::theme_classic() +
    # theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"), 
      axis.text.x = ggplot2::element_text(size = 10, family = "Helvetica"),
      axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"), 
      axis.text.y = ggplot2::element_text(size = 8, family = "Helvetica")
    )
  # individual.heterozygosity.boxplot
  
  ## Step 2: Blacklist outlier individuals -------------------------------------
  # Blacklist individuals based a threshold of mean heterozygosity
  if (!is.null(ind.heterozygosity.threshold)) {
    # ind.heterozygosity.threshold <- c(0.035, 0.10)
    threshold.min <- ind.heterozygosity.threshold[1]
    threshold.max <- ind.heterozygosity.threshold[2]
    
    blacklist.ind.het  <- dplyr::ungroup(het.ind) %>%
      dplyr::filter(HET_PROP > threshold.max | HET_PROP < threshold.min) %>% 
      dplyr::distinct(INDIVIDUALS)
    
    message(stringi::stri_join("Filter individual's heterozygosity: ", length(blacklist.ind.het$INDIVIDUALS), " individual(s) blacklisted"))
    
  } else {
    blacklist.ind.het <- "ind.heterozygosity.threshold is necessary to get a blacklist of individuals"
  }
  message(stringi::stri_join("Computation time: ", round((proc.time() - timing)[[3]]), " sec"))
  cat("############################## completed ##############################\n")
  res <- list(
    individual.heterozygosity = het.ind,
    heterozygosity.statistics = het.ind.stats,
    blacklist.ind.het = blacklist.ind.het,
    individual.heterozygosity.boxplot = individual.heterozygosity.boxplot,
    individual.heterozygosity.manhattan.plot = individual.heterozygosity.manhattan.plot
  )
  return(res)
}
