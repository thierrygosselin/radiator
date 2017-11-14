# Observed heterozygosity
#' @title Heterozygosity filter
#' @description Observed Heterozygosity based filtering.
#' The filter arguments of \code{filter_het} allows you to test rapidly
#' if departure from realistic expectations of heterozygosity statistics
#' are a problem in downstream analysis.

#' \enumerate{
#' \item Highlight outliers individual's observed heterozygosity for a quick
#' diagnostic of mixed samples or poor polymorphism discovery. The statistic is
#' also contrasted with missing data to help differentiate problems.
#' (also computed in \link{detect_mixed_genomes})
#'
#' \item The observed heterozygosity in the dataset: an assembly artefact,
#' a genotyping problem, a problem of population groupings
#' or a reliable signal of biological polymorphism?
#' Detect assembly artifact or genotyping problem (e.g. under/over-splitting loci)
#' by looking at marker's observed heterozygosity statistics by population
#' or overall.
#'
#' \item Use haplotype or snp level statistics. When the haplotype approach
#' is selected, consistencies of heterozygosity statistics along the read are
#' highlighted.
#'
#' \item Interactive approach help by visualizing the data before making
#' a decision on thresholds.
#' }

# Most arguments are inherited from tidy_genomic_data
#' @inheritParams tidy_genomic_data

#' @param interactive.filter (optional, logical) Do you want the filtering session to
#' be interactive. With the default, the user is asked to see figures of
#' distribution before making decisions for filtering with heterozygosity statistics.
#' Default: \code{interactive.filter == TRUE}.

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

#' @param het.approach (Character string). First value, \code{"SNP"} or
#' \code{"haplotype"}. The haplotype approach considers the statistic
#' consistency on the read (locus/haplotype). The major difference:
#' the haplotype approach results in blacklisting the entire locus/haplotype
#' with all the SNPs on the read.
#' With the SNP approach, SNPs are independently analyzed and blacklisted.
#' The second value, will use the statistics by population \code{"pop"} or
#' will consider the data overall, as 1 large group \code{"overall"}.
#' Default: \code{het.approach = c("SNP", "overall")}.

#' @param het.threshold Number Biallelic markers usually max 0.5. But departure
#' from this value is common. The higher the proportion threshold, the more relaxed
#' the filter is. With default there is no filtering.
#' Default: \code{het.threshold = 1}.

#' @param het.dif.threshold Number (0 - 1). For \code{het.approach = "haplotype"} only.
#' You can set a threshold for the difference in het along your read.
#' Set the number your willing to tolerate on the same read/haplotype.
#' e.g. if you have 2 SNP on a read/haplotype and on as a het of 0.9 and the other
#' 0.1 and you set \code{het.dif.threshold = 0.3}, this markers will be blacklisted.
#' You should strive to have similar statistics along short read like RAD.
#' The higher the proportion threshold, the more relaxed the filter is.
#' With default, there is no filtering and all the range of differences are allowed.
#' Default: \code{het.dif.threshold = 1}.

#' @param outlier.pop.threshold (integer, optional) Useful to incorporate problematic
#' populations dragging down polymorphism discovery, but still wanted for analysis.
#' Use this threshold to allow variance in the number of populations passing
#' the thresholds described above.
#' e.g. with \code{outlier.pop.threshold = 2},
#' you tolerate a maximum of 2 populations failing the
#' \code{het.threshold} and/or \code{het.dif.threshold}.
#' Manage outlier markers, individuals and populations
#' downstream with blacklists and whitelists produced by the function.
#' Default: \code{outlier.pop.threshold = 1}. See details for more info.

#' @param coverage.info (optional, logical) Use
#' \code{coverage.info = TRUE}, if you want to visualize
#' the relationshio between locus coverage and locus observed heterozygosity
#' statistics. Coverage information is required (e.g. in a vcf file...).
#' Default: {coverage.info = FALSE}.

#' @param helper.tables (logical) Output tables that show
#' the number of markers blacklisted or whitelisted based on a series of
#' automatic thresholds to guide decisions. When \code{interactive.filter == TRUE},
#' helper tables are written to the directory.
#' Default: \code{helper.tables = FALSE}.

#' @param filename (optional) The function uses \code{\link[fst]{write.fst}},
#' to write the tidy data frame in
#' the folder created in the working directory. The file extension appended to
#' the \code{filename} provided is \code{.rad}.
#' Default: \code{filename = NULL}.

#' @details
#' \strong{Interactive version}
#'
#' There are 4 steps in the interactive version to visualize and filter
#' the data based heterozygosity statistics:
#'
#' Step 1. Individual's observed heterozygosity: outliers that might represent mixed samples
#' Step 2. Blacklist outliers based on a proportion threshold of mean observed heterozygosity
#' Step 3. Observed heterozygosity statistics per populations and overall
#' Step 4: Blacklist markers based on observed heterozygosity
#'
#'
#'
#' \strong{outlier.pop.threshold}
#'
#' If your a regular radiator user, you've seen the \code{pop.num.threshold}.
#' \code{outlier.pop.threshold}, is different and requires more thinking,
#' because the number of populations genotyped potentially vary across markers,
#' which makes the use of \code{pop.num.threshold} less optimal.
#' e.g. If only 4 populations out of 10 are genotyped, for a marker you want
#' to keep, using \code{pop.num.threshold = 0.6} will lead to unwanted results
#' and inconsistensis.
#' It's easier to tolerate outliers with this new approach:
#' \code{outlier.pop.threshold = 2}.
#'
#'
#' \strong{Individual observed heterozygosity (averaged across markers):}
#' To help discard an individual based on his observed heterozygosity
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
#' 15 objects, the objects names are found by using \code{names(your.list.name)}:
#'
#' \enumerate{
#' \item filtered tidy data frame: \code{$tidy.filtered.het}
#' \item whitelist of markers:\code{$whitelist.markers}
#' \item the strata:\code{$strata}
#' \item the filters parameters used:\code{$filters.parameters}
#' \item the individual's heterozigosity:\code{$individual.heterozigosity}
#' \item the heterozygosity statistics per populations and overall:\code{$heterozygosity.statistics}
#' \item the blacklisted individuals based on the individual's heterozigosity:\code{$blacklist.ind.het}
#' \item a list containing the helper tables:\code{$helper.table.het}
#' \item the boxplot of individual observed heterozygosity:\code{$individual.heterozygosity.boxplot}
#' \item the manhattan plot of individual heterozygosity:\code{$individual.heterozygosity.manhattan.plot}
#' \item the boxplot of observed heterozygosity averaged across markers and pop:\code{$markers.pop.heterozygosity.boxplot}
#' \item the density plot of observed heterozygosity averaged across markers and pop:\code{$markers.pop.heterozygosity.density.plot}
#' \item the manhattan plot of observed heterozygosity averaged across markers and pop:\code{$markers.pop.heterozygosity.manhattan.plot}
#' \item the relationship between the number of SNPs on a locus and locus observed heterozygosity statistics:\code{$snp.per.locus.het.plot}
#' \item the relationship between locus coverage (read depth) and locus observed heterozygosity statistics:\code{$het.cov.plot}
#' \item whitelist of markers:\code{$whitelist.markers}
#' \item whitelist of markers:\code{$whitelist.markers}
#' }
#'
#'
#' In the working directory, output is conditional to \code{interactive.filter} argument

#' @rdname filter_het
#' @export

#' @importFrom stringi stri_join stri_replace_all_fixed stri_sub
#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join funs
#' @importFrom readr write_tsv
#' @importFrom tibble data_frame
#' @importFrom tidyr complete gather unite spread nesting
#' @importFrom ggplot2 ggplot geom_jitter labs aes scale_y_continuous as_labeller waiver scale_color_discrete scale_size_continuous theme element_blank element_text geom_hline facet_grid labeller geom_boxplot theme_classic geom_line expand_limits ggsave theme_minimal element_line

#' @seealso \link{plot_density_distribution_het}


filter_het <- function(
  data,
  interactive.filter = TRUE,
  ind.heterozygosity.threshold = NULL,
  het.approach = c("SNP", "overall"),
  het.threshold = 1,
  het.dif.threshold = 1,
  outlier.pop.threshold = 1,
  helper.tables = FALSE,
  coverage.info = FALSE,
  filename = NULL,
  vcf.metadata = TRUE,
  blacklist.id = NULL,
  blacklist.genotype = NULL,
  whitelist.markers = NULL,
  monomorphic.out = TRUE,
  max.marker = NULL,
  snp.ld = NULL,
  common.markers = FALSE,
  strata = NULL,
  pop.levels = NULL,
  pop.labels = NULL,
  pop.select = NULL,
  parallel.core = parallel::detectCores() - 1
) {
  cat("#######################################################################\n")
  cat("######################### radiator::filter_het ##########################\n")
  cat("#######################################################################\n")
  opt.change <- getOption("width")
  options(width = 70)
  timing <- proc.time()

  # manage missing arguments -----------------------------------------------------
  if (missing(data)) stop("Input file missing")
  if (!is.null(pop.levels) & is.null(pop.labels)) {
    pop.levels <- stringi::stri_replace_all_fixed(pop.levels, pattern = " ", replacement = "_", vectorize_all = FALSE)
    pop.labels <- pop.levels
  }

  if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")

  if (!is.null(pop.labels)) {
    if (length(pop.labels) != length(pop.levels)) stop("pop.labels and pop.levels must have the same length (number of groups)")
    pop.labels <- stringi::stri_replace_all_fixed(pop.labels, pattern = " ", replacement = "_", vectorize_all = FALSE)
  }

  if (!is.null(pop.select)) {
    pop.select <- stringi::stri_replace_all_fixed(pop.select, pattern = " ", replacement = "_", vectorize_all = FALSE)
  }
  # Message about steps taken during the process ---------------------------------
  if (interactive.filter) {
    message("Interactive mode: on, 4 steps of data visualization and filtering based on heterozygosity:\n")
    message("Step 1. Individual's heterozygosity: outliers that might represent mixed samples")
    message("Step 2. Blacklist outliers based on a proportion threshold of mean heterozygosity")
    message("Step 3. Observed heterozygosity statistics per populations and overall")
    message("Step 4: Blacklist markers based on observed heterozygosity\n\n")
  }
  # Folder -------------------------------------------------------------------
  # Date and time
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  folder.extension <- stringi::stri_join("filter_het_", file.date, sep = "")
  path.folder <- stringi::stri_join(getwd(),"/", folder.extension, sep = "")
  dir.create(file.path(path.folder))
  message("\nFolder created: ", folder.extension)
  file.date <- NULL #unused object

  # Filter parameter file ------------------------------------------------------
  message("Parameters used in this run are stored in a file")
  filters.parameters <- list.files(path = getwd(), pattern = "filters_parameters.tsv", full.names = TRUE)
  if (length(filters.parameters) == 0) {
    filters.parameters <- tibble::data_frame(FILTERS = as.character(), PARAMETERS = as.character(), VALUES = as.integer(), BEFORE = as.character(), AFTER = as.character(), BLACKLIST = as.integer(), UNITS = as.character(), COMMENTS = as.character())
    readr::write_tsv(x = filters.parameters, path = "filters_parameters.tsv", append = FALSE, col_names = TRUE)
    message("    Created a parameter file: filters_parameters.tsv")
  } else {
    message("    Using the filters parameters file: filters_parameters.tsv")
  }
  # File type detection----------------------------------------------------------
  data.type <- radiator::detect_genomic_format(data)

  if (data.type == "haplo.file") {
    message("With stacks haplotype file the approach is automatically set to: haplotype")
    # counter intuitive, but there is no snp or haplotype info in a stacks haplotypes file
    het.approach[1] <- "SNP"
  }

  if (coverage.info) {
    vcf.metadata <- TRUE
  }

  # import data ----------------------------------------------------------------
  message("Importing data ...")
  input <- radiator::tidy_genomic_data(
    data = data,
    vcf.metadata = vcf.metadata,
    blacklist.id = blacklist.id,
    blacklist.genotype = blacklist.genotype,
    whitelist.markers = whitelist.markers,
    monomorphic.out = monomorphic.out,
    max.marker = max.marker,
    snp.ld = snp.ld,
    common.markers = common.markers,
    strata = strata,
    pop.levels = pop.levels,
    pop.labels = pop.labels,
    pop.select = pop.select,
    filename = NULL,
    parallel.core = parallel.core
  )
  # create a strata.df
  strata.df <- input %>%
    dplyr::select(INDIVIDUALS, POP_ID) %>%
    dplyr::distinct(INDIVIDUALS, .keep_all = TRUE)

  pop.levels <- levels(input$POP_ID)
  pop.labels <- pop.levels

  pop.number <- dplyr::n_distinct(input$POP_ID) # number of pop

  # Biallelic data ?
  biallelic <- radiator::detect_biallelic_markers(data = input, verbose = FALSE)

  # scan for the columnn CHROM and keep the info to include it after imputations
  if (tibble::has_name(input, "CHROM")) {
    marker.meta <- dplyr::distinct(.data = input, MARKERS, CHROM, LOCUS, POS)
  } else {
    marker.meta <- NULL
  }

  # double check the approach vs the file used
  if (!tibble::has_name(input, "CHROM") & het.approach[1] == "haplotype") {
    stop("The haplotype approach during HET filtering requires LOCUS and SNP
information")
  }

  # Step 1. Highlight individual's heterozygosity  -----------------------------
  # Heterozygosity at the individual level before looking at the markers level per population
  # It's a good way to do outlier diagnostic ... mixed individuals

  # highlight heterozygote and missing (optimized for speed depending on input)
  # you see the difference with > 30K SNP

  n.markers.pop <- dplyr::filter(input, GT != "000000") %>%
    dplyr::distinct(MARKERS, POP_ID) %>%
    dplyr::count(x = ., POP_ID)

  # n.markers.overall <- dplyr::n_distinct(input$MARKERS[input$GT != "000000"])
  n.markers.overall <- dplyr::n_distinct(input$MARKERS)

  if (tibble::has_name(input, "GT_BIN")) {
    het.summary <- dplyr::mutate(
      .data = input,
      HET = dplyr::if_else(GT_BIN == 1, 1, 0, missing = 0)
    ) %>%
      dplyr::group_by(INDIVIDUALS) %>%
      dplyr::mutate(GENOTYPED = length(GT_BIN[!is.na(GT_BIN)])) %>%
      dplyr::ungroup(.)

  } else if (biallelic && tibble::has_name(input, "GT_VCF")) {
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

  # Get stats...
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

  message("Generating plots")

  rounder <- function(x, accuracy, f = round) {
    f(x / accuracy) * accuracy
  }
  y.breaks.by <- rounder(max(het.ind$HET_PROP)/10, 0.001, ceiling)
  y.breaks.max <- rounder(max(het.ind$HET_PROP), 0.001, ceiling)
  y.breaks <- seq(0, y.breaks.max + y.breaks.by , by = y.breaks.by)

  # labeller to rename in the facet_grid or facet_wrap call:
  facet_names <- ggplot2::as_labeller(
    c(`MISSING_PROP_OVERALL` = "Missing (overall)", `MISSING_PROP_POP` = "Missing (populations)"))

  individual.heterozygosity.manhattan.plot <- ggplot2::ggplot(
    data = het.ind.overall,
    ggplot2::aes(x = POP_ID, y = HET_PROP, size = MISSING_PROP, colour = POP_ID)) +
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

  individual.heterozygosity.boxplot <- ggplot2::ggplot(
    data = het.ind.overall,
    ggplot2::aes(x = POP_ID, y = HET_PROP, colour = POP_ID)) +
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

  if (interactive.filter) {
    message("\nStep 1. Individual's heterozygosity visualization")
    message("    Inspect the plot for outliers that might represent
    mixed or pooled samples, on the higher end of the spectrum, or
    poorly sequenced individuals, on the lower end.\n")
    print(individual.heterozygosity.manhattan.plot)
    message("\n    If pooled samples with higher heterozygosity are
    problematic, use one of the replicate instead.")
  }

  # save
  ggplot2::ggsave(stringi::stri_join(path.folder, "/individual.heterozygosity.manhattan.plot.pdf"), width = pop.number * 4, height = 15, dpi = 600, units = "cm", useDingbats = FALSE)
  ggplot2::ggsave(stringi::stri_join(path.folder, "/individual.heterozygosity.manhattan.plot.png"), width = pop.number * 4, height = 15, dpi = 300, units = "cm")
  message("\n    2 versions (pdf and png) of the plot (individual.heterozygosity.manhattan.plot)
    were written in the folder")

  if (interactive.filter) {
    message("\n\n    Same data, but shown with a box plot of individual's heterozygosity:")
    # boxplot <- as.character(readLines(n = 1))
    # if (boxplot == "y") {
    print(individual.heterozygosity.boxplot)
  }
  # save
  ggplot2::ggsave(stringi::stri_join(path.folder, "/individual.heterozygosity.boxplot.pdf"), width = pop.number * 4, height = 10, dpi = 600, units = "cm", useDingbats = FALSE)
  ggplot2::ggsave(stringi::stri_join(path.folder, "/individual.heterozygosity.boxplot.png"), width = pop.number * 4, height = 10, dpi = 300, units = "cm")
  message("    2 versions (pdf and png) of the plot (individual.heterozygosity.boxplot)
    were written in the folder")

  message("\nNote: in Rstudio it's easy to visualize the 2 plot produced and going back and forth in the plot window")

  ## Step 2: Blacklist outlier individuals -------------------------------------
  if (!is.null(ind.heterozygosity.threshold)) {
    threshold.min <- ind.heterozygosity.threshold[1]
    threshold.max <- ind.heterozygosity.threshold[2]
  }

  # Blacklist individuals based a threshold of mean heterozygosity
  if (interactive.filter) {
    threshold.min <- 2#to make sure we go through selection of threshold
    while (isTRUE(threshold.min > 1)) {
      message("\n\nStep 2. Blacklist outliers based on a proportion threshold of mean observed heterozygosity\n")
      message("If you want to blacklist individuals based on the mean observed heterozygosity,
use this ind.heterozygosity.threshold filter. We have 2 values (min and max) to select.
First, select the min threshold, individuals with mean observed heterozygosity
lower (<) than the threshold will be blacklisted. 0 turns off the min filter.
Enter the values (proportion, e.g. 0.05 for the min threshold:")
      threshold.min <- as.numeric(readLines(n = 1))
    }
  }

  if (interactive.filter) {
    threshold.max <- 2#to make sure we go through selection of threshold
    while (isTRUE(threshold.max > 1)) {
      message("Second, the max threshold...
Individuals with mean observed heterozygosity higher (>) than the thresholds
will be blacklisted. 1 turns off the max filter.
Enter the values (proportion, e.g. 0.34 for max threshold:")
      threshold.max <- as.numeric(readLines(n = 1))
    }
    ind.heterozygosity.threshold <- c(threshold.min, threshold.max)
  }

  if (!is.null(ind.heterozygosity.threshold)) {
    blacklist.ind.het  <- dplyr::ungroup(het.ind) %>%
      dplyr::filter(HET_PROP > threshold.max | HET_PROP < threshold.min) %>%
      dplyr::distinct(INDIVIDUALS)
    message("Filter individual's heterozygosity: ", length(blacklist.ind.het$INDIVIDUALS), " individual(s) blacklisted")
  }


  # Remove the individuals from the dataset
  if (length(blacklist.ind.het$INDIVIDUALS > 0)) {
    readr::write_tsv(
      x = blacklist.ind.het,
      path = stringi::stri_join(path.folder, "/blacklist.individuals.heterozygosity.tsv"),
      col_names = TRUE
    )
    message("Blacklist (blacklist.individuals.heterozygosity) written in the folder")
    het.summary <- dplyr::anti_join(het.summary, blacklist.ind.het, by = "INDIVIDUALS")
  }

  # Step 3. Markers observed heterozygosity statistics per populations and overall----
  # input <- NULL # no longer needed
  # decide approach
  # Haplotype or SNP approach ...
  if (interactive.filter) {
    message("\n\nStep 3. Markers observed heterozygosity statistics per populations and overall\n")

    het.approach[1] <- "not defined"
    while (isTRUE(het.approach[1] != "haplotype" & het.approach[1] != "SNP")) {# to make sure the answer is ok
      message("het.approach argument: haplotype or SNP ?\n")
      message("Decide on the best approach to filter the data based on oserved heterozygosity\n")
      message("The haplotype approach considers the statistic consistency on the read (locus/haplotype).
The major difference: the haplotype approach results in blacklisting the entire
locus/haplotype with all the SNPs on the read. With the SNP approach, SNPs are
independently analyzed and blacklisted.\n")
      message("Enter if you want the haplotype or SNP approach (haplotype/SNP):")
      het.approach[1] <- as.character(readLines(n = 1))
    }
  }
  # overall or by populations...
  # het.approach[2] <- "test" #test
  if (interactive.filter) {
    het.approach[2] <- "not defined"
    while (isTRUE(het.approach[2] != "overall" & het.approach[2] != "pop")) {# to make sure the answer is ok
      message("het.approach argument: overall or by pop ?\n")
      message("Decide on the best approach to filter the data based on observed heterozygosity\n")
      message("The overall approach doesn't look at the observed heterozygosity by populations.
It's just 1 large population of sample.
If you're not sure about your population structure and/or population sampling,
use the overall approach.\n")
      message("Enter if you want the overall or pop approach (overall/pop):")
      het.approach[2] <- as.character(readLines(n = 1))
    }
  }
  message("Computing heterozygosity summary statistics...")

  # For the next steps we keep only genotyped markers
  het.summary <- dplyr::filter(.data = het.summary, GT != "000000")

  if (tibble::has_name(het.summary, "LOCUS") & het.approach[1] == "haplotype") {
    # By pop
    het.missing.pop <- dplyr::distinct(input, LOCUS, POP_ID, INDIVIDUALS, GT) %>%
      dplyr::group_by(LOCUS, POP_ID) %>%
      dplyr::summarise(
        GENOTYPED = length(INDIVIDUALS[GT != "000000"]),
        TOTAL = n(),
        MISSING_PROP = (TOTAL - GENOTYPED)/TOTAL
      ) %>%
      dplyr::ungroup(.) %>%
      dplyr::filter(MISSING_PROP != 1) %>%
      dplyr::select(LOCUS, POP_ID, MISSING_PROP)


    het.missing.overall <- dplyr::distinct(input, LOCUS, INDIVIDUALS, GT) %>%
      dplyr::group_by(LOCUS) %>%
      dplyr::summarise(
        GENOTYPED = length(INDIVIDUALS[GT != "000000"]),
        TOTAL = n(),
        MISSING_PROP = (TOTAL - GENOTYPED)/TOTAL
      ) %>%
      dplyr::select(LOCUS, MISSING_PROP) %>%
      dplyr::mutate(POP_ID = rep("OVERALL", n())) %>%
      dplyr::bind_rows(dplyr::mutate(.data = het.missing.pop, POP_ID = as.character(POP_ID))) %>%
      dplyr::mutate(POP_ID = factor(POP_ID, levels = c(levels(het.missing.pop$POP_ID), "OVERALL")))

    het.missing.pop <- NULL

    het.summary.pop <- het.summary %>%
      dplyr::group_by(MARKERS, LOCUS, POP_ID) %>%
      dplyr::summarise(HET_O = as.numeric(length(HET[HET == 1]) / n())) %>%
      dplyr::group_by(LOCUS, POP_ID) %>%
      dplyr::summarise(
        HET_MEAN = mean(HET_O),
        HET_MAX = max(HET_O),
        HET_MIN = min(HET_O),
        HET_DIF = HET_MAX - HET_MIN
      )

    # overall
    het.summary.overall <- het.summary %>%
      dplyr::group_by(MARKERS, LOCUS) %>%
      dplyr::summarise(HET_O = as.numeric(length(HET[HET == 1]) / n())) %>%
      dplyr::group_by(LOCUS) %>%
      dplyr::summarise(
        HET_MEAN = mean(HET_O),
        HET_MAX = max(HET_O),
        HET_MIN = min(HET_O),
        HET_DIF = HET_MAX - HET_MIN
      ) %>%
      dplyr::mutate(POP_ID = rep("OVERALL", n()))



    het.summary.tidy <- suppressWarnings(
      dplyr::bind_rows(het.summary.pop, het.summary.overall) %>%
        dplyr::full_join(het.missing.overall, by = c("LOCUS", "POP_ID")) %>%
        tidyr::gather(
          data = .,
          key = HET_GROUP,
          value = VALUE,
          -c(LOCUS, POP_ID, MISSING_PROP)
        ) %>%
        dplyr::rename(MARKERS = LOCUS) %>%
        dplyr::mutate(
          HET_GROUP = factor(
            HET_GROUP,
            levels = c("HET_MEAN", "HET_MIN", "HET_MAX", "HET_DIF"),
            ordered = TRUE
          )))
  } else {
    # By pop
    het.missing.pop <- input %>%
      dplyr::group_by(MARKERS, POP_ID) %>%
      dplyr::summarise(
        GENOTYPED = length(INDIVIDUALS[GT != "000000"]),
        TOTAL = n(),
        MISSING_PROP = (TOTAL - GENOTYPED)/TOTAL
      ) %>%
      dplyr::ungroup(.) %>%
      # next step removes pop with 100% missing (not genotyped at all)
      # no impact if common.markers was used.
      dplyr::filter(MISSING_PROP != 1) %>%
      dplyr::select(MARKERS, POP_ID, MISSING_PROP)

    het.summary.pop <- het.summary %>%
      dplyr::group_by(MARKERS, POP_ID) %>%
      dplyr::summarise(HET_MEAN = as.numeric(length(HET[HET == 1]) / n())) %>%
      dplyr::full_join(het.missing.pop, by = c("MARKERS", "POP_ID"))

    # Overall
    het.missing.overall <- input %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::summarise(
        GENOTYPED = length(INDIVIDUALS[GT != "000000"]),
        TOTAL = n(),
        MISSING_PROP = (TOTAL - GENOTYPED)/TOTAL
      ) %>%
      dplyr::select(MARKERS, MISSING_PROP) %>%
      dplyr::mutate(POP_ID = rep("OVERALL", n()))

    het.summary.overall <- het.summary %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::summarise(HET_MEAN = as.numeric(length(HET[HET == 1]) / n())) %>%
      dplyr::mutate(POP_ID = rep("OVERALL", n())) %>%
      dplyr::full_join(het.missing.overall, by = c("MARKERS", "POP_ID"))


    het.summary.tidy <- suppressWarnings(
      dplyr::bind_rows(het.summary.pop, het.summary.overall) %>%
        tidyr::gather(
          data = .,
          key = HET_GROUP,
          value = VALUE,
          -c(MARKERS, POP_ID, MISSING_PROP)
        )
    )
  }
  # Tidy use for figures
  # # unused object
  het.summary <- het.summary.pop
  het.summary.pop <- NULL


  if (!is.null(pop.levels)) {
    het.summary.tidy <- het.summary.tidy %>%
      dplyr::mutate(POP_ID = factor(POP_ID, levels = c(pop.labels, "OVERALL"), ordered = TRUE)) %>%
      dplyr::arrange(POP_ID)
  }

  # Visualization --------------------------------------------------------------
  # Density plot markers het obs -------------------------------------------------
  markers.pop.heterozygosity.density.plot <- ggplot2::ggplot(
    het.summary.tidy, ggplot2::aes(x = VALUE, na.rm = FALSE)) +
    ggplot2::geom_line(ggplot2::aes(y = ..scaled..), stat = "density", adjust = 0.3) +
    ggplot2::labs(x = "Markers Observed Heterozygosity") +
    ggplot2::labs(y = "Density of SNP (scaled)") +
    ggplot2::expand_limits(y = 0) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
      legend.title = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
      legend.text = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
      strip.text.y = ggplot2::element_text(angle = 0, size = 12, family = "Helvetica", face = "bold"),
      strip.text.x = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold")
    ) +
    ggplot2::facet_grid(POP_ID ~ HET_GROUP)

  # manhattan plot markers het obs -------------------------------------------------
  markers.pop.heterozygosity.manhattan.plot <- ggplot2::ggplot(
    data = dplyr::filter(het.summary.tidy, HET_GROUP == "HET_MEAN"),
    ggplot2::aes(x = POP_ID, y = VALUE, colour = POP_ID, size = MISSING_PROP)) +
    ggplot2::geom_jitter(alpha = 0.3) +
    ggplot2::labs(y = "Markers Observed Heterozygosity (mean)") +
    ggplot2::labs(x = "Populations") +
    ggplot2::labs(colour = "Populations") +
    ggplot2::scale_color_discrete(guide = "none") +
    ggplot2::scale_size_continuous(name = "Missing proportion") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(linetype = "solid", colour = "black"),
      panel.grid.minor.y = ggplot2::element_line(linetype = "dotted", colour = "blue"),
      panel.grid.major.x = ggplot2::element_blank(),
      # legend.position = "none",
      axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.text.x = ggplot2::element_text(size = 10, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5),
      axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.text.y = ggplot2::element_text(size = 8, family = "Helvetica")
    ) #+ ggplot2::facet_grid(~ HET_GROUP)
  # markers.pop.heterozygosity.manhattan.plot

  # box plot markers het obs -------------------------------------------------
  markers.pop.heterozygosity.boxplot <- ggplot2::ggplot(
    data = het.summary.tidy, ggplot2::aes(x = POP_ID, y = VALUE, colour = POP_ID)) +
    ggplot2::geom_boxplot() +
    ggplot2::labs(y = "Markers Observed Heterozygosity") +
    ggplot2::labs(x = "Populations") +
    ggplot2::labs(colour = "Populations") +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.text.x = ggplot2::element_text(size = 10, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5),
      axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.text.y = ggplot2::element_text(size = 8, family = "Helvetica")
    ) +
    ggplot2::facet_grid(~ HET_GROUP, scales = "free_y")

  het.summary.tidy <- NULL #unused object

  # helper tables ----------------------------------------------------------------

  if (helper.tables) {
    message("Generating several helper table(s)...")
    if (tibble::has_name(het.summary, "LOCUS") &  het.approach[1] == "haplotype") {# by Haplotype
      n.markers <- dplyr::n_distinct(het.summary$LOCUS)

      # by pop
      het.helper <- dplyr::select(het.summary, LOCUS, POP_ID, HET_MAX, HET_DIF) %>%
        dplyr::filter(POP_ID != "OVERALL") %>%
        dplyr::mutate(
          MAX_0.4 = dplyr::if_else(HET_MAX <= 0.4, "whitelist", "blacklist"),
          MAX_0.5 = dplyr::if_else(HET_MAX <= 0.5, "whitelist", "blacklist"),
          MAX_0.6 = dplyr::if_else(HET_MAX <= 0.6, "whitelist", "blacklist"),
          MAX_0.7 = dplyr::if_else(HET_MAX <= 0.7, "whitelist", "blacklist"),
          MAX_0.8 = dplyr::if_else(HET_MAX <= 0.8, "whitelist", "blacklist"),
          MAX_0.9 = dplyr::if_else(HET_MAX <= 0.9, "whitelist", "blacklist")
        ) %>%
        tidyr::gather(
          data = .,
          key = MAX_THRESHOLD,
          value = MAX_OUTLIERS,
          MAX_0.4:MAX_0.9
        ) %>%
        dplyr::mutate(
          DIF_0.1 = dplyr::if_else(HET_DIF <= 0.1, "whitelist", "blacklist"),
          DIF_0.2 = dplyr::if_else(HET_DIF <= 0.2, "whitelist", "blacklist"),
          DIF_0.3 = dplyr::if_else(HET_DIF <= 0.3, "whitelist", "blacklist"),
          DIF_0.4 = dplyr::if_else(HET_DIF <= 0.4, "whitelist", "blacklist"),
          DIF_0.5 = dplyr::if_else(HET_DIF <= 0.5, "whitelist", "blacklist"),
          DIF_0.6 = dplyr::if_else(HET_DIF <= 0.6, "whitelist", "blacklist"),
          DIF_0.7 = dplyr::if_else(HET_DIF <= 0.7, "whitelist", "blacklist"),
          DIF_0.8 = dplyr::if_else(HET_DIF <= 0.8, "whitelist", "blacklist"),
          DIF_0.9 = dplyr::if_else(HET_DIF <= 0.9, "whitelist", "blacklist")
        ) %>%
        tidyr::gather(
          data = .,
          key = DIF_THRESHOLD,
          value = DIF_OUTLIERS,
          DIF_0.1:DIF_0.9
        ) %>%
        tidyr::unite(
          data = .,
          col = MAX_DIF_THRESHOLD, MAX_THRESHOLD, DIF_THRESHOLD,
          sep = "_",
          remove = FALSE
        ) %>%
        dplyr::mutate(
          MAX_DIF_OUTLIERS = dplyr::if_else(
            MAX_OUTLIERS == "whitelist" & DIF_OUTLIERS == "whitelist",
            "whitelist",
            "blacklist")
        ) %>%
        dplyr::group_by(LOCUS, MAX_DIF_THRESHOLD, MAX_THRESHOLD, DIF_THRESHOLD) %>%
        dplyr::summarise(
          MAX_DIF_OUTLIERS = length(MAX_DIF_OUTLIERS[MAX_DIF_OUTLIERS == "whitelist"]),
          MAX_OUTLIERS = length(MAX_OUTLIERS[MAX_OUTLIERS == "whitelist"]),
          DIF_OUTLIERS = length(DIF_OUTLIERS[DIF_OUTLIERS == "whitelist"])
        ) %>%
        dplyr::ungroup(.) %>%
        dplyr::mutate(
          MAX_DIF_OUTLIERS = pop.number - MAX_DIF_OUTLIERS,
          MAX_OUTLIERS = pop.number - MAX_OUTLIERS,
          DIF_OUTLIERS = pop.number - DIF_OUTLIERS
        )

      max.threshold <- dplyr::distinct(
        het.helper, LOCUS, MAX_THRESHOLD, MAX_OUTLIERS
      ) %>%
        dplyr::group_by(MAX_THRESHOLD, MAX_OUTLIERS) %>%
        dplyr::tally(.) %>%
        dplyr::ungroup(.) %>%
        tidyr::complete(
          data = .,
          MAX_OUTLIERS,
          tidyr::nesting(MAX_THRESHOLD),
          fill = list(n = 0)
        ) %>%
        dplyr::group_by(MAX_THRESHOLD) %>%
        dplyr::mutate(
          WHITELIST = cumsum(n),
          BLACKLIST = n.markers - WHITELIST,
          MARKERS = stringi::stri_join(WHITELIST, " (", BLACKLIST, ")")
        ) %>%
        dplyr::select(-c(n, WHITELIST, BLACKLIST)) %>%
        dplyr::group_by(MAX_THRESHOLD, MAX_OUTLIERS) %>%
        tidyr::spread(data = ., key = MAX_THRESHOLD, value = MARKERS, fill = 0)

      dif.threshold <- dplyr::distinct(
        het.helper, LOCUS, DIF_THRESHOLD, DIF_OUTLIERS
      ) %>%
        dplyr::group_by(DIF_THRESHOLD, DIF_OUTLIERS) %>%
        dplyr::tally(.) %>%
        dplyr::ungroup(.) %>%
        tidyr::complete(
          data = .,
          DIF_OUTLIERS,
          tidyr::nesting(DIF_THRESHOLD),
          fill = list(n = 0)
        ) %>%
        dplyr::group_by(DIF_THRESHOLD) %>%
        dplyr::mutate(
          WHITELIST = cumsum(n),
          BLACKLIST = n.markers - WHITELIST,
          MARKERS = stringi::stri_join(WHITELIST, " (", BLACKLIST, ")")
        ) %>%
        dplyr::select(-c(n, WHITELIST, BLACKLIST)) %>%
        dplyr::group_by(DIF_THRESHOLD, DIF_OUTLIERS) %>%
        tidyr::spread(data = ., key = DIF_THRESHOLD, value = MARKERS, fill = 0)


      max.dif.threshold.combined <- dplyr::distinct(
        het.helper, LOCUS, MAX_DIF_THRESHOLD, MAX_DIF_OUTLIERS
      ) %>%
        dplyr::group_by(MAX_DIF_THRESHOLD, MAX_DIF_OUTLIERS) %>%
        dplyr::tally(.) %>%
        dplyr::ungroup(.) %>%
        tidyr::complete(
          data = .,
          MAX_DIF_OUTLIERS,
          tidyr::nesting(MAX_DIF_THRESHOLD),
          fill = list(n = 0)
        ) %>%
        dplyr::group_by(MAX_DIF_THRESHOLD) %>%
        dplyr::mutate(
          WHITELIST = cumsum(n),
          BLACKLIST = n.markers - WHITELIST,
          MARKERS = stringi::stri_join(WHITELIST, " (", BLACKLIST, ")")
        ) %>%
        dplyr::select(-c(n, WHITELIST, BLACKLIST)) %>%
        dplyr::group_by(MAX_DIF_THRESHOLD, MAX_DIF_OUTLIERS) %>%
        tidyr::spread(data = ., key = MAX_DIF_THRESHOLD, value = MARKERS, fill = 0)

      # overall
      het.helper.overall <- dplyr::select(het.summary.overall, LOCUS, HET_MAX, HET_DIF) %>%
        dplyr::mutate(
          MAX_0.4 = dplyr::if_else(HET_MAX <= 0.4, "whitelist", "blacklist"),
          MAX_0.5 = dplyr::if_else(HET_MAX <= 0.5, "whitelist", "blacklist"),
          MAX_0.6 = dplyr::if_else(HET_MAX <= 0.6, "whitelist", "blacklist"),
          MAX_0.7 = dplyr::if_else(HET_MAX <= 0.7, "whitelist", "blacklist"),
          MAX_0.8 = dplyr::if_else(HET_MAX <= 0.8, "whitelist", "blacklist"),
          MAX_0.9 = dplyr::if_else(HET_MAX <= 0.9, "whitelist", "blacklist")
        ) %>%
        tidyr::gather(
          data = .,
          key = MAX_THRESHOLD,
          value = MAX_OUTLIERS,
          MAX_0.4:MAX_0.9
        ) %>%
        dplyr::mutate(
          DIF_0.1 = dplyr::if_else(HET_DIF <= 0.1, "whitelist", "blacklist"),
          DIF_0.2 = dplyr::if_else(HET_DIF <= 0.2, "whitelist", "blacklist"),
          DIF_0.3 = dplyr::if_else(HET_DIF <= 0.3, "whitelist", "blacklist"),
          DIF_0.4 = dplyr::if_else(HET_DIF <= 0.4, "whitelist", "blacklist"),
          DIF_0.5 = dplyr::if_else(HET_DIF <= 0.5, "whitelist", "blacklist"),
          DIF_0.6 = dplyr::if_else(HET_DIF <= 0.6, "whitelist", "blacklist"),
          DIF_0.7 = dplyr::if_else(HET_DIF <= 0.7, "whitelist", "blacklist"),
          DIF_0.8 = dplyr::if_else(HET_DIF <= 0.8, "whitelist", "blacklist"),
          DIF_0.9 = dplyr::if_else(HET_DIF <= 0.9, "whitelist", "blacklist")
        ) %>%
        tidyr::gather(
          data = .,
          key = DIF_THRESHOLD,
          value = DIF_OUTLIERS,
          DIF_0.1:DIF_0.9
        ) %>%
        tidyr::unite(
          data = .,
          col = MAX_DIF_THRESHOLD, MAX_THRESHOLD, DIF_THRESHOLD,
          sep = "_",
          remove = FALSE
        ) %>%
        dplyr::mutate(
          MAX_DIF_OUTLIERS = dplyr::if_else(
            MAX_OUTLIERS == "whitelist" & DIF_OUTLIERS == "whitelist",
            "whitelist",
            "blacklist")
        ) %>%
        dplyr::group_by(LOCUS, MAX_DIF_THRESHOLD, MAX_THRESHOLD, DIF_THRESHOLD) %>%
        dplyr::summarise(
          MAX_DIF_OUTLIERS = length(MAX_DIF_OUTLIERS[MAX_DIF_OUTLIERS == "whitelist"]),
          MAX_OUTLIERS = length(MAX_OUTLIERS[MAX_OUTLIERS == "whitelist"]),
          DIF_OUTLIERS = length(DIF_OUTLIERS[DIF_OUTLIERS == "whitelist"])
        ) %>%
        dplyr::ungroup(.) %>%
        dplyr::mutate(
          MAX_DIF_OUTLIERS = 1 - MAX_DIF_OUTLIERS,
          MAX_OUTLIERS = 1 - MAX_OUTLIERS,
          DIF_OUTLIERS = 1 - DIF_OUTLIERS
        )

      max.threshold.overall <- dplyr::distinct(
        het.helper.overall, LOCUS, MAX_THRESHOLD, MAX_OUTLIERS
      ) %>%
        dplyr::group_by(MAX_THRESHOLD, MAX_OUTLIERS) %>%
        dplyr::tally(.) %>%
        dplyr::ungroup(.) %>%
        tidyr::complete(
          data = .,
          MAX_OUTLIERS,
          tidyr::nesting(MAX_THRESHOLD),
          fill = list(n = 0)
        ) %>%
        dplyr::group_by(MAX_THRESHOLD) %>%
        dplyr::mutate(
          WHITELIST = cumsum(n),
          BLACKLIST = n.markers - WHITELIST,
          MARKERS = stringi::stri_join(WHITELIST, " (", BLACKLIST, ")")
        ) %>%
        dplyr::select(-c(n, WHITELIST, BLACKLIST)) %>%
        dplyr::group_by(MAX_THRESHOLD, MAX_OUTLIERS) %>%
        tidyr::spread(data = ., key = MAX_THRESHOLD, value = MARKERS, fill = 0) %>%
        dplyr::filter(MAX_OUTLIERS  == 0) %>%
        dplyr::ungroup(.) %>%
        dplyr::select(-MAX_OUTLIERS)

      dif.threshold.overall <- dplyr::distinct(
        het.helper.overall, LOCUS, DIF_THRESHOLD, DIF_OUTLIERS
      ) %>%
        dplyr::group_by(DIF_THRESHOLD, DIF_OUTLIERS) %>%
        dplyr::tally(.) %>%
        dplyr::ungroup(.) %>%
        tidyr::complete(
          data = .,
          DIF_OUTLIERS,
          tidyr::nesting(DIF_THRESHOLD),
          fill = list(n = 0)
        ) %>%
        dplyr::group_by(DIF_THRESHOLD) %>%
        dplyr::mutate(
          WHITELIST = cumsum(n),
          BLACKLIST = n.markers - WHITELIST,
          MARKERS = stringi::stri_join(WHITELIST, " (", BLACKLIST, ")")
        ) %>%
        dplyr::select(-c(n, WHITELIST, BLACKLIST)) %>%
        dplyr::group_by(DIF_THRESHOLD, DIF_OUTLIERS) %>%
        tidyr::spread(data = ., key = DIF_THRESHOLD, value = MARKERS, fill = 0) %>%
        dplyr::filter(DIF_OUTLIERS  == 0) %>%
        dplyr::ungroup(.) %>%
        dplyr::select(-DIF_OUTLIERS)


      max.dif.threshold.combined.overall <- dplyr::distinct(
        het.helper.overall, LOCUS, MAX_DIF_THRESHOLD, MAX_DIF_OUTLIERS
      ) %>%
        dplyr::group_by(MAX_DIF_THRESHOLD, MAX_DIF_OUTLIERS) %>%
        dplyr::tally(.) %>%
        dplyr::ungroup(.) %>%
        tidyr::complete(
          data = .,
          MAX_DIF_OUTLIERS,
          tidyr::nesting(MAX_DIF_THRESHOLD),
          fill = list(n = 0)
        ) %>%
        dplyr::group_by(MAX_DIF_THRESHOLD) %>%
        dplyr::mutate(
          WHITELIST = cumsum(n),
          BLACKLIST = n.markers - WHITELIST,
          MARKERS = stringi::stri_join(WHITELIST, " (", BLACKLIST, ")")
        ) %>%
        dplyr::select(-c(n, WHITELIST, BLACKLIST)) %>%
        dplyr::group_by(MAX_DIF_THRESHOLD, MAX_DIF_OUTLIERS) %>%
        tidyr::spread(data = ., key = MAX_DIF_THRESHOLD, value = MARKERS, fill = 0) %>%
        dplyr::filter(MAX_DIF_OUTLIERS  == 0) %>%
        dplyr::ungroup(.) %>%
        dplyr::select(-MAX_DIF_OUTLIERS)

      helper.table.het <- list(
        helper.table.het.max.threshold = max.threshold,
        helper.table.het.max.threshold.overall = max.threshold.overall,
        helper.table.het.dif.threshold = dif.threshold,
        helper.table.het.dif.threshold.overall = dif.threshold.overall,
        helper.table.het.max.dif.threshold.combined = max.dif.threshold.combined,
        helper.table.het.max.dif.threshold.combined.overall = max.dif.threshold.combined.overall

      )
      readr::write_tsv(
        x = max.threshold,
        path = stringi::stri_join(path.folder, "/helper.table.het.max.threshold.tsv")
      )
      readr::write_tsv(
        x = max.threshold.overall,
        path = stringi::stri_join(path.folder, "/helper.table.het.max.threshold.overall.tsv")
      )

      readr::write_tsv(
        x = dif.threshold,
        path = stringi::stri_join(path.folder, "/helper.table.het.dif.threshold.tsv")
      )
      readr::write_tsv(
        x = dif.threshold.overall,
        path = stringi::stri_join(path.folder, "/helper.table.het.dif.threshold.overall.tsv")
      )
      readr::write_tsv(
        x = max.dif.threshold.combined,
        path = stringi::stri_join(path.folder, "/helper.table.het.max.dif.threshold.tsv")
      )
      readr::write_tsv(
        x = max.dif.threshold.combined.overall,
        path = stringi::stri_join(path.folder, "/helper.table.het.max.dif.threshold.overall.tsv")
      )
      n.markers <- het.helper <- max.threshold <- max.threshold.overall <- dif.threshold <- dif.threshold.overall <- max.dif.threshold.combined <- max.dif.threshold.combined.overall <- NULL
    } else {# by SNP
      n.markers <- dplyr::n_distinct(het.summary$MARKERS)

      # by pop
      helper.table.het.pop <- dplyr::select(het.summary, MARKERS, POP_ID, HET_MEAN) %>%
        dplyr::filter(POP_ID != "OVERALL") %>%
        dplyr::mutate(
          `0.4` = dplyr::if_else(HET_MEAN <= 0.4, "whitelist", "blacklist"),
          `0.5` = dplyr::if_else(HET_MEAN <= 0.5, "whitelist", "blacklist"),
          `0.6` = dplyr::if_else(HET_MEAN <= 0.6, "whitelist", "blacklist"),
          `0.7` = dplyr::if_else(HET_MEAN <= 0.7, "whitelist", "blacklist"),
          `0.8` = dplyr::if_else(HET_MEAN <= 0.8, "whitelist", "blacklist"),
          `0.9` = dplyr::if_else(HET_MEAN <= 0.9, "whitelist", "blacklist")
        ) %>%
        tidyr::gather(data = ., key = THRESHOLD, value = OUTLIERS, `0.4`:`0.9`) %>%
        dplyr::group_by(MARKERS, THRESHOLD, THRESHOLD) %>%
        dplyr::summarise(OUTLIERS = length(OUTLIERS[OUTLIERS == "whitelist"])) %>%
        dplyr::ungroup(.) %>%
        dplyr::mutate(OUTLIERS = pop.number - OUTLIERS) %>%
        dplyr::group_by(THRESHOLD, OUTLIERS) %>%
        dplyr::tally(.) %>%
        dplyr::ungroup(.) %>%
        tidyr::complete(data = ., OUTLIERS, tidyr::nesting(THRESHOLD), fill = list(n = 0)) %>%
        dplyr::group_by(THRESHOLD) %>%
        dplyr::mutate(
          WHITELIST = cumsum(n),
          BLACKLIST = n.markers - WHITELIST,
          MARKERS = stringi::stri_join(WHITELIST, " (", BLACKLIST, ")")
        ) %>%
        dplyr::select(-c(n, WHITELIST, BLACKLIST)) %>%
        dplyr::group_by(THRESHOLD, OUTLIERS) %>%
        tidyr::spread(data = ., key = THRESHOLD, value = MARKERS, fill = 0)

      # by overall
      helper.table.het.overall <- dplyr::select(het.summary.overall, MARKERS, HET_MEAN) %>%
        dplyr::mutate(
          `0.4` = dplyr::if_else(HET_MEAN <= 0.4, "whitelist", "blacklist"),
          `0.5` = dplyr::if_else(HET_MEAN <= 0.5, "whitelist", "blacklist"),
          `0.6` = dplyr::if_else(HET_MEAN <= 0.6, "whitelist", "blacklist"),
          `0.7` = dplyr::if_else(HET_MEAN <= 0.7, "whitelist", "blacklist"),
          `0.8` = dplyr::if_else(HET_MEAN <= 0.8, "whitelist", "blacklist"),
          `0.9` = dplyr::if_else(HET_MEAN <= 0.9, "whitelist", "blacklist")
        ) %>%
        tidyr::gather(data = ., key = THRESHOLD, value = OUTLIERS, `0.4`:`0.9`) %>%
        dplyr::group_by(MARKERS, THRESHOLD, THRESHOLD) %>%
        dplyr::summarise(OUTLIERS = length(OUTLIERS[OUTLIERS == "whitelist"])) %>%
        dplyr::ungroup(.) %>%
        dplyr::mutate(OUTLIERS = pop.number - OUTLIERS) %>%
        dplyr::group_by(THRESHOLD, OUTLIERS) %>%
        dplyr::tally(.) %>%
        dplyr::ungroup(.) %>%
        tidyr::complete(data = ., OUTLIERS, tidyr::nesting(THRESHOLD), fill = list(n = 0)) %>%
        dplyr::group_by(THRESHOLD) %>%
        dplyr::mutate(
          WHITELIST = cumsum(n),
          BLACKLIST = 1 - WHITELIST,
          MARKERS = stringi::stri_join(WHITELIST, " (", BLACKLIST, ")")
        ) %>%
        dplyr::select(-c(n, WHITELIST, BLACKLIST)) %>%
        dplyr::group_by(THRESHOLD, OUTLIERS) %>%
        tidyr::spread(data = ., key = THRESHOLD, value = MARKERS, fill = 0) %>%
        dplyr::filter(OUTLIERS  == 0) %>%
        dplyr::ungroup(.) %>%
        dplyr::select(-OUTLIERS)

      helper.table.het <- list(
        helper.table.het.pop = helper.table.het.pop,
        helper.table.het.overall = helper.table.het.overall
      )

      readr::write_tsv(
        x = helper.table.het.pop,
        path = stringi::stri_join(path.folder, "/helper.table.het.threshold.pop.tsv")
      )
      readr::write_tsv(
        x = helper.table.het.overall,
        path = stringi::stri_join(path.folder, "/helper.table.het.threshold.overall.tsv")
      )
      n.markers <- NULL
    }# End by SNP approach
  } else {
    helper.table.het <- "not selected"
  }# End helper table

  message("Generating ", het.approach[1], " heterozygosity statistics plots")
  if (interactive.filter) {
    print(markers.pop.heterozygosity.density.plot)
    message("\n\nGenerating the density plot of heterozygosity statistics: ")
    message("    HET_MEAN: mean heterozygosity")
    message("    HET_MIN: min heterozygosity")
    message("    HET_MAX: max heterozygosity")
    message("    HET_DIF: difference between MAX and MIN heterozygosity of SNPs on the same LOCUS")
  }
  # save
  ggplot2::ggsave(filename = stringi::stri_join(path.folder, "/markers.pop.heterozygosity.density.plot.pdf"), plot = markers.pop.heterozygosity.density.plot, width = pop.number * 4, height = pop.number * 2, dpi = 600, units = "cm", useDingbats = FALSE)
  ggplot2::ggsave(filename = stringi::stri_join(path.folder, "/markers.pop.heterozygosity.density.plot.png"), plot = markers.pop.heterozygosity.density.plot, width = pop.number * 4, height = pop.number * 2, dpi = 300, units = "cm")
  message("2 versions (pdf and png) of the plot (markers.pop.heterozygosity.density.plot) were written in the folder")


  if (interactive.filter) {
    message("\n\nNext, manhattan plot of markers mean observed heterozygosity")
    message("   note: rendering the plot may take some time...")
    print(markers.pop.heterozygosity.manhattan.plot)
  }

  # save
  ggplot2::ggsave(filename = stringi::stri_join(path.folder, "/markers.pop.heterozygosity.manhattan.plot.pdf"), plot = markers.pop.heterozygosity.manhattan.plot, width = pop.number * 4, height = 10, dpi = 600, units = "cm", useDingbats = FALSE)
  ggplot2::ggsave(filename = stringi::stri_join(path.folder, "/markers.pop.heterozygosity.manhattan.plot.png"), plot = markers.pop.heterozygosity.manhattan.plot, width = pop.number * 4, height = 10, dpi = 300, units = "cm")
  message("2 versions (pdf and png) of the plot (markers.pop.heterozygosity.manhattan.plot) were saved in the folder")

  if (interactive.filter) {
    message("\n\nNext plot is a boxplot of markers observed heterozygosity statistics")
    message("   note: rendering the plot may take some time...")
    print(markers.pop.heterozygosity.boxplot)
  }
  # save
  ggplot2::ggsave(filename = stringi::stri_join(path.folder, "/markers.pop.heterozygosity.boxplot.pdf"), plot = markers.pop.heterozygosity.boxplot, width = pop.number * 2, height = 10, dpi = 600, units = "cm", useDingbats = FALSE)
  ggplot2::ggsave(filename = stringi::stri_join(path.folder, "/markers.pop.heterozygosity.boxplot.png"), plot = markers.pop.heterozygosity.boxplot, width = pop.number * 2, height = 10, dpi = 300, units = "cm")
  message("2 versions (pdf and png) of the plot (markers.pop.heterozygosity.boxplot) were saved in the folder")


  if (tibble::has_name(input, "LOCUS") && tibble::has_name(input, "POS")) {
    if (interactive.filter) {
      message("\nThe next plot shows the relationship between the number of SNPs per locus
    and the observed heterozygosity statistics.")
      message("\nNote: more SNPs on the same locus usually = to higher heterozygosity for the locus")
    }
    snp.per.locus.het <- dplyr::distinct(input, LOCUS, POS) %>%
      dplyr::group_by(LOCUS) %>%
      dplyr::tally(.) %>%
      dplyr::rename(SNP_NUM = n) %>%
      dplyr::left_join(het.summary.overall, by = "LOCUS") %>%
      dplyr::select(-POP_ID) %>%
      tidyr::gather(
        data = .,
        key = HET_GROUP,
        value = VALUE,
        -c(LOCUS, SNP_NUM)
      ) %>%
      dplyr::mutate(
        HET_GROUP = factor(
          HET_GROUP,
          levels = c("HET_MEAN", "HET_MIN", "HET_MAX", "HET_DIF"),
          ordered = TRUE))

    snp.per.locus.het.plot <- ggplot2::ggplot(snp.per.locus.het, ggplot2::aes(y = VALUE, x = SNP_NUM)) +
      ggplot2::geom_point() +
      ggplot2::stat_smooth(method = stats::lm, level = 0.99, fullrange = FALSE) +
      # labs(title = "Correlation between missingness and inbreeding coefficient") +
      ggplot2::labs(x = "Number of SNPs per Locus") +
      ggplot2::labs(y = "Observed Heterozygosity Statistics") +
      ggplot2::theme(
        axis.title.x = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
        axis.title.y = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
        legend.title = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
        legend.text = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
        strip.text.x = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold")
      ) +
      ggplot2::facet_grid(~HET_GROUP)
    # snp.per.locus.het.plot

    if (interactive.filter) print(snp.per.locus.het.plot)
    ggplot2::ggsave(filename = stringi::stri_join(path.folder, "/snp.per.locus.het.plot.pdf"), plot = snp.per.locus.het.plot, width = 30, height = 15, dpi = 600, units = "cm", useDingbats = FALSE)
    ggplot2::ggsave(filename = stringi::stri_join(path.folder, "/snp.per.locus.het.plot.png"), plot = snp.per.locus.het.plot, width = 30, height = 15, dpi = 300, units = "cm")
    message("2 versions (pdf and png) of the plot (snp.per.locus.het.plot) were saved in the folder")
  } else {
    snp.per.locus.het.plot <- "not available for the data"
  }

  if (coverage.info && tibble::has_name(input, "READ_DEPTH")) {
    if (interactive.filter) {
      message("\nThe next plot shows the relationship between locus coverage and
    the observed heterozygosity statistics.")
    }
    het.cov <- dplyr::group_by(input, LOCUS) %>%
      dplyr::summarise(COVERAGE = mean(READ_DEPTH, na.rm = TRUE)) %>%
      dplyr::left_join(het.summary.overall, by = "LOCUS") %>%
      dplyr::select(-POP_ID) %>%
      tidyr::gather(
        data = .,
        key = HET_GROUP,
        value = VALUE,
        -c(LOCUS, COVERAGE)
      ) %>%
      dplyr::mutate(
        HET_GROUP = factor(
          HET_GROUP,
          levels = c("HET_MEAN", "HET_MIN", "HET_MAX", "HET_DIF"),
          ordered = TRUE))

    het.cov.plot <- ggplot2::ggplot(het.cov, ggplot2::aes(y = VALUE, x = COVERAGE)) +
      ggplot2::geom_point() +
      ggplot2::stat_smooth(method = stats::lm, level = 0.99, fullrange = FALSE) +
      # labs(title = "Correlation between missingness and inbreeding coefficient") +
      ggplot2::labs(x = "Locus coverage (read depth)") +
      ggplot2::labs(y = "Observed Heterozygosity Statistics") +
      ggplot2::theme(
        axis.title.x = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
        axis.title.y = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
        legend.title = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
        legend.text = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
        strip.text.x = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold")
      ) +
      ggplot2::facet_grid(~HET_GROUP)
    # het.cov.plot

    if (interactive.filter) print(het.cov.plot)
    ggplot2::ggsave(filename = stringi::stri_join(path.folder, "/het.cov.plot.plot.pdf"), plot = het.cov.plot, width = 30, height = 15, dpi = 600, units = "cm", useDingbats = FALSE)
    ggplot2::ggsave(filename = stringi::stri_join(path.folder, "/het.cov.plot.png"), plot = het.cov.plot, width = 30, height = 15, dpi = 300, units = "cm")
    message("2 versions (pdf and png) of the plot (het.cov.plot) were saved in the folder")
  } else {
    het.cov.plot <- "not available for the data"
  }

  # Step 4: Blacklist markers based on observed heterozygosity------------------
  if (interactive.filter) {
    message("\n\nStep 4: Blacklist markers based on observed heterozygosity")
  }

  if (het.approach[1] == "haplotype") {
    message("Approach selected: haplotype")
    # het.threshold
    if (interactive.filter) {
      het.threshold <- 2
      while (isTRUE(het.threshold > 1)) {
        message("het.threshold argument:\n
With the haplotype approach, the het.threshold argument is independent of the
number of SNP/locus and is for filtering the mean maximum observed heterozygosity
of the haplotype. Anything higher (>) than the threshold will be blacklisted.\n")
        message("Example 1: LOCUS 9865 as 2 SNP on the 100 bp read, snp1 HET_OBS = 0.2 and snp2 HET_OBS = 0.9.
The max observed heterozygosity is 0.9, consequently, if you choose a threshold
of 0.5, all SNPs on this locus will be blacklisted, including snp1.\n")
        message("Example 2: LOCUS 2377 as 1 SNP on the 100 bp read, snp1 HET_OBS = 0.6.
The max observed heterozygosity is 0.6, consequently, if you choose a threshold
of 0.5, the SNP/locus/haplotype are blacklisted.\n")

        message("Markers will be discarded from the dataset if they don't pass the 2 next filters
(het.dif.threshold and het.pop.threshold).
Enter the het.threshold value (0 to 1, where 0.1 is very strict and 1 turns off the filter):")
        het.threshold <- as.double(readLines(n = 1))
      }
    }
    # het.dif.threshold
    if (interactive.filter) {
      het.dif.threshold <- 2
      while (isTRUE(het.dif.threshold > 1)) {
        message("het.dif.threshold argument:\n
Markers with 1 SNP/read are not affected by this threshold.
Locus/haplotypes with > 1 SNP are inspected for consistencies of observed heterozygosity along the read.
Locus (read/haplotypes and all SNP on it) are blacklisted if the difference in observed het > threshold.\n")
        message("Example: LOCUS 7654 as 2 SNP, snp1 HET_OBS = 0.1 and snp2 HET_OBS = 0.5.
If you choose a threshold of 0.2 all SNP on this locus will be blacklisted.
Choosing a threshold of 0.6 will not blacklist the SNPs and it's locus.
The locus and it's SNPs will be discarded if they don't pass the remaining filter (het.pop.threshold, with het.approach by pop is selected)\n")
        message("Enter the het.dif.threshold value (0 to 1, where 0.1 is very strict and 1 turn off the filter):")
        het.dif.threshold <- as.double(readLines(n = 1))
      }
    }
    #outlier.pop.threshold
    if (interactive.filter) {
      if (het.approach[2] == "pop") {
        outlier.pop.threshold <- pop.number + 10
        while (isTRUE(outlier.pop.threshold > pop.number)) {
          message("outlier.pop.threshold argument:\n
This filter works by counting the number of populations that didn't pass
the previous 2 filters (het.threshold, het.dif.threshold).
If the number of populations is higher (>) than the threshold, the SNPs and it's locus are discarded.\n")
          message("Useful to incorporate problematic populations dragging down
polymorphism discovery, but still wanted for analysis.
Use this threshold to allow variance in the number of populations passing
the previous thresholds. Blacklist and whitelist produced by the filter
allows to manage outlier markers, individuals and populations.\n")
          message("Example: with a outlier.pop.threshold = 2, you tolerate a
maximum of 2 outlier populations (failing het.threshold and/or het.dif.threshold).
The lower the number, the more severe the filtering, while entering the
number of populations in the dataset turns off the filter.\n")
          message("Enter the outlier.pop.threshold:")

          outlier.pop.threshold <- as.numeric(readLines(n = 1))
        }
      }
    }

    # het.threshold <- 0.5 # test
    # het.dif.threshold <- 0.5# test
    # outlier.pop.threshold <- 2# test
    if (het.approach[2] == "pop") {
      filter <- dplyr::ungroup(het.summary) %>%
        dplyr::group_by(LOCUS) %>%
        dplyr::summarise(OUTLIERS = length(POP_ID[HET_DIF > het.dif.threshold | HET_MAX > het.threshold])) %>%
        dplyr::filter(OUTLIERS <= outlier.pop.threshold) %>%
        dplyr::select(LOCUS) %>%
        dplyr::left_join(input, by = "LOCUS") %>%
        dplyr::arrange(LOCUS, POP_ID) %>%
        dplyr::ungroup(.)
    } else {
      filter <- dplyr::ungroup(het.summary.overall) %>%
        dplyr::group_by(LOCUS) %>%
        dplyr::filter(HET_DIF <= het.dif.threshold) %>%
        dplyr::filter(HET_MAX <= het.threshold) %>%
        dplyr::select(LOCUS) %>%
        dplyr::left_join(input, by = "LOCUS") %>%
        dplyr::arrange(LOCUS, POP_ID) %>%
        dplyr::ungroup(.)
    }
  } else {# by SNP
    message("Approach selected: SNP")
    # het.threshold
    if (interactive.filter) {
      het.threshold <- 2
      while (isTRUE(het.threshold > 1)) {
        message("\nhet.threshold argument:\n
With the SNP approach, the het.threshold argument is independent of the
number of SNP/locus and is for filtering the observed heterozygosity
of each SNP independently of it's haplotype
(there is no het.dif.threshold with this approach, this argument is set to 1 (= off)).\n")
        message("Example 1: LOCUS 9865 as 2 SNP on the 100 bp read, snp1 HET_OBS = 0.2 and snp2 HET_OBS = 0.9.
If you choose a threshold of 0.5, only snp2 is blacklisted.\n")
        message("Markers will be discarded from the dataset if they don't pass
the next filter (het.pop.threshold).
Enter the het.threshold value (0 to 1, where 0.1 is very strict and 1 turns off the filter):")
        het.threshold <- as.double(readLines(n = 1))
      }
    }
    #outlier.pop.threshold
    if (interactive.filter) {
      if (het.approach[2] == "pop") {
        outlier.pop.threshold <- pop.number + 10
        while (isTRUE(outlier.pop.threshold > pop.number)) {
          message("outlier.pop.threshold argument:\n
This filter works by counting the number of populations that didn't pass
the previous filter (het.threshold).
If the number of populations is higher (>) than the threshold, the SNPs and it's locus are discarded.\n")
          message("Useful to incorporate problematic populations dragging down
polymorphism discovery, but still wanted for analysis.
Use this threshold to allow variance in the number of populations passing
the previous thresholds. Blacklist and whitelist produced by the filter
allows to manage outlier markers, individuals and populations.\n")
          message("Example: with a outlier.pop.threshold = 2, you tolerate a
maximum of 2 outlier populations (failing het.threshold and/or het.dif.threshold).
The lower the number, the more severe the filtering, while entering the
number of populations in the dataset turns off the filter.\n")
          message("Enter the outlier.pop.threshold:")
          outlier.pop.threshold <- as.numeric(readLines(n = 1))
        }
      }
    }
    # het.threshold <- 0.5 # test
    # outlier.pop.threshold <- 2# test
    if (het.approach[2] == "pop") {
      filter <- dplyr::ungroup(het.summary) %>%
        dplyr::group_by(MARKERS) %>%
        dplyr::summarise(OUTLIERS = length(POP_ID[HET_MEAN > het.threshold])) %>% #using HET_MEAN is the same with snp approach
        dplyr::filter(OUTLIERS <= outlier.pop.threshold) %>%
        dplyr::select(MARKERS) %>%
        dplyr::left_join(input, by = "MARKERS") %>%
        dplyr::arrange(MARKERS, POP_ID) %>%
        dplyr::ungroup(.)
    } else {
      filter <- dplyr::ungroup(het.summary.overall) %>%
        dplyr::group_by(MARKERS) %>%
        dplyr::filter(HET_MEAN <= het.threshold) %>%
        dplyr::select(MARKERS) %>%
        dplyr::left_join(input, by = "MARKERS") %>%
        dplyr::arrange(MARKERS, POP_ID) %>%
        dplyr::ungroup(.)
    }
  }# end snp approach to filtering

  # Update filters.parameters SNP ----------------------------------------------
  # Prepare a list of markers and number of markers before filtering

  if (data.type == "haplo.file") {
    snp.before <- as.character("NA")
    snp.after <- as.character("NA")
    snp.blacklist <- as.character("NA")
  } else {
    snp.before <- as.integer(dplyr::n_distinct(input$MARKERS))
    snp.after <- as.integer(dplyr::n_distinct(filter$MARKERS))
    snp.blacklist <- as.integer(dplyr::n_distinct(input$MARKERS) - dplyr::n_distinct(filter$MARKERS))
  }

  if (tibble::has_name(input, "LOCUS")) {
    locus.before <- as.integer(dplyr::n_distinct(input$LOCUS))
    locus.after <- as.integer(dplyr::n_distinct(filter$LOCUS))
    locus.blacklist <- as.integer(dplyr::n_distinct(input$LOCUS) - dplyr::n_distinct(filter$LOCUS))
  } else {
    locus.before <- as.character("NA")
    locus.after <- as.character("NA")
    locus.blacklist <- as.character("NA")
  }

  ind.before <- dplyr::n_distinct(input$INDIVIDUALS)
  ind.blacklisted <- length(blacklist.ind.het$INDIVIDUALS)
  ind.after <- ind.before - ind.blacklisted

  markers.before <- stringi::stri_join(snp.before, locus.before, sep = "/")
  markers.after <- stringi::stri_join(snp.after, locus.after, sep = "/")
  markers.blacklist <- stringi::stri_join(snp.blacklist, locus.blacklist, sep = "/")


  if (tibble::has_name(het.summary, "LOCUS")) {
    markers.df <- dplyr::distinct(input, CHROM, LOCUS, POS)
  } else {
    markers.df <- dplyr::distinct(input, MARKERS)
  }


  if (het.approach[2] == "overall") outlier.pop.threshold <- "using overall"

  filters.parameters <- tibble::data_frame(
    FILTERS = c("Observed Heterozygosity", rep(as.character(""), 4)),
    PARAMETERS = c("ind.mean.het", "het.approach", "het.threshold", "het.dif.threshold", "outlier.pop.threshold"),
    VALUES = c(paste(ind.heterozygosity.threshold, collapse = "/"), paste(het.approach, collapse = " and "), paste("<=", het.threshold), paste("<=", het.dif.threshold), paste("<=", outlier.pop.threshold)),
    BEFORE = c(ind.before, "", "", "", markers.before),
    AFTER = c(ind.after, "", "", "", markers.after),
    BLACKLIST = c(ind.blacklisted, "", "", "", markers.blacklist),
    UNITS = c("individual" , "", "", "", "SNP/LOCUS"),
    COMMENTS = c("min/max values", "", "", "", "")
  )
  readr::write_tsv(x = filters.parameters, path = "filters_parameters.tsv", append = TRUE, col_names = FALSE)

  # saving filtered tidy data --------------------------------------------------
  if (!is.null(filename)) {
    tidy.name <- stringi::stri_join(filename, ".rad")
    message("Writing the filtered tidy data set: ", tidy.name)
    fst::write.fst(x = filter, path = stringi::stri_join(path.folder, "/", tidy.name), compress = 85)
  }

  # saving whitelist -----------------------------------------------------------
  message("Writing the whitelist of markers: whitelist.markers.pop.tsv")

  if (tibble::has_name(het.summary, "LOCUS")) {
    whitelist.markers <- dplyr::ungroup(filter) %>%
      dplyr::distinct(MARKERS, CHROM, LOCUS, POS)
  } else {
    whitelist.markers <- dplyr::ungroup(filter) %>%
      dplyr::distinct(MARKERS)
  }

  # Integrate marker.meta columns (double check from above)
  # if (!is.null(marker.meta)) {
  #   whitelist.markers <- dplyr::left_join(whitelist.markers, marker.meta, by = "MARKERS")
  # }

  readr::write_tsv(whitelist.markers, stringi::stri_join(path.folder, "/", "whitelist.markers.het.tsv"), append = FALSE, col_names = TRUE)


  # saving blacklist -----------------------------------------------------------
  message("Writing the blacklist of markers: blacklist.markers.pop.tsv")
  if (tibble::has_name(het.summary, "LOCUS")) {
    blacklist.markers <- dplyr::anti_join(markers.df, whitelist.markers, by = c("CHROM", "LOCUS", "POS"))
  } else {
    blacklist.markers <- dplyr::anti_join(markers.df, whitelist.markers, by = "MARKERS")
  }

  # Integrate marker.meta columns (double check from above)
  # if (!is.null(marker.meta)) {
  #   blacklist.markers <- dplyr::left_join(blacklist.markers, marker.meta, by = "MARKERS")
  # }
  readr::write_tsv(blacklist.markers, paste0(path.folder,"/blacklist.markers.het.tsv"), append = FALSE, col_names = TRUE)

  # results --------------------------------------------------------------------
  cat("############################### RESULTS ###############################\n")
  message("ind.heterozygosity.threshold (min/max): ", paste(ind.heterozygosity.threshold, collapse = "/"))
  message("Blacklisted individuals: ", ind.blacklisted)
  message("het.approach: ", paste(het.approach, collapse = " and "))
  message("het.threshold: ", het.threshold)
  message("het.dif.threshold: ", het.dif.threshold)
  message("outlier.pop.threshold: ", outlier.pop.threshold)
  message("The number of markers (SNP/LOCUS) removed by the HET filter:\n", markers.blacklist)
  message("The number of markers (SNP/LOCUS) before -> after the HET filter:\n", markers.before, " -> ", markers.after)
  if (!interactive.filter) {
    message("Computation time: ", round((proc.time() - timing)[[3]]), " sec")
  }
  cat("############################## completed ##############################\n")
  res <- list(
    tidy.filtered.het = filter,
    whitelist.markers = whitelist.markers,
    blacklist.markers = blacklist.markers,
    strata = strata.df,
    filters.parameters = filters.parameters,
    individual.heterozygosity = het.ind,
    heterozygosity.statistics = het.ind.stats,
    blacklist.ind.het = blacklist.ind.het,
    helper.table.het = helper.table.het,
    individual.heterozygosity.boxplot = individual.heterozygosity.boxplot,
    individual.heterozygosity.manhattan.plot = individual.heterozygosity.manhattan.plot,
    markers.pop.heterozygosity.boxplot = markers.pop.heterozygosity.boxplot,
    markers.pop.heterozygosity.density.plot = markers.pop.heterozygosity.density.plot,
    markers.pop.heterozygosity.manhattan.plot = markers.pop.heterozygosity.manhattan.plot,
    snp.per.locus.het.plot = snp.per.locus.het.plot,
    het.cov.plot = het.cov.plot
  )
  return(res)
}

# param ind.heterozygosity.threshold (string, double, optional)
# Blacklist individuals based on observed heterozygosity (averaged across markers).
#
#
# The string contains 2 thresholds values (min and max).
# The values are proportions (0 to 1), where 0 turns off the min threshold and
# 1 turns off the max threshold.
# Individuals with mean observed heterozygosity higher (>) or lower (<)
# than the thresholds will be blacklisted.
#
# Default: \code{ind.heterozygosity.threshold = NULL} will turn off completely
# the filter and the function will only output the plots and table of heterozygosity.



#  \strong{Individual observed heterozygosity (averaged across markers):}
# To help discard an individual based on his observed heterozygosity
# (averaged across markers),
# use the manhanttan plot to:
# \enumerate{
# \item contrast the individual with population and overall samples.
# \item visualize the impact of missigness information (based on population or
# overall number of markers) and the individual observed heterozygosity. The
# larger the point, the more missing genotypes.
# }
# \strong{Outlier above average:}
# \itemize{
# \item potentially represent two samples mixed together (action: blacklist), or...
# \item a sample with more sequecing effort (point size small): did you merge your replicates fq files ? (action : keep and monitor)
# \item a sample with poor sequencing effort (point size large) where the genotyped markers are
# all heterozygotes, verify this with missingness (action: discard)
# }
# In all cases, if there is no bias in the population sequencing effort,
# the size of the point will usually be "average" based on the population or
# overall number of markers.




#
# \strong{Outlier below average:}
# \itemize{
# \item A point with a size larger than the population or overall average (= lots of missing):
# the poor polymorphism discovery of the sample is probably the result of bad
# DNA quality, a bias in sequencing effort, etc. (action: blacklist)
# \item A point with a size that looks average (not much missing):
# this sample requires more attention (action: blacklist) and conduct more tests.
# e.g. for biallelic data, look for coverage imbalance between ALT/REF allele.
# At this point you need to distinguish between an artifact of poor polymorphism discovery
# or a biological reason (highly inbred individual, etc.).
# }
#

