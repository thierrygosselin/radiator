#### Visualization



# count ind_total_reads -------------------------------------------------------
#' @title ind_total_reads
#' @description Counts the total number of reads per samples
#' @rdname ind_total_reads
#' @keywords internal
#' @export
ind_total_reads <- function(x, path.folder = NULL) {
  # x <- unfiltered
  # path.folder <- path.folder.coverage
  if (is.null(path.folder)) path.folder <- getwd()
  # generate the stats
  read.info <- dplyr::group_by(x, INDIVIDUALS, POP_ID) %>%
    dplyr::summarise(TOTAL_READ_COUNTS = sum(READ_DEPTH, na.rm = TRUE))

  if (is.factor(read.info$POP_ID)) {
    read.pop.levels <- levels(read.info$POP_ID)
  } else {
    read.pop.levels <- unique(read.info$POP_ID)
  }
  read.pop.levels <- c(read.pop.levels, "OVERALL")

  overall.read <- dplyr::mutate(read.info, POP_ID = "OVERALL")


  read.info <- suppressWarnings(dplyr::bind_rows(read.info, overall.read)) %>%
    dplyr::mutate(POP_ID = factor(POP_ID, levels = read.pop.levels))

  read.counts.stats <- read.info %>%
    dplyr::group_by(POP_ID) %>%
    dplyr::summarise(
      MEAN = mean(TOTAL_READ_COUNTS, na.rm = TRUE),
      MEDIAN = stats::median(TOTAL_READ_COUNTS, na.rm = TRUE),
      SD = stats::sd(TOTAL_READ_COUNTS, na.rm = TRUE),
      Q25 = stats::quantile(TOTAL_READ_COUNTS, 0.25, na.rm = TRUE),
      Q75 = stats::quantile(TOTAL_READ_COUNTS, 0.75, na.rm = TRUE),
      IQR = stats::IQR(TOTAL_READ_COUNTS, na.rm = TRUE),
      MIN = min(TOTAL_READ_COUNTS, na.rm = TRUE),
      MAX = max(TOTAL_READ_COUNTS, na.rm = TRUE),
      OUTLIERS_LOW = Q25 - (1.5 * IQR),
      OUTLIERS_HIGH = Q75 + (1.5 * IQR),
      OUTLIERS_LOW_N = length(TOTAL_READ_COUNTS[TOTAL_READ_COUNTS < OUTLIERS_LOW]),
      OUTLIERS_HIGH_N = length(TOTAL_READ_COUNTS[TOTAL_READ_COUNTS > OUTLIERS_HIGH]),
      OUTLIERS_TOTAL = OUTLIERS_HIGH_N + OUTLIERS_LOW_N,
      OUTLIERS_PROP = round(OUTLIERS_TOTAL / length(x), 3)
    )

  # plots
  element.text <- ggplot2::element_text(size = 10,
                                        face = "bold")
  n.pop <- dplyr::n_distinct(x$POP_ID)
  ind.plot <- suppressWarnings(
    ggplot2::ggplot(
      read.info, ggplot2::aes(x = POP_ID, y = TOTAL_READ_COUNTS, na.rm = TRUE)) +
      ggplot2::geom_boxplot() +
      ggplot2::labs(x = "Populations",
                    y = "Individuals total read counts",
                    title = "Individuals total read counts") +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
        # plot.subtitle = ggplot2::element_text(size = 10, hjust = 0.5),
        axis.title.y = element.text,
        axis.title.x = element.text,
        axis.text.x = element.text))

  suppressWarnings(
    ggplot2::ggsave(
      plot = ind.plot,
      filename = file.path(path.folder, "plot.ind.total.reads.pdf"),
      width = n.pop * 2, height = 10, dpi = 300, units = "cm",
      useDingbats = FALSE)
  )

}#End ind_coverage

# markers_genotyped_helper------------------------------------------------------
#' @title markers_genotyped_helper
#' @description Help individual's genotyped threshold
#' @rdname markers_genotyped_helper
#' @export
#' @keywords internal
markers_genotyped_helper <- function(x, y, overall.only = FALSE) {
  # x <- res$missing.genotypes.markers.pop
  # Set the breaks for the figure
  max.markers <- dplyr::n_distinct(y$MARKERS)

  threshold.helper.overall <- y %>%
    dplyr::ungroup(.) %>%
    dplyr::summarise(
      `0` = length(PERCENT[PERCENT == 0]),
      `10` = length(PERCENT[PERCENT <= 10]),
      `20` = length(PERCENT[PERCENT <= 20]),
      `30` = length(PERCENT[PERCENT <= 30]),
      `40` = length(PERCENT[PERCENT <= 40]),
      `50` = length(PERCENT[PERCENT <= 50]),
      `60` = length(PERCENT[PERCENT <= 60]),
      `70` = length(PERCENT[PERCENT <= 70]),
      `80` = length(PERCENT[PERCENT <= 80]),
      `90` = length(PERCENT[PERCENT <= 90]),
      `100` = length(PERCENT[PERCENT <= 100])
    ) %>%
    tidyr::pivot_longer(
      data = .,
      cols = tidyselect::everything(),
      names_to = "GENOTYPED_THRESHOLD",
      values_to = "NUMBER_MARKERS"
    ) %>%
    dplyr::mutate(POP_ID = rep("OVERALL", n()))

  if (!overall.only) {
    threshold.helper.pop <- x %>%
      dplyr::group_by(POP_ID) %>%
      dplyr::summarise(
        `0` = length(PERCENT[PERCENT == 0]),
        `10` = length(PERCENT[PERCENT <= 10]),
        `20` = length(PERCENT[PERCENT <= 20]),
        `30` = length(PERCENT[PERCENT <= 30]),
        `40` = length(PERCENT[PERCENT <= 40]),
        `50` = length(PERCENT[PERCENT <= 50]),
        `60` = length(PERCENT[PERCENT <= 60]),
        `70` = length(PERCENT[PERCENT <= 70]),
        `80` = length(PERCENT[PERCENT <= 80]),
        `90` = length(PERCENT[PERCENT <= 90]),
        `100` = length(PERCENT[PERCENT <= 100])
      ) %>%
      tidyr::pivot_longer(
        data = .,
        cols = -POP_ID,
        names_to = "GENOTYPED_THRESHOLD",
        values_to = "NUMBER_MARKERS"
      )

    mean.pop <- threshold.helper.pop %>%
      dplyr::group_by(GENOTYPED_THRESHOLD) %>%
      dplyr::summarise(
        NUMBER_MARKERS = round(mean(NUMBER_MARKERS), 0)
      ) %>%
      dplyr::mutate(POP_ID = rep("MEAN_POP", n()))

    # Check if x$POP_ID is a factor
    if (is.factor(x$POP_ID)) {
      x.pop.levels <- levels(x$POP_ID)
    } else {
      x.pop.levels <- unique(x$POP_ID)
    }

    threshold.helper <- suppressWarnings(
      dplyr::bind_rows(threshold.helper.pop, mean.pop, threshold.helper.overall) %>%
        dplyr::mutate(
          GENOTYPED_THRESHOLD = as.numeric(GENOTYPED_THRESHOLD),
          POP_ID = factor(POP_ID, levels = c(x.pop.levels, "MEAN_POP", "OVERALL"), ordered = TRUE)
        ))
    threshold.helper.pop <- mean.pop <- threshold.helper.overall <- x <- y <- NULL
  } else {
    threshold.helper <- threshold.helper.overall %>%
      dplyr::mutate(
        GENOTYPED_THRESHOLD = as.numeric(GENOTYPED_THRESHOLD)
      )
    threshold.helper.overall <- x <- y <- NULL
  }
  #Function to replace package plyr round_any
  rounder <- function(x, accuracy, f = round) {
    f(x / accuracy) * accuracy
  }

  if (max.markers >= 1000) {
    y.breaks.by <- rounder(max.markers / 10, 100, ceiling)
    y.breaks.max <- rounder(max.markers, 1000, ceiling)
    y.breaks <- seq(0, y.breaks.max, by = y.breaks.by)
  } else {
    y.breaks.by <- rounder(max.markers / 10, 10, ceiling)
    y.breaks.max <- rounder(max.markers, 100, ceiling)
    y.breaks <- seq(0, y.breaks.max, by = y.breaks.by)
  }


  plot.markers.geno.threshold <- ggplot2::ggplot(
    threshold.helper,
    ggplot2::aes(x = GENOTYPED_THRESHOLD, y = NUMBER_MARKERS)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2, shape = 21, fill = "white") +
    ggplot2::scale_x_continuous(name = "Marker's missing genotype threshold (percent)", breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
    ggplot2::scale_y_continuous(name = "Markers\n(whitelisted number)", breaks = y.breaks, limits = c(0, y.breaks.max)) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 10, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 10, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 8),#, angle = 90, hjust = 1, vjust = 0.5),
      strip.text.x = ggplot2::element_text(size = 10, face = "bold")
    ) +
    ggplot2::theme_bw() +
    ggplot2::facet_grid(~POP_ID)
  # plot.markers.geno.threshold
  return(plot.markers.geno.threshold)
}#End markers_genotyped_helper

# tibble_stats-----------------------------------------------------------------
#' @title tibble_stats
#' @description Generate a tibble of statistics
#' @rdname tibble_stats
#' @keywords internal
#' @export
tibble_stats <- function(x, group, subsample = NULL) {
  if (is.null(subsample)) subsample <- 1L
  Q <- stats::quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  res <- tibble::tibble(
    GROUP = group,
    MIN = min(x, na.rm = TRUE),
    # Q25 = stats::quantile(x, 0.25, na.rm = TRUE),
    Q25 = Q[1],
    MEAN = mean(x, na.rm = TRUE),
    MEDIAN = stats::median(x, na.rm = TRUE),
    Q75 = Q[2],
    MAX = max(x, na.rm = TRUE),
    IQR = abs(diff(Q)),
    # IQR = stats::IQR(depth, na.rm = TRUE),
    OUTLIERS_LOW = Q25 - (1.5 * IQR),
    OUTLIERS_HIGH =  Q75 + (1.5 * IQR)) %>%
    dplyr::mutate(dplyr::across(where(is.integer), .fns = as.numeric)) %>%
    dplyr::mutate(
      OUTLIERS_LOW = dplyr::if_else(OUTLIERS_LOW < MIN, MIN, OUTLIERS_LOW),
      OUTLIERS_HIGH = dplyr::if_else(OUTLIERS_HIGH > MAX, MAX, OUTLIERS_HIGH),
      SUBSAMPLE = subsample
    )
  Q <- NULL
  return(res)
}#End tibble_stats

# boxplot_stats-----------------------------------------------------------------
# boxplot of stats
#' @title boxplot_stats
#' @description Generate a boxplot
#' @rdname boxplot_stats
#' @keywords internal
#' @export
boxplot_stats <- function(
  data,
  title,
  subtitle = NULL,
  x.axis.title = NULL,
  y.axis.title,
  facet.columns = FALSE,
  facet.rows = FALSE,
  bp.filename = NULL,
  path.folder = NULL
) {
  # data <- test
  # x.axis.title = NULL
  # x.axis.title <- "SNP position on the read groupings"
  # title <- "Individual's QC stats"
  # subtitle = "testing"
  # y.axis.title <- "Statistics"
  # y.axis.title <- "SNP position (base pair)"
  # bp.filename <- "vcf.snp.position.read.pdf"
  # bp.filename <- "test.pdf"
  # facet.columns = TRUE
  # facet.rows = FALSE
  # path.folder = NULL

  if (is.null(path.folder)) path.folder <- getwd()

  n.group <- dplyr::n_distinct(data$GROUP)
  element.text <- ggplot2::element_text(size = 10,
                                        face = "bold")

  if (facet.columns) {
    data <- dplyr::mutate(data, X = "1")
    fig.boxplot <- ggplot2::ggplot(data = data, ggplot2::aes(X))
  } else {
    fig.boxplot <- ggplot2::ggplot(data = data, ggplot2::aes(GROUP))
  }


  fig.boxplot <- fig.boxplot +
    ggplot2::geom_boxplot(
      ggplot2::aes(ymin = OUTLIERS_LOW, lower = Q25, middle = MEDIAN, upper = Q75,
                   ymax = OUTLIERS_HIGH), stat = "identity") +
    ggplot2::labs(y = y.axis.title, title = title)

  if (!is.null(subtitle)) fig.boxplot <- fig.boxplot + ggplot2::labs(subtitle = subtitle)

  # Draw upper outliers
  if (facet.columns) {
    fig.boxplot <- fig.boxplot +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = "1", xend = "1",
          y = OUTLIERS_HIGH, yend = MAX),
        linetype = "dashed")
  } else {
    fig.boxplot <- fig.boxplot +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = GROUP, xend = GROUP,
          y = OUTLIERS_HIGH, yend = MAX),
        linetype = "dashed")
  }

  # Draw lower outliers
  if (facet.columns) {
    fig.boxplot <- fig.boxplot +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = "1", xend = "1",
          y = OUTLIERS_LOW, yend = MIN),
        linetype = "dashed")
  } else {
    fig.boxplot <- fig.boxplot +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = GROUP, xend = GROUP,
          y = OUTLIERS_LOW, yend = MIN),
        linetype = "dashed")
  }

  fig.boxplot <- fig.boxplot +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
      legend.position = "none",
      axis.title.y = element.text,
      axis.text.y = element.text
    )

  if (is.null(x.axis.title)) {
    fig.boxplot <- fig.boxplot +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  } else {
    fig.boxplot <- fig.boxplot +
      ggplot2::xlab(x.axis.title) +
      ggplot2::theme(
        axis.title.x = element.text,
        axis.ticks.x = ggplot2::element_blank()
      )
  }

  if (!is.null(subtitle)) {
    fig.boxplot <- fig.boxplot +
      ggplot2::theme(
        plot.subtitle = ggplot2::element_text(size = 10))
  }

  if (facet.columns) {
    fig.boxplot <- fig.boxplot + ggplot2::facet_grid(GROUP ~ ., scales = "free")
    n.facet <- n.group * 2
    width <- 15
    height <- 5 + (4 * n.group)
  }

  if (facet.rows) {
    fig.boxplot <- fig.boxplot + ggplot2::facet_grid(FACET_ROWS ~ ., scales = "free")
    n.facet <- n.group * 2
    width <- 10
    height <- 5 + (4 * n.group)
  }

  if (!facet.rows && !facet.columns) {
    width <-  13 + (5 * n.group) + 1
    height <-  8
  }

  print(fig.boxplot)
  if (!is.null(bp.filename)) {
    suppressMessages(
      ggplot2::ggsave(
        filename = file.path(path.folder, bp.filename),
        plot = fig.boxplot,
        width = width,
        height = height,
        dpi = 300,
        units = "cm",
        limitsize = FALSE,
        useDingbats = FALSE
      )
    )
  }
  return(fig.boxplot)
}#Endboxplot_stats




# plot_density_distribution_coverage -------------------------------------------
#' @title Figure density distribution of coverage summary statistics
#' @description Create density distribution of coverage summary statistics.
#' Use the coverage summary file created with coverage_summary function.
#' @param data Coverage summary file.
#' @param aes.colour GGPLOT2 aesthetics colour,
#' e.g. aes(y = ..scaled.., color = COVERAGE_GROUP).
#' @param adjust.bin Adjust GGPLOT2 bin size (0 to 1).
#' @export
#' @rdname plot_density_distribution_coverage
#' @keywords internal


plot_density_distribution_coverage <- function(data, aes.colour, adjust.bin) {

  VALUE <- NULL



  ggplot2::ggplot(data, ggplot2::aes(x = VALUE, na.rm = TRUE)) +
    ggplot2::geom_line(aes.colour, stat = "density", size = 0.5, adjust = adjust.bin) +
    ggplot2::labs(x = "Depth of coverage(read number)") +
    ggplot2::labs(y = "Loci (scaled density)") +
    ggplot2::expand_limits(x = 0) +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_text(size = 10, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 10, face = "bold"),
      legend.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.text = ggplot2::element_text(size = 10, face = "bold"),
      strip.text.x = ggplot2::element_text(size = 10, face = "bold"))
}


#' @title Figure box plot of coverage summary statistics
#' @description Create box plots of coverage summary statistics.
#' Use the coverage summary file created with coverage_summary function.
#' @param data Coverage summary file.
#' @export
#' @rdname plot_boxplot_coverage
#' @keywords internal

plot_boxplot_coverage <- function(data) {

  POP_ID <- NULL
  VALUE <- NULL
  POP_ID <- NULL
  VALUE <- NULL


  ggplot2::ggplot(data, ggplot2::aes(x = factor(POP_ID), y = VALUE, na.rm = TRUE)) +
    ggplot2::geom_violin(trim = FALSE) +
    ggplot2::geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA) +
    ggplot2::stat_summary(fun = "mean", geom = "point", shape = 21, size = 2.5, fill = "white") +
    ggplot2::labs(x = "Sampling sites") +
    ggplot2::labs(y = "Read depth coverage") +
    ggplot2::facet_wrap(facets = ~COVERAGE_GROUP, scales = "free") +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_text(size = 10, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 10, face = "bold"),
      legend.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.text = ggplot2::element_text(size = 10, face = "bold"),
      strip.text.x = ggplot2::element_text(size = 10, face = "bold")
    )
}




#' @title Visual diagnostic of coverage imbalance
#' @description GBS data and STACKS pipeline sometimes output REF and ALT
#' alleles that have coverage imbalance, i.e. the coverage is not equal
#' and skewed towards the REF or ALT alleles.
#' Thw density distribution figure of coverage imbalance between REF and ALT
#' alleles will highlight the problem in you data.
#' @param tidy.vcf.file The tidy VCF file created with tidy_genomic_data.
#' @param pop.levels Character string defining your ordered populations.
#' @param read.depth.threshold Define the threshold you wish to analyse.
#' @param aes.colour GGPLOT2 aesthetics,
#' e.g. aes(y = ..count..).
#' @param adjust.bin Adjust GGPLOT2 bin size (0 to 1).
#' @return 4-plots highlighting the different combination of under
#'  or over the coverage threshold and mean genotype likelihood.
#'  Y- axis show the distribution of genotypes.
#'  X- axis show the coverage imbalance the
#'  Negative ratio (left x axis) : REF > ALT.
#'  Positive ratio (right x axis) : ALT > REF.
#' @details The figures shows the results of the of coverage threshold
#' selected and mean genotype likelihood.
#' You can test different threshold to inspect your data.
#' Ideally the lower left pannel of the 4-plot should be empty. If it is, this
#' shows that setting the threshold of the genotype likelihood filter
#' to the mean or close to it take care of the allelic coverage imbalance.
#' #' e.g. fig <- plot_coverage_imbalance_diagnostic(
#' tidy.vcf.file, pop.levels, read.depth.threshold, aes.colour, adjust.bin)
#' Use ( fig + facet_grid(GROUP_GL ~ GROUP_COVERAGE)). The ratio is calculated :
#' (read depth ALT allele - read depth REF allele)/(read depth ALT allele + read depth REF allele).
#' @export
#' @rdname plot_coverage_imbalance_diagnostic
#' @keywords internal


plot_coverage_imbalance_diagnostic <- function(tidy.vcf.file, pop.levels, read.depth.threshold, aes.colour, adjust.bin) {

  INDIVIDUALS <- NULL
  POP_ID <- NULL
  READ_DEPTH <- NULL
  GL <- NULL
  GROUP_COVERAGE <- NULL
  GROUP_GL <- NULL
  ALLELE_COVERAGE_RATIO <- NULL



  if (is.vector(tidy.vcf.file)) {
    data <- readr::read_tsv(tidy.vcf.file, col_names = T, col_types = "diidccddccccdddddc") %>%
      dplyr::mutate(INDIVIDUALS = factor(INDIVIDUALS))
  } else {
    data <- tidy.vcf.file
  }

  if (missing(pop.levels)) {
    data <- data
  } else {
    data <- suppressWarnings(
      data %>%
        dplyr::mutate(POP_ID = factor(POP_ID, levels = pop.levels, ordered = T))
    )
  }

  below.threshold <- stringi::stri_join("READ DEPTH <=", read.depth.threshold, sep = " ")
  over.threshold <- stringi::stri_join("READ DEPTH >", read.depth.threshold, sep = " ")

  imbalance.coverage <- data %>%
    dplyr::mutate(
      #       GROUP_COVERAGE = ifelse(READ_DEPTH <= read.depth.threshold, "READ DEPTH <= 8", "READ DEPTH > 8"),
      GROUP_COVERAGE = ifelse(READ_DEPTH <= read.depth.threshold, below.threshold, over.threshold),
      GROUP_GL = ifelse(GL < (mean(data$GL, na.rm = T)), "< mean GL", ">= mean GL"),
      GROUP_COVERAGE = factor(GROUP_COVERAGE, levels = c(below.threshold, over.threshold), ordered = T),
      #       GROUP_COVERAGE = factor(GROUP_COVERAGE, levels = c("READ DEPTH <= 8", "READ DEPTH > 8"), ordered = T),
      GROUP_GL = factor(GROUP_GL, levels = c("< mean GL", ">= mean GL"), ordered = T)
    ) %>%
    dplyr::filter(GROUP_COVERAGE != "NA" & GROUP_GL != "NA")

  ggplot2::ggplot(imbalance.coverage, ggplot2::aes(x = ALLELE_COVERAGE_RATIO, na.rm = T)) +
    #   geom_bar() +
    ggplot2::geom_line(aes.colour, stat = "density", adjust = adjust.bin) +
    ggplot2::labs(x = "Coverage imbalance between REF and ALT alleles (ratio)") +
    ggplot2::labs(y = "Distribution of genotypes") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 12, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 12, face = "bold"),
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.text = ggplot2::element_text(size = 12, face = "bold"),
      strip.text.x = ggplot2::element_text(size = 12, face = "bold"),
      strip.text.y = ggplot2::element_text(size = 12, face = "bold")
    )

}










#' @title Figure density distribution of minor allele frequency (MAF)
#' summary statistics.
#' @description Create density distribution of MAF summary statistics.
#' @param data sumstats or tidy vcf summarised files.
#' @param maf.group The GGPLOT2 aes (e.g. aes(x = FREQ_ALT, na.rm = F)).
#' @param aes.colour GGPLOT2 aesthetics colour,
#' e.g. aes(y = ..scaled.., color = POP_ID) or aes(y = ..count.., color = POP_ID)
#' @param adjust.bin Adjust GGPLOT2 bin size (0 to 1).
#' @param x.title Title of the x-axis. e.g. "MAF distribution"
#' @export
#' @rdname plot_density_distribution_maf
#' @details turn off the legend using \code{fig + theme(legend.position = "none")} or
#' zoom in a section with 'fig + coord_cartesian(xlim = c(0, 0.1), ylim = c(0, 1))'.
#' Save the figure with : \code{ggsave("figure name.pdf", width = 40, height = 20, dpi = 600, units = "cm", useDingbats = F)}.
#' @seealso \link{tidy_genomic_data}
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @keywords internal

plot_density_distribution_maf <- function(data, maf.group, aes.colour = ggplot2::aes(y = ..scaled.., color = POP_ID), adjust.bin = 1, x.title) {

  ..scaled.. <- NULL #get rid of R CMD check note

  if (is.vector(data)) {
    data <- readr::read_tsv(data, col_names = T)
  } else {
    data <- data
  }
  #   font.group <-
  graph <- ggplot2::ggplot(data, maf.group) +
    ggplot2::geom_line(aes.colour, stat = "density", adjust = adjust.bin) + # pop colored
    #   scale_colour_manual(name ="Sampling sites", values = colour_palette_sites.pink) +
    ggplot2::scale_x_continuous(breaks = c(0, 0.05, 0.1, 0.2, 0.5, 1),
                                labels = c("0", "0.05", "0.1", "0.2", "0.5", "1.0")) +
    ggplot2::labs(x = x.title) +
    ggplot2::labs(y = "Density of SNP (scaled)") +
    ggplot2::expand_limits(y = 0) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 12, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 12, face = "bold"),
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.text = ggplot2::element_text(size = 12, face = "bold"),
      strip.text.y = ggplot2::element_text(angle = 0, size = 12, face = "bold"),
      strip.text.x = ggplot2::element_text(size = 12, face = "bold")
    )
}





# HET Figure
#' @title Figure density distribution of the observed heterozygosity
#' summary statistics.
#' @description Create density distribution of the observed heterozygosity
#' summary statistics.
#' @param data sumstats or tidy vcf files.
#' @param pop.levels Character string defining your ordered populations.
#' @param het.group = aes(x = HET_MAX, na.rm = F)
#' @param aes.colour GGPLOT2 aesthetics colour,
#' e.g. aes(y = ..scaled.., color = GROUP).
#' @param adjust.bin Adjust GGPLOT2 bin size (0 to 1).
#' @param x.title Title of the x-axis.
#' @export
#' @rdname plot_density_distribution_het
#' @keywords internal

plot_density_distribution_het <- function(data, pop.levels, het.group, aes.colour, adjust.bin, x.title){

  POP_ID <- NULL
  HET_O <- NULL
  HET_MAX <- NULL
  HET_MIN <- NULL
  VALUE <- NULL

  if (is.vector(data)) {
    data <- readr::read_tsv(data, col_names = T)
  } else {
    data = data
  }

  data.summary <- data %>%
    dplyr::group_by(LOCUS, POP_ID) %>%
    dplyr::summarise(
      HET_MEAN = mean(HET_O),
      HET_MAX = max(HET_O),
      HET_MIN = min(HET_O),
      HET_DIFF = HET_MAX - HET_MIN
    ) %>%
    dplyr::group_by(LOCUS, POP_ID) %>%
    tidyr::pivot_longer(
      data = .,
      cols = tidyselect::everything(),
      names_to = "HET_GROUP",
      values_to = "VALUE"
    )

  if (missing(pop.levels) == "TRUE") {
    data.summary <- data.summary
  } else {
    data.summary <- data.summary %>%
      dplyr::mutate(POP_ID = factor(POP_ID, levels = pop.levels, ordered = T)) %>%
      dplyr::arrange(POP_ID)
  }

  graph <- ggplot2::ggplot(data.summary, ggplot2::aes(x = VALUE, na.rm = F)) +
    ggplot2::geom_line(aes.colour, stat = "density", adjust = adjust.bin) + # pop colored
    #   scale_colour_manual(name ="Sampling sites", values = colour_palette_sites.pink) +
    #     scale_x_continuous(breaks=c(0, 0.05, 0.1, 0.2, 0.5, 1),
    #                        labels = c("0", "0.05", "0.1", "0.2", "0.5", "1.0")) +
    ggplot2::labs(x = x.title) +
    ggplot2::labs(y = "Density of SNP (scaled)") +
    ggplot2::expand_limits(y = 0) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 12, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 12, face = "bold"),
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.text = ggplot2::element_text(size = 12, face = "bold"),
      strip.text.y = ggplot2::element_text(angle = 0, size = 12, face = "bold"),
      strip.text.x = ggplot2::element_text(size = 12, face = "bold")
    )

}



#' @title Figure of the distribution of SNP per locus before and after filters
#' @description Create density distribution of SNP per locus before and after filters.
#' @param before.filter.data Data set before filter.
#' @param after.filter.data Data set after filter.
#' @export
#' @rdname plot_snp_number_loci
#' @keywords internal

plot_snp_number_loci <- function(before.filter.data, after.filter.data) {

  GROUP <- NULL
  SNP_N <- NULL
  POP_ID <- NULL


  if (is.vector(before.filter.data)) {
    before.filter.data <- readr::read_tsv(before.filter.data, col_names = T)
  } else {
    before.filter.data <- before.filter.data
  }

  if (is.vector(after.filter.data)) {
    after.filter.data <- readr::read_tsv(after.filter.data, col_names = T)
  } else {
    after.filter.data <- after.filter.data
  }

  number.snp.loci <- before.filter.data %>% # Before
    dplyr::group_by(LOCUS) %>%
    dplyr::summarise(SNP_N = dplyr::n_distinct(POS)) %>%
    dplyr::mutate(GROUP = rep("pre-filters", n())) %>%
    dplyr::bind_rows(
      after.filter.data %>% # After
        dplyr::group_by(LOCUS) %>%
        dplyr::summarise(SNP_N = dplyr::n_distinct(POS)) %>%
        dplyr::mutate(GROUP = rep("post-filters", n()))) %>%
    dplyr::mutate(
      GROUP = factor(GROUP, levels = c("pre-filters", "post-filters"),
                     ordered = TRUE)
    )

  graph <- ggplot2::ggplot(number.snp.loci, ggplot2::aes(factor(SNP_N))) +
    ggplot2::geom_bar() +
    ggplot2::labs(x = "Number of SNP per haplotypes") +
    ggplot2::labs(y = "Distribution (number)") +
    ggplot2::facet_wrap(~GROUP, nrow = 1, ncol = 2) +
    ggplot2::theme(axis.title.x = ggplot2::element_text(size = 12, face = "bold"),
                   axis.title.y = ggplot2::element_text(size = 12, face = "bold"),
                   legend.title = ggplot2::element_text(size = 12, face = "bold"),
                   legend.text = ggplot2::element_text(size = 12, face = "bold"),
                   strip.text.x = ggplot2::element_text(size = 12, face = "bold"))

  graph
}




#' @title Figure of the distribution of SNP nucleotide position alond the read
#' @description Distribution of SNP nucleotide position alond the read.
#' @param data Data for the figure.
#' @param aes.colour GGPLOT2 aesthetic.
#' @param y.title Title of the Y-axis.
#' @export
#' @rdname plot_snp_position_read
#' @keywords internal



plot_snp_position_read <- function(data, aes.colour, y.title) {

  COL <- NULL


  ggplot2::ggplot(data, ggplot2::aes(x = COL, na.rm = TRUE)) +
    ggplot2::geom_line((aes.colour), stat = "density", size = 1, adjust = 0.7) +
    ggplot2::scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90)) +
    ggplot2::labs(x = "Nucleotide position on the read") +
    ggplot2::labs(y = y.title) +
    ggplot2::expand_limits(x = 0) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 12, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 12, face = "bold"),
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.text = ggplot2::element_text(size = 12, face = "bold"),
      strip.text.x = ggplot2::element_text(size = 12, face = "bold")
    )
}




#' @title Density distribution of diversity (Gene and Haplotypes)
#' @description GGPLOT2 Density distribution of diversity (Gene and Haplotypes).
#' @param data The hapstats summary file or object.
#' @param aes.x GGPLOT2 aesthetics,
#' e.g. aes.x = aes(x = GENE_DIVERSITY, na.rm = T).
#' @param aes.colour GGPLOT2 aesthetics colour,
#' e.g. aes.colour = aes(y = ..scaled.., colour = POP_ID).
#' @param x.title Title of the x-axis.
#' @param y.title Title of the y-axis.
#' @export
#' @rdname plot_distribution_diversity
#' @keywords internal

plot_distribution_diversity <- function(data, aes.x, aes.colour, x.title, y.title) {

  hapstats.summary <- NULL


  ggplot2::ggplot(hapstats.summary, aes.x) +
    ggplot2::geom_line(aes.colour, stat = "density", adjust = 0.8) +
    #   scale_colour_manual(name = "Populations", values = colour_palette_sites.pink, breaks = c("BUR", "GRA", "GUL", "LLI", "ANG", "WEI", "HAY", "GOD")) +
    #   geom_density(aes(fill=POP_ID, color=NA), alpha=0.4) +
    #   scale_fill_manual(name="Populations", values=colour_palette_sites.pink, breaks = c("BUR", "GRA", "GUL", "LLI", "ANG", "WEI", "HAY", "GOD")) +
    #   scale_color_manual(name="Populations", values=colour_palette_sites.pink, breaks = c("BUR", "GRA", "GUL", "LLI", "ANG", "WEI", "HAY", "GOD")) +
    ggplot2::labs(x = x.title) +
    ggplot2::labs(y = y.title) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_text(size = 12,face = "bold"),
      axis.title.y = ggplot2::element_text(size = 12,face = "bold"),
      legend.title = ggplot2::element_text(size = 12,face = "bold"),
      legend.text = ggplot2::element_text(size = 12,face = "bold"),
      strip.text.y = ggplot2::element_text(angle = 0,size = 12,face = "bold"),
      strip.text.x = ggplot2::element_text(size = 12,face = "bold")
    )
}




#' @title Box plot of the diversity (Gene and Haplotypes)
#' @description GGPLOT2 Box plot of the diversity (Gene and Haplotypes).
#' @param data The hapstats summary file or object.
#' @param aes.x.y The GGPLOT2 aesthetics,
#' e.g. aes.x.y = aes(x = factor(POP_ID), y = GENE_DIVERSITY, na.rm = T).
#' @param y.title Title of the y-axis.
#' @export
#' @rdname plot_boxplot_diversity
#' @keywords internal
plot_boxplot_diversity <- function(data, aes.x.y, y.title) {
  ggplot2::ggplot(data, aes.x.y) +
    ggplot2::geom_violin(trim = F) +
    ggplot2::geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA) +
    ggplot2::stat_summary(fun = "mean", geom = "point", shape = 21, size = 2.5, fill = "white") +
    ggplot2::labs(x = "Sampling sites") +
    ggplot2::labs(y = y.title) +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_text(size = 12, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 12, face = "bold"),
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.text = ggplot2::element_text(size = 12, face = "bold"),
      strip.text.x = ggplot2::element_text(size = 12, face = "bold")
    )
}

