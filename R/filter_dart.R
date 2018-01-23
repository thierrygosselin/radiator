# Import, filter and transform a dart output file to different formats

#' @name filter_dart

#' @title Swiss Army knife tool to prepare \href{http://www.diversityarrays.com}{DArT}
#' output file for population genetics analysis.
#' The function uses \code{\link[radiator]{tidy_dart}} to import and tidy DArT input file.
#' Currently 3 formats are recognized: 1- and 2- row format (also called binary),
#' and count data.

#' @description Import, filter and generate imputed dataset of DArT output file.

#' @param interactive.filter (optional, logical) Do you want the filtering session to
#' be interactive. The user is asked to see figures of distribution before
#' making decisions for filtering.
#' Default: \code{interactive.filter = TRUE}.

#' @inheritParams tidy_dart
#' @inheritParams genomic_converter
#' @inheritParams tidy_genomic_data

#' @param strata A tab delimited file with columns header:
#' \code{INDIVIDUALS} and \code{STRATA}.
#' Note: the column \code{STRATA} refers to any grouping of individuals. If a third column
#' named \code{NEW_ID} is used, this column will be used to replace the
#' \code{INDIVIDUALS} column in the main data file.

#' @param filter.reproducibility (optional, numerical) Filter the \code{RepAvg}
#' column in the data set. Default: \code{filter.reproducibility = NULL}.
#' e.g to keep markers with reproducibility >= 99%,
#' use: \code{filter.reproducibility = 0.99}.

#' @param filter.call.rate (optional, numerical) Filter the \code{CallRate}
#' column in the data set. Default: \code{filter.call.rate = NULL}. e.g to keep
#' markers genotyped in more than 95% of the individuals use :
#' \code{filter.call.rate = 0.95}

#' @param filter.coverage (optional, string, numerical) Filter the lower and
#' upper bound of locus/read coverage. The locus/read coverage combines the markers
#' average count for REF and ALT allele (respectively the \code{AvgCountRef} and
#' \code{AvgCountSnp} info).
#' Default: \code{filter.coverage = NULL}.
#' e.g to keep markers with coverage inbetween 7 and 200,
#' use : \code{filter.coverage = c(7, 200)}.

#' @param filter.ind.missing.geno (optional, string) Similar to call rate, but
#' more adapted to the data. 3 values are required in the string, corresponding
#' to the \code{\link[radiator]{filter_individual}} module of radiator.
#' First value is the approach to count genotyped individuals, \code{"overall"}
#' or by \code{"pop"}. Second value is the percent threshold for the marker, with
#' \code{70}, 70 percent of genotyped individuals are required to keep the marker.
#' The last threshold is the number of problematic population that are allowed to skip
#' the threshold. In doubt, use the interactive mode that take step by step these
#' arguments. e.g to keep individuals genotyped at >= 70 percent for the markers,
#' without considering the population info and allowing 1 population to be problematic for the
#' threshold, use: \code{c("overall", 70, 1)}.
#' Default: \code{filter.ind.missing.geno = NULL}.

#' @param number.snp.reads (optional, integer) This filter removes outlier markers
#' with too many SNP number per locus/read.
#' Having a higher than "normal" SNP number is usually the results of
#' assembly artifacts or bad assembly parameters.
#' This filter is population-agnostic. This is best decide after viewing the figures,
#' with the interactive mode.
#' If the argument is set to \code{number.snp.reads = 2},
#' locus with 3 and more SNPs will be blacklisted.
#' Default: \code{number.snp.reads = NULL}.

#' @param mixed.genomes.analysis (optional, logical) Highlight outliers individual's
#' observed heterozygosity for a quick
#' diagnostic of mixed samples or poor polymorphism discovery due to DNA quality,
#' sequencing effort, etc.
#' See this function for more info: \code{\link[radiator]{detect_mixed_genomes}}.
#' Default: \code{detect_mixed_genomes = TRUE}.

#' @inheritParams detect_mixed_genomes

#' @param duplicate.genomes.analysis (optional, string) Detect duplicate individuals.
#' The function can compute two methods (distance or genome pairwise similarity)
#' to highligh potential duplicate individuals.
#' See this function for more info: \code{\link[radiator]{detect_duplicate_genomes}}.
#' The string required to run the analysis as 2 values:
#' \enumerate{
#' \item TRUE/FALSE to run the analysis;
#' \item Computes pairwise genome similarity (TRUE/FALSE),
#' with FALSE just the distance measure is used.
#' The pairwise genome similarity is longer to run, but is better because it
#' integrates markers in common/missing data.
#' Using \code{interactive.filter = TRUE}, can overide this value,
#' you can opt in for the pairwise genome similarity after viewing the figures
#' used with distance measure... handy!
#' }
#' Default: \code{detect_duplicate_genomes = c(TRUE, FALSE)}.


#' @param filename (optional) The filename prefix for the objet in the global environment
#' or the working directory. Default: \code{filename = NULL}. A default name will be used,
#' customized with the output file(s) selected.


#' @return The function returns an object (list). The content of the object
#' can be listed with \code{names(object)} and use \code{$} to isolate specific
#' object (see examples). Some output format will write the output file in the
#' working directory. The tidy genomic data frame is generated automatically.

#' @export
#' @rdname filter_dart
#' @importFrom dplyr group_by select rename filter mutate summarise distinct n_distinct arrange left_join semi_join anti_join inner_join full_join tally bind_rows
#' @importFrom parallel detectCores
#' @importFrom stringi stri_replace_all_fixed stri_join stri_sub stri_replace_na stri_pad_left
#' @importFrom purrr discard
#' @importFrom readr read_tsv write_tsv
#' @importFrom tibble as_data_frame data_frame
#' @importFrom tidyr spread gather unite separate
# @importFrom fst write.fst read.fst

#' @examples
#' \dontrun{
#' testing
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com} and Peter Grewe \email{peter.grewe@csiro.au}

filter_dart <- function(
  interactive.filter = TRUE,
  data,
  strata,
  output = NULL,
  pop.levels = NULL,
  blacklist.id = NULL,
  pop.select = NULL,
  monomorphic.out = TRUE,
  common.markers = TRUE,
  filter.reproducibility = NULL,
  filter.coverage = NULL,
  filter.call.rate = NULL,
  filter.ind.missing.geno = NULL,
  number.snp.reads = NULL,
  maf.thresholds = NULL,
  mixed.genomes.analysis = TRUE,
  ind.heterozygosity.threshold = NULL,
  duplicate.genomes.analysis = c(TRUE, FALSE),
  snp.ld = NULL,
  imputation.method = NULL,
  hierarchical.levels = "populations",
  num.tree = 50,
  filename = NULL,
  verbose = TRUE,
  parallel.core = parallel::detectCores() - 1
) {
  if (verbose) {
    cat("#######################################################################\n")
    cat("######################## radiator::filter_dart ########################\n")
    cat("#######################################################################\n")
  }
  opt.change <- getOption("width")
  options(width = 70)
  timing <- proc.time()# for timing
  res <- list()

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file missing")
  if (missing(output)) stop("At least 1 output format is required")

  # Filename -------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")

  if (is.null(filename)) {
    folder.extension <- filename <- stringi::stri_join("filter_dart_", file.date)

    # if (!is.null(imputation.method)) {
    #   filename.imp <- stringi::stri_join("filter_dart_imputed_", file.date)
    # }
  } else {
    # if (!is.null(imputation.method)) {
    #   filename.imp <- stringi::stri_join(filename, "imputed", file.date, sep = "_")
    # }
    folder.extension <- stringi::stri_join("filter_dart", filename, file.date, sep = "_")
    filename <- stringi::stri_join(filename, file.date, sep = "_")
  }

  working.dir <- getwd()
  path.folder <- stringi::stri_join(working.dir,"/", folder.extension, sep = "")
  dir.create(file.path(path.folder))

  message("\nFolder created: \n", folder.extension)
  filename <- file.path(path.folder, filename)

  # Filter parameter file ------------------------------------------------------
  filters.parameters <- tibble::data_frame(
    FILTERS = as.character(),
    PARAMETERS = as.character(),
    VALUES = as.integer(),
    BEFORE = as.character(),
    AFTER = as.character(),
    BLACKLIST = as.integer(),
    UNITS = as.character(),
    COMMENTS = as.character())
  filters.parameters.path <- stringi::stri_join(
    path.folder, "/filters_parameters_dart.tsv")
  readr::write_tsv(x = filters.parameters,
                   path = filters.parameters.path,
                   append = FALSE,
                   col_names = TRUE)
  message("Generated a filters parameters file: filters_parameters_dart.tsv")

  # Markers meta----------------------------------------------------------------
  message("Importing markers metadata")
  metadata <- radiator::tidy_dart_metadata(
    data = data,
    filename = filename,
    parallel.core = parallel.core,
    verbose = TRUE)

  metadata.file <- list.files(path = path.folder, pattern = "metadata")

  # create 2 data.info
  data.info <- first.data.info <- data_info(metadata, print.info = FALSE)


  # message("interactive.filter: ", interactive.filter)
  # message("filter.reproducibility: ", filter.reproducibility)
  # Filtering reproducibility  -------------------------------------------------
  if (interactive.filter || !is.null(filter.reproducibility)) {
    message("Filtering reproducibility")

    folder.extension <- stringi::stri_join("filter_reproducibility_dart_", file.date, sep = "")
    path.folder.reproducibility <- file.path(path.folder, folder.extension)
    dir.create(path.folder.reproducibility)
    message("Folder created: \n", folder.extension)

    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "REP_AVG")
    reproducibility <- dplyr::select(metadata, dplyr::one_of(want)) %>%
      dplyr::mutate(Markers = rep("markers", n()))

    plot.reproducibility.violinplot <- ggplot2::ggplot(
      reproducibility, ggplot2::aes(x = Markers, y = REP_AVG, na.rm = TRUE)) +
      ggplot2::geom_violin(trim = TRUE) +
      ggplot2::geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA) +
      ggplot2::stat_summary(fun.y = "mean", geom = "point", shape = 21,
                            size = 2.5, fill = "white") +
      ggplot2::labs(x = "Markers", y = "Markers reproducibility averaged") +
      ggplot2::theme(
        legend.position = "none",
        axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.x = ggplot2::element_blank()
      )
    # plot.reproducibility.violinplot
    ggplot2::ggsave(
      filename = file.path(path.folder.reproducibility, "plot.reproducibility.violinplot.pdf"),
      plot = plot.reproducibility.violinplot,
      width = 20, height = 15, dpi = 600, units = "cm", useDingbats = FALSE)

    plot.reproducibility.histo <- ggplot2::ggplot(
      data = reproducibility,
      ggplot2::aes(x = REP_AVG)) +
      ggplot2::geom_histogram() +
      ggplot2::labs(x = "Markers", y = "Markers reproducibility averaged") +
      ggplot2::theme(
        legend.position = "none",
        axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5),
        strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
      )
    # plot.reproducibility.histo
    suppressMessages(ggplot2::ggsave(
      filename = file.path(path.folder.reproducibility, "plot.reproducibility.histo.pdf"),
      plot = plot.reproducibility.histo,
      width = 20, height = 15, dpi = 600, units = "cm", useDingbats = FALSE))

    if (interactive.filter) {
      message("    Inspect plots in folder created to help choose reproducibility threshold...")
      message("    Enter the value (between 0 and 1) for filter.reproducibility threshold \n    (below threshold < markers are discarded): ")
      filter.reproducibility <- as.numeric(readLines(n = 1))
    }

    if (!is.null(filter.reproducibility)) {
      want <- c("MARKERS", "CHROM", "LOCUS", "POS")

      blacklist.reproducibility.markers <- reproducibility %>%
        dplyr::filter(REP_AVG < filter.reproducibility) %>%
        dplyr::select(dplyr::one_of(want)) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE)

      if (nrow(blacklist.reproducibility.markers) > 0) {
        n.snp.before <- data.info$n.snp
        if (verbose) message("    Number of markers before = ", n.snp.before)
        n.snp.blacklist <- nrow(blacklist.reproducibility.markers)
        if (verbose) message("    Number of markers removed = ", n.snp.blacklist)
        if (verbose) message("    Number of markers after = ", n.snp.before - n.snp.blacklist)
        # input <- dplyr::anti_join(input, blacklist.reproducibility.markers, by = "MARKERS")
        # whitelist.markers <- dplyr::select(input, dplyr::one_of(want)) %>%
        #   dplyr::distinct(MARKERS, .keep_all = TRUE)

        metadata <- metadata %>%
          dplyr::filter(!MARKERS %in% blacklist.reproducibility.markers$MARKERS)

        new.data.info <- data_info(metadata) # updating parameters

        filters.parameters <- tibble::data_frame(
          FILTERS = "reproducibility",
          PARAMETERS = "",
          VALUES = filter.reproducibility,
          BEFORE = stringi::stri_join(data.info$n.chrom, data.info$n.locus, data.info$n.snp, sep = "/"),
          AFTER = stringi::stri_join(new.data.info$n.chrom, new.data.info$n.locus, new.data.info$n.snp, sep = "/"),
          BLACKLIST = stringi::stri_join(data.info$n.chrom - new.data.info$n.chrom, data.info$n.locus - new.data.info$n.locus, data.info$n.snp - new.data.info$n.snp, sep = "/"),
          UNITS = "CHROM/LOCUS/SNP",
          COMMENTS = ""
        )
        readr::write_tsv(x = filters.parameters,
                         path = filters.parameters.path, append = TRUE,
                         col_names = FALSE)
        # update data.info
        data.info <- new.data.info

        blacklist.markers <- blacklist.reproducibility.markers
        whitelist.markers <- dplyr::select(metadata, MARKERS, CHROM, LOCUS, POS)

        # write blacklist and whitelist
        readr::write_tsv(x = blacklist.reproducibility.markers, path = file.path(path.folder.reproducibility, "blacklist.reproducibility.markers.tsv"))
        readr::write_tsv(x = whitelist.markers, path = file.path(path.folder.reproducibility, "whitelist.reproducibility.markers.tsv"))

      }
      blacklist.reproducibility.markers <- NULL
      whitelist.markers <- dplyr::select(metadata, MARKERS, CHROM, LOCUS, POS)
    } else {
      stop("A filter.reproducibility threshold value is required...")
    }

    reproducibility <- plot.reproducibility.histo <- plot.reproducibility.violinplot <- NULL
  }#End filter.reproducibility

  # Filtering call rate ---------------------------------------------------------
  if (interactive.filter || !is.null(filter.call.rate)) {
    message("Filtering markers call rate")

    folder.extension <- stringi::stri_join("filter_call_rate_dart_", file.date, sep = "")
    path.folder.call.rate <- file.path(path.folder, folder.extension)
    dir.create(path.folder.call.rate)
    message("Folder created: \n", folder.extension)

    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "CALL_RATE")
    call.rate <- dplyr::select(metadata, dplyr::one_of(want)) %>%
      dplyr::mutate(Markers = rep("markers", n()))

    plot.call.rate.violinplot <- ggplot2::ggplot(
      call.rate, ggplot2::aes(x = Markers, y = CALL_RATE, na.rm = TRUE)) +
      ggplot2::geom_violin(trim = TRUE) +
      ggplot2::geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA) +
      ggplot2::stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white") +
      ggplot2::labs(x = "Markers") +
      ggplot2::labs(y = "Markers Call Rate (proportion)") +
      ggplot2::theme(
        legend.position = "none",
        # axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.x = ggplot2::element_blank()
        # axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5),
        # legend.title = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        # legend.text = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
        # strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
      ) #+ ggplot2::facet_grid(~POP_ID)
    # plot.call.rate.violinplot
    ggplot2::ggsave(
      filename = file.path(path.folder.call.rate, "plot.call.rate.violinplot.pdf"),
      plot = plot.call.rate.violinplot,
      width = 20, height = 15, dpi = 600, units = "cm", useDingbats = FALSE)

    plot.call.rate.histo <- ggplot2::ggplot(
      data = call.rate,
      ggplot2::aes(x = CALL_RATE)) +
      ggplot2::geom_histogram() +
      ggplot2::labs(x = "Markers") +
      ggplot2::labs(y = "Markers Call Rate (proportion)") +
      ggplot2::theme(
        legend.position = "none",
        axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5),
        strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
      )
    # plot.call.rate.histo
    suppressMessages(ggplot2::ggsave(
      filename = file.path(path.folder.call.rate, "plot.call.rate.histo.pdf"),
      plot = plot.call.rate.histo,
      width = 20, height = 15, dpi = 600, units = "cm", useDingbats = FALSE))

    if (interactive.filter) {
      message("    Inspect plots in folder created to help choose call rate threshold...")
      message("    Enter the value (between 0 and 1) for filter.call.rate threshold \n    (below threshold < markers are discarded): ")
      filter.call.rate <- as.numeric(readLines(n = 1))
    }

    if (!is.null(filter.call.rate)) {
      want <- c("MARKERS", "CHROM", "LOCUS", "POS")
      blacklist.call.rate.markers <- call.rate %>%
        dplyr::filter(CALL_RATE < filter.call.rate) %>%
        dplyr::select(dplyr::one_of(want)) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE)

      if (nrow(blacklist.call.rate.markers) > 0) {
        n.snp.before <- data.info$n.snp
        if (verbose) message("    Number of markers before = ", n.snp.before)
        n.snp.blacklist <- nrow(blacklist.call.rate.markers)
        if (verbose) message("    Number of markers removed = ", n.snp.blacklist)
        if (verbose) message("    Number of markers after = ", n.snp.before - n.snp.blacklist)

        # input <- dplyr::anti_join(input, blacklist.call.rate.markers, by = "MARKERS")

        metadata <- metadata %>%
          dplyr::filter(!MARKERS %in% blacklist.call.rate.markers$MARKERS)

        new.data.info <- data_info(metadata) # updating parameters
        filters.parameters <- tibble::data_frame(
          FILTERS = "call rate",
          PARAMETERS = "",
          VALUES = filter.call.rate,
          BEFORE = stringi::stri_join(data.info$n.chrom, data.info$n.locus, data.info$n.snp, sep = "/"),
          AFTER = stringi::stri_join(new.data.info$n.chrom, new.data.info$n.locus, new.data.info$n.snp, sep = "/"),
          BLACKLIST = stringi::stri_join(data.info$n.chrom - new.data.info$n.chrom, data.info$n.locus - new.data.info$n.locus, data.info$n.snp - new.data.info$n.snp, sep = "/"),
          UNITS = "CHROM/LOCUS/SNP",
          COMMENTS = ""
        )
        readr::write_tsv(x = filters.parameters,
                         path = filters.parameters.path, append = TRUE,
                         col_names = FALSE)
        # update data.info
        data.info <- new.data.info

        if (!is.null(blacklist.markers)) {
          blacklist.markers <- dplyr::bind_rows(blacklist.markers, blacklist.call.rate.markers)
        } else {
          blacklist.markers <- blacklist.call.rate.markers
        }

        whitelist.markers <- dplyr::select(metadata, dplyr::one_of(want)) %>%
          dplyr::distinct(MARKERS, .keep_all = TRUE)

        # write blacklist and whitelist
        readr::write_tsv(x = blacklist.call.rate.markers, path = file.path(path.folder.call.rate, "blacklist.call.rate.markers.tsv"))
        readr::write_tsv(x = whitelist.markers, path = file.path(path.folder.call.rate, "whitelist.call.rate.markers.tsv"))
      }
      blacklist.call.rate.markers <- NULL
    } else {
      stop("A filter.call.rate threshold value is required...")
    }
    call.rate <- plot.call.rate.histo <- plot.call.rate.violinplot <- NULL
  }#End filter.call.rate

  # Filtering coverage --------------------------------------------------------
  if (interactive.filter || !is.null(filter.coverage)) {
    message("Filtering markers mean coverage")

    folder.extension <- stringi::stri_join("filter_coverage_dart_", file.date, sep = "")
    path.folder.coverage <- file.path(path.folder, folder.extension)
    dir.create(path.folder.coverage)
    message("Folder created: \n", folder.extension)

    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "AVG_COUNT_REF", "AVG_COUNT_SNP")

    coverage <- dplyr::select(metadata, dplyr::one_of(want)) %>%
      dplyr::rename(REF = AVG_COUNT_REF, ALT = AVG_COUNT_SNP) %>%
      dplyr::distinct(MARKERS, REF, ALT, .keep_all = TRUE) %>%
      dplyr::mutate(READS = REF + ALT)

    coverage.long <- tidyr::gather(data = coverage, key = GROUP, value = COVERAGE, -c(MARKERS, CHROM, LOCUS, POS)) %>%
      dplyr::mutate(GROUP = factor(GROUP, levels = c("READS", "REF", "ALT"), ordered = TRUE))
    plot.coverage <- ggplot2::ggplot(
      data = coverage.long,
      ggplot2::aes(x = COVERAGE)) +
      ggplot2::geom_histogram() +
      ggplot2::labs(x = "Markers") +
      ggplot2::labs(y = "Coverage (count)") +
      ggplot2::theme(
        legend.position = "none",
        axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5),
        strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
      ) +
      ggplot2::facet_grid(~ GROUP)
    # plot.coverage
    suppressMessages(ggplot2::ggsave(
      filename = file.path(path.folder.coverage, "plot.coverage.pdf"),
      plot = plot.coverage,
      width = 30, height = 15, dpi = 600, units = "cm", useDingbats = FALSE))

    plot.coverage.boxplot <- ggplot2::ggplot(
      coverage.long,
      ggplot2::aes(x = factor(GROUP), y = COVERAGE, na.rm = TRUE)) +
      ggplot2::geom_violin(trim = TRUE) +
      ggplot2::geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA) +
      ggplot2::stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white") +
      ggplot2::labs(x = "Markers") +
      ggplot2::labs(y = "Coverage (count)") +
      ggplot2::theme(
        legend.position = "none",
        axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica"),
        legend.title = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        legend.text = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
      )
    # plot.coverage.boxplot
    suppressMessages(ggplot2::ggsave(
      filename = file.path(path.folder.coverage, "plot.coverage.boxplot.pdf"),
      plot = plot.coverage.boxplot,
      width = 15, height = 15, dpi = 600, units = "cm", useDingbats = FALSE))

    if (interactive.filter) {
      filter.coverage <- c(1, 1000000)
      message("    Inspect plots in folder created to help choose coverage thresholds (min and max)...")
      message("    Enter the minimum coverage allowed to keep a marker (e.g. 7): ")
      filter.coverage[1] <- as.numeric(readLines(n = 1))
    }

    if (interactive.filter) {
      message("    Enter the maximum coverage allowed to keep a marker (e.g. 150): ")
      filter.coverage[2] <- as.numeric(readLines(n = 1))
      filter.coverage <- as.numeric(filter.coverage)
    }

    if (!is.null(filter.coverage)) {
      want <- c("MARKERS", "CHROM", "LOCUS", "POS")
      blacklist.coverage.markers <- coverage %>%
        dplyr::filter(READS < filter.coverage[1] | READS > filter.coverage[2]) %>%
        dplyr::select(dplyr::one_of(want)) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE)

      if (nrow(blacklist.coverage.markers) > 0) {
        n.snp.before <- data.info$n.snp
        if (verbose) message("    Number of markers before = ", n.snp.before)
        n.snp.blacklist <- nrow(blacklist.coverage.markers)
        if (verbose) message("    Number of markers removed = ", n.snp.blacklist)
        if (verbose) message("    Number of markers after = ", n.snp.before - n.snp.blacklist)
        # input <- dplyr::anti_join(input, blacklist.coverage.markers, by = "MARKERS")

        metadata <- metadata %>%
          dplyr::filter(!MARKERS %in% blacklist.coverage.markers$MARKERS)

        whitelist.markers <- dplyr::select(metadata, dplyr::one_of(want)) %>%
          dplyr::distinct(MARKERS, .keep_all = TRUE)

        new.data.info <- data_info(metadata) # updating parameters
        filters.parameters <- tibble::data_frame(
          FILTERS = "coverage",
          PARAMETERS = "min/max",
          VALUES = stringi::stri_join(filter.coverage, collapse = "/"),
          BEFORE = stringi::stri_join(data.info$n.chrom, data.info$n.locus, data.info$n.snp, sep = "/"),
          AFTER = stringi::stri_join(new.data.info$n.chrom, new.data.info$n.locus, new.data.info$n.snp, sep = "/"),
          BLACKLIST = stringi::stri_join(data.info$n.chrom - new.data.info$n.chrom, data.info$n.locus - new.data.info$n.locus, data.info$n.snp - new.data.info$n.snp, sep = "/"),
          UNITS = "CHROM/LOCUS/SNP",
          COMMENTS = ""
        )
        readr::write_tsv(x = filters.parameters,
                         path = filters.parameters.path, append = TRUE,
                         col_names = FALSE)
        # update data.info
        data.info <- new.data.info

        if (!is.null(blacklist.markers)) {
          blacklist.markers <- dplyr::bind_rows(blacklist.markers, blacklist.coverage.markers)
        } else {
          blacklist.markers <- blacklist.coverage.markers
        }

        # write blacklist and whitelist
        readr::write_tsv(x = blacklist.coverage.markers, path = file.path(path.folder.coverage, "blacklist.coverage.markers.tsv"))
        readr::write_tsv(x = whitelist.markers, path = file.path(path.folder.coverage, "whitelist.coverage.markers.tsv"))
      }

      blacklist.coverage.markers <- NULL
    } else {
      stop("A filter.coverage thresholds values are required...")
    }
    coverage <- coverage.long <- plot.coverage <- plot.coverage.boxplot <- NULL
  }#End filter.coverage

  # update metadata file -------------------------------------------------------
  fst::write.fst(
    x = metadata,
    path = file.path(path.folder, stringi::stri_replace_all_fixed(metadata.file, ".rad", "_filtered.rad", vectorize_all = FALSE)),
    compress =85)

  # readr::write_tsv(
  # x = metadata,
  # path = file.path(path.folder, stringi::stri_replace_all_fixed(metadata.file, ".rad", "_filtered.rad", vectorize_all = FALSE))
  # )

  # Tidy DArT ------------------------------------------------------------------
  message("\nNext step requires the genotypes")

  input <- radiator::tidy_dart(data = data, strata = strata,
                               whitelist.markers = whitelist.markers,
                               filename = filename, verbose = TRUE,
                               parallel.core = parallel.core)

  # Change populations names or order/levels -----------------------------------
  input <- change_pop_names(data = input, pop.levels = pop.levels)

  # Blacklist id ---------------------------------------------------------------
  if (is.null(blacklist.id)) { # No blacklist of ID
    if (verbose) message("\nBlacklisted individuals: no")
    res$blacklist.id <- tibble::data_frame(INDIVIDUALS = character(0))
  } else {# With blacklist of ID
    res$blacklist.id <- suppressMessages(readr::read_tsv(blacklist.id, col_names = TRUE))
    n.ind.blacklisted <- length(res$blacklist.id$INDIVIDUALS)
    if (verbose) message("\nBlacklisted individuals: ", n.ind.blacklisted, " ind.")
    if (verbose) message("    Filtering with blacklist of individuals")

    n.id.before <- dplyr::n_distinct(input$INDIVIDUALS)
    n.pop.before <- dplyr::n_distinct(input$POP_ID)
    # system.time(input2 <- suppressWarnings(dplyr::anti_join(input, res$blacklist.id, by = "INDIVIDUALS")))
    input <- dplyr::filter(input, !INDIVIDUALS %in% res$blacklist.id$INDIVIDUALS)
    n.id.after <- dplyr::n_distinct(input$INDIVIDUALS)
    n.pop.after <- dplyr::n_distinct(input$POP_ID)

    n.ind.blacklisted <- n.id.after - n.id.before
    if (n.ind.blacklisted > 0) {
      # updating parameters
      filters.parameters <- tibble::data_frame(
        FILTERS = "blacklisted ids",
        PARAMETERS = "",
        VALUES = "",
        BEFORE = n.id.before,
        AFTER = n.id.after,
        BLACKLIST = n.id.before - n.id.after,
        UNITS = "individuals",
        COMMENTS = ""
      )
      readr::write_tsv(x = filters.parameters, path = filters.parameters.path,
                       append = TRUE, col_names = FALSE)
    } else {
      message("    None of the individuals in the blacklist were used to filter the data")
    }
    # update data.info
    data.info$n.ind <- n.id.after
    data.info$n.pop <- n.pop.after
  }

  # pop.select -----------------------------------------------------------------
  if (!is.null(pop.select)) {
    if (verbose) {
      n.pop.new <- length(pop.select)
      message(
        n.pop.new, " population(s) selected: ",
        stringi::stri_join(pop.select, collapse = ", "), sep = " ")
    }
    input <- suppressWarnings(input %>% dplyr::filter(POP_ID %in% pop.select))

    # updating parameters
    filters.parameters <- tibble::data_frame(
      FILTERS = "pop selected",
      PARAMETERS = "",
      VALUES = stringi::stri_join(pop.select, collapse = ", "),
      BEFORE = data.info$n.pop,
      AFTER = n.pop.new,
      BLACKLIST = data.info$n.pop - n.pop.new,
      UNITS = "populations",
      COMMENTS = ""
    )
    readr::write_tsv(x = filters.parameters, path = filters.parameters.path,
                     append = TRUE, col_names = FALSE)
    # update data.info
    data.info$n.pop <- n.pop.new
  }
  input$POP_ID <- droplevels(input$POP_ID)
  pop.levels <- levels(input$POP_ID)

  # Filter monomorphic markers  ------------------------------------------------
  if (monomorphic.out) {
    input <- discard_monomorphic_markers(data = input, verbose = verbose)
    blacklist.monomorphic.markers <- input$blacklist.monomorphic.markers
    whitelist.markers <- input$whitelist.polymorphic.markers
    input <- input$input

    # update blacklist.markers
    if (nrow(blacklist.monomorphic.markers) > 0) {
      readr::write_tsv(x = blacklist.monomorphic.markers,
                       path = file.path(path.folder, "blacklist.monomorphic.markers.tsv"))
      if (verbose) message("    Blacklist of monomorphic markers written in: blacklist.monomorphic.markers.tsv")

      if (is.null(blacklist.markers)) {
        blacklist.markers<- blacklist.monomorphic.markers
      } else {
        blacklist.markers <- dplyr::bind_rows(blacklist.markers, blacklist.monomorphic.markers)
      }

      new.data.info <- data_info(input)

      # updating parameters
      filters.parameters <- tibble::data_frame(
        FILTERS = "removing monomorphic markers",
        PARAMETERS = "",
        VALUES = "",
        BEFORE = stringi::stri_join(data.info$n.chrom, data.info$n.locus, data.info$n.snp, sep = "/"),
        AFTER = stringi::stri_join(new.data.info$n.chrom, new.data.info$n.locus, new.data.info$n.snp, sep = "/"),
        BLACKLIST = stringi::stri_join(data.info$n.chrom - new.data.info$n.chrom, data.info$n.locus - new.data.info$n.locus, data.info$n.snp - new.data.info$n.snp, sep = "/"),
        UNITS = "CHROM/LOCUS/SNP",
        COMMENTS = ""
      )
      readr::write_tsv(x = filters.parameters,
                       path = filters.parameters.path, append = TRUE,
                       col_names = FALSE)
      # update data.info
      data.info <- new.data.info
    }
    blacklist.monomorphic.markers <- NULL
  }# End monomorphic.out

  # Filter common markers between all populations  ------------------------------
  if (common.markers) {
    input <- keep_common_markers(data = input, plot = FALSE,
                                 verbose = verbose)
    blacklist.not.in.common.markers <- input$blacklist.not.in.common.markers
    whitelist.markers <- input$whitelist.common.markers
    input <- input$input

    if (nrow(blacklist.not.in.common.markers) > 0) {
      readr::write_tsv(x = blacklist.not.in.common.markers,
                       path = file.path(path.folder, "blacklist.not.in.common.markers.tsv"))
      if (verbose) message("    Blacklist of markers not in common written in: blacklist.not.in.common.markers.tsv")

      new.data.info <- data_info(input) # updating parameters
      filters.parameters <- tibble::data_frame(
        FILTERS = "keeping common markers",
        PARAMETERS = "",
        VALUES = "",
        BEFORE = stringi::stri_join(data.info$n.chrom, data.info$n.locus, data.info$n.snp, sep = "/"),
        AFTER = stringi::stri_join(new.data.info$n.chrom, new.data.info$n.locus, new.data.info$n.snp, sep = "/"),
        BLACKLIST = stringi::stri_join(data.info$n.chrom - new.data.info$n.chrom, data.info$n.locus - new.data.info$n.locus, data.info$n.snp - new.data.info$n.snp, sep = "/"),
        UNITS = "CHROM/LOCUS/SNP",
        COMMENTS = ""
      )
      readr::write_tsv(x = filters.parameters,
                       path = filters.parameters.path, append = TRUE,
                       col_names = FALSE)
      # update data.info
      data.info <- new.data.info
      if (!is.null(blacklist.markers)) {
        blacklist.markers <- dplyr::bind_rows(blacklist.markers,
                                              blacklist.not.in.common.markers)
      } else {
        blacklist.markers <- blacklist.not.in.common.markers
      }
    }
    blacklist.not.in.common.markers <- NULL
  } #End common.markers

  # Filtering genotyped individuals --------------------------------------------
  if (!is.null(filter.ind.missing.geno) || interactive.filter) {
    ind.approach <- as.character(filter.ind.missing.geno[1])
    ind.threshold <- as.numeric(filter.ind.missing.geno[2])
    prob.pop.threshold <- as.integer(filter.ind.missing.geno[3])

    if (interactive.filter) {
      message("2 steps to visualize and filter the data based on the number of genotyped individuals:")
      message("Step 1. Impact of individual threshold on marker discovery")
      message("Step 2. Choose the filtering approach and thresholds")
    }

    # Folder
    folder.extension <- stringi::stri_join("filter_individual_dart_", file.date, sep = "")
    path.folder.ind.filter <- file.path(path.folder, folder.extension)
    dir.create(path.folder.ind.filter)
    message("Folder created: \n", folder.extension)

    # prepare filter, table and figure
    ind.total <- dplyr::n_distinct(input$INDIVIDUALS) # total number of individuals
    pop.number <- dplyr::n_distinct(input$POP_ID) # number of pop

    # individuals per pop
    ind.pop <- input %>%
      dplyr::distinct(POP_ID, INDIVIDUALS) %>%
      dplyr::group_by(POP_ID) %>%
      dplyr::tally(.) %>%
      dplyr::rename(N_IND = n)

    # input genotyped
    input.genotyped <- input %>%
      dplyr::filter(GT != "000000")

    # overall genotyped individuals
    overall <- input.genotyped %>%
      dplyr::select(MARKERS, INDIVIDUALS) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::tally(.) %>%
      dplyr::rename(GENOTYPED = n) %>%
      dplyr::mutate(PERCENT = ceiling(GENOTYPED/ind.total*100))

    # number of pop. genotyped per marker
    pop.genotyped.marker <- input.genotyped %>%
      dplyr::distinct(MARKERS, POP_ID) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::tally(.) %>%
      dplyr::rename(POP_GENOTYPED = n)

    # genotyped individuals per pop
    pop <- input.genotyped %>%
      dplyr::select(MARKERS, INDIVIDUALS, POP_ID) %>%
      dplyr::group_by(MARKERS, POP_ID) %>%
      dplyr::tally(.) %>%
      dplyr::rename(GENOTYPED = n) %>%
      dplyr::inner_join(ind.pop, by = "POP_ID") %>%
      dplyr::mutate(PERCENT = ceiling(GENOTYPED/N_IND*100))

    input.genotyped <- NULL # unused object

    # Step 1. Impact of individual threshold on marker discovery
    threshold.helper.overall <- overall %>%
      dplyr::ungroup(.) %>%
      dplyr::summarise(
        `10` = length(PERCENT[PERCENT >= 10]),
        `20` = length(PERCENT[PERCENT >= 20]),
        `30` = length(PERCENT[PERCENT >= 30]),
        `40` = length(PERCENT[PERCENT >= 40]),
        `50` = length(PERCENT[PERCENT >= 50]),
        `60` = length(PERCENT[PERCENT >= 60]),
        `70` = length(PERCENT[PERCENT >= 70]),
        `80` = length(PERCENT[PERCENT >= 80]),
        `90` = length(PERCENT[PERCENT >= 90]),
        `100` = length(PERCENT[PERCENT == 100])
      ) %>%
      tidyr::gather(data = ., key = IND_THRESHOLD, value = MARKER_NUMBER) %>%
      dplyr::mutate(POP_ID = rep("OVERALL", n()))

    threshold.helper.pop <- pop %>%
      dplyr::group_by(POP_ID) %>%
      dplyr::summarise(
        `10` = length(PERCENT[PERCENT >= 10]),
        `20` = length(PERCENT[PERCENT >= 20]),
        `30` = length(PERCENT[PERCENT >= 30]),
        `40` = length(PERCENT[PERCENT >= 40]),
        `50` = length(PERCENT[PERCENT >= 50]),
        `60` = length(PERCENT[PERCENT >= 60]),
        `70` = length(PERCENT[PERCENT >= 70]),
        `80` = length(PERCENT[PERCENT >= 80]),
        `90` = length(PERCENT[PERCENT >= 90]),
        `100` = length(PERCENT[PERCENT == 100])
      ) %>%
      tidyr::gather(data = ., key = IND_THRESHOLD, value = MARKER_NUMBER, -POP_ID)

    mean.pop <- threshold.helper.pop %>%
      dplyr::group_by(IND_THRESHOLD) %>%
      dplyr::summarise(
        MARKER_NUMBER = round(mean(MARKER_NUMBER), 0)
      ) %>%
      dplyr::mutate(POP_ID = rep("MEAN_POP", n()))

    # Make sure POP_ID is factor
    if (!is.factor(input$POP_ID)) {
      input$POP_ID <- factor(input$POP_ID)
    }

    threshold.helper <- suppressWarnings(
      dplyr::bind_rows(threshold.helper.pop, mean.pop, threshold.helper.overall) %>%
        dplyr::mutate(
          IND_THRESHOLD = as.numeric(IND_THRESHOLD),
          POP_ID = factor(POP_ID, levels = c(levels(input$POP_ID), "MEAN_POP", "OVERALL"), ordered = TRUE)
        ))

    # Set the breaks for the figure
    max.markers <- dplyr::n_distinct(input$MARKERS)

    #Function to replace plyr::round_any
    rounder <- function(x, accuracy, f = round) {
      f(x / accuracy) * accuracy
    }

    # max.markers <- 658
    if (max.markers >= 1000) {
      y.breaks.by <- rounder(max.markers/10, 100, ceiling)
      y.breaks.max <- rounder(max.markers, 1000, ceiling)
      y.breaks <- seq(0, y.breaks.max, by = y.breaks.by)
    } else {
      y.breaks.by <- rounder(max.markers/10, 10, ceiling)
      y.breaks.max <- rounder(max.markers, 100, ceiling)
      y.breaks <- seq(0, y.breaks.max, by = y.breaks.by)
    }

    plot.ind.threshold <- ggplot2::ggplot(threshold.helper, ggplot2::aes(x = IND_THRESHOLD, y = MARKER_NUMBER)) +
      ggplot2::geom_line() +
      ggplot2::geom_point(size = 2, shape = 21, fill = "white") +
      ggplot2::scale_x_continuous(name = "Individual threshold (percent)", breaks = seq(10, 100, by = 10)) +
      ggplot2::scale_y_continuous(name = "Number of markers", breaks = y.breaks, limits = c(0, y.breaks.max)) +
      ggplot2::theme(
        axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5),
        strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
      ) +
      ggplot2::facet_grid(~POP_ID)
    ggplot2::ggsave(file.path(path.folder.ind.filter, "plot.ind.threshold.pdf"), width = pop.number*4, height = 10, dpi = 600, units = "cm", useDingbats = FALSE)

    if (interactive.filter) {
      message("Step 1. Impact of individual threshold on marker discovery")
      print(plot.ind.threshold)
      message("Look at the plot and inspect the change in the number of markers
in relation to the individual percentage thresholds\n")
    }

    # Helper table for individual thresholds
    ind.threshold.helper.table <- threshold.helper %>%
      dplyr::group_by(POP_ID) %>%
      tidyr::spread(data = ., key = IND_THRESHOLD, MARKER_NUMBER)

    readr::write_tsv(
      x = ind.threshold.helper.table,
      path = file.path(path.folder.ind.filter, "ind.threshold.helper.table.tsv"))

    if (interactive.filter) {
      message("\nInspect the table (ind.threshold.helper.table) to help you view
the relation between individual thresholds and marker discovery:")
      message("First column: POP_ID")
      message("Remaining columns: the individual thresholds in percent
with the value = the number of markers discovered")
      message("The last 2 rows: MEAN_POP is the mean across your populations
and OVERALL is if you had 1 large population")
    }

    # Step 2. Choose the filtering approach and thresholds
    # 2 approach: filtering with the overall n. ind. ("overall") or by pop ("pop")
    if (interactive.filter) {
      message("Step 2. Choose the filtering approach and thresholds")
      message("The approach to filter a marker: do you want it based on the overall
number of genotyped individuals or
on the number of genotyped individuals per pop ? (overall or pop):")
      ind.approach <- as.character(readLines(n = 1))

      message("Enter the individual threshold percentage: ")
      ind.threshold <- as.numeric(readLines(n = 1))
      if (ind.approach == "pop") {
        message("Tolerance for deviation: look at the plot produced ealier and if you see some populations dragging down
                the number of markers for certain percentage thresholds, you have 3 options:\n
                1. remove the population (use pop.select argument to keep the desired populations)
                2. remove individuals with too many missing genotypes (use blacklist.id argument)
                3. use the next threshold (below) to allow variance and then
                manage the missing values with blacklist of individuals and/or
                missing data imputations.\n
                Enter the number of problematic population that you allow to deviate from the threshold:")
        prob.pop.threshold <- as.numeric(readLines(n = 1))
      }
    }
    if (verbose) message("Filtering data")
    # some discrepencies need to be highllighted here. If you have entire pop not genotyped for a markers
    # this will compute them when doing the filtering:
    # summarise(GENOTYPED = length(INDIVIDUALS[GT != "000000"]))
    # if we first remove the individuals not genotyped with :
    # filter(GT != "000000")
    # the pop not genotyped are not accounted for. And this is what we want here.
    # filter_populations take care the ungenotyped pop and common.markers make sure that for certain analysis
    # you can have common markers or not.
    # so here we focus on when pop got a marker is it at 50% 60% 70%  ... genotyped?
    if (ind.approach == "overall") {
      threshold.id <- "(percent)"
      prob.pop.threshold <- "NA"
      ind.approach <- "overall individuals (no pop)"
      filter <- overall %>%
        dplyr::filter(PERCENT >= ind.threshold) %>%
        dplyr::distinct(MARKERS)
    } else {# approach by pop
      threshold.id <- "(percent)"
      ind.approach <- "individuals by pop"
      # ind.threshold <- 0.6
      # prob.pop.threshold <- 3
      filter <- pop %>%
        dplyr::ungroup(.) %>%
        dplyr::filter(PERCENT >= ind.threshold) %>%
        dplyr::group_by(MARKERS) %>%
        dplyr::tally(.) %>%
        dplyr::inner_join(pop.genotyped.marker, by = "MARKERS") %>%
        dplyr::mutate(PROB_POP = POP_GENOTYPED - n) %>%
        dplyr::filter(PROB_POP <= prob.pop.threshold) %>%
        dplyr::distinct(MARKERS)
    }

    # Apply the filter to the tidy data
    input <- dplyr::left_join(x = filter, input, by = "MARKERS")
    new.data.info <- data_info(input) # updating parameters

    n.snp.before <- data.info$n.snp
    if (verbose) message("    Number of markers before = ", n.snp.before)
    n.snp.blacklist <- n.snp.before - new.data.info$n.snp
    if (verbose) message("    Number of markers removed = ", n.snp.blacklist)
    if (verbose) message("    Number of markers after = ", n.snp.before - n.snp.blacklist)

    filters.parameters <- tibble::data_frame(
      FILTERS = "genotyped individuals",
      PARAMETERS = "ind.approach/ind.threshold/prob.pop.threshold",
      VALUES = stringi::stri_join(ind.approach, ind.threshold, prob.pop.threshold, sep = "/"),
      BEFORE = stringi::stri_join(data.info$n.chrom, data.info$n.locus, n.snp.before, sep = "/"),
      AFTER = stringi::stri_join(new.data.info$n.chrom, new.data.info$n.locus, new.data.info$n.snp, sep = "/"),
      BLACKLIST = stringi::stri_join(data.info$n.chrom - new.data.info$n.chrom, data.info$n.locus - new.data.info$n.locus, n.snp.blacklist, sep = "/"),
      UNITS = "CHROM/LOCUS/SNP",
      COMMENTS = ""
    )
    # update data.info
    data.info <- new.data.info
    readr::write_tsv(x = filters.parameters,
                     path = filters.parameters.path, append = TRUE,
                     col_names = FALSE)
    new.whitelist.markers <- dplyr::select(input, dplyr::one_of(want)) %>%
      dplyr::distinct(MARKERS, .keep_all = TRUE)

    blacklist.markers.geno <- dplyr::setdiff(whitelist.markers, new.whitelist.markers)
    if (nrow(blacklist.markers.geno) > 0) {
      readr::write_tsv(x = blacklist.markers.geno, path = file.path(path.folder.ind.filter, "blacklist.markers.geno.tsv"))
      if (!is.null(blacklist.markers)) {
        blacklist.markers <- dplyr::bind_rows(blacklist.markers, blacklist.markers.geno)
      } else {
        blacklist.markers <- blacklist.markers.geno
      }
    }
    whitelist.markers <- new.whitelist.markers
    new.whitelist.markers <- blacklist.markers.geno <- ind.total <- NULL
    pop.number <- ind.pop <- input.genotyped <- overall <- pop.genotyped.marker <- NULL
    pop <- threshold.helper.overall <- threshold.helper.pop <- NULL
    mean.pop <- threshold.helper <- max.markers <- filter <- ind.threshold.helper.table <- NULL
    plot.ind.threshold <- NULL
  }#End filter.ind.missing.geno

  # change to the new directory
  old.dir <- getwd()
  setwd(path.folder)

  # Filter Minor Allele Frequency  ---------------------------------------------
  if (!is.null(maf.thresholds) || interactive.filter) {

    if (interactive.filter) {
      maf.info <- radiator::filter_maf(
        data = input,
        maf.thresholds = c("SNP", 1, "OR", 1, 1),
        parallel.core = parallel.core,
        interactive.filter = TRUE)
    } else {
      maf.info <- radiator::filter_maf(
        data = input,
        maf.thresholds = maf.thresholds,
        parallel.core = parallel.core,
        interactive.filter = FALSE)
    }

    param.delete <- list.files(path = path.folder, pattern = "filters_parameters.tsv", full.names = TRUE)
    suppressWarnings(param.delete <- file.remove(param.delete))
    param.delete <- NULL

    param <- maf.info$filters.parameters

    blacklist.markers.maf <- maf.info$blacklist.markers
    whitelist.markers <- maf.info$whitelist.markers
    input <- maf.info$tidy.filtered.maf
    maf.info <- NULL
    if (nrow(blacklist.markers.maf) > 0) {
      readr::write_tsv(x = param,
                       path = filters.parameters.path, append = TRUE,
                       col_names = FALSE)
      n.snp.before <- data.info$n.snp
      n.snp.blacklist <- nrow(blacklist.markers.maf)
      data.info <- data_info(input)
      if (!is.null(blacklist.markers)) {
        blacklist.markers <- dplyr::bind_rows(blacklist.markers, blacklist.markers.maf)
      } else {
        blacklist.markers <- blacklist.markers.maf
      }
    }
    blacklist.markers.maf <- param <- NULL
    # change to the new directory
    old.dir <- getwd()
    setwd(path.folder)
  }# End maf.thresholds

  # Filter snp number  ---------------------------------------------------------
  setwd(old.dir)
  if (!is.null(number.snp.reads) || interactive.filter) {
    message("Filtering markers based on the number of SNP on the read/locus")
    folder.extension <- stringi::stri_join("filter_snp_number_dart_", file.date, sep = "")
    path.folder.snp.number <- file.path(path.folder, folder.extension)
    dir.create(path.folder.snp.number)
    message("Folder created: \n", folder.extension)

    # get the number of SNP per reads/haplotypes/locus
    number.snp <- dplyr::distinct(input, MARKERS, CHROM, LOCUS, POS) %>%
      dplyr::group_by(LOCUS) %>%
      dplyr::mutate(SNP_N = dplyr::n_distinct(POS))

    number.snp.plot.data <- number.snp %>% dplyr::distinct(LOCUS, SNP_N)

    number.snp.reads.plot <- ggplot2::ggplot(number.snp.plot.data, ggplot2::aes(factor(SNP_N))) +
      ggplot2::geom_bar() +
      ggplot2::labs(x = "Number of SNP per haplotypes (reads)") +
      ggplot2::labs(y = "Distribution (number)") +
      ggplot2::theme(
        axis.title.x = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
        axis.title.y = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
        legend.title = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
        legend.text = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
        strip.text.x = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"))

    # save
    ggplot2::ggsave(
      filename = file.path(path.folder.snp.number, "number.snp.locus.plot.pdf"),
      plot = number.snp.reads.plot,
      width = 20, height = 10, dpi = 600, units = "cm", useDingbats = FALSE)

    if (interactive.filter) {
      print(number.snp.reads.plot)
      message("Based on the plot, choose the threshold in maximum number of SNP per locus allowed (an integer): ")
      number.snp.reads <- as.integer(readLines(n = 1))
    }

    blacklist.snp.number.markers <- number.snp %>%
      dplyr::filter(SNP_N > number.snp.reads) %>%
      dplyr::distinct(LOCUS, .keep_all = TRUE) %>%
      dplyr::select(-SNP_N)

    if (nrow(blacklist.snp.number.markers) > 0) {
      n.snp.before <- data.info$n.snp
      readr::write_tsv(blacklist.snp.number.markers, path = file.path(path.folder.snp.number, "blacklist.snp.number.markers.tsv"))
      input <- dplyr::anti_join(input, blacklist.snp.number.markers, by = "LOCUS")
      whitelist.markers <- dplyr::select(input, dplyr::one_of(c("MARKERS", "CHROM", "LOCUS", "POS"))) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE)
      new.data.info <- data_info(input) # updating parameters

      if (verbose) message("    Number of markers before = ", n.snp.before)
      if (verbose) message("    Number of markers removed = ", data.info$n.snp - new.data.info$n.snp)
      if (verbose) message("    Number of markers after = ", new.data.info$n.snp)

      filters.parameters <- tibble::data_frame(
        FILTERS = "SNP number per reads/locus",
        PARAMETERS = "",
        VALUES = number.snp.reads,
        BEFORE = stringi::stri_join(data.info$n.chrom, data.info$n.locus, data.info$n.snp, sep = "/"),
        AFTER = stringi::stri_join(new.data.info$n.chrom, new.data.info$n.locus, new.data.info$n.snp, sep = "/"),
        BLACKLIST = stringi::stri_join(data.info$n.chrom - new.data.info$n.chrom, data.info$n.locus - new.data.info$n.locus, data.info$n.snp - new.data.info$n.snp, sep = "/"),
        UNITS = "CHROM/LOCUS/SNP",
        COMMENTS = ""
      )
      readr::write_tsv(x = filters.parameters,
                       path = filters.parameters.path, append = TRUE,
                       col_names = FALSE)
      # update data.info
      data.info <- new.data.info
      if (!is.null(blacklist.markers)) {
        blacklist.markers <- dplyr::bind_rows(blacklist.markers, blacklist.snp.number.markers)
      } else {
        blacklist.markers <- blacklist.snp.number.markers
      }
    }
    blacklist.snp.number.markers <- number.snp <- number.snp.reads.plot <- NULL
  }#End number.snp.reads

  # update metadata info
  metadata <- dplyr::filter(metadata, !MARKERS %in% blacklist.markers$MARKERS)
  fst::write.fst(
    x = metadata,
    path = file.path(path.folder, stringi::stri_replace_all_fixed(metadata.file, ".rad", "_filtered.rad", vectorize_all = FALSE)),
    compress =85)

  # readr::write_tsv(
  # x = metadata,
  # path = file.path(path.folder, stringi::stri_replace_all_fixed(metadata.file, ".rad", "_filtered.rad", vectorize_all = FALSE)),
  # )

  # Data quality AFTER filters --------------------------------------------------
  setwd(path.folder)

  # Detect mixed genomes -------------------------------------------------------
  if (interactive.filter || mixed.genomes.analysis) {
    if (verbose) message("Mixed genomes analysis ...")
    mixed.genomes.analysis <- radiator::detect_mixed_genomes(
      data = input,
      ind.heterozygosity.threshold = ind.heterozygosity.threshold)

    blacklist.ind.het <- mixed.genomes.analysis$blacklist.ind.het

    if (interactive.filter) {
      message("\n\nInspect plots and tables in folder created...")
      message("    Do you want to exclude individuals based on heterozygosity ? (y/n): ")
      mixed.gen.analysis <- as.character(readLines(n = 1))
      if (mixed.gen.analysis == "y") {
        message("    Enter the min value for ind.heterozygosity.threshold argument (0 turns off): ")
        threshold.min <- as.numeric(readLines(n = 1))
        message("    Enter the max value for ind.heterozygosity.threshold argument (1 turns off): ")
        threshold.max <- as.numeric(readLines(n = 1))
        blacklist.ind.het  <- dplyr::ungroup(mixed.genomes.analysis$individual.heterozygosity) %>%
          dplyr::filter(HET_PROP > threshold.max | HET_PROP < threshold.min) %>%
          dplyr::distinct(INDIVIDUALS)
        ind.heterozygosity.threshold <- as.numeric(c(threshold.min, threshold.max))
      }
    }


    if (!is.null(nrow(blacklist.ind.het)) &&
        nrow(blacklist.ind.het > 0)) {
      n.ind.blacklisted <- length(blacklist.ind.het$INDIVIDUALS)
      message("Filter individual's heterozygosity: ", n.ind.blacklisted, " individual(s) blacklisted")
      mixed.genome.folder <- list.files(path = path.folder, pattern = "detect_mixed_genomes", full.names = TRUE)

      if (length(mixed.genome.folder) > 1) {
        mixed.genome.folder <- file.info(mixed.genome.folder) %>%
          tibble::rownames_to_column(df = ., var = "FILE") %>%
          dplyr::filter(mtime == max(mtime))
        mixed.genome.folder <- mixed.genome.folder$FILE
      }

      readr::write_tsv(
        x = blacklist.ind.het,
        path = stringi::stri_join(mixed.genome.folder, "/blacklist.ind.het.tsv"),
        col_names = TRUE)

      res$blacklist.id <- res$blacklist.id %>% dplyr::bind_rows(blacklist.ind.het)

      input <- dplyr::anti_join(
        input,
        blacklist.ind.het, by = "INDIVIDUALS")

      # updating parameters
      filters.parameters <- tibble::data_frame(
        FILTERS = "detect mixed genomes",
        PARAMETERS = "ind.heterozygosity.threshold (min/max)",
        VALUES = stringi::stri_join(ind.heterozygosity.threshold, collapse = "/"),
        BEFORE = data.info$n.ind,
        AFTER = data.info$n.ind - n.ind.blacklisted,
        BLACKLIST = n.ind.blacklisted,
        UNITS = "individuals",
        COMMENTS = ""
      )
      readr::write_tsv(x = filters.parameters, path = filters.parameters.path,
                       append = TRUE, col_names = FALSE)
      # update data.info
      data.info$n.ind <- data.info$n.ind - n.ind.blacklisted
    }
    mixed.genomes.analysis <- NULL
  }# End mixed genomes

  # Detect duplicate genomes ---------------------------------------------------
  genome <- duplicate.genomes.analysis[2]
  duplicate.genomes.analysis <- duplicate.genomes.analysis[1]

  if (interactive.filter || duplicate.genomes.analysis) {
    if (verbose) message("Duplicate genomes analysis...")
    duplicate.genomes <- radiator::detect_duplicate_genomes(
      data = input,
      distance.method = "manhattan",
      genome = genome,
      parallel.core = parallel.core)

    if (interactive.filter) {
      message("\n\n    Inspect plots and tables")
      # filtering by distance or pairwise genome similarity ?

      message("    Suspicious about the data?")
      message("    Run the full pairwise genome comparisons")
      message("    This approach integrates markers in common & missing data\n")

      message("    Do you want to run the pairwise genome comparison (y/n): ")
      genome <- as.character(readLines(n = 1))
      if (genome == "y") {
        duplicate.genomes <- radiator::detect_duplicate_genomes(
          data = input,
          distance.method = "manhattan",
          genome = TRUE,
          parallel.core = parallel.core)
      }
      message("    Inspect tables and decide which individual(s) to blacklist(s)")
      message("    Do you want to remove individuals ? (y/n): ")
      remove.id <- as.character(readLines(n = 1))
      if (remove.id == "y") {
        blacklist.id.similar <- tibble::data_frame(INDIVIDUALS = as.character())
        dup.genome.folder <- list.files(
          path = path.folder,
          pattern = "detect_duplicate_genomes", full.names = TRUE)

        if (length(dup.genome.folder) > 1) {
          dup.genome.folder <- file.info(dup.genome.folder) %>%
            tibble::rownames_to_column(df = ., var = "FILE") %>%
            dplyr::filter(mtime == max(mtime))
          dup.genome.folder <- dup.genome.folder$FILE
        }
        blacklist.id.similar.path <- stringi::stri_join(dup.genome.folder, "/blacklist.id.similar.tsv")
        readr::write_tsv(x = blacklist.id.similar, path = blacklist.id.similar.path, append = FALSE, col_names = TRUE)
        message("    An empty blacklist file was generated: \n", blacklist.id.similar.path)
        message("    Keep column name, just add the individual(s) to blacklist(s), save and close the file")
        message("    Type (y) when ready to import blacklist of id(s): ")
        import.blacklist <- as.character(readLines(n = 1))
        if (import.blacklist == "y") {
          blacklist.id.similar <- suppressMessages(readr::read_tsv(blacklist.id.similar.path, col_names = TRUE))
          n.ind.blacklisted <- length(blacklist.id.similar$INDIVIDUALS)
          if (verbose) message("Blacklisted individuals: ", n.ind.blacklisted, " ind.")
          if (verbose) message("    Filtering with blacklist of individuals")
          input <- suppressWarnings(dplyr::anti_join(input, blacklist.id.similar, by = "INDIVIDUALS"))

          res$blacklist.id <- res$blacklist.id %>% dplyr::bind_rows(blacklist.id.similar)
          blacklist.id.similar <- NULL

          # updating parameters
          filters.parameters <- tibble::data_frame(
            FILTERS = "detect duplicate genomes",
            PARAMETERS = "",
            VALUES = "",
            BEFORE = data.info$n.ind,
            AFTER = data.info$n.ind - n.ind.blacklisted,
            BLACKLIST = n.ind.blacklisted,
            UNITS = "individuals",
            COMMENTS = ""
          )
          readr::write_tsv(x = filters.parameters, path = filters.parameters.path,
                           append = TRUE, col_names = FALSE)
          # update data.info
          data.info$n.ind <- data.info$n.ind - n.ind.blacklisted

          # After removing individuals => check for monomorphic markers
          message("\nScan and remove monomorphic markers...")
          input <- discard_monomorphic_markers(data = input, verbose = verbose)
          blacklist.monomorphic.markers <- input$blacklist.monomorphic.markers
          whitelist.markers <- input$whitelist.polymorphic.markers
          input <- input$input

          if (nrow(blacklist.monomorphic.markers) > 0) {
            new.data.info <- data_info(input)

            # updating parameters
            filters.parameters <- tibble::data_frame(
              FILTERS = "removing monomorphic markers",
              PARAMETERS = "",
              VALUES = "",
              BEFORE = stringi::stri_join(data.info$n.chrom, data.info$n.locus, data.info$n.snp, sep = "/"),
              AFTER = stringi::stri_join(new.data.info$n.chrom, new.data.info$n.locus, new.data.info$n.snp, sep = "/"),
              BLACKLIST = stringi::stri_join(data.info$n.chrom - new.data.info$n.chrom, data.info$n.locus - new.data.info$n.locus, data.info$n.snp - new.data.info$n.snp, sep = "/"),
              UNITS = "CHROM/LOCUS/SNP",
              COMMENTS = ""
            )
            readr::write_tsv(x = filters.parameters,
                             path = filters.parameters.path, append = TRUE,
                             col_names = FALSE)
            # update data.info
            data.info <- new.data.info
            if (!is.null(blacklist.markers)) {
              blacklist.markers <- dplyr::bind_rows(blacklist.markers, blacklist.monomorphic.markers)
            } else {
              blacklist.markers <- blacklist.monomorphic.markers
            }
          }
        }
      }
    }
    duplicate.genomes <- blacklist.monomorphic.markers <- NULL
  }#End duplicate.genomes.analysis

  # Missing visualization analysis before filters------------------------------
  # if (missing.analysis) {
  #   if (verbose) message("Missing data analysis: after filters")
  #   missing.visualization <- grur::missing_visualization(data = input, write.plot = TRUE)
  # }

  # Writing to working directory the filtered data frame -----------------------
  # Whitelist
  res$whitelist.markers <- whitelist.markers
  readr::write_tsv(x = res$whitelist.markers, path = "whitelist.markers.tsv", col_names = TRUE)
  if (verbose) message("Writing the whitelist of markers: whitelist.markers.tsv")

  if (nrow(blacklist.markers) > 0) {
    res$blacklist.markers <- blacklist.markers
    readr::write_tsv(x = res$blacklist.markers, path = "blacklist.markers.tsv", col_names = TRUE)
    if (verbose) message("Writing the blacklist of markers: blacklist.markers.tsv")
  }

  # writing the blacklist of id
  if (nrow(res$blacklist.id) > 0) {
    readr::write_tsv(x = res$blacklist.id, path = "blacklist.id.tsv", col_names = TRUE)
    if (verbose) message("Writing the blacklist of ids: blacklist.id.tsv")
  }
  # tidy data
  want <- c("MARKERS", "CHROM", "LOCUS", "POS", "INDIVIDUALS", "POP_ID", "REF", "ALT", "GT", "GT_VCF", "GT_VCF_NUC", "GT_BIN")
  res$tidy.data <- dplyr::select(input, dplyr::one_of(want)) %>%
    dplyr::arrange(MARKERS, POP_ID, INDIVIDUALS)
  input <- NULL

  # write tidy to working directory
  if (!is.null(filename)) {
    tidy.name <- stringi::stri_join(filename, ".filtered.rad")
    fst::write.fst(x = res$tidy.data, path = tidy.name, compress = 85)
    message("Tidy DArT data, filtered, written to folder: \n", tidy.name)
  }


  # Import back the filter parameter file
  res$filters.parameters <- readr::read_tsv(file = filters.parameters.path, col_types = "cccccccc")

  # Generate new strata --------------------------------------------------------
  res$strata <-res$tidy.data %>%
    dplyr::distinct(INDIVIDUALS, POP_ID) %>%
    dplyr::rename(STRATA = POP_ID) %>%
    readr::write_tsv(x = ., path = "new_filtered_strata.tsv")

  # genomic_converter & Imputations --------------------------------------------
  if (!is.null(output)) {
    res$output <- radiator::genomic_converter(
      data = res$tidy.data,
      output = output,
      snp.ld = snp.ld,
      imputation.method = imputation.method,
      hierarchical.levels = hierarchical.levels,
      num.tree = num.tree,
      parallel.core = parallel.core,
      verbose = verbose)
  }

  if (verbose) {
    last.data.info <- data_info(res$tidy.data)
    cat("\n\n\n############################### RESULTS ###############################\n")
    message("DArT data info (before -> after) filters: ")
    message("Number of populations: ", first.data.info$n.pop, " -> ", last.data.info$n.pop)
    message("Number of individuals: ", first.data.info$n.ind, " -> ", last.data.info$n.ind)
    message("Number of chrom: ", first.data.info$n.chrom, " -> ", last.data.info$n.chrom)
    message("Number of locus: ", first.data.info$n.locus, " -> ", last.data.info$n.locus)
    message("Number of SNPs: ", first.data.info$n.snp, " -> ", last.data.info$n.snp)
    timing <- proc.time() - timing
    message("\nComputation time: ", round(timing[[3]]), " sec")
    cat("############################ completed ################################\n")
  }
  setwd(working.dir) #back to the original working directory
  return(res)
}

# update data.info
#' @title data_info
#' @description function generate tidy data main info
#' @rdname data_info
#' @keywords internal
#' @export
data_info <- function(x, print.info = FALSE) {

  if (tibble::has_name(x, "POP_ID")) {
    x.pop.ind <- dplyr::distinct(x, POP_ID, INDIVIDUALS)
    n.pop <- dplyr::n_distinct(x.pop.ind$POP_ID)
    n.ind <- dplyr::n_distinct(x.pop.ind$INDIVIDUALS)
  } else {
    n.pop <- 0
    n.ind <- 0
  }
  if (tibble::has_name(x, "MARKERS")) {
    x <- dplyr::distinct(x, MARKERS, CHROM, LOCUS)
    n.chrom <- dplyr::n_distinct(x$CHROM)
    n.locus <- dplyr::n_distinct(x$LOCUS)
    n.snp <- dplyr::n_distinct(x$MARKERS)
  } else {
    n.chrom <- 0
    n.locus <- 0
    n.snp <- 0
  }
  res <- list(
    n.pop = n.pop,
    n.ind = n.ind,
    n.chrom = n.chrom,
    n.locus = n.locus,
    n.snp = n.snp
  )
  if (print.info) {
    message("Number of populations: ", res$n.pop)
    message("Number of individuals: ", res$n.ind)
    message("Number of chrom: ", res$n.chrom)
    message("Number of locus: ", res$n.locus)
    message("Number of SNPs: ", res$n.snp)
  }
  return(res)
}
