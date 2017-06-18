# Import, filter and transform a dart output file to different formats

#' @name tidy_dart

#' @title Swiss Army knife tool to prepare \href{http://www.diversityarrays.com}{DArT}
#' output file for population genetics analysis.

#' @description Import, filter and generate imputed dataset of DArT output file.

#' @param data DArT output file in \emph{wide} or \emph{binary} format
#' used in CSIRO genomic projects.

#' @inheritParams genomic_converter

#' @param strata A tab delimited file with columns header:
#' \code{INDIVIDUALS} and \code{POP_ID}.
#' Note: the column \code{POP_ID} refers to any grouping of individuals. If a third column
#' named \code{NEW_ID} is used, this column will be used to replace the
#' \code{INDIVIDUALS} column in the main data file.

#' @param duplicate.genomes.analysis (optional) Detect duplicate individuals
#' before and after filters.
#' Default: \code{detect_duplicate_genomes = TRUE}.

#' @param mixed.genomes.analysis (optional) Detect mixed genomes before and
#' after filters.
#' Default: \code{detect_mixed_genomes = TRUE}.
#' @inheritParams detect_mixed_genomes

#' @param missing.analysis (optional) Visualize pattern of missing data before
#' and after filters.
#' Default: \code{missing.analysis = TRUE}

#' @param filter.reproducibility (optional, numerical) Filter the \code{RepAvg}
#' column in the data set. Default: \code{filter.reproducibility = NULL}.
#' e.g to keep markers with reproducibility >= 99%,
#' use: \code{filter.reproducibility = 0.99}.

#' @param plot.reproducibility (optional, logical) Plot the distribution
#' (violin plot) of reproducibility. Default: \code{plot.reproducibility = FALSE}.

#' @param filter.coverage.high (optional, numerical) Filter the upper bound of the
#' \code{AvgCountSnp} column in the data set. Default: \code{filter.coverage.high = NULL}.
#' e.g to keep markers with coverage <= 150 (depth of coverage),
#' use : \code{filter.coverage.high = 150}.

#' @param filter.coverage.low (optional, numerical) Filter the lower bound of the
#' \code{AvgCountSnp} column in the data set. Default: \code{filter.coverage.low = NULL}.
#' e.g to keep markers with coverage >= 10 (depth of coverage),
#' use : \code{filter.coverage.low = 10}.
#' @param plot.coverage (optional, logical) Plot the coverage distribution
#' (violin plot). Default: \code{plot.coverage = FALSE}.

#' @param filter.call.rate (optional, numerical) Filter the \code{CallRate}
#' column in the data set. Default: \code{filter.call.rate = NULL}. e.g to keep
#' markers genotyped in more than 95% of the individuals use :
#' \code{filter.call.rate = 0.95}
#' @param plot.call.rate (optional, logical) Plot the distribution
#' (violin plot) of call rate. Default: \code{plot.call.rate = FALSE}.

#' @param filter.ind.missing.geno (optional, numerical) Filter with the individual's
#' proportion of missing genotype. Similar to call rate.
#' e.g to keep individuals genotyped at >= 0.90 of the markers, use:
#' \code{filter.ind.missing.geno = 0.90}.
#' Default: \code{filter.ind.missing.geno = NULL}.
#' @param plot.ind.missing.geno (optional, logical) Plot the distribution
#' (violin plot) of an individual's genotype proportion.
#' Default: \code{plot.ind.missing.geno = FALSE}.

#' @param filter.markers.missing.ind (optional, numerical) Filter markers based
#' on the proportion of individuals genotyped.
#' e.g to keep markers with >= 0.95 of the genotyped individuals, use:
#' \code{filter.markers.missing.ind = 0.95}.
#' Default: \code{filter.markers.missing.ind = NULL}.
#' @param plot.markers.missing.ind (optional, logical) Plot the distribution
#' (violin plot) of markers genotype individuals proportion.
#' Default: \code{plot.markers.missing.ind = FALSE}.

#' @inheritParams snp_ld

#' @param plot.number.snp.reads (optional, logical) Plot the distribution of SNP
#' per read.
#' Default: \code{plot.number.snp.reads = FALSE}.

#' @inheritParams tidy_genomic_data

#' @param filename (optional) The filename prefix for the objet in the global environment
#' or the working directory. Default: \code{filename = NULL}. A default name will be used,
#' customized with the output file(s) selected.

#' @inheritParams radiator_imputations_module

#' @return The function returns an object (list). The content of the object
#' can be listed with \code{names(object)} and use \code{$} to isolate specific
#' object (see examples). Some output format will write the output file in the
#' working directory. The tidy genomic data frame is generated automatically.

#' @export
#' @rdname tidy_dart
#' @importFrom dplyr group_by select rename filter mutate summarise distinct n_distinct arrange left_join semi_join anti_join inner_join full_join tally bind_rows
#' @importFrom parallel detectCores
#' @importFrom stringi stri_replace_all_fixed stri_join stri_sub stri_replace_na stri_pad_left
#' @importFrom purrr discard
#' @importFrom data.table fread melt.data.table as.data.table
#' @importFrom readr read_tsv write_tsv
#' @importFrom tibble as_data_frame data_frame
#' @importFrom tidyr spread gather unite separate
#importFrom grur missing visualization

#' @examples
#' \dontrun{
#' testing
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

tidy_dart <- function(
  data,
  strata,
  output,
  pop.levels = NULL,
  blacklist.id = NULL,
  pop.select = NULL,
  monomorphic.out = TRUE,
  common.markers = TRUE,
  mixed.genomes.analysis = TRUE,
  ind.heterozygosity.threshold = NULL,
  duplicate.genomes.analysis = TRUE,
  missing.analysis = TRUE,
  filter.reproducibility = NULL,
  plot.reproducibility = FALSE,
  filter.coverage.high = NULL,
  filter.coverage.low = NULL,
  plot.coverage = FALSE,
  filter.call.rate = NULL,
  plot.call.rate = FALSE,
  filter.ind.missing.geno = NULL,
  plot.ind.missing.geno = FALSE,
  filter.markers.missing.ind = NULL,
  plot.markers.missing.ind = FALSE,
  snp.ld = NULL,
  plot.number.snp.reads = FALSE,
  maf.thresholds = NULL,
  maf.pop.num.threshold = 1,
  maf.approach = "SNP",
  maf.operator = "OR",
  filename = NULL,
  imputation.method = NULL,
  hierarchical.levels = "populations",
  verbose = FALSE,
  parallel.core = parallel::detectCores() - 1
) {
  # for timing
  timing <- proc.time()

  if (verbose) {
    cat("#######################################################################\n")
    cat("######################### radiator::tidy_dart ###########################\n")
    cat("#######################################################################\n")
  }
  if (verbose) message("Importing data ...")
  # Checking for missing and/or default arguments ******************************
  if (missing(data)) stop("Input file missing")
  if (missing(strata)) stop("strata file missing")
  if (missing(output)) stop("At least 1 output format is required")

  # Filename -------------------------------------------------------------------
  # Get date and time to have unique filenaming
  if (is.null(filename)) {
    file.date <- stringi::stri_replace_all_fixed(
      Sys.time(),
      pattern = " EDT",
      replacement = "",
      vectorize_all = FALSE
    )
    file.date <- stringi::stri_replace_all_fixed(
      file.date,
      pattern = c("-", " ", ":"),
      replacement = c("", "@", ""),
      vectorize_all = FALSE
    )
    file.date <- stringi::stri_sub(file.date, from = 1, to = 13)

    filename <- stringi::stri_join("radiator_data_", file.date)

    if (!is.null(imputation.method)) {
      filename.imp <- stringi::stri_join("radiator_data_imputed_", file.date)
    }
  } else {
    if (!is.null(imputation.method)) {
      filename.imp <- stringi::stri_join(filename, "_imputed")
    }
  }

  # Strata file ------------------------------------------------------------------
  strata.df <- suppressMessages(readr::read_tsv(file = strata, col_names = TRUE))

  # Import data ---------------------------------------------------------------
  colnames.keeper <- c(c("AlleleID", "SNP", "SnpPosition", "CallRate", "AvgCountRef", "AvgCountSnp", "RepAvg"), strata.df$INDIVIDUALS)

  input <- suppressWarnings(
    data.table::fread(
      input = data,
      sep = "\t",
      stringsAsFactors = FALSE,
      header = TRUE,
      na.strings = "-",
      strip.white = TRUE,
      select = colnames.keeper,
      showProgress = TRUE,
      verbose = FALSE
    ) %>%
      tibble::as_data_frame(.) %>%
      dplyr::rename(LOCUS = AlleleID, POS = SnpPosition, CALL_RATE = CallRate, AVG_COUNT_REF = AvgCountRef, AVG_COUNT_SNP = AvgCountSnp, REP_AVG = RepAvg) %>%
      dplyr::arrange(LOCUS, POS)
  )

  # Screen for duplicate names -------------------------------------------------
  remove.list <- c("LOCUS", "SNP", "POS", "CALL_RATE", "AVG_COUNT_REF", "AVG_COUNT_SNP", "REP_AVG")
  individuals.df <- tibble::data_frame(INDIVIDUALS = purrr::discard(.x = colnames(input), .p = colnames(input) %in% remove.list))
  duplicate.individuals <- length(individuals.df$INDIVIDUALS) - dplyr::n_distinct(individuals.df$INDIVIDUALS)
  if (duplicate.individuals == 0) {
    if (verbose) message("Duplicate individual names in the data: no")
  } else {
    stop(stringi::stri_join("Duplicated individuals names found in the data set.\nNumber of duplicate names = ", duplicate.individuals))
  }
  # removing unused object
  remove.list <- individuals.df <- duplicate.individuals <- NULL

  # Tidying data ---------------------------------------------------------------
  input <- suppressWarnings(
    input %>%
      tidyr::separate(col = LOCUS, into = c("LOCUS", "NOT_USEFUL"), sep = "\\|", extra = "drop") %>%
      dplyr::select(-NOT_USEFUL) %>%
      tidyr::separate(col = SNP, into = c("NOT_USEFUL", "KEEPER"), sep = ":", extra = "drop") %>%
      dplyr::select(-NOT_USEFUL) %>%
      tidyr::separate(col = KEEPER, into = c("REF", "ALT"), sep = ">") %>%
      dplyr::mutate(
        CHROM = rep("CHROM_1", n()),
        MARKERS = stringi::stri_join(CHROM, LOCUS, POS, sep = "__"))
  )

  # Determine the type of DArT file
  binary <- anyDuplicated(input$LOCUS)


  if (binary != 2) {
    if (verbose) message("Tidying the dataset...")
    input <- data.table::melt.data.table(
      data = data.table::as.data.table(input),
      id.vars = c("MARKERS", "CHROM", "LOCUS", "POS", "REF", "ALT", "CALL_RATE", "AVG_COUNT_REF", "AVG_COUNT_SNP", "REP_AVG"),
      variable.name = "INDIVIDUALS",
      variable.factor = FALSE,
      value.name = "GT"
    ) %>%
      tibble::as_data_frame(.)

    n.row <- nrow(input)
    # as.integer is usually twice as light as numeric vector...
    split.vec <- as.integer(floor((parallel.core * 3 * (1:n.row - 1) / n.row) + 1))
    n.row <- NULL

    dart2gt <- function (x) {
      res <- x %>%
        dplyr::mutate(
          REF = stringi::stri_replace_all_fixed(
            str = REF, pattern = c("A", "C", "G", "T"),
            replacement = c("001", "002", "003", "004"), vectorize_all = FALSE),
          ALT = stringi::stri_replace_all_fixed(
            str = ALT, pattern = c("A", "C", "G", "T"),
            replacement = c("001", "002", "003", "004"), vectorize_all = FALSE),
          GT = as.character(GT),
          GT = stringi::stri_replace_all_fixed(
            str = GT, pattern = c("0", "1", "2"),
            replacement = c("REF_REF", "ALT_ALT", "REF_ALT"),
            vectorize_all = FALSE),
          GT = stringi::stri_replace_na(str= GT, replacement = "000_000"),
          GENOTYPE = dplyr::if_else(
            GT == "REF_REF", stringi::stri_join(REF, REF, sep = ""),
            dplyr::if_else(GT == "ALT_ALT", stringi::stri_join(ALT, ALT, sep = ""),
                           dplyr::if_else(GT == "REF_ALT", stringi::stri_join(REF, ALT, sep = ""),
                                          dplyr::if_else(GT == "ALT_REF",  stringi::stri_join(ALT, REF, sep = ""),"000000"))))
        ) %>%
        dplyr::select(-GT, GT = GENOTYPE)
      return(res)
    }

    input <- dplyr::bind_cols(
      dplyr::select(input, -c(REF, ALT, GT)),
      dplyr::select(input, GT, REF, ALT) %>%
        split(x = ., f = split.vec) %>%
        .radiator_parallel(
          X = ., FUN = dart2gt, mc.cores = parallel.core) %>%
        dplyr::bind_rows(.))
    split.vec <- NULL
  }
  if (binary == 2) {
    if (verbose) message("Tidying DArT binary data set")
    # necessary to deal with the duplication of lines because of the GT in 2 lines
    grouping.column <- input %>%
      ungroup() %>%
      dplyr::select(MARKERS, CHROM, LOCUS, POS, REF, ALT, CALL_RATE, AVG_COUNT_REF,
                    AVG_COUNT_SNP, REP_AVG) %>%
      dplyr::filter(!is.na(REF) | !is.na(ALT)) %>%
      dplyr::distinct(MARKERS, CHROM, LOCUS, POS, REF, ALT, CALL_RATE, AVG_COUNT_REF,
                      AVG_COUNT_SNP, REP_AVG, .keep_all = TRUE) %>%
      dplyr::mutate(
        REF = stringi::stri_replace_all_fixed(
          str = REF, pattern = c("A", "C", "G", "T"),
          replacement = c("01", "02", "03", "04"),
          vectorize_all = FALSE), # replace nucleotide with numbers
        ALT = stringi::stri_replace_all_fixed(
          str = ALT, pattern = c("A", "C", "G", "T"),
          replacement = c("01", "02", "03", "04"),
          vectorize_all = FALSE)# replace nucleotide with numbers
      )

    # Data tidying
    input <- input %>%
      dplyr::select(-c(CHROM, LOCUS, POS, REF, ALT, CALL_RATE, AVG_COUNT_REF,
                       AVG_COUNT_SNP, REP_AVG)) %>%
      dplyr::mutate(
        ALLELE_NAME = rep("A", n()),
        ALLELE_NUMBER = rep(1:2, each = 1, times = n()/2)
      ) %>%
      tidyr::unite(ALLELE, c(ALLELE_NAME, ALLELE_NUMBER), sep = "") %>%
      tidyr::gather(INDIVIDUALS, GENOTYPE, -c(MARKERS, ALLELE)) %>%
      tidyr::spread(data = ., ALLELE, GENOTYPE) %>%
      tidyr::unite(GENOTYPE, c(A1, A2)) %>%
      dplyr::inner_join(grouping.column, by = c("MARKERS")) %>%
      dplyr::mutate(
        GT = ifelse(GENOTYPE == "1_0", "REF_REF",
                    ifelse(GENOTYPE == "0_1", "ALT_ALT",
                           ifelse(GENOTYPE == "1_1", "REF_ALT", "0_0")))
      ) %>%
      dplyr::select(MARKERS, CHROM, LOCUS, POS, REF, ALT, INDIVIDUALS, GT, CALL_RATE,
                    AVG_COUNT_REF, AVG_COUNT_SNP, REP_AVG)

    grouping.column <- NULL # remove unused object
  }

  # get the number of markers before filters
  locus.before.filters <- dplyr::n_distinct(input$LOCUS)
  snp.before.filters <- dplyr::n_distinct(input$MARKERS)

  # Strata file ----------------------------------------------------------------
  input <- dplyr::left_join(input, strata.df, by = "INDIVIDUALS")

  if (ncol(strata.df) == 3) {
    input <- input %>%
      dplyr::select(-INDIVIDUALS) %>%
      dplyr::rename(INDIVIDUALS = NEW_ID)
  }

  if (tibble::has_name(input, "STRATA")) {
    input <- dplyr::rename(input, POP_ID = STRATA)
    # colnames(input) <- stringi::stri_replace_all_fixed(str = colnames(input), pattern = "STRATA", replacement = "POP_ID", vectorize_all = FALSE)
  }

  # pop.levels -------------------------------------------------------------------
  input <- change_pop_names(data = input, pop.levels = pop.levels)

  # Prepare to store results in list -------------------------------------------
  res <- list()

  # Import blacklist id --------------------------------------------------------
  if (is.null(blacklist.id)) { # No blacklist of ID
    if (verbose) message("Blacklisted individuals: no")
  } else {# With blacklist of ID
    blacklist.id <- suppressMessages(readr::read_tsv(blacklist.id, col_names = TRUE))
    if (verbose) message("Blacklisted individuals: yes (", length(blacklist.id$INDIVIDUALS), " ind.)")
    if (verbose) message("    Filtering with blacklist of individuals")
    input <- suppressWarnings(dplyr::anti_join(input, blacklist.id, by = "INDIVIDUALS"))
    res$blacklist.id <- blacklist.id
    blacklist.id <- NULL
  }

  # pop.select -----------------------------------------------------------------
  if (!is.null(pop.select)) {
    if (verbose) message(stringi::stri_join(length(pop.select), "population(s) selected", sep = " "))
    input <- suppressWarnings(input %>% dplyr::filter(POP_ID %in% pop.select))
  }
  input$POP_ID <- droplevels(input$POP_ID)
  pop.levels <- levels(input$POP_ID)

  # Filter monomorphic markers  ------------------------------------------------
  if (monomorphic.out) {
    mono <- discard_monomorphic_markers(data = input, verbose = verbose)

    input <- mono$input
    blacklist.monomorphic <- mono$blacklist.momorphic.markers
    mono <- NULL

    if (nrow(blacklist.monomorphic) > 0) {
      res$blacklist.monomorphic <- blacklist.monomorphic
      readr::write_tsv(x = blacklist.monomorphic, path = "blacklist.monomorphic.tsv", col_names = TRUE)
    } else {
      res$blacklist.monomorphic <- "no monomorphic markers blacklisted"
    }
    blacklist.monomorphic <- NULL
  }

  # Filter common markers between all populations  ------------------------------
  if (common.markers) {
    input <- keep_common_markers(data = input, plot = TRUE,
                                 verbose = verbose)
  }

  # Detect mixed genomes -------------------------------------------------------
  if (mixed.genomes.analysis) {
    if (verbose) {
      message("Mixed genomes analysis before filters...")
    }
    res$mixed.genomes.before.filters <- detect_mixed_genomes(
        data = input,
        ind.heterozygosity.threshold = ind.heterozygosity.threshold)

    if (!is.null(nrow(res$mixed.genomes.before.filters$blacklist.ind.het)) && nrow(res$mixed.genomes.before.filters$blacklist.ind.het > 0)) {
      readr::write_tsv(x = res$mixed.genomes.before.filters$blacklist.ind.het, path = "blacklist.individuals.heterozygosity.tsv", col_names = TRUE)
      input <- dplyr::anti_join(input, res$mixed.genomes.before.filters$blacklist.ind.het, by = "INDIVIDUALS")
    }
  }
  # Detect duplicate genomes ---------------------------------------------------
  if (duplicate.genomes.analysis) {
    if (verbose) message("Duplicate genomes analysis before filters...")
    res$duplicate.genomes.before.filters <- detect_duplicate_genomes(
      data = input, distance.method = "manhattan", parallel.core = parallel.core)
  }

  # Missing visualization analysis before filters-------------------------------
  # if (missing.analysis) {
  #   if (verbose) message("Missing data analysis: before filters")
  #   res$missing.before.filters <- grur::missing_visualization(data = input)
  #   res$missing.before.filters$tidy.data <- NULL
  #   res$missing.before.filters$tidy.data.binary <- NULL
  # }

  # Filtering reproducibility  -------------------------------------------------
  if (!is.null(filter.reproducibility)) {
    filter <- input %>%
      dplyr::filter(REP_AVG >= filter.reproducibility)

    whitelist.filter <- filter %>%
      dplyr::select(MARKERS, LOCUS, POS) %>%
      dplyr::distinct(MARKERS, LOCUS, POS)

    blacklist.markers.reproducibility <- input %>%
      dplyr::distinct(MARKERS, LOCUS, POS) %>%
      dplyr::filter(!MARKERS %in% whitelist.filter$MARKERS)

    if (verbose) message("Filter reproducibility: ", dplyr::n_distinct(input$MARKERS) - dplyr::n_distinct(filter$MARKERS), " markers deleted")
    if (length(blacklist.markers.reproducibility$MARKERS > 0)) {
      readr::write_tsv(
        x = blacklist.markers.reproducibility,
        path = "blacklist.markers.reproducibility.tsv", col_names = TRUE)
      res$blacklist.markers.reproducibility <- blacklist.markers.reproducibility
    }
    blacklist.markers.reproducibility <- NULL

    if (plot.reproducibility) {
      data.combined <- dplyr::bind_rows(
        data.before <- input %>%
          dplyr::select(POP_ID, INDIVIDUALS, REP_AVG) %>%
          dplyr::mutate(GROUP = rep("before", n())),
        data.after <- filter %>%
          dplyr::select(POP_ID, INDIVIDUALS, REP_AVG) %>%
          dplyr::mutate(GROUP = rep("after", n()))
      ) %>%
        dplyr::mutate(GROUP = factor(GROUP, levels = c("before", "after"), ordered = TRUE))

      res$plot.reproducibility <- ggplot2::ggplot(data.combined, ggplot2::aes(x = factor(POP_ID), y = REP_AVG, na.rm = TRUE)) +
        ggplot2::geom_violin(trim = TRUE) +
        ggplot2::geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA) +
        ggplot2::stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white") +
        ggplot2::labs(x = "Sampling sites") +
        ggplot2::labs(y = "Markers reproducibility average") +
        ggplot2::theme(
          legend.position = "none",
          axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
          axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
          axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5),
          legend.title = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
          legend.text = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
          strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
        ) +
        ggplot2::facet_grid(~GROUP)
    }

    if (!(plot.reproducibility)) {
      res$plot.reproducibility <- "not selected"
    }
    input <- filter
  }

  if (is.null(filter.reproducibility) & (plot.reproducibility)) {
    data.combined <- input %>%
      dplyr::select(POP_ID, INDIVIDUALS, REP_AVG) %>%
      dplyr::mutate(GROUP = rep("before filter", n()))

    res$plot.reproducibility <- ggplot2::ggplot(data.combined, ggplot2::aes(x = factor(POP_ID), y = REP_AVG, na.rm = TRUE)) +
      ggplot2::geom_violin(trim = TRUE) +
      ggplot2::geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA) +
      ggplot2::stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white") +
      ggplot2::labs(x = "Sampling sites") +
      ggplot2::labs(y = "Markers reproducibility average before filter") +
      ggplot2::theme(
        legend.position = "none",
        axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5),
        legend.title = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        legend.text = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
      ) +
      ggplot2::facet_grid(~GROUP)
  }

  filter <- whitelist.filter <- data.combined <- data.before <- data.after <- NULL

  # Filtering coverage --------------------------------------------------------
  input.before.filter <- input

  # high bound
  if (!is.null(filter.coverage.high)) {
    filter <- input %>%
      dplyr::filter(AVG_COUNT_SNP <= filter.coverage.high)

    whitelist.filter <- filter %>% dplyr::distinct(MARKERS, LOCUS, POS)

    blacklist.markers.coverage.high <- input %>%
      dplyr::distinct(MARKERS, LOCUS, POS) %>%
      dplyr::filter(!MARKERS %in% whitelist.filter$MARKERS)

    if (verbose) message("Filter coverage high: ", dplyr::n_distinct(input$MARKERS) - dplyr::n_distinct(filter$MARKERS), " markers deleted")
    if (length(blacklist.markers.coverage.high$MARKERS > 0)) {
      readr::write_tsv(x = blacklist.markers.coverage.high, path = "blacklist.markers.coverage.high.tsv", col_names = TRUE)
      input <- filter
      res$blacklist.markers.coverage.high <- blacklist.markers.coverage.high
    }
  }
  # lower bound
  if (!is.null(filter.coverage.low)) {
    filter <- input %>%
      dplyr::filter(AVG_COUNT_SNP >= filter.coverage.low)

    whitelist.filter <- filter %>% dplyr::distinct(MARKERS, LOCUS, POS)

    blacklist.markers.coverage.low <- input %>%
      dplyr::distinct(MARKERS, LOCUS, POS) %>%
      dplyr::filter(!MARKERS %in% whitelist.filter$MARKERS)

    if (verbose) message("Filter coverage low: ", dplyr::n_distinct(input$MARKERS) - dplyr::n_distinct(filter$MARKERS), " markers deleted")
    if (length(blacklist.markers.coverage.low$MARKERS > 0)) {
      readr::write_tsv(x = blacklist.markers.coverage.low, path = "blacklist.markers.coverage.low.tsv", col_names = TRUE)
      input <- filter
      res$blacklist.markers.coverage.low <- blacklist.markers.coverage.low
    }
  }

  if (!is.null(filter.coverage.high) | !is.null(filter.coverage.low) & (plot.coverage)) {
    data.combined <- dplyr::bind_rows(
      data.before <- input.before.filter %>%
        dplyr::select(POP_ID, INDIVIDUALS, AVG_COUNT_SNP) %>%
        dplyr::mutate(GROUP = rep("before", n())),
      data.after <- input %>%
        dplyr::select(POP_ID, INDIVIDUALS, AVG_COUNT_SNP) %>%
        dplyr::mutate(GROUP = rep("after", n()))
    ) %>%
      dplyr::mutate(GROUP = factor(GROUP, levels = c("before", "after"), ordered = TRUE))

    res$plot.coverage <- ggplot2::ggplot(data.combined, ggplot2::aes(x = factor(POP_ID), y = AVG_COUNT_SNP, na.rm = TRUE)) +
      ggplot2::geom_violin(trim = TRUE) +
      ggplot2::geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA) +
      ggplot2::stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white") +
      ggplot2::labs(x = "Sampling sites") +
      ggplot2::labs(y = "Markers coverage") +
      ggplot2::theme(
        legend.position = "none",
        axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5),
        legend.title = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        legend.text = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
      ) +
      ggplot2::facet_grid(~GROUP)
  }
  if (is.null(filter.coverage.high) & is.null(filter.coverage.low) & (plot.coverage)) {
    data.combined <- input.before.filter %>%
      dplyr::select(POP_ID, INDIVIDUALS, AVG_COUNT_SNP) %>%
      dplyr::mutate(GROUP = rep("before filter", n()))

    res$plot.coverage <- ggplot2::ggplot(data.combined, ggplot2::aes(x = factor(POP_ID), y = AVG_COUNT_SNP, na.rm = TRUE)) +
      ggplot2::geom_violin(trim = TRUE) +
      ggplot2::geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA) +
      ggplot2::stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white") +
      ggplot2::labs(x = "Sampling sites") +
      ggplot2::labs(y = "Markers coverage before filter") +
      ggplot2::theme(
        legend.position = "none",
        axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5),
        legend.title = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        legend.text = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
      ) +
      ggplot2::facet_grid(~GROUP)
  }
  if (!(plot.coverage)) {
    res$plot.coverage <- "not selected"
  }
  filter <- whitelist.filter <- data.combined <- data.before <- data.after <- input.before.filter <- NULL
  blacklist.markers.coverage.low <- blacklist.markers.coverage.high <- NULL

  # Filtering call rate ---------------------------------------------------------
  if (!is.null(filter.call.rate)) {
    filter <- input %>%
      dplyr::filter(CALL_RATE >= filter.call.rate)

    whitelist.filter <- filter %>% dplyr::distinct(MARKERS, LOCUS, POS)

    blacklist.call.rate <- input %>%
      dplyr::distinct(MARKERS, LOCUS, POS) %>%
      dplyr::filter(!MARKERS %in% whitelist.filter$MARKERS)

    if (verbose) message("Filter call rate: ", dplyr::n_distinct(input$MARKERS) - dplyr::n_distinct(filter$MARKERS), " markers deleted")
    if (length(blacklist.call.rate$MARKERS > 0)) {
      readr::write_tsv(x = blacklist.call.rate, path = "blacklist.call.rate.tsv", col_names = TRUE)
      res$blacklist.call.rate <- blacklist.call.rate
    }
    if (plot.call.rate) {
      data.combined <- dplyr::bind_rows(
        data.before <- input %>%
          dplyr::select(POP_ID, INDIVIDUALS, CALL_RATE) %>%
          dplyr::mutate(GROUP = rep("before", n())),
        data.after <- filter %>%
          dplyr::select(POP_ID, INDIVIDUALS, CALL_RATE) %>%
          dplyr::mutate(GROUP = rep("after", n()))
      ) %>%
        dplyr::mutate(GROUP = factor(GROUP, levels = c("before", "after"), ordered = TRUE))

      res$plot.call.rate <- ggplot2::ggplot(data.combined, ggplot2::aes(x = factor(POP_ID), y = CALL_RATE, na.rm = TRUE)) +
        ggplot2::geom_violin(trim = TRUE) +
        ggplot2::geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA) +
        ggplot2::stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white") +
        ggplot2::labs(x = "Sampling sites") +
        ggplot2::labs(y = "Markers call rate") +
        ggplot2::theme(
          legend.position = "none",
          axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
          axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
          axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5),
          legend.title = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
          legend.text = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
          strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
        ) +
        ggplot2::facet_grid(~GROUP)

    }
    if (!(plot.call.rate)) {
      res$plot.call.rate <- "not selected"
    }
    input <- filter
  }
  if (is.null(filter.call.rate) & (plot.call.rate)) {
    data.combined <- input %>%
      dplyr::select(POP_ID, INDIVIDUALS, CALL_RATE) %>%
      dplyr::mutate(GROUP = rep("before filter", n()))

    res$plot.call.rate <- ggplot2::ggplot(data.combined, ggplot2::aes(x = factor(POP_ID), y = CALL_RATE, na.rm = TRUE)) +
      ggplot2::geom_violin(trim = TRUE) +
      ggplot2::geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA) +
      ggplot2::stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white") +
      ggplot2::labs(x = "Sampling sites") +
      ggplot2::labs(y = "Markers call rate before filter") +
      ggplot2::theme(
        legend.position = "none",
        axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5),
        legend.title = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        legend.text = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
      ) +
      ggplot2::facet_grid(~GROUP)
  }

  filter <- whitelist.filter <- data.combined <- data.before <- data.after <- input.before.filter <- NULL
  blacklist.call.rate <- NULL

  # Filtering genotyped individuals --------------------------------------------
  # filter.ind.missing.geno = NULL,
  # plot.ind.missing.geno = FALSE,

  if (!is.null(filter.ind.missing.geno)) {
    # filter.ind.missing.geno <- 0.90 # test

    # Create a new df with genotyped prop.
    ind.geno.prop <- input %>%
      dplyr::group_by(INDIVIDUALS) %>%
      dplyr::summarise(
        GENOTYPED_PROP = length(GT[GT != "000000"])/length(GT)
      )

    # merge with dataset
    input <- input %>%
      dplyr::full_join(ind.geno.prop, by = "INDIVIDUALS")

    # filter
    filter <- input %>%
      dplyr::filter(GENOTYPED_PROP >= filter.ind.missing.geno)


    whitelist.filter <- filter %>% dplyr::distinct(INDIVIDUALS)

    blacklist.ind.missing.geno <- input %>%
      dplyr::distinct(INDIVIDUALS) %>%
      dplyr::filter(!INDIVIDUALS %in% whitelist.filter$INDIVIDUALS)

    if (verbose) message("Filter individuals with less than ", filter.ind.missing.geno, " missing genotype prop: ", dplyr::n_distinct(input$INDIVIDUALS) - dplyr::n_distinct(filter$INDIVIDUALS), " individuals removed")
    if (length(blacklist.ind.missing.geno$INDIVIDUALS > 0)) {
      readr::write_tsv(x = blacklist.ind.missing.geno, path = "blacklist.ind.missing.geno.tsv", col_names = TRUE)
      res$blacklist.ind.missing.geno <- blacklist.ind.missing.geno
    }

    if (plot.ind.missing.geno) {
      data.combined <- dplyr::bind_rows(
        data.before <- input %>%
          dplyr::select(POP_ID, INDIVIDUALS, GENOTYPED_PROP) %>%
          dplyr::mutate(GROUP = rep("before", n())),
        data.after <- filter %>%
          dplyr::select(POP_ID, INDIVIDUALS, GENOTYPED_PROP) %>%
          dplyr::mutate(GROUP = rep("after", n()))
      ) %>%
        dplyr::mutate(GROUP = factor(GROUP, levels = c("before", "after"), ordered = TRUE))

      res$plot.ind.missing.geno <- ggplot2::ggplot(data.combined, ggplot2::aes(x = factor(POP_ID), y = GENOTYPED_PROP, na.rm = TRUE)) +
        ggplot2::geom_violin(trim = TRUE) +
        ggplot2::geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA) +
        ggplot2::stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white") +
        ggplot2::labs(x = "Sampling sites") +
        ggplot2::labs(y = "Individual's genotyped proportion") +
        ggplot2::theme(
          legend.position = "none",
          axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
          axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
          axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5),
          legend.title = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
          legend.text = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
          strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
        ) +
        ggplot2::facet_grid(~GROUP)
    }

    if(!(plot.ind.missing.geno)) {
      res$plot.ind.missing.geno <- "not selected"
    }
    input <- filter
  }

  if (is.null(filter.ind.missing.geno) & (plot.ind.missing.geno)) {
    # Create a new df with genotyped prop.
    ind.geno.prop <- input %>%
      dplyr::group_by(INDIVIDUALS) %>%
      dplyr::summarise(
        GENOTYPED_PROP = length(GT[GT != "000000"])/length(GT)
      )

    # merge with dataset
    input <- input %>%
      dplyr::full_join(ind.geno.prop, by = "INDIVIDUALS")

    data.combined <- input %>%
      dplyr::select(POP_ID, INDIVIDUALS, GENOTYPED_PROP) %>%
      dplyr::mutate(GROUP = rep("before filter", n()))

    res$plot.ind.missing.geno <- ggplot2::ggplot(data.combined, ggplot2::aes(x = factor(POP_ID), y = GENOTYPED_PROP, na.rm = TRUE)) +
      ggplot2::geom_violin(trim = TRUE) +
      ggplot2::geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA) +
      ggplot2::stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white") +
      ggplot2::labs(x = "Sampling sites") +
      ggplot2::labs(y = "Individual's genotyped proportion") +
      ggplot2::theme(
        legend.position = "none",
        axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5),
        legend.title = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        legend.text = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
      ) +
      ggplot2::facet_grid(~GROUP)
  }

  filter <- whitelist.filter <- data.combined <- data.before <- data.after <- input.before.filter <- NULL
  blacklist.ind.missing.geno <- ind.geno.prop <- NULL

  # Filtering filter.markers.missing.ind  --------------------------------------
  # filter.markers.missing.ind = NULL,
  # plot.markers.missing.ind = FALSE,
  if (!is.null(filter.markers.missing.ind)) {
    # filter.markers.missing.ind <- 0.95 # test

    # Create a new df with marker.missing.ind prop.
    marker.missing.ind.prop <- input %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::summarise(
        MISSING_IND_PROP = length(GT[GT != "000000"])/length(GT)
      )

    # merge with dataset
    input <- input %>%
      dplyr::full_join(marker.missing.ind.prop, by = "MARKERS")

    # filter
    filter <- input %>%
      dplyr::filter(MISSING_IND_PROP >= filter.markers.missing.ind)

    whitelist.filter <- filter %>% dplyr::distinct(MARKERS)

    blacklist.marker.missing.ind.prop <- input %>%
      dplyr::distinct(MARKERS) %>%
      dplyr::filter(!MARKERS %in% whitelist.filter$MARKERS)

    if (verbose) message("Filter markers with less than ", filter.markers.missing.ind, " missing genotype ind: ", dplyr::n_distinct(input$MARKERS) - dplyr::n_distinct(filter$MARKERS), " markers removed")
    if (length(blacklist.marker.missing.ind.prop$MARKERS > 0)) {
      readr::write_tsv(x = blacklist.marker.missing.ind.prop, path = "blacklist.marker.missing.ind.tsv", col_names = TRUE)
      res$blacklist.marker.missing.ind.prop <- blacklist.marker.missing.ind.prop
    }

    if (plot.markers.missing.ind) {
      data.combined <- dplyr::bind_rows(
        data.before <- input %>%
          dplyr::select(POP_ID, MARKERS, MISSING_IND_PROP) %>%
          dplyr::mutate(GROUP = rep("before", n())),
        data.after <- filter %>%
          dplyr::select(POP_ID, MARKERS, MISSING_IND_PROP) %>%
          dplyr::mutate(GROUP = rep("after", n()))
      ) %>%
        dplyr::mutate(GROUP = factor(GROUP, levels = c("before", "after"), ordered = TRUE))

      res$plot.markers.missing.ind <- ggplot2::ggplot(data.combined, ggplot2::aes(x = factor(POP_ID), y = MISSING_IND_PROP, na.rm = TRUE)) +
        ggplot2::geom_violin(trim = TRUE) +
        ggplot2::geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA) +
        ggplot2::stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white") +
        ggplot2::labs(x = "Sampling sites") +
        ggplot2::labs(y = "Marker's genotyped proportion") +
        ggplot2::theme(
          legend.position = "none",
          axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
          axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
          axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5),
          legend.title = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
          legend.text = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
          strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
        ) +
        ggplot2::facet_grid(~GROUP)
    }

    if(!(plot.markers.missing.ind)) {
      res$plot.markers.missing.ind <- "not selected"
    }
    input <- filter
  }

  if (is.null(filter.markers.missing.ind) & (plot.markers.missing.ind)) {
    data.combined <- input %>%
      dplyr::select(POP_ID, MARKERS, MISSING_IND_PROP) %>%
      dplyr::mutate(GROUP = rep("before filter", n()))

    res$plot.markers.missing.ind <- ggplot2::ggplot(data.combined, ggplot2::aes(x = factor(POP_ID), y = MISSING_IND_PROP, na.rm = TRUE)) +
      ggplot2::geom_violin(trim = TRUE) +
      ggplot2::geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA) +
      ggplot2::stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white") +
      ggplot2::labs(x = "Sampling sites") +
      ggplot2::labs(y = "Marker's genotyped proportion") +
      ggplot2::theme(
        legend.position = "none",
        axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5),
        legend.title = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        legend.text = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
      ) +
      ggplot2::facet_grid(~GROUP)
  }

  filter <- whitelist.filter <- data.combined <- data.before <- data.after <- input.before.filter <- NULL
  blacklist.marker.missing.ind.prop <- marker.missing.ind.prop <- NULL

  # snp.ld  --------------------------------------------------------------------

  if (plot.number.snp.reads) {
    # get the number of SNP per reads/haplotypes/locus
    res$number.snp.reads <- input %>%
      dplyr::group_by (LOCUS) %>%
      dplyr::summarise(SNP_N = dplyr::n_distinct(POS))

    res$number.snp.reads.plot <- ggplot2::ggplot(res$number.snp.reads, ggplot2::aes(factor(SNP_N))) +
      ggplot2::geom_bar() +
      ggplot2::labs(x="Number of SNP per haplotypes (reads)") +
      ggplot2::labs(y="Distribution (number)") +
      ggplot2::theme(axis.title.x = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
                     axis.title.y = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
                     legend.title = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
                     legend.text = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
                     strip.text.x = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"))
  }

  if (!(plot.number.snp.reads)) {
    res$number.snp.reads.plot <- "not selected"
  }

  # filter snp per reads
  if (!is.null(snp.ld)) {
    input <- snp_ld(data = input, snp.ld = snp.ld)
  }

  # Filter Minor Allele Frequency  ---------------------------------------------
  if (!is.null(maf.thresholds)) {
    maf.info <- radiator_maf_module(
      data = input,
      maf.thresholds = maf.thresholds,
      maf.pop.num.threshold = maf.pop.num.threshold,
      maf.approach = maf.approach,
      maf.operator = maf.operator,
      parallel.core = parallel.core
    )

    input <- maf.info$input
    res$maf.data <- maf.info$maf.data
    maf.info <- NULL
  } # End of MAF filters

  # Data quality AFTER filter --------------------------------------------------

  # Detect mixed genomes
  if (mixed.genomes.analysis) {
    if (verbose) message("Mixed genomes analysis after filters...")
    res$mixed.genomes.after.filters <- detect_mixed_genomes(
      data = input, ind.heterozygosity.threshold = NULL)
  }

  # Detect duplicate genomes
  if (duplicate.genomes.analysis) {
    if (verbose) message("Duplicate genomes analysis after filters...")
    res$duplicate.genomes.after.filters <- detect_duplicate_genomes(
      data = input,
      distance.method = "manhattan",
      # genome = TRUE,
      parallel.core = parallel.core)
  }

  # Missing visualization analysis before filters
  # if (missing.analysis) {
  #   if (verbose) message("Missing data analysis: after filters")
  #   res$missing.after.filters <- grur::missing_visualization(data = input)
  #   res$missing.after.filters$tidy.data <- NULL
  #   res$missing.after.filters$tidy.data.binary <- NULL
  # }


  # Writing to working directory the filtered data frame -----------------------

  # Whitelist
  res$whitelist.markers <- dplyr::distinct(input, MARKERS, CHROM, LOCUS, POS)
  readr::write_tsv(x = res$whitelist.markers, path = "whitelist.markers.tsv", col_names = TRUE)
  if (verbose) message("Writing the whitelist of markers: whitelist.markers.tsv")

  # filtered tidy data
  res$tidy.data <- input %>%
    dplyr::arrange(MARKERS, POP_ID, INDIVIDUALS)
  input <- NULL

  # write to working directory
  filtered.data.name <- stringi::stri_join(filename, "_tidy.tsv", sep = "")
  if (verbose) message("Writing the tidy and filtered data set: ",
                       filtered.data.name, "\nWorking directory: ", getwd())
  readr::write_tsv(x = res$tidy.data, path = filtered.data.name, col_names = TRUE)

  # genomic_converter & Imputations --------------------------------------------
  res$output <- genomic_converter(data = res$tidy.data,
                              output = output,
                              imputation.method = imputation.method,
                              hierarchical.levels = hierarchical.levels,
                              parallel.core = parallel.core,
                              verbose = verbose)

  # Results --------------------------------------------------------------------
  snp.after.filters <- dplyr::n_distinct(res$tidy.data$MARKERS)
  locus.after.filters <- dplyr::n_distinct(res$tidy.data$LOCUS)

  if (verbose) {
    cat("############################### RESULTS ###############################\n")
    message("The number of markers removed by the filters:\nSNP: ",
            snp.before.filters - snp.after.filters, "\nLOCUS: ",
            locus.before.filters - locus.after.filters)
    message("The number of markers before -> after the filters")
    message("SNP: ", snp.before.filters, " -> ", snp.after.filters)
    message("LOCUS: ", locus.before.filters, " -> ", locus.after.filters)
    timing <- proc.time() - timing
    message("\nComputation time: ", round(timing[[3]]), " sec")
    cat("############################ completed ################################\n")
  }
  return(res)
}
