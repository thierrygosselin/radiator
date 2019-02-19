## VCF
#' @title Summary statistics of a tidy VCF by population and markers
#' @description Summarise and prepare the tidy VCF.
#' Summary, by population and markers (SNP), of frequency of the REF
#' and the ALT alleles, the observed and the expected heterozygosity
#' and the inbreeding coefficient. The Global MAF of Loci,
#' with STACKS GBS/RAD loci = read or de novo haplotypes,
#' is included and repeated over SNP.
#' @param data The tidy VCF file created with \link{tidy_genomic_data}.
#' @param filename (optional) Name of the file written to the working directory.
#' @rdname summary_stats_vcf_tidy
#' @export
#' @keywords internal

summary_stats_vcf_tidy <- function(data, filename = NULL) {
  cat("#######################################################################\n")
  cat("################## radiator::summary_stats_vcf_tidy ###################\n")
  cat("#######################################################################\n")
  timing <- proc.time()

  # Checking for missing and/or default arguments-------------------------------
  if (missing(data)) rlang::abort("Input file missing")

  vcf.summary <- data %>%
    dplyr::filter(GT_VCF != "./.") %>%
    dplyr::group_by(LOCUS, POS, POP_ID) %>%
    dplyr::summarise(
      N = as.numeric(n()),
      PP = as.numeric(length(GT_VCF[GT_VCF == "0/0"])),
      PQ = as.numeric(length(GT_VCF[GT_VCF == "1/0" | GT_VCF == "0/1"])),
      QQ = as.numeric(length(GT_VCF[GT_VCF == "1/1"]))
    ) %>%
    dplyr::mutate(
      FREQ_REF = ((PP*2) + PQ)/(2*N),
      FREQ_ALT = ((QQ*2) + PQ)/(2*N),
      HET_O = PQ/N,
      HET_E = 2 * FREQ_REF * FREQ_ALT,
      FIS = dplyr::if_else(HET_O == 0, 1, round(((HET_E - HET_O) / HET_E), 6))
    )

  global.maf <- vcf.summary %>%
    dplyr::group_by(LOCUS, POS) %>%
    dplyr::summarise_at(.tbl = ., .vars = c("N", "PP", "PQ", "QQ"), .funs = sum) %>%
    dplyr::mutate(GLOBAL_MAF = (PQ + (2 * QQ)) / (2*N)) %>%
    dplyr::select(LOCUS, POS, GLOBAL_MAF)

  vcf.summary <- dplyr::left_join(global.maf, vcf.summary, by = c("LOCUS", "POS"))

  vcf.summary <- vcf.summary[c("LOCUS", "POS", "POP_ID", "N", "PP", "PQ", "QQ", "FREQ_REF", "FREQ_ALT", "GLOBAL_MAF", "HET_O", "HET_E", "FIS")]

  if (!is.null(filename)) {
    readr::write_tsv(vcf.summary, filename, append = FALSE, col_names = TRUE)
    message("writing the summary file in your working directory: \n", getwd())
  } else {
    message("writing the summary file to the working directory: not selected")
  }
  timing <- proc.time() - timing
  message("Computation time: ", round(timing[[3]]), " sec")
  cat("############################## completed ##############################\n")
  return(vcf.summary)
}#End summary_stats_vcf_tidy

#' @title Table of low coverage genotypes
#' @description This function create a table summary of the genotypes
#' below a user-define threshold.
#' coverage statistics by populations.

#' @param data A tidy vcf object or file (using ".tsv"),
#' created with tidy_genomic_data.
#' @param pop.levels Character string defining your ordered populations.
#' @param read.depth.threshold The read depth threshold to evaluate.
#' @param filename.low.coverage Filename of the low coverage table written
#' in the working directory.
#' @param filename.low.coverage.imbalance Filename of ...

#' @return a list of 2 tables (accessed with $). The values in the tables
#' represent percentage of samples.

#' @details work in progress....
#' Table 1: low coverage summary $low.coverage.summary (homo- and
#' hetero- zygotes genotypes).
#' Table 2: summary of coverage imbalance between alleles in the heterozygotes.
#' 0/0 : homozygote REF allele.
#' 1/1 : homozygote ALT allele.
#' 0/1 : heterozygote with coverage REF > ALT allele.
#' 1/0 : heterozygote with coverage REF < ALT allele.

#' @importFrom tidyr spread
#' @rdname table_low_coverage_summary
#' @export
#' @keywords internal


table_low_coverage_summary <- function(
  data,
  pop.levels = NULL,
  read.depth.threshold,
  filename.low.coverage,
  filename.low.coverage.imbalance
) {

  ALLELE_COVERAGE_RATIO <- SAMPLES_NUMBER <- TOTAL_NUMBER <- IMBALANCE_NUMBER <- NULL

  if (is.vector(data)) {
    data <- readr::read_tsv(
      data,
      col_names = TRUE,
      col_types = "diidccddccccdddddc"
    ) %>%
      dplyr::mutate(INDIVIDUALS = factor(INDIVIDUALS))
    message("Using the file in your directory")

  } else {
    data <- data
    message("Using the file from your global environment")

  }

  if (missing(pop.levels) == "TRUE") {
    data <- data

  } else {
    data <- data %>%
      dplyr::mutate(POP_ID = factor(POP_ID, levels = pop.levels, ordered = T))
  }

  low.coverage.summary <- data %>%
    dplyr::filter(GT != "./.") %>%
    dplyr::select(GT, POP_ID) %>%
    dplyr::group_by(GT, POP_ID) %>%
    dplyr::summarise(
      TOTAL_NUMBER = n()
    ) %>%
    dplyr::full_join(
      data %>%
        dplyr::filter(READ_DEPTH < read.depth.threshold & GT != "./.") %>%
        dplyr::group_by(GT, POP_ID) %>%
        dplyr::summarise(
          SAMPLES_NUMBER = n()
        ),
      by = c("GT", "POP_ID")
    ) %>%
    dplyr::full_join(
      data %>%
        dplyr::filter(READ_DEPTH < read.depth.threshold & GT != "./.") %>%
        dplyr::filter(ALLELE_COVERAGE_RATIO != "NA" & ALLELE_COVERAGE_RATIO != 0 ) %>%
        dplyr::group_by(GT, POP_ID) %>%
        dplyr::summarise(
          IMBALANCE_NUMBER = n()
        ),
      by = c("GT", "POP_ID")
    ) %>%
    dplyr::mutate(
      LOW_COVERAGE_PERCENT = round(SAMPLES_NUMBER / TOTAL_NUMBER * 100, 2),
      IMBALANCE_PERCENT = round(IMBALANCE_NUMBER / TOTAL_NUMBER * 100, 2)
    )

  low.coverage.summary.table <- low.coverage.summary %>%
    tidyr::spread(data = ., key = GT, value = LOW_COVERAGE_PERCENT, -POP_ID)

  utils::write.table(low.coverage.summary.table, filename.low.coverage, sep = "\t", row.names = F, col.names = T, quote = F)

  low.coverage.imbalance.summary.table <- low.coverage.summary %>%
    dplyr::filter(GT != "0/0" & GT != "1/1") %>%
    tidyr::spread(data = ., key = GT, value = IMBALANCE_PERCENT, -POP_ID)

  utils::write.table(low.coverage.imbalance.summary.table, filename.low.coverage.imbalance, sep = "\t", row.names = F, col.names = T, quote = F)

  invisible(
    cat(
      sprintf(
        "2 files:
%s
%s\n
Written in the directory:
%s",
        filename.low.coverage, filename.low.coverage.imbalance, getwd()
      )))

  res <- list()
  res$low.coverage.summary <- low.coverage.summary.table
  res$heterozygote.imbalance <- low.coverage.imbalance.summary.table

  return(res)
}# End table_low_coverage_summary


## Genotype likelihood ###

#' @title Genotype likelihood summary

#' @description This function create 3 summary tables with
#' genotype likelihood statistics. Individuals, populations and markers information is provided.
#' The input data is a tidy vcf created.
#' with \link{tidy_genomic_data}.

#' @param data A tidy vcf object or file (using ".tsv"),
#' created with tidy_genomic_data.
#' @param pop.levels (optional) Character string defining your ordered populations.
#' Default: \code{pop.levels = NULL}
#' @param gl.approach Character. By \code{"SNP"} or by \code{"haplotype"}.
#' The function will consider the SNP or haplotype GL statistics to filter the marker.
#' Default: \code{gl.approach = "haplotype"}.
#' @param folder (optional) Path to folder where results will be saved. With the
#' default, \code{folder = NULL}, the working directory is used.

#' @return A list with 3 object (data frames): the statistics at the
#' individual level ($gl.individuals,
#' saved as: "genotype.likelihood.individual.summary.tsv"), a long format with information by
#' marker and population ($gl.summary.marker.pop) and the overall summary by
#' population ($gl.summary.pop). Those 3 data frames are also written in your
#' working directory.

#' @details The table contains summary statistics: mean, median, min, max and
#' diff (max-min), of genotype likelihood for individuals, markers and populations. To access
#' the two tables, use $. The table that summarise by populations was created
#' using average nested: loci -> individuals -> populations.
#' Using the haplotype approach: The summary statistics are averaged
#' and nested SNP -> individuals -> population -> loci. e.g. the mean GL is the average
#' genotype likelihood for all individuals of pop x for loci x.

#' @importFrom stats median
#' @importFrom dplyr group_by mutate summarise ungroup arrange summarise_at
#' @importFrom readr write_tsv
#' @importFrom stringi stri_join
#' @importFrom tidyr gather
#'
#' @export
#' @rdname summary_genotype_likelihood
#' @keywords internal

summary_genotype_likelihood <- function(
  data,
  pop.levels = NULL,
  gl.approach = "haplotype",
  folder = NULL
){
  if (missing(data)) rlang::abort("missing input file")

  if (gl.approach == "haplotype") {
    message("Approach: haplotype")
  } else {
    message("Approach: SNP")
  }

  # set the working directory to save results
  if (is.null(folder)) folder <- getwd()

  # import data
  message("Importing data")

  input <- radiator::tidy_wide(data = data, import.metadata = TRUE)

  if (!tibble::has_name(input, "GL") & !tibble::has_name(input, "PL") ) {
    rlang::abort("GL or PL information is required")
  }

  pop.levels <- levels(input$POP_ID)
  pop.labels <- pop.levels

  # summary individuals
  message("Summarising GL info: individual")
  gl.individuals <- input %>%
    dplyr::group_by(POP_ID, INDIVIDUALS) %>%
    dplyr::summarise(
      GL_MEAN = mean(GL, na.rm = TRUE),
      GL_MEDIAN = stats::median(GL, na.rm = TRUE),
      GL_MIN = min(GL, na.rm = TRUE),
      GL_MAX = max(GL, na.rm = TRUE),
      GL_DIFF = GL_MAX - GL_MIN
    ) %>%
    dplyr::arrange(POP_ID, INDIVIDUALS)
  readr::write_tsv(x = gl.individuals, path = stringi::stri_join(folder,"/genotype.likelihood.individual.summary.tsv"), col_names = TRUE)

  # summary by locus and pop
  message("Summarising GL info: marker and pop")
  if (gl.approach == "haplotype") { # by haplotype
    gl.summary.marker.pop <- input %>%
      dplyr::ungroup(.) %>%
      dplyr::group_by(LOCUS, POP_ID) %>%
      dplyr::summarise(
        GL_MEAN = mean(GL, na.rm = T),
        GL_MEDIAN = stats::median(GL, na.rm = T),
        GL_MIN = min(GL, na.rm = T),
        GL_MAX = max(GL, na.rm = T),
        GL_DIFF = GL_MAX - GL_MIN
      ) %>%
      tidyr::gather(data = ., key = "GENOTYPE_LIKELIHOOD_GROUP", value = "VALUE", -c(LOCUS, POP_ID)) %>%
      dplyr::arrange(LOCUS, POP_ID, GENOTYPE_LIKELIHOOD_GROUP)

  } else {# by SNP
    gl.summary.marker.pop <- input %>%
      dplyr::ungroup(.) %>%
      dplyr::group_by(LOCUS, POS, POP_ID) %>%
      dplyr::summarise(
        GL_MEAN = mean(GL, na.rm = T),
        GL_MEDIAN = stats::median(GL, na.rm = T),
        GL_MIN = min(GL, na.rm = T),
        GL_MAX = max(GL, na.rm = T),
        GL_DIFF = GL_MAX - GL_MIN
      ) %>%
      tidyr::gather(data = ., key = "GENOTYPE_LIKELIHOOD_GROUP", value = "VALUE", -c(LOCUS, POS, POP_ID)) %>%
      dplyr::arrange(LOCUS, POS, POP_ID, GENOTYPE_LIKELIHOOD_GROUP)
  }
  readr::write_tsv(x = gl.summary.marker.pop, path = stringi::stri_join(folder,"/genotype.likelihood.summary.marker.pop.tsv"), col_names = TRUE)

  # Summary by pop overall locus
  message("Summarising GL info: pop")
  gl.summary.pop <- gl.individuals %>%
    dplyr::group_by(POP_ID) %>%
    dplyr::summarise_at(.tbl = ., .vars = c("GL_MEAN", "GL_MEDIAN", "GL_MIN", "GL_MAX", "GL_DIFF"), .funs = mean)

  readr::write_tsv(x = gl.summary.pop, path = stringi::stri_join(folder,"/genotype.likelihood.summary.pop.tsv"), col_names = TRUE)

  # results
  results <- list()
  results$gl.individuals <- gl.individuals
  results$gl.summary.marker.pop <- gl.summary.marker.pop
  results$gl.summary.pop <- gl.summary.pop
  return(results)
}# End summary_genotype_likelihood


## FIS
#' @title Fis summary
#' @description Fis summary (per markers) by populations and overall.

#' @inheritParams tidy_genomic_data

#' @details The tables contains summary statistics (mean, median, min, max).

#' @return A list with 2 tables: the long format of loci and populations
#' coverage statistics and the short format by populations.


#' @rdname fis_summary
#' @keywords internal
#' @export

#' @importFrom stats median
#' @importFrom dplyr select distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join funs summarise_at
#' @importFrom stringi stri_replace_all_fixed stri_join stri_sub
#' @importFrom tibble has_name as_data_frame
#' @importFrom tidyr gather spread complete separate nesting

fis_summary <- function(
  data,
  vcf.metadata = TRUE,
  whitelist.markers = NULL,
  monomorphic.out = TRUE,
  snp.ld = NULL,
  common.markers = TRUE,
  maf.thresholds = NULL,
  blacklist.id = NULL,
  pop.levels = NULL,
  pop.labels = NULL,
  strata = NULL,
  pop.select = NULL,
  filename = NULL,
  verbose = FALSE
) {
  timing <- proc.time()

  # Checking for missing and/or default arguments-------------------------------
  if (missing(data)) rlang::abort("Input file missing")
  if (!is.null(pop.levels) & is.null(pop.labels)) pop.labels <- pop.levels
  if (!is.null(pop.labels) & is.null(pop.levels)) rlang::abort("pop.levels is required if you use pop.labels")


  # Import data ---------------------------------------------------------------
  input <- suppressMessages(
    radiator::tidy_genomic_data(
      data = data,
      vcf.metadata = TRUE,
      blacklist.id = blacklist.id,
      whitelist.markers = whitelist.markers,
      monomorphic.out = monomorphic.out,
      snp.ld = snp.ld,
      common.markers = common.markers,
      maf.thresholds = maf.thresholds,
      strata = strata,
      pop.levels = pop.levels,
      pop.labels = pop.labels,
      pop.select = pop.select,
      filename = NULL,
      verbose = FALSE
    )
  )

  input$GT <- stringi::stri_replace_all_fixed(
    str = input$GT,
    pattern = c("/", ":", "_", "-", "."),
    replacement = c("", "", "", "", ""),
    vectorize_all = FALSE
  )

  pop.levels <- levels(input$POP_ID)
  pop.labels <- pop.levels

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

  # Detect if biallelic --------------------------------------------------------
  biallelic <- radiator::detect_biallelic_markers(input)


  if (tibble::has_name(input, "GT_VCF") & biallelic) {
    fis <- input %>%
      dplyr::filter(GT_VCF != "./.") %>%
      dplyr::group_by(MARKERS, POP_ID) %>%
      dplyr::summarise(
        N = n(),
        HET = length(GT_VCF[GT_VCF == "1/0" | GT_VCF == "0/1"]),
        HOM_ALT = length(GT_VCF[GT_VCF == "1/1"])
      ) %>%
      dplyr::mutate(
        FREQ_ALT = ((HOM_ALT * 2) + HET) / (2 * N),
        FREQ_REF = 1 - FREQ_ALT,
        HET_O = HET / N,
        HET_E = 2 * FREQ_REF * FREQ_ALT,
        FIS = (HET_E - HET_O) / HET_E,
        FIS = replace(FIS, which(FIS == "NaN"), 1)
      )

    fis.overall <- input %>%
      dplyr::filter(GT_VCF != "./.") %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::summarise(
        N = n(),
        HET = length(GT_VCF[GT_VCF == "1/0" | GT_VCF == "0/1"]),
        HOM_ALT = length(GT_VCF[GT_VCF == "1/1"])
      ) %>%
      dplyr::mutate(
        FREQ_ALT = ((HOM_ALT * 2) + HET) / (2 * N),
        FREQ_REF = 1 - FREQ_ALT,
        HET_O = HET / N,
        HET_E = 2 * FREQ_REF * FREQ_ALT,
        FIS = (HET_E - HET_O) / HET_E,
        FIS = replace(FIS, which(FIS == "NaN"), 1)
      )
  } else {
    input.alleles <- dplyr::select(.data = input, MARKERS, POP_ID, INDIVIDUALS, GT) %>%
      dplyr::filter(GT != "000000") %>%
      dplyr::mutate(
        A1 = stringi::stri_sub(GT, 1, 3),
        A2 = stringi::stri_sub(GT, 4,6)
      ) %>%
      dplyr::select(-GT)

    alleles.count <- input.alleles %>%
      tidyr::gather(data = ., key = ALLELE_GROUP, value = ALLELES, -c(MARKERS, INDIVIDUALS, POP_ID)) %>%
      dplyr::group_by(MARKERS, ALLELES, POP_ID) %>%
      dplyr::tally(.) %>%
      dplyr::ungroup() %>%
      tidyr::complete(data = ., POP_ID, tidyr::nesting(MARKERS, ALLELES), fill = list(n = 0)) %>%
      dplyr::group_by(MARKERS, POP_ID) %>%
      dplyr::mutate(
        FREQ = n/sum(n),
        HOM_E = FREQ^2
      ) %>%
      dplyr::summarise(HOM_E = sum(HOM_E)) %>%
      dplyr::mutate(HET_E = 1 - HOM_E)

    # Observed Het et Hom per pop and markers
    fis <- input.alleles %>%
      dplyr::mutate(IND_LEVEL_POLYMORPHISM = dplyr::if_else(A1 != A2, "HET", "HOM")) %>%
      dplyr::group_by(MARKERS, POP_ID, IND_LEVEL_POLYMORPHISM) %>%
      dplyr::tally(.) %>%
      dplyr::group_by(MARKERS, POP_ID) %>%
      tidyr::spread(data = ., IND_LEVEL_POLYMORPHISM, n, fill = 0) %>%
      dplyr::group_by(MARKERS, POP_ID) %>%
      dplyr::summarise(
        N = HOM + HET,
        HOM_O = HOM / N,
        HET_O = HET / N
      ) %>%
      dplyr::full_join(alleles.count, by = c("MARKERS", "POP_ID")) %>%
      dplyr::mutate(
        FIS = 1 - (HET_O/HET_E),
        FIS = replace(FIS, which(FIS == "NaN"), 1)
      ) %>%
      dplyr::arrange(MARKERS, POP_ID)

    # overall
    alleles.count.overall <- dplyr::select(.data = input.alleles, -POP_ID) %>%
      tidyr::gather(data = ., key = ALLELE_GROUP, value = ALLELES, -c(MARKERS, INDIVIDUALS)) %>%
      dplyr::group_by(MARKERS, ALLELES) %>%
      dplyr::tally(.) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(
        FREQ = n/sum(n),
        HOM_E = FREQ^2
      ) %>%
      dplyr::summarise(HOM_E = sum(HOM_E)) %>%
      dplyr::mutate(HET_E = 1 - HOM_E)

    # Observed Het et Hom per pop and markers
    fis.overall <- input.alleles %>%
      dplyr::mutate(IND_LEVEL_POLYMORPHISM = dplyr::if_else(A1 != A2, "HET", "HOM")) %>%
      dplyr::group_by(MARKERS, IND_LEVEL_POLYMORPHISM) %>%
      dplyr::tally(.) %>%
      dplyr::group_by(MARKERS) %>%
      tidyr::spread(data = ., IND_LEVEL_POLYMORPHISM, n, fill = 0) %>%
      dplyr::summarise(
        N = HOM + HET,
        HOM_O = HOM / N,
        HET_O = HET / N
      ) %>%
      dplyr::full_join(alleles.count.overall, by = "MARKERS") %>%
      dplyr::mutate(
        FIS = 1 - (HET_O/HET_E),
        FIS = replace(FIS, which(FIS == "NaN"), 1)
      ) %>%
      dplyr::arrange(MARKERS)
  }

  # summary
  fis.summary <- fis %>%
    dplyr::summarise(
      FIS_MEAN = mean(FIS, na.rm = TRUE),
      FIS_MEDIAN = stats::median(FIS, na.rm = TRUE),
      FIS_MIN = min(FIS, na.rm = TRUE),
      FIS_MAX = max(FIS, na.rm = TRUE),
      FIS_DIFF = FIS_MAX - FIS_MIN
    ) %>%
    dplyr::full_join(
      dplyr::select(.data = fis.overall, MARKERS, FIS_OVERALL = FIS)
      , by = "MARKERS"
    )
  res = list(fis.pop = fis, fis.overall = fis.overall, fis.summary = fis.summary)
  message(stringi::stri_join("Computation time: ", round((proc.time() - timing)[[3]]), " sec"))
  return(res)
}# End fis_summary
