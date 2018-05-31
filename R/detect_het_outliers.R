#' @title Detect heterozygotes outliers and estimate miscall rate
#' @description Explore departure from H-W equilibrium in bi-allelic RADseq data.
#' Highlight excess of homozygotes present in numeros RADseq studies.
#' The function estimate the genotyping error rate and heterozygote miscall rate.
#' The model focus on heterozygotes being
#' incorrectly called as homozygotes. See details below for more info.
#'
#' @param nreps (integer, optional) The number of MCMC sweeps to do.
#' Default: \code{nreps = 2000}.

#' @param burn.in (integer, optional) The number of MCMC burn-in reps.
#' With default, during execution, you will be asked to enter the nuber of burn-in.
#' For this, a plot showing the heterozygote miscall rate for all
#' the MCMC sweeps will be printed. This plot will help pinpoint the
#' number of burn-in. The remaining MCMC sweeps will be used
#' to average the heterozygote miscall rate.
#' e.g. of common value \code{burn.in = 500}.
#' With default: \code{burn.in = NULL}.

#' @inheritParams tidy_genomic_data

#' @param blacklist.markers (optional) Path to a file with markers to blacklist
#' before generating the miscall rate. Usefull to test the impact of different
#' HWE thresholds rapidly...
#' Default: \code{blacklist.markers = NULL}.

#' @details
#' \strong{Before using the function:}
#' \enumerate{
#' \item Don't use raw RADseq data, this function will work best with filtered data
#' \item Remove duplicate \code{\link[radiator]{detect_duplicate_genomes}}.
#' \item Remove mixed samples \code{\link[radiator]{detect_mixed_genomes}}.
#' \item Look at other filters in radiator package...
#' }
#'
#' \strong{During import:}
#'
#' By default the function will keep only polymorphic markers and markers common
#' between all populations. If you supply a tidy data frame or a \code{.rad} file,
#' the function skip all the filters, pop selection, etc. It will however scan and
#' remove monomorphic markers automatically.
#'
#' \strong{Keep track of the data:}
#'
#' Use the argument filename to write the imported (and maybe further filtered)
#' tidy genomic data set inside the folder. The filename will be automatically
#' appended \code{.rad} to it. This file can be used again directly inside this
#' function and other radiator functions. See \code{\link[radiator]{read_rad}}.


#' @return A folder generated automatically with date and time,
#' the file \code{het.summary.tsv} contains the summary statistics. The file
#' \code{markers.genotypes.boundaries.pdf} is the plot with boundaries.
#' The function also returns a list inside the global environment with
#' 8 objects:
#'
#' \enumerate{
#' \item input the input data, cleaned if filters were used during import.
#' \item outlier.summary a list with a tibble and plot of genotypes frequencies
#' and boundaries (also written in the folder).
#' \item summary.alt.allele a tibble summarizing the number of markers with:
#' \itemize{
#' \item no homozygote for the alternate allele (NO_HOM_ALT)
#' \item no heterozygote genotype (NO_HET)
#' \item one homozygote for the alternate allele(ONE_HOM_ALT)
#' \item one heterozygote genotype (ONE_HET)
#' \item one homozygote for the alternate allele only (ONE_HOM_ALT_ONLY)
#' \item one heterozygote genotype only (ONE_HET_ONLY)
#' \item one homozygote for the alternate allele and one heterozygote genotype only (ONE_HOM_ALT_ONE_HET_ONLY)
#' }
#' \item m.nreps A tibble with the heterozygote miscall rate for each MCMC replicate
#' \item overall.genotyping.error.rate The overall genotyping error rate
#' \item overall.m The overall heterozygote miscall rate
#' \item simmed_genos The simulated genotypes
#' }
#'
#' The statistics are summarized per population and overall,
#' the grouping is found in the last column called \code{POP_ID}.
#'

#' @examples
#' \dontrun{
#' het.prob <- radiator::detect_het_outliers(
#' data = "tuna.vcf", strata = "tuna.strata.tsv", nreps = 2000)
#' }


#' @rdname detect_het_outliers
#' @export
#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange
#' tally filter if_else mutate summarise left_join inner_join right_join anti_join
#' semi_join full_join funs
#' @importFrom ggplot2 as_labeller ggplot theme_classic theme aes geom_jitter
#' scale_y_continuous scale_y_continuous scale_color_discrete
#' scale_size_continuous theme element_blank element_text geom_hline labeller
#' facet_grid ggsave geom_boxplot labs geom_polygon geom_abline theme_bw
#' @importFrom tibble tibble
#' @importFrom tidyr complete gather unite spread nesting
#' @importFrom parallel detectCores
#' @importFrom stats rnorm rbeta runif rbinom
#' @importFrom fst read.fst


#' @author Eric Anderson \email{eric.anderson@noaa.gov} and
#' Thierry Gosselin \email{thierrygosselin@@icloud.com}

detect_het_outliers <- function (
  data,
  nreps = 2000,
  burn.in = NULL,
  blacklist.id = NULL,
  whitelist.markers = NULL,
  blacklist.markers = NULL,
  monomorphic.out = TRUE,
  max.marker = NULL,
  snp.ld = NULL,
  common.markers = TRUE,
  strata = NULL,
  pop.levels = NULL,
  pop.labels = NULL,
  pop.select = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1
) {
  opt.change <- getOption("width")
  options(width = 70)
  timing <- proc.time()
  cat("###############################################################################\n")
  cat("######################### radiator::detect_het_outliers #######################\n")
  cat("###############################################################################\n")
  res <- list() # to store results

  # manage missing arguments
  if (missing(data)) stop("missing data argument")

  # folder ---------------------------------------------------------------------
  # Get date and time to have unique filenaming
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  folder.extension <- stringi::stri_join("detect_het_outliers_", file.date, sep = "")
  path.folder <- file.path(getwd(), folder.extension)
  dir.create(path.folder)
  message(stringi::stri_join("Folder created: \n", folder.extension))
  file.date <- NULL #unused object

  # import data ----------------------------------------------------------------
  message("Importing data ...")

  # if filename argument present, add path.folder to it
  if (!is.null(filename)) filename <- file.path(path.folder, filename)

  data.type <- radiator::detect_genomic_format(data)

  if (data.type %in% c("tbl_df", "fst.file")) {
    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "INDIVIDUALS", "POP_ID", "REF", "ALT", "GT", "GT_BIN")
    message("    using tidy data frame of genotypes as input")
    message("    skipping all filters except removal of monomorphic markers")

    if (data.type == "tbl_df") {
      res$input <- suppressWarnings(dplyr::select(data, dplyr::one_of(want)))
    }
    if (data.type == "fst.file") {
      import.col <- colnames(fst::read.fst(path = data, from = 1, to = 1))
      import.col <- purrr::discard(.x = import.col, .p = !import.col %in% want)
      res$input <- fst::read.fst(path = data, columns = import.col)
      import.col <- want <- NULL
    }

    # Remove blacklisted markers -------------------------------------------------
    # Note to myself: I think this argument is useful because it allows to test
    # quickly the impact of different HWE thresholds.
    if (!is.null(blacklist.markers)) {
      message("Removing markers in the blacklist")
      blacklist.markers <- readr::read_tsv(
        file = blacklist.markers,
        col_types = readr::cols(.default = readr::col_character()))
      res$input <- dplyr::filter(res$input, !MARKERS %in% blacklist.markers$MARKERS)
    }

    # because this step skip the import and filter process we include monomorphic filter
    mono.out <- discard_monomorphic_markers(res$input, verbose = TRUE)
    mono.markers <- mono.out$blacklist.monomorphic.markers
    if (nrow(mono.markers) > 0) {
      res$input <- mono.out$input
      readr::write_tsv(mono.markers, file.path(path.folder, "blacklist.monomorphic.markers.tsv"))
    }
    mono.out <- mono.markers <- NULL
  } else {
    res$input <- radiator::tidy_genomic_data(
      data = data,
      vcf.metadata = FALSE,
      blacklist.id = blacklist.id,
      blacklist.genotype = NULL,
      whitelist.markers = whitelist.markers,
      monomorphic.out = monomorphic.out,
      max.marker = max.marker,
      snp.ld = snp.ld,
      common.markers = common.markers,
      strata = strata,
      pop.levels = pop.levels,
      pop.labels = pop.labels,
      pop.select = pop.select,
      filename = filename,
      parallel.core = parallel.core,
      verbose = FALSE
    )
  }
  data <- NULL # no need to store this extra data

  # Quick check that dataset is biallelic --------------------------------------
  biallelic <- detect_biallelic_markers(res$input)
  if (!biallelic) stop("Analysis requires a biallelic dataset")

  # Generate the summary statistics and plot -----------------------------------
  message("\nGenerating genotypes summary statistics and plot with boundaries...")
  res$outlier.summary <- plot_het_outliers(data = res$input, path.folder = path.folder)

  res$summary.alt.allele <-  dplyr::ungroup(res$outlier.summary$het.summary) %>%
    dplyr::filter(POP_ID == "OVERALL") %>%
    dplyr::summarise(
      TOTAL = n(),
      NO_HOM_ALT = length(MARKERS[HOM_ALT == 0]),
      NO_HET = length(MARKERS[HET == 0]),
      ONE_HOM_ALT = length(MARKERS[HOM_ALT == 1]),
      ONE_HET = length(MARKERS[HET == 1]),
      ONE_HOM_ALT_ONLY = length(MARKERS[HOM_ALT == 1 & HET == 0]),
      ONE_HET_ONLY = length(MARKERS[HOM_ALT == 0 & HET == 1]),
      ONE_HOM_ALT_ONE_HET_ONLY = length(MARKERS[HOM_ALT == 1 & HET == 1])
    ) %>%
    tidyr::gather(data = ., key = MARKERS, value = NUMBERS) %>%
    dplyr::mutate(PROPORTION = NUMBERS / (NUMBERS[MARKERS == "TOTAL"]))

  # Estimate heterozygotes miscall rate -------------------------------------------
  message("\nCalculating heterozygotes miscall rate...")

  mest <- estimate_m(data = res$input, nreps = nreps, m_init = 0.1)
  res$m.nreps <- tibble::tibble(iter = 1:nreps, m = mest$m)
  res$simmed_genos <- mest$simmed_genos
  res$overall.genotyping.error.rate <- tibble::tibble(overall.genotyping.error.rate = mest$overall_geno_err_est)
  mest <- NULL


  res$trace.mcmc.plot <- ggplot2::ggplot(res$m.nreps, ggplot2::aes(x = iter, y = m)) +
    ggplot2::geom_line() +
    ggplot2::labs(y= "Heterozygote miscall rate", x = "Number of MCMC sweeps") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.text.x = ggplot2::element_text(size = 10, family = "Helvetica"),
      axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.text.y = ggplot2::element_text(size = 10, family = "Helvetica")
    )
  if (!is.null(path.folder)) {
    ggplot2::ggsave(
      filename = file.path(path.folder, "trace.mcmc.plot.pdf"),
      plot = res$trace.mcmc.plot,
      width = 20, height = 10,
      dpi = 600, units = "cm",
      useDingbats = FALSE, limitsize = FALSE)
  }

  print(res$trace.mcmc.plot)

  if (is.null(burn.in)) {
    message("    The plot shows the heterozygote miscall rate for all the MCMC sweeps")
    message("    after the burn-in, the remaining MCMC sweeps will be used
    to average the heterozygote miscall rate")
    message("\n    enter the max number of burn-in reps:")

    burn.in <- as.integer(readLines(n = 1))
  }

  res$m.post.means <- dplyr::filter(res$m.nreps, iter > burn.in) %>%
    dplyr::summarise(POSTERIOR_MEAN = mean(m))

  message("\n    Overall genotyping error rate = ", round(res$overall.genotyping.error.rate, digits = 4))
  message("    Overall heterozygotes miscall rate = ", round(res$m.post.means, digits = 4))
  message("\nComputation time: ", round((proc.time() - timing)[[3]]), " sec")
  cat("################################## completed ##################################\n")
  options(width = opt.change)
  return(res)
}#End detect_het_outliers


# Internal nested functions ----------------------------------------------------
#' @title summarise genotypes
#' @description Function that summarise genotypes info
#' @rdname summarise_genotypes
#' @keywords internal
#' @export
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

summarise_genotypes <- function(data, path.folder = NULL) {
  if(is.null(path.folder)) path.folder <- getwd()

  # data.bk <- data
  # data <- data.bk
  want <- c("MARKERS", "POP_ID", "INDIVIDUALS", "GT_BIN", "READ_DEPTH")
  data <- suppressWarnings(dplyr::select(data, dplyr::one_of(want))) #%>% dplyr::filter(!is.na(GT_BIN))

  n.pop <- dplyr::n_distinct(data$POP_ID)
  if (is.factor(data$POP_ID)) {
    pop.levels <- c(levels(data$POP_ID), "OVERALL")
  } else {
    pop.levels <- c(sort(unique(data$POP_ID)), "OVERALL")
  }
  data$POP_ID <- as.character(data$POP_ID)
  replace_zero <- function(x) replace(x = x, list = which(is.na(x)), 0)

  # sum of read depth per pop and overall
  # scaled separately for pop and overall before merging
  if (tibble::has_name(data, "READ_DEPTH")) {
    rd.pop <- dplyr::select(data, MARKERS, POP_ID, READ_DEPTH) %>%
      dplyr::group_by(MARKERS, POP_ID) %>%
      dplyr::summarise(READ_DEPTH = sum(READ_DEPTH, na.rm = TRUE)) %>%
      dplyr::ungroup(.) %>%
      dplyr::mutate(READ_DEPTH_SCALED = READ_DEPTH / max(READ_DEPTH, na.rm = TRUE))

    rd <- dplyr::bind_rows(
      dplyr::select(rd.pop, MARKERS, POP_ID, READ_DEPTH = READ_DEPTH_SCALED),
      dplyr::mutate(rd.pop, POP_ID = "OVERALL") %>%
        dplyr::group_by(MARKERS, POP_ID) %>%
        dplyr::summarise(READ_DEPTH = sum(READ_DEPTH, na.rm = TRUE)) %>%
        dplyr::ungroup(.) %>%
        dplyr::mutate(READ_DEPTH = READ_DEPTH / max(READ_DEPTH, na.rm = TRUE))
    ) %>%
      dplyr::arrange(MARKERS, POP_ID) %>%
      dplyr::select(-MARKERS, -POP_ID)
    rd.pop <- NULL
  } else {
    rd <- NULL
  }

  pop <- data %>%
    dplyr::mutate(
      GT_BIN = dplyr::case_when(
        GT_BIN == 0 ~ "HOM_REF",
        GT_BIN == 1 ~ "HET",
        GT_BIN == 2 ~ "HOM_ALT",
          is.na(GT_BIN) ~ "MISSING")
    ) %>%
    dplyr::group_by(MARKERS, POP_ID, GT_BIN) %>%
    dplyr::tally(.) %>%
    data.table::as.data.table(.) %>%
    data.table::dcast.data.table(
      data = .,
      formula = MARKERS + POP_ID ~ GT_BIN,
      value.var = "n"
    ) %>%
    tibble::as_data_frame(.) %>%
    dplyr::mutate_if(.tbl = ., .predicate = is.integer, .funs = replace_zero) %>%
    dplyr::mutate(N = HOM_REF + HET + HOM_ALT)

  if (!tibble::has_name(pop, "MISSING")) {
    pop <- dplyr::mutate(pop, MISSING = as.integer("0"))
  }


  data <- dplyr::bind_rows(
    pop,
    dplyr::mutate(pop, POP_ID = "OVERALL") %>%
      dplyr::group_by(MARKERS, POP_ID) %>%
      dplyr::summarise_all(.tbl = ., .funs = sum)) %>%
    dplyr::arrange(MARKERS, POP_ID)
  pop <- NULL

  data <- data %>%
    dplyr::mutate(
      FREQ_ALT = ((HOM_ALT * 2) + HET) / (2 * N),
      FREQ_REF = 1 - FREQ_ALT,
      FREQ_HET = HET / (2 * N),
      FREQ_HOM_REF_O = HOM_REF / N,
      FREQ_HET_O = HET / N,
      FREQ_HOM_ALT_O = HOM_ALT / N,
      FREQ_HOM_REF_E = FREQ_REF^2,
      FREQ_HET_E = 2 * FREQ_REF * FREQ_ALT,
      FREQ_HOM_ALT_E = FREQ_ALT^2,
      N_HOM_REF_EXP = N * FREQ_HOM_REF_E,
      N_HET_EXP = N * FREQ_HET_E,
      N_HOM_ALT_EXP = N * FREQ_HOM_ALT_E,
      HOM_REF_Z_SCORE = (HOM_REF - N_HOM_REF_EXP) / sqrt(N * FREQ_HOM_REF_E * (1 - FREQ_HOM_REF_E)),
      HOM_HET_Z_SCORE = (HET - N_HET_EXP) / sqrt(N * FREQ_HET_E * (1 - FREQ_HET_E)),
      HOM_ALT_Z_SCORE = (HOM_ALT - N_HOM_ALT_EXP) / sqrt(N * FREQ_HOM_ALT_E * (1 - FREQ_HOM_ALT_E))
    ) %>%
    dplyr::bind_cols(rd) %>%
    dplyr::mutate(POP_ID = factor(POP_ID, levels = pop.levels, ordered = TRUE)) %>%
    readr::write_tsv(x = ., path = file.path(path.folder, "genotypes.summary.tsv"))
  rd <- NULL
  return(data)
}#End summarise_genotypes

# Generate genotypes frequencies boundaries

#' @title plot_het_outliers
#' @description Function that calculates genotypes observed and expected frequencies
#' returns a tibble with all the summary info and a plot with boundaries.
#' @rdname plot_het_outliers
#' @keywords internal
#' @export
#' @author Eric Anderson \email{eric.anderson@noaa.gov} and
#' Thierry Gosselin \email{thierrygosselin@@icloud.com}

plot_het_outliers <- function(data, path.folder = NULL) {
  res <- list() # to store results
  n.pop <- dplyr::n_distinct(data$POP_ID)
  res$het.summary <- summarise_genotypes(data, path.folder = path.folder)

  # prepare data for figure
  freq.summary <- dplyr::bind_cols(
    res$het.summary %>%
      dplyr::select(MARKERS, POP_ID, HOM_REF = FREQ_HOM_REF_O, HET = FREQ_HET_O, HOM_ALT = FREQ_HOM_ALT_O) %>%
      tidyr::gather(data = ., key = GENOTYPES, value = OBSERVED, -c(MARKERS, POP_ID)) %>%
      dplyr::arrange(MARKERS, POP_ID),
    res$het.summary %>%
      dplyr::select(MARKERS, POP_ID, HOM_REF = FREQ_HOM_REF_E, HET = FREQ_HET_E, HOM_ALT = FREQ_HOM_ALT_E) %>%
      tidyr::gather(data = ., key = GENOTYPES, value = EXPECTED, -c(MARKERS, POP_ID)) %>%
      dplyr::arrange(MARKERS, POP_ID) %>%
      dplyr::select(EXPECTED)
  ) %>%
    dplyr::mutate(GENOTYPES = factor(
      GENOTYPES,
      levels = c("HOM_REF", "HET", "HOM_ALT"),
      labels = c("Homozygote REF allele", "Heterozygote", "Homozygote ALT allele"),
      ordered = TRUE))

  # generate boundaries
  boundaries <- generate_geno_freq_boundaries() %>%
    dplyr::mutate(
      GENOTYPES = factor(
        GENOTYPES,
        levels = c("Homozygote REF allele", "Heterozygote", "Homozygote ALT allele")))

  res$gt.boundaries.plot <- ggplot2::ggplot(freq.summary , ggplot2::aes(x = EXPECTED, y = OBSERVED, colour = GENOTYPES)) +
    ggplot2::geom_jitter(alpha = 0.1, width = 0.01, height = 0.01) +
    ggplot2::geom_polygon(data = boundaries, fill = NA, linetype = "dashed", colour = "black") +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "solid") +
    ggplot2::labs(x = "Genotypes (expected frequency) ", y = "Genotypes (observed frequency)") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.text.x = ggplot2::element_text(size = 10, family = "Helvetica"),
      axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.text.y = ggplot2::element_text(size = 10, family = "Helvetica")
    ) +
    ggplot2::facet_grid(POP_ID ~ GENOTYPES)

  if (!is.null(path.folder)) {
    ggplot2::ggsave(
      filename = file.path(path.folder, "markers.genotypes.boundaries.pdf"),
      plot = res$gt.boundaries.plot,
      width = 20, height = (n.pop + 1) * 5,# + 1 for overall always present
      dpi = 300, units = "cm",
      useDingbats = FALSE, limitsize = FALSE)
  }
  return(res)
}#End plot_het_outliers


#' @title generate_geno_freq_boundaries
#' @description Function that returns a tibble with the min/max values of
#' genotype freqs possible.
#'
#' These mins and maxes occur because the genotypes are used to estimate
#' the allele frequencie.
#' @rdname generate_geno_freq_boundaries
#' @keywords internal
#' @export
#' @author Eric Anderson \email{eric.anderson@noaa.gov}

generate_geno_freq_boundaries <- function() {
  # first do it for the homozygote category
  Poe <- seq(0,1, by = 0.005)
  phat <- 1 - sqrt(Poe)
  minPo <- pmax(0, 1 - 2 * phat)
  maxPo <- 1 - phat

  # these are the values for the two homozygote categories
  homo_tib <- tibble::tibble(
    EXPECTED = rep(Poe, 4),
    OBSERVED = rep(c(minPo, maxPo), 2),
    GENOTYPES = as.character(rep(c("Homozygote REF allele", "Homozygote ALT allele"), each = length(Poe) * 2)))

  # now, it should be easy to get the max/min values for heterozygotes.
  # They will occur where one of the homozygotes is min or max.
  P1e <- 2 * phat * (1 - phat)
  maxP1 <- 2 * (1 - phat - minPo)
  minP1 <- 2 * (1 - phat - maxPo)

  het_tib <- tibble::tibble(
    EXPECTED = rep(P1e, 2),
    OBSERVED = c(minP1, maxP1),
    GENOTYPES = "Heterozygote")

  dplyr::bind_rows(homo_tib, het_tib)
}#End generate_geno_freq_boundaries

#' @title heterozygotes miscall rate
#' @description estimate the heterozygotes miscall rate from a simple model
#' @param data A tidy data frame (biallelic only).
#' @param nreps number of MCMC sweeps to do.
#' Default: \code{nreps = 200}.
#' @param m_init initial starting value for m must be between 0 and 1
#' @param a0 beta parameter for reference alleles
#' @param a1 beta parameter for alternate alleles
#' @param sm standard devation of proposal distribution for m
#' @rdname estimate_m
#' @keywords internal
#' @export
#' @author Eric Anderson \email{eric.anderson@noaa.gov}

estimate_m <- function(
  data,
  nreps = 200,
  m_init = stats::runif(1),
  a0 = 0.5,
  a1 = 0.5,
  sm = 0.005
) {
  D <- dplyr::select(data, INDIVIDUALS, MARKERS, GT_BIN) %>%
    dplyr::group_by(INDIVIDUALS) %>%
    tidyr::spread(data = ., key = MARKERS, value = GT_BIN) %>%
    dplyr::ungroup(.) %>%
    dplyr::select(-INDIVIDUALS) %>%
    as.matrix(.)

  stopifnot(m_init > 0 & m_init < 1)

  D[is.na(D)] <- -1

  # get the N variables
  N0 <- colSums(D == 0)
  N1 <- colSums(D == 1)
  N2 <- colSums(D == 2)

  # initialize the Zs to the Ns
  Z0 <- N0
  Z1 <- N1
  Z2 <- N2

  # make some place to return the m values visited
  m <- rep(NA, nreps)
  m[1] <- m_init

  # then do the sweeps
  for (r in 2:nreps) {

    # new estimate of frequency of the "1" allele from Gibbs sampling
    p <- stats::rbeta(n = length(Z0),
                      shape1 = a1 + 2 * Z2 + Z1,
                      shape2 = a0 + 2 * Z0 + Z1)

    # propose then accept or reject a new value for m
    mprop <- m[r - 1] + stats::rnorm(1, 0, sm)
    reject <- TRUE  # reject it unless we don't
    if (mprop > 0 & mprop < 1) {
      numer <- sum(N0 * log((1 - p)^2 + mprop * p * (1 - p)) +
                     N1 * log((1 - mprop) * 2 * p * (1 - p)) +
                     N2 * log(p ^ 2 + mprop * p * (1 - p)))
      denom <- sum(N0 * log((1 - p)^2 + m[r - 1] * p * (1 - p)) +
                     N1 * log((1 - m[r - 1]) * 2 * p * (1 - p)) +
                     N2 * log(p ^ 2 + m[r - 1] * p * (1 - p)))
      if (log(stats::runif(1)) < numer - denom) {
        reject <- FALSE
      }
    }
    if (reject == FALSE) {
      m[r] <- mprop
    } else {
      m[r] <- m[r - 1]
    }

    # new values for Z from Gibbs sampling
    A0 <- stats::rbinom(n = length(N0), size = N0, prob = (m[r] * p) / (1 - p + m[r] * p))
    A2 <- stats::rbinom(n = length(N2), size = N2, prob = (m[r] * (1 - p)) / (p + m[r] * (1 - p)))

    Z0 <- N0 - A0
    Z1 <- N1 + A0 + A2
    Z2 <- N2 - A2

  }
  # return m, and eventually I need to also return the final Zs and the Ns
  # and I may as well return a new 012 file with "corrected" genotypes, which
  # I can make by broadcasting the Zs around, for example...

  # inferring/realizing/simulating genotypes. I can simulate these from their posterior
  # given the estimated allele freq and the observed genotype.  To do this I will cycle
  # over the columns (the snps) in D, and for each one, I will compute the posterior of the
  # the genotype given the observed genotype (only have to for 0's and 2's) and then I will
  # sample from those posteriors.  We have a separate function that does this
  ret <- list()
  ret$simmed_genos <- simulate_genos_from_posterior(D, p, m[nreps])

  # compute an overall genotyping error rate
  diff <- ret$simmed_genos != D
  diff[D == -1] <- NA
  ret$overall_geno_err_est <- mean(diff, na.rm = TRUE)

  # return the trace of m values
  ret$m <- m

  ret
}#End estimate_m

#' @title simulate_genos_from_posterior
#' @description simulate values for the genotypes given the observed genotype,
#' then estimate allele frequencies, and the genotyping error rate.
#'
#' This is a helper function for the estimate_m function
#' @param D an 012,-1 matrix of observed genotypes
#' @param p the estimated allele freqs
#' @param m the genotyping error rate (must be a scalar)
#' @rdname simulate_genos_from_posterior
#' @keywords internal
#' @export
#' @author Eric Anderson \email{eric.anderson@noaa.gov}

simulate_genos_from_posterior <- function(D, p, m) {
  stopifnot(length(m) == 1)

  glist <- lapply(1:ncol(D), function(i) {
    obs <- D[, i] # the observed genotypes
    pl <- p[i]  # the alle freq at locus i
    post0 <- c(
      (1 - pl) / (1 - pl + m * pl),  # posterior that observed 0 is truly a 0
      (m * pl) / (1 - pl + m * pl)   # posterior that observed 0 is truly a 1
    )
    post2 <- c(
      (m * (1 - pl)) / (pl + m * (1 - pl)),  # posterior that observed 2 is truly a 1
      pl / (pl + m * (1 - pl))               # posterior that observed 2 is truly a 2
    )
    obs[obs == 0] <- sample(x = c(0, 1), size = sum(obs == 0), replace = TRUE, prob = post0)
    obs[obs == 2] <- sample(x = c(1, 2), size = sum(obs == 2), replace = TRUE, prob = post2)
    obs
  })

  # then turn it into a matrix with the same dimensions and dimnames as D
  ret <- matrix(unlist(glist), nrow = nrow(D))
  dimnames(ret) <- dimnames(D)
  ret
}#End simulate_genos_from_posterior
