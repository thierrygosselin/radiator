# write a SeqArray object from a tidy data frame

#' @name write_seqarray
#' @title Write a SeqArray GDS file from a vcf file and generate a connection object.
#' @description Write a SeqArray \href{http://zhengxwen.github.io/SeqArray/}{SeqArray}
#' file (Zheng et al. 2017) from a vcf file and generate a connection object.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param vcf The path to the vcf file.

#' @param filename (optional) The file name of the Genomic Data Structure (GDS) file.
#' radiator will append \code{.gds} to the filename.
#' If filename chosen is already present in the
#' working directory, the default \code{radiator_datetime.gds} is chosen.
#' Default: \code{filename = NULL}.

#' @inheritParams tidy_genomic_data

#' @param vcf.stats (logical) When \code{vcf.stats = TRUE}, individual's missing
#' genotype proportion, averaged heterozygosity, total coverage, mean genotype
#' coverage and marker's metadata along count for ref and alt alleles and mean
#' coverage is generated and written in the working directory.
#' Default: \code{vcf.stats = FALSE}.

#' @param ... (optional) To pass further argument for fine-tuning the function.

#' @details
#' A vcf file of 35 GB with ~4 millions SNPs take about ~7 min with 8 CPU.
#' A vcf file of 21 GB with ~2 millions SNPs take about ~5 min with 7 CPU.
#'
#' After the file is generated, it's a matter of sec to open a connection.

#' @export
#' @rdname write_seqarray

#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join
#' @importFrom stringi stri_replace_all_fixed
#' @importFrom tibble has_name
#' @importFrom tidyr spread

#' @references Zheng X, Gogarten S, Lawrence M, Stilp A, Conomos M, Weir BS,
#' Laurie C, Levine D (2017). SeqArray -- A storage-efficient high-performance
#' data format for WGS variant calls.
#' Bioinformatics.
#'
#' @references Danecek P, Auton A, Abecasis G et al. (2011)
#' The variant call format and VCFtools.
#' Bioinformatics, 27, 2156-2158.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_seqarray <- function(
  vcf,
  filename = NULL,
  vcf.stats = FALSE,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
  ...
) {

  #Test
  # vcf = populations.snps.vcf
  # filename = NULL
  # vcf.stats <- TRUE
  # parallel.core = parallel::detectCores() - 1
  # verbose = TRUE
  # snp.read.position.filter <- c("outliers", "q75", "iqr")

  res <- list() #store the results
  timing.import <- proc.time()
  # Check that SeqArray is installed
  if (!"SeqArray" %in% utils::installed.packages()[,"Package"]) {
    stop('Please install SeqArray for this option:\n
         devtools::install_github("zhengxwen/SeqArray")')
  }

  # Checking for missing and/or default arguments ------------------------------
  if (missing(vcf)) stop("vcf file missing")

  # dotslist -------------------------------------------------------------------
  dotslist <- list(...)
  want <- c("snp.read.position.filter", "mac.threshold")
  unknowned_param <- setdiff(names(dotslist), want)

  if (length(unknowned_param) > 0) {
    stop("Unknowned \"...\" parameters ",
         stringi::stri_join(unknowned_param, collapse = " "))
  }

  radiator.dots <- dotslist[names(dotslist) %in% want]
  snp.read.position.filter <- radiator.dots[["snp.read.position.filter"]]
  mac.threshold <- radiator.dots[["mac.threshold"]]

  snp.read.position.filter <- match.arg(
    arg = snp.read.position.filter,
    choices = c("outliers", "iqr", "q75"),
    several.ok = TRUE)

  mac.threshold <- 4

  # Start --------
  message("\nReading VCF...")

  # Get file size
  big.vcf <- file.size(vcf)

  if (verbose) {
    if (big.vcf > 500000000) message("Large vcf file may take several minutes...")
    if (big.vcf > 5000000000) message("Actually, you have time for a coffee...")
  }

  # Filename -------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (is.null(filename)) {
    ind.file <- stringi::stri_join("vcf_individuals_info_", file.date, ".tsv")
    markers.file <- stringi::stri_join("vcf_markers_metadata_", file.date, ".tsv")
    filename <- stringi::stri_join("radiator_", file.date, ".gds")
  } else {
    ind.file <- stringi::stri_join(filename, "_vcf_individuals_info_", file.date, ".tsv")
    markers.file <- stringi::stri_join(filename, "_vcf_markers_metadata_", file.date, ".tsv")
    filename.problem <- file.exists(filename)
    if (filename.problem) {
      filename <- stringi::stri_join(filename, "_", file.date, ".gds")
    } else {
      filename <- stringi::stri_join(filename, ".gds")
    }
  }

  # Read vcf -------------------------------------------------------------------
  timing.vcf <- proc.time()
  res$vcf.connection <- SeqArray::seqVCF2GDS(
    vcf.fn = vcf,
    out.fn = filename,
    parallel = parallel.core,
    storage.option = "ZIP_RA",
    verbose = FALSE) %>%
    SeqArray::seqOpen(gds.fn = .)
  if (verbose) message("\nconversion timing: ", round((proc.time() - timing.vcf)[[3]]), " sec")

  if (verbose) {
    message("\nSeqArray GDS file generated: ", filename)
    message("To close the connection use SeqArray::seqClose(OBJECT_NAME$vcf.connection)")
  }

  # Strata ---------------------------------------------------------------------
  strata <- radiator::read_strata(strata) %$% strata
  res$individuals <- tibble::tibble(
    INDIVIDUALS_VCF = SeqArray::seqGetData(res$vcf.connection, "sample.id")) %>%
    dplyr::mutate(INDIVIDUALS = radiator::clean_ind_names(INDIVIDUALS_VCF))

  # test <- res$individuals
  gdsfmt::add.gdsn(
    node = res$vcf.connection,
    name = "sample.id",
    val = res$individuals$INDIVIDUALS,
    replace = TRUE,
    compress = "ZIP_RA",
    closezip = TRUE)

  res$individuals %<>% dplyr::left_join(strata, by = "INDIVIDUALS")

  SeqArray::seqSetFilter(object = res$vcf.connection,
                         sample.id = res$individuals$INDIVIDUALS,
                         action = "set",
                         verbose = FALSE)
  # SeqArray::seqResetFilter(object = res$vcf.connection, sample = TRUE, variant = TRUE, verbose = TRUE)

  # Markers metadata  -----------------------------------------------------------
  if (vcf.stats) {
    if (verbose) message("radiator is working on the file ...")
    res$markers.meta <- tibble::tibble(
      VARIANT_ID = SeqArray::seqGetData(res$vcf.connection, "variant.id"),
      CHROM = SeqArray::seqGetData(res$vcf.connection, "chromosome"),
      LOCUS = SeqArray::seqGetData(res$vcf.connection, "annotation/id"),
      POS = SeqArray::seqGetData(res$vcf.connection, "position"))


    # detect ":" or "_" in locus id
    if (stringi::stri_detect_regex(str = res$markers.meta[1,3], pattern = "[^[:alnum:]]+")) {
      res$markers.meta <- tidyr::separate(
        data = res$markers.meta,
        col = LOCUS, into = c("LOCUS", "COL"),
        # sep = ":", # To extract rom different stacks version...
        extra = "drop", remove = TRUE, convert = TRUE)
      # SeqArray::seqClose(object = res$vcf.connection)
      # SeqVarTools::setVariantID(filename, res$markers.meta$LOCUS)
      # res$vcf.connection <- SeqArray::seqOpen(gds.fn = filename)
      # check <- tibble::tibble(LOCUS = SeqArray::seqGetData(res$vcf.connection, "variant.id"))
    }

    # Markers stats ------------------------------------------------------------
    if (verbose) message("Updating markers metadata and stats")
    n.markers <- dplyr::n_distinct(res$markers.meta$VARIANT_ID)
    want <- c("VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS",
              "COL","MISSING_PROP", "REF_COUNT", "MAC", "SNP_PER_LOCUS")
    res$markers.meta  <- suppressWarnings(
      dplyr::bind_cols(
        dplyr::mutate(res$markers.meta,
                      MARKERS = stringi::stri_join(CHROM, LOCUS, POS, sep = "__")),
        SeqArray::seqAlleleCount(
          gdsfile = res$vcf.connection,
          ref.allele = NULL,
          .progress = TRUE,
          parallel = parallel.core) %>%
          unlist(.) %>%
          matrix(
            data = .,
            nrow = n.markers, ncol = 2, byrow = TRUE,
            dimnames = list(rownames = res$markers.stats$VARIANT_ID,
                            colnames = c("REF_COUNT", "ALT_COUNT"))) %>%
          tibble::as_tibble(.),
        tibble::tibble(
          MISSING_PROP = SeqArray::seqMissing(
            gdsfile = res$vcf.connection,
            per.variant = TRUE, .progress = TRUE, parallel = parallel.core))) %>%
        dplyr::mutate(
          MAC = dplyr::if_else(ALT_COUNT < REF_COUNT, ALT_COUNT, REF_COUNT),
          ALT_COUNT = NULL) %>%
        dplyr::group_by(LOCUS) %>%
        dplyr::mutate(SNP_PER_LOCUS = n()) %>%
        dplyr::ungroup(.) %>%
        dplyr::select(want)
    )
    # test <- res$markers.meta
    # when low MAC is remove that number drops...
    if (!is.null(mac.threshold)) {
      res$markers.meta %<>%
        dplyr::mutate(
          MAC_FILTER = dplyr::if_else(MAC >= 4, TRUE, FALSE))
    } else {
      res$markers.meta %<>%
        dplyr::mutate(MAC_FILTER = TRUE)
    }
    # test <- res$markers.meta

    # position of the SNPs on the read -----------------------------------------
    if (tibble::has_name(res$markers.meta, "COL")) {
      if (verbose) message("Generating SNP position on read stats")
      res$snp.col.stats <- tibble_stats(
        x = dplyr::filter(res$markers.meta, MAC_FILTER) %>% # using filtered dataset
          dplyr::distinct(MARKERS,COL) %$%
          COL,
        group = "snp position on read") #%>% dplyr::select(-OUTLIERS_LOW)
      snp.col.iqr.threshold <- c(res$snp.col.stats$Q25, res$snp.col.stats$Q75)

      # outliers
      if ("outliers" %in% snp.read.position.filter) {
        res$markers.meta %<>%
          dplyr::mutate(
            SNP_POS_READ_OUTLIERS = dplyr::if_else(
              COL < res$snp.col.stats$OUTLIERS_HIGH, TRUE, FALSE))
      } else {
        res$markers.meta %<>% dplyr::mutate(SNP_POS_READ_OUTLIERS = TRUE)
      }

      # Q75
      if ("q75" %in% snp.read.position.filter) {
        res$markers.meta %<>%
          dplyr::mutate(
            SNP_POS_READ_Q75 = dplyr::if_else(
              COL <= snp.col.iqr.threshold[2], TRUE, FALSE))
      } else {
        res$markers.meta %<>% dplyr::mutate(SNP_POS_READ_Q75 = TRUE)
      }

      # IQR
      if ("iqr" %in% snp.read.position.filter) {
        res$markers.meta %<>%
          dplyr::mutate(
            SNP_POS_READ_IQR = dplyr::if_else(
              COL >= snp.col.iqr.threshold[1] & COL <= snp.col.iqr.threshold[2], TRUE, FALSE))
      } else {
        res$markers.meta %<>% dplyr::mutate(SNP_POS_READ_IQR = TRUE)
      }

      # Update
      res$snp.col.stats %<>% dplyr::bind_rows(
        tibble_stats(
          x = unique.col <- res$markers.meta %>%
            dplyr::filter(MAC_FILTER) %>% # mac filter
            dplyr::filter(SNP_POS_READ_OUTLIERS) %>% # outlier snp on read filter
            dplyr::filter(SNP_POS_READ_Q75) %>% # snp below Q75
            dplyr::filter(SNP_POS_READ_IQR) %>% # focus on snp within IQR range
            dplyr::distinct(MARKERS, COL) %$% COL,
          group = "snp position on read filtered")) #%>% dplyr::select(-OUTLIERS_LOW))
      snp.col.iqr.threshold <- NULL
      # test <- res$snp.col.stats
      # Generate box plot
      res$snp.pos.read.fig <- boxplot_stats(
        data = res$snp.col.stats,
        title = "Impact of filter on the SNP position on the read",
        x.axis.title = "SNP position on the read groupings",
        y.axis.title = "SNP position (base pair)",
        bp.filename = "vcf.snp.position.read.pdf")
    }
    # number of SNPs per locus -------------------------------------------------
    # before and after filtering
    res$snp.locus.stats <- tibble_stats(
      x = res$markers.meta$SNP_PER_LOCUS,
      group =  "unfiltered") %>%
      dplyr::bind_rows(
        tibble_stats(
          x = res$markers.meta %>%
            dplyr::select(LOCUS, MAC_FILTER) %>%
            dplyr::filter(MAC_FILTER) %>%
            dplyr::group_by(LOCUS) %>%
            dplyr::tally(.) %>%
            dplyr::rename(SNP_PER_LOCUS_MAC = n) %$% SNP_PER_LOCUS_MAC,
          group = "after MAC filter"),
        tibble_stats(
          x = res$markers.meta %>%
            dplyr::select(LOCUS, MAC_FILTER, SNP_POS_READ_OUTLIERS,
                          SNP_POS_READ_Q75, SNP_POS_READ_IQR) %>%
            dplyr::filter(MAC_FILTER) %>% # mac filter
            dplyr::filter(SNP_POS_READ_OUTLIERS) %>% # outlier snp on read filter
            dplyr::filter(SNP_POS_READ_Q75) %>% # snp below Q75
            dplyr::filter(SNP_POS_READ_IQR) %>%
            dplyr::group_by(LOCUS) %>%
            dplyr::tally(.) %$% n,
          group = "after SNP read position filters")) %>%
      dplyr::mutate(GROUP = factor(
        x = GROUP,
        levels = c("unfiltered",
                   "after MAC filter",
                   "after SNP read position filters"))) %>%
      dplyr::arrange(GROUP)

    # Generate box plot
    res$snp.per.locus.fig <- boxplot_stats(
      data = res$snp.locus.stats,
      title = "Impact of filters on the number of SNPs per locus",
      x.axis.title = "groupings",
      y.axis.title = "Number of SNPs per locus",
      bp.filename = "vcf.number.snp.per.locus.pdf")

    # get unique locus id and generate a variant.locus.id vector ---------------
    variant.locus.id <- res$markers.meta %>%
      dplyr::filter(MAC_FILTER) %>% # mac filter
      dplyr::filter(SNP_POS_READ_OUTLIERS) %>% # outlier snp on read filter
      dplyr::filter(SNP_POS_READ_Q75) %>% # snp below Q75
      dplyr::filter(SNP_POS_READ_IQR) %>%
      dplyr::select(LOCUS, VARIANT_ID) %>%
      dplyr::distinct(LOCUS, .keep_all = TRUE) %>%
      dplyr::select(VARIANT_ID) %>%
      purrr::flatten_int(.)

    # SeqArray::seqSetFilter(object = res$vcf.connection, variant.id = variant.locus.id, verbose = FALSE)
    # SeqArray::seqResetFilter(object = res$vcf.connection, verbose = FALSE)

    if (verbose) message("Generating coverage stats")
    coverage.info <- extract_coverage(data = res$vcf.connection,
                                      variant.locus.id = variant.locus.id)

    #Coverage: MARKERS
    res$markers.meta %<>% dplyr::mutate(COVERAGE_MEAN = coverage.info$markers.mean) %>%
      readr::write_tsv(x = ., path = markers.file)

    # Individuals stats --------------------------------------------------------
    if (verbose) message("Generating individual stats")
    res$individuals %<>% dplyr::mutate(
      # INDIVIDUALS = SeqArray::seqGetData(res$vcf.connection, "sample.id"),
      MISSING_PROP = round(SeqArray::seqMissing(
        gdsfile = res$vcf.connection, per.variant = FALSE,
        .progress = TRUE,
        parallel = parallel.core), 6),
      HETEROZYGOSITY = round(SeqVarTools::heterozygosity(
        gdsobj = res$vcf.connection, margin = "by.sample", use.names = FALSE), 6),
      COVERAGE_TOTAL = coverage.info$ind.cov.tot,
      COVERAGE_MEAN = coverage.info$ind.cov.mean
    ) %>%
      readr::write_tsv(x = ., path = ind.file)
    coverage.info <- NULL

    # test <- res$individuals
    res$ind.stats <- tibble_stats(
      x = res$individuals$MISSING_PROP,
      group = "individual's missing genotypes") %>%
      dplyr::bind_rows(
        tibble_stats(
          x = res$individuals$HETEROZYGOSITY,
          group = "individual's heterozygosity"),
        tibble_stats(
          x = as.numeric(res$individuals$COVERAGE_TOTAL),
          group = "individual's total coverage"),
        tibble_stats(
          x = as.numeric(res$individuals$COVERAGE_MEAN),
          group = "individual's mean coverage")
      )
    # test <- res$ind.stats

    ind.stats.fig <- boxplot_stats(
      data = res$ind.stats,
      title = "Individual's QC stats",
      x.axis.title = NULL,
      y.axis.title = "Statistics",
      facet.columns = TRUE,
      bp.filename = "vcf.individuals.qc.pdf")

    # Stats ----------------------------------------------------------------------
    res$n.chromosome <- dplyr::n_distinct(res$markers.meta$CHROM)
    res$n.locus <- dplyr::n_distinct(res$markers.meta$LOCUS)
    res$n.markers <- dplyr::n_distinct(res$markers.meta$MARKERS)
    res$n.individuals <- nrow(res$individuals)
    res$ind.missing <- round(mean(res$individuals$MISSING_PROP, na.rm = TRUE), 2)
    res$ind.cov.total <- round(mean(res$individuals$COVERAGE_TOTAL, na.rm = TRUE), 0)
    res$ind.cov.mean <- round(mean(res$individuals$COVERAGE_MEAN, na.rm = TRUE), 0)
    res$markers.missing <- round(mean(res$markers.meta$MISSING_PROP, na.rm = TRUE), 2)
    res$markers.cov <- round(mean(res$markers.meta$COVERAGE_MEAN, na.rm = TRUE), 0)


    if (verbose) {
      message("\n\nNumber of chromosome/contig/scaffold: ", res$n.chromosome)
      message("Number of locus: ", res$n.locus)
      message("Number of markers: ", res$n.markers)
      message("Number of individuals: ", res$n.individuals)
      message("\n\nMissing data (averaged): ")
      message("    markers: ", res$markers.missing)
      message("    individuals: ", res$ind.missing)
      message("\n\nCoverage info:")
      message("    individuals mean read depth: ", res$ind.cov.total)
      message("    individuals mean genotype coverage: ", res$ind.cov.mean)
      message("    markers mean coverage: ", res$markers.cov)
    }


  }#End stats

  # prepare data for the SeqVarTools package
  # seqOptimize("tmp.gds", target="by.sample")
  timing.import <- proc.time() - timing.import
  if (verbose) message("\nWorking time: ", round(timing.import[[3]]), " sec")
  return(res)
} # End write_SeqArray

# Internal nested Function -----------------------------------------------------
#' @title extract_info
#' @description Extract vcf information
#' @rdname extract_info
#' @keywords internal
#' @export
extract_info <- function(vcf) {
  res <- list()
  # print(vcf, all=TRUE, attribute=TRUE)
  # tic()
  vcf.info <- SeqArray::seqVCF_Header(vcf.fn = vcf, getnum = TRUE)
  # toc()
  res$vcf.source <- vcf.info$header$value[2]
  res$n.ind <- vcf.info$num.sample
  res$n.markers <- vcf.info$num.variant
  res$sample.id <- vcf.info$sample.id
  return(res)
}#End extract_info

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
    MEDIAN = stats::median(x, na.rm = TRUE),
    Q75 = Q[2],
    MAX = max(x, na.rm = TRUE),
    IQR = abs(diff(Q)),
    # IQR = stats::IQR(depth, na.rm = TRUE),
    OUTLIERS_LOW = Q25 - (1.5 * IQR),
    OUTLIERS_HIGH =  Q75 + (1.5 * IQR)) %>%
    dplyr::mutate_if(.tbl = ., .predicate = is.integer, .funs = as.numeric) %>%
    dplyr::mutate(
      OUTLIERS_LOW = dplyr::if_else(OUTLIERS_LOW < MIN, MIN, OUTLIERS_LOW),
      OUTLIERS_HIGH = dplyr::if_else(OUTLIERS_HIGH > MAX, MAX, OUTLIERS_HIGH),
      SUBSAMPLE = subsample
    )
  Q <- NULL
  return(res)
}#End tibble_stats

#' @title extract_coverage
#' @description Extract coverage information from a GDS file
#' @rdname extract_coverage
#' @keywords internal
#' @export
extract_coverage <- function(data = NULL, variant.locus.id = NULL) {
  # data <- res$vcf.connection
  coverage.info <- list()
  depth <- SeqArray::seqGetData(data, "annotation/format/DP")
  coverage.info$ind.cov.tot <- as.integer(round(rowSums(x = depth$data[,variant.locus.id], na.rm = TRUE, dims = 1L), 0))
  coverage.info$ind.cov.mean <- as.integer(round(rowMeans(x = depth$data[,variant.locus.id], na.rm = TRUE, dims = 1L), 0))
  coverage.info$markers.mean <- as.integer(round(colMeans(x = depth$data, na.rm = TRUE, dims = 1L), 0))
  return(coverage.info)
}#End extract_coverage

# boxplot of stats
#' @title boxplot_stats
#' @description Generate a boxplot
#' @rdname boxplot_stats
#' @keywords internal
#' @export
boxplot_stats <- function(data, title,
                          x.axis.title = NULL,
                          y.axis.title,
                          facet.columns = FALSE,
                          facet.rows = FALSE,
                          bp.filename = NULL) {
  # data <- test
  # x.axis.title <- "SNP position on the read groupings"
  # title <- "Impact of filter on the SNP position in the read"
  # y.axis.title <- "SNP position (base pair)"
  # bp.filename <- "vcf.snp.position.read.pdf"

  n.group <- dplyr::n_distinct(data$GROUP)
  element.text <- ggplot2::element_text(size = 10,
                                        family = "Helvetica", face = "bold")

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
    ggplot2::labs(
      y = y.axis.title,
      title = title) +
    ggplot2::theme_bw()

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
        # axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  }


  fig.boxplot <- fig.boxplot +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, family = "Helvetica",
                                         face = "bold", hjust = 0.5),
      legend.position = "none",
      axis.title.y = element.text,
      axis.text.y = element.text
    )

  if (facet.columns) {
    fig.boxplot <- fig.boxplot + ggplot2::facet_grid(GROUP ~ ., scales = "free")
    n.facet <- n.group * 2
    width <- 10
    height <- 10 + (2 * n.group)
  } else {
    width <-  10 + (5 * n.group) + 1
    height <-  10
  }

  if (facet.rows) {
    fig.boxplot <- fig.boxplot + ggplot2::facet_grid(FACET_ROWS ~ ., scales = "free")
    n.facet <- n.group * 2
    width <- 10
    height <- 10 + (2 * n.group)
  } else {
    width <-  10 + (5 * n.group) + 1
    height <-  10
  }

  print(fig.boxplot)
  if (!is.null(bp.filename)) {
    suppressMessages(ggplot2::ggsave(
      filename = bp.filename,
      plot = fig.boxplot,
      width = width,
      height = height,
      dpi = 300, units = "cm", useDingbats = FALSE))
  }
  return(fig.boxplot)
}#Endboxplot_stats
