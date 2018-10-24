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

#' @param ... (optional) Advance mode that allows to pass further arguments
#' for fine-tuning the function (see details).

#' @seealso \code{\link{snp_ld}}

#' @details
#' A vcf file of 35 GB with ~4 millions SNPs take about ~7 min with 8 CPU.
#' A vcf file of 21 GB with ~2 millions SNPs take about ~5 min with 7 CPU.
#'
#' After the file is generated, it's a matter of sec to open a connection.
#'
#' \strong{Advance mode, using \emph{dots-dots-dots}}
#' \enumerate{
#' \item \code{whitelist.markers}: detailed in \code{\link[radiator]{tidy_genomic_data}}.
#' \item \code{blacklist.id}: detailed in \code{\link[radiator]{tidy_genomic_data}}.
#' \item \code{pop.select}: detailed in \code{\link[radiator]{tidy_genomic_data}}.
#' \item \code{pop.levels}: detailed in \code{\link[radiator]{tidy_genomic_data}}.
#' \item \code{pop.labels}: detailed in \code{\link[radiator]{tidy_genomic_data}}.
#' \item \code{keep.both.strands}: (logical) Radiator removes
#' duplicate SNPs found on different strands,
#' by default (\code{keep.both.strands = FALSE}).
#' \item \code{common.markers}: (logical) Default: code{common.markers = TRUE}.
#' Detailed in \code{\link[radiator]{tidy_genomic_data}}.
#' \item \code{filter.mac}: (integer) To blacklist markers below a specific MAC
#' (calculated overall).
#' \item \code{filter.coverage.outliers}: (logical) To blacklist markers with low and
#' high coverage based on outlier statistics.
#' \item \code{filter.markers.missing}: (integer) To blacklist markers with too
#' many missing data. e.g. \code{filter.markers.missing = 10}, will only keep
#' markers with missing rate <= to 10%.
#' \item \code{filter.snp.read.position}: 3 options available, \code{"outliers", "q75", "iqr"}.
#' This argument will blacklist markers based on it's position on the read.
#' \code{filter.snp.read.position = "outliers"}, will remove markers based
#' on outlier statistics of the position on the reads. e.g. if majority of SNPs
#' are found between 10 and 90 pb, and very few above and below, those markers are
#' discarded. Use this function argument with dataset with problematic assembly,
#' or \emph{de novo} assembly with undocumented or poorly selected
#' mismatch threshold.
#' \item \code{filter.short.ld}: this is the \code{snp.ld} argument
#' in \code{\link[radiator]{snp_ld}}
#' \item \code{filter.long.ld}: this is the \code{ld.threshold} argument
#' in \code{\link[radiator]{snp_ld}}
#' \item \code{long.ld.missing}: this is the \code{long.ld.missing} argument
#' in \code{\link[radiator]{snp_ld}}
#' \item \code{filter.individuals.missing}: (double) Use this argument to
#' blacklist individuals with too many missing data.
#' e.g. \code{filter.individuals.missing = 0.7}, will remove individuals with >
#' 0.3 or 30% missing genotypes. This can help discover more polymorphic markers
#' with some dataset.
#' \item \code{path.folder}: to write ouput in a specific path
#' (used internally in radiator).
#' }

#' @export
#' @rdname write_seqarray

#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join
#' @importFrom stringi stri_replace_all_fixed
#' @importFrom tibble has_name
#' @importFrom tidyr spread
# @importFrom SeqVarTools heterozygosity
# @importFrom gdsfmt add.gdsn

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
  strata = NULL,
  filename = NULL,
  vcf.stats = FALSE,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
  ...
) {

  ##Test
  # vcf = "populations.snps.vcf"
  # strata <- "spis-popmap-448samples.tsv"
  # filename <- NULL
  # vcf.stats <- TRUE
  # parallel.core <- parallel::detectCores() - 1
  # verbose <- TRUE
  # keep.both.strands = FALSE
  # common.markers = TRUE
  # filter.mac <- 4
  # filter.coverage.outliers = TRUE
  # filter.markers.missing <- 10
  # filter.snp.read.position <- c("outliers", "q75", "iqr")
  # filter.short.ld <- "maf"
  # filter.long.ld <- 0.8
  # long.ld.missing <- TRUE
  # filter.individuals.missing <- 0.7
  # blacklist.id = NULL
  # pop.select = NULL
  # pop.levels = NULL
  # pop.labels = NULL
  # whitelist.markers = NULL
  # keep.gds <- TRUE
  # path.folder = NULL

  res <- list() #store the results
  timing.import <- proc.time()
  # Check that SeqArray is installed
  if (!"SeqArray" %in% utils::installed.packages()[,"Package"]) {
    stop('Please install SeqArray for this option:\n
         devtools::install_github("zhengxwen/SeqArray")
         or the bioconductor version:
         source("https://bioconductor.org/biocLite.R")
         biocLite("SeqArray")')
  }

  if (!"SeqVarTools" %in% utils::installed.packages()[,"Package"]) {
    stop('Please install SeqVarTools for this option:\n
         source("https://bioconductor.org/biocLite.R")
         biocLite("SeqVarTools")')
  }

  if (!"gdsfmt" %in% utils::installed.packages()[,"Package"]) {
    stop('Please install gdsfmt for this option:\n
         source("https://bioconductor.org/biocLite.R")
         biocLite("gdsfmt")')
  }

  # Checking for missing and/or default arguments ------------------------------
  if (missing(vcf)) stop("vcf file missing")

  # dotslist -------------------------------------------------------------------
  dotslist <- list(...)
  want <- c("whitelist.markers",
            "filter.snp.read.position", "filter.mac",
            "filter.coverage.outliers", "filter.markers.missing", "filter.short.ld",
            "filter.long.ld", "filter.individuals.missing", "common.markers",
            "keep.both.strands",
            "blacklist.id", "pop.select", "pop.levels", "pop.labels", "keep.gds",
            "path.folder", "long.ld.missing")
  unknowned_param <- setdiff(names(dotslist), want)

  if (length(unknowned_param) > 0) {
    stop("Unknowned \"...\" parameters ",
         stringi::stri_join(unknowned_param, collapse = " "))
  }

  radiator.dots <- dotslist[names(dotslist) %in% want]
  whitelist.markers <- radiator.dots[["whitelist.markers"]]

  filter.snp.read.position <- radiator.dots[["filter.snp.read.position"]]
  filter.mac <- radiator.dots[["filter.mac"]]
  filter.coverage.outliers <- radiator.dots[["filter.coverage.outliers"]]
  filter.markers.missing <- radiator.dots[["filter.markers.missing"]]
  filter.short.ld <- radiator.dots[["filter.short.ld"]]
  filter.long.ld <- radiator.dots[["filter.long.ld"]]
  filter.individuals.missing <- radiator.dots[["filter.individuals.missing"]]
  common.markers <- radiator.dots[["common.markers"]]
  keep.both.strands <- radiator.dots[["keep.both.strands"]]
  long.ld.missing <- radiator.dots[["long.ld.missing"]]

  blacklist.id <- radiator.dots[["blacklist.id"]]
  pop.select <- radiator.dots[["pop.select"]]
  pop.levels <- radiator.dots[["pop.levels"]]
  pop.labels <- radiator.dots[["pop.labels"]]
  keep.gds <- radiator.dots[["keep.gds"]]
  path.folder <- radiator.dots[["path.folder"]]

  # useful outside this function
  if (is.null(keep.gds)) keep.gds <- TRUE
  if (is.null(path.folder)) path.folder <- getwd()
  if (is.null(keep.both.strands)) keep.both.strands <- FALSE
  if (is.null(long.ld.missing)) long.ld.missing <- FALSE
  if (is.null(filter.coverage.outliers)) filter.coverage.outliers <- FALSE
  if (is.null(common.markers)) common.markers <- TRUE

  if (!is.null(filter.snp.read.position)) {
    filter.snp.read.position <- match.arg(
      arg = filter.snp.read.position,
      choices = c("outliers", "iqr", "q75"),
      several.ok = TRUE)
  }


  # LD
  # currently: will only filter long ld if short ld as been taken care of first...
  if (!is.null(filter.long.ld) && is.null(filter.short.ld)) {
    message("\nfilter.short.ld argument set by default to: maf")
    filter.short.ld <- "maf"
  }

  # Start --------
  message("\nReading VCF...")

  # Get file size
  big.vcf <- file.size(vcf)

  if (verbose) {
    if (big.vcf > 500000000) message("Large vcf file may take several minutes...")
    if (big.vcf > 5000000000) message("    you have time for a coffee!")
  }

  # Filename -------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (is.null(filename)) {
    ind.file <- stringi::stri_join("vcf_individuals_info_", file.date, ".tsv")
    markers.file <- stringi::stri_join("vcf_markers_metadata_", file.date, ".tsv")
    filename <- stringi::stri_join("radiator_", file.date, ".gds")
    blacklist.markers <- stringi::stri_join("blacklist.markers_", file.date, ".tsv")
    blacklist.id.filename <- stringi::stri_join("blacklist.individuals_", file.date, ".tsv")
  } else {
    ind.file <- stringi::stri_join(filename, "_vcf_individuals_info_", file.date, ".tsv")
    markers.file <- stringi::stri_join(filename, "_vcf_markers_metadata_", file.date, ".tsv")
    blacklist.markers <- stringi::stri_join(filename, "_blacklist.markers_", file.date, ".tsv")
    blacklist.id.filename <- stringi::stri_join(filename, "_blacklist.individuals_", file.date, ".tsv")
    filename.problem <- file.exists(filename)
    if (filename.problem) {
      filename <- stringi::stri_join(filename, "_", file.date, ".gds")
    } else {
      filename <- stringi::stri_join(filename, ".gds")
    }
  }

  filename.short <- filename
  filename <- file.path(path.folder, filename)

  # Read vcf -------------------------------------------------------------------
  timing.vcf <- proc.time()

  # Check for bad header generated by stacks
  detect.source <- check_header_source(vcf)
  stacks.2 <- detect.source$source
  check.header <- detect.source$check.header

  res$vcf.connection <- SeqArray::seqVCF2GDS(
    vcf.fn = vcf,
    out.fn = filename,
    parallel = parallel.core,
    storage.option = "ZIP_RA",
    verbose = FALSE,
    header = check.header
  ) %>%
    SeqArray::seqOpen(gds.fn = ., readonly = FALSE)

  check.header <- detect.source <- NULL
  # Summary --------------------------------------------------------------------
  # vcf.sum <- SeqArray::seqSummary(res$vcf.connection, verbose = FALSE)
  # this might be safer, because using seqSummary seems to forget if filters were used
  n.ind <- length(SeqArray::seqGetData(res$vcf.connection, "sample.id"))
  n.markers <- length(SeqArray::seqGetData(res$vcf.connection, "variant.id"))
  message("\nNumber of SNPs: ", n.markers)
  message("Number of samples: ", n.ind)
  # vcf.sum <- NULL

  if (verbose) message("\nconversion timing: ", round((proc.time() - timing.vcf)[[3]]), " sec")

  if (verbose && keep.gds) {
    message("\nGDS file generated: \n", filename.short)
    message("To close the connection use SeqArray::seqClose(OBJECT_NAME$vcf.connection)")
  }

  # Strata ---------------------------------------------------------------------
  # import strata and filter with blacklist of id if present...
  strata <- radiator::read_strata(
    strata = strata,
    pop.levels = pop.levels,
    pop.labels = pop.labels,
    pop.select = pop.select,
    blacklist.id = blacklist.id) %$% strata

  # clean sample id in VCF -----------------------------------------------------
  res$individuals <- tibble::tibble(
    INDIVIDUALS_VCF = SeqArray::seqGetData(res$vcf.connection, "sample.id")) %>%
    dplyr::mutate(INDIVIDUALS = radiator::clean_ind_names(INDIVIDUALS_VCF))

  # replace id in VCF
  gdsfmt::add.gdsn(
    node = res$vcf.connection,
    name = "sample.id",
    val = res$individuals$INDIVIDUALS,
    replace = TRUE,
    compress = "ZIP_RA",
    closezip = TRUE)

  # radiator gds folder --------------------------------------------------------
  radiator.gds <- gdsfmt::addfolder.gdsn(
    node = res$vcf.connection,
    name = "radiator",
    replace = TRUE)

  # Sync id in STRATA and VCF --------------------------------------------------
  if (!is.null(strata)) {
    if (verbose) message("\nSynchronizing sample IDs in VCF and strata...")
    res$individuals %<>%
      dplyr::filter(INDIVIDUALS %in% strata$INDIVIDUALS) %>%
      dplyr::left_join(strata, by = "INDIVIDUALS")

    SeqArray::seqSetFilter(object = res$vcf.connection,
                           sample.id = res$individuals$INDIVIDUALS,
                           # action = "set",
                           verbose = FALSE)
    # SeqArray::seqResetFilter(object = res$vcf.connection, sample = TRUE, variant = TRUE, verbose = TRUE)

    # Add STRATA to GDS
    gdsfmt::add.gdsn(
      node = radiator.gds,
      name = "STRATA",
      val = strata$STRATA,
      replace = TRUE,
      compress = "ZIP_RA",
      closezip = TRUE)
  }

  # Markers metadata  -----------------------------------------------------------
  if (vcf.stats) {
    if (verbose) message("\nWorking on the vcf ...")
    if (stacks.2) {
      vcf.locus <- SeqArray::seqGetData(res$vcf.connection, "annotation/id")
      sample.size <- min(length(unique(vcf.locus)), 100)
      ref.genome <- sample(x = unique(vcf.locus), size = sample.size, replace = FALSE) %>%
        stringi::stri_detect_regex(str = ., pattern = "[^[:alnum:]]+") %>%
        unique
      sample.size <- NULL

      if (ref.genome) {
        res$markers.meta <- tibble::tibble(
          VARIANT_ID = SeqArray::seqGetData(res$vcf.connection, "variant.id"),
          CHROM = SeqArray::seqGetData(res$vcf.connection, "chromosome"),
          LOCUS = vcf.locus,
          POS = SeqArray::seqGetData(res$vcf.connection, "position"))
      } else {
        res$markers.meta <- tibble::tibble(
          VARIANT_ID = SeqArray::seqGetData(res$vcf.connection, "variant.id"),
          CHROM = "1",
          LOCUS = SeqArray::seqGetData(res$vcf.connection, "chromosome"),
          POS = SeqArray::seqGetData(res$vcf.connection, "position")) %>%
          dplyr::mutate(COL = POS - 1)
      }
      vcf.locus <- NULL
    } else {
      res$markers.meta <- tibble::tibble(
        VARIANT_ID = SeqArray::seqGetData(res$vcf.connection, "variant.id"),
        CHROM = SeqArray::seqGetData(res$vcf.connection, "chromosome"),
        LOCUS = SeqArray::seqGetData(res$vcf.connection, "annotation/id"),
        POS = SeqArray::seqGetData(res$vcf.connection, "position"))

      ref.genome <- sample(
        x = unique(res$markers.meta$CHROM),
        size = min(length(unique(res$markers.meta$CHROM)), 100), replace = FALSE) %>%
        stringi::stri_detect_regex(str = ., pattern = "[^[:alnum:]]+") %>%
        unique

      if (!ref.genome) {
        res$markers.meta <- res$markers.meta %>%
          dplyr::mutate(CHROM = "1")
      }
    }

    # LOCUS cleaning and Strands detection -------------------------------------
    if (stringi::stri_detect_regex(str = res$markers.meta[1,3], pattern = "[^[:alnum:]]+")) {
      locus.sep <- unique(stringi::stri_extract_all_regex(
        str = res$markers.meta[1,3],
        pattern = "[^a-zA-Z0-9-+-]",
        omit_no_match = TRUE)[[1]])

      res$markers.meta <- suppressWarnings(
        tidyr::separate(
          data = res$markers.meta,
          col = LOCUS, into = c("LOCUS", "COL", "STRANDS"),
          sep = locus.sep,
          # sep = "",
          extra = "drop", fill = "warn",
          remove = TRUE, convert = TRUE)
      )
      if(anyNA(res$markers.meta$STRANDS)) res$markers.meta$STRANDS <- NA_character_
      locus.sep <- NULL

      detect.strand <- any(stringi::stri_detect_fixed(str = unique(res$markers.meta$STRANDS), pattern = "+"))
      if(anyNA(detect.strand)) detect.strand <- FALSE

      if(detect.strand) {
        blacklist.strands <- dplyr::distinct(res$markers.meta, CHROM, LOCUS, POS) %>%
          dplyr::group_by(CHROM, POS) %>%
          dplyr::tally(.) %>%
          dplyr::filter(n > 1) %>%
          dplyr::ungroup(.) %>%
          dplyr::select(CHROM, POS)

        if (nrow(blacklist.strands) == 0) {
          blacklist.strands <- NULL
        } else {
          blacklist.strands %<>%
            dplyr::mutate_at(.tbl = .,
                             .vars = c("CHROM", "POS"),
                             .funs = radiator::clean_markers_names)
          if (verbose) {
            message("\n\nDetected ", nrow(blacklist.strands)," duplicate SNPs on different strands (+/-)")
            message("    By default radiator prune SNPs on one of the strand (-)")
            message("    To change this behavior, use argument: keep.both.strands = TRUE\n\n")
          }
        }

      } else {
        blacklist.strands <- NULL
      }
      detect.strand <- NULL
    } else {
      blacklist.strands <- NULL
    }

    # Generate MARKERS column and fix types
    res$markers.meta %<>%
      dplyr::mutate_at(
        .tbl = .,
        .vars = c("CHROM", "LOCUS", "POS"),
        .funs = radiator::clean_markers_names) %>%
      dplyr::mutate(
        MARKERS = stringi::stri_join(CHROM, LOCUS, POS, sep = "__"),
        REF = SeqArray::seqGetData(gdsfile = res$vcf.connection, var.name = "$ref"),
        ALT = SeqArray::seqGetData(gdsfile = res$vcf.connection, var.name = "$alt")
      )

    # # ADD MARKERS META to GDS --------------------------------------------------
    gdsfmt::add.gdsn(
      node = radiator.gds,
      name = "markers.meta",
      val = res$markers.meta,
      replace = TRUE,
      compress = "ZIP_RA",
      closezip = TRUE)

    # Add to GDS if we're dealing with de novo data or not
    if (ref.genome) {
      gdsfmt::add.gdsn(
        node = radiator.gds,
        name = "reference.genome",
        val = "TRUE",
        replace = FALSE,
        compress = "ZIP_RA",
        closezip = TRUE)
    } else {
      gdsfmt::add.gdsn(
        node = radiator.gds,
        name = "reference.genome",
        val = "FALSE",
        replace = FALSE,
        compress = "ZIP_RA",
        closezip = TRUE)
    }
    # Scan and filter with FILTER column ---------------------------------------
    res$markers.meta$FILTER <- SeqArray::seqGetData(
      res$vcf.connection, "annotation/filter")
    filter.check.unique <- unique(res$markers.meta$FILTER)

    if (length(filter.check.unique) > 1) {
      message("Filtering markers based on VCF FILTER column")

      variant.select <- dplyr::filter(res$markers.meta, FILTER == "PASS") %>%
        dplyr::select(VARIANT_ID) %>%
        purrr::flatten_int(.)

      n.markers.before <- nrow(res$markers.meta)
      n.markers.after <- length(variant.select)
      message("    Number of markers before = ", n.markers.before)
      message("    Number of markers removed = ", n.markers.before - n.markers.after)
      message("    Number of markers after = ", n.markers.after)

      SeqArray::seqSetFilter(object = res$vcf.connection,
                             variant.id = variant.select,
                             verbose = FALSE)
    }
    filter.check.unique <- NULL
    res$markers.meta %<>% dplyr::select(-FILTER)

    # bi- or multi-alllelic VCF ------------------------------------------------
    if (max(unique(SeqArray::seqNumAllele(gdsfile = res$vcf.connection))) - 1 > 1) {
      res$biallelic <- FALSE
      message("VCF is multi-allelic")
    } else {
      res$biallelic <- TRUE
      message("VCF is biallelic")
    }

    # add biallelic info to GDS
    gdsfmt::add.gdsn(
      node = radiator.gds,
      name = "biallelic",
      val = res$biallelic,
      replace = FALSE,
      compress = "ZIP_RA",
      closezip = TRUE)

    # Filter with whitelist of markers------------------------------------------
    if (!is.null(whitelist.markers)) {

      # Import whitelist of markers
      if (is.vector(whitelist.markers)) {
        whitelist.markers <- suppressMessages(
          readr::read_tsv(whitelist.markers, col_names = TRUE) %>%
            dplyr::mutate_all(.tbl = ., .funs = as.character))
      }
      columns.names.whitelist <- colnames(whitelist.markers)
      nrow.before <- nrow(whitelist.markers)
      whitelist.markers <- dplyr::distinct(whitelist.markers)
      nrow.after <- nrow(whitelist.markers)
      duplicate.whitelist.markers <- nrow.before - nrow.after
      if (duplicate.whitelist.markers > 0) {
        message("Whitelist of markers with ", duplicate.whitelist.markers, " duplicated identifiers...")
        message("    Creating unique whitelist")
        message("    Warning: downstream results might be impacted by this, check how you made your VCF file...")
      }
      nrow.before <- nrow.after <- duplicate.whitelist.markers <- NULL

      whitelist.markers <- dplyr::mutate_all(
        .tbl = whitelist.markers, .funs = clean_markers_names)

      if (!biallelic) {
        if (ncol(whitelist.markers) >= 3) {
          message("Note: whitelist with CHROM LOCUS POS columns and VCF haplotype:
                If the whitelist was not created from this VCF,
                the filtering could result in losing all the markers.
                The POS column is different in biallelic and multiallelic file...\n")

          message("Discarding the POS column in the whitelist")
          whitelist.markers <- dplyr::select(whitelist.markers, -POS)
        }

        if (ncol(whitelist.markers) == 1 && tibble::has_name(whitelist.markers, "MARKERS")) {
          message("Note: whitelist MARKERS column and VCF haplotype:
                If the whitelist was not created from this VCF,
                the filtering could result in losing all the markers.
                The POS column used in the MARKERS column is different in biallelic and multiallelic file...\n")
        }
      }
      message("Filtering: ", nrow(whitelist.markers), " markers in whitelist")
      columns.names.whitelist <- colnames(whitelist.markers)

      variant.select <- suppressWarnings(
        dplyr::semi_join(
          res$markers.meta, whitelist.markers, by = columns.names.whitelist)) %>%
        dplyr::select(VARIANT_ID) %>%
        purrr::flatten_int(.)

      SeqArray::seqSetFilter(object = res$vcf.connection,
                             variant.id = variant.select,
                             verbose = FALSE)

      if (length(SeqArray::seqGetData(res$vcf.connection, "variant.id")) == 0) {
        stop("No markers left in the dataset, check whitelist...")
      }
    }

    # Individuals stats --------------------------------------------------------
    if (verbose) message("Generating individual stats")
    id.stats.before.filter <- extract_coverage(data = res$vcf.connection)

    res$individuals %<>% dplyr::mutate(
      MISSING_PROP = round(SeqArray::seqMissing(
        gdsfile = res$vcf.connection, per.variant = FALSE,
        .progress = TRUE,
        parallel = parallel.core), 6),
      HETEROZYGOSITY = round(SeqVarTools::heterozygosity(
        gdsobj = res$vcf.connection, margin = "by.sample", use.names = FALSE), 6),
      COVERAGE_TOTAL = id.stats.before.filter$ind.cov.tot,
      COVERAGE_MEAN = id.stats.before.filter$ind.cov.mean
    ) %>%
      readr::write_tsv(x = ., path = file.path(path.folder, ind.file))
    id.stats.before.filter <- NULL

    # test <- res$individuals
    res$stats$ind.stats <- tibble_stats(
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
    # test <- res$stats$ind.stats

    res$figures$ind.stats.fig <- boxplot_stats(
      data = res$stats$ind.stats,
      title = "Individual's QC stats",
      x.axis.title = NULL,
      y.axis.title = "Statistics",
      facet.columns = TRUE,
      bp.filename = "vcf.individuals.qc.pdf",
      path.folder = path.folder)

    # test <- res$stats$ind.stats

    if (!is.null(filter.individuals.missing)) {

      if (!purrr::is_double(filter.individuals.missing)) {
        outlier.id.missing <- 1 - res$stats$ind.stats$OUTLIERS_HIGH[1]
        message("Removing outlier individuals based on genotyping statistics: ", outlier.id.missing)
        filter.individuals.missing <- outlier.id.missing
      }

      message("Removing individuals with too many missing genotypes")
      blacklist.id <- res$individuals %>%
        dplyr::filter(MISSING_PROP > 1 - filter.individuals.missing) %>%
        dplyr::ungroup(.) %>%
        dplyr::distinct(INDIVIDUALS)
      blacklisted.id <- nrow(blacklist.id)
      if (blacklisted.id > 0) {
        if (verbose) message("Number of individuals blacklisted: ", blacklisted.id)
        res$blacklist.id <- blacklist.id
        res$individuals %<>% dplyr::mutate(
          FILTER_INDIVIDUALS_MISSING = dplyr::if_else(
            INDIVIDUALS %in% blacklist.id$INDIVIDUALS, FALSE, TRUE)) %>%
          readr::write_tsv(x = ., path = file.path(path.folder, ind.file))
        readr::write_tsv(x = blacklist.id, path = file.path(path.folder, blacklist.id.filename))

        # update the strata and the GDS
        strata %<>% dplyr::filter(!INDIVIDUALS %in% blacklist.id$INDIVIDUALS)
        SeqArray::seqSetFilter(object = res$vcf.connection,
                               sample.id = strata$INDIVIDUALS,
                               # action = "set",
                               verbose = FALSE)
        gdsfmt::add.gdsn(
          node = radiator.gds,
          name = "STRATA",
          val = strata$STRATA,
          replace = TRUE,
          compress = "ZIP_RA",
          closezip = TRUE)
      } else {
        res$individuals$FILTER_INDIVIDUALS_MISSING = TRUE
      }
    } else {
      res$individuals$FILTER_INDIVIDUALS_MISSING = TRUE
    }#filter.individuals.missing
    blacklist.id <- NULL

    # Markers stats ------------------------------------------------------------
    if (verbose) message("Updating markers metadata and stats")
    n.markers <- dplyr::n_distinct(res$markers.meta$VARIANT_ID)
    want <- c("VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS",
              "COL", "STRANDS", "REF", "ALT", "MISSING_PROP", "MISSING_POP",
              "REF_COUNT", "MAC", "SNP_PER_LOCUS")
    res$markers.meta <- suppressWarnings(
      dplyr::bind_cols(
        res$markers.meta,
        SeqArray::seqAlleleCount(
          gdsfile = res$vcf.connection,
          ref.allele = NULL,
          .progress = TRUE,
          parallel = parallel.core) %>%
          unlist(.) %>%
          matrix(
            data = .,
            nrow = n.markers, ncol = 2, byrow = TRUE,
            dimnames = list(rownames = res$markers.meta$VARIANT_ID,
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
        dplyr::ungroup(.))

    if (is.null(strata)) {
      suppressWarnings(res$markers.meta %<>% dplyr::select(dplyr::one_of(want)))
    } else {
      suppressWarnings(res$markers.meta %<>%
                         dplyr::left_join(
                           missing_per_pop(vcf.connection = res$vcf.connection,
                                           strata = strata,
                                           parallel.core = parallel.core)
                           , by = "MARKERS") %>%
                         dplyr::select(dplyr::one_of(want)))
    }
    # test <- res$markers.meta

    # make sure GDS is updated with sample -------------------------------------
    SeqArray::seqSetFilter(object = res$vcf.connection,
                           sample.id = strata$INDIVIDUALS,
                           # action = "set",
                           verbose = FALSE)

    # Generate a blacklist of markers to store info ----------------------------
    res$blacklist.markers <- tibble::tibble(MARKERS = character(0), FILTER = character(0))

    # COMMON MARKERS -----------------------------------------------------------
    # Remove markers not in common if option selected
    if (common.markers && !is.null(strata)) {
      message("Scanning for markers in common...")
      not.common.markers <- dplyr::select(res$markers.meta, MARKERS, MISSING_POP) %>%
        tidyr::unnest(MISSING_POP) %>%
        dplyr::filter(MISSING_PROP == 1) %>%
        dplyr::distinct(MARKERS)
      not.in.common <- nrow(not.common.markers)

      if (not.in.common > 0) {
        message("    number of markers not in common between strata: ", not.in.common)
        # update the blacklist
        res$blacklist.markers <- dplyr::bind_rows(
          res$blacklist.markers,
          dplyr::mutate(not.common.markers, FILTER = "common.markers"))
        res$markers.meta %<>% dplyr::filter(!MARKERS %in% not.common.markers$MARKERS)
        # Update GDS
        gdsfmt::add.gdsn(
          node = radiator.gds,
          name = "markers.meta",
          val = res$markers.meta,
          replace = TRUE,
          compress = "ZIP_RA",
          closezip = TRUE)
        SeqArray::seqSetFilter(object = res$vcf.connection,
                               variant.id = res$markers.meta$VARIANT_ID,
                               verbose = FALSE)
      } else {
        message("    all markers in common between strata")
      }
      not.in.common <- not.common.markers <- NULL
    }

    # FILTER STRANDS -----------------------------------------------------------
    if (!is.null(blacklist.strands)) {
      blacklist.strands <- dplyr::right_join(res$markers.meta, blacklist.strands,
                                             by = c("CHROM", "POS")) %>%
        dplyr::mutate_at(.tbl = .,
                         .vars = c("MISSING_PROP", "MAC", "SNP_PER_LOCUS"),
                         .funs = as.numeric) %>%
        dplyr::group_by(CHROM, POS) %>%
        dplyr::filter(MISSING_PROP == max(MISSING_PROP)) %>%
        dplyr::filter(SNP_PER_LOCUS == max(SNP_PER_LOCUS)) %>%
        dplyr::filter(MAC == min(MAC)) %>%
        dplyr::ungroup(.) %>%
        dplyr::arrange(CHROM, POS) %>%
        dplyr::distinct(CHROM, POS, .keep_all = TRUE) %>%
        dplyr::select(MARKERS)

      if (!keep.both.strands) {
        if (nrow(blacklist.strands) > 0) {
          message("Number of duplicated markers on different strands removed: ", nrow(blacklist.strands))

          # update the blacklist
          res$blacklist.markers <- dplyr::bind_rows(
            res$blacklist.markers,
            dplyr::mutate(blacklist.strands, FILTER = "keep.both.strands"))

          # Update markers.meta
          res$markers.meta %<>% dplyr::filter(!MARKERS %in% blacklist.strands$MARKERS)
          # check <- res$markers.meta

          # Update GDS
          gdsfmt::add.gdsn(
            node = radiator.gds,
            name = "markers.meta",
            val = res$markers.meta,
            replace = TRUE,
            compress = "ZIP_RA",
            closezip = TRUE)

          SeqArray::seqSetFilter(object = res$vcf.connection,
                                 variant.id = res$markers.meta$VARIANT_ID,
                                 verbose = FALSE)
        }
      }
    }

    blacklist.strands <- NULL

    # FILTER MAC----------------------------------------------------------------
    if (!is.null(filter.mac)) {
      blacklist.mac <- dplyr::filter(res$markers.meta, MAC < filter.mac) %>%
        dplyr::select(MARKERS)

      if (nrow(blacklist.mac) > 0) {
        message("Number of markers with low MAC removed: ", nrow(blacklist.mac))
        # update the blacklist
        res$blacklist.markers <- dplyr::bind_rows(
          res$blacklist.markers,
          dplyr::mutate(blacklist.mac, FILTER = "filter.mac"))
        # check <- res$blacklist.markers

        # Update markers.meta
        res$markers.meta %<>% dplyr::filter(!MARKERS %in% blacklist.mac$MARKERS)
        # check <- res$markers.meta

        # Update GDS
        gdsfmt::add.gdsn(
          node = radiator.gds,
          name = "markers.meta",
          val = res$markers.meta,
          replace = TRUE,
          compress = "ZIP_RA",
          closezip = TRUE)

        SeqArray::seqSetFilter(object = res$vcf.connection,
                               variant.id = res$markers.meta$VARIANT_ID,
                               verbose = FALSE)
      }
      blacklist.mac <- NULL
    }
    # test <- res$markers.meta
    # Coverage Stats -----------------------------------------------------------
    if (verbose) message("Generating coverage stats")
    coverage.info <- extract_coverage(data = res$vcf.connection)

    res$stats$coverage.stats <- tibble_stats(
      x = coverage.info$markers.mean,
      group =  "unfiltered")


    res$stats$coverage.stats <- dplyr::bind_rows(
      res$stats$coverage.stats,
      tibble_stats(
        x = coverage.info$markers.mean[coverage.info$markers.mean >= res$stats$coverage.stats$OUTLIERS_LOW &
                                         coverage.info$markers.mean <= res$stats$coverage.stats$OUTLIERS_HIGH],
        group = "filtered outliers")
    ) %>%
      dplyr::mutate(GROUP = factor(
        x = GROUP,
        levels = c("unfiltered", "filtered outliers"))) %>%
      dplyr::arrange(GROUP)

    # Generate box plot
    res$figures$coverage.stats.fig <- boxplot_stats(
      data = res$stats$coverage.stats,
      title = "Impact of filtering outliers\n on coverage statistics",
      x.axis.title = NULL,
      y.axis.title = "Markers coverage (mean read depth)",
      facet.columns = TRUE,
      bp.filename = "coverage.statistics.pdf",
      path.folder = path.folder)

    res$markers.meta$COVERAGE_MEAN <- coverage.info$markers.mean

    if (filter.coverage.outliers) {
      blacklist.coverage <- res$markers.meta %>%
        dplyr::filter(
          COVERAGE_MEAN < res$stats$coverage.stats$OUTLIERS_LOW[1] |
            COVERAGE_MEAN > res$stats$coverage.stats$OUTLIERS_HIGH[1]) %>%
        dplyr::select(MARKERS)

      if (nrow(blacklist.coverage) > 0) {
        message("Number of markers with outlier coverage removed: ", nrow(blacklist.coverage))
        # update the blacklist
        res$blacklist.markers <- dplyr::bind_rows(
          res$blacklist.markers,
          dplyr::mutate(blacklist.coverage, FILTER = "filter.coverage.outliers"))
        # check <- res$blacklist.markers

        # Update markers.meta
        res$markers.meta %<>% dplyr::filter(!MARKERS %in% blacklist.coverage$MARKERS)
        # check <- res$markers.meta

        # Update GDS
        gdsfmt::add.gdsn(
          node = radiator.gds,
          name = "markers.meta",
          val = res$markers.meta,
          replace = TRUE,
          compress = "ZIP_RA",
          closezip = TRUE)

        SeqArray::seqSetFilter(object = res$vcf.connection,
                               variant.id = res$markers.meta$VARIANT_ID,
                               verbose = FALSE)
      }
      blacklist.coverage <- NULL
    }
    coverage.info <- NULL

    # Filter markers genotyping/missing ----------------------------------------
    if (!is.null(strata)) {
      res$figures$missing.markers.fig <- markers_genotyped_helper(
        x = dplyr::select(res$markers.meta, MARKERS, MISSING_POP) %>%
          tidyr::unnest(MISSING_POP) %>%
          dplyr::rename(PERCENT = MISSING_PROP, POP_ID = STRATA) %>%
          dplyr::mutate(PERCENT = PERCENT * 100),
        y = dplyr::select(res$markers.meta, MARKERS, PERCENT = MISSING_PROP) %>%
          dplyr::mutate(PERCENT = PERCENT * 100),
        overall.only = FALSE
      )
      n.pop <- dplyr::n_distinct(strata$STRATA) + 2
      ggplot2::ggsave(
        plot = res$figures$missing.markers.fig,
        filename = file.path(path.folder, "plot.markers.genotyping.rate.pdf"),
        width = n.pop * 10, height = 10,
        dpi = 300, units = "cm", useDingbats = FALSE, limitsize = FALSE)
    } else {#overall only
      res$figures$missing.markers.fig <- markers_genotyped_helper(
        x = NULL,
        y = dplyr::select(res$markers.meta, MARKERS, PERCENT = MISSING_PROP) %>%
          dplyr::mutate(PERCENT = PERCENT * 100),
        overall.only = TRUE
      )
      ggplot2::ggsave(
        plot = res$figures$missing.markers.fig,
        filename = file.path(path.folder, "plot.markers.genotyping.rate.pdf"),
        width = 15, height = 10,
        dpi = 300, units = "cm", useDingbats = FALSE)
    }

    if (!is.null(filter.markers.missing)) {
      res$markers.meta %<>%
        dplyr::mutate(
          MARKERS_MISSING = dplyr::if_else(MISSING_PROP <= filter.markers.missing / 100, TRUE, FALSE))

      blacklist.markers.missing <- res$markers.meta %>%
        dplyr::filter(MISSING_PROP > filter.markers.missing / 100) %>%
        dplyr::select(MARKERS)

      if (nrow(blacklist.markers.missing) > 0) {
        message("Number of markers with too many missing genotypes: ", nrow(blacklist.markers.missing))
        # update the blacklist
        res$blacklist.markers <- dplyr::bind_rows(
          res$blacklist.markers,
          dplyr::mutate(blacklist.markers.missing, FILTER = "filter.markers.missing"))
        # check <- res$blacklist.markers

        # Update markers.meta
        res$markers.meta %<>% dplyr::filter(!MARKERS %in% blacklist.markers.missing$MARKERS)
        # check <- res$markers.meta

        # Update GDS
        gdsfmt::add.gdsn(
          node = radiator.gds,
          name = "markers.meta",
          val = res$markers.meta,
          replace = TRUE,
          compress = "ZIP_RA",
          closezip = TRUE)

        SeqArray::seqSetFilter(object = res$vcf.connection,
                               variant.id = res$markers.meta$VARIANT_ID,
                               verbose = FALSE)
      }
      blacklist.markers.missing <- NULL
    }
    # test <- res$markers.meta

    # position of the SNPs on the read -----------------------------------------
    if (tibble::has_name(res$markers.meta, "COL")) {
      if (verbose) message("Generating SNP position on read stats")
      res$stats$snp.col.stats <- tibble_stats(
        x = dplyr::filter(res$markers.meta) %>% # using unfiltered dataset
          dplyr::distinct(MARKERS,COL) %$%
          COL,
        group = "snp position on read")
      snp.col.iqr.threshold <- c(res$stats$snp.col.stats$Q25, res$stats$snp.col.stats$Q75)

      # outliers
      if ("outliers" %in% filter.snp.read.position) {
        res$markers.meta %<>%
          dplyr::mutate(
            SNP_POS_READ_OUTLIERS = dplyr::if_else(
              COL < res$stats$snp.col.stats$OUTLIERS_HIGH, TRUE, FALSE))
      } else {
        res$markers.meta %<>% dplyr::mutate(SNP_POS_READ_OUTLIERS = TRUE)
      }

      # Q75
      if ("q75" %in% filter.snp.read.position) {
        res$markers.meta %<>%
          dplyr::mutate(
            SNP_POS_READ_Q75 = dplyr::if_else(
              COL <= snp.col.iqr.threshold[2], TRUE, FALSE))
      } else {
        res$markers.meta %<>% dplyr::mutate(SNP_POS_READ_Q75 = TRUE)
      }

      # IQR
      if ("iqr" %in% filter.snp.read.position) {
        res$markers.meta %<>%
          dplyr::mutate(
            SNP_POS_READ_IQR = dplyr::if_else(
              COL >= snp.col.iqr.threshold[1] & COL <= snp.col.iqr.threshold[2], TRUE, FALSE))
      } else {
        res$markers.meta %<>% dplyr::mutate(SNP_POS_READ_IQR = TRUE)
      }

      # Update with ALL filters available
      res$stats$snp.col.stats %<>% dplyr::bind_rows(
        tibble_stats(
          x = unique.col <- res$markers.meta %>%
            dplyr::filter(SNP_POS_READ_OUTLIERS) %>% # outlier snp on read filter
            dplyr::filter(SNP_POS_READ_Q75) %>% # snp below Q75
            dplyr::filter(SNP_POS_READ_IQR) %>% # focus on snp within IQR range
            dplyr::distinct(MARKERS, COL) %$% COL,
          group = "snp position on read filtered")) #%>% dplyr::select(-OUTLIERS_LOW))
      snp.col.iqr.threshold <- NULL
      # Generate box plot
      res$figures$snp.pos.read.fig <- boxplot_stats(
        data = res$stats$snp.col.stats,
        title = "Impact of filter\non the SNP position on the read",
        x.axis.title = NULL,
        y.axis.title = "SNP position (base pair)",
        facet.columns = TRUE,
        bp.filename = "vcf.snp.position.read.pdf",
        path.folder = path.folder)

      # Filtering ...
      if (!is.null(filter.snp.read.position)) {
        blacklist.snp.read.position <- res$markers.meta %>%
          dplyr::filter(!SNP_POS_READ_OUTLIERS | !SNP_POS_READ_Q75 | !SNP_POS_READ_IQR) %>%
          dplyr::select(MARKERS)

        if (nrow(blacklist.snp.read.position) > 0) {
          message("Number of markers blacklisted based on position on the read: ", nrow(blacklist.snp.read.position))
          # update the blacklist
          res$blacklist.markers <- dplyr::bind_rows(
            res$blacklist.markers,
            dplyr::mutate(blacklist.snp.read.position, FILTER = "filter.snp.read.position"))
          # check <- res$blacklist.markers

          # Update markers.meta
          res$markers.meta %<>% dplyr::filter(SNP_POS_READ_OUTLIERS) %>%
            dplyr::filter(SNP_POS_READ_Q75) %>%
            dplyr::filter(SNP_POS_READ_IQR) %>%
            dplyr::select(-c(SNP_POS_READ_OUTLIERS, SNP_POS_READ_Q75, SNP_POS_READ_IQR))
          # check <- res$markers.meta

          # Update GDS
          gdsfmt::add.gdsn(
            node = radiator.gds,
            name = "markers.meta",
            val = res$markers.meta,
            replace = TRUE,
            compress = "ZIP_RA",
            closezip = TRUE)

          SeqArray::seqSetFilter(object = res$vcf.connection,
                                 variant.id = res$markers.meta$VARIANT_ID,
                                 verbose = FALSE)
        }
        blacklist.snp.read.position <- NULL
      }
      # test <- res$markers.meta
    }
    # number of SNPs per locus -------------------------------------------------

    # Note to myself: there is no filtering here, just an output of the stats...
    # before and after filtering
    res$stats$snp.locus.stats <- tibble_stats(
      x = dplyr::distinct(res$markers.meta, LOCUS, SNP_PER_LOCUS) %$%
        SNP_PER_LOCUS,
      group =  "filtered")

    # read.length <- max(res$markers.meta$COL)
    # Generate box plot
    res$figures$snp.per.locus.fig <- boxplot_stats(
      data = res$stats$snp.locus.stats,
      title = "Number of SNPs per locus",
      x.axis.title = NULL,
      y.axis.title = "Number of SNPs per locus",
      bp.filename = "vcf.number.snp.per.locus.pdf",
      path.folder = path.folder)

    # generate a variant.id vector ---------------------------------------------
    variant.id.select <- res$markers.meta$VARIANT_ID
    # Check filter
    # SeqArray::seqGetFilter(gdsfile = res$vcf.connection)
    # SeqArray::seqResetFilter(object = res$vcf.connection, verbose = FALSE)


    # SNP LD -------------------------------------------------------------------
    if (!is.null(filter.short.ld)) {
      filter.ld <- snp_ld(
        data = res$vcf.connection,
        snp.ld = filter.short.ld,
        maf.data = NULL,
        ld.threshold = filter.long.ld,
        parallel.core = parallel.core,
        filename = NULL,
        long.ld.missing = long.ld.missing
      )
      # names(filter.ld)
      blacklist.snp.ld <- res$markers.meta %>%
        dplyr::filter(!MARKERS %in% filter.ld$whitelist.snp.ld$MARKERS) %>%
        dplyr::select(MARKERS)

      if (nrow(blacklist.snp.ld) > 0) {
        message("Number of markers blacklisted based LD: ", nrow(blacklist.snp.ld))
        # update the blacklist
        res$blacklist.markers <- dplyr::bind_rows(
          res$blacklist.markers,
          dplyr::mutate(blacklist.snp.ld, FILTER = "filter LD (short/long)"))
        # check <- res$blacklist.markers

        # Update markers.meta
        res$markers.meta %<>%
          dplyr::filter(MARKERS %in% filter.ld$whitelist.snp.ld$MARKERS)
        # check <- res$markers.meta

        # Update GDS
        gdsfmt::add.gdsn(
          node = radiator.gds,
          name = "markers.meta",
          val = res$markers.meta,
          replace = TRUE,
          compress = "ZIP_RA",
          closezip = TRUE)

        SeqArray::seqSetFilter(object = res$vcf.connection,
                               variant.id = res$markers.meta$VARIANT_ID,
                               verbose = FALSE)
      }
      blacklist.snp.ld <- filter.ld <- NULL

      # make sure GDS is updated with sample -------------------------------------
      SeqArray::seqSetFilter(object = res$vcf.connection,
                             sample.id = strata$INDIVIDUALS,
                             # action = "set",
                             verbose = FALSE)
    }

    # blacklisted markers ------------------------------------------------------

    # Stats --------------------------------------------------------------------
    res$n.individuals <- length(SeqArray::seqGetData(res$vcf.connection, "sample.id"))
    res$n.markers <- length(SeqArray::seqGetData(res$vcf.connection, "variant.id"))

    # # vcf.sum <- SeqArray::seqSummary(res$vcf.connection, verbose = FALSE)
    # res$n.markers <-  vcf.sum$num.variant
    # res$n.individuals <- vcf.sum$num.sample
    # vcf.sum <- NULL

    if (vcf.stats) {
      res$n.chromosome <- dplyr::n_distinct(res$markers.meta$CHROM)
      res$n.locus <- dplyr::n_distinct(res$markers.meta$LOCUS)
      res$stats$ind.missing <- round(mean(res$individuals$MISSING_PROP[res$individuals$FILTER_INDIVIDUALS_MISSING], na.rm = TRUE), 2)
      res$stats$ind.cov.total <- round(mean(res$individuals$COVERAGE_TOTAL[res$individuals$FILTER_INDIVIDUALS_MISSING], na.rm = TRUE), 0)
      res$stats$ind.cov.mean <- round(mean(res$individuals$COVERAGE_MEAN[res$individuals$FILTER_INDIVIDUALS_MISSING], na.rm = TRUE), 0)
      res$stats$markers.missing <- round(mean(res$markers.meta$MISSING_PROP, na.rm = TRUE), 2)
      res$stats$markers.cov <- round(mean(res$markers.meta$COVERAGE_MEAN, na.rm = TRUE), 0) # same as above because NA...
    }
    res$filename <- filename

    if (verbose) {
      if (vcf.stats) {
        message("\n\nMissing data (averaged): ")
        message("    markers: ", res$stats$markers.missing)
        message("    individuals: ", res$stats$ind.missing)
        message("\n\nCoverage info:")
        message("    individuals mean read depth: ", res$stats$ind.cov.total)
        message("    individuals mean genotype coverage: ", res$stats$ind.cov.mean)
        message("    markers mean coverage: ", res$stats$markers.cov)
        message("\n\nNumber of chromosome/contig/scaffold: ", res$n.chromosome)
        message("Number of locus: ", res$n.locus)
      }
      message("Number of markers: ", res$n.markers)
      message("Number of individuals: ", res$n.individuals)
    }


  }#End stats

  # prepare data for the SeqVarTools package
  # seqOptimize("tmp.gds", target="by.sample")
  timing.import <- proc.time() - timing.import
  if (verbose) message("\nWorking time: ", round(timing.import[[3]]), " sec\n")
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
extract_coverage <- function(data = NULL) {#, variant.id.select = NULL) {
  # data <- res$vcf.connection
  coverage.info <- list()
  depth <- SeqArray::seqGetData(data, "annotation/format/DP")
  # test <- depth$data
  coverage.info$ind.cov.tot <- as.integer(round(rowSums(x = depth$data, na.rm = TRUE, dims = 1L), 0))
  coverage.info$ind.cov.mean <- as.integer(round(rowMeans(x = depth$data, na.rm = TRUE, dims = 1L), 0))
  coverage.info$markers.mean <- as.integer(round(colMeans(x = depth$data, na.rm = TRUE, dims = 1L), 0))
  # coverage.info$ind.cov.tot <- as.integer(round(rowSums(x = depth$data[,variant.id.select], na.rm = TRUE, dims = 1L), 0))
  # coverage.info$ind.cov.mean <- as.integer(round(rowMeans(x = depth$data[,variant.id.select], na.rm = TRUE, dims = 1L), 0))
  return(coverage.info)
}#End extract_coverage

# boxplot of stats
#' @title boxplot_stats
#' @description Generate a boxplot
#' @rdname boxplot_stats
#' @keywords internal
#' @export
boxplot_stats <- function(data,
                          title,
                          x.axis.title = NULL,
                          y.axis.title,
                          facet.columns = FALSE,
                          facet.rows = FALSE,
                          bp.filename = NULL,
                          path.folder = NULL) {
  # data <- test
  # x.axis.title = NULL
  # x.axis.title <- "SNP position on the read groupings"
  # title <- "Individual's QC stats"
  # title <- "Impact of filter on the SNP position in the read"
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
      title = title)

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
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, family = "Helvetica",
                                         face = "bold", hjust = 0.5),
      legend.position = "none",
      axis.title.y = element.text,
      axis.text.y = element.text
    ) +
    ggplot2::theme_bw()

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

  if (facet.columns) {
    fig.boxplot <- fig.boxplot + ggplot2::facet_grid(GROUP ~ ., scales = "free")
    n.facet <- n.group * 2
    width <- 10
    height <- 10 + (4 * n.group)
  }

  # else {
  # width <-  10 + (5 * n.group) + 1
  # height <-  10
  # }

  if (facet.rows) {
    fig.boxplot <- fig.boxplot + ggplot2::facet_grid(FACET_ROWS ~ ., scales = "free")
    n.facet <- n.group * 2
    width <- 10
    height <- 10 + (4 * n.group)
  }

  if (!facet.rows && !facet.columns) {
    width <-  10 + (5 * n.group) + 1
    height <-  10
  }

  print(fig.boxplot)
  if (!is.null(bp.filename)) {
    suppressMessages(ggplot2::ggsave(
      filename = file.path(path.folder, bp.filename),
      plot = fig.boxplot,
      width = width,
      height = height,
      dpi = 300, units = "cm", useDingbats = FALSE))
  }
  return(fig.boxplot)
}#Endboxplot_stats

# missingness per markers per pop
#' @title missing_per_pop
#' @description Generate missingness per markers per pop
#' @rdname missing_per_pop
#' @keywords internal
#' @export
missing_per_pop <- function(
  vcf.connection,
  strata,
  parallel.core = parallel::detectCores() - 1
) {
  # vcf.connection <- res$vcf.connection

  missing_pop <- function(
    id.select, vcf.connection,
    parallel.core = parallel::detectCores() - 1
  ) {
    SeqArray::seqSetFilter(object = vcf.connection,
                           sample.id = id.select$INDIVIDUALS,
                           action = "set",
                           verbose = FALSE)
    res <- tibble::tibble(
      MARKERS = gdsfmt::read.gdsn(gdsfmt::index.gdsn(node = vcf.connection, path = "radiator/markers.meta/MARKERS")),
      MISSING_PROP = SeqArray::seqMissing(
        gdsfile = vcf.connection,
        per.variant = TRUE, .progress = TRUE, parallel = parallel.core))
    SeqArray::seqResetFilter(
      object = vcf.connection, sample = TRUE, variant = FALSE, verbose = FALSE)
    return(res)
  }#End missing_pop

  res <- strata %>%
    dplyr::group_by(STRATA) %>%
    tidyr::nest(data = ., .key = id.select) %>%
    dplyr::mutate(MISSING_POP = purrr::map(
      .x = .$id.select,
      .f = missing_pop,
      vcf.connection = vcf.connection,
      parallel.core = parallel.core)) %>% #102
    tidyr::unnest(data = ., MISSING_POP) %>%#114
    #406 MB ... not very storage efficient by handy:
    dplyr::mutate(STRATA = as.character(STRATA)) %>% # less space
    dplyr::group_by(MARKERS) %>%
    tidyr::nest(data = ., .key = MISSING_POP)
  return(res)
}#End missing_per_pop
