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
#'
#'
#' \strong{... :dot dot dot arguments, to pass further argument to the function}
#' x arguments are available:
#' \enumerate{
#' \item whitelist.markers
#' \item blacklist.id
#' \item pop.select
#' \item pop.levels
#' \item pop.labels
#' \item snp.read.position.filter
#' \item mac.threshold
#' }
#' Documentation for these arguments is detailed
#' in \code{\link[radiator]{tidy_genomic_data}}.

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
  # vcf = "test.vcf" # data = "populations.snps.vcf"
  # strata <- "StLaw_popmap_thierry.tsv"
  # filename <- NULL
  # vcf.stats <- TRUE
  # parallel.core <- parallel::detectCores() - 1
  # verbose <- TRUE
  # snp.read.position.filter <- c("outliers", "q75", "iqr")
  # mac.threshold <- 4
  # blacklist.id = "blacklist.id.missing.50.tsv"
  # pop.select = NULL
  # pop.levels = NULL
  # pop.labels = NULL
  # whitelist.markers = NULL
  # keep.gds <- TRUE

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
  want <- c("whitelist.markers", "snp.read.position.filter", "mac.threshold",
            "blacklist.id", "pop.select", "pop.levels", "pop.labels", "keep.gds")
  unknowned_param <- setdiff(names(dotslist), want)

  if (length(unknowned_param) > 0) {
    stop("Unknowned \"...\" parameters ",
         stringi::stri_join(unknowned_param, collapse = " "))
  }

  radiator.dots <- dotslist[names(dotslist) %in% want]
  whitelist.markers <- radiator.dots[["whitelist.markers"]]
  snp.read.position.filter <- radiator.dots[["snp.read.position.filter"]]
  mac.threshold <- radiator.dots[["mac.threshold"]]
  blacklist.id <- radiator.dots[["blacklist.id"]]
  pop.select <- radiator.dots[["pop.select"]]
  pop.levels <- radiator.dots[["pop.levels"]]
  pop.labels <- radiator.dots[["pop.labels"]]
  keep.gds <- radiator.dots[["keep.gds"]]

  # useful outside this function
  if (is.null(keep.gds)) keep.gds <- TRUE

  if (!is.null(snp.read.position.filter)) {
    snp.read.position.filter <- match.arg(
      arg = snp.read.position.filter,
      choices = c("outliers", "iqr", "q75"),
      several.ok = TRUE)
  }

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
    blacklist.markers <- stringi::stri_join("blacklist.markers_", file.date, ".tsv")
  } else {
    ind.file <- stringi::stri_join(filename, "_vcf_individuals_info_", file.date, ".tsv")
    markers.file <- stringi::stri_join(filename, "_vcf_markers_metadata_", file.date, ".tsv")
    blacklist.markers <- stringi::stri_join(filename, "_blacklist.markers_", file.date, ".tsv")
    filename.problem <- file.exists(filename)
    if (filename.problem) {
      filename <- stringi::stri_join(filename, "_", file.date, ".gds")
    } else {
      filename <- stringi::stri_join(filename, ".gds")
    }
  }

  # Read vcf -------------------------------------------------------------------
  timing.vcf <- proc.time()

  # Check for bad header generated by stacks
  check.header <- SeqArray::seqVCF_Header(vcf)

  if (check.header$format$Number[check.header$format$ID == "AD"] == 1) {
    check.header$format$Number[check.header$format$ID == "AD"] <- "."
  }

  # remove GL info in stacks <v.2
  # Check stacks version (used below with metadata)
  check.source <- check.header$header$value[check.header$header$id == "source"]
  is.stacks <- stringi::stri_detect_fixed(str = check.source, pattern = "Stacks")
  if (is.stacks) {
    stacks.2 <- keep.stacks.gl <- stringi::stri_detect_fixed(str = check.source, pattern = "Stacks v2")
    if (!keep.stacks.gl) {
      check.header$format <- dplyr::filter(check.header$format, ID != "GL")
    }
  } else {
    stacks.2 <- FALSE
  }

  res$vcf.connection <- SeqArray::seqVCF2GDS(
    vcf.fn = vcf,
    out.fn = filename,
    parallel = parallel.core,
    storage.option = "ZIP_RA",
    verbose = FALSE, header = check.header
    ) %>%
    SeqArray::seqOpen(gds.fn = ., readonly = FALSE)
  check.header <- NULL

  if (verbose) message("\nconversion timing: ", round((proc.time() - timing.vcf)[[3]]), " sec")

  if (verbose && keep.gds) {
    message("\nSeqArray GDS file generated: ", filename)
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

  # keep id in VCF and clean
  res$individuals <- tibble::tibble(
    INDIVIDUALS_VCF = SeqArray::seqGetData(res$vcf.connection, "sample.id")) %>%
    dplyr::mutate(INDIVIDUALS = radiator::clean_ind_names(INDIVIDUALS_VCF))

  # test <- res$individuals

  # replace id in VCF
  gdsfmt::add.gdsn(
    node = res$vcf.connection,
    name = "sample.id",
    val = res$individuals$INDIVIDUALS,
    replace = TRUE,
    compress = "ZIP_RA",
    closezip = TRUE)

  # potentially filter ids based on strata
  res$individuals %<>%
    dplyr::filter(INDIVIDUALS %in% strata$INDIVIDUALS) %>%
    dplyr::left_join(strata, by = "INDIVIDUALS")

  SeqArray::seqSetFilter(object = res$vcf.connection,
                         sample.id = res$individuals$INDIVIDUALS,
                         action = "set",
                         verbose = FALSE)
  # SeqArray::seqResetFilter(object = res$vcf.connection, sample = TRUE, variant = TRUE, verbose = TRUE)

  # Markers metadata  -----------------------------------------------------------
  if (vcf.stats) {
    if (verbose) message("\nradiator is working on the file ...")
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
    }


    # test <- res$markers.meta

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

    # Generate MARKERS column and fix types
    res$markers.meta %<>%
      dplyr::mutate_at(.tbl = .,
                       .vars = c("CHROM", "LOCUS", "POS"),
                       .funs = radiator::clean_markers_names) %>%
      dplyr::mutate(
        MARKERS = stringi::stri_join(CHROM, LOCUS, POS, sep = "__"),
        REF = SeqArray::seqGetData(gdsfile = res$vcf.connection, var.name = "$ref"),
        ALT = SeqArray::seqGetData(gdsfile = res$vcf.connection, var.name = "$alt")
      )

    # test <- res$markers.meta

    # Scan and filter with FILTER column ---------------------------------------
    res$markers.meta$FILTER <- SeqArray::seqGetData(res$vcf.connection, "annotation/filter")
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

    # bi- or multi-alllelic VCF --------------------------------------------------
    alt.num <- max(unique(SeqArray::seqNumAllele(gdsfile = res$vcf.connection))) - 1

    if (alt.num > 1) {
      res$biallelic <- FALSE
      message("VCF is multi-allelic")
    } else {
      res$biallelic <- TRUE
      message("VCF is biallelic")
    }

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

    # Markers stats ------------------------------------------------------------
    if (verbose) message("Updating markers metadata and stats")
    n.markers <- dplyr::n_distinct(res$markers.meta$VARIANT_ID)
    want <- c("VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS",
              "COL", "REF", "ALT", "MISSING_PROP", "REF_COUNT", "MAC", "SNP_PER_LOCUS")
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
        dplyr::select(dplyr::one_of(want))
    )
    # test <- res$markers.meta

    # when low MAC is remove that number drops...
    if (!is.null(mac.threshold)) {
      res$markers.meta %<>%
        dplyr::mutate(
          MAC_FILTER = dplyr::if_else(MAC >= mac.threshold, TRUE, FALSE))
    } else {
      res$markers.meta %<>%
        dplyr::mutate(MAC_FILTER = TRUE)
    }
    # test <- res$markers.meta

    # position of the SNPs on the read -----------------------------------------
    if (tibble::has_name(res$markers.meta, "COL")) {
      if (verbose) message("Generating SNP position on read stats")
      res$stats$snp.col.stats <- tibble_stats(
        x = dplyr::filter(res$markers.meta, MAC_FILTER) %>% # using filtered dataset
          dplyr::distinct(MARKERS,COL) %$%
          COL,
        group = "snp position on read") #%>% dplyr::select(-OUTLIERS_LOW)
      snp.col.iqr.threshold <- c(res$stats$snp.col.stats$Q25, res$stats$snp.col.stats$Q75)

      # outliers
      if ("outliers" %in% snp.read.position.filter) {
        res$markers.meta %<>%
          dplyr::mutate(
            SNP_POS_READ_OUTLIERS = dplyr::if_else(
              COL < res$stats$snp.col.stats$OUTLIERS_HIGH, TRUE, FALSE))
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
      res$stats$snp.col.stats %<>% dplyr::bind_rows(
        tibble_stats(
          x = unique.col <- res$markers.meta %>%
            dplyr::filter(MAC_FILTER) %>% # mac filter
            dplyr::filter(SNP_POS_READ_OUTLIERS) %>% # outlier snp on read filter
            dplyr::filter(SNP_POS_READ_Q75) %>% # snp below Q75
            dplyr::filter(SNP_POS_READ_IQR) %>% # focus on snp within IQR range
            dplyr::distinct(MARKERS, COL) %$% COL,
          group = "snp position on read filtered")) #%>% dplyr::select(-OUTLIERS_LOW))
      snp.col.iqr.threshold <- NULL
      # test <- res$stats$snp.col.stats
      # Generate box plot
      res$figures$snp.pos.read.fig <- boxplot_stats(
        data = res$stats$snp.col.stats,
        title = "Impact of filter on the SNP position on the read",
        x.axis.title = "SNP position on the read groupings",
        y.axis.title = "SNP position (base pair)",
        bp.filename = "vcf.snp.position.read.pdf")
    }
    # number of SNPs per locus -------------------------------------------------
    # before and after filtering
    res$stats$snp.locus.stats <- tibble_stats(
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
    res$figures$snp.per.locus.fig <- boxplot_stats(
      data = res$stats$snp.locus.stats,
      title = "Impact of filters on the number of SNPs per locus",
      x.axis.title = "groupings",
      y.axis.title = "Number of SNPs per locus",
      bp.filename = "vcf.number.snp.per.locus.pdf")

    # generate a variant.id vector ---------------------------------------------
    variant.id.select <- res$markers.meta %>%
      dplyr::filter(MAC_FILTER) %>% # mac filter
      dplyr::filter(SNP_POS_READ_OUTLIERS) %>% # outlier snp on read filter
      dplyr::filter(SNP_POS_READ_Q75) %>% # snp below Q75
      dplyr::filter(SNP_POS_READ_IQR) %>%
      # dplyr::select(LOCUS, VARIANT_ID) %>%
      # dplyr::distinct(LOCUS, .keep_all = TRUE) %>%
      dplyr::select(VARIANT_ID) %>%
      purrr::flatten_int(.)

    # Apply the filter on the gds ...
    SeqArray::seqSetFilter(object = res$vcf.connection,
                           variant.id = variant.id.select, verbose = FALSE)
    # Check filter
    # SeqArray::seqGetFilter(gdsfile = res$vcf.connection)
    # SeqArray::seqResetFilter(object = res$vcf.connection, verbose = FALSE)

    if (verbose) message("Generating coverage stats")
    coverage.info <- extract_coverage(data = res$vcf.connection)
    # variant.id.select = variant.id.select) # old version discard if no longer necessary

    # blacklisted markers ------------------------------------------------------
    blacklisted.markers <- nrow(res$markers.meta) - length(variant.id.select)

    if (blacklisted.markers > 0) {
      if (verbose) message("Number of blacklisted markers: ", blacklisted.markers)
      res$blacklist.markers <- dplyr::filter(
        res$markers.meta, !VARIANT_ID %in% variant.id.select) %>%
        readr::write_tsv(x = ., path = blacklist.markers)
      res$markers.meta %<>% dplyr::filter(VARIANT_ID %in% variant.id.select)
    }

    #Coverage: MARKERS
    res$markers.meta %<>% dplyr::mutate(COVERAGE_MEAN = coverage.info$markers.mean) %>%
      readr::write_tsv(x = ., path = markers.file)
    # dplyr::left_join(
    #   tibble::tibble(VARIANT_ID = variant.id.select,
    #                  COVERAGE_MEAN = coverage.info$markers.mean
    #   ), by = "VARIANT_ID")

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
      bp.filename = "vcf.individuals.qc.pdf")

    # Stats ----------------------------------------------------------------------
    res$n.chromosome <- dplyr::n_distinct(res$markers.meta$CHROM)
    res$n.locus <- dplyr::n_distinct(res$markers.meta$LOCUS)
    res$n.markers <- dplyr::n_distinct(res$markers.meta$MARKERS)
    res$n.individuals <- nrow(res$individuals)
    res$stats$ind.missing <- round(mean(res$individuals$MISSING_PROP, na.rm = TRUE), 2)
    res$stats$ind.cov.total <- round(mean(res$individuals$COVERAGE_TOTAL, na.rm = TRUE), 0)
    res$stats$ind.cov.mean <- round(mean(res$individuals$COVERAGE_MEAN, na.rm = TRUE), 0)
    res$stats$markers.missing <- round(mean(res$markers.meta$MISSING_PROP, na.rm = TRUE), 2)
    res$stats$markers.cov <- round(mean(res$markers.meta$COVERAGE_MEAN, na.rm = TRUE), 0) # same as above because NA...

    res$filename <- filename

    if (verbose) {
      message("\n\nMissing data (averaged): ")
      message("    markers: ", res$stats$markers.missing)
      message("    individuals: ", res$stats$ind.missing)
      message("\n\nCoverage info:")
      message("    individuals mean read depth: ", res$stats$ind.cov.total)
      message("    individuals mean genotype coverage: ", res$stats$ind.cov.mean)
      message("    markers mean coverage: ", res$stats$markers.cov)
      message("\n\nNumber of chromosome/contig/scaffold: ", res$n.chromosome)
      message("Number of locus: ", res$n.locus)
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
                          bp.filename = NULL) {
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


  # fig.boxplot <- fig.boxplot +
  #   ggplot2::theme(
  #     plot.title = ggplot2::element_text(size = 12, family = "Helvetica",
  #                                        face = "bold", hjust = 0.5),
  #     legend.position = "none",
  #     axis.title.y = element.text,
  #     axis.text.y = element.text
  #   ) +
  #   ggplot2::theme_bw()

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
      filename = bp.filename,
      plot = fig.boxplot,
      width = width,
      height = height,
      dpi = 300, units = "cm", useDingbats = FALSE))
  }
  return(fig.boxplot)
}#Endboxplot_stats
