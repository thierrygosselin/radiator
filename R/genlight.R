# tidy_genlight ----------------------------------------------------------------
#' @name tidy_genlight
#' @title Tidy a genlight object to a tidy dataframe and/or GDS object/file
#' @description Tidy genlight object to a tidy dataframe and/or GDS object/file.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data (path or object) A genlight object in the global environment or
#' path to a genlight file that will be open with \code{readRDS}.

#' @inheritParams radiator_common_arguments

#' @param tidy (logical) Generate a tidy dataset.
#' Default: \code{tidy = TRUE}.

#' @param gds (optional, logical) To write a radiator gds object.
#' Default: \code{gds = TRUE}.

#' @param write (optional, logical) To write in the working directory the tidy
#' data. The file is written with \code{radiator_genlight_DATE@TIME.rad}.
#' Default: \code{write = FALSE}.

#' @export
#' @rdname tidy_genlight

#' @references Jombart T (2008) adegenet: a R package for the multivariate
#' analysis of genetic markers. Bioinformatics, 24, 1403-1405.
#' @references Jombart T, Ahmed I (2011) adegenet 1.3-1:
#' new tools for the analysis of genome-wide SNP data.
#' Bioinformatics, 27, 3070-3071.

#' @details
#' A string of the same dimension is generated when genlight:
#' \enumerate{
#' \item \code{is.null(genlight@pop)}: pop will be integrated
#' in the tidy dataset.
#' \item \code{is.null(data@chromosome)}: CHROM1 will be integrated
#' in the tidy dataset.
#' \item \code{is.null(data@loc.names)}: LOCUS1 to dim(genlight)[2]
#' will be integrated in the tidy dataset.
#' \item \code{is.null(data@position)}: an integer string of
#' length = dim(genlight)[2] will be integrated in the tidy dataset.
#' }
#'
#' \strong{Note: that if all CHROM, LOCUS and POS is missing the function will be terminated}


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


tidy_genlight <- function(
  data,
  tidy = TRUE,
  gds = TRUE,
  write = FALSE,
  verbose = FALSE,
  parallel.core = parallel::detectCores() - 1
) {
  # Test
  # data = "radiator_genlight_20191211@1836.RData"
  # tidy = TRUE
  # gds = TRUE
  # write = FALSE
  # verbose = TRUE
  # parallel.core = 12L


  # Package requirement --------------------------------------------------------
  radiator_packages_dep(package = "adegenet")

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("data argument required")

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) data <- readRDS(data)
  if (class(data)[1] != "genlight") rlang::abort("Input is not a genlight object")

  if (verbose) message("genlight info:")
  # strata ?
  strata <- tibble::tibble(INDIVIDUALS = data@ind.names)
  if (is.null(data@pop)) {
    if (verbose) message("    strata: no")
    if (verbose) message("    'pop' will be added")
    strata %<>% dplyr::mutate(POP_ID = "pop")
  } else {
    if (verbose) message("    strata: yes")
    strata$POP_ID = data@pop
  }

  n.markers <- dim(data)[2]
  n.ind <- nrow(strata)

  # Chromosome ?
  if (is.null(data@chromosome)) {
    if (verbose) message("    Chromosome/contig/scaffold: no")
    data@chromosome <- factor(rep("CHROM1", n.markers))
    chrom.info <- FALSE
  } else {
    if (verbose) message("    Chromosome/contig/scaffold: yes")
    chrom.info <- TRUE
  }

  # Locus ?
  if (is.null(data@loc.names)) {
    if (verbose) message("    Locus: no")
    locus.info <- FALSE
    data@loc.names <- stringi::stri_join("LOCUS", seq(from = 1, to = n.markers, by = 1))
  } else {
    if (verbose) message("    Locus: yes")
    locus.info <- TRUE
  }

  # POS ?
  if (is.null(data@position)) {
    if (verbose) message("    POS: no")
    pos.info <- FALSE
    data@position <- rlang::as_integer(seq(from = 1, to = n.markers, by = 1))
  } else {
    if (verbose) message("    POS: yes")
    pos.info <- TRUE
  }


  if (!chrom.info && !locus.info && !pos.info) {
    rlang::abort("Tidying the genlight requires at least one of these 3 markers metadata:
       CHROM (genlight@chromosome), LOCUS (genlight@loc.names) or POS (genlight@position)")
  }


  # markers
  markers <- tibble::tibble(
    CHROM = data@chromosome,#adegenet::chromosome(data),
    LOCUS = data@loc.names,#adegenet::locNames(data),
    POS = data@position#adegenet::position(data)
  ) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), .fns = as.character)) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), .fns = radiator::clean_markers_names)) %>%
    tidyr::unite(data = ., col = MARKERS, CHROM, LOCUS, POS, sep = "__", remove = FALSE) %>%
    dplyr::select(MARKERS, CHROM, LOCUS, POS)

  # Nuc info
  nuc.data <- data@loc.all
  if (!is.null(nuc.data)) {
    markers %<>%
      dplyr::mutate(
        REF = stringi::stri_sub(str = nuc.data, from = 1, to = 1),
        ALT = stringi::stri_sub(str = nuc.data, from = 3, to = 3)
      )
  }

  if (gds) tidy <- TRUE

  if (tidy) {
    if (write) {
      filename.temp <- generate_filename(extension = "rad")
      filename.short <- filename.temp$filename.short
      filename.genlight <- filename.temp$filename
    }

    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "REF", "ALT","POP_ID", "INDIVIDUALS",
              "GT_VCF", "GT_BIN", "GT")

    if (verbose) message("Generating tidy data...")
    tidy.data <- suppressWarnings(
      data.frame(data) %>%
        magrittr::set_colnames(x = ., value = markers$MARKERS) %>%
        tibble::add_column(.data = ., INDIVIDUALS = rownames(.), .before = 1) %>%
        data.table::as.data.table(.) %>%
        data.table::melt.data.table(
          data = .,
          id.vars = "INDIVIDUALS",
          variable.name = "MARKERS",
          value.name = "GT_BIN"
        ) %>%
        tibble::as_tibble(.) %>%
        dplyr::full_join(markers, by = "MARKERS") %>%
        dplyr::full_join(strata, by =  "INDIVIDUALS") %>%
        dplyr::mutate(
          INDIVIDUALS = radiator::clean_ind_names(INDIVIDUALS),
          POP_ID = radiator::clean_pop_names(POP_ID)
        ) %>%
        dplyr::arrange(MARKERS, POP_ID, INDIVIDUALS) %>%
        dplyr::select(dplyr::one_of(want))
    ) %>%
      radiator::calibrate_alleles(
        data = .,
        biallelic = TRUE,
        parallel.core = parallel.core,
        verbose = verbose
      ) %$%
      input
    if (write) {
      radiator::write_rad(data = tidy.data, path = filename.genlight)
      if (verbose) message("File written: ", filename.short)
    }
  }#End tidy genlight

  if (gds) {
    # generate the GDS --------------------------------------------------------------
    # markers %<>% dplyr::mutate(VARIANT_ID = as.integer(factor(MARKERS)))
    gds.filename <- radiator_gds(
      data.source = "genlight",
      genotypes = gt2array(
        gt.bin = tidy.data$GT_BIN,
        n.ind = n.ind,
        n.snp = n.markers
      ),
      # genotypes = tibble::as_tibble(data.frame(data) %>% t) %>%
      #   tibble::add_column(.data = ., MARKERS = markers$MARKERS, .before = 1) %>%
      #   dplyr::arrange(MARKERS) %>%
      #   tibble::column_to_rownames(.data = ., var = "MARKERS"),
      strata = dplyr::rename(strata, STRATA = POP_ID),
      biallelic = TRUE,
      markers.meta = markers,
      filename = NULL,
      verbose = verbose
    )
    # if (verbose) message("Written: GDS filename: ", gds.filename)
  }# End gds genlight

  if (tidy) {
    return(tidy.data)
  } else {
    message("returning GDS filename")
    return(gds.filename)
  }
} # End tidy_genlight

# write_genlight ----------------------------------------------------------------
#' @name write_genlight
#' @title Write a genlight object from a tidy data frame or GDS file or object.
#' @description Write a genlight object from a tidy data frame or GDS file or object.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @inheritParams radiator_common_arguments

#' @param biallelic (logical, optional) If you already know that the data is
#' biallelic use this argument to speed up the function.
#' Default: \code{biallelic = TRUE}.

#' @param write (logical, optional) To write in the working directory the genlight
#' object. The file is written with \code{radiator_genlight_DATE@TIME.RData} and
#' can be open with load or readRDS.
#' Default: \code{write = FALSE}.

#' @param dartr (logical, optional) For non-dartR users who wants to have a genlight
#' object ready for the package. This option transfer or generates:
#' \code{CALL_RATE, AVG_COUNT_REF, AVG_COUNT_SNP, REP_AVG,
#' ONE_RATIO_REF, ONE_RATIO_SNP}. These markers metadata are stored into
#' the genlight slot: \code{genlight.obj@other$loc.metrics}.
#' \strong{Use the radiator generated GDS data for best result}.
#' Default: \code{dartr = FALSE}.

#' @export
#' @rdname write_genlight

#' @references Jombart T (2008) adegenet: a R package for the multivariate
#' analysis of genetic markers. Bioinformatics, 24, 1403-1405.
#' @references Jombart T, Ahmed I (2011) adegenet 1.3-1:
#' new tools for the analysis of genome-wide SNP data.
#' Bioinformatics, 27, 3070-3071.

#' @examples
#' \dontrun{
#' # With defaults:
#' turtle <- radiator::write_genlight(data = "my.radiator.gds.rad")
#'
#' # Write gl object in directory:
#' turtle <- radiator::write_genlight(data = "my.radiator.gds.rad", write = TRUE)
#'
#' # Generate a dartR ready genlight and verbose = TRUE:
#' turtle <- radiator::write_genlight(
#'     data = "my.radiator.gds.rad",
#'     write = TRUE,
#'     dartr = TRUE,
#'     verbose = TRUE
#'  )
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_genlight <- function(
  data,
  biallelic = TRUE,
  write = FALSE,
  dartr = FALSE,
  verbose = FALSE,
  parallel.core = parallel::detectCores() - 2
) {


  # NOTE: Make sure it's the same levels.... id, markers etc...
  # TEST
  # biallelic = TRUE
  # write = TRUE
  # dartr = TRUE
  # verbose = TRUE
  # parallel.core = parallel::detectCores() - 1

  # Checking for missing and/or default arguments ------------------------------
  radiator_packages_dep(package = "adegenet")
  if (missing(data)) rlang::abort("Input file missing")
  if (verbose) message("Generating genlight...")

  # File type detection---------------------------------------------------------
  data.type <- radiator::detect_genomic_format(data)

  # Import data ---------------------------------------------------------------
  if (data.type %in% c("SeqVarGDSClass", "gds.file")) {
    # Check that SeqVarTools is installed (it requires automatically: SeqArray and gdsfmt)
    radiator_packages_dep(package = "SeqVarTools", cran = FALSE, bioc = TRUE)

    if (data.type == "gds.file") {
      data <- radiator::read_rad(data, verbose = verbose)
    }
    biallelic <- radiator::detect_biallelic_markers(data)# faster with GDS
    if (!biallelic) rlang::abort("genlight object requires biallelic genotypes")

    data.source <- radiator::extract_data_source(data)
    # markers.meta <- extract_markers_metadata(gds = data, whitelist = TRUE)
    data.bk <- data
    # data.bk -> data
    if (dartr && !any(unique(c("1row", "2rows") %in% data.source))) {
      # counts data and data with read depth alleles depth info...
      if (!"counts" %in% data.source) {
        genotypes.meta <- radiator::extract_genotypes_metadata(gds = data, whitelist = TRUE)
      } else {
        genotypes.meta <- NULL
      }
      data <- gds2tidy(gds = data.bk,
                       wide = FALSE,
                       # markers.meta.select = "MARKERS",
                       parallel.core = parallel.core,
                       calibrate.alleles = FALSE)
      markers.levels <- unique(data$MARKERS)
      ind.levels <- unique(data$INDIVIDUALS)

      if (!"counts" %in% data.source) {
        genotypes.meta <- extract_coverage(
          gds = data.bk,
          ind = FALSE,
          depth.tibble = TRUE,
          parallel.core = parallel.core
        )
        suppressWarnings(
          data %<>%
            dplyr::left_join(genotypes.meta, by = c("INDIVIDUALS", "MARKERS"))
        )
      }

      genotypes.meta <- NULL
      want <- c("MARKERS", "CHROM", "LOCUS", "POS",
                "REF", "ALT",
                "CALL_RATE", "AVG_COUNT_REF", "AVG_COUNT_SNP", "REP_AVG",
                "ONE_RATIO_REF", "ONE_RATIO_SNP")
      suppressWarnings(
        markers.meta <- data %>%
          dplyr::select(dplyr::one_of(want)) %>%
          dplyr::distinct(MARKERS, .keep_all = TRUE)
      )
    } else {
      data <- radiator::extract_individuals_metadata(
        gds = data, ind.field.select = c("INDIVIDUALS", "STRATA"), whitelist = TRUE) %>%
        dplyr::bind_cols(
          SeqArray::seqGetData(
            gdsfile = data, var.name = "$dosage_alt") %>%
            magrittr::set_colnames(x = ., value = markers.meta$MARKERS) %>%
            tibble::as_tibble(.)
        ) %>%
        dplyr::rename(POP_ID = STRATA)
    }

    # dartR-----------------------------------------------------------------------
    if (dartr) {
      # That bits of code below generate what's nenessary for dartR
      if (verbose) message("Calculating read depth for each SNP\n")
      if ("dart" %in% data.source) {
        # if (any(unique(c("1row", "2rows") %in% data.source))) {
        markers.meta %<>%
          dplyr::mutate(
            N_IND = SeqArray::seqApply(
              gdsfile = data.bk,
              var.name = "$dosage_alt",
              FUN = function(g) length(g[!is.na(g)]),
              margin = "by.variant", as.is = "integer",
              parallel = parallel.core),
            rdepth = round(
              ((N_IND * ONE_RATIO_REF * AVG_COUNT_REF) +
                 (N_IND * ONE_RATIO_SNP * AVG_COUNT_SNP)) / N_IND
            )
          ) %>%
          dplyr::ungroup(.)
        data.bk <- NULL
      } else {
        if (rlang::has_name(data, "READ_DEPTH")) {
          # dart 2 rows counts...
          if (rlang::has_name(markers.meta, "AVG_COUNT_REF")) {

            #NOTE TO MYSELF: will have to keep REP_AVG from count somehow... to do
            # for now, I expect people will use this with filtered data so not important
            # to fill with 1
            if ("counts" %in% data.source) rep.avg <- markers.meta$REP_AVG
            not.wanted <- c("CALL_RATE", "AVG_COUNT_REF", "AVG_COUNT_SNP",
                            "REP_AVG", "ONE_RATIO_REF", "ONE_RATIO_SNP")
            markers.meta %<>% dplyr::select(-dplyr::one_of(not.wanted))
          }

          suppressWarnings(
            markers.meta %<>%
              dplyr::left_join(
                data %>%
                  dplyr::group_by(MARKERS) %>%
                  dplyr::summarise(
                    rdepth = mean(READ_DEPTH, na.rm = TRUE),
                    AVG_COUNT_REF = mean(ALLELE_REF_DEPTH, na.rm = TRUE),
                    AVG_COUNT_SNP = mean(ALLELE_ALT_DEPTH, na.rm = TRUE),
                    CALL_RATE = length(INDIVIDUALS[!is.na(GT_BIN)]) / length(INDIVIDUALS),
                    ONE_RATIO_REF = length(INDIVIDUALS[GT_BIN == 0]) + length(INDIVIDUALS[GT_BIN == 1]) / length(INDIVIDUALS),
                    ONE_RATIO_SNP = length(INDIVIDUALS[GT_BIN == 2]) + length(INDIVIDUALS[GT_BIN == 1]) / length(INDIVIDUALS)
                  ) %>%
                  dplyr::mutate(REP_AVG = 1L)
                , by = "MARKERS"
              ) %>%
              dplyr::ungroup(.)
          )
          genotypes.meta <- NULL
          #Note to myself: will have to remove duplicated info in code between genotypes.meta and data...
          if ("counts" %in% data.source) {
            markers.meta$REP_AVG <- rep.avg
          }
        }
      }

      # ~25 times slower
      # Calculate Read Depth (from Arthur Georges)
      # gl.obj@other$loc.metrics$rdepth <- array(NA, adegenet::nLoc(gl.obj))
      # for (i in 1:adegenet::nLoc(gl.obj)){
      #   called.ind <- round(adegenet::nInd(gl.obj) * gl.obj@other$loc.metrics$CallRate[i],0)
      #   ref.count <- called.ind * gl.obj@other$loc.metrics$OneRatioRef[i]
      #   alt.count <- called.ind * gl.obj@other$loc.metrics$OneRatioSnp[i]
      #   sum.count.ref <- ref.count * gl.obj@other$loc.metrics$AvgCountRef[i]
      #   sum.count.alt <- ref.count * gl.obj@other$loc.metrics$AvgCountSnp[i]
      #   gl.obj@other$loc.metrics$rdepth[i] <- round((sum.count.alt + sum.count.ref) / called.ind, 1)
      # }

      # gl.obj@other$loc.metrics$rdepth
      data %<>%
        dplyr::select(MARKERS, POP_ID, INDIVIDUALS, GT_BIN) %>%
        dplyr::mutate(GT_BIN = as.integer(GT_BIN)) %>%
        data.table::as.data.table(.) %>%
        data.table::dcast.data.table(
          data = .,
          formula = POP_ID + INDIVIDUALS ~ MARKERS,
          value.var = "GT_BIN"
        ) %>%
        tibble::as_tibble(.)
    }#End dartr

    data.type <- "tbl_df"
  } else {
    # Tidy data
    if (rlang::has_name(data, "STRATA")) data %<>% dplyr::rename(POP_ID = STRATA)
    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "POP_ID", "INDIVIDUALS",
              "REF", "ALT", "GT_VCF", "GT_BIN",
              "CALL_RATE", "AVG_COUNT_REF", "AVG_COUNT_SNP", "REP_AVG",
              "ONE_RATIO_REF", "ONE_RATIO_SNP")
    data <- suppressWarnings(
      radiator::tidy_wide(data = data, import.metadata = TRUE) %>%
        dplyr::select(dplyr::one_of(want)) %>%
        dplyr::arrange(POP_ID, INDIVIDUALS)
    )

    # Detect if biallelic data ---------------------------------------------------
    if (is.null(biallelic)) biallelic <- radiator::detect_biallelic_markers(data)
    if (!biallelic) rlang::abort("genlight object requires biallelic genotypes")
    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "REF", "ALT",
              "CALL_RATE", "AVG_COUNT_REF", "AVG_COUNT_SNP", "REP_AVG",
              "ONE_RATIO_REF", "ONE_RATIO_SNP")
    data %<>% dplyr::arrange(MARKERS, POP_ID, INDIVIDUALS)
    markers.meta <- suppressWarnings(
      dplyr::select(data, dplyr::one_of(want)) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
        separate_markers(
          data = .,
          sep = "__",
          markers.meta.all.only = TRUE,
          biallelic = TRUE,
          verbose = verbose)
    )

    if (!rlang::has_name(data, "GT_BIN") && rlang::has_name(data, "GT_VCF")) {
      data$GT_BIN <- stringi::stri_replace_all_fixed(
        str = data$GT_VCF,
        pattern = c("0/0", "1/1", "0/1", "1/0", "./."),
        replacement = c("0", "2", "1", "1", NA),
        vectorize_all = FALSE
      )
    }

    data %<>%
      dplyr::select(MARKERS, POP_ID, INDIVIDUALS, GT_BIN) %>%
      dplyr::mutate(GT_BIN = as.integer(GT_BIN)) %>%
      data.table::as.data.table(.) %>%
      data.table::dcast.data.table(
        data = .,
        formula = POP_ID + INDIVIDUALS ~ MARKERS,
        value.var = "GT_BIN"
      ) %>%
      tibble::as_tibble(.)
  }# End tidy data

  # Generate genlight
  if (length(markers.meta$MARKERS) > 10000) {
    parallel.core.temp <- parallel.core
  } else {
    parallel.core.temp <- FALSE
  }

  # longer
  # gl.obj2 <- methods::new(
  #   "genlight",
  #   gen = data[,-(1:2)],
  #   ploidy = 2,
  #   ind.names = data$INDIVIDUALS,
  #   chromosome = markers.meta$CHROM,
  #   loc.names = markers.meta$LOCUS,
  #   position = markers.meta$POS,
  #   pop = data$POP_ID,
  #   parallel = parallel.core.temp)
  # tictoc::toc()


  gl.obj <- methods::new(
    "genlight",
    data[,-(1:2)],
    parallel = parallel.core.temp
  )
  adegenet::indNames(gl.obj)   <- data$INDIVIDUALS
  adegenet::pop(gl.obj)        <- data$POP_ID
  adegenet::chromosome(gl.obj) <- markers.meta$CHROM
  adegenet::locNames(gl.obj)   <- markers.meta$LOCUS
  adegenet::position(gl.obj)   <- markers.meta$POS

  if (dartr) {
    gl.obj@other$loc.metrics <- markers.meta %>%
      dplyr::select(
        OneRatioRef = ONE_RATIO_REF,
        OneRatioSnp = ONE_RATIO_SNP,
        AvgCountRef = AVG_COUNT_REF,
        AvgCountSnp = AVG_COUNT_SNP,
        CallRate = CALL_RATE,
        RepAvg = REP_AVG,
        rdepth
      )
  }#End dartr


  # Check
  # gl.obj@n.loc
  # gl.obj@ind.names
  # gl.obj@chromosome
  # gl.obj@position
  # length(gl.obj@position)
  # gl.obj@loc.names
  # length(gl.obj@loc.names)
  # gl.obj@pop
  # gl.obj@strata
  # adegenet::nLoc(gl.obj)
  # adegenet::popNames(gl.obj)
  # adegenet::indNames(gl.obj)
  # adegenet::nPop(gl.obj)
  # adegenet::NA.posi(gl.obj)
  # names(gl.obj@other$loc.metrics)


  if (write) {
    filename.temp <- generate_filename(extension = "genlight")
    filename.short <- filename.temp$filename.short
    filename.genlight <- filename.temp$filename
    saveRDS(object = gl.obj, file = filename.genlight)
    if (verbose) message("File written: ", filename.short)
  }
  return(gl.obj)
} # End write_genlight
