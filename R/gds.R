# radiator gds constructor------------------------------------------------------
#' @title radiator gds constructor
#' @description helper function to construct a radiator gds
#' @name radiator_gds
#' @rdname radiator_gds
#' @keywords internal
#' @export
radiator_gds <- function(
  genotypes.df,
  strata = NULL,
  biallelic = TRUE,
  markers.meta = NULL,
  genotypes.meta = NULL,
  dp = NULL,
  ad = NULL,
  filename = NULL,
  source = NULL,
  open = FALSE,
  verbose = TRUE
) {
  # genotypes.df <- data

  if (!biallelic) rlang::abort("Biallelic data required")

  # Filename
  if (is.null(filename)) {
    filename.temp.gds <- generate_filename(
      name.shortcut = filename,
      extension = "gds.rad")
    filename.short.gds <- filename.temp.gds$filename.short
    filename.gds <- filename.temp.gds$filename
    filename.temp.gds <- NULL
    radiator.temp.file <- stringi::stri_join(filename.gds, "_radiator_temp")
  } else {
    filename.short.gds <- folder_short(filename)
    radiator.temp.file <- stringi::stri_join(filename, "_radiator_temp")
    filename.gds <- filename
  }

  if (rlang::has_name(genotypes.df, "MARKERS")) {
    variant.id <- as.integer(factor(order(unique(genotypes.df$MARKERS))))
    n.markers <- length(variant.id)
  } else {
    n.markers <- nrow(genotypes.df)
    variant.id <- rownames(genotypes.df)
  }

  # check if data is long or wide format...
  if (rlang::has_name(genotypes.df, "MARKERS") &&
      rlang::has_name(genotypes.df, "INDIVIDUALS")) {
    genotypes.df %<>%
      dplyr::ungroup(.) %>%
      dplyr::select(MARKERS, INDIVIDUALS, GT_BIN) %>%
      dplyr::arrange(MARKERS, INDIVIDUALS) %>%
      data.table::as.data.table(.) %>%
      data.table::dcast.data.table(
        data = .,
        formula = MARKERS ~ INDIVIDUALS,
        value.var = "GT_BIN"
      ) %>%
      dplyr::arrange(MARKERS) %>%
      dplyr::select(-MARKERS) %>%
      data.matrix(.) %>%
      magrittr::set_rownames(x = ., value = variant.id)
  }

  if (rlang::has_name(genotypes.df, "MARKERS")) {
    genotypes.df %<>%
      dplyr::select(-MARKERS)
    data.matrix(.) %>%
      magrittr::set_rownames(x = ., value = variant.id)
  }

  # Genotypes coding
  genotypes.df %<>%
    # 4 steps for SNPRelate genotype coding (change from ALT to REF dosage)
    magrittr::inset(is.na(.), 3L) %>%
    magrittr::inset(. == 0L, 9L) %>%
    magrittr::inset(. == 2L, 0L) %>%
    magrittr::inset(. == 9L, 2L)


  # STRATA
  if (is.null(strata)) {
    strata <- tibble::tibble(
      STRATA = rep("pop", length(colnames(genotypes.df))),
      INDIVIDUALS = colnames(genotypes.df)
    )
  } else {
    if (!rlang::has_name(strata, "POP_ID") && !rlang::has_name(strata, "STRATA")) {
      strata %<>% dplyr::mutate(STRATA = "pop")
    }
    if (!rlang::has_name(strata, "INDIVIDUALS")) {
      strata %<>% dplyr::mutate(INDIVIDUALS = colnames(genotypes.df))
    }
  }

  strata %<>% dplyr::distinct(INDIVIDUALS, .keep_all = TRUE)

  # MARKERS META
  if (is.null(markers.meta)) {
    if (rlang::has_name(genotypes.df, "MARKERS")) {
      markers.meta <- unique(genotypes.df$MARKERS)
    } else{
      markers.meta <- rownames(genotypes.df)
    }
  }

  markers.meta  %<>%
    dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
    separate_markers(
      data = .,
      sep = "__",
      # markers.meta.all.only = TRUE,
      biallelic = TRUE,
      verbose = verbose)

  if (!rlang::has_name(markers.meta, "VARIANT_ID") &&
      nrow(markers.meta) == n.markers) {
    markers.meta %<>% dplyr::mutate(VARIANT_ID = variant.id)
  }

  if (is.null(source)) {
    if (rlang::has_name(markers.meta, "CALL_RATE") ||
        rlang::has_name(markers.meta, "REP_AVG")) {
      source <- "dart"
    } else {
      source <- "radiator"
    }
  }

  markers.meta %<>% dplyr::arrange(MARKERS)

  if (!is.null(strata)) {
    if ("dart" %in% source) {
      df.order <- colnames(genotypes.df)
      strata %<>%
        dplyr::mutate(
          TARGET_ID = factor(x = TARGET_ID, levels = df.order, ordered = TRUE)
        ) %>%
        dplyr::arrange(TARGET_ID) %>%
        dplyr::mutate(TARGET_ID = as.character(TARGET_ID))
      check <- strata %>%
        dplyr::mutate(
          CHECK = df.order,
          CHECK = dplyr::if_else(CHECK == TARGET_ID, TRUE, FALSE)
        ) %>%
        dplyr::select(CHECK) %>%
        purrr::flatten_lgl(.) %>%
        unique
    } else {
      df.order <- colnames(genotypes.df)
      strata %<>%
        dplyr::mutate(
          INDIVIDUALS = factor(x = INDIVIDUALS, levels = df.order, ordered = TRUE)
        ) %>%
        dplyr::arrange(INDIVIDUALS) %>%
        dplyr::mutate(INDIVIDUALS = as.character(INDIVIDUALS))
      check <- strata %>%
        dplyr::mutate(
          CHECK = df.order,
          CHECK = dplyr::if_else(CHECK == INDIVIDUALS, TRUE, FALSE)
        ) %>%
        dplyr::select(CHECK) %>%
        purrr::flatten_lgl(.) %>%
        unique
    }
    if (length(check) > 1 || !check) {
      rlang::abort("Problem with samples in strata and data...")
    }
    check <- df.order <- NULL
  }

  # Generate GDS format
  # 1. SNPRelate

  SNPRelate::snpgdsCreateGeno(
    gds.fn = radiator.temp.file,
    genmat = genotypes.df,
    sample.id = strata$INDIVIDUALS,
    snp.id = markers.meta$VARIANT_ID,
    snp.rs.id = markers.meta$LOCUS,
    snp.chromosome = markers.meta$CHROM,
    snp.position = markers.meta$POS,
    snp.allele = stringi::stri_join(markers.meta$REF, markers.meta$ALT, sep = "/"),
    snpfirstdim = TRUE,
    compress.annotation = "ZIP_RA",
    compress.geno = ""
  )

  # 2. SeqArray
  data.gds <- SeqArray::seqSNP2GDS(
    gds.fn = radiator.temp.file,
    out.fn = filename.gds,
    storage.option = "ZIP_RA",
    major.ref = TRUE, verbose = FALSE) %>%
    SeqArray::seqOpen(gds.fn = ., readonly = FALSE)
  file.remove(radiator.temp.file)

  # radiator skeleton folder
  radiator.gds <- radiator_gds_skeleton(data.gds)

  # source
  update_radiator_gds(data.gds, node.name = "source", value = source)

  # reference genome or de novo
  update_radiator_gds(
    data.gds,
    node.name = "reference.genome",
    value = detect_ref_genome(chromosome = markers.meta$CHROM, verbose = FALSE)
  )


  # bi- or multi-alllelic VCF
  update_radiator_gds(data.gds, node.name = "biallelic", value = biallelic)

  # Add STRATA to GDS
  update_radiator_gds(data.gds, node.name = "individuals", value = strata)

  # ADD MARKERS META to GDS
  update_radiator_gds(data.gds, node.name = "markers.meta", value = markers.meta)

  # Add genotypes metadata
  if (!is.null(genotypes.meta)) {
    if ("dart" %in% source && "counts" %in% source) {
      suppressWarnings(
        genotypes.meta %<>%
          dplyr::left_join(dplyr::select(strata, TARGET_ID, INDIVIDUALS), by = "TARGET_ID") %>%
          dplyr::select(-TARGET_ID) %>%
          dplyr::select(MARKERS, INDIVIDUALS, dplyr::everything(.))
      )
    }


    update_radiator_gds(
      data.gds,
      node.name = "genotypes.meta",
      value = genotypes.meta
    )
  }


  if (!is.null(dp)) {
    update_radiator_gds(data.gds, node.name = "DP", value = dp)
  }

  if (!is.null(ad)) {
    update_radiator_gds(data.gds, node.name = "AD", value = ad)
  }

  if (verbose) message("File written: ", filename.short.gds)

  if (open) {
    return (data.gds)
  } else {
    radiator::write_rad(data = data.gds)
    return(filename.gds)
  }
} # End rad_gds

# radiator_gds_skeleton---------------------------------------------------------
#' @title radiator_gds_skeleton
#' @description Generate a radiator.gds skeleton inside GDS
#' @name radiator_gds_skeleton
#' @rdname radiator_gds_skeleton
#' @keywords internal
#' @export
radiator_gds_skeleton <- function(gds) {
  radiator.gds <- gdsfmt::addfolder.gdsn(
    node = gds,
    name = "radiator",
    replace = TRUE)

  purrr::walk(
    .x = c(
      "source",
      "reference.genome",
      "biallelic",
      "id.clean",
      "individuals",
      "markers.meta",
      "genotypes.meta",
      "blacklist.individuals",
      "blacklist.markers"
    ),
    .f = update_radiator_gds,
    gds = gds,
    value = NULL,
    replace = FALSE,
    sync = FALSE,
    summary = FALSE,
    verbose = FALSE
  )
  return(radiator.gds)
}#End radiator_gds_skeleton

# update_radiator_gds-----------------------------------------------------------
#' @title update_radiator_gds
#' @description Update radiator gds node
#' @name update_radiator_gds
#' @rdname update_radiator_gds
#' @keywords internal
#' @export
update_radiator_gds <- function(
  gds,
  node.name,
  value,
  replace = TRUE,
  sync = FALSE,
  summary = FALSE,
  verbose = FALSE
) {

  if (node.name %in% c("AD", "DP")) {
    if (node.name == "AD") {
      format <- gdsfmt::index.gdsn(
        node = gds, path = "annotation/format", silent = TRUE)

      gdsfmt::add.gdsn(
        node = format,
        name = "AD",
        val = value,
        # storage = "vl_int",
        replace = TRUE,
        compress = "ZIP_RA",
        closezip = TRUE)
    }

    if (node.name == "DP") {
      format <- gdsfmt::index.gdsn(
        node = gds, path = "annotation/format", silent = TRUE)

      gdsfmt::add.gdsn(
        node = format,
        name = "DP",
        val = value,
        # storage = "vl_int",
        replace = TRUE,
        compress = "ZIP_RA",
        closezip = TRUE)
    }
  } else {
    radiator.gds <- gdsfmt::index.gdsn(
      node = gds, path = "radiator", silent = TRUE)

    gdsfmt::add.gdsn(
      node = radiator.gds,
      name = node.name,
      val = value,
      replace = replace,
      compress = "ZIP_RA",
      closezip = TRUE)
  }



  if (sync) {
    message("Synchronizing markers.meta")
    if (node.name == "markers.meta") {
      sync_gds(gds = gds, markers = as.integer(value$VARIANT_ID))
    }

    if (node.name == "individuals") {
      message("Synchronizing individuals")
      sync_gds(gds = gds, samples = value$INDIVIDUALS)
    }
  }

  if (summary) sum <- summary_gds(gds, verbose = verbose)
}#End update_radiator_gds

# extract_source----------------------------------------------------------------
#' @title extract_source
#' @description Extract the source of the radiator gds file
#' @name extract_source
#' @rdname extract_source
#' @keywords internal
#' @export
extract_source <- function(gds) {
  source <- gdsfmt::index.gdsn(node = gds, path = "radiator/source", silent = TRUE)
  if (!is.null(source)) source <- gdsfmt::read.gdsn(source)
  return(source)
}# End extract_source

# extract_markers_metadata------------------------------------------------------
#' @title extract_markers_metadata
#' @description Import gds or radiator markers meta node
#' @name extract_markers_metadata
#' @rdname extract_markers_metadata
#' @keywords internal
#' @export
extract_markers_metadata <- function(
  gds,
  markers.meta.select = NULL,
  radiator.node = TRUE,
  verbose = FALSE
) {

  # will switch radiator.mode to FALSE if returns null

  if (radiator.node) {
    markers.index <- gdsfmt::ls.gdsn(gdsfmt::index.gdsn(
      node = gds, path = "radiator/markers.meta", silent = TRUE))
    if (is.null(markers.index) || length(markers.index) == 0) radiator.node <- FALSE
  }

  if (!is.null(markers.meta.select)) {
    if (is.null(markers.index) || length(markers.index) == 0) {
      markers.index <- markers.meta.select
    } else {
      markers.index %<>% intersect(markers.meta.select)
    }
  }

  if (radiator.node) {
    markers.index %<>% magrittr::set_names(x = ., value = .)
    integrate_meta <- function(x, gds, markers.df) {
      markers.df[[x]] <- gdsfmt::read.gdsn(
        gdsfmt::index.gdsn(
          node = gds,
          path = stringi::stri_join("radiator/markers.meta/", x),
          silent = TRUE))
    }
  } else {
    if (!is.null(markers.meta.select)) {
      markers.meta.select %<>% intersect(
        c("VARIANT_ID", "CHROM", "LOCUS", "POS")) %>%
        stringi::stri_replace_all_fixed(
          str = .,
          pattern = c("VARIANT_ID", "CHROM", "LOCUS", "POS"),
          replacement =  c("variant.id", "chromosome", "annotation/id", "position"),
          vectorize_all = FALSE)
      markers.index <- markers.meta.select
    } else {
      markers.index <- c("variant.id", "chromosome", "annotation/id", "position")
    }

    markers.index %<>% magrittr::set_names(x = ., value = .)
    integrate_meta <- function(x, gds, markers.df) {
      markers.df[[x]] <- SeqArray::seqGetData(gds, x)
    }

    # markers.meta <- tibble::tibble(
    #   VARIANT_ID = SeqArray::seqGetData(gds, "variant.id"),
    #   CHROM = SeqArray::seqGetData(gds, "chromosome"),
    #   LOCUS = SeqArray::seqGetData(gds, "annotation/id"),
    #   POS = SeqArray::seqGetData(gds, "position")
    # )
  }
  markers.df <- list()
  markers.meta <- purrr::map_df(
    .x = markers.index,
    .f = integrate_meta,
    gds, markers.df
  )

  colnames(markers.meta)  %<>%
    stringi::stri_replace_all_fixed(
      str = .,
      pattern = c("variant.id", "chromosome", "annotation/id", "position"),
      replacement =  c("VARIANT_ID", "CHROM", "LOCUS", "POS"),
      vectorize_all = FALSE)

  if (rlang::has_name(markers.meta, "CHROM")) markers.meta$CHROM <- as.character(markers.meta$CHROM)
  if (rlang::has_name(markers.meta, "LOCUS")) markers.meta$LOCUS <- as.character(markers.meta$LOCUS)
  if (rlang::has_name(markers.meta, "POS")) markers.meta$POS <- as.character(markers.meta$POS)

  return(markers.meta)
}#End import_metadata

# extract_individuals-----------------------------------------------------------
#' @title extract_individuals
#' @description Import gds or radiator individuals node
#' @name extract_individuals
#' @rdname extract_individuals
#' @keywords internal
#' @export
extract_individuals <- function(
  gds,
  ind.field.select = NULL,
  radiator.node = TRUE,
  verbose = FALSE
) {

  # For SNPRelate data
  snprelate <- "SNPGDSFileClass" %in% class(gds)[1]

  if (snprelate) {
    individuals <- tibble::tibble(
      INDIVIDUALS = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "sample.id")))
  } else {
    # will switch radiator.mode to FALSE if returns null
    if (radiator.node) {
      id.index <- gdsfmt::ls.gdsn(gdsfmt::index.gdsn(
        node = gds, path = "radiator/individuals", silent = TRUE))
      if (is.null(id.index)) radiator.node <- FALSE
    }

    if (radiator.node) {
      if (!is.null(ind.field.select)) {
        id.index %<>% intersect(ind.field.select)
      }

      id.index %<>% magrittr::set_names(x = ., value = .)

      id.df <- list()
      id_fields <- function(x, gds, id.df) {
        id.df[[x]] <- gdsfmt::read.gdsn(
          gdsfmt::index.gdsn(
            node = gds,
            path = stringi::stri_join("radiator/individuals/", x),
            silent = TRUE))
      }

      individuals <- purrr::map_df(.x = id.index, .f = id_fields, gds, id.df)
    }

    if (!radiator.node) {
      individuals <- tibble::tibble(
        INDIVIDUALS = SeqArray::seqGetData(gds, "sample.id")
      )
    }
  }

  return(individuals)
}#End extract_individuals


# extract_genotypes_metadata----------------------------------------------------
#' @title extract_genotypes_metadata
#' @description Import gds or radiator genotypes meta node
#' @name extract_genotypes_metadata
#' @rdname extract_genotypes_metadata
#' @keywords internal
#' @export
extract_genotypes_metadata <- function(
  gds,
  genotypes.meta.select = NULL,
  radiator.node = TRUE,
  index.only = FALSE,
  verbose = FALSE
) {
  # gds = input
  if (!radiator.node) return(NULL)

  geno.index <- gdsfmt::ls.gdsn(gdsfmt::index.gdsn(
    node = gds, path = "radiator/genotypes.meta", silent = TRUE))

  if (!index.only) {
    if (!is.null(genotypes.meta.select)) {
      geno.index %<>% intersect(genotypes.meta.select)
    }

    geno.index %<>% magrittr::set_names(x = ., value = .)

    genotypes.df <- list()
    integrate_meta <- function(x, gds, genotypes.df) {
      genotypes.df[[x]] <- gdsfmt::read.gdsn(
        gdsfmt::index.gdsn(
          node = gds,
          path = stringi::stri_join("radiator/genotypes.meta/", x),
          silent = TRUE))
    }

    genotypes.meta <- purrr::map_df(.x = geno.index, .f = integrate_meta, gds, genotypes.df)
    sum <- summary_gds(gds, verbose = FALSE)
    n.obs <- sum$n.ind * sum$n.markers
    if (nrow(genotypes.meta) != n.obs) {
      markers <- extract_markers_metadata(gds = gds, markers.meta.select = "MARKERS") %$% MARKERS
      individuals <- extract_individuals(gds = gds, ind.field.select = "INDIVIDUALS") %$% INDIVIDUALS
      genotypes.meta %<>% dplyr::filter(
        MARKERS %in% markers,
        INDIVIDUALS %in% individuals
      )
    }

    return(genotypes.meta)
  } else {
    return(geno.index)
  }
}#End extract_genotypes_metadata

# extract_coverage--------------------------------------------------------------
#' @title extract_coverage
#' @description Extract coverage information from a GDS file
#' @rdname extract_coverage
#' @keywords internal
#' @export
extract_coverage <- function(gds, markers = TRUE, ind = TRUE, depth = NULL) {
  coverage.info <- list()
  source <- extract_source(gds)
  if ("dart" %in% source) {
    # 1 and 2 rows
    if ("counts" %in% source) {
      dart.cov <- extract_genotypes_metadata(gds) %>%
        tidy2wide(x = ., gds = gds, wide = FALSE) %$%
        data.tidy
      want <- c("READ_DEPTH", "ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH")
      m <- dplyr::group_by(dart.cov, MARKERS) %>%
        dplyr::summarise_at(.tbl = ., .vars = want, .funs = mean, na.rm = TRUE) %>%
        dplyr::mutate_at(.tbl = ., .vars = want, .funs = round, digits = 0) %>%
        dplyr::mutate_at(.tbl = ., .vars = want, .funs = as.integer) %>%
        dplyr::ungroup(.)
      coverage.info$markers.tot <- dplyr::group_by(dart.cov, MARKERS) %>%
        dplyr::summarise(READ_DEPTH = sum(READ_DEPTH, na.rm = TRUE)) %$%
        READ_DEPTH
      coverage.info$markers.mean <- m$READ_DEPTH
      coverage.info$ref.mean <- m$ALLELE_REF_DEPTH
      coverage.info$alt.mean <- m$ALLELE_ALT_DEPTH
      m <- NULL
      if (ind) {
        coverage.info$ind.cov.tot <- dplyr::group_by(dart.cov, INDIVIDUALS) %>%
          dplyr::summarise(READ_DEPTH = sum(READ_DEPTH, na.rm = TRUE)) %$%
          READ_DEPTH
        coverage.info$ind.cov.mean <- dplyr::group_by(dart.cov, INDIVIDUALS) %>%
          dplyr::summarise(READ_DEPTH = as.integer(round(mean(READ_DEPTH, na.rm = TRUE), 0))) %$%
          READ_DEPTH
      }
      dart.cov <- NULL
    } else {
      dart.cov <- extract_markers_metadata(
        gds,
        markers.meta.select = c("AVG_COUNT_REF", "AVG_COUNT_SNP"))
      if (!is.null(dart.cov)) {
        coverage.info$markers.mean <- as.integer(
          round(dart.cov$AVG_COUNT_REF + dart.cov$AVG_COUNT_SNP, 0))
        coverage.info$ref.mean <- as.integer(round(dart.cov$AVG_COUNT_REF, 0))
        coverage.info$alt.mean <- as.integer(round(dart.cov$AVG_COUNT_SNP, 0))
        dart.cov <- NULL
      }
    }

  } else {
    depth.info <- gdsfmt::index.gdsn(gds, "annotation/format/DP", silent = TRUE)
    if (!is.null(depth.info)) {
      if (is.null(depth)) {
        depth <- SeqArray::seqGetData(gds, "annotation/format/DP")
      }
      if (ind) {
        coverage.info$ind.cov.tot <- as.integer(round(rowSums(x = depth$data, na.rm = TRUE, dims = 1L), 0))
        coverage.info$ind.cov.mean <- as.integer(round(rowMeans(x = depth$data, na.rm = TRUE, dims = 1L), 0))
      }

      if (markers) {
        coverage.info$markers.mean <- as.integer(round(colMeans(x = depth$data, na.rm = TRUE, dims = 1L), 0))
        coverage.info$markers.tot <- as.integer(round(colSums(x = depth$data, na.rm = TRUE, dims = 1L), 0))
      }
    }
  }
  return(coverage.info)
}#End extract_coverage

# sync GDS----------------------------------------------------------------------
#' @title sync_gds
#' @description Synchronize gds with samples and markers. If left NULL, the info
#' is first search in the radiator node, if not found, it goes in the level above.
#' An argument also allows to reset the filters.
#' @rdname sync_gds
#' @keywords internal
#' @export
sync_gds <- function(gds, samples = NULL, markers = NULL, reset = FALSE, verbose = FALSE) {

  if (reset) {
    SeqArray::seqResetFilter(
      object = gds,
      sample = TRUE,
      variant = TRUE,
      verbose = verbose
    )
  } else {
    if (is.null(markers)) {
      if (verbose) message("synchronizing GDS with current markers")
      markers <- extract_markers_metadata(
        gds = gds, markers.meta.select = "VARIANT_ID", verbose = verbose) %$%
        VARIANT_ID
    } else {
      if (verbose) message("synchronizing GDS with provided markers")
      markers <- as.integer(markers)
    }

    if (is.null(samples)) {
      if (verbose) message("synchronizing GDS with current samples")
      samples <- extract_individuals(
        gds = gds, ind.field.select = "INDIVIDUALS", verbose = verbose) %$%
        INDIVIDUALS
    } else {
      if (verbose) message("synchronizing GDS with provided samples")
    }

    SeqArray::seqSetFilter(
      object = gds, sample.id = samples, variant.id = markers, verbose = verbose)
  }
}#End sync_gds


# summarize GDS-----------------------------------------------------------------

#' @title summary_gds
#' @description Summary of gds object or file: number of samples and markers
#' @rdname summary_gds
#' @keywords internal
#' @export
summary_gds <- function(gds, verbose = TRUE) {
  data.type <- radiator::detect_genomic_format(gds)
  if (data.type == "gds.file") {
    gds <- radiator::read_rad(gds, verbose = FALSE)
    data.type <- "SeqVarGDSClass"
  }
  if (verbose) message("Data summary: ")
  check <- SeqArray::seqGetFilter(gds)
  n.ind <- length(check$sample.sel[check$sample.sel])
  n.markers <- length(check$variant.sel[check$variant.sel])
  if (verbose) message("    number of samples: ", n.ind)
  if (verbose) message("    number of markers: ", n.markers)
  return(res = list(n.ind = n.ind, n.markers = n.markers))
}# End summary_gds

# update blacklist of markers---------------------------------------------------
#' @title update_bl_markers
#' @description Generate and update blacklist of markers on radiator GDS node
#' @rdname update_bl_markers
#' @keywords internal
#' @export
update_bl_markers <- function(
  gds,
  generate = FALSE,
  update = NULL,
  extract = FALSE,
  bl.gds = NULL
) {

  # Check conditions...
  # if(!generate && !is.null(update) && !extract && is.null(bl.gds)) generate <- TRUE
  if (is.null(update) && !extract && is.null(bl.gds)) generate <- TRUE
  if (!generate && is.null(bl.gds)) extract <- TRUE

  radiator.gds <- gdsfmt::index.gdsn(node = gds, path = "radiator", silent = TRUE)
  if (is.null(radiator.gds)) radiator.gds <- gds

  #Run
  if (!generate) {
    test <- gdsfmt::get.attr.gdsn(gdsfmt::index.gdsn(
      node = radiator.gds, path = "blacklist.markers", silent = TRUE))$R.class[1]

    if (is.null(test)) generate <- TRUE
  }

  if (generate) {
    bl.gds <- tibble::tibble(MARKERS = character(0), FILTER = character(0))
    gdsfmt::add.gdsn(
      node = radiator.gds,
      name = "blacklist.markers",
      val = bl.gds,
      replace = TRUE,
      compress = "ZIP_RA",
      closezip = TRUE)
  }

  if (extract) {
    bl.gds <- tibble::tibble(
      MARKERS = gdsfmt::read.gdsn(gdsfmt::index.gdsn(
        node = radiator.gds, path = "blacklist.markers/MARKERS", silent = TRUE)),
      FILTER = gdsfmt::read.gdsn(gdsfmt::index.gdsn(
        node = radiator.gds, path = "blacklist.markers/FILTER", silent = TRUE))
    )
  }


  if (!is.null(update)) {
    bl.gds %<>% dplyr::bind_rows(update)
    suppressWarnings(gdsfmt::add.gdsn(
      node = radiator.gds,
      name = "blacklist.markers",
      val = bl.gds,
      replace = TRUE,
      compress = "ZIP_RA",
      closezip = TRUE))
  }
  return(bl.gds)
}#End update_bl_markers

# update blacklist of individuals-----------------------------------------------
#' @title update_bl_individuals
#' @description Generate and update blacklist of individuals on radiator GDS node
#' @rdname update_bl_individuals
#' @keywords internal
#' @export
update_bl_individuals <- function(
  gds,
  generate = FALSE,
  update = NULL,
  extract = FALSE,
  bl.i.gds = NULL
) {

  # Check conditions...
  # if(!generate && !is.null(update) && !extract && is.null(bl.gds)) generate <- TRUE
  if (is.null(update) && !extract && is.null(bl.i.gds)) generate <- TRUE
  if (!generate && is.null(bl.i.gds)) extract <- TRUE

  radiator.gds <- gdsfmt::index.gdsn(node = gds, path = "radiator", silent = TRUE)
  if (is.null(radiator.gds)) radiator.gds <- radiator_gds_skeleton(gds)

  #Run
  if (!generate) {
    test <- gdsfmt::get.attr.gdsn(gdsfmt::index.gdsn(
      node = radiator.gds, path = "blacklist.individuals", silent = TRUE))$R.class[1]

    if (is.null(test)) generate <- TRUE
  }

  if (generate) {
    bl.i.gds <- tibble::tibble(INDIVIDUALS = character(0), FILTER = character(0))
    gdsfmt::add.gdsn(
      node = radiator.gds,
      name = "blacklist.individuals",
      val = bl.i.gds,
      replace = TRUE,
      compress = "ZIP_RA",
      closezip = TRUE)
  }

  if (extract) {
    bl.i.gds <- tibble::tibble(
      INDIVIDUALS = gdsfmt::read.gdsn(gdsfmt::index.gdsn(
        node = radiator.gds, path = "blacklist.individuals/INDIVIDUALS", silent = TRUE)),
      FILTER = gdsfmt::read.gdsn(gdsfmt::index.gdsn(
        node = radiator.gds, path = "blacklist.individuals/FILTER", silent = TRUE))
    )
  }


  if (!is.null(update)) {
    bl.i.gds %<>% dplyr::bind_rows(update)
    gdsfmt::add.gdsn(
      node = radiator.gds,
      name = "blacklist.individuals",
      val = bl.i.gds,
      replace = TRUE,
      compress = "ZIP_RA",
      closezip = TRUE)
  }
  return(bl.i.gds)
}#End update_bl_individuals


# generate id/sample statistics-------------------------------------------------
#' @title generate_id_stats
#' @description Generate id/sample statistics
#' @rdname generate_id_stats
#' @keywords internal
#' @export
generate_id_stats <- function (
  gds,
  missing = TRUE,
  heterozygosity = TRUE,
  coverage = TRUE,
  subsample = NULL,
  depth.info = FALSE,
  path.folder = NULL,
  plot = TRUE,
  file.date = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE
) {

  res <- list() # return result in this list
  if (is.null(file.date)) file.date <- format(Sys.time(), "%Y%m%d@%H%M")

  if (!is.null(subsample)) {
    SeqArray::seqSetFilter(object = gds,
                           variant.id = subsample,
                           # sample.id = id.info$INDIVIDUALS,
                           action = "push+set",
                           verbose = FALSE)
  }

  radiator.gds <- gdsfmt::index.gdsn(node = gds, path = "radiator", silent = TRUE)
  if (is.null(radiator.gds)) radiator.gds <- radiator_gds_skeleton(gds)

  id.info.class <- gdsfmt::get.attr.gdsn(gdsfmt::index.gdsn(
    node = radiator.gds, path = "individuals", silent = TRUE))$R.class[1]

  if (is.null(id.info.class)) {
    id.info <- tibble::tibble(
      INDIVIDUALS = SeqArray::seqGetData(gdsfile = gds, var.name = "sample.id"))
  } else {
    id.info <- tibble::tibble(INDIVIDUALS = gdsfmt::read.gdsn(
      gdsfmt::index.gdsn(
        node = radiator.gds, path = "individuals/INDIVIDUALS", silent = TRUE)))
    strata <- gdsfmt::index.gdsn(
      node = radiator.gds, path = "individuals/STRATA", silent = TRUE)
    if (!is.null(strata)) {
      id.info %<>% dplyr::mutate(STRATA = gdsfmt::read.gdsn(strata))
    }
  }

  # missing
  if (missing) {
    # info
    id.info %<>%
      dplyr::mutate(
        MISSING_PROP = round(
          SeqArray::seqMissing(
            gdsfile = gds, per.variant = FALSE,
            .progress = TRUE, parallel = parallel.core
          ) , 6)
      )
    # stats
    id.stats.m <- tibble_stats(
      x = id.info$MISSING_PROP,
      group = "missing genotypes")
  } else {
    id.stats.m <- NULL
  }


  # heterozygosity
  if (heterozygosity) {
    # info
    id.info %<>%
      dplyr::mutate(
        HETEROZYGOSITY = round(SeqVarTools::heterozygosity(
          gdsobj = gds, margin = "by.sample", use.names = FALSE
        ), 6)
      )

    # stats
    id.stats.h <- tibble_stats(
      x = id.info$HETEROZYGOSITY,
      group = "heterozygosity")
  } else {
    id.stats.h <- NULL
  }

  # coverage
  if (coverage && !depth.info) {
    check_dp <- function(gds) {
      SeqArray::seqGetData(
        gdsfile = gds,
        var.name = "annotation/format/DP")
    }#check_dp
    check_dp_safe <- purrr::safely(.f = check_dp)
    depth <- check_dp_safe(gds = gds) %$% result
    if (is.null(depth)) {
      if (verbose) message("No coverage info: skipping for individual's stats")
      depth.info <- FALSE
      depth <- NULL
      coverage <- FALSE
    } else {
      depth.info <- TRUE
    }
    source <- extract_source(gds)
    if ("dart" %in% source && "counts" %in% source) {
      depth.info <- TRUE
      coverage <- TRUE
      depth <- NULL
    }
  } else {
    depth <- NULL
  }

  if (coverage) {
    # info
    dp <- extract_coverage(gds, markers = FALSE, depth = depth)

    id.info %<>%
      dplyr::mutate(
        COVERAGE_TOTAL = dp$ind.cov.tot,
        COVERAGE_MEAN = dp$ind.cov.mean
      )

    # stats
    id.stats.c <- dplyr::bind_rows(
      tibble_stats(
        x = as.numeric(id.info$COVERAGE_TOTAL),
        group = "total coverage"),
      tibble_stats(
        x = as.numeric(id.info$COVERAGE_MEAN),
        group = "mean coverage")
    )
  } else {
    id.stats.c <- NULL
  }

  # Generate the stats
  id.stats <- dplyr::bind_rows(id.stats.m, id.stats.h, id.stats.c)
  id.levels <- c("total coverage", "mean coverage", "missing genotypes", "heterozygosity")
  id.stats$GROUP <- factor(x = id.stats$GROUP, levels = id.levels, ordered = TRUE)
  id.stats$GROUP <- droplevels(x = id.stats$GROUP)

  id.stats.filename <- stringi::stri_join("individuals.qc.stats_", file.date, ".tsv")
  readr::write_tsv(x = id.info, path = file.path(path.folder, id.stats.filename))
  id.stats.filename <- stringi::stri_join("individuals.qc.stats.summary_", file.date, ".tsv")
  readr::write_tsv(x = id.stats, path = file.path(path.folder, id.stats.filename))
  if (verbose) message("File written: individuals qc info and stats summary")

  # Generate plots
  if (plot) {
    # correlations info

    # Note to myself: here we need to add conditions when some stats are not requested...
    corr.info <- stringi::stri_join("Correlations:\n")

    if (coverage && missing) {
      cm <- floor(stats::cor(id.info$COVERAGE_TOTAL, id.info$MISSING_PROP, use = "pairwise.complete.obs") * 100) / 100
      cmt <- stringi::stri_join("    total coverage & missing = ", cm)
      corr.info <- stringi::stri_join(corr.info, cmt)
    }
    if (coverage) {
      cc <- ceiling(stats::cor(id.info$COVERAGE_TOTAL, id.info$COVERAGE_MEAN,use = "pairwise.complete.obs") * 100) / 100
      cct <- stringi::stri_join("\n    total coverage & mean coverage = ", cc)
      corr.info <- stringi::stri_join(corr.info, cct)
    }
    if (coverage && heterozygosity) {
      ch <- ceiling(stats::cor(id.info$COVERAGE_TOTAL, id.info$HETEROZYGOSITY, use = "pairwise.complete.obs") * 100) / 100
      cht <- stringi::stri_join("\n    total coverage & heterozygosity = ", ch)
      corr.info <- stringi::stri_join(corr.info, cht)
    }
    if (missing && heterozygosity) {
      mh <- floor(stats::cor(id.info$HETEROZYGOSITY, id.info$MISSING_PROP, use = "pairwise.complete.obs") * 100) / 100
      mht <- stringi::stri_join("\n    missing & heterozygosity = ", mh)
      corr.info <- stringi::stri_join(corr.info, mht)
    }

    # heterozygosity.info
    if (heterozygosity) {
      n.markers <- summary_gds(gds, verbose = FALSE) %$% n.markers
      n.markers.range <- ceiling((id.stats[[2,6]] - id.stats[[2,2]] ) * n.markers)
      n.markers.iqr <- ceiling(id.stats[[2,7]] * n.markers)

      corr.info <- stringi::stri_join(
        "n. het markers in the bp range = ", n.markers.range,
        "\nn. het markers in the bp IQR = ", n.markers.iqr,
        "\n\n", corr.info
      )
    }


    if (!is.null(subsample)) {
      subtitle.ind.stats <- stringi::stri_join(
        "Markers subsampled: ", length(subsample), "\n\n", corr.info)
    } else {
      subtitle.ind.stats <- corr.info
    }

    res$fig <- boxplot_stats(
      data = id.stats,
      title = "Individual's QC stats",
      subtitle = subtitle.ind.stats,
      x.axis.title = NULL,
      y.axis.title = "Statistics",
      facet.columns = TRUE,
      bp.filename = stringi::stri_join("individuals.qc_", file.date, ".pdf"),
      path.folder = path.folder)
    if (verbose) message("File written: individuals qc plot")

  }

  if (!is.null(subsample)) {
    SeqArray::seqSetFilter(gds, action = "pop", verbose = FALSE)
  }
  res$stats <- id.stats
  res$info <- id.info
  return(res)
}#End generate_id_stats

# generate markers statistics---------------------------------------------------
#' @title generate_markers_stats
#' @description Generate markers statistics
#' @rdname generate_markers_stats
#' @keywords internal
#' @export
generate_markers_stats <- function (
  gds,
  coverage = TRUE,
  allele.coverage = TRUE,
  missing = TRUE,
  mac = TRUE,
  heterozygosity = TRUE,
  snp.position.read = TRUE,
  snp.per.locus = TRUE,
  path.folder = NULL,
  filename = NULL,
  fig.filename = NULL,
  plot = TRUE,
  force.stats = FALSE,
  subsample = NULL,
  file.date = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE
) {

  # coverage <- TRUE
  # allele.coverage <- TRUE
  # mac <- TRUE
  # missing <- TRUE
  # heterozygosity <- TRUE
  # snp.per.locus <- TRUE
  # snp.position.read <- TRUE
  # fig.filename = NULL
  # plot = TRUE

  if (is.null(file.date)) file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  source <- extract_source(gds)

  # Filenames
  # Stats
  if (is.null(filename)) {
    markers.stats.file <- stringi::stri_join("markers_metadata_stats_", file.date, ".tsv")
  } else {
    markers.stats.file <- stringi::stri_join(filename, "_stats_", file.date, ".tsv")
  }

  # tables
  if (is.null(filename)) {
    filename <- stringi::stri_join("markers_metadata_", file.date, ".tsv")
  } else {
    filename <- stringi::stri_join(filename, "_", file.date, ".tsv")
  }

  # Figure
  if (is.null(fig.filename)) {
    fig.filename <- stringi::stri_join("markers_qc_", file.date, ".pdf")
  } else {
    fig.filename <- stringi::stri_join(fig.filename, "_", file.date, ".pdf")
  }

  info <- extract_markers_metadata(gds = gds)

  if (!is.null(subsample)) {
    SeqArray::seqSetFilter(object = gds,
                           variant.id = subsample,
                           action = "push+set",
                           verbose = FALSE)
    info %<>% dplyr::filter(VARIANT_ID %in% subsample)
  }
  n.markers <- length(info$VARIANT_ID)

  if (allele.coverage) mac <- TRUE # by default the mac info is needed...
  if (!rlang::has_name(info, "COL")) snp.position.read <- FALSE


  # Coverage
  if (coverage) {
    if (!rlang::has_name(info, "COVERAGE_TOTAL") || force.stats) {
      if ("dart" %in% source) {
        if ("counts" %in% source) {
          mc <- extract_coverage(gds, ind = FALSE)#coverage
        } else {
          if (is.null(subsample)) {
            mc <- extract_coverage(gds, ind = FALSE)#coverage
            # mc$markers.tot <- NULL
          } else {
            mc <- list()
            # mc$markers.tot <- NULL
            mc$markers.mean <- info$COVERAGE_MEAN
            mc$ref.mean <- info$REF_DEPTH_MEAN
            mc$alt.mean <- info$ALT_DEPTH_MEAN
          }
        }
      } else {
        mc <- extract_coverage(gds, ind = FALSE)#coverage
      }

      if (!is.null(mc$markers.tot)) {
        info$COVERAGE_TOTAL <- mc$markers.tot
        stats.ct <- tibble_stats(
          x = info$COVERAGE_TOTAL,
          group = "total coverage")
      } else {
        stats.ct <- NULL
      }

      if (!is.null(mc$markers.mean)) {
        info$COVERAGE_MEAN <- mc$markers.mean
        stats.cm <- tibble_stats(
          x = info$COVERAGE_MEAN,
          group = "mean coverage")
      } else {
        stats.cm <- NULL
      }

      if (!is.null(mc$ref.mean)) {
        info$REF_DEPTH_MEAN <- mc$ref.mean
        allele.coverage <- mac <- TRUE
      }
      if (!is.null(mc$alt.mean)){
        info$ALT_DEPTH_MEAN <- mc$alt.mean
      }
      mc <- NULL
    } else {
      stats.ct <- tibble_stats(
        x = info$COVERAGE_TOTAL,
        group = "total coverage")
      stats.cm <- tibble_stats(
        x = info$COVERAGE_MEAN,
        group = "mean coverage")
    }
  } else {
    stats.cm <- stats.ct <- NULL
  }

  if (mac) {
    if (!rlang::has_name(info, "MAC_GLOBAL") || force.stats) {
      info %<>%
        dplyr::bind_cols(
          SeqArray::seqAlleleCount(
            gdsfile = gds,
            ref.allele = NULL,
            .progress = TRUE,
            parallel = parallel.core) %>%
            unlist(.) %>%
            matrix(
              data = .,
              nrow = n.markers, ncol = 2, byrow = TRUE,
              dimnames = list(rownames = info$VARIANT_ID,
                              colnames = c("REF_COUNT", "ALT_COUNT"))) %>%
            tibble::as_tibble(.)
        ) %>%
        dplyr::mutate(
          MAC_GLOBAL = dplyr::if_else(ALT_COUNT < REF_COUNT, ALT_COUNT, REF_COUNT),
          MAF_GLOBAL = MAC_GLOBAL / (REF_COUNT + ALT_COUNT),
          ALT_CHECK = dplyr::if_else(MAC_GLOBAL == ALT_COUNT, "ALT", "REF"),
          ALT_COUNT = NULL
        )
    }

    stats.mac <- tibble_stats(
      x = info$MAC_GLOBAL,
      group = "MAC")

  } else {
    stats.mac <- NULL
  }


  if (allele.coverage) {

    if (!rlang::has_name(info, "REF_DEPTH_TOTAL") || force.stats) {

      ad.info <- gdsfmt::index.gdsn(gds, "annotation/format/AD", silent = TRUE)
      if (!is.null(ad.info)) {
        ad <- SeqArray::seqGetData(gds, "annotation/format/AD") %$% data

        info %<>%
          dplyr::bind_cols(
            #REF and ALT Total read depth
            colSums(x = ad, na.rm = TRUE, dims = 1L) %>%
              unlist(.) %>%
              matrix(
                data = .,
                nrow = n.markers, ncol = 2, byrow = TRUE,
                dimnames = list(rownames = info$VARIANT_ID,
                                colnames = c("RT", "AT"))) %>%
              tibble::as_tibble(.),
            #REF and ALT Mean read depth
            colMeans(x = ad, na.rm = TRUE, dims = 1L) %>%
              unlist(.) %>%
              matrix(
                data = .,
                nrow = n.markers, ncol = 2, byrow = TRUE,
                dimnames = list(rownames = info$VARIANT_ID,
                                colnames = c("RM", "AM"))) %>%
              tibble::as_tibble(.)
          ) %>%
          dplyr::mutate(
            REF_DEPTH_TOTAL = dplyr::if_else(ALT_CHECK == "ALT", RT, AT),
            ALT_DEPTH_TOTAL = dplyr::if_else(ALT_CHECK == "ALT", AT, RT),
            RT = NULL, AT = NULL,
            REF_DEPTH_MEAN = dplyr::if_else(ALT_CHECK == "ALT", RM, AM),
            ALT_DEPTH_MEAN = dplyr::if_else(ALT_CHECK == "ALT", AM, RM),
            ALT_CHECK = NULL, RM = NULL, AM = NULL
          )
        ad <- NULL
      }
      ad.info <- NULL
    }

    if (rlang::has_name(info, "REF_DEPTH_TOTAL")) {
      stats.rt <- tibble_stats(
        x = info$REF_DEPTH_TOTAL,
        group = "total ref depth")
      stats.at <- tibble_stats(
        x = info$ALT_DEPTH_TOTAL,
        group = "total alt depth")
    } else {
      stats.rt <- stats.at <-  NULL
    }

    if (rlang::has_name(info, "REF_DEPTH_MEAN")) {
      stats.rm <- tibble_stats(
        x = info$REF_DEPTH_MEAN,
        group = "mean ref depth")
      stats.am <- tibble_stats(
        x = info$ALT_DEPTH_MEAN,
        group = "mean alt depth")
    } else {
      stats.rm <- stats.am <-  NULL
    }


  } else {
    if (mac) info %<>% dplyr::select(-ALT_CHECK)
    stats.rt <- stats.at <- stats.am <- stats.rm <- NULL

  }

  if (missing) {
    if (!rlang::has_name(info, "MISSING_PROP") || force.stats) {
      info %<>%
        dplyr::mutate(
          MISSING_PROP = SeqArray::seqMissing(
            gdsfile = gds,
            per.variant = TRUE, .progress = TRUE, parallel = parallel.core)
        )
    }
    stats.m <- tibble_stats(
      x = info$MISSING_PROP,
      group = "missing genotypes")

  } else {
    stats.m <- NULL
  }

  if (heterozygosity) {
    if (!rlang::has_name(info, "HET_OBS") || force.stats) {
      info %<>%
        dplyr::mutate(
          HET_OBS = round(SeqVarTools::heterozygosity(
            gdsobj = gds, margin = "by.variant", use.names = FALSE
          ), 6),
          FIS = SeqVarTools::inbreedCoeff(gdsobj = gds, margin = "by.variant")
        )
    }
    stats.h <- tibble_stats(
      x = info$HET_OBS,
      group = "observed heterozygosity")

    stats.f <- tibble_stats(
      x = info$FIS,
      group = "inbreeding coefficient (Fis)")

  } else {
    stats.h <- stats.f <- NULL
  }

  if (snp.per.locus) {
    if (!rlang::has_name(info, "SNP_PER_LOCUS") || force.stats) {
      biallelic <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(
        node = gds,
        path = "radiator/biallelic", silent = TRUE))
      if (!is.logical(biallelic)) biallelic <- FALSE
      if (biallelic) {
        info %<>%
          dplyr::group_by(LOCUS) %>%
          dplyr::mutate(SNP_PER_LOCUS = n()) %>%
          dplyr::ungroup(.)
      } else {
        if (verbose) message("With haplotype vcf, number of SNP/locus is counted based on the REF allele")
        if (verbose) message("    this stat is good only if nuc length is equal between REF and ALT haplotypes")
        if (verbose) message("    stacks haplotype vcf: ok")
        if (verbose) message("    dDocent/freebayes haplotype vcf: be careful")

        info %<>%
          dplyr::mutate(
            SNP_PER_LOCUS = stringi::stri_length(
              str = SeqArray::seqGetData(gdsfile = gds, var.name = "$ref"))
          )
      }
    }
    stats.s <- tibble_stats(
      x = info$SNP_PER_LOCUS,
      group = "SNPs per locus")
  } else {
    stats.s <- NULL
  }

  if (snp.position.read) {
    stats.sp <- tibble_stats(
      x = dplyr::distinct(info, COL) %$% COL,
      group = "SNP position on read")
  } else {
    stats.sp <- NULL
  }

  # Generate the stats
  stats <- dplyr::bind_rows(
    stats.cm, stats.ct,#coverage
    stats.mac, #MAC
    stats.rt, stats.at, stats.rm, stats.am, # allele depth
    stats.m, # missing
    stats.h, stats.f,#het obs and fis
    stats.sp, # snp position on read
    stats.s#snps per locus
  )
  stats.cm <- stats.ct <- stats.mac <- stats.rt <- stats.at <- stats.rm <- stats.am <- NULL
  stats.m <- stats.h <- stats.f <- stats.s <- stats.sp <- NULL

  markers.levels <- c(
    "total coverage", "total ref depth", "total alt depth",
    "mean coverage", "mean ref depth", "mean alt depth",
    "missing genotypes",
    "MAC",  "observed heterozygosity",
    "inbreeding coefficient (Fis)",
    "SNP position on read",
    "SNPs per locus"
  )
  stats$GROUP <- factor(x = stats$GROUP, levels = markers.levels,
                        ordered = TRUE)
  stats$GROUP <- droplevels(x = stats$GROUP)

  info %<>% dplyr::arrange(MARKERS) %>%
    readr::write_tsv(x = ., path = file.path(path.folder, filename))

  # Generate plots
  if (plot) {
    # correlations info

    # Note to myself: here we need to add conditions when some stats are not requested...
    corr.info <- stringi::stri_join("Correlations:\n")

    if (rlang::has_name(info, "COVERAGE_TOTAL") && missing) {
      cm <- floor(stats::cor(info$COVERAGE_TOTAL, info$MISSING_PROP, use = "pairwise.complete.obs") * 100) / 100
      cmt <- stringi::stri_join("    total coverage & missing = ", cm)
      corr.info <- stringi::stri_join(corr.info, cmt)
    }
    if (rlang::has_name(info, "COVERAGE_TOTAL") && rlang::has_name(info, "COVERAGE_MEAN")) {
      cc <- ceiling(stats::cor(info$COVERAGE_TOTAL, info$COVERAGE_MEAN, use = "pairwise.complete.obs") * 100) / 100
      cct <- stringi::stri_join("\n    total coverage & mean coverage = ", cc)
      corr.info <- stringi::stri_join(corr.info, cct)
    }
    if (rlang::has_name(info, "COVERAGE_TOTAL") && heterozygosity) {
      ch <- ceiling(stats::cor(info$COVERAGE_TOTAL, info$HET_OBS, use = "pairwise.complete.obs") * 100) / 100
      cht <- stringi::stri_join("\n    total coverage & heterozygosity = ", ch)
      corr.info <- stringi::stri_join(corr.info, cht)
    }
    if (missing && heterozygosity) {
      mh <- floor(stats::cor(info$HET_OBS, info$MISSING_PROP, use = "pairwise.complete.obs") * 100) / 100
      mht <- stringi::stri_join("\n    missing & heterozygosity = ", mh)
      corr.info <- stringi::stri_join(corr.info, mht)
    }

    # corr.info <- stringi::stri_join(
    #   "Correlations:\n",
    #   "    total coverage & missing = ", cm,
    #   "\n    total coverage & mean coverage = ", cc,
    #   "\n    total coverage & heterozygosity = ", ch,
    #   "\n    missing & heterozygosity = ", mh
    # )


    cm <- cc <- ch <- mh <- mht <- cmt <- cct <- cht <- NULL

    if (!is.null(subsample)) {
      subtitle.stats <- stringi::stri_join(
        "Markers subsampled: ", length(subsample), "\n\n", corr.info)
    } else {
      subtitle.stats <- corr.info
    }

    if (!heterozygosity && !rlang::has_name(info, "COVERAGE_TOTAL")) {
      subtitle.stats <- NULL
    }


    fig <- boxplot_stats(
      data = stats,
      title = "Marker's QC stats",
      subtitle = subtitle.stats,
      x.axis.title = NULL,
      y.axis.title = "Statistics",
      facet.columns = TRUE,
      bp.filename = fig.filename,
      path.folder = path.folder)


    dplyr::mutate_if(
      .tbl = stats, .predicate = is.numeric, .funs = format, scientific = FALSE) %>%
      readr::write_tsv(x = ., path = file.path(path.folder, markers.stats.file))
  } else {
    fig <- NULL
  }

  # restore original filters (subsampling...)
  if (!is.null(subsample)) {
    SeqArray::seqSetFilter(gds, action = "pop", verbose = TRUE)
  } else {
    radiator.gds <- gdsfmt::index.gdsn(node = gds, path = "radiator", silent = TRUE)
    # Update GDS
    if (!is.null(radiator.gds)) {
      gdsfmt::add.gdsn(
        node = radiator.gds,
        name = "markers.meta",
        val = info,
        replace = TRUE,
        compress = "ZIP_RA",
        closezip = TRUE)
    }
  }

  return(list(info = info, stats = stats, fig = fig))
}#End generate_markers_stats

# missingness per markers per pop-----------------------------------------------
#' @title missing_per_pop
#' @description Generate missingness per markers per pop helper table
#' @rdname missing_per_pop
#' @keywords internal
#' @export
missing_per_pop <- function(
  gds,
  strata,
  parallel.core = parallel::detectCores() - 1
) {
  missing_pop <- function(
    id.select, gds, markers.df,
    parallel.core = parallel::detectCores() - 1
  ) {
    SeqArray::seqSetFilter(object = gds,
                           sample.id = id.select$INDIVIDUALS,
                           action = "set",
                           verbose = FALSE)

    res <- markers.df %>%
      dplyr::mutate(
        MISSING_PROP = SeqArray::seqMissing(
          gdsfile = gds,
          per.variant = TRUE, .progress = TRUE, parallel = parallel.core - 1))

    mis_many_markers <- function(threshold, x) {
      nrow(dplyr::filter(x, MISSING_PROP <= threshold))
    }#End how_many_markers

    helper.table <- tibble::tibble(MISSING_PROP = seq(0.1, 1, 0.1)) %>%
      dplyr::mutate(
        WHITELISTED_MARKERS = purrr::map_int(
          .x = seq(0.1, 1, 0.1), .f = mis_many_markers, x = res)
      )
    SeqArray::seqResetFilter(
      object = gds, sample = TRUE, variant = FALSE, verbose = FALSE)
    return(helper.table)
  }#End missing_pop
  markers.df <- extract_markers_metadata(gds = gds, markers.meta.select = "MARKERS")
  n.markers <- nrow(markers.df)

  res <- strata %>%
    dplyr::group_by(STRATA) %>%
    tidyr::nest(data = ., .key = id.select) %>%
    dplyr::mutate(MISSING_POP = purrr::map(
      .x = .$id.select,
      .f = missing_pop,
      gds = gds,
      markers.df = markers.df,
      parallel.core = parallel.core)) %>%
    tidyr::unnest(data = ., MISSING_POP) %>%
    dplyr::mutate(BLACKLISTED_MARKERS = n.markers - WHITELISTED_MARKERS)
  return(res)
}#End missing_per_pop


# gds2tidy ---------------------------------------------------------------------
#' @title gds2tidy
#' @description GDS to tidy...
#' @rdname gds2tidy
#' @keywords internal
#' @export
gds2tidy <- function(
  gds,
  markers.meta = NULL,
  individuals = NULL,
  parallel.core = parallel::detectCores() - 1,
  ...
) {
  if (is.null(individuals)) {
    individuals <- extract_individuals(gds = gds, verbose = FALSE)
  }

  if (is.null(markers.meta)) {
    markers.meta <- extract_markers_metadata(gds = gds, verbose = TRUE)
  }
  want <- intersect(c("MARKERS", "CHROM", "LOCUS", "POS", "COL", "REF", "ALT"),
                    names(markers.meta))

  # summary_gds(gds)
  tidy.data <- suppressWarnings(
    SeqArray::seqGetData(
      gdsfile = gds, var.name = "$dosage_alt") %>%
    # magrittr::set_colnames(x = ., value = markers.meta$VARIANT_ID) %>%
    magrittr::set_colnames(x = ., value = markers.meta$MARKERS) %>%
    magrittr::set_rownames(x = ., value = individuals$INDIVIDUALS) %>%
      data.table::as.data.table(x = ., keep.rownames = "INDIVIDUALS") %>%
      data.table::melt.data.table(
        data = .,
        id.vars = "INDIVIDUALS",
        variable.name = "MARKERS",
        value.name = "GT_BIN",
        variable.factor = FALSE) %>%
      tibble::as_tibble(.) %>%
      dplyr::left_join(dplyr::select(markers.meta, dplyr::one_of(want)), by = "MARKERS") %>%
      dplyr::mutate(
        MARKERS = factor(x = MARKERS,
                         levels = markers.meta$MARKERS, ordered = TRUE),
        INDIVIDUALS = factor(x = INDIVIDUALS,
                             levels = individuals$INDIVIDUALS,
                             ordered = TRUE)) %>%
      dplyr::arrange(MARKERS, INDIVIDUALS)
  )


  # re-calibration of ref/alt alleles ------------------------------------------
  # if (verbose) message("\nCalculating REF/ALT alleles...")
  tidy.data <- radiator::calibrate_alleles(
    data = tidy.data,
    # biallelic = TRUE,
    parallel.core = parallel.core,
    verbose = FALSE,
    gt = FALSE, gt.vcf = FALSE
  ) %$% input

  # include strata
  colnames(individuals) <- stringi::stri_replace_all_fixed(
    str = colnames(individuals),
    pattern = "STRATA",
    replacement = "POP_ID",
    vectorize_all = FALSE)
  if (rlang::has_name(individuals, "POP_ID")) {
    suppressWarnings(
      tidy.data %<>%
        dplyr::left_join(dplyr::select(individuals, INDIVIDUALS, POP_ID), by = "INDIVIDUALS")
    )
  } else {
    tidy.data %<>% dplyr::mutate(POP_ID = 1L)
  }

  tidy.data %<>% dplyr::arrange(POP_ID, INDIVIDUALS)
  return(tidy.data)
} #End tidy gds
