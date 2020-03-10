# radiator gds constructor------------------------------------------------------
#' @title radiator gds constructor
#' @description helper function to construct a radiator gds
#' @name radiator_gds
#' @rdname radiator_gds
#' @keywords internal
#' @export
radiator_gds <- function(
  genotypes = NULL,
  biallelic = TRUE,
  data.source = NULL,
  geno.coding = c("alt.dos", "ref.dos"),
  strata = NULL,
  markers.meta = NULL,
  genotypes.meta = NULL,
  dp = NULL, # read depth
  ad = NULL, # allele depth
  gl = NULL, # genotype likelihood
  filename = NULL,
  open = FALSE,
  verbose = TRUE
) {
  # genotypes <- data
  message("Generating GDS...")
  if (!biallelic) rlang::abort("Biallelic data required")

  # Types of genotypes coding:
  # SNPRelate is coded by ref allele dosage...
  geno.coding <- match.arg(
    arg = geno.coding,
    choices = c("alt.dos", "ref.dos"),
    several.ok = FALSE
  )

  # Filename -------------------------------------------------------------------
  if (is.null(filename)) {
    filename <- generate_filename(
      name.shortcut = filename,
      extension = "gds.rad"
    ) %$% filename
  }

  if (length(filename) > 1) filename <- filename$filename

  temp.file <- stringi::stri_join("radiator_temp_", format(Sys.time(), "%Y%m%d@%H%M"))

  #Empty vcf -------------------------------------------------------------------
  write_vcf(filename = temp.file, source = NULL, empty = TRUE)

  # Empty GDS ------------------------------------------------------------------
  data.gds <- suppressMessages(
    suppressWarnings(
      SeqArray::seqVCF2GDS(
        vcf.fn = paste0(temp.file, ".vcf"),
        out.fn = filename,
        parallel = 1L,
        storage.option = "ZIP_RA",
        verbose = FALSE
      ) %>%
        SeqArray::seqOpen(gds.fn = ., readonly = FALSE)
    )
  )

  # Removing VCF file
  file.remove(paste0(temp.file, ".vcf"))

  # Add genotype node now...
  genotype.node <- gdsfmt::addfolder.gdsn(
    node = data.gds,
    visible = TRUE,
    name = "genotype",
    replace = TRUE
  )

  # radiator skeleton folder ---------------------------------------------------
  radiator.gds <- radiator_gds_skeleton(data.gds)

  # data.source ----------------------------------------------------------------
  if (is.null(data.source)) data.source <- "radiator"
  update_radiator_gds(data.gds, node.name = "data.source", value = data.source)

  # bi- or multi-alllelic VCF --------------------------------------------------
  update_radiator_gds(data.gds, node.name = "biallelic", value = biallelic)

  # Coverage -------------------------------------------------------------------
  update_radiator_gds(data.gds, node.name = "DP", value = dp)
  update_radiator_gds(data.gds, node.name = "AD", value = ad)
  update_radiator_gds(data.gds, node.name = "GL", value = gl)

  # STRATA ---------------------------------------------------------------------
  if (!is.null(strata)) {
    if (!rlang::has_name(strata, "POP_ID") && !rlang::has_name(strata, "STRATA")) {
      strata %<>% dplyr::mutate(STRATA = "pop")
    }
    # rare case #3
    if (!rlang::has_name(strata, "INDIVIDUALS")) {
      strata %<>% dplyr::mutate(INDIVIDUALS = colnames(genotypes))
    }
    strata %<>% dplyr::distinct(INDIVIDUALS, .keep_all = TRUE)
  } else {
    if (rlang::has_name(genotypes, "STRATA")) {
      strata <- dplyr::distinct(.data = genotypes, INDIVIDUALS, STRATA)
    } else {
      strata <- dplyr::distinct(.data = genotypes, INDIVIDUALS) %>%
        dplyr::mutate(STRATA = 1L)
    }
  }
  n.ind <- nrow(strata)

  # MARKERS META ---------------------------------------------------------------
  # rare case (usually DArT) where some IDs are still in the metadata...
  if (!is.null(strata)) {
    suppressWarnings(
      markers.meta %<>% dplyr::select(-dplyr::one_of(strata$INDIVIDUALS))
    )
  }

  # check <- markers.meta$MARKERS
  # generate all the markers metadata
  markers.meta  %<>%
    dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
    separate_markers(
      data = .,
      sep = "__",
      markers.meta.all.only = TRUE,
      biallelic = TRUE,
      verbose = verbose)
  n.markers <- nrow(markers.meta)

  # reference genome or de novo
  update_radiator_gds(
    data.gds,
    node.name = "reference.genome",
    value = detect_ref_genome(chromosome = markers.meta$CHROM, verbose = FALSE)
  )

  # ADD MARKERS META to GDS
  update_radiator_gds(data.gds, node.name = "markers.meta", value = markers.meta)
  update_radiator_gds(
    data.gds,
    node.name = "variant.id",
    value = markers.meta$VARIANT_ID,
    radiator.gds = FALSE,
    replace = TRUE)
  update_radiator_gds(
    data.gds,
    node.name = "chromosome",
    value = markers.meta$CHROM,
    radiator.gds = FALSE,
    replace = TRUE)
  update_radiator_gds(
    data.gds,
    node.name = "position",
    value = markers.meta$POS,
    radiator.gds = FALSE,
    replace = TRUE)
  update_radiator_gds(
    data.gds,
    node.name = "allele",
    value = stringi::stri_join(markers.meta$REF, markers.meta$ALT, sep = ","),
    radiator.gds = FALSE,
    replace = TRUE)

  annotation <- gdsfmt::index.gdsn(
    node = data.gds, path = "annotation", silent = TRUE)
  update_radiator_gds(
    annotation,
    node.name = "id",
    value = markers.meta$LOCUS,
    radiator.gds = FALSE,
    replace = TRUE)

  if ("dart" %in% data.source) {
    annotation <- gdsfmt::index.gdsn(
      node = data.gds, path = "annotation/info", silent = TRUE)
    update_radiator_gds(
      annotation,
      node.name = "NS",
      value = NULL,
      radiator.gds = FALSE,
      replace = TRUE)
  }

  # variant.id <- markers.meta$VARIANT_ID
  n.snp <- nrow(markers.meta)
  markers.meta <- NULL

  # Genotypes  -----------------------------------------------------------------
  # Add individuals.meta
  update_radiator_gds(gds = data.gds, node.name = "individuals.meta", value = strata)
  update_radiator_gds(
    gds = data.gds,
    node.name = "sample.id",
    value = strata$INDIVIDUALS,
    radiator.gds = FALSE,
    replace = TRUE)

  suppressWarnings(
    gdsfmt::add.gdsn(
      node = genotype.node,
      name = "data",
      val = genotypes,
      valdim = c(2L, n.ind, n.snp),
      replace = TRUE,
      compress = "ZIP_RA",
      storage = "bit2",
      closezip = TRUE)
  )
  # gdsfmt::readmode.gdsn(node = genotype.node)
  # gdsfmt::compression.gdsn(node = genotype.node, compress="ZIP_RA")
  # data.node <- gdsfmt::index.gdsn(node = genotype.node, path = "data", silent = TRUE)
  # gdsfmt::write.gdsn(node = genotype.node, val = genotypes)


  gdsfmt::put.attr.gdsn(genotype.node, "VariableName", "GT")
  gdsfmt::put.attr.gdsn(genotype.node, "Description", "Genotype")

  gdsfmt::add.gdsn(
    node = genotype.node,
    name = "@data",
    val = rep(1L, n.snp),
    replace = TRUE,
    compress = "ZIP_RA",
    storage = "uint8",
    visible = FALSE,
    closezip = TRUE)

  genotypes <- NULL

  # Add genotypes metadata
  if (!is.null(genotypes.meta)) {
    update_radiator_gds(
      data.gds,
      node.name = "genotypes.meta",
      value = genotypes.meta
    )
  }

  genotypes.meta <- NULL
  message("File written: ", folder_short(filename))
  data.gds <- radiator::write_rad(data = data.gds)
  if (open) data.gds <- radiator::read_rad(data = data.gds)
  return (data.gds)
} # End rad_gds


# gt2array----------------------------------------------------------------------
#' @title gt2array
#' @description Alternate allele dosage gt (GT_BIN) to
#' presence/absence array for the genotypes in GDS
#' @rdname gt2array
#' @keywords internal
#' @export
gt2array <- function(gt.bin, n.ind, n.snp) {
  # The genotypes coding and the array...
  genotypes <- cbind(
    dplyr::case_when(
      gt.bin == 0L ~ 0L,
      gt.bin == 1L ~ 0L,
      gt.bin == 2L ~ 1L
    ),
    dplyr::if_else(gt.bin == 2L, 1L, gt.bin)
  )

  #NA...
  genotypes[is.na(genotypes)] <- 0x0F

  # dimensions
  dim(genotypes) <- c(n.snp, n.ind, 2)

  # permute the array: alleles, samples, markers
  genotypes <- aperm(a = genotypes, c(3,2,1))

  # dimension names
  dimnames(genotypes) <- list(allele=NULL, sample=NULL, variant=NULL)

  return(genotypes)
}# End gt2array

# tidy2gds ---------------------------------------------------------------------
#' @title tidy2gds
#' @description Generate a GDS file/object from a tidy dataset requires (GT_BIN)
#' @rdname tidy2gds
#' @keywords internal
#' @export
tidy2gds <- function(x) {

  # markers.meta will always be sorted by variant id
  if (!rlang::has_name(x, "VARIANT_ID")) {
    mk.col <- dplyr::intersect(colnames(x), c("CHROM", "LOCUS", "POS"))
    x %<>%
      dplyr::mutate_at(.tbl = ., .vars = mk.col, .funs = as.character) %>%
      dplyr::arrange(MARKERS) %>%
      dplyr::mutate(VARIANT_ID = as.integer(factor(MARKERS)))
  }
  want <- c("VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS", "COL", "REF",
            "ALT")
  markers.meta <- suppressWarnings(
    x %>%
      dplyr::select(dplyr::one_of(want)) %>%
      dplyr::distinct(.) %>%
      dplyr::arrange(VARIANT_ID)
  )
  notwanted <- c("MARKERS", "CHROM", "LOCUS", "POS", "COL", "REF","ALT")
  x <- suppressWarnings(x %>% dplyr::select(-dplyr::one_of(notwanted)))
  if (rlang::has_name(x, "POP_ID")) x %<>% dplyr::rename(STRATA = POP_ID)

  # the genotypes array is sorted by individuals and then by variant...
  x %<>% dplyr::arrange(STRATA, INDIVIDUALS, VARIANT_ID)

  x <- radiator_gds(
    genotypes = gt2array(
      gt.bin = x$GT_BIN,
      n.ind = length(unique(x$INDIVIDUALS)),
      n.snp = length(unique(x$VARIANT_ID))
    ),
    strata = generate_strata(data = x),
    markers.meta = markers.meta,
    open = TRUE
  )
  return(x)
} #End tidy2gds

# gds2tidy ---------------------------------------------------------------------
#' @title gds2tidy
#' @description GDS to tidy...
#' @rdname gds2tidy
# @param pop.id (logical) When \code{pop.id = TRUE}, the strata returns
# the stratification colname \code{POP_ID}.
# Default: \code{pop.id = FALSE}, returns \code{STRATA}.
#' @keywords internal
#' @export
gds2tidy <- function(
  gds,
  markers.meta = NULL,
  markers.meta.select = NULL,
  wide = FALSE,
  individuals = NULL,
  pop.id = TRUE,
  calibrate.alleles = TRUE,
  parallel.core = parallel::detectCores() - 1,
  ...
) {
  if (is.null(individuals)) {
    individuals <- extract_individuals_metadata(gds = gds, whitelist = TRUE, verbose = FALSE)
  }

  if (wide) markers.meta.select <- c("MARKERS")
  if (is.null(markers.meta)) {
    if (is.null(markers.meta.select)) {
      markers.meta <- extract_markers_metadata(
        gds = gds,
        whitelist = TRUE,
        verbose = TRUE
      )
    } else {
      markers.meta <- extract_markers_metadata(
        gds = gds,
        markers.meta.select = markers.meta.select,
        whitelist = TRUE,
        verbose = TRUE
      )
    }
  }
  want <- intersect(
    c("MARKERS", "CHROM", "LOCUS", "POS", "COL", "REF", "ALT", "COL",
      "CALL_RATE", "AVG_COUNT_REF", "AVG_COUNT_SNP", "REP_AVG",
      "ONE_RATIO_REF", "ONE_RATIO_SNP"),
    names(markers.meta))

  # summary_gds(gds)
  tidy.data <-
    SeqArray::seqGetData(
      gdsfile = gds, var.name = "$dosage_alt") %>%
    magrittr::set_colnames(x = ., value = markers.meta$MARKERS) %>%
    magrittr::set_rownames(x = ., value = individuals$INDIVIDUALS) %>%
    tibble::as_tibble(x = ., rownames = "INDIVIDUALS")

  if (!wide) {
    tidy.data <- suppressWarnings(
      data.table::as.data.table(x = tidy.data) %>%
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
  }


  # should make this optional --------------------------------------------------
  if (pop.id) {
    # include strata
    colnames(individuals) <- stringi::stri_replace_all_fixed(
      str = colnames(individuals),
      pattern = "STRATA",
      replacement = "POP_ID",
      vectorize_all = FALSE)
    if (rlang::has_name(individuals, "POP_ID")) {
      suppressWarnings(
        tidy.data %<>%
          dplyr::left_join(
            dplyr::select(individuals, INDIVIDUALS, POP_ID)
            , by = "INDIVIDUALS"
          )
      )
    }
    tidy.data %<>% dplyr::arrange(POP_ID, INDIVIDUALS)
  } else {
    # include strata
    colnames(individuals) <- stringi::stri_replace_all_fixed(
      str = colnames(individuals),
      pattern = "POP_ID",
      replacement = "STRATA",
      vectorize_all = FALSE)
    if (rlang::has_name(individuals, "STRATA")) {
      suppressWarnings(
        tidy.data %<>%
          dplyr::left_join(dplyr::select(individuals, INDIVIDUALS, STRATA), by = "INDIVIDUALS")
      )
    }
    if (!rlang::has_name(tidy.data, "STRATA") && !rlang::has_name(tidy.data, "POP_ID")) {
      tidy.data %<>% dplyr::mutate(POP_ID = 1L)
    }
    # re-calibration of ref/alt alleles ------------------------------------------
    # if (verbose) message("\nCalculating REF/ALT alleles...")
    if (calibrate.alleles && !wide) {
      tidy.data <- radiator::calibrate_alleles(
        data = tidy.data,
        # biallelic = TRUE,
        parallel.core = parallel.core,
        verbose = FALSE,
        gt = FALSE, gt.vcf = FALSE
      ) %$% input
    }
    tidy.data %<>% dplyr::arrange(STRATA, INDIVIDUALS)
  }
  return(tidy.data)
} #End tidy gds

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
    replace = TRUE
    )

  purrr::walk(
    .x = c(
      "data.source",
      "reference.genome",
      "biallelic",
      "id.clean",
      "individuals.meta",
      "markers.meta",
      "genotypes.meta"
    ),
    .f = update_radiator_gds,
    gds = gds,
    radiator.gds = TRUE,
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
  radiator.gds = TRUE,
  node.name,
  value,
  replace = TRUE,
  sync = FALSE,
  summary = FALSE,
  verbose = FALSE
) {

  if (node.name %in% c("AD", "DP", "GL")) {
    format <- gdsfmt::index.gdsn(
      node = gds, path = "annotation/format", silent = TRUE)

    if (node.name == "AD" && !is.null(value)) {
      gdsfmt::add.gdsn(
        node = format,
        name = "AD",
        val = value,
        replace = TRUE,
        compress = "ZIP_RA",
        closezip = TRUE)
    }

    if (node.name == "DP" && !is.null(value)) {
      gdsfmt::add.gdsn(
        node = format,
        name = "DP",
        val = value,
        replace = TRUE,
        compress = "ZIP_RA",
        closezip = TRUE)
    }

    if (node.name == "GL" && !is.null(value)) {
      gdsfmt::add.gdsn(
        node = format,
        name = "GL",
        val = value,
        replace = TRUE,
        compress = "ZIP_RA",
        closezip = TRUE)
    }

  } else {

    if (radiator.gds) {
      radiator.gds <- gdsfmt::index.gdsn(
        node = gds, path = "radiator", silent = TRUE)

      suppressWarnings(
        gdsfmt::add.gdsn(
          node = radiator.gds,
          name = node.name,
          val = value,
          replace = replace,
          compress = "ZIP_RA",
          closezip = TRUE)
      )
    } else {
      suppressWarnings(
        gdsfmt::add.gdsn(
          node = gds,
          name = node.name,
          val = value,
          replace = replace,
          compress = "ZIP_RA",
          closezip = TRUE)
      )
    }
  }


  if (sync) {
    if (verbose) message("Synchronizing markers.meta")
    if (node.name == "markers.meta") {
      if (rlang::has_name(value, "FILTERS")) {
        sync_gds(gds = gds, variant.id = as.integer(value$VARIANT_ID[value$FILTERS == "whitelist"]))
      } else {
        sync_gds(gds = gds, variant.id = as.integer(value$VARIANT_ID))
      }
    }

    if (node.name == "individuals.meta") {
      if (verbose) message("Synchronizing individuals")
      if (rlang::has_name(value, "FILTERS")) {
        sync_gds(gds = gds, samples = value$INDIVIDUALS[value$FILTERS == "whitelist"])
      } else {
        sync_gds(gds = gds, samples = value$INDIVIDUALS)
      }
    }
  }

  if (summary) sum <- summary_gds(gds, verbose = verbose)
}#End update_radiator_gds

# extract_data_source----------------------------------------------------------------
#' @title extract_data_source
#' @description Extract the data.source of the radiator gds file
#' @name extract_data_source
#' @rdname extract_data_source
#' @keywords internal
#' @export
extract_data_source <- function(gds) {
  data.source <- gdsfmt::index.gdsn(node = gds, path = "radiator/data.source", silent = TRUE)
  if (!is.null(data.source)) data.source <- gdsfmt::read.gdsn(data.source)
  return(data.source)
}# End extract_data_source

# extract_individuals_metadata-----------------------------------------------------------
#' @title extract_individuals_metadata
#' @description Import gds or radiator individuals node
#' @name extract_individuals_metadata
#' @rdname extract_individuals_metadata
#' @param gds The gds object.
#' @param ind.field.select (optional, character) Default:\code{ind.field.select = NULL}.
#' @param radiator.node (optional, logical) Default:\code{radiator.node = TRUE}.
#' @param whitelist (optional, logical) Default:\code{whitelist = FALSE}.
#' @param blacklist (optional, logical) Default:\code{blacklist = FALSE}.
#' @inheritParams radiator_common_arguments
# @keywords internal
#' @export
extract_individuals_metadata <- function(
  gds,
  ind.field.select = NULL,
  radiator.node = TRUE,
  whitelist = FALSE,
  blacklist = FALSE,
  verbose = FALSE
) {

  # For SNPRelate data
  snprelate <- "SNPGDSFileClass" %in% class(gds)[1]
  keep.one <- FALSE
  if (whitelist) blacklist <- FALSE
  if (blacklist) whitelist <- FALSE
  if (snprelate) {
    individuals <- tibble::tibble(
      INDIVIDUALS = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "sample.id")))
  } else {
    # will switch radiator.mode to FALSE if returns null
    if (radiator.node) {
      id.index <- gdsfmt::ls.gdsn(gdsfmt::index.gdsn(
        node = gds, path = "radiator/individuals.meta", silent = TRUE))
      if (is.null(id.index)) radiator.node <- FALSE
    }

    if (radiator.node) {
      if (!is.null(ind.field.select)) {
        if (length(ind.field.select) == 1) keep.one <- TRUE

        if (whitelist || blacklist) {
          if (!"FILTERS" %in% ind.field.select) ind.field.select  %<>% c("FILTERS")
        }
        id.index %<>% intersect(ind.field.select)
      }

      id.index %<>% magrittr::set_names(x = ., value = .)

      id.df <- list()
      id_fields <- function(x, gds, id.df) {
        id.df[[x]] <- gdsfmt::read.gdsn(
          gdsfmt::index.gdsn(
            node = gds,
            path = stringi::stri_join("radiator/individuals.meta/", x),
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

  if (!whitelist && !blacklist && !rlang::has_name(individuals, "FILTERS")) {
    individuals %<>% dplyr::mutate(FILTERS = "whitelist")
  }

  if (whitelist && rlang::has_name(individuals, "FILTERS")){
    individuals %<>%
      dplyr::filter(FILTERS == "whitelist") %>%
      dplyr::select(-FILTERS)
  }
  if (blacklist && rlang::has_name(individuals, "FILTERS")){
    individuals %<>%
      dplyr::filter(FILTERS != "whitelist") %>%
      dplyr::select(-FILTERS)
  }

  return(individuals)
}#End extract_individuals_metadata


# extract_markers_metadata------------------------------------------------------
#' @title extract_markers_metadata
#' @description Import gds or radiator markers meta node
#' @name extract_markers_metadata
#' @rdname extract_markers_metadata
#' @param gds The gds object.
#' @param markers.meta.select (optional, character) Default:\code{markers.meta.select = NULL}.
#' @param radiator.node (optional, logical) Default:\code{radiator.node = TRUE}.
#' @param whitelist (optional, logical) Default:\code{whitelist = FALSE}.
#' @param blacklist (optional, logical) Default:\code{blacklist = FALSE}.
#' @inheritParams radiator_common_arguments

# @keywords internal
#' @export
extract_markers_metadata <- function(
  gds,
  markers.meta.select = NULL,
  radiator.node = TRUE,
  whitelist = FALSE,
  blacklist = FALSE,
  verbose = FALSE
) {

  # will switch radiator.mode to FALSE if returns null
  keep.one <- FALSE

  if (whitelist) blacklist <- FALSE
  if (blacklist) whitelist <- FALSE

  if (radiator.node) {
    markers.index <- gdsfmt::ls.gdsn(gdsfmt::index.gdsn(
      node = gds, path = "radiator/markers.meta", silent = TRUE))
    if (is.null(markers.index) || length(markers.index) == 0) radiator.node <- FALSE
  }

  if (!is.null(markers.meta.select)) {

    if (length(markers.meta.select) == 1) keep.one <- TRUE

    if (whitelist || blacklist) {
      if (!"FILTERS" %in% markers.meta.select) markers.meta.select %<>% c("FILTERS")
    }
    if (is.null(markers.index) || length(markers.index) == 0) {
      markers.index <- markers.meta.select
    } else {
      markers.index %<>% intersect(markers.meta.select)
    }
  }

  if (radiator.node) {
    markers.index %<>% magrittr::set_names(x = ., value = .)
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
  }


  integrate_meta <- function(x, gds, markers.df, radiator.node) {
    if (radiator.node) {
      markers.df[[x]] <- gdsfmt::read.gdsn(
        gdsfmt::index.gdsn(
          node = gds,
          path = stringi::stri_join("radiator/markers.meta/", x),
          silent = TRUE))
    } else {
      markers.df[[x]] <- SeqArray::seqGetData(gds, x)
    }
  }

  markers.df <- list()
  markers.meta <- purrr::map_df(
    .x = markers.index,
    .f = integrate_meta,
    gds,
    markers.df,
    radiator.node
  )

  colnames(markers.meta) %<>%
    stringi::stri_replace_all_fixed(
      str = .,
      pattern = c("variant.id", "chromosome", "annotation/id", "position"),
      replacement =  c("VARIANT_ID", "CHROM", "LOCUS", "POS"),
      vectorize_all = FALSE
    )

  mk.col <- dplyr::intersect(colnames(markers.meta), c("CHROM", "LOCUS", "POS"))
  markers.meta  %<>%
    dplyr::mutate_at(.tbl = ., .vars = mk.col, .funs = as.character)

  if (!whitelist && !blacklist && !rlang::has_name(markers.meta, "FILTERS")) {
    markers.meta %<>% dplyr::mutate(FILTERS = "whitelist")
  }
  if (whitelist && rlang::has_name(markers.meta, "FILTERS")){
    markers.meta %<>%
      dplyr::filter(FILTERS == "whitelist") %>%
      dplyr::select(-FILTERS)
  }
  if (blacklist && rlang::has_name(markers.meta, "FILTERS")){
    markers.meta %<>%
      dplyr::filter(FILTERS != "whitelist")# %>%
    # dplyr::select(-FILTERS)
  }
  return(markers.meta)
}#End import_metadata


# extract_genotypes_metadata----------------------------------------------------
#' @title extract_genotypes_metadata
#' @description Import gds or radiator genotypes meta node
#' @name extract_genotypes_metadata
#' @rdname extract_genotypes_metadata
#' @param gds The gds object.
#' @param genotypes.meta.select (optional, character) Default:\code{genotypes.meta.select = NULL}.
#' @param genotypes (optional, character) Default: \code{genotypes = FALSE}.
#' @param radiator.node (optional, logical) Default: \code{radiator.node = TRUE}.
#' @param index.only (optional, logical) Default: \code{index.only = FALSE}.
#' @param sync.markers.individuals (optional, logical) Default: \code{sync.markers.individuals = TRUE}.
#' @param whitelist (optional, logical) Default: \code{whitelist = FALSE}.
#' @param blacklist (optional, logical) Default: \code{blacklist = FALSE}.
#' @inheritParams radiator_common_arguments
# @keywords internal
#' @export
extract_genotypes_metadata <- function(
  gds,
  genotypes.meta.select = NULL,
  genotypes = FALSE,
  radiator.node = TRUE,
  index.only = FALSE,
  sync.markers.individuals = TRUE,
  whitelist = FALSE,
  blacklist = FALSE,
  verbose = FALSE
) {
  # gds = input
  keep.one <- FALSE
  if (!radiator.node) return(NULL)
  if (whitelist) blacklist <- FALSE
  if (blacklist) whitelist <- FALSE

  geno.index <- gdsfmt::ls.gdsn(gdsfmt::index.gdsn(
    node = gds, path = "radiator/genotypes.meta", silent = TRUE))

  if (!index.only) {
    if (!is.null(genotypes.meta.select)) {
      if (length(genotypes.meta.select) == 1) keep.one <- TRUE
      if (whitelist || blacklist) {
        if (!"FILTERS" %in% genotypes.meta.select) genotypes.meta.select %<>% c("FILTERS")
      }
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
    if (sync.markers.individuals) {
      sum <- summary_gds(gds, verbose = FALSE)
      n.obs <- sum$n.ind * sum$n.markers
      if (nrow(genotypes.meta) != n.obs) {
        markers <- extract_markers_metadata(
          gds = gds,
          markers.meta.select = "MARKERS",
          whitelist = TRUE
        ) %$% MARKERS
        individuals <- extract_individuals_metadata(
          gds = gds,
          ind.field.select = "INDIVIDUALS",
          whitelist = TRUE
        ) %$% INDIVIDUALS
        genotypes.meta %<>%
          dplyr::filter(MARKERS %in% markers, INDIVIDUALS %in% individuals)
      }
    }
    if (!whitelist && !blacklist && !rlang::has_name(genotypes.meta, "FILTERS")) {
      genotypes.meta %<>% dplyr::mutate(FILTERS = "whitelist")
    }
    if (whitelist && rlang::has_name(genotypes.meta, "FILTERS")){
      genotypes.meta %<>%
        dplyr::filter(FILTERS == "whitelist") %>%
        dplyr::select(-FILTERS)
    }
    if (blacklist && rlang::has_name(genotypes.meta, "FILTERS")){
      genotypes.meta %<>%
        dplyr::filter(FILTERS != "whitelist") %>%
        dplyr::select(-FILTERS)
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
#' @param gds The gds object.
#' @param markers (optional, logical) Default:\code{markers = TRUE}.
#' @param ind (optional, logical) Default:\code{ind = TRUE}.
#' @param update.gds (optional, logical) Default:\code{update.gds = FALSE}.
#' @param depth.tibble (optional, logical) Returns the depth info in a tibble instead of list
#' with total and mean coverage info.
#' Used internally.
#' Default:\code{depth.tibble = FALSE}.
#' @inheritParams radiator_common_arguments

# @keywords internal
#' @export
extract_coverage <- function(
  gds,
  markers = TRUE,
  ind = TRUE,
  update.gds = FALSE,
  depth.tibble = FALSE,
  parallel.core = parallel::detectCores() - 2,
  verbose = FALSE
) {
  #Test
  # markers = TRUE
  # ind = TRUE
  # update.gds = FALSE
  # parallel.core = parallel::detectCores() - 2
  # verbose = FALSE


  coverage.info <- list()
  data.source <- extract_data_source(gds)



  n.markers<- summary_gds(gds, verbose = FALSE)$n.markers

  # DArT counts and VCFs -------------------------------------------------------
  if (!"dart" %in% data.source || "counts" %in% data.source) {
    depth <- extract_genotypes_metadata(
      gds,
      genotypes.meta.select = c("MARKERS", "INDIVIDUALS", "READ_DEPTH",
                                "ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH"),
      sync.markers.individuals = TRUE,
      whitelist = TRUE
    )
    if (length(depth) == 0 || is.null(depth)) {

      if (stringi::stri_detect_fixed(str = data.source, pattern = "Stacks") &&
          !gdsfmt::read.gdsn(gdsfmt::index.gdsn(
            node = gds, path = "radiator/biallelic", silent = TRUE))) {
        have <- integer(0)
      } else {
        # detect FORMAT fields available
        have <-  SeqArray::seqSummary(
          gdsfile = gds,
          varname = "annotation/format",
          check = "none", verbose = FALSE)$ID
      }

      if (length(have) > 0) {
        want <- c("DP", "AD", "CATG")
        parse.format.list <- purrr::keep(.x = have, .p = have %in% want)
        if (verbose) message("Extracting ", paste0(parse.format.list, collapse = ", "), " information...")
        # work on parallelization of this part
        depth <- tidy2wide(
          x = purrr::map(
            .x = parse.format.list,
            .f = parse_gds_metadata,
            gds = gds,
            verbose = FALSE,
            parallel.core = parallel.core
          ) %>%
            purrr::flatten(.) %>%
            purrr::flatten_df(.),
          gds = gds,
          tidy = TRUE,
          wide = FALSE) %$% data.tidy
        parse.format.list <- want <- NULL
      } else {
        depth <- NULL
      }
      have <- NULL
    }

    if (!is.null(depth)) {
      # detect.catg
      catg.depth <- rlang::has_name(depth, "C_DEPTH")
      if (catg.depth) {
        depth %<>%
          dplyr::mutate(MARKERS = as.character(MARKERS)) %>%
          dplyr::left_join(
            extract_markers_metadata(
              gds = gds,
              markers.meta.select = c("MARKERS", "REF", "ALT"),
              whitelist = TRUE
            ), by = "MARKERS") %>%
          dplyr::mutate(
            ALLELE_REF_DEPTH = dplyr::case_when(
              REF == "C" ~ C_DEPTH,
              REF == "A" ~ A_DEPTH,
              REF == "T" ~ T_DEPTH,
              REF == "G" ~ G_DEPTH
            ),
            ALLELE_ALT_DEPTH = dplyr::case_when(
              ALT == "C" ~ C_DEPTH,
              ALT == "A" ~ A_DEPTH,
              ALT == "T" ~ T_DEPTH,
              ALT == "G" ~ G_DEPTH
            )
          )
      }#catg.depth

      if (depth.tibble) return(depth)

      want <- c("READ_DEPTH", "ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH")
      have <- colnames(depth)
      want <- purrr::keep(.x = have, .p = have %in% want)
      colnames_rep <- function(x, total = FALSE, mean = FALSE) {
        if (total) {
          x <- stringi::stri_replace_all_fixed(
            str = x,
            pattern = c("READ_DEPTH", "ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH"),
            replacement = c("COVERAGE_TOTAL", "REF_DEPTH_TOTAL", "ALT_DEPTH_TOTAL"),
            vectorize_all = FALSE
          )
        }
        if (mean) {
          x <- stringi::stri_replace_all_fixed(
            str = x,
            pattern = c("READ_DEPTH", "ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH"),
            replacement = c("COVERAGE_MEAN", "REF_DEPTH_MEAN", "ALT_DEPTH_MEAN"),
            vectorize_all = FALSE
          )
        }
        return(x)
      }#End colnames_rep



      if (markers) {
        m <- dplyr::group_by(depth, MARKERS) %>%
          dplyr::summarise_at(.tbl = ., .vars = want, .funs = sum, na.rm = TRUE) %>%
          dplyr::mutate_at(.tbl = ., .vars = want, .funs = round, digits = 0) %>%
          dplyr::mutate_at(.tbl = ., .vars = want, .funs = as.integer) %>%
          dplyr::ungroup(.) %>%
          dplyr::rename_all(.funs = list(colnames_rep), total = TRUE) %>%
          dplyr::bind_cols(
            dplyr::group_by(depth, MARKERS) %>%
              dplyr::summarise_at(.tbl = ., .vars = want, .funs = mean, na.rm = TRUE) %>%
              dplyr::mutate_at(.tbl = ., .vars = want, .funs = round, digits = 0) %>%
              dplyr::mutate_at(.tbl = ., .vars = want, .funs = as.integer) %>%
              dplyr::ungroup(.) %>%
              dplyr::rename_all(.funs = list(colnames_rep), mean = TRUE) %>%
              dplyr::select(-MARKERS)
          )
      }
      if (ind) {
        i <- dplyr::group_by(depth, INDIVIDUALS) %>%
          dplyr::summarise_at(.tbl = ., .vars = want, .funs = sum, na.rm = TRUE) %>%
          dplyr::mutate_at(.tbl = ., .vars = want, .funs = round, digits = 0) %>%
          dplyr::mutate_at(.tbl = ., .vars = want, .funs = as.integer) %>%
          dplyr::ungroup(.) %>%
          dplyr::rename_all(.funs = list(colnames_rep), total = TRUE) %>%
          dplyr::bind_cols(
            dplyr::group_by(depth, INDIVIDUALS) %>%
              dplyr::summarise_at(.tbl = ., .vars = want, .funs = mean, na.rm = TRUE) %>%
              dplyr::mutate_at(.tbl = ., .vars = want, .funs = round, digits = 0) %>%
              dplyr::mutate_at(.tbl = ., .vars = want, .funs = as.integer) %>%
              dplyr::ungroup(.) %>%
              dplyr::rename_all(.funs = list(colnames_rep), mean = TRUE) %>%
              dplyr::select(-INDIVIDUALS)
          )
      }
    } else {
      m <- i <- coverage.info <- NULL
    }

  }#End DART counts

  #DArT 1 row and 2rows
  if ("dart" %in% data.source && any(c("2rows", "1row") %in% data.source)) {
    depth <- extract_markers_metadata(
      gds,
      markers.meta.select = c("AVG_COUNT_REF", "AVG_COUNT_SNP"),
      whitelist = TRUE
    )
    if (!is.null(depth)) {
      coverage.info$markers.mean <- as.integer(
        round(depth$AVG_COUNT_REF + depth$AVG_COUNT_SNP, 0))
      coverage.info$ref.mean <- as.integer(round(depth$AVG_COUNT_REF, 0))
      coverage.info$alt.mean <- as.integer(round(depth$AVG_COUNT_SNP, 0))
      depth <- NULL
      markers <- ind <- FALSE
    } else {
      coverage.info <- NULL
      markers <- ind <- FALSE
    }

  }#End DART 1row and 2 rows


  # Update
  if (update.gds && markers) {
    markers.meta <- extract_markers_metadata(gds, whitelist = FALSE)
    m.levels <- markers.meta$MARKERS
    not.wanted <- c("COVERAGE_MEAN", "REF_DEPTH_MEAN", "ALT_DEPTH_MEAN", "COVERAGE_TOTAL", "REF_DEPTH_TOTAL", "ALT_DEPTH_TOTAL")
    if (any(not.wanted %in% colnames(markers.meta))) {
      suppressWarnings(markers.meta %<>% dplyr::select(-dplyr::one_of(not.wanted)))
    }
    common.col <- intersect(colnames(m), colnames(markers.meta))
    suppressWarnings(markers.meta %<>% dplyr::left_join(m, by = common.col))
    not.wanted <- NULL
    if (!identical(markers.meta$MARKERS, m.levels)) {
      rlang::abort("Problem with order of markers: contact author")
    }
    m.levels <- NULL
    update_radiator_gds(gds = gds,
                        node.name = "markers.meta",
                        value = markers.meta)
    markers.meta <- NULL
  }

  if (markers) {
    coverage.info$markers.tot <- if (rlang::has_name(m, "COVERAGE_TOTAL")) m$COVERAGE_TOTAL
    coverage.info$markers.mean <- if (rlang::has_name(m, "COVERAGE_MEAN")) m$COVERAGE_MEAN
    coverage.info$ref.tot <- if (rlang::has_name(m, "REF_DEPTH_TOTAL")) m$REF_DEPTH_TOTAL
    coverage.info$ref.mean <- if (rlang::has_name(m, "REF_DEPTH_MEAN")) m$REF_DEPTH_MEAN
    coverage.info$alt.tot <- if (rlang::has_name(m, "ALT_DEPTH_TOTAL")) m$ALT_DEPTH_TOTAL
    coverage.info$alt.mean <- if (rlang::has_name(m, "ALT_DEPTH_MEAN")) m$ALT_DEPTH_MEAN
    m <- NULL
  }

  if (update.gds && ind) {
    strata <- extract_individuals_metadata(gds = gds, whitelist = FALSE)
    i.levels <- strata$INDIVIDUALS
    not.wanted <- c("COVERAGE_MEAN", "REF_DEPTH_MEAN", "ALT_DEPTH_MEAN", "COVERAGE_TOTAL", "REF_DEPTH_TOTAL", "ALT_DEPTH_TOTAL")
    if (any(not.wanted %in% colnames(strata))) {
      suppressWarnings(strata %<>% dplyr::select(-dplyr::one_of(not.wanted)))
    }
    common.col <- intersect(colnames(i), colnames(strata))
    suppressWarnings(strata %<>% dplyr::left_join(i, by = common.col))
    not.wanted <- NULL
    if (!identical(strata$INDIVIDUALS, i.levels)) {
      rlang::abort("Problem with order of individuals: contact author")
    }
    i.levels <- NULL
    update_radiator_gds(gds = gds,
                        node.name = "individuals.meta",
                        value = strata)
    strata <- NULL
  }
  if (ind) {
    coverage.info$ind.cov.tot <- if (rlang::has_name(i, "COVERAGE_TOTAL")) i$COVERAGE_TOTAL
    coverage.info$ind.cov.mean <- if (rlang::has_name(i, "COVERAGE_MEAN")) i$COVERAGE_MEAN
    coverage.info$ind.ref.tot <- if (rlang::has_name(i, "REF_DEPTH_TOTAL")) i$REF_DEPTH_TOTAL
    coverage.info$ind.ref.mean <- if (rlang::has_name(i, "REF_DEPTH_MEAN")) i$REF_DEPTH_MEAN
    coverage.info$ind.alt.tot <- if (rlang::has_name(i, "ALT_DEPTH_TOTAL")) i$ALT_DEPTH_TOTAL
    coverage.info$ind.alt.mean <- if (rlang::has_name(i, "ALT_DEPTH_MEAN")) i$ALT_DEPTH_MEAN
    i <- NULL
  }
  return(coverage.info)
}#End extract_coverage

# sync GDS----------------------------------------------------------------------
#' @title sync_gds
#' @description Synchronize gds with samples and markers. If left NULL, the info
#' is first search in the radiator node, if not found, it goes in the level above.
#' An argument also allows to reset the filters.
#' @param gds The gds object.
#' @param samples (optional, character string). Will sync the gds object/file with
#' these samples. With default, uses the individuals in the radiator node.
#' If not found, goes a level above and uses the individuals in the main GDS.
#' Default: \code{samples = NULL}.
#' @param variant.id (optional, integer string). Will sync the gds object/file with
#' these variant.id With default, uses the variant.id in the radiator node.
#' If not found, goes a level above and uses the variant.id in the main GDS.
#' Default: \code{variant.id = NULL}.
#' @param reset.gds (optional, logical) Default: \code{reset.gds = FALSE}.
#' @param reset.filters.m (optional, logical) To reset only markers/variant. Default: \code{reset.filters.m = FALSE}.
#' @param reset.filters.i (optional, logical) To reset only individuals. Default: \code{reset.filters.i = FALSE}.
#' @param verbose (optional, logical) Default: \code{verbose = FALSE}.
#' @rdname sync_gds
# @keywords internal
#' @seealso \code{\link{sync_gds}}, \code{\link{list_filters}}.
#' @export
sync_gds <- function(
  gds,
  samples = NULL,
  variant.id = NULL,
  reset.gds = FALSE,
  reset.filters.m = FALSE,
  reset.filters.i = FALSE,
  verbose = FALSE) {

  if (reset.gds || reset.filters.m || reset.filters.i) {
    if (reset.gds) {
      SeqArray::seqResetFilter(
        object = gds,
        sample = TRUE,
        variant = TRUE,
        verbose = verbose
      )
      reset.filters.m <- reset.filters.i <- TRUE
    }

    if (reset.filters.m) {
      m <- extract_markers_metadata(gds = gds, whitelist = FALSE) %>%
        dplyr::mutate(FILTERS = "whitelist")
      update_radiator_gds(gds = gds, node.name = "markers.meta", value = m, sync = FALSE)
      SeqArray::seqSetFilter(
        object = gds, variant.id = m$VARIANT_ID, verbose = FALSE)
    }

    if (reset.filters.i) {
      i <- extract_individuals_metadata(gds = gds, whitelist = FALSE) %>%
        dplyr::mutate(FILTERS = "whitelist")
      update_radiator_gds(gds = gds, node.name = "individuals.meta", value = i, sync = FALSE)
      SeqArray::seqSetFilter(
        object = gds, sample.id = i$INDIVIDUALS, verbose = FALSE)
    }

  } else {
    if (is.null(variant.id)) {
      if (verbose) message("synchronizing GDS with current variant.id")
      variant.id <- extract_markers_metadata(
        gds = gds,
        markers.meta.select = "VARIANT_ID",
        whitelist = TRUE,
        verbose = verbose
      ) %$% VARIANT_ID
    } else {
      if (verbose) message("synchronizing GDS with provided variant.id")
      variant.id <- as.integer(variant.id)
    }

    if (is.null(samples)) {
      if (verbose) message("synchronizing GDS with current samples")
      samples <- extract_individuals_metadata(
        gds = gds,
        ind.field.select = "INDIVIDUALS",
        whitelist = TRUE,
        verbose = verbose
      ) %$% INDIVIDUALS
    } else {
      if (verbose) message("synchronizing GDS with provided samples")
    }

    SeqArray::seqSetFilter(
      object = gds, sample.id = samples, variant.id = variant.id, verbose = verbose)
  }
  # return(gds)
}#End sync_gds

# list_filters------------------------------------------------------------------
#' @name list_filters
#' @rdname list_filters
#' @title List current active filters (individuals and markers) in radiator GDS object.
#' @description List current active filters (individuals and markers) in radiator GDS object.
#' @inheritParams radiator_common_arguments
#' @export
#' @examples
#' \dontrun{
#' # List active filters for individuals and markers
#' list_filters(gds)
#' }
#' @seealso \code{\link{sync_gds}}, \code{\link{list_filters}}.
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

list_filters <- function(gds) {
  radiator_packages_dep("SeqVarTools", cran = FALSE, bioc = TRUE)

  data.type <- radiator::detect_genomic_format(gds)

  if (!data.type %in% c("SeqVarGDSClass", "gds.file")) {
    rlang::abort("Input not supported for this function: read function documentation")
  }
  if (data.type == "gds.file") {
    gds <- radiator::read_rad(gds)
    data.type <- "SeqVarGDSClass"
  }
  i <- extract_individuals_metadata(gds = gds, ind.field.select = "FILTERS", blacklist = FALSE) %>%
    dplyr::count(FILTERS, FILTERS) #%>% readr::write_tsv(x = i, path = "filters.individuals.tsv")
  i %<>% dplyr::filter(FILTERS != "whitelist")
  message("Number of filters for individuals: ", nrow(i))
  if (nrow(i) > 0) {
    message("Filter(s): ")
    message(stringi::stri_join(i$FILTERS, collapse = "\n"))
  }
  m <- extract_markers_metadata(gds = gds, markers.meta.select = "FILTERS", blacklist = FALSE) %>%
    dplyr::count(FILTERS, FILTERS) #%>% readr::write_tsv(x = i, path = "filters.individuals.tsv")
  m %<>% dplyr::filter(FILTERS != "whitelist")
  message("\nNumber of filters for markers: ", nrow(m))
  if (nrow(m) > 0) {
    message("Filter(s): ")
    message(stringi::stri_join(m$FILTERS, collapse = "\n"))
  }
}#list_filters

# reset_filters----------------------------------------------------------------------
#' @name reset_filters
#' @rdname reset_filters
#' @title Reset filters (individuals and markers) in radiator GDS object.
#' @description List current active filters in radiator GDS object.
#' Reset specific filters or all filters at once.
#' @inheritParams radiator_common_arguments
#' @param list.filters filters (logical, optional) List current active filters for individuals and markers.
#' Default: \code{list.filters = TRUE}.
#' @param reset.all (logical, optional) Reset all individuals and markers filters.
#' Default: \code{reset.all = FALSE}.
#' @param filter.individuals.missing (logical, optional)
#' Default: \code{filter.individuals.missing = FALSE}.
#' @param filter.individuals.heterozygosity (logical, optional)
#' Default: \code{filter.individuals.heterozygosity = FALSE}.
#' @param filter.individuals.coverage.total (logical, optional)
#' Default: \code{filter.individuals.coverage.total = FALSE}.
#' @param detect.mixed.genomes (logical, optional)
#' Default: \code{detect.mixed.genomes = FALSE}.
#' @param detect.duplicate.genomes (logical, optional)
#' Default: \code{detect.duplicate.genomes = FALSE}.

#' @param filter.reproducibility (logical, optional)
#' Default: \code{filter.reproducibility = FALSE}.
#' @param filter.monomorphic (logical, optional)
#' Default: \code{filter.monomorphic = FALSE}.
#' @param filter.common.markers (logical, optional)
#' Default: \code{filter.common.markers = FALSE}.
#' @param filter.mac (logical, optional)
#' Default: \code{filter.mac = FALSE}.
#' @param filter.mean.coverage (logical, optional)
#' Default: \code{filter.mean.coverage = FALSE}.
#' @param filter.genotyping (logical, optional)
#' Default: \code{filter.genotyping = FALSE}.
#' @param filter.snp.position.read (logical, optional)
#' Default: \code{filter.snp.position.read = FALSE}.
#' @param filter.snp.number (logical, optional)
#' Default: \code{filter.snp.number = FALSE}.
#' @param filter.short.ld (logical, optional)
#' Default: \code{filter.short.ld = FALSE}.
#' @param filter.long.ld (logical, optional)
#' Default: \code{filter.long.ld = FALSE}.
#' @param filter.hwe (logical, optional)
#' Default: \code{filter.hwe = FALSE}.
#' @param filter.whitelist (logical, optional)
#' Default: \code{filter.whitelist = FALSE}.
#' @keywords internal
#' @export
#'
#'
#' @examples
#' \dontrun{
#' # List active filters for individuals and markers
#' reset_filters(gds)
#'
#' # You changed your decision concerning the genotyping threshold or
#' # entered the wrong one:
#' reset_filters(gds, filter.genotyping = TRUE)
#' }
#' @seealso \code{\link{sync_gds}}, \code{\link{list_filters}}.
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

reset_filters <- function(
  gds,
  list.filters = TRUE,
  reset.all = FALSE,
  filter.individuals.missing = FALSE,
  filter.individuals.heterozygosity = FALSE,
  filter.individuals.coverage.total = FALSE,
  detect.mixed.genomes = FALSE,
  detect.duplicate.genomes = FALSE,
  filter.reproducibility = FALSE,
  filter.monomorphic = FALSE,
  filter.common.markers = FALSE,
  filter.mac = FALSE,
  filter.mean.coverage = FALSE,
  filter.genotyping = FALSE,
  filter.snp.position.read = FALSE,
  filter.snp.number = FALSE,
  filter.short.ld = FALSE,
  filter.long.ld = FALSE,
  filter.hwe = FALSE,
  filter.whitelist = FALSE
) {
  radiator_packages_dep(package = "SeqVarTools", cran = FALSE, bioc = TRUE)
  data.type <- radiator::detect_genomic_format(gds)

  if (!data.type %in% c("SeqVarGDSClass", "gds.file")) {
    rlang::abort("Input not supported for this function: read function documentation")
  }
  if (data.type == "gds.file") {
    gds <- radiator::read_rad(gds)
    data.type <- "SeqVarGDSClass"
  }


  if (list.filters) list_filters(gds)

  # reset all ------------------------------------------------------------------
  if (reset.all) {
    markers.meta <- extract_markers_metadata(gds = gds, whitelist = FALSE) %>%
      dplyr::mutate(FILTERS = "whitelist")
    update_radiator_gds(gds = gds, node.name = "markers.meta", value = markers.meta, sync = FALSE)
    SeqArray::seqSetFilter(object = gds, variant.id = markers.meta$VARIANT_ID, verbose = FALSE)

    individuals <- extract_individuals_metadata(gds = gds, whitelist = FALSE) %>%
      dplyr::mutate(FILTERS = "whitelist")
    update_radiator_gds(gds = gds, node.name = "individuals.meta", value = individuals, sync = FALSE)
    SeqArray::seqSetFilter(object = gds, sample.id = individuals$INDIVIDUALS, verbose = FALSE)
  }
  # reset individuals ----------------------------------------------------------
  reset.i <- list()
  if (filter.individuals.missing) {
    reset.i$filter.individuals.missing <- "filter.individuals.missing"
  }

  if (filter.individuals.heterozygosity) {
    reset.i$filter.individuals.heterozygosity <- "filter.individuals.heterozygosity"
  }
  if (filter.individuals.coverage.total) {
    reset.i$filter.individuals.coverage.total <- "filter.individuals.coverage.total"
  }
  if (detect.mixed.genomes) {
    reset.i$detect.mixed.genomes <- "detect.mixed.genomes"
  }

  if (detect.duplicate.genomes) {
    reset.i$detect.duplicate.genomes <- "detect.duplicate.genomes"
  }

  if (length(reset.i) > 0) {
    message("Resetting ", length(reset.i), " filters for individuals")
    individuals <- extract_individuals_metadata(gds = gds, whitelist = FALSE) %>%
      dplyr::mutate(FILTERS = dplyr::if_else(FILTERS %in% reset.i, "whitelist", FILTERS))
    # update
    update_radiator_gds(gds = gds, node.name = "individuals.meta", value = individuals, sync = FALSE)
    #reset
    SeqArray::seqSetFilter(object = gds, sample.id = individuals$INDIVIDUALS, verbose = FALSE)
  }

  # reset markers: -------------------------------------------------------------
  reset.m <- list()

  if (filter.reproducibility) {
    reset.m$filter.dart.reproducibility <- "filter.dart.reproducibility"
  }

  if (filter.monomorphic) {
    reset.m$filter.monomorphic <- "filter.monomorphic"
  }

  if (filter.common.markers) {
    reset.m$filter.common.markers <- "filter.common.markers"
  }

  if (filter.mac) {
    reset.m$filter.mac <- "filter.mac"
  }

  if (filter.mean.coverage) {
    reset.m$filter.mean.coverage <- "filter.mean.coverage"
  }

  if (filter.genotyping) {
    reset.m$filter.genotyping <- "filter.genotyping"
  }

  if (filter.snp.position.read) {
    reset.m$filter.snp.position.read <- "filter.snp.position.read"
  }

  if (filter.snp.number) {
    reset.m$filter.snp.number <- "filter.snp.number"
  }

  if (filter.short.ld) {
    reset.m$filter.short.ld <- "filter.short.ld"
  }

  if (filter.long.ld) {
    reset.m$filter.long.ld <- "filter.long.ld"
  }

  if (filter.hwe) {
    reset.m$filter.hwe <- "filter.hwe"
  }

  if (filter.whitelist) {
    reset.m$filter.whitelist <- "filter.whitelist"
  }

  if (length(reset.m) > 0) {
    message("Resetting ", length(reset.m), " filters for markers")
    markers.meta <- extract_markers_metadata(gds = gds, whitelist = FALSE) %>%
      dplyr::mutate(FILTERS = dplyr::if_else(FILTERS %in% reset.m, "whitelist", FILTERS))
    # update
    update_radiator_gds(gds = gds, node.name = "markers.meta", value = markers.meta, sync = FALSE)
    #reset
    SeqArray::seqSetFilter(object = gds, variant.id = markers.meta$VARIANT_ID, verbose = FALSE)
  }

}#reset_filters

# summary_gds ------------------------------------------------------------------
#' @name summary_gds
#' @title summary_gds
#' @description Summary of gds object or file: number of samples and markers
#' @inheritParams radiator_common_arguments
#' @param check.sync (logical, optional) Check if GDS object is in sync with
#' radiator metada.
#' @rdname summary_gds
# @keywords internal
#' @export
summary_gds <- function(gds, check.sync = FALSE, verbose = TRUE) {
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

  if (check.sync) {
    id.filtered <- nrow(
      extract_individuals_metadata(
        gds = gds,
        ind.field.select = "INDIVIDUALS",
        whitelist = TRUE
      )
    )

    markers.filtered <- nrow(
      extract_markers_metadata(
        gds = gds,
        markers.meta.select = "VARIANT_ID",
        whitelist = TRUE
      )
    )

    if (!identical(n.ind, id.filtered) || !identical(n.markers, markers.filtered)) {
      if (verbose) message("\nGDS file not in sync with data")
    } else {
      if (verbose) message("\nGDS file is in sync with data")
    }
  }#End check sync

  invisible(x = list(n.ind = n.ind, n.markers = n.markers))
  # return(res = list(n.ind = n.ind, n.markers = n.markers))
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
  path.folder = NULL,
  plot = TRUE,
  digits = 6,
  file.date = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE
) {

  ## TEST
  # missing = TRUE
  # heterozygosity = TRUE
  # coverage = TRUE
  # plot = TRUE
  # digits = 6

  if (is.null(path.folder)) path.folder <- getwd()
  res <- list() # return result in this list
  if (is.null(file.date)) file.date <- format(Sys.time(), "%Y%m%d@%H%M")

  # subsampling ----------------------------------------------------------------
  if (!is.null(subsample)) {
    SeqArray::seqSetFilter(object = gds,
                           variant.id = subsample,
                           action = "push+set",
                           verbose = FALSE)
  }

  id.info <- extract_individuals_metadata(gds, whitelist = TRUE)

  # missing --------------------------------------------------------------------
  if (missing) {
    # info
    id.info %<>%
      dplyr::mutate(
        MISSING_PROP = round(
          SeqArray::seqMissing(
            gdsfile = gds, per.variant = FALSE,
            .progress = TRUE, parallel = parallel.core
          ) , digits)
      )
    # stats
    id.stats.m <- tibble_stats(
      x = id.info$MISSING_PROP,
      group = "missing genotypes")
  } else {
    missing <- FALSE
    id.stats.m <- NULL
  }# End missing


  # heterozygosity -------------------------------------------------------------
  if (heterozygosity) {
    # info
    id.info %<>%
      dplyr::mutate(
        HETEROZYGOSITY = round(SeqVarTools::heterozygosity(
          gdsobj = gds, margin = "by.sample", use.names = FALSE
        ), digits)
      )

    id.info$HETEROZYGOSITY[is.na(id.info$HETEROZYGOSITY)] <- 0

    # stats
    id.stats.h <- tibble_stats(
      x = id.info$HETEROZYGOSITY,
      group = "heterozygosity")
  } else {
    heterozygosity <- FALSE
    id.stats.h <- NULL
  }

  # coverage -------------------------------------------------------------------
  if (coverage) {
    # info
    dp <- extract_coverage(gds, markers = FALSE)

    if (!is.null(dp) && "ind.cov.tot" %in% names(dp)) {
      id.info %<>%
        dplyr::mutate(
          COVERAGE_TOTAL = dp$ind.cov.tot,
          COVERAGE_MEAN = dp$ind.cov.mean
        )

      id.info$COVERAGE_MEAN[is.na(id.info$COVERAGE_MEAN)] <- 0


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
      coverage <- FALSE
      id.stats.c <- NULL
    }
  } else {
    coverage <- FALSE
    id.stats.c <- NULL
  } # End coverage

  # Generate the stats
  id.stats <- dplyr::bind_rows(id.stats.m, id.stats.h, id.stats.c)
  id.levels <- c("total coverage", "mean coverage", "missing genotypes", "heterozygosity")
  id.stats$GROUP <- factor(x = id.stats$GROUP, levels = id.levels, ordered = TRUE)
  id.stats$GROUP <- droplevels(x = id.stats$GROUP)

  id.stats.filename <- stringi::stri_join("individuals.qc.stats_", file.date, ".tsv")
  if (!is.null(id.info)) readr::write_tsv(x = id.info, path = file.path(path.folder, id.stats.filename))
  id.stats.filename <- stringi::stri_join("individuals.qc.stats.summary_", file.date, ".tsv")
  if (!is.null(id.stats))readr::write_tsv(x = id.stats, path = file.path(path.folder, id.stats.filename))
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
      if (stats::sd(id.info$COVERAGE_MEAN) != 0) {
        cc <- ceiling(stats::cor(id.info$COVERAGE_TOTAL, id.info$COVERAGE_MEAN,use = "pairwise.complete.obs") * 100) / 100
      } else {
        cc <- "NA"
        message("Note: mean individual coverage SD = 0")
        message("correlation with total coverage is not possible")
      }
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

    if (missing) {
      missing.out <- stringi::stri_join("Missing genotypes outlier: ", id.stats[[1, 9]])
      corr.info <- stringi::stri_join(missing.out, "\n", corr.info)
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
  force.stats = TRUE,
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

  if (is.null(path.folder)) path.folder <- getwd()
  if (is.null(file.date)) file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  data.source <- extract_data_source(gds)

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

  info <- extract_markers_metadata(gds = gds, whitelist = TRUE)

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
      if ("dart" %in% data.source) {
        if ("counts" %in% data.source) {
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
  # test <- SeqArray::seqMissing(
  #   gdsfile = data,
  #   per.variant = TRUE, .progress = TRUE, parallel = parallel.core)


  if (missing) {
    if (!rlang::has_name(info, "MISSING_PROP") || force.stats) {
      # if (rlang::has_name(info, "MISSING_PROP")) info  %<>% dplyr::select(-MISSING_PROP)
      info$MISSING_PROP <- SeqArray::seqMissing(
        gdsfile = gds,
        per.variant = TRUE,
        .progress = TRUE,
        parallel = parallel.core
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
  }
  # else {
  # radiator.gds <- gdsfmt::index.gdsn(node = gds, path = "radiator", silent = TRUE)
  # Update GDS
  # if (!is.null(radiator.gds)) {
  # gdsfmt::add.gdsn(
  # node = radiator.gds,
  # name = "markers.meta",
  # val = info,
  # replace = TRUE,
  # compress = "ZIP_RA",
  # closezip = TRUE)
  # }
  # }

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
    id.select,
    gds,
    markers.meta,
    parallel.core = parallel::detectCores() - 1
  ) {
    SeqArray::seqSetFilter(
      object = gds,
      sample.id = id.select$INDIVIDUALS,
      action = "set",
      verbose = FALSE)
    n.ind <- length(id.select$INDIVIDUALS)
    if (n.ind < parallel.core) {
      parallel.core.temp <- n.ind
    } else {
      parallel.core.temp = parallel.core
    }
    res <- markers.meta %>%
      dplyr::mutate(
        MISSING_PROP = SeqArray::seqMissing(
          gdsfile = gds,
          per.variant = TRUE,
          .progress = TRUE,
          parallel = parallel.core.temp
        )
      )

    mis_many_markers <- function(threshold, x) {
      nrow(dplyr::filter(x, MISSING_PROP <= threshold))
    }#End how_many_markers

    helper.table <- tibble::tibble(MISSING_PROP = seq(0.1, 1, 0.1)) %>%
      dplyr::mutate(
        WHITELISTED_MARKERS = purrr::map_int(
          .x = seq(0.1, 1, 0.1), .f = mis_many_markers, x = res)
      )
    return(helper.table)
  }#End missing_pop

  sample.bk <- extract_individuals_metadata(
    gds = gds,
    ind.field.select = "INDIVIDUALS",
    whitelist = TRUE
  ) %$% INDIVIDUALS

  markers.meta <- extract_markers_metadata(
    gds = gds,
    markers.meta.select = "MARKERS",
    whitelist = TRUE)
  n.markers <- nrow(markers.meta)

  res <- strata %>%
    # dplyr::group_by(STRATA) %>%
    # tidyr::nest(data = ., .key = id.select) %>%
    tidyr::nest(.data = ., id.select = "INDIVIDUALS") %>%
    dplyr::mutate(
      MISSING_POP = purrr::map(
        .x = .$id.select,
        .f = missing_pop,
        gds = gds,
        markers.meta = markers.meta,
        parallel.core = parallel.core),
      id.select = NULL
    ) %>%
    tidyr::unnest(data = ., MISSING_POP) %>%
    dplyr::mutate(BLACKLISTED_MARKERS = n.markers - WHITELISTED_MARKERS)

  # reset filter
  SeqArray::seqSetFilter(
    object = gds,
    sample.id = sample.bk,
    action = "set",
    verbose = FALSE)

  return(res)
}#End missing_per_pop


# write_gds --------------------------------------------------------------------
# write a GDS object from a tidy data frame

#' @name write_gds
#' @title Write a GDS object from a tidy data frame
#' @description Write a Genomic Data Structure (GDS) file format
#' \href{https://github.com/zhengxwen/gdsfmt}{gdsfmt}
#' and object of class \code{SeqVarGDSClass} from a tidy data frame.
#'
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.

#' @param data.source (optional, character) The name of the software that
#' generated the data. e.g. \code{data.source = "Stacks v.2.4"}.
#' Default: \code{data.source = NULL}.

#' @param filename (optional) The file name of the Genomic Data Structure (GDS) file.
#' radiator will append \code{.gds.rad} to the filename.
#' If filename chosen is already present in the
#' working directory, the default \code{radiator_datetime.gds.rad} is chosen.
#' Default: \code{filename = NULL}.

#' @param open (optional, logical) Open or not the radiator connection.
#' Default: \code{open = TRUE}.

#' @inheritParams tidy_genomic_data

#' @return An object in the global environment of class \code{SeqVarGDSClass} and
#' a file in the working directory.

#' @export
#' @rdname write_gds
#' @references Zheng X, Levine D, Shen J, Gogarten SM, Laurie C, Weir BS.
#' A high-performance computing toolset for relatedness and principal component
#' analysis of SNP data. Bioinformatics. 2012;28: 3326-3328.
#' doi:10.1093/bioinformatics/bts606
#' @references Zheng X, Gogarten S, Lawrence M, Stilp A, Conomos M, Weir BS,
#' Laurie C, Levine D (2017). SeqArray -- A storage-efficient high-performance
#' data format for WGS variant calls.
#' Bioinformatics.

#' @examples
#' \dontrun{
#' require(SeqVarTools)
#' data.gds <- radiator::write_gds(data = "shark.rad")
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_gds <- function(
  data,
  data.source = NULL,
  filename = NULL,
  open = TRUE,
  verbose = TRUE
) {
  timing <- proc.time()

  ## testing
  # filename = NULL
  # verbose = TRUE

  # Check that SeqVarTools is installed (it requires automatically: SeqArray and gdsfmt)
  radiator_packages_dep(package = "SeqVarTools", cran = FALSE, bioc = TRUE)

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("Input file missing")

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    data <- radiator::tidy_wide(data = data, import.metadata = TRUE)
  }

  biallelic <- radiator::detect_biallelic_markers(data)
  if (!biallelic) rlang::abort("Biallelic data required")

  # GDS prep -------------------------------------------------------------------
  if (verbose) message("Generating radiator GDS format...")
  n.markers <- length(unique(data$MARKERS))
  n.ind <- length(unique(data$INDIVIDUALS))

  strata <-
    suppressWarnings(
      data %>%
        dplyr::select(INDIVIDUALS, dplyr::one_of(c("STRATA", "POP_ID", "TARGET_ID")))
    )
  if (tibble::has_name(strata, "POP_ID")) strata  %<>% dplyr::rename(STRATA = POP_ID)

  want <- c("VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS", "COL", "REF",
            "ALT", "CALL_RATE", "REP_AVG", "AVG_COUNT_REF",
            "AVG_COUNT_SNP", "ONE_RATIO_REF", "ONE_RATIO_SNP", "SEQUENCE")
  markers.meta <-
    suppressWarnings(
      data %>%
        dplyr::select(
          dplyr::one_of(want))
    )
  if (ncol(markers.meta) == 0) markers.meta <- NULL

  want <- c(
    "GT", "GT_VCF", "GT_VCF_NUC",
    "READ_DEPTH", "ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH",
    "GQ",
    "GL_HOM_REF", "GL_HET", "GL_HOM_ALT",
    "DP", "AD", "GL", "PL", "HQ", "GOF", "NR", "NV", "RO", "QR", "AO", "QA")
  genotypes.meta <-
    suppressWarnings(
      data %>%
        dplyr::select(
          dplyr::one_of(want))
    )
  want <- NULL
  if (ncol(genotypes.meta) == 0) genotypes.meta <- NULL

  data.gds <- radiator_gds(
    genotypes = data,
    strata = strata,
    biallelic = TRUE,
    markers.meta = markers.meta,
    genotypes.meta = genotypes.meta,
    filename = filename,
    data.source = data.source,
    verbose = verbose
  )
  markers.meta <- genotypes.meta <- data <- strata <- NULL

  # RETURN ---------------------------------------------------------------------
  if (verbose) message("Number of individuals: ", n.ind)
  if (verbose) message("Number of markers: ", n.markers)

  timing <- proc.time() - timing
  if (verbose) message("\nConversion time: ", round(timing[[3]]), " sec\n")

  if (open) data.gds <- read_rad(data.gds, verbose = FALSE)
  return(data.gds)
} # End write_gds

# write_snprelate --------------------------------------------------------------
#' @name write_snprelate
#' @title Write a SNPRelate object from a tidy data frame
#' @description Write a \href{https://github.com/zhengxwen/SNPRelate}{SNPRelate}
#' object from a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.
#' \strong{The genotypes are biallelic.}

#' @param biallelic (logical, optional) If you already know that the data is
#' biallelic use this argument to speed up the function.
#' Default: \code{biallelic = TRUE}.

#' @param filename (optional) The file name of the Genomic Data Structure (GDS) file.
#' radiator will append \code{.gds} to the filename.
#' If filename chosen is already present in the
#' working directory, the default \code{radiator_snprelate_datetime.gds} is chosen.
#' Default: \code{filename = NULL}.

#' @inheritParams tidy_genomic_data

#' @export
#' @rdname write_snprelate

#' @references Zheng X, Levine D, Shen J, Gogarten SM, Laurie C, Weir BS.
#' A high-performance computing toolset for relatedness and principal component
#' analysis of SNP data. Bioinformatics. 2012;28: 3326-3328.
#' doi:10.1093/bioinformatics/bts606

#' @seealso \href{https://github.com/zhengxwen/SNPRelate}{SNPRelate}

#' @return An object in the global environment of class
#' \code{"SNPGDSFileClass", "gds.class"} and
#' a file in the working directory.

#' @examples
#' \dontrun{
#' require(SNPRelate)
#' data.gds <- radiator::write_snprelate(data = "shark.rad")
#' }


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_snprelate <- function(data, biallelic = TRUE, filename = NULL, verbose = TRUE) {

  #TEST
  # biallelic = TRUE
  # filename = NULL
  # verbose = TRUE

  # Check that snprelate is installed
  radiator_packages_dep(package = "SNPRelate", cran = FALSE, bioc = TRUE)

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("Input file missing")
  if (verbose) message("Generating SNPRelate object/file...")

  # Filename -------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")

  if (!is.null(filename)) filename <- stringi::stri_join(filename, "_snprelate")
  filename <- generate_filename(name.shortcut = filename, extension = "gds.rad")

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    data <- radiator::tidy_wide(data = data, import.metadata = TRUE)
  }

  biallelic <- radiator::detect_biallelic_markers(data)
  if (!biallelic) rlang::abort("Biallelic data required")

  # Generate GDS format --------------------------------------------------------
  if (!rlang::has_name(data, "POP_ID")) data %<>% dplyr::mutate(POP_ID = "pop")
  data %<>% dplyr::arrange(MARKERS, INDIVIDUALS, POP_ID) %>%
    dplyr::mutate(VARIANT_ID = as.integer(factor(MARKERS))) %>%
    dplyr::ungroup(.)

  # MARKERS META WORK ----------------------------------------------------------
  markers.meta <- separate_markers(
    data = data,
    sep = "__",
    markers.meta.all.only = TRUE,
    biallelic = biallelic,
    verbose = verbose) %>%
    dplyr::arrange(VARIANT_ID)
  n.markers <- nrow(markers.meta)

  notwanted <- c("MARKERS", "CHROM", "LOCUS", "POS", "COL", "REF", "ALT")
  suppressWarnings(data %<>% dplyr::select(-dplyr::one_of(notwanted)))

  # genotypes metadata
  want <- c("ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH", "READ_DEPTH", "GL_HOM_REF",
            "GL_HET", "GL_HOM_ALT", "GQ", "DP", "AD", "GL", "PL", "GQ", "HQ", "GOF", "NR", "NV")
  if (TRUE %in% rlang::has_name(data, want)) {
    genotypes.meta <- suppressWarnings(data %>% dplyr::select(dplyr::one_of(want)))
    suppressWarnings(data %<>% dplyr::select(-dplyr::one_of(want)))
  } else {
    genotypes.meta <- NULL
  }

  # strata
  strata <- dplyr::distinct(data, POP_ID, INDIVIDUALS) %>%
    dplyr::rename(STRATA = POP_ID) %>%
    dplyr::arrange(INDIVIDUALS)

  # prep data
  data %<>% dplyr::select(VARIANT_ID, INDIVIDUALS, GT_BIN) %>%
    data.table::as.data.table(.) %>%
    data.table::dcast.data.table(
      data = .,
      formula = VARIANT_ID ~ INDIVIDUALS,
      value.var = "GT_BIN"
    ) %>%
    dplyr::arrange(VARIANT_ID) %>%
    magrittr::set_rownames(x = ., value = .$VARIANT_ID) %>%
    dplyr::select(-VARIANT_ID) %>%
    data.matrix(.) %>%
    # 4 steps for SNPRelate genotype coding (change from ALT to REF dosage)
    magrittr::inset(is.na(.), 3L) %>%
    magrittr::inset(. == 0L, 9L) %>%
    magrittr::inset(. == 2L, 0L) %>%
    magrittr::inset(. == 9L, 2L)


  SNPRelate::snpgdsCreateGeno(
    gds.fn = filename$filename,
    genmat = data,
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
  data <- SNPRelate::snpgdsOpen(filename$filename, readonly = FALSE, allow.fork = TRUE)

  if (verbose) message("SNPRelate GDS: ", filename$filename.short)

  # timing <- proc.time() - timing
  # if (verbose) message("\nConversion time: ", round(timing[[3]]), " sec\n")
  return(data)
} # End write_snprelate
