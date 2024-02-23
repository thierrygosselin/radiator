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
    markers.meta %<>% dplyr::select(-dplyr::any_of(strata$INDIVIDUALS))
  }

  # check <- markers.meta$MARKERS
  # generate all the markers metadata


    # test <- markers.meta  %>%
  #   dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
  #   separate_markers(
  #     data = .,
  #     sep = "__",
  #     markers.meta.all.only = TRUE,
  #     biallelic = TRUE,
  #     verbose = verbose)

  # data.bk <- data
  # data <- data.bk
  # data <- markers.meta  %>%
  #   dplyr::distinct(MARKERS, .keep_all = TRUE)
  # sep = "__"
  # markers.meta.all.only = TRUE
  # markers.meta.lists.only = FALSE
  # biallelic = TRUE
  # generate.ref.alt = FALSE
  # parallel.core = parallel::detectCores() - 1

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
  # gdsfmt::closefn.gds(gdsfile = data.gds)
  # data.gds <- gdsfmt::openfn.gds(filename = filename.gds$filename, readonly = FALSE)
  # ?gdsfmt::unload.gdsn
  # gdsfmt::read.gdsn(node = genotype.node)
  # gdsfmt::readmode.gdsn(node = genotype.node)
  # gdsfmt::compression.gdsn(node = genotype.node, compress="ZIP_RA")
  # gdsfmt::compression.gdsn(node = genotype.node, compress="")
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
  return(data.gds)
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
  dimnames(genotypes) <- list(allele = NULL, sample = NULL, variant = NULL)

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
      dplyr::mutate(dplyr::across(.cols = mk.col, .fns = as.character)) %>%
      dplyr::arrange(MARKERS) %>%
      dplyr::mutate(VARIANT_ID = as.integer(factor(MARKERS)))
  }
  want <- c("VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS", "COL", "REF",
            "ALT")
  markers.meta <- x %>%
    dplyr::select(dplyr::any_of(want)) %>%
    dplyr::distinct(.) %>%
    dplyr::arrange(VARIANT_ID)

  notwanted <- c("MARKERS", "CHROM", "LOCUS", "POS", "COL", "REF","ALT")
  x %<>%
    dplyr::select(-dplyr::any_of(notwanted)) %>%
    dplyr::rename(STRATA = tidyselect::any_of("POP_ID"))

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
    strip.rad = FALSE,
    parallel.core = parallel::detectCores() - 1,
    close.gds = FALSE,
    ...
) {
  if (is.null(individuals)) {
    individuals <- extract_individuals_metadata(gds = gds, whitelist = TRUE, verbose = FALSE)
  }
  if (!rlang::has_name(individuals, "STRATA") && !rlang::has_name(individuals, "POP_ID")) {
    if (strip.rad) {
      individuals %<>% dplyr::mutate(STRATA_SEQ = 1L)
    } else {
      individuals %<>% dplyr::mutate(STRATA = 1L)
    }
  }
  if (pop.id) individuals %<>% dplyr::rename(POP_ID = tidyselect::any_of("STRATA"))


  if (wide) markers.meta.select <- "MARKERS"
  if (is.null(markers.meta.select) && strip.rad) markers.meta.select <- "M_SEQ"

  if (is.null(markers.meta)) {
    markers.meta <- extract_markers_metadata(
      gds = gds,
      markers.meta.select = markers.meta.select,
      whitelist = TRUE,
      verbose = TRUE
    )
  }

  # added 20210616

  weird.locus <- suppressWarnings(length(unique(markers.meta$LOCUS)) <= 1)
  if (weird.locus) {
    message("LOCUS field empty... adding unique id instead")
    markers.meta$LOCUS <- suppressWarnings(markers.meta$VARIANT_ID)
  }

  # colnames and "id" tidy data...
  if (rlang::has_name(markers.meta, "M_SEQ")) {
    cn <- markers.meta %>% dplyr::pull(M_SEQ)
    want.m <- "M_SEQ"
  } else {
    if (!rlang::has_name(markers.meta, "MARKERS")) {
      markers.meta %<>%
        dplyr::mutate(
          dplyr::across(
            .cols = c(CHROM, LOCUS, POS),
            .fns = radiator::clean_markers_names
          )
        ) %>%
        dplyr::mutate(
          MARKERS = stringi::stri_join(CHROM, LOCUS, POS, sep = "__")
        )
    }
    cn <- markers.meta %>% dplyr::pull(MARKERS)
    want.m <- "MARKERS"
  }
  if (rlang::has_name(individuals, "ID_SEQ")) {
    want.id <- "ID_SEQ"
  } else {
    want.id <- "INDIVIDUALS"
  }
  id <- individuals %>% dplyr::select(tidyselect::any_of(want.id))


  # summary_gds(gds)
  tidy.data <- SeqArray::seqGetData(gdsfile = gds, var.name = "$dosage_alt") %>%
    magrittr::set_colnames(x = ., value = cn) %>%
    tibble::as_tibble(x = .) %>%
    tibble::add_column(id, .before = 1)

  if (!wide) {
    tidy.data %<>%
      radiator::rad_long(
        x = .,
        cols = want.id,
        names_to = want.m,
        values_to = "GT_BIN",
        variable_factor = FALSE
      )

    if (strip.rad) {
      tidy.data %<>% dplyr::left_join(individuals %>% dplyr::select(ID_SEQ, STRATA_SEQ), by = "ID_SEQ")
    } else {
      want <- intersect(colnames(tidy.data), colnames(markers.meta))
      tidy.data %<>% dplyr::left_join(markers.meta, by = want)
      want <- intersect(colnames(tidy.data), colnames(individuals))
      tidy.data %<>% dplyr::left_join(individuals, by = want)
    }
  }

  # re-calibration of ref/alt alleles
  # Note to myself: maybe this should be done on the GDS before conversion...
  # for small dataset, it won't matter, but for large ones, this could be the bottleneck
  if (calibrate.alleles && !wide) {
    tidy.data %<>%
      radiator::calibrate_alleles(
        data = .,
        verbose = FALSE
      ) %$% input
  }

  #  close the connection with GDS file
  if (close.gds) {
    SeqArray::seqClose(gds)
    gds <- NULL
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
      id.index <- gdsfmt::index.gdsn(node = gds, path = "radiator/individuals.meta", silent = TRUE)
      if (is.null(id.index)) {
        radiator.node <- FALSE
      } else {
        id.index <- gdsfmt::ls.gdsn(id.index)
      }
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

  # will switch radiator.node to FALSE if returns null
  keep.one <- FALSE

  if (whitelist) blacklist <- FALSE
  if (blacklist) whitelist <- FALSE

  if (radiator.node) {
    markers.index <- gdsfmt::index.gdsn(node = gds, path = "radiator/markers.meta", silent = TRUE)

    if (is.null(markers.index) || length(markers.index) == 0 || length(gdsfmt::ls.gdsn(markers.index)) == 0) {
      radiator.node <- FALSE
    } else {
      markers.index <- gdsfmt::ls.gdsn(markers.index)
    }
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
      markers.meta.select %<>%
        intersect(c("VARIANT_ID", "CHROM", "LOCUS", "POS")) %>%
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
  markers.meta  %<>% dplyr::mutate(dplyr::across(.cols = mk.col, .fns = as.character))

  if (!whitelist && !blacklist && !rlang::has_name(markers.meta, "FILTERS")) {
    markers.meta %<>% dplyr::mutate(FILTERS = "whitelist")
  }
  if (whitelist && rlang::has_name(markers.meta, "FILTERS")) {
    markers.meta %<>%
      dplyr::filter(FILTERS == "whitelist") %>%
      dplyr::select(-FILTERS)
  }
  if (blacklist && rlang::has_name(markers.meta, "FILTERS")) {
    markers.meta %<>% dplyr::filter(FILTERS != "whitelist")# %>%
  }
  return(markers.meta)
}#End extract_markers_metadata


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
  ## TEST
  # genotypes = FALSE
  # radiator.node = TRUE
  # index.only = FALSE
  # sync.markers.individuals = TRUE
  # whitelist = FALSE
  # blacklist = FALSE
  # verbose = TRUE

  keep.one <- FALSE
  if (!radiator.node) return(NULL)
  if (whitelist) blacklist <- FALSE
  if (blacklist) whitelist <- FALSE

  geno.index <- gdsfmt::ls.gdsn(gdsfmt::index.gdsn(
    node = gds, path = "radiator/genotypes.meta", silent = TRUE))
  if (length(geno.index) == 0L) return(geno.index)

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
    if (whitelist && rlang::has_name(genotypes.meta, "FILTERS")) {
      genotypes.meta %<>%
        dplyr::filter(FILTERS == "whitelist") %>%
        dplyr::select(-FILTERS)
    }
    if (blacklist && rlang::has_name(genotypes.meta, "FILTERS")) {
      genotypes.meta %<>%
        dplyr::filter(FILTERS != "whitelist") %>%
        dplyr::select(-FILTERS)
    }

    return(genotypes.meta)
  } else {
    return(geno.index)
  }
}#End extract_genotypes_metadata

# check_coverage--------------------------------------------------------------
#' @title check_coverage
#' @description Check that the coverage info is in the GDS. By default, it will
#' look for the DP and AD info in the FORMAT field.
#' @rdname extract_coverage
#' @param gds The gds object.
#' @param genotypes.metadata.check (optional, logical) Look for already extracted coverage information in the
#' radiator genotypes_metadata field of the GDS.
#' Default: \code{genotypes.metadata.check = FALSE}.
#' @param stacks.haplo.check (optional, logical) stacks haplotypes VCF header is
#' baddly generated. It will say you have Read and allele Depth info, but you don't.
#' Default: \code{stacks.haplo.check = FALSE}.
#' @param dart.check (optional, logical) DArT have different reporting for coverage
#' information.
#' Default: \code{dart.check = FALSE}.
# @keywords internal
#' @export

check_coverage <- function(gds, genotypes.metadata.check = FALSE, stacks.haplo.check = FALSE, dart.check = FALSE){

  got.coverage <- NULL # default

  # stacks haplotype vcf have the info fields for depth in the VCF header
  # but they do not have the info with genotypes...
  # it's laziness from stacks...
  data.source <- extract_data_source(gds)

  if (stacks.haplo.check) {
    biallelic <- radiator::detect_biallelic_markers(data = gds)
    biallelic <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(
      node = gds, path = "radiator/biallelic", silent = TRUE))
    if (!biallelic) {
      # check again to be 100% sure...
      biallelic <- radiator::detect_biallelic_markers(data = gds)
    }
    if (!biallelic && stringi::stri_detect_fixed(str = data.source, pattern = "Stacks"))  got.coverage <- NULL
  }

  # DART
  # DArT count and 1 and 2 rows have different information...
  if (dart.check && "dart" %in% data.source) {
    if (any(c("2rows", "1row") %in% data.source)) {
      got.coverage <- extract_markers_metadata(
        gds,
        markers.meta.select = c("AVG_COUNT_REF", "AVG_COUNT_SNP"),
        whitelist = TRUE
      )
      if (!is.null(got.coverage)) {
        got.coverage <- c("AVG_COUNT_REF", "AVG_COUNT_SNP")
        return(got.coverage)
      }
    }#End DART 1row and 2 rows
  }

  # check if coverage info already extracted and in the genotypes_metadata field
  if (genotypes.metadata.check) {
    geno.index <- extract_genotypes_metadata(
      gds = gds,
      genotypes.meta.select = c("READ_DEPTH", "ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH"),
      index.only = TRUE
    )
    if (length(geno.index) == 0L) {
      got.coverage <- NULL
    } else {
      # this part might generate an error if you actually have genotypes metadata...
      # need to run tests...
      got.coverage <- geno.index
      return(got.coverage)
    }
    geno.index <- NULL
  }

  # Check GDS
  # detect FORMAT fields available
  have <-  SeqArray::seqSummary(
    gdsfile = gds,
    varname = "annotation/format",
    check = "none", verbose = FALSE)$ID

  if (length(have) > 0) {
    want <- c("DP", "AD", "CATG")
    got.coverage <- purrr::keep(.x = have, .p = have %in% want)
    return(got.coverage)
  } else {
    got.coverage <- NULL
  }
  have <- NULL
  return(got.coverage)
}#End check_coverate

# extract_coverage--------------------------------------------------------------
#' @title extract_coverage
#' @description Extract coverage information from a GDS file
#' @rdname extract_coverage
#' @param gds The gds object.
#' @param markers (optional, logical) Default: \code{markers = TRUE}.
#' @param individuals (optional, logical) Default: \code{individuals = TRUE}.
#' @param dp (optional, logical) Default: \code{dp = TRUE}.
#' @param ad (optional, logical) Default: \code{ad = TRUE}.
#' @param coverage.stats (optional, character string). Choice of stats to use with
#' coverage.
#' Default: \code{coverage.stats = c("sum", "mean", "median", "iqr")}.
#' @param subsample.info (optional, double) Default: \code{subsample.info = 1}.
#' The subsample proportion used (e.g. 0.3 or none the default).

# @param verbose (optional, logical) Default: \code{verbose = FALSE}.
# @param exhaustive (optional, logical) With default, will output the total, mean
# median and IQR coverage information for all depth info (AD, DP).
# With \code{exhaustive = FALSE}, only the total and mean coverage for DP is extracted.
# Default: \code{exhaustive = TRUE}.
# @param total.cov (optional, logical) For individual or markers total coverage.
# Default: \code{total.cov = FALSE}.
# @param mean.cov (optional, logical) For individual or markers mean coverage.
# Default: \code{mean.cov = TTRUE}.
# @param median.cov (optional, logical) For individual or markers median coverage.
# Default: \code{median.cov = FALSE}.
# @param iqr.cov (optional, logical) For individual or Inter-Quartile Range (iqr) coverage.
# Default: \code{iqr.cov = FALSE}.
# @param update.gds (optional, logical) Default: \code{update.gds = FALSE}.
# @param depth.tibble (optional, logical) Returns the depth info in a tibble instead of list
# with total and mean coverage info.
# Used internally.
# Default:\code{depth.tibble = FALSE}.
#' @inheritParams radiator_common_arguments

# @keywords internal
#' @export
extract_coverage <- function(
    gds,
    individuals = TRUE,
    markers = TRUE,
    dp = TRUE,
    ad = TRUE,
    coverage.stats = c("sum", "mean", "median", "iqr"),
    subsample.info = 1,
    verbose = TRUE
    ) {

  if (verbose) cli::cli_progress_step("Coverage ...")
  i.info <- i.stats <- m.info <- m.stats <- NULL

  # internal function
  coverage_stats <- function(
    gds,
    coverage.stats = c("sum", "mean", "median", "iqr"),
    dp = TRUE,
    ad = TRUE,
    individuals = TRUE,
    markers = TRUE
  ) {

    coverage.stats <- match.arg(
      arg = coverage.stats,
      choices = c("sum", "mean", "median", "iqr"),
      several.ok = TRUE
    )
    dp.i <- dp.m <- ad.m <- ad.i <- NULL


    if (dp) {
      coverage.stats.l <- as.list(coverage.stats)
      names(coverage.stats.l) <- stringi::stri_replace_all_fixed(
        str = coverage.stats,
        pattern = c("sum", "mean", "median", "iqr"),
        replacement = c("COVERAGE_TOTAL", "COVERAGE_MEAN", "COVERAGE_MEDIAN", "COVERAGE_IQR"),
        vectorize_all = FALSE
      )

      if (markers) {
        dp_f_m <- function(gds, coverage.stats) {

          # Using switch instead was not optimal for additional options in the func...
          if (coverage.stats == "sum") rad_cov_stats <- function(x) round(sum(x, na.rm = TRUE))
          if (coverage.stats == "mean") rad_cov_stats <- function(x) round(mean(x, na.rm = TRUE))
          if (coverage.stats == "median") rad_cov_stats <- function(x) round(stats::median(x, na.rm = TRUE))
          # if (coverage.stats == "iqr") rad_cov_stats <- function(x) round(abs(diff(stats::quantile(x, probs = c(0.25, 0.75), na.rm = TRUE))))
          if (coverage.stats == "iqr") rad_cov_stats <- function(x) round(matrixStats::iqr(x, na.rm = TRUE)) # faster

          x <- SeqArray::seqApply(
            gdsfile = gds,
            var.name = "annotation/format/DP",
            FUN = rad_cov_stats,
            as.is = "integer",
            margin = "by.variant",
            parallel = TRUE
          )
        }

        dp.m <- purrr::map_dfc(.x = coverage.stats.l, .f = dp_f_m, gds = gds)
      }

      if (individuals) {
        # changing the margin doesnt work with seqarray, the GDS needs to be optimized by sample
        # this operation is very costly in both time and disk space.
        # faster to do matrix calculations by rows and sums
        # Note to myself: the huge speed gain by using other packages robustbse, Rfast, etc.
        # is not worth the unreliability of the results check your testing files...

        dp.temp <- SeqArray::seqGetData(
          gdsfile = gds,
          var.name = "annotation/format/DP"
        )

        dp_f_i <- function(coverage.stats, x) {
          if ("sum" %in% coverage.stats) x <- rowSums(x, na.rm = TRUE)
          if ("mean" %in% coverage.stats) x <- rowMeans(x, na.rm = TRUE)
          if ("median" %in% coverage.stats) x <- matrixStats::rowMedians(x, na.rm = TRUE)
          if ("iqr" %in% coverage.stats) x <- matrixStats::rowIQRs(x, na.rm = TRUE)
          x <- as.integer(round(x))
          return(x)
        }

        dp.i <- purrr::map_dfc(.x = coverage.stats.l, .f = dp_f_i, x = dp.temp)
        dp.temp <- NULL
      }
    }


    if (ad) {
      #temp object contains AD for REF and ALT
      ref <- SeqArray::seqGetData(
        gdsfile = gds,
        var.name = "annotation/format/AD"
      )$data


      # to extract the REF and ALT
      column.vec <- seq_len(length.out = dim(ref)[2])
      alt <- ref[, column.vec %% 2 == 0]
      alt[alt == 0] <- NA
      ref <- ref[, column.vec %% 2 == 1]
      ref[ref == 0] <- NA
      column.vec <- NULL

      ad_f <- function(coverage.stats, x, margin = c("markers", "individuals")) {

        margin <- match.arg(
          arg = margin,
          choices = c("markers", "individuals"),
          several.ok = FALSE
        )

        if (margin == "markers") {
          if ("sum" %in% coverage.stats) x <- colSums(x, na.rm = TRUE)
          if ("mean" %in% coverage.stats) x <- colMeans(x, na.rm = TRUE)
          if ("median" %in% coverage.stats) x <- matrixStats::colMedians(x, na.rm = TRUE)
          if ("iqr" %in% coverage.stats) x <- matrixStats::colIQRs(x, na.rm = TRUE)
          x <- as.integer(round(x))
          return(x)
        }
        if (margin == "individuals") {
          if ("sum" %in% coverage.stats) x <- rowSums(x, na.rm = TRUE)
          if ("mean" %in% coverage.stats) x <- rowMeans(x, na.rm = TRUE)
          if ("median" %in% coverage.stats) x <- matrixStats::rowMedians(x, na.rm = TRUE)
          if ("iqr" %in% coverage.stats) x <- matrixStats::rowIQRs(x, na.rm = TRUE)
          x <- as.integer(round(x))
          return(x)
        }
      }

      # for ref and alt
      coverage.stats.ref <- coverage.stats.alt <- as.list(coverage.stats)

      names(coverage.stats.ref) <- stringi::stri_replace_all_fixed(
        str = coverage.stats,
        pattern = c("sum", "mean", "median", "iqr"),
        replacement = c("REF_DEPTH_TOTAL", "REF_DEPTH_MEAN", "REF_DEPTH_MEDIAN", "REF_DEPTH_IQR"),
        vectorize_all = FALSE
      )
      names(coverage.stats.alt) <- stringi::stri_replace_all_fixed(
        str = coverage.stats,
        pattern = c("sum", "mean", "median", "iqr"),
        replacement = c("ALT_DEPTH_TOTAL", "ALT_DEPTH_MEAN", "ALT_DEPTH_MEDIAN", "ALT_DEPTH_IQR"),
        vectorize_all = FALSE
      )

      if (markers) {
        ad.m <- dplyr::bind_cols(
          purrr::map_dfc(.x = coverage.stats.ref, .f = ad_f, x = ref, margin = "markers"),
          purrr::map_dfc(.x = coverage.stats.alt, .f = ad_f, x = alt, margin = "markers")
        )
      }

      if (individuals) {
        ad.i <- dplyr::bind_cols(
          purrr::map_dfc(.x = coverage.stats.ref, .f = ad_f, x = ref, margin = "individuals"),
          purrr::map_dfc(.x = coverage.stats.alt, .f = ad_f, x = alt, margin = "individuals")
        )
      }
      ref <- alt <- NULL
    }

    cov.m <- dplyr::bind_cols(dp.m, ad.m)
    cov.i <- dplyr::bind_cols(dp.i, ad.i)

    cov.stats <- list(markers = cov.m, individuals = cov.i)

    return(cov.stats)
  } # END dp_stats


  c.s <- coverage_stats(
    gds = gds,
    coverage.stats = coverage.stats,
    dp = dp,
    ad = ad,
    individuals = individuals,
    markers = markers
  )

  # required for individuals and markers
  cov_tibble_stats <- function(have, tibble.group, data, subsample.info) {
    if (have %in% names(data)) {
      cov.tib <- tibble_stats(x = as.numeric(data[[have]]), group = tibble.group, subsample = subsample.info)
    } else {
      cov.tib <- NULL
    }
    return(cov.tib)
  }#End cov_tibble_stats



  if (individuals) {
    i.info %<>%
      dplyr::bind_cols(c.s$individuals)

    have <- names(c.s$individuals)
    want <- c("COVERAGE_TOTAL", "COVERAGE_MEAN", "COVERAGE_MEDIAN", "COVERAGE_IQR", "REF_DEPTH_TOTAL", "REF_DEPTH_MEAN", "REF_DEPTH_MEDIAN", "REF_DEPTH_IQR", "ALT_DEPTH_TOTAL", "ALT_DEPTH_MEAN", "ALT_DEPTH_MEDIAN", "ALT_DEPTH_IQR")
    tibble.group <- c("total coverage", "mean coverage", "median coverage", "iqr coverage", "total ref depth", "mean ref depth", "median ref depth", "iqr ref depth", "total alt depth", "mean alt depth", "median alt depth", "iqr alt depth")

    have <- purrr::keep(.x = have, .p = have %in% want)

    tibble.group <- stringi::stri_replace_all_fixed(
      str = have,
      pattern = want,
      replacement = tibble.group,
      vectorize_all = FALSE
    )

    i.stats %<>%
      dplyr::bind_rows(
        purrr::map2_dfr(
          .x = as.list(have),
          .y = as.list(tibble.group),
          .f = cov_tibble_stats,
          data = i.info,
          subsample.info = subsample.info
        )
      )
  }#End ind

  if (markers) {
    m.info %<>%
      dplyr::bind_cols(c.s$markers)

    have <- names(c.s$markers)
    want <- c("COVERAGE_TOTAL", "COVERAGE_MEAN", "COVERAGE_MEDIAN", "COVERAGE_IQR", "REF_DEPTH_TOTAL", "REF_DEPTH_MEAN", "REF_DEPTH_MEDIAN", "REF_DEPTH_IQR", "ALT_DEPTH_TOTAL", "ALT_DEPTH_MEAN", "ALT_DEPTH_MEDIAN", "ALT_DEPTH_IQR")
    tibble.group <- c("total coverage", "mean coverage", "median coverage", "iqr coverage", "total ref depth", "mean ref depth", "median ref depth", "iqr ref depth", "total alt depth", "mean alt depth", "median alt depth", "iqr alt depth")

    have <- purrr::keep(.x = have, .p = have %in% want)

    tibble.group <- stringi::stri_replace_all_fixed(
      str = have,
      pattern = want,
      replacement = tibble.group,
      vectorize_all = FALSE
    )


    m.stats %<>%
      dplyr::bind_rows(
        purrr::map2_dfr(
          .x = as.list(have),
          .y = as.list(tibble.group),
          .f = cov_tibble_stats,
          data = m.info,
          subsample.info = subsample.info
        )
      )
  }#End markers
  return(list(i.info = i.info, i.stats = i.stats, m.info = m.info, m.stats = m.stats))
}#End extract_coverage


# parse_gds_metadata -----------------------------------------------------------
#' @title parse_gds_metadata
#' @description function to parse the format field and tidy the results of VCF
#' @rdname parse_gds_metadata
#' @keywords internal
#' @export
parse_gds_metadata <- function(
    x,
    gds = NULL,
    verbose = TRUE,
    strip.rad = FALSE
) {

  # x <- parse.format.list
  # format.name <- x <- "DP"
  # format.name <- x <- "AD"
  # format.name <- x <- "GL"
  # format.name <- x <- "PL"
  # format.name <- x <- "HQ"
  # format.name <- x <- "GQ"
  # format.name <- x <- "GOF"
  # format.name <- x <- "NR"
  # format.name <- x <- "NV"
  # format.name <- x <- "RO" # not yet implemented
  # format.name <- x <- "QR" # not yet implemented
  # format.name <- x <- "AO" # not yet implemented
  # format.name <- x <- "QA" # not yet implemented
  res <- list()
  format.name <- as.list(x)
  names(format.name) <- x
  if (verbose) message("\nParsing and tidying genotypes metadata: ", stringi::stri_join(x, collapse = ", "))

  # sample and markers in long format...
  if (strip.rad) {
    i <- "ID_SEQ"
    m <- "M_SEQ"
  } else {
    i <- "INDIVIDUALS"
    m <- "VARIANT_ID"
  }

  # When subsampling is used, different numbers of markers ....
  sum.gds <- summary_gds(gds, verbose = FALSE)
  n.ind <- sum.gds$n.ind
  n.markers <- sum.gds$n.markers
  sum.gds <- NULL


  i <- radiator::extract_individuals_metadata(
    gds = gds,
    ind.field.select = i,
    whitelist = TRUE
  ) %>%
    dplyr::pull()

  if (length(i) != n.ind) {
    rlang::abort("Problem with GDS sync, contact author and mention: parsing problem")
  }

  m.temp <- radiator::extract_markers_metadata(
    gds = gds,
    markers.meta.select = m,
    whitelist = TRUE
  ) %>%
    dplyr::pull()

  if (length(m.temp) != n.markers) {
    v.select <- radiator::extract_markers_metadata(gds = gds, radiator.node = FALSE) %>% dplyr::pull("VARIANT_ID")
    m.temp <- radiator::extract_markers_metadata(gds = gds, radiator.node = TRUE, markers.meta.select = c("VARIANT_ID","M_SEQ"), whitelist = TRUE) %>%
      dplyr::filter(VARIANT_ID %in% v.select) %>%
      dplyr::pull(m)
  }
  m <- m.temp
  n.markers <- length(m)
  m.temp <- NULL

  if (strip.rad) {
    res$info <- tibble::tibble(
      ID_SEQ = rep(i, n.markers),
      M_SEQ = sort(rep(m, n.ind)) # faster than dplyr::arrange(M_SEQ, ID_SEQ)
    )
  } else {
    res$info <- tibble::tibble(
      INDIVIDUALS = rep(i, n.markers),
      VARIANT_ID = sort(rep(m, n.ind)) # faster than dplyr::arrange(M_SEQ, ID_SEQ)
    ) %>%
      dplyr::mutate(INDIVIDUALS = factor(x = INDIVIDUALS, levels = i))
  }
  m <- i <- NULL

  gds_metadata <- function(gds, format.name, verbose = TRUE) {
    # Allele Depth
    if ("AD" %in% format.name) {
      # NOTE TO MYSELF: THIS SHOULD BE DONE FOR BIALLELIC DATA ONLY...
      if (verbose) message("AD column: splitting into ALLELE_REF_DEPTH and ALLELE_ALT_DEPTH")
      temp.ad <- SeqArray::seqGetData(
        gdsfile = gds,
        var.name = "annotation/format/AD"
      )$data

      column.vec <- seq_len(length.out = dim(temp.ad)[2])

      make_tibble <- function(...) tibble::tibble(...)
      safe_make_tibble <- purrr::safely(.f = make_tibble)
      safe.ad <- safe_make_tibble(
        ALLELE_REF_DEPTH = temp.ad[, column.vec %% 2 == 1] %>%
          as.vector(.) %>% replace_by_na(data = ., what = 0L),
        ALLELE_ALT_DEPTH = temp.ad[, column.vec %% 2 == 0] %>%
          as.vector(.) %>% replace_by_na(data = ., what = 0L)
      )

      if (is.null(safe.ad$error)) {
        AD <- safe.ad$result
      } else {
        message("\n\nCaution: something is wrong with the dataset AD")
        message("Likely missing information for some sample(s)")
      }
      make_tibble <- safe_make_tibble <- safe.ad <- temp.ad <- column.vec <- NULL
      return(AD)
    }# End AD

    # Read depth
    if ("DP" %in% format.name) {
      if (verbose) message("DP column: cleaning and renaming to READ_DEPTH")
      DP <- tibble::tibble(
        READ_DEPTH = SeqArray::seqGetData(
          gdsfile = gds,
          var.name = "annotation/format/DP"
        ) %>%
          as.vector(.)
      )
      return(DP)
    } # End DP

    # Cleaning HQ: Haplotype quality as phred score
    if ("HQ" %in% format.name) {
      HQ <- tibble::tibble(HQ = SeqArray::seqGetData(
        gdsfile = gds,
        var.name = "annotation/format/HQ")$data %>% as.vector(.))

      # test <- res$HQ

      # check HQ and new stacks version with no HQ
      all.missing <- nrow(HQ)
      if (all.missing != 0) {
        if (verbose) message("HQ column: Haplotype Quality")
      } else {
        message("HQ values are all missing: removing column")
        HQ <- NULL
      }
      return(HQ)
    } # End HQ

    # Cleaning GQ: Genotype quality as phred score
    if ("GQ" %in% format.name) {
      if (verbose) message("GQ column: Genotype Quality")
      GQ <- tibble::tibble(
        GQ = SeqArray::seqGetData(
          gdsfile = gds,
          var.name = "annotation/format/GQ"
        ) %>%
          as.vector(.)
      )
      return(GQ)
    } # End GQ

    # GL cleaning
    if ("GL" %in% format.name) {
      if (verbose) message("GL column: cleaning Genotype Likelihood column")
      gl <- unique(SeqArray::seqGetData(gdsfile = gds,
                                        var.name = "annotation/format/GL")$length)
      if (gl > 0) {
        GL <- SeqArray::seqGetData(
          gdsfile = gds,
          var.name = "annotation/format/GL")$data

        column.vec <- seq_len(length.out = dim(GL)[2])
        GL <- tibble::tibble(
          GL_HOM_REF = GL[, column.vec %% 3 == 1] %>%
            as.vector(.),
          GL_HET = GL[, column.vec %% 3 == 2] %>%
            as.vector(.),
          GL_HOM_ALT = GL[, column.vec %% 3 == 0] %>%
            as.vector(.)
        )
        GL[GL == "NaN"] <- NA
        gl <- column.vec <- NULL
      }
      return(GL)
    } # End GL

    # Cleaning GOF: Goodness of fit value
    if ("GOF" %in% format.name) {
      if (verbose) message("GOF column: Goodness of fit value")
      GOF <- tibble::tibble(
        GOF = SeqArray::seqGetData(
          gdsfile = gds,
          var.name = "annotation/format/GOF"
        ) %>%
          as.vector(.)
      )
      return(GOF)
    } # End GOF

    # Cleaning NR: Number of reads covering variant location in this sample
    if ("NR" %in% format.name) {
      if (verbose) message("NR column: splitting column into the number of variant")
      NR <- tibble::tibble(
        NR = SeqArray::seqGetData(
          gdsfile = gds,
          var.name = "annotation/format/NR"
        ) %>%
          as.vector(.)
      )
      return(NR)
    }#End cleaning NR column

    # Cleaning NV: Number of reads containing variant in this sample
    if ("NV" %in% format.name) {
      if (verbose) message("NV column: splitting column into the number of variant")
      NV <- tibble::tibble(
        NV = SeqArray::seqGetData(
          gdsfile = gds,
          var.name = "annotation/format/NV"
        ) %>%
          as.vector(.)
      )
      return(NV)
    }#End cleaning NV column

    # CATG from ipyrad...
    if ("CATG" %in% format.name) {
      if (verbose) message("CATG columns: cleaning and renaming C_DEPTH, A_DEPTH, T_DEPTH, G_DEPTH")
      CATG <- SeqArray::seqGetData(gdsfile = gds,var.name = "annotation/format/CATG")$data

      column.vec <- seq_len(length.out = dim(CATG)[2])

      CATG <- tibble::tibble(
        C_DEPTH =  CATG[, column.vec %% 4 == 1] %>%
          as.integer(.),
        A_DEPTH =  CATG[, column.vec %% 4 == 2] %>%
          as.integer(.),
        T_DEPTH =  CATG[, column.vec %% 4 == 3] %>%
          as.integer(.),
        G_DEPTH =  CATG[, column.vec %% 4 == 0] %>%
          as.integer(.)
      )
      column.vec <- NULL
      return(CATG)
    } # End CATG
  }# END gds_metadata

  res$metadata <- purrr::map_dfc(.x = format.name, .f = gds_metadata, gds = gds)

  # Flatten the tibble for results
  res %<>% purrr::flatten_dfc(.)
  return(res)
}#End parse_gds_metadata


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
      if (verbose) message("synchronizing GDS with samples provided")
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
  radiator_packages_dep("SeqArray", cran = FALSE, bioc = TRUE)

  data.type <- radiator::detect_genomic_format(gds)

  if (!data.type %in% c("SeqVarGDSClass", "gds.file")) {
    rlang::abort("Input not supported for this function: read function documentation")
  }
  if (data.type == "gds.file") {
    gds <- radiator::read_rad(gds)
    data.type <- "SeqVarGDSClass"
  }
  i <- extract_individuals_metadata(gds = gds, ind.field.select = "FILTERS", blacklist = FALSE) %>%
    dplyr::count(FILTERS, FILTERS) #%>% readr::write_tsv(x = i, file = "filters.individuals.tsv")
  i %<>% dplyr::filter(FILTERS != "whitelist")
  message("Number of filters for individuals: ", nrow(i))
  if (nrow(i) > 0) {
    message("Filter(s): ")
    message(stringi::stri_join(i$FILTERS, collapse = "\n"))
  }
  m <- extract_markers_metadata(gds = gds, markers.meta.select = "FILTERS", blacklist = FALSE) %>%
    dplyr::count(FILTERS, FILTERS) #%>% readr::write_tsv(x = i, file = "filters.individuals.tsv")
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
#' @param filter.ma (logical, optional)
#' Default: \code{filter.ma = FALSE}.
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
    filter.ma = FALSE,
    filter.mean.coverage = FALSE,
    filter.genotyping = FALSE,
    filter.snp.position.read = FALSE,
    filter.snp.number = FALSE,
    filter.short.ld = FALSE,
    filter.long.ld = FALSE,
    filter.hwe = FALSE,
    filter.whitelist = FALSE
) {
  radiator_packages_dep(package = "SeqArray", cran = FALSE, bioc = TRUE)
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

  if (filter.ma) {
    reset.m$filter.ma <- "filter.ma"
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
  if (verbose) message("\nData summary: ")
  check <- SeqArray::seqGetFilter(gds)
  n.ind <- length(check$sample.sel[check$sample.sel])
  n.markers <- length(check$variant.sel[check$variant.sel])

  if (verbose) message("    Number of individuals: ", n.ind)
  if (verbose) message("    Number of markers: ", n.markers)

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


# generate statistics--------------------------------------------------------
#' @title generate_stats
#' @description Generate individuals and markers statistics
#' @rdname generate_stats
#' @keywords internal
#' @export
generate_stats <- function(
    gds,
    individuals = TRUE,
    markers = TRUE,
    missing = TRUE,
    heterozygosity = TRUE,
    coverage = TRUE,
    allele.coverage = TRUE,
    mac = TRUE,
    snp.position.read = TRUE,
    snp.per.locus = TRUE,
    subsample = NULL,
    exhaustive = TRUE,
    force.stats = TRUE,
    path.folder = NULL,
    plot = TRUE,
    digits = 6,
    file.date = NULL,
    parallel.core = parallel::detectCores() - 1,
    verbose = TRUE
) {

  ## TEST
  # individuals = TRUE
  # markers = FALSE
  # missing = TRUE
  # heterozygosity = TRUE
  # coverage = TRUE
  # allele.coverage = TRUE
  # mac = TRUE
  # snp.position.read = TRUE
  # snp.per.locus = TRUE
  # subsample = variant.select
  # exhaustive = FALSE
  # force.stats = TRUE
  # # path.folder = NULL
  # plot = TRUE
  # digits = 6
  # # file.date = NULL
  # parallel.core = parallel::detectCores() - 1
  # verbose = TRUE

  if (verbose) message("Generating statistics")

  if (is.null(path.folder)) path.folder <- getwd()
  if (is.null(file.date)) file.date <- format(Sys.time(), "%Y%m%d@%H%M")

  # Defaults -------------------------------------------------------------------
  m.stats <- i.stats <- i.info <- m.info <- NULL

  if (allele.coverage) mac <- TRUE # by default the mac info is needed...
  if (!markers) {
    snp.position.read <- snp.per.locus <- mac <- FALSE
  }

  # required info --------------------------------------------------------------
  data.source <- extract_data_source(gds)
  if (individuals) {
    i.info <- extract_individuals_metadata(gds, ind.field.select = c("INDIVIDUALS", "STRATA"), whitelist = TRUE)
  }

  # need this info even if markers = FALSE
  m.info <- extract_markers_metadata(gds = gds, whitelist = TRUE)
  n.markers.tot <- nrow(m.info)
  if (!rlang::has_name(m.info, "COL")) snp.position.read <- FALSE

  subsample.info <- NULL # default for reporting in the tibble

  # snp.per.locus and snp.position.read are very fast with huge dataset, no need for subsampling...

  # snp.per.locus --------------------------------------------------------------
  if (snp.per.locus) {
    if (verbose) cli::cli_progress_step("SNPs per locus")
    if (!rlang::has_name(m.info, "SNP_PER_LOCUS") || force.stats) {
      biallelic <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(
        node = gds,
        path = "radiator/biallelic", silent = TRUE))
      if (!is.logical(biallelic)) biallelic <- FALSE
      if (biallelic) {
        m.info %<>%
          dplyr::group_by(LOCUS) %>%
          dplyr::mutate(SNP_PER_LOCUS = n()) %>%
          dplyr::ungroup(.)
      } else {
        if (verbose) message("With haplotype vcf, number of SNP/locus is counted based on the REF allele")
        if (verbose) message("    this stat is good only if nuc length is equal between REF and ALT haplotypes")
        if (verbose) message("    stacks haplotype vcf: ok")
        if (verbose) message("    dDocent/freebayes haplotype vcf: be careful")

        m.info %<>%
          dplyr::mutate(
            SNP_PER_LOCUS = stringi::stri_length(
              str = SeqArray::seqGetData(gdsfile = gds, var.name = "$ref"))
          )
      }
    }
    m.stats %<>%
      dplyr::bind_rows(
        tibble_stats(x = m.info$SNP_PER_LOCUS, group = "SNPs per locus", subsample = subsample.info)
      )
  } # snp.per.locus

  # snp.position.read -----------------------------------------------------------
  if (snp.position.read) {
    if (verbose) cli::cli_progress_step("SNP position on the read")
    m.stats %<>%
      dplyr::bind_rows(
        tibble_stats(
          x = dplyr::distinct(m.info, COL) %$% COL,
          group = "SNP position on read",
          subsample = subsample.info
        )
      )
  } # snp.position.read

  # missing --------------------------------------------------------------------
  if (missing) {
    if (verbose) cli::cli_progress_step("Missing genotypes")

    if (individuals) {
      # info
      i.info %<>%
        dplyr::mutate(
          MISSING_PROP = round(
            SeqArray::seqMissing(
              gdsfile = gds,
              per.variant = FALSE,
              parallel = parallel.core
            )
            , digits
          )
        )
      # stats
      i.stats %<>%
        dplyr::bind_rows(
          tibble_stats(
            x = i.info$MISSING_PROP,
            group = "missing genotypes",
            subsample = subsample.info
          )
        )
    } # End ind
    if (markers) {
      if (!rlang::has_name(m.info, "MISSING_PROP") || force.stats) {
        m.info$MISSING_PROP <- SeqArray::seqMissing(
          gdsfile = gds,
          per.variant = TRUE,
          parallel = parallel.core
        )
      }
      m.stats %<>%
        dplyr::bind_rows(
          tibble_stats(
            x = m.info$MISSING_PROP,
            group = "missing genotypes",
            subsample = subsample.info
          )
        )
    }#End markers
  }# End missing

  # heterozygosity -------------------------------------------------------------
  if (heterozygosity) {
    if (verbose) cli::cli_progress_step("Heterozygosity")

    if (individuals) {
      # info
      i.info %<>%
        dplyr::mutate(HETEROZYGOSITY = round(individual_het(gds), digits))

      i.info$HETEROZYGOSITY[is.na(i.info$HETEROZYGOSITY)] <- 0

      # stats
      i.stats %<>%
        dplyr::bind_rows(
          tibble_stats(
            x = i.info$HETEROZYGOSITY,
            group = "heterozygosity",
            subsample = subsample.info
          )
        )
    }#end ind
    if (markers) {
      if (!rlang::has_name(m.info, "HET_OBS") || force.stats) {
        m.info %<>%
          dplyr::mutate(
            HET_OBS = round(markers_het(gds), 6),
            FIS = round(markers_fis(gds), 6)
          )
      }
      m.stats %<>%
        dplyr::bind_rows(
          tibble_stats(
            x = m.info$HET_OBS,
            group = "observed heterozygosity",
            subsample = subsample.info
          )
        ) %>%
        dplyr::bind_rows(
          tibble_stats(
            x = m.info$FIS,
            group = "inbreeding coefficient (Fis)",
            subsample = subsample.info
          )
        )

    }#End markers
  }#End het

  # MAC ------------------------------------------------------------------------
  if (mac) {
    if (verbose) cli::cli_progress_step("MAC")

    if (!rlang::has_name(m.info, "MAC_GLOBAL") || force.stats) {

      ac <- SeqArray::seqAlleleCount(
        gdsfile = gds,
        ref.allele = NULL,
        parallel = parallel.core)

      # check if more than 2 alternate alleles...
      m.info %<>%
        dplyr::mutate(
          N_ALLELES = purrr::map_int(.x = ac, .f = length),
          REF_COUNT = purrr::map_int(.x = ac, .f = max),
          A_SUM = purrr::map_int(.x = ac, .f = sum),
          MAC_GLOBAL = A_SUM - REF_COUNT,
          MAF_GLOBAL = MAC_GLOBAL / A_SUM,
          A_SUM = NULL
        )
      ac <- NULL
      n.al.max <- max(m.info$N_ALLELES, na.rm = TRUE)

      # when more than 2 alternate alleles...
      if (n.al.max > 2) {
        m.info %<>% dplyr::mutate(ALT_CHECK = "ALT")
        message("\n\nCaution: More than 2 alleles detected in the dataset")
        wanted.info <- c("VARIANT_ID", "CHROM", "LOCUS", "POS", "MARKERS", "REF", "ALT", "N_ALLELES")
        m.info %>%
          dplyr::filter(N_ALLELES > 2) %>%
          dplyr::select(dplyr::any_of(wanted.info)) %>%
          write_rad(
            data = .,
            path = path.folder,
            filename = "markers_number_alleles_problem.tsv",
            tsv = TRUE,
            write.message = "standard",
            verbose = TRUE
          )
      } else {
        m.info %<>%
          dplyr::mutate(ALT_CHECK = dplyr::if_else(MAC_GLOBAL <= REF_COUNT, "ALT", "REF"))
      }

      if (!allele.coverage) m.info %<>% dplyr::select(-ALT_CHECK)

    }
    m.stats %<>%
      dplyr::bind_rows(
        tibble_stats(
          x = m.info$MAC_GLOBAL,
          group = "MAC",
          subsample = subsample.info
        )
      )
  } #mac

  # subsampling ----------------------------------------------------------------
  if (!is.null(subsample)) {
    SeqArray::seqSetFilter(object = gds,
                           variant.id = subsample,
                           action = "push+set",
                           verbose = FALSE)
    m.info %<>% dplyr::filter(VARIANT_ID %in% subsample)
    n.markers <- length(m.info$VARIANT_ID)
    subsample.info <- round(n.markers / n.markers.tot, 2)
  }# End subsampling


  # # checking for REAL ALTERNATE ALLELE using REF and ALT READ DEPTH
  # if (allele.coverage) {
  #
  #   if (!rlang::has_name(m.info, "REF_DEPTH_TOTAL") || force.stats) {
  #
  #     ad.info <- gdsfmt::index.gdsn(gds, "annotation/format/AD", silent = TRUE)
  #     if (!is.null(ad.info)) {
  #       ad <- SeqArray::seqGetData(gds, "annotation/format/AD") %$% data
  #
  #       m.info %<>%
  #         dplyr::bind_cols(
  #           #REF and ALT Total read depth
  #           colSums(x = ad, na.rm = TRUE, dims = 1L) %>%
  #             unlist(.) %>%
  #             matrix(
  #               data = .,
  #               nrow = n.markers, ncol = 2, byrow = TRUE,
  #               dimnames = list(rownames = m.info$VARIANT_ID,
  #                               colnames = c("RT", "AT"))) %>%
  #             tibble::as_tibble(.),
  #           #REF and ALT Mean read depth
  #           colMeans(x = ad, na.rm = TRUE, dims = 1L) %>%
  #             unlist(.) %>%
  #             matrix(
  #               data = .,
  #               nrow = n.markers, ncol = 2, byrow = TRUE,
  #               dimnames = list(rownames = m.info$VARIANT_ID,
  #                               colnames = c("RM", "AM"))) %>%
  #             tibble::as_tibble(.)
  #         ) %>%
  #         dplyr::mutate(
  #           REF_DEPTH_TOTAL = dplyr::if_else(ALT_CHECK == "ALT", RT, AT),
  #           ALT_DEPTH_TOTAL = dplyr::if_else(ALT_CHECK == "ALT", AT, RT),
  #           RT = NULL, AT = NULL,
  #           REF_DEPTH_MEAN = dplyr::if_else(ALT_CHECK == "ALT", RM, AM),
  #           ALT_DEPTH_MEAN = dplyr::if_else(ALT_CHECK == "ALT", AM, RM),
  #           ALT_CHECK = NULL, RM = NULL, AM = NULL
  #         )
  #       ad <- ad.info <- NULL
  #     }
  #   }
  #
  #   if (rlang::has_name(m.info, "REF_DEPTH_TOTAL")) {
  #     stats.rt <- tibble_stats(
  #       x = m.info$REF_DEPTH_TOTAL,
  #       group = "total ref depth")
  #     stats.at <- tibble_stats(
  #       x = m.info$ALT_DEPTH_TOTAL,
  #       group = "total alt depth")
  #   }
  #
  #   if (rlang::has_name(m.info, "REF_DEPTH_MEAN")) {
  #     stats.rm <- tibble_stats(
  #       x = m.info$REF_DEPTH_MEAN,
  #       group = "mean ref depth")
  #     stats.am <- tibble_stats(
  #       x = m.info$ALT_DEPTH_MEAN,
  #       group = "mean alt depth")
  #   }
  #
  #
  # }
  # coverage -------------------------------------------------------------------
  # Confirm coverage information is available
  if (coverage || allele.coverage) {
    got.coverage <- check_coverage(
      gds = gds,
      stacks.haplo.check = TRUE,
      genotypes.metadata.check = TRUE,
      dart.check = TRUE
    )
    if (!"DP" %in% got.coverage && !"READ_DEPTH" %in% got.coverage) coverage <- FALSE
    if (!"AD" %in% got.coverage) allele.coverage <- FALSE
    if (!exhaustive) allele.coverage <- FALSE
    got.coverage <- NULL
  }

  if (coverage || allele.coverage) {
    if (verbose) cli::cli_progress_step("Coverage ...")

    # internal function
    coverage_stats <- function(
    gds,
    coverage.stats = c("sum", "mean", "median", "iqr"),
    dp = TRUE,
    ad = TRUE,
    individuals = TRUE,
    markers = TRUE
    ) {

      coverage.stats <- match.arg(
        arg = coverage.stats,
        choices = c("sum", "mean", "median", "iqr"),
        several.ok = TRUE
      )
      dp.i <- dp.m <- ad.m <- ad.i <- NULL


      if (dp) {
        coverage.stats.l <- as.list(coverage.stats)
        names(coverage.stats.l) <- stringi::stri_replace_all_fixed(
          str = coverage.stats,
          pattern = c("sum", "mean", "median", "iqr"),
          replacement = c("COVERAGE_TOTAL", "COVERAGE_MEAN", "COVERAGE_MEDIAN", "COVERAGE_IQR"),
          vectorize_all = FALSE
        )
        data.source <- radiator::extract_data_source(gds = gds)

        if ("dart" %in% data.source) {
          dart.data <- radiator::extract_genotypes_metadata(
            gds = gds,
            genotypes.meta.select = c("M_SEQ", "ID_SEQ", "READ_DEPTH"),
            whitelist = TRUE
          ) %>%
            radiator::rad_wide(
              x = .,
              formula = "ID_SEQ ~ M_SEQ",
              values_from = "READ_DEPTH"
            ) %>%
            dplyr::select(-ID_SEQ)
          colnames(dart.data) <- NULL
          dart.data <- as.matrix(dart.data)
        } else {
          dart.data <- NULL
        }

        if (markers) {
          dp_f_m <- function(gds, coverage.stats, dart.data) {
            # Using switch instead was not optimal for additional options in the func...
            if (coverage.stats == "sum") rad_cov_stats <- function(x) round(sum(x, na.rm = TRUE))
            if (coverage.stats == "mean") rad_cov_stats <- function(x) round(mean(x, na.rm = TRUE))
            if (coverage.stats == "median") rad_cov_stats <- function(x) round(stats::median(x, na.rm = TRUE))
            # if (coverage.stats == "iqr") rad_cov_stats <- function(x) round(abs(diff(stats::quantile(x, probs = c(0.25, 0.75), na.rm = TRUE))))
            if (coverage.stats == "iqr") rad_cov_stats <- function(x) round(matrixStats::iqr(x, na.rm = TRUE)) # faster

            if (!is.null(dart.data)) {
              x <- as.integer(apply(X = dart.data, MARGIN = 2, FUN = rad_cov_stats))
              dart.data <- NULL
            } else {
              x <- SeqArray::seqApply(
                gdsfile = gds,
                var.name = "annotation/format/DP",
                FUN = rad_cov_stats,
                as.is = "integer",
                margin = "by.variant",
                parallel = TRUE
              )
            }
            return(x)
          }

          dp.m <- purrr::map_dfc(.x = coverage.stats.l, .f = dp_f_m, gds = gds, dart.data = dart.data)
        }

        if (individuals) {
          # changing the margin doesnt work with seqarray, the GDS needs to be optimized by sample
          # this operation is very costly in both time and disk space.
          # faster to do matrix calculations by rows and sums
          # Note to myself: the huge speed gain by using other packages robustbse, Rfast, etc.
          # is not worth the unreliability of the results check your testing files...

          if ("dart" %in% data.source) {
            dp.temp <- dart.data
            dart.data <- NULL
          } else {
            dp.temp <- SeqArray::seqGetData(
              gdsfile = gds,
              var.name = "annotation/format/DP"
            )
          }


          dp_f_i <- function(coverage.stats, x) {
            if ("sum" %in% coverage.stats) x <- rowSums(x, na.rm = TRUE)
            if ("mean" %in% coverage.stats) x <- rowMeans(x, na.rm = TRUE)
            if ("median" %in% coverage.stats) x <- matrixStats::rowMedians(x, na.rm = TRUE)
            if ("iqr" %in% coverage.stats) x <- matrixStats::rowIQRs(x, na.rm = TRUE)
            x <- as.integer(round(x))
            return(x)
          }

          dp.i <- purrr::map_dfc(.x = coverage.stats.l, .f = dp_f_i, x = dp.temp)
          dp.temp <- NULL
        }
      }


      if (ad) {
        #temp object contains AD for REF and ALT
        ref <- SeqArray::seqGetData(
          gdsfile = gds,
          var.name = "annotation/format/AD"
        )$data


        # to extract the REF and ALT
        column.vec <- seq_len(length.out = dim(ref)[2])
        alt <- ref[, column.vec %% 2 == 0]
        alt[alt == 0] <- NA
        ref <- ref[, column.vec %% 2 == 1]
        ref[ref == 0] <- NA
        column.vec <- NULL

        ad_f <- function(coverage.stats, x, margin = c("markers", "individuals")) {

          margin <- match.arg(
            arg = margin,
            choices = c("markers", "individuals"),
            several.ok = FALSE
          )

          if (margin == "markers") {
            if ("sum" %in% coverage.stats) x <- colSums(x, na.rm = TRUE)
            if ("mean" %in% coverage.stats) x <- colMeans(x, na.rm = TRUE)
            if ("median" %in% coverage.stats) x <- matrixStats::colMedians(x, na.rm = TRUE)
            if ("iqr" %in% coverage.stats) x <- matrixStats::colIQRs(x, na.rm = TRUE)
            x <- as.integer(round(x))
            return(x)
          }
          if (margin == "individuals") {
            if ("sum" %in% coverage.stats) x <- rowSums(x, na.rm = TRUE)
            if ("mean" %in% coverage.stats) x <- rowMeans(x, na.rm = TRUE)
            if ("median" %in% coverage.stats) x <- matrixStats::rowMedians(x, na.rm = TRUE)
            if ("iqr" %in% coverage.stats) x <- matrixStats::rowIQRs(x, na.rm = TRUE)
            x <- as.integer(round(x))
            return(x)
          }
        }

        # for ref and alt
        coverage.stats.ref <- coverage.stats.alt <- as.list(coverage.stats)

        names(coverage.stats.ref) <- stringi::stri_replace_all_fixed(
          str = coverage.stats,
          pattern = c("sum", "mean", "median", "iqr"),
          replacement = c("REF_DEPTH_TOTAL", "REF_DEPTH_MEAN", "REF_DEPTH_MEDIAN", "REF_DEPTH_IQR"),
          vectorize_all = FALSE
        )
        names(coverage.stats.alt) <- stringi::stri_replace_all_fixed(
          str = coverage.stats,
          pattern = c("sum", "mean", "median", "iqr"),
          replacement = c("ALT_DEPTH_TOTAL", "ALT_DEPTH_MEAN", "ALT_DEPTH_MEDIAN", "ALT_DEPTH_IQR"),
          vectorize_all = FALSE
        )

        if (markers) {
          ad.m <- dplyr::bind_cols(
            purrr::map_dfc(.x = coverage.stats.ref, .f = ad_f, x = ref, margin = "markers"),
            purrr::map_dfc(.x = coverage.stats.alt, .f = ad_f, x = alt, margin = "markers")
          )
        }

        if (individuals) {
          ad.i <- dplyr::bind_cols(
            purrr::map_dfc(.x = coverage.stats.ref, .f = ad_f, x = ref, margin = "individuals"),
            purrr::map_dfc(.x = coverage.stats.alt, .f = ad_f, x = alt, margin = "individuals")
          )
        }
        ref <- alt <- NULL
      }

      cov.m <- dplyr::bind_cols(dp.m, ad.m)
      cov.i <- dplyr::bind_cols(dp.i, ad.i)

      cov.stats <- list(markers = cov.m, individuals = cov.i)

      return(cov.stats)
    } # END dp_stats

    c.s <- coverage_stats(
      gds = gds,
      coverage.stats = c("sum", "mean", "median", "iqr"),
      dp = coverage,
      ad = allele.coverage,
      individuals = individuals,
      markers = markers
    )

    # required for individuals and markers
    cov_tibble_stats <- function(have, tibble.group, data, subsample.info) {
      if (have %in% names(data)) {
        cov.tib <- tibble_stats(x = as.numeric(data[[have]]), group = tibble.group, subsample = subsample.info)
      } else {
        cov.tib <- NULL
      }
      return(cov.tib)
    }#End cov_tibble_stats



    if (individuals) {
      i.info %<>%
        dplyr::bind_cols(c.s$individuals)

      have <- names(c.s$individuals)
      want <- c("COVERAGE_TOTAL", "COVERAGE_MEAN", "COVERAGE_MEDIAN", "COVERAGE_IQR", "REF_DEPTH_TOTAL", "REF_DEPTH_MEAN", "REF_DEPTH_MEDIAN", "REF_DEPTH_IQR", "ALT_DEPTH_TOTAL", "ALT_DEPTH_MEAN", "ALT_DEPTH_MEDIAN", "ALT_DEPTH_IQR")
      tibble.group <- c("total coverage", "mean coverage", "median coverage", "iqr coverage", "total ref depth", "mean ref depth", "median ref depth", "iqr ref depth", "total alt depth", "mean alt depth", "median alt depth", "iqr alt depth")

      have <- purrr::keep(.x = have, .p = have %in% want)

      tibble.group <- stringi::stri_replace_all_fixed(
        str = have,
        pattern = want,
        replacement = tibble.group,
        vectorize_all = FALSE
      )

      i.stats %<>%
        dplyr::bind_rows(
          purrr::map2_dfr(
            .x = as.list(have),
            .y = as.list(tibble.group),
            .f = cov_tibble_stats,
            data = i.info,
            subsample.info = subsample.info
          )
        )
    }#End ind

    if (markers) {
      m.info %<>%
        dplyr::bind_cols(c.s$markers)

      have <- names(c.s$markers)
      want <- c("COVERAGE_TOTAL", "COVERAGE_MEAN", "COVERAGE_MEDIAN", "COVERAGE_IQR", "REF_DEPTH_TOTAL", "REF_DEPTH_MEAN", "REF_DEPTH_MEDIAN", "REF_DEPTH_IQR", "ALT_DEPTH_TOTAL", "ALT_DEPTH_MEAN", "ALT_DEPTH_MEDIAN", "ALT_DEPTH_IQR")
      tibble.group <- c("total coverage", "mean coverage", "median coverage", "iqr coverage", "total ref depth", "mean ref depth", "median ref depth", "iqr ref depth", "total alt depth", "mean alt depth", "median alt depth", "iqr alt depth")

      have <- purrr::keep(.x = have, .p = have %in% want)

      tibble.group <- stringi::stri_replace_all_fixed(
        str = have,
        pattern = want,
        replacement = tibble.group,
        vectorize_all = FALSE
      )


      m.stats %<>%
        dplyr::bind_rows(
          purrr::map2_dfr(
            .x = as.list(have),
            .y = as.list(tibble.group),
            .f = cov_tibble_stats,
            data = m.info,
            subsample.info = subsample.info
          )
        )
    }#End markers
  } # End coverage

  # NOTE TO MYSELF need to work on including that one with DArT files...--------
  #DArT 1 row and 2rows --------------------------------------------------------
  # if ("dart" %in% data.source && any(c("2rows", "1row") %in% data.source)) {
  #   depth <- extract_markers_metadata(
  #     gds,
  #     markers.meta.select = c("AVG_COUNT_REF", "AVG_COUNT_SNP"),
  #     whitelist = TRUE
  #   )
  #   markers <- ind <- FALSE
  #
  #   if (is.null(depth)) return(NULL)
  #
  #   coverage.info$markers.mean <- as.integer(
  #     round(depth$AVG_COUNT_REF + depth$AVG_COUNT_SNP, 0)
  #   )
  #   coverage.info$ref.mean <- as.integer(round(depth$AVG_COUNT_REF))
  #   coverage.info$alt.mean <- as.integer(round(depth$AVG_COUNT_SNP))
  #   depth <- NULL
  # }#End DART 1row and 2 rows
  #


  # work on the stats ----------------------------------------------------------
  markers.levels <- c(
    "missing genotypes",
    "MAC",  "observed heterozygosity",
    "inbreeding coefficient (Fis)",
    "SNP position on read",
    "SNPs per locus",
    "total coverage", "total ref depth", "total alt depth",
    "mean coverage", "mean ref depth", "mean alt depth",
    "median coverage", "median ref depth", "median alt depth",
    "iqr coverage", "iqr ref depth", "iqr alt depth"
  )
  m.stats$GROUP <- factor(x = m.stats$GROUP, levels = markers.levels, ordered = TRUE)
  m.stats$GROUP <- droplevels(x = m.stats$GROUP)

  id.levels <- c("missing genotypes", "heterozygosity", "total coverage",
                 "mean coverage", "median coverage", "iqr coverage",
                 "total ref depth", "mean ref depth", "median ref depth",
                 "iqr ref depth",
                 "total alt depth", "mean alt depth", "median alt depth",
                 "iqr alt depth")
  i.stats$GROUP <- factor(x = i.stats$GROUP, levels = id.levels, ordered = TRUE)
  i.stats$GROUP <- droplevels(x = i.stats$GROUP)


  # Generate plots -------------------------------------------------------------
  if (plot) {
    if (verbose) cli::cli_progress_step("Generating figures")
    i.fig.filename <- stringi::stri_join("individuals_qc_", file.date, ".pdf") # Figure
    m.fig.filename <- stringi::stri_join("markers_qc_", file.date, ".pdf") # Figure


    if (individuals) {
      # conditions for figures ---------------------------------------------------

      # info printed on top of figures
      corr.info <- stringi::stri_join("Correlations:\n")

      if (coverage && missing) {
        cm <- floor(stats::cor(i.info$COVERAGE_TOTAL, i.info$MISSING_PROP, use = "pairwise.complete.obs") * 100) / 100
        cmt <- stringi::stri_join("    total coverage & missing = ", cm)
        corr.info <- stringi::stri_join(corr.info, cmt)
      }
      if (coverage) {
        if (stats::sd(i.info$COVERAGE_MEAN) != 0) {
          cc <- ceiling(stats::cor(i.info$COVERAGE_TOTAL, i.info$COVERAGE_MEAN, use = "pairwise.complete.obs") * 100) / 100
        } else {
          cc <- "NA"
          message("Note: mean individual coverage SD = 0")
          message("correlation with total coverage is not possible")
        }
        cct <- stringi::stri_join("\n    total coverage & mean coverage = ", cc)
        corr.info <- stringi::stri_join(corr.info, cct)
      }
      if (coverage && heterozygosity) {
        ch <- ceiling(stats::cor(i.info$COVERAGE_TOTAL, i.info$HETEROZYGOSITY, use = "pairwise.complete.obs") * 100) / 100
        cht <- stringi::stri_join("\n    total coverage & heterozygosity = ", ch)
        corr.info <- stringi::stri_join(corr.info, cht)
      }
      if (missing && heterozygosity) {
        mh <- floor(stats::cor(i.info$HETEROZYGOSITY, i.info$MISSING_PROP, use = "pairwise.complete.obs") * 100) / 100
        mht <- stringi::stri_join("\n    missing & heterozygosity = ", mh)
        corr.info <- stringi::stri_join(corr.info, mht)
      }
      if (heterozygosity) {
        n.markers <- summary_gds(gds, verbose = FALSE) %$% n.markers
        n.markers.range <- ceiling((i.stats[[2,6]] - i.stats[[2,2]] ) * n.markers)
        n.markers.iqr <- ceiling(i.stats[[2,7]] * n.markers)

        corr.info <- stringi::stri_join(
          "n.markers = ", n.markers,
          "\nn. het markers in the bp range = ", n.markers.range,
          "\nn. het markers in the bp IQR = ", n.markers.iqr,
          "\n\n", corr.info
        )
      }
      if (missing) {
        missing.out <- stringi::stri_join("Missing genotypes outlier: ", i.stats[[1, 9]])
        corr.info <- stringi::stri_join(missing.out, "\n", corr.info)
      }

      cm <- cmt <- cc <- cct <- ch <- cht <- mh <- mht <- NULL


      if (!is.null(subsample)) {
        subtitle.ind.stats <- stringi::stri_join(
          "Markers subsampled: ", length(subsample), "\n\n", corr.info)
      } else {
        subtitle.ind.stats <- corr.info
      }

      i.fig <- boxplot_stats(
        data = i.stats,
        title = "Individual's QC stats",
        subtitle = subtitle.ind.stats,
        x.axis.title = NULL,
        y.axis.title = "Statistics",
        facet.columns = TRUE,
        scale.log = TRUE,
        scientific.format = FALSE,
        bp.filename = i.fig.filename,
        path.folder = path.folder)
      # if (verbose) message("File written: individuals qc plot")

    }#End ind

    if (markers) {
      # correlations info
      corr.info <- stringi::stri_join("Correlations:\n")

      if (rlang::has_name(m.info, "COVERAGE_TOTAL") && missing) {
        cm <- floor(stats::cor(m.info$COVERAGE_TOTAL, m.info$MISSING_PROP, use = "pairwise.complete.obs") * 100) / 100
        cmt <- stringi::stri_join("    total coverage & missing = ", cm)
        corr.info <- stringi::stri_join(corr.info, cmt)
      }
      if (rlang::has_name(m.info, "COVERAGE_TOTAL") && rlang::has_name(m.info, "COVERAGE_MEAN")) {
        cc <- ceiling(stats::cor(m.info$COVERAGE_TOTAL, m.info$COVERAGE_MEAN, use = "pairwise.complete.obs") * 100) / 100
        cct <- stringi::stri_join("\n    total coverage & mean coverage = ", cc)
        corr.info <- stringi::stri_join(corr.info, cct)
      }
      if (rlang::has_name(m.info, "COVERAGE_TOTAL") && heterozygosity) {
        ch <- ceiling(stats::cor(m.info$COVERAGE_TOTAL, m.info$HET_OBS, use = "pairwise.complete.obs") * 100) / 100
        cht <- stringi::stri_join("\n    total coverage & heterozygosity = ", ch)
        corr.info <- stringi::stri_join(corr.info, cht)
      }
      if (missing && heterozygosity) {
        mh <- floor(stats::cor(m.info$HET_OBS, m.info$MISSING_PROP, use = "pairwise.complete.obs") * 100) / 100
        mht <- stringi::stri_join("\n    missing & heterozygosity = ", mh)
        corr.info <- stringi::stri_join(corr.info, mht)
      }

      cm <- cc <- ch <- mh <- mht <- cmt <- cct <- cht <- NULL

      if (!is.null(subsample)) {
        subtitle.stats <- stringi::stri_join(
          "Markers subsampled: ", length(subsample), "\n\n", corr.info)
      } else {
        subtitle.stats <- corr.info
      }

      if (!heterozygosity && !rlang::has_name(m.info, "COVERAGE_TOTAL")) {
        subtitle.stats <- NULL
      }

      m.fig <- boxplot_stats(
        data = m.stats,
        title = "Marker's QC stats",
        subtitle = subtitle.stats,
        x.axis.title = NULL,
        scale.log = TRUE,
        scientific.format = FALSE,
        y.axis.title = "Statistics",
        facet.columns = TRUE,
        bp.filename = m.fig.filename,
        path.folder = path.folder)
    }#End markers
  }#End fig

  # Subsampling back to normal -------------------------------------------------
  if (!is.null(subsample)) SeqArray::seqSetFilter(gds, action = "pop", verbose = FALSE)

  # write files ----------------------------------------------------------------
  if (verbose) cli::cli_progress_step("Writing files")
  i.stats.f <- stringi::stri_join("individuals.qc.stats_", file.date, ".tsv")
  m.stats.f <- stringi::stri_join("markers.qc.stats_", file.date, ".tsv")
  i.stats.f.sum <- stringi::stri_join("individuals.qc.stats.summary_", file.date, ".tsv")# tibble summary stats
  m.stats.f.sum <- stringi::stri_join("markers.qc.stats.summary_", file.date, ".tsv")# tibble summary stats

  if (!markers) m.stats <- m.info <- m.fig <- NULL
  if (!individuals) i.stats <- i.info <- i.fig <- NULL
  if (!plot) m.fig <- i.fig <- NULL


  if (!is.null(m.stats)) {
    dplyr::mutate(.data = m.stats,
                  dplyr::across(.cols = where(is.numeric), .fns = round, digits = digits),
                  dplyr::across(.cols = where(is.numeric), .fns = format, scientific = FALSE),
    ) %>%
      readr::write_tsv(x = ., file = file.path(path.folder, m.stats.f.sum))
  }
  if (!is.null(i.stats)) {
    dplyr::mutate(.data = i.stats,
                  dplyr::across(.cols = where(is.numeric), .fns = round, digits = digits),
                  dplyr::across(.cols = where(is.numeric), .fns = format, scientific = FALSE),
    ) %>%
      readr::write_tsv(x = ., file = file.path(path.folder, i.stats.f.sum))
  }
  if (!is.null(i.info)) readr::write_tsv(x = i.info, file = file.path(path.folder, i.stats.f))
  if (!is.null(m.info)) readr::write_tsv(x = m.info, file = file.path(path.folder, m.stats.f))
  if (verbose) cli::cli_process_done()

  # return stats ---------------------------------------------------------------
  return(list(i.info = i.info, m.info = m.info, i.stats = i.stats, m.stats = m.stats, i.fig = i.fig, m.fig = m.fig))
}#End generate_stats


# Calculate individual het------------------------------------------------------
#' @title individual_het
#' @description Calculate individual heterozygosity from GDS
#' @rdname individual_het
#' @keywords internal
#' @export
individual_het <- function(gds) {

  # PLAN A

  # Longer with big datasets
  # info <- summary_gds(gds = gds, check.sync = TRUE, verbose = FALSE)
  # nonmiss <- het <- integer(info$n.ind)
  #
  # SeqArray::seqApply(
  #   gdsfile = gds,
  #   var.name = "genotype",
  #   FUN = function(x) {
  #     nm <- !is.na(x[1,]) & !is.na(x[2,])
  #     het <<- het + (x[1,] != x[2,] & nm)
  #     nonmiss <<- nonmiss + nm
  #   },
  #   margin = "by.variant",
  #   as.is = "none",
  #   parallel = FALSE
  # )
  #

  # return(het / nonmiss)


  # PLAN B
  # da <- SeqArray::seqGetData(gdsfile = gds, var.name = "$dosage_alt")
  # het <- rowSums(da == 1L, na.rm = TRUE) / rowSums(!is.na(da), na.rm = FALSE)

  # faster
  het <- SeqArray::seqGetData(gdsfile = gds, var.name = "$dosage_alt")
  rs1 <- rowSums(het == 1L, na.rm = TRUE)
  rs2 <- rowSums(!is.na(het), na.rm = FALSE)
  het <- rs1 / rs2
  return(het)
}#End individual_het

# Calculate markers het------------------------------------------------------
#' @title markers_het
#' @description Calculate markers het
#' @rdname markers_het
#' @keywords internal
#' @export
markers_het <- function(gds) {
  # PLAN A
  SeqArray::seqApply(
    gdsfile = gds,
    var.name = "$dosage_alt",
    FUN = function(x) sum(x == 1, na.rm = TRUE) / sum(!is.na(x)),
    margin = "by.variant",
    as.is = "double",
    parallel = TRUE
  )
  # PLAN B
  # not faster... strange because for sample it is faster...
  # het <- SeqArray::seqGetData(gdsfile = gds, var.name = "$dosage_alt")
  # rs1 <- colSums(het == 1L, na.rm = TRUE)
  # rs2 <- colSums(!is.na(het), na.rm = FALSE)
  # het <- rs1 / rs2

}#End markers_het

# Calculate markers FIS------------------------------------------------------
#' @title markers_fis
#' @description Calculate markers FIS
#' @rdname markers_fis
#' @keywords internal
#' @export
markers_fis <- function(gds) {
  fis <- SeqArray::seqApply(
    gdsfile = gds,
    var.name = "$dosage_alt",
    FUN = function(x) {
      c(
        REF = sum(x == 0L, na.rm = TRUE),
        HET = sum(x == 1L, na.rm = TRUE),
        ALT = sum(x > 1L, na.rm = TRUE)
      )
    },
    margin = "by.variant",
    as.is = "list",
    parallel = FALSE
  ) %>%
    dplyr::bind_rows(.) %>%
    dplyr::mutate(
      N = REF + HET + ALT,
      ALT = NULL,
      FREQ_REF = ((2 * REF) + HET) / (2 * N),
      HET_EXP = 2 * FREQ_REF * (1 - FREQ_REF) * N,
      FREQ_REF = NULL,
      FIS = 1 - (HET / HET_EXP)
    ) %>%
    dplyr::pull(FIS)
  return(fis)
}# End markers_fis


# Generate gt_vcf_nuc ----------------------------------------------------------
#' @title generate_gt_vcf_nuc
#' @description Generate gt_vcf_nuc
#' @rdname generate_gt_vcf_nuc
#' @keywords internal
#' @export
generate_gt_vcf_nuc <- function(gds) {
  gt.vcf.nuc <- SeqArray::seqApply(
    gdsfile = gds,
    var.name = c(geno = "genotype", phase = "phase", allele = "allele"),
    FUN = function(x) {
      alleles <- unlist(strsplit(x$allele, ",", fixed = TRUE), use.names = FALSE)
      names(alleles) <- 0:(length(alleles) - 1)
      a <- alleles[as.character(x$geno[1,])]
      b <- alleles[as.character(x$geno[2,])]
      sep = ifelse(x$phase, "|", "/")
      paste0(a, sep, b)
    },
    margin = "by.variant",
    as.is = "list",
    parallel = FALSE
  )

  gt.vcf.nuc <- matrix(unlist(gt.vcf.nuc, use.names = FALSE), ncol = length(gt.vcf.nuc))
  gt.vcf.nuc[gt.vcf.nuc == "NA/NA"] <- NA
}#End generate_gt_vcf_nuc



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

  # Check that SeqArray is installed
  radiator_packages_dep(package = "SeqArray", cran = FALSE, bioc = TRUE)

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

  strata <- data %>%
    dplyr::select(INDIVIDUALS, dplyr::any_of(c("STRATA", "POP_ID", "TARGET_ID"))) %>%
    dplyr::rename(STRATA = tidyselect::any_of("POP_ID"))

  want <- c("VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS", "COL", "REF",
            "ALT", "CALL_RATE", "REP_AVG", "AVG_COUNT_REF",
            "AVG_COUNT_SNP", "ONE_RATIO_REF", "ONE_RATIO_SNP", "SEQUENCE")
  markers.meta <- dplyr::select(data, dplyr::any_of(want))
  if (ncol(markers.meta) == 0) markers.meta <- NULL

  want <- c(
    "GT", "GT_VCF", "GT_VCF_NUC",
    "READ_DEPTH", "ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH",
    "GQ",
    "GL_HOM_REF", "GL_HET", "GL_HOM_ALT",
    "DP", "AD", "GL", "PL", "HQ", "GOF", "NR", "NV", "RO", "QR", "AO", "QA")
  genotypes.meta <- dplyr::select(data, dplyr::any_of(want))
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


