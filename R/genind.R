# tidy_genind ------------------------------------------------------------------
#' @name tidy_genind
#' @title Tidy a genind object to a tidy dataframe
#' @description Tidy genind object from
#' \href{https://github.com/thibautjombart/adegenet}{adegenet} to a tidy dataframe.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data (path or object) A genind object in the global environment or
#' path to a genind file that will be open with \code{readRDS}.

#' @param tidy (logical) Generate a tidy dataset.
#' Default: \code{tidy = TRUE}.

#' @param gds (optional, logical) To write a radiator gds object.
#' Currently, for biallelic datasets only.
#' Default: \code{gds = TRUE}.

#' @param write (optional, logical) To write in the working directory the tidy
#' data. The file is written with \code{radiator_genind_DATE@TIME.rad}.
#' Default: \code{write = FALSE}.

#' @inheritParams radiator_common_arguments

#' @note \href{https://github.com/thibautjombart/adegenet}{genind} objects, like
#' genepop, are not optimal genomic format for RADseq datasets,
#' they lack important genotypes and markers metadata: chromosome, locus, snp,
#' position, read depth, allele depth, etc.
#' \href{https://github.com/thibautjombart/adegenet}{genlight} object is a more
#' interesting container and is memory efficient, see \code{\link{tidy_genlight}}.
#'
#'
#' By default allele names will be kept for the tidy dataset,
#' if the alleles is numeric and length < 3.
#'
#'
#' In the unlikely event that the genind object as no stratification/population,
#' \emph{pop} will be added to the strata column.


#' @export
#' @rdname tidy_genind
#' @references Jombart T (2008) adegenet: a R package for the multivariate
#' analysis of genetic markers. Bioinformatics, 24, 1403-1405.
#' @references Jombart T, Ahmed I (2011) adegenet 1.3-1:
#' new tools for the analysis of genome-wide SNP data.
#' Bioinformatics, 27, 3070-3071.


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


tidy_genind <- function(
    data,
    tidy = TRUE,
    gds = TRUE,
    write = FALSE,
    verbose = FALSE
) {

  # TEST
  # tidy = TRUE
  # gds = TRUE
  # write = FALSE
  # verbose = TRUE
  if (verbose) cli::cli_progress_step("Reading genind")

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("data argument required")
  if (is.vector(data)) data <- readRDS(data)
  if (class(data)[1] != "genind") rlang::abort("Input is not a genind object")


  # Working on individuals and strata ------------------------------------------
  strata <- tibble::tibble(INDIVIDUALS = rownames(data@tab))
  if (is.null(data@pop)) {
    strata %<>% dplyr::mutate(STRATA = "pop")
  } else {
    strata$STRATA = data@pop
  }

  n.ind <- length(strata$INDIVIDUALS)
  n.snp <- length(unique(data@loc.fac))

  biallelic <- max(unique(data@loc.n.all)) == 2
  if (!biallelic) gds <- FALSE
  if (gds) tidy <- TRUE

  if (write) {
    filename.temp <- generate_filename(extension = "rad")
    filename.short <- filename.temp$filename.short
    filename.genind <- filename.temp$filename
  }

  if (tidy) {
    # detect A2... facilitate the conversion if made by radiator...
    A2 <- FALSE
    codom <- FALSE
    colnames.coding <- sample(colnames(data@tab), min(ncol(data@tab), 100))
    A2 <- TRUE %in% unique(stringi::stri_detect_fixed(str = colnames.coding, pattern = ".A2"))
    if (data@type == "codom") codom <- TRUE

    # data.temp <- tibble::as_tibble(data@tab) %>%
    #   tibble::add_column(.data = ., INDIVIDUALS = rownames(data@tab), .before = 1)
    # if (nostrata) {
    #   data.temp %<>% tibble::add_column(STRATA = "pop")
    # } else {
    #   data.temp %<>% tibble::add_column(STRATA = data@pop)
    # }
    # data <- data.temp
    # data.temp <- NULL

    if (biallelic) {
      if (verbose) cli::cli_progress_step("Preparing biallelic tidy data")
      markers.meta <- colnames(data@tab)
      n.snp <- length(markers.meta)/2

      alt.alleles <- tibble::tibble(
        MARKERS_ALLELES = markers.meta,
        COUNT = colSums(x = data@tab, na.rm = TRUE)
      ) %>%
        dplyr::mutate(
          MARKERS = stringi::stri_extract_first_regex(str = markers.meta, pattern = "^[^.]+"),
          ALLELES = stringi::stri_extract_last_regex(str = markers.meta, pattern = "(?<=\\.).*"),
          ALLELES = stringi::stri_replace_all_fixed(str = ALLELES, pattern = c("A1__", "A2__"), replacement = c("", ""), vectorize_all = FALSE),
        ) %>%
        dplyr::arrange(MARKERS, COUNT, ALLELES) %>%
        dplyr::mutate(REF = rep(c("ALT", "REF"), n() / 2))

      # check that REF/ALT are A, C, G, T before going further .....
      discard.alleles <- FALSE
      unique.alleles <- unique(alt.alleles$ALLELES)
      if (length(unique.alleles) > 4) discard.alleles <- TRUE
      if (any(stringi::stri_detect_regex(str = unique.alleles, pattern = "[0-9]"))) discard.alleles <- TRUE

      if (!discard.alleles) {
        ref.alt <- dplyr::select(alt.alleles, MARKERS, REF, ALLELES) %>%
          radiator::rad_wide(x = ., formula = "MARKERS ~ REF", values_from = "ALLELES")
      }

      alt.alleles %<>% dplyr::filter(REF == "ALT") %$% MARKERS_ALLELES
      if (length(alt.alleles) != n.snp) rlang::abort("Contact author problem with tidying the genind")

      # test1 <- tibble::as_tibble(t(data@tab), rownames = "MARKERS") %>%
      #   dplyr::filter(MARKERS %in% alt.alleles) %>%
      #   dplyr::arrange(MARKERS)
      # test2 <- tibble::as_tibble(t(data@tab[, alt.alleles]), rownames = "MARKERS") %>% dplyr::arrange(MARKERS)
      # test3 <- tibble::as_tibble(t(data@tab[, seq(2, length(markers.meta), by = 2)]), rownames = "MARKERS") %>% dplyr::arrange(MARKERS)

      data <- tibble::as_tibble(t(data@tab), rownames = "MARKERS") %>%
        dplyr::filter(MARKERS %in% alt.alleles) %>%
        dplyr::mutate(
          MARKERS = stringi::stri_extract_first_regex(str = MARKERS, pattern = "^[^.]+"),
          MARKERS = stringi::stri_replace_all_fixed(
            str = MARKERS,
            pattern = c("__A1", "__A2"),
            replacement = c("", ""),
            vectorize_all = FALSE
          ),
          VARIANT_ID = as.integer(factor(MARKERS))) %>%
        dplyr::arrange(VARIANT_ID)

      alt.alleles <- NULL
      markers.meta <- dplyr::distinct(data, VARIANT_ID, MARKERS)

      if (!discard.alleles) {
        markers.meta %<>%
          dplyr::left_join(ref.alt, by = "MARKERS")
      }
      ref.alt <- NULL

      # generate markers metadata
      markers.meta %<>% separate_markers(data = ., generate.ref.alt = TRUE, biallelic = TRUE)

      data <- radiator::rad_long(
        x = data,
        cols = "MARKERS",
        names_to = "INDIVIDUALS",
        values_to = "GT_BIN"
      ) %>%
        dplyr::left_join(strata, by = "INDIVIDUALS") %>%
        dplyr::left_join(markers.meta, by = "MARKERS") %>%
        dplyr::mutate(
          INDIVIDUALS = radiator::clean_ind_names(INDIVIDUALS),
          STRATA = radiator::clean_pop_names(STRATA)
        ) %>%
        dplyr::arrange(MARKERS, STRATA, INDIVIDUALS)

    } else {
      if (verbose) cli::cli_progress_step("Preparing multi-allelic tidy data")
      data <- radiator::rad_long(
          x = tibble::as_tibble(data@tab, rownames = "INDIVIDUALS"),
          cols = "INDIVIDUALS",
          names_to = c("MARKERS", "ALLELES"),
          values_to = "COUNT",
          names_sep = "\\.",
          tidy = TRUE
        ) %>%
        dplyr::left_join(strata, by = "INDIVIDUALS")

      check.alleles <- unique(stringi::stri_detect_regex(str = unique(data$ALLELES), pattern = "[0-9]"))
      if (!check.alleles) data %<>% dplyr::mutate(ALLELES = as.numeric(factor(ALLELES)))
      data %<>%
        dplyr::mutate(ALLELES = stringi::stri_pad_left(str = ALLELES, pad = "0", width = 3)) %>%
        dplyr::filter(is.na(COUNT) | COUNT != 0)

      multi_genind <- function(x) {
        count.type <- unique(x$COUNT)
        if (is.na(count.type)) {
          x %<>%
            dplyr::distinct(STRATA, INDIVIDUALS, MARKERS) %>%
            dplyr::mutate(GT = rep("000000", dplyr::n())) %>%
            dplyr::ungroup(.)
        } else if (count.type == 2L) {
          x %<>%
            dplyr::group_by(STRATA, INDIVIDUALS, MARKERS) %>%
            dplyr::summarise(GT = stringi::stri_join(ALLELES, ALLELES, sep = ""), .groups = "drop")
        } else {
          x %<>%
            dplyr::group_by(STRATA, INDIVIDUALS, MARKERS) %>%
            dplyr::summarise(GT = stringi::stri_join(ALLELES, collapse = ""), .groups = "drop")
        }
      }

      if (n.ind * n.snp >= 5000000) {
        data <- radiator_future(
          .x = data,
          .f = multi_genind,
          flat.future = "dfr",
          split.with = "COUNT"
        )
      } else {
        data %<>%
          dplyr::group_split(COUNT) %>%
          purrr::map_dfr(.x = ., .f = multi_genind)
      }
    }#End for multi-allelic or weird genind

    if (write) radiator::write_rad(data = data, path = filename.genind, verbose = verbose)

  }# End tidy genind
  if (gds) {
    if (verbose) cli::cli_progress_step("Working on the GDS dataset")
    if (!rlang::has_name(x = data, "GT_BIN")) rlang::abort("Missing GT_BIN format to generate GDS, contact author")
    n.snp <- dplyr::n_distinct(data$MARKERS)
    n.ind <- dplyr::n_distinct(data$INDIVIDUALS)

    gds.filename <- radiator_gds(
      data.source = "genind",
      geno.coding = "alt.dos",
      genotypes = gt2array(
        gt.bin = data$GT_BIN,
        n.ind = n.ind,
        n.snp = n.snp
      ),
      strata = strata,
      biallelic = TRUE,
      markers.meta = markers.meta,
      filename = NULL,
      verbose = verbose
    )
  }# End gds genind

  return(data)
} # End tidy_genind


# write_genind ------------------------------------------------------------------

#' @name write_genind
#' @title Write a genind object from a tidy data frame or GDS file or object.

#' @description Write a genind object from a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @inheritParams radiator_common_arguments

#' @param write (logical, optional) To write in the working directory the genind
#' object. The file is written with \code{radiator_genind_DATE@TIME.RData} and
#' can be open with load or readRDS.
#' Default: \code{write = FALSE}.

#' @export
#' @rdname write_genind

#' @references Jombart T (2008) adegenet: a R package for the multivariate
#' analysis of genetic markers. Bioinformatics, 24, 1403-1405.
#' @references Jombart T, Ahmed I (2011) adegenet 1.3-1:
#' new tools for the analysis of genome-wide SNP data.
#' Bioinformatics, 27, 3070-3071.


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_genind <- function(data, write = FALSE, verbose = FALSE) {


  # TEST
  # write = TRUE
  # verbose = TRUE


  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("Input file missing")

  # File type detection----------------------------------------------------------
  data.type <- radiator::detect_genomic_format(data)


  # Import data ---------------------------------------------------------------

  if (data.type %in% c("SeqVarGDSClass", "gds.file")) {
    if (data.type == "gds.file") data %<>% radiator::read_rad(data = ., verbose = verbose)
    data <- gds2tidy(gds = data, parallel.core = parallel::detectCores() - 1)
    data.type <- "tbl_df"
  } else {
    want <- c("MARKERS", "STRATA", "INDIVIDUALS", "REF", "ALT", "GT", "GT_BIN")
    data %<>%
      radiator::tidy_wide(data = ., import.metadata = TRUE) %>%
      dplyr::select(tidyselect::any_of(want))
  }

  if (is.factor(data$STRATA)) {
    pop.levels <- levels(data$STRATA)
  } else {
    pop.levels <- unique(data$STRATA)
  }
  # Make sure that STRATA and INDIVIDUALS are character
  data$INDIVIDUALS <- as.character(data$INDIVIDUALS)
  data$STRATA <- as.character(data$STRATA)

  # Isolate the strata
  # we convert STRATA to factor because adegenet does it automatically...
  pop.num <- unique(stringi::stri_detect_regex(str = pop.levels, pattern = "^[0-9]+$"))
  if (length(pop.num) == 1 && pop.num) pop.levels <- as.character(sort(as.numeric(pop.levels)))

  # When GT_BIN available
  if (rlang::has_name(data, "GT_BIN")) {
    if (rlang::has_name(data, "REF")) {
      data <- dplyr::bind_rows(
        dplyr::select(data, MARKERS, STRATA, INDIVIDUALS, REF, n = GT_BIN) %>%
          dplyr::mutate(
            REF = stringi::stri_join("A1", REF, sep = "__"),
            MARKERS_ALLELES = stringi::stri_join(MARKERS, REF, sep = "."),
            MARKERS = NULL,
            REF = NULL,
            n = as.integer(abs(n - 2))
          ),
        dplyr::select(data, MARKERS, STRATA, INDIVIDUALS, ALT, n = GT_BIN) %>%
          dplyr::mutate(
            ALT = stringi::stri_join("A2", ALT, sep = "__"),
            MARKERS_ALLELES = stringi::stri_join(MARKERS, ALT, sep = "."),
            MARKERS = NULL,
            ALT = NULL
          )
      ) %>%
        radiator::rad_wide(
          x = .,
          formula = "STRATA + INDIVIDUALS ~ MARKERS_ALLELES",
          values_from = "n"
        ) %>%
        dplyr::mutate(STRATA = factor(as.character(STRATA), levels = pop.levels)) %>%# xvalDapc doesn't accept pop as ordered factor
        dplyr::arrange(STRATA, INDIVIDUALS)
    } else {
      data <- dplyr::bind_rows(
        dplyr::select(data, MARKERS, STRATA, INDIVIDUALS, n = GT_BIN) %>%
          dplyr::mutate(
            MARKERS_ALLELES = stringi::stri_join(MARKERS, "A1", sep = "."),
            MARKERS = NULL,
            REF = NULL,
            n = as.integer(abs(n - 2))
          ),
        dplyr::select(data, MARKERS, STRATA, INDIVIDUALS, n = GT_BIN) %>%
          dplyr::mutate(
            MARKERS_ALLELES = stringi::stri_join(MARKERS, "A2", sep = "."),
            MARKERS = NULL,
            ALT = NULL
          )
      ) %>%
        radiator::rad_wide(
          x = .,
          formula = "STRATA + INDIVIDUALS ~ MARKERS_ALLELES",
          values_from = "n"
        ) %>%
        dplyr::mutate(STRATA = factor(as.character(STRATA), levels = pop.levels)) %>%# xvalDapc doesn't accept pop as ordered factor
        dplyr::arrange(STRATA, INDIVIDUALS)
    }
  } else {
    missing.geno <- dplyr::ungroup(data) %>%
      dplyr::select(MARKERS, INDIVIDUALS, GT) %>%
      dplyr::filter(GT == "000000") %>%
      dplyr::select(MARKERS, INDIVIDUALS) #%>%  dplyr::mutate(MISSING = rep("blacklist", n()))

    data <- suppressWarnings(
      dplyr::ungroup(data) %>%
        dplyr::select(MARKERS, INDIVIDUALS, GT, STRATA) %>%
        dplyr::filter(GT != "000000") %>%
        dplyr::mutate(
          A1 = stringi::stri_sub(str = GT, from = 1, to = 3),
          A2 = stringi::stri_sub(str = GT, from = 4, to = 6),
          GT = NULL
        ) %>%
        radiator::rad_long(
          x = .,
          cols = c("INDIVIDUALS", "STRATA", "MARKERS"),
          names_to = "ALLELES",
          values_to = "GT"
        ) %>%
        dplyr::arrange(MARKERS, STRATA, INDIVIDUALS, GT) %>%
        dplyr::count(x = ., STRATA, INDIVIDUALS, MARKERS, GT) %>%
        dplyr::ungroup(.) %>%
        tidyr::complete(data = ., tidyr::nesting(INDIVIDUALS, STRATA), tidyr::nesting(MARKERS, GT), fill = list(n = 0)) %>%
        dplyr::mutate(MARKERS_ALLELES = stringi::stri_join(MARKERS, GT, sep = ".")) %>%
        dplyr::anti_join(missing.geno, by = c("MARKERS", "INDIVIDUALS")) %>%
        dplyr::select(-MARKERS, -GT) %>%
        dplyr::mutate(STRATA = factor(as.character(STRATA), levels = pop.levels)) %>%# xvalDapc doesn't accept pop as ordered factor
        dplyr::arrange(MARKERS_ALLELES, INDIVIDUALS) %>%
        radiator::rad_wide(x = ., formula = "STRATA + INDIVIDUALS ~ MARKERS_ALLELES", values_from = "n") %>%
        dplyr::arrange(STRATA, INDIVIDUALS))
  }

  strata.genind <- dplyr::distinct(.data = data, INDIVIDUALS, STRATA) %>%
    dplyr::mutate(INDIVIDUALS = factor(INDIVIDUALS, levels = unique(data$INDIVIDUALS)))

  # genind arguments common to all data.type
  ind <- data$INDIVIDUALS
  pop <- data$STRATA
  data <-  dplyr::ungroup(data) %>%
    dplyr::select(-c(INDIVIDUALS, STRATA))
  suppressWarnings(rownames(data) <- ind)

  # genind constructor
  prevcall <- match.call()
  res <- adegenet::genind(
    tab = data,
    pop = pop,
    prevcall = prevcall,
    ploidy = 2,
    type = "codom",
    strata = strata.genind,
    hierarchy = NULL
  )
  data <- strata.genind <- NULL

  if (write) {
    filename.temp <- generate_filename(extension = "genind")
    filename.short <- filename.temp$filename.short
    filename.genind <- filename.temp$filename
    saveRDS(object = res, file = filename.genind)
    if (verbose) message("File written: ", filename.short)
  }

  return(res)
} # End write_genind
