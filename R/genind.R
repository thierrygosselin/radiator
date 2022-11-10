# tidy_genind ------------------------------------------------------------------
#' @name tidy_genind
#' @title Tidy a genind object to a tidy dataframe
#' @description Tidy genind object from
#' \href{https://github.com/thibautjombart/adegenet}{adegenet} to a tidy dataframe.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data (path or object) A genind object in the global environment or
#' path to a genind file that will be open with \code{readRDS}.

#' @param keep.allele.names Allows to keep allele names for the tidy dataset.
#' Requires the alleles to be numeric. To have this argument in
#' \pkg{radiator} \code{\link{tidy_genomic_data}} or
#' \pkg{radiator} \code{\link{genomic_converter}}
#' use it at the end. \code{...} in those function looks for it.
#' Default: \code{keep.allele.names = FALSE}.

#' @param tidy (logical) Generate a tidy dataset.
#' Default: \code{tidy = TRUE}.


#' @param gds (optional, logical) To write a radiator gds object.
#' Currently, for biallelic datasets only.
#' Default: \code{gds = TRUE}.

#' @param write (optional, logical) To write in the working directory the tidy
#' data. The file is written with \code{radiator_genind_DATE@TIME.rad}.
#' Default: \code{write = FALSE}.

#' @inheritParams radiator_common_arguments


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
    keep.allele.names = FALSE,
    tidy = TRUE,
    gds = TRUE,
    write = FALSE,
    verbose = FALSE
) {

  # TEST
  # keep.allele.names = TRUE
  # tidy = TRUE
  # gds = FALSE
  # write = FALSE
  # verbose = TRUE

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("data argument required")
  if (is.vector(data)) data <- readRDS(data)
  if (class(data)[1] != "genind") rlang::abort("Input is not a genind object")
  if (verbose) message("genind info:")
  strata <- tibble::tibble(INDIVIDUALS = rownames(data@tab))
  if (is.null(data@pop)) {
    if (verbose) message("    strata: no")
    if (verbose) message("    'pop' will be added")
    strata %<>% dplyr::mutate(STRATA = "pop")
    nostrata <- TRUE
  } else {
    if (verbose) message("    strata: yes")
    strata$STRATA = data@pop
    nostrata <- FALSE
  }

  biallelic <- max(unique(data@loc.n.all)) == 2
  if (!biallelic) gds <- FALSE
  if (gds) tidy <- TRUE

  if (tidy) {
    if (write) {
      filename.temp <- generate_filename(extension = "rad")
      filename.short <- filename.temp$filename.short
      filename.genind <- filename.temp$filename
    }

    # detect A2... facilitate the conversion if made by radiator...
    A2 <- TRUE %in% unique(stringi::stri_detect_fixed(str = sample(colnames(data@tab), min(ncol(data@tab), 100)), pattern = ".A2"))
    if (isTRUE(TRUE %in% unique(stringi::stri_detect_fixed(str = sample(colnames(data@tab), min(ncol(data@tab), 100)), pattern = ".A2__")))) {
      A2 <- FALSE
    }
    codom <- FALSE
    if (data@type == "codom") codom <- TRUE

    data.temp <- tibble::as_tibble(data@tab) %>%
      tibble::add_column(.data = ., INDIVIDUALS = rownames(data@tab), .before = 1)
    if (nostrata) {
      data.temp %<>% tibble::add_column(STRATA = "pop")
    } else {
      data.temp %<>% tibble::add_column(STRATA = data@pop)
    }
    data <- data.temp
    data.temp <- NULL

    if (biallelic) {
      if (A2) {
        # changed adegenet::indNames to rownames(data@tab) to lower dependencies
        data %<>%
          dplyr::select(STRATA, INDIVIDUALS, dplyr::ends_with(match = ".A2")) %>%
          radiator::rad_long(
            x = .,
            cols = c("INDIVIDUALS", "STRATA"),
            names_to = "MARKERS",
            values_to = "GT_BIN"
          ) %>%
          dplyr::mutate(
            MARKERS = stringi::stri_replace_all_fixed(
              str = MARKERS, pattern = ".A2", replacement = "", vectorize_all = FALSE),
            GT_VCF = dplyr::if_else(GT_BIN == 0, "0/0",
                                    dplyr::if_else(GT_BIN == 1, "0/1", "1/1"), missing = "./."),
            GT = dplyr::if_else(GT_BIN == 0, "001001", dplyr::if_else(GT_BIN == 1, "001002", "002002") , missing = "000000")
          )
      } else {
        if (codom) {
          data %<>%
            radiator::rad_long(
              x = .,
              cols = c("INDIVIDUALS", "STRATA"),
              names_to = "MARKERS_ALLELES",
              values_to = "GT_BIN"
            ) %>%
            tidyr::separate(data = ., col = MARKERS_ALLELES, into = c("MARKERS", "ALLELES"), sep = "\\.") %>%
            dplyr::mutate(
              ALLELES = NULL,
              GT_VCF = dplyr::if_else(GT_BIN == 0, "0/0",
                                      dplyr::if_else(GT_BIN == 1, "0/1", "1/1"), missing = "./."),
              GT = dplyr::if_else(GT_BIN == 0, "001001", dplyr::if_else(GT_BIN == 1, "001002", "002002") , missing = "000000")
            )
        }
      }
    } else {

      if (keep.allele.names) {
        message("Alleles names are kept if all numeric and padded with 0 if length < 3")
        data %<>%
          radiator::rad_long(
            x = .,
            cols = c("INDIVIDUALS", "STRATA"),
            names_to = "MARKERS_ALLELES",
            values_to = "COUNT"
          ) %>%
          dplyr::filter(COUNT > 0 | is.na(COUNT)) %>%
          tidyr::separate(data = ., col = MARKERS_ALLELES, into = c("MARKERS", "ALLELES"), sep = "\\.") %>%
          dplyr::mutate(
            ALLELES = stringi::stri_replace_all_fixed(str = ALLELES, pattern = c("A1__", "A2__"), replacement = c("", ""), vectorize_all = FALSE),
          )

        check.alleles <- unique(stringi::stri_detect_regex(str = unique(data$ALLELES), pattern = "[0-9]"))
        check.alleles <- length(check.alleles) == 1 && check.alleles

        if (check.alleles) {
          data %<>%
            dplyr::mutate(
              ALLELES = stringi::stri_pad_left(str = ALLELES, pad = "0", width = 3)
            )
        } else {
          data %<>%
            dplyr::mutate(
              ALLELES = as.numeric(factor(ALLELES)),
              ALLELES = stringi::stri_pad_left(str = ALLELES, pad = "0", width = 3)
            )
        }

      } else {
        message("Alleles names for each markers will be converted to factors and padded with 0")
        data %<>%
          radiator::rad_long(
            x = .,
            cols = c("INDIVIDUALS", "STRATA"),
            names_to = "MARKERS_ALLELES",
            values_to = "COUNT"
          ) %>%
          dplyr::filter(COUNT > 0 | is.na(COUNT)) %>%
          tidyr::separate(data = ., col = MARKERS_ALLELES, into = c("MARKERS", "ALLELES"), sep = "\\.") %>%
          dplyr::mutate(
            ALLELES = stringi::stri_replace_all_fixed(str = ALLELES, pattern = c("A1__", "A2__"), replacement = c("", ""), vectorize_all = FALSE),
            ALLELES = as.numeric(factor(ALLELES)),
            ALLELES = stringi::stri_pad_left(str = ALLELES, pad = "0", width = 3)
          )
      }

      # #If the genind was coded with allele 0, this will generate missing data with this code
      # allele.zero <- TRUE %in% stringi::stri_detect_regex(str = unique(data3$ALLELES), pattern = "^0$")
      # if (allele.zero) stop("alleles in this multiallelic genind were coded as 0 and won't work with this script")

      #Isolate missing genotype
      missing.hom <- dplyr::filter(data, is.na(COUNT)) %>%
        dplyr::distinct(STRATA, INDIVIDUALS, MARKERS) %>%
        dplyr::mutate(GT = rep("000000", n())) %>%
        dplyr::ungroup(.)

      #Isolate all Homozygote genotypes and combine with the missings
      missing.hom <- dplyr::filter(data, COUNT == 2) %>%
        dplyr::group_by(STRATA, INDIVIDUALS, MARKERS) %>%
        dplyr::summarise(GT = stringi::stri_join(ALLELES, ALLELES, sep = "")) %>%
        dplyr::ungroup(.) %>%
        dplyr::bind_rows(missing.hom)

      #Isolate all Het genotypes and combine
      data <- dplyr::filter(data, COUNT != 2) %>%#this will also remove the missing
        dplyr::group_by(STRATA, INDIVIDUALS, MARKERS) %>%
        dplyr::summarise(GT = stringi::stri_join(ALLELES, collapse = "")) %>%
        dplyr::ungroup(.) %>%
        dplyr::bind_rows(missing.hom)
      missing.hom <- NULL
    }#End for multi-allelic or weird genind

    if (write) radiator::write_rad(data = data, path = filename.genind, verbose = verbose)

  }# End tidy genind
  if (gds) {

    chrom.info <- FALSE
    locus.info <- FALSE
    pos.info <- FALSE
    n.markers <- dplyr::n_distinct(data$MARKERS)
    n.ind <- dplyr::n_distinct(data$INDIVIDUALS)

    markers <- data %>%
      # dplyr::group_by(MARKERS) %>%
      dplyr::distinct(MARKERS) %>%
      dplyr::mutate(
        VARIANT_ID = as.integer(factor(MARKERS)),
        CHROM = factor(rep("CHROM1", n.markers)),
        LOCUS = MARKERS # this is usually the adegenet behavior
        # POS = vctrs::vec_cast(seq(from = 1, to = n.markers, by = 1), integer()) # need to remove the dependecy
      )


    data %<>%
      dplyr::mutate(
        INDIVIDUALS = radiator::clean_ind_names(INDIVIDUALS),
        STRATA = radiator::clean_pop_names(STRATA)
      ) %>%
      dplyr::arrange(MARKERS, STRATA, INDIVIDUALS)

    gds.filename <- radiator_gds(
      data.source = "genind",
      genotypes = gt2array(
        gt.bin = data$GT_BIN,
        n.ind = n.ind,
        n.snp = n.markers
      ),
      strata = strata,
      biallelic = TRUE,
      markers.meta = markers,
      filename = NULL,
      verbose = verbose
    )





    # doesnt work with all genind object + GDS requires an array...
    # # prepare genind
    # alt.alleles <- tibble::tibble(
    #   MARKERS_ALLELES = colnames(data@tab),
    #   COUNT = colSums(x = data@tab, na.rm = TRUE)
    # ) %>%
    #   dplyr::mutate(
    #     MARKERS = stringi::stri_extract_first_regex(str = colnames(data@tab), pattern = "^[^.]+"),
    #     ALLELES = stringi::stri_extract_last_regex(str = colnames(data@tab), pattern = "(?<=\\.).*"),
    #     ALLELES = stringi::stri_replace_all_fixed(str = ALLELES, pattern = c("A1__", "A2__"), replacement = c("", ""), vectorize_all = FALSE),
    #   ) %>%
    #   dplyr::arrange(MARKERS, COUNT, ALLELES) %>%
    #   dplyr::mutate(REF = rep(c("ALT", "REF"), n() / 2))
    #
    # # check that REF/ALT are A, C, G, T before going further
    # # any(unique(alt.alleles$ALLELES) %in% c("A", "C", "G", "T"))
    # # any(unique(c("A", "DD", "G")) %in% c("A", "C", "G", "T"))
    #
    # ref.alt <- dplyr::select(alt.alleles, MARKERS, REF, ALLELES) %>%
    #   radiator::rad_wide(x = ., formula = "MARKERS ~ REF", values_from = "ALLELES")
    #
    # alt.alleles %<>% dplyr::filter(REF == "ALT") %$% MARKERS_ALLELES
    #
    # geno <- tibble::as_tibble(t(data@tab), rownames = "MARKERS") %>%
    #   dplyr::filter(MARKERS %in% alt.alleles) %>%
    #   dplyr::mutate(
    #     MARKERS = stringi::stri_extract_first_regex(str = MARKERS, pattern = "^[^.]+"),
    #     MARKERS = stringi::stri_replace_all_fixed(
    #       str = MARKERS,
    #       pattern = c("__A1", "__A2"),
    #       replacement = c("", ""),
    #       vectorize_all = FALSE
    #     ),
    #     VARIANT_ID = as.integer(factor(MARKERS))) %>%
    #   dplyr::arrange(VARIANT_ID)
    #
    # alt.alleles <- NULL
    # markers.meta <- dplyr::select(geno, VARIANT_ID, MARKERS) %>%
    #   dplyr::left_join(ref.alt, by = "MARKERS")
    # ref.alt <- NULL
    #
    # suppressWarnings(
    #   geno %<>%
    #     dplyr::select(-MARKERS) %>%
    #     tibble::column_to_rownames(.data = ., var = "VARIANT_ID")
    # )
    #
    # gds.filename <- radiator_gds(
    #   data.source = "genind",
    #   genotypes = geno,
    #   strata = strata,
    #   biallelic = TRUE,
    #   markers.meta = markers.meta,
    #   filename = NULL,
    #   open = FALSE,
    #   verbose = verbose
    # )
    if (verbose) message("Written: GDS filename: ", gds.filename)
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
  data <- NULL

  if (write) {
    filename.temp <- generate_filename(extension = "genind")
    filename.short <- filename.temp$filename.short
    filename.genind <- filename.temp$filename
    saveRDS(object = res, file = filename.genind)
    if (verbose) message("File written: ", filename.short)
  }

  return(res)
} # End write_genind
