# write_rubias ----------------------------------------------------------------

#' @name write_rubias

#' @title Write a rubias object

#' @description Write a rubias object from a tidy data frame or GDS file/object.

#' @inheritParams radiator_common_arguments

#' @param strata (optional, tibble file or object) This tibble of individual's
#' metadata must contain
#' four columns: \code{SAMPLE_TYPE, REPUNIT, COLLECTION, INDIVIDUALS}. Those columns
#' are described in \strong{rubias}.
#' Default:\code{strata = NULL}. With default, \code{SAMPLE_TYPE} is filled
#' with \code{reference}. \code{REPUNIT} and \code{COLLECTION} will be filled by
#' the \code{STRATA} or \code{POP_ID} column found in the data.

#' @param filename The prefix for the name of the file written to the working directory.
#' Default: \code{filename = NULL}.
#' With default, only the rubias object is generated. The filename will be appended
#' \code{_rubias.tsv}.

#' @seealso \href{https://github.com/eriqande/rubias}{rubias}: genetic stock
#' identification (GSI) in the tidyverse.

#' @return A rubias object in the global environment and a file is written in the
#' working directory if \code{filename} argument was used.
#' @export
#' @rdname write_rubias
#' @references Moran BM, and Anderson, E. C. 2018 Bayesian inference from the
#' conditional genetic stock identification model.
#' Can. J. Fish. Aquat. Sci. 76: 551-560.
#' @references Anderson, Eric C., Robin S. Waples, and Steven T. Kalinowski. (2008)
#' An improved method for predicting the accuracy of genetic stock identification.
#' Canadian Journal of Fisheries and Aquatic Sciences 65, 7:1475-1486.
#' @references Anderson, E. C. (2010) Assessing the power of informative subsets of
#' loci for population assignment: standard methods are upwardly biased.
#' Molecular ecology resources 10, 4:701-710.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_rubias <- function(
  data,
  strata = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1
) {

  ## TEST
  # strata = NULL
  # filename = NULL
  # parallel.core = parallel::detectCores() - 1

  # Checking for missing and/or default arguments
  if (missing(data)) rlang::abort("Input file necessary to write the rubias file is missing")

  # Date and time --------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")

  # Filename -------------------------------------------------------------------
  if (!is.null(filename)) {
    filename.problem <- file.exists(filename)
    if (filename.problem) {
      filename <- stringi::stri_join(filename, "_rubias_", file.date, ".tsv")
    } else {
      filename <- stringi::stri_join(filename, "_rubias.tsv")
    }
  }

  # File type detection----------------------------------------------------------
  data.type <- radiator::detect_genomic_format(data)

  # Import data ---------------------------------------------------------------

  if (data.type %in% c("SeqVarGDSClass", "gds.file")) {
    if (data.type == "gds.file") data %<>% radiator::read_rad(data = ., verbose = FALSE)
    strata.id <- radiator::extract_individuals_metadata(
      gds = data,
      ind.field.select = c("STRATA", "INDIVIDUALS"),
      whitelist = TRUE
    )

    markers <- radiator::extract_markers_metadata(
      gds = data,
      # markers.meta.select = c("VARIANT_ID", "MARKERS"),
      markers.meta.select = "MARKERS",
      whitelist = TRUE
    ) %$% MARKERS


    id <- SeqArray::seqGetData(gdsfile = data, var.name = "sample.id")
    # checks
    if (!identical(id, strata.id$INDIVIDUALS)) {
      rlang::abort("Problem with GDS, contact author")
    }
    id <- NULL

    data <- generate_gt_vcf_nuc(data)
    if (length(markers) != ncol(data)) {
      rlang::abort("Problem with function conversion, contact author")
    }

    # generate the appropriate colnames with .A1 and .A2 appended
    markers %<>%
      purrr::map(.x = ., .f = function(x) paste0(x, c(".A1", ".A2"))) %>%
      purrr::flatten_chr(.)

    data <- dplyr::bind_cols(
      strata.id,
      data  %<>%
        radiator_split_tibble(x = ., parallel.core = parallel.core) %>%
        magrittr::set_colnames(markers)
    )

    markers <- strata.id <- NULL

  } else {
    if (is.vector(data)) data %<>% radiator::tidy_wide(data = ., import.metadata = TRUE)
    split.chunks <- 1L
    if (parallel.core > 1) split.chunks <- 3L

    if (rlang::has_name(data, "GT_BIN")) {
      want <- c("MARKERS", "POP_ID", "INDIVIDUALS", "GT_BIN")
      gt.format <- "GT_BIN"
    } else {
      if (!rlang::has_name(data, "GT")) {
        data %<>% calibrate_alleles(data = ., verbose = FALSE, gt = TRUE) %$% input
        want <- c("MARKERS", "POP_ID", "INDIVIDUALS", "GT")
        gt.format <- "GT"
      }
    }

    data %<>%
      dplyr::select(tidyselect::any_of(want)) %>%
      separate_gt(
        x = .,
        gt = gt.format,
        gather = TRUE,
        exclude = c("MARKERS", "POP_ID", "INDIVIDUALS"),
        remove = TRUE,
        split.chunks = split.chunks
      ) %>%
      dplyr::mutate(
        MARKERS_ALLELES = stringi::stri_join(MARKERS, ALLELES_GROUP, sep = "."),
        MARKERS = NULL, ALLELES_GROUP = NULL
      ) %>%
      {if (gt.format == "GT") dplyr::mutate(.data = ., ALLELES = replace(x = ALLELES, which(ALLELES == "000"), NA)) else .} %>%
      rad_wide(x = ., formula = "POP_ID + INDIVIDUALS ~ MARKERS_ALLELES", values_from = "ALLELES")
  }

  # individuals metadata -------------------------------------------------------
  if (!is.null(strata)) {
    strata <- radiator::read_strata(
      strata = strata,
      pop.id = TRUE,
      verbose = FALSE, keep.two = FALSE
    ) %$%
      strata

    # check strata -------------------------------------------------------------
    if (ncol(strata) < 4 || FALSE %in% (c("SAMPLE_TYPE", "REPUNIT", "COLLECTION", "INDIVIDUALS") %in% colnames(strata))) {
      message("The strata file/object requires 4 columns: SAMPLE_TYPE, REPUNIT, COLLECTION, INDIVIDUALS")
      message("SAMPLE_TYPE filled with reference, REPUNIT and COLLECTION filled with STRATA/POP_ID column")

      strata %<>%
        dplyr::distinct(POP_ID, INDIVIDUALS) %>%
        dplyr::mutate(
          sample_type = "reference",
          repunit = as.character(POP_ID),
          collection = as.character(POP_ID),
          indiv = INDIVIDUALS,
          POP_ID = NULL
        )

    } else {
      strata %<>%
        dplyr::rename(
          sample_type = SAMPLE_TYPE,
          repunit = REPUNIT,
          collection = COLLECTION
        ) %>%
        dplyr::mutate(indiv = INDIVIDUALS)
    }



  } else {
    if (!rlang::has_name(data, "POP_ID") && !rlang::has_name(data, "STRATA")) {
      rlang::abort("function requires data with STRATA or POP_ID information")
    }
    if (rlang::has_name(data, "POP_ID")) data %<>% dplyr::rename(STRATA = POP_ID)

    strata <- dplyr::distinct(data, STRATA, INDIVIDUALS) %>%
      dplyr::mutate(
        sample_type = "reference",
        repunit = as.character(STRATA),
        collection = as.character(STRATA),
        indiv = INDIVIDUALS,
        STRATA = NULL
      )
  }

  # check
  if (!identical(strata$INDIVIDUALS, data$INDIVIDUALS)) {
    rlang::abort("Problem with individuals metadata, contact author")
  }
  if (rlang::has_name(data, "POP_ID")) data$POP_ID <- NULL
  if (rlang::has_name(data, "STRATA")) data$STRATA <- NULL
  data$INDIVIDUALS <- NULL
  strata$INDIVIDUALS <- NULL

  # merge data and metadata-----------------------------------------------------
  data <- dplyr::bind_cols(strata, data)
  strata <- NULL

  # checks for mixture----------------------------------------------------------
  if ("mixture" %in% unique(data$sample_type)) {
    data %<>%
      dplyr::mutate(
        repunit = dplyr::if_else(sample_type == "mixture", NA_character_, repunit)
      )
  }

  # write to disk---------------------------------------------------------------
  if (!is.null(filename)) readr::write_tsv(x = data, file = filename)
  return(data)
} # End write_rubias function
