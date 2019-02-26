# private_alleles

#' @name private_alleles
#' @title Find private alleles
#' @description The function highlight private alleles by strata, using a GDS or tidy file or object.

#' @inheritParams read_strata
#' @inheritParams radiator_common_arguments

#' @return A list with an object highlighting private alleles by markers and strata and
#' a second object with private alleles summarized by strata.

#' @examples
#' \dontrun{
#' corals.private.alleles.by.pop <- radiator::private_alleles(data = tidy, strata = strata.pop)
#' }



#' @rdname private_alleles
#' @export

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

private_alleles <- function(data, strata = NULL, verbose = TRUE) {
  if (verbose) {
    cat("################################################################################\n")
    cat("########################## radiator::private_alleles ###########################\n")
    cat("################################################################################\n")
  }


  # Cleanup---------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (verbose) message("Execution date/time: ", file.date)
  old.dir <- getwd()
  opt.change <- getOption("width")
  options(width = 70)
  timing <- proc.time()# for timing
  res <- list()
  #back to the original directory and options
  on.exit(setwd(old.dir), add = TRUE)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(timing <- proc.time() - timing, add = TRUE)
  on.exit(if (verbose) message("\nComputation time, overall: ", round(timing[[3]]), " sec"), add = TRUE)
  on.exit(if (verbose) cat("########################## completed private_alleles ###########################\n"), add = TRUE)

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("Input file missing")

  # Detect format --------------------------------------------------------------
  data.type <- radiator::detect_genomic_format(data)
  if (!data.type %in% c("tbl_df", "fst.file", "SeqVarGDSClass", "gds.file")) {
    rlang::abort("Input not supported for this function: read function documentation")
  }

  # Import data ---------------------------------------------------------------
  if (verbose) message("Importing data ...")
  if (data.type %in% c("tbl_df", "fst.file")) {
    if (is.vector(data)) data <- radiator::tidy_wide(data, import.metadata = TRUE)
    data.type <- "tbl_df"

    if (is.null(strata)) {
      if (rlang::has_name(data, "POP_ID") || rlang::has_name(data, "STRATA")) {
        if (rlang::has_name(data, "POP_ID")) data %<>% dplyr::rename(STRATA = POP_ID)
        strata <- generate_strata(data = data)
      } else {
        rlang::abort("strata required and if not provided a STRATA or POP_ID column in the data...")
      }
    } else {
      strata <- read_strata(strata = strata) %$% strata
      data %<>% join_strata(data = ., strata = strata)
    }

    if (rlang::has_name(data, "GT_VCF_NUC")) {
      private <- dplyr::filter(data, GT_VCF_NUC != "./.") %>%
        dplyr::distinct(STRATA, MARKERS, GT_VCF_NUC) %>%
        separate_gt(x = ., exclude = c("MARKERS", "STRATA")) %>%
        dplyr::select(-ALLELE_GROUP, ALLELE = HAPLOTYPES) %>%
        dplyr::distinct(MARKERS, STRATA, ALLELE)
    } else if (rlang::has_name(data, "GT")) {
      private <- dplyr::filter(data, GT != "000000") %>%
        dplyr::distinct(STRATA, MARKERS, GT) %>%
        separate_gt(x = ., sep = 3, gt = "GT",exclude = c("MARKERS", "STRATA")) %>%
        dplyr::select(-ALLELE_GROUP, ALLELE = HAPLOTYPES) %>%
        dplyr::distinct(MARKERS, STRATA, ALLELE)
    } else {
      # work with GT_BIN
      private <- dplyr::filter(data, !is.na(GT_BIN)) %>%
        dplyr::distinct(STRATA, MARKERS, GT_BIN) %>%
        dplyr::mutate(
          ALLELE = dplyr::case_when(
            GT_BIN == 0 ~ "REFREF",
            GT_BIN == 1 ~ "REFALT",
            GT_BIN == 2 ~ "ALTALT"
          ),
          GT_BIN = NULL
        ) %>%
        separate_gt(x = ., sep = 3, gt = "ALLELE", exclude = c("MARKERS", "STRATA")) %>%
        dplyr::select(-ALLELE_GROUP, ALLELE = HAPLOTYPES) %>%
        dplyr::distinct(MARKERS, STRATA, ALLELE)
    }
  }

  if (data.type %in% c("SeqVarGDSClass", "gds.file")) {
    if (!"SeqVarTools" %in% utils::installed.packages()[,"Package"]) {
      rlang::abort('Please install SeqVarTools for this option:\n
                   install.packages("BiocManager")
                   BiocManager::install("SeqVarTools")')
    }
    if (data.type == "gds.file") {
      data <- radiator::read_rad(data, verbose = verbose)
      data.type <- "SeqVarGDSClass"
    }
    strata <- extract_individuals(
      gds = data, ind.field.select = c("INDIVIDUALS", "STRATA")
    )
    private <- SeqVarTools::getGenotypeAlleles(
      gdsobj = data, use.names = TRUE) %>%
      magrittr::set_colnames(
        x = .,
        value = extract_markers_metadata(
          gds = data, markers.meta.select = "MARKERS"
        ) %$%
          MARKERS
      ) %>%
      tibble::as_tibble(x = ., rownames = "INDIVIDUALS") %>%
      data.table::as.data.table(.) %>%
      data.table::melt.data.table(
        data = .,
        id.vars = "INDIVIDUALS",
        variable.name = "MARKERS",
        value.name = "GT_VCF_NUC",
        variable.factor = FALSE) %>%
      tibble::as_tibble(.) %>%
      separate_gt(x = ., exclude = c("MARKERS", "INDIVIDUALS")) %>%
      dplyr::select(-ALLELE_GROUP, ALLELE = HAPLOTYPES)

    private %<>%
      join_strata(data = ., strata = strata, verbose = FALSE) %>%
      dplyr::distinct(MARKERS, STRATA, ALLELE)
  }

  data <- NULL
  private.search <- private %>%
    dplyr::group_by(MARKERS, ALLELE) %>%
    dplyr::tally(.) %>%
    dplyr::filter(n == 1) %>%
    dplyr::distinct(MARKERS, ALLELE) %>%
    dplyr::left_join(private, by = c("MARKERS", "ALLELE")) %>%
    dplyr::ungroup(.) %>%
    dplyr::select(STRATA, MARKERS, ALLELE) %>%
    dplyr::arrange(STRATA, MARKERS, ALLELE) %>%
    readr::write_tsv(x = ., path = "private.alleles.tsv")

  private.summary <- private.search %>%
    dplyr::group_by(STRATA) %>%
    dplyr::tally(.) %>%
    dplyr::ungroup(.) %>%
    tibble::add_row(.data = ., STRATA = "OVERALL", n = sum(.$n)) %>%
    dplyr::rename(PRIVATE_ALLELES = n) %>%
    readr::write_tsv(x = ., path = "private.alleles.summary.tsv")

  if(nrow(private.summary) > 0) {
    message("Number of private alleles per strata:")
    priv.message <- dplyr::mutate(private.summary, PRIVATE = stringi::stri_join(STRATA, PRIVATE_ALLELES, sep = " = "))
    message(stringi::stri_join(priv.message$PRIVATE, collapse = "\n"))
    res <- list(private.alleles = private.search, private.alleles.summary = private.summary)
  } else {
    message("Private allele(s) found: 0")
    res <- "0 private allele"
  }

  return(res)
}#End private_alleles
