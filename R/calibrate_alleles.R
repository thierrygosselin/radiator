# Compute the ref and alt alleles of a tidy dataset

#' @name calibrate_alleles

#' @title Calibrate REF and ALT alleles based on count

#' @description Calibrate REF and ALT alleles based on counts. The REF allele is
#' designated as the allele with more counts in the dataset.
#' The function will generate a REF and ALT columns.
#'
#' \strong{reference genome}: for people using a reference genome, the reference
#' allele terminology is different and is not based on counts...
#'
#'
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.
#'
#' @param data A genomic data set in the global environment tidy formats.
#' See details for more info.

#' @param biallelic (optional) If \code{biallelic = TRUE/FALSE} will be use
#' during multiallelic REF/ALT decision and speed up computations.
#' Used internally in
#' \href{https://github.com/thierrygosselin/radiator}{radiator}.
#' Default: \code{biallelic = NULL}.

#' @param parallel.core (optional) The number of core used for parallel
#' execution.
#' Default: \code{parallel.core = parallel::detectCores() - 1}.


#' @param verbose (optional, logical) \code{verbose = TRUE} to be chatty
#' during execution.
#' Default: \code{verbose = FALSE}.

#' @param ... (optional) To pass further argument for fine-tuning the tidying
#' (details below).


#' @return
#' Depending if the input file is biallelic or multiallelic,
#' the function will output additional to REF and ALT column several genotype codings:
#' \itemize{
#' \item \code{GT}: the genotype in 6 digits format with integers.
#' \item \code{GT_VCF}: the genotype in VCF format with integers.
#' \item \code{GT_VCF_NUC}: the genotype in VCF format with letters corresponding to nucleotide.
#' \item \code{GT_BIN}: biallelic coding similar to PLINK,
#' the coding \code{0, 1, 2, NA} correspond to the number of ALT allele in the
#' genotype and \code{NA} for missing genotypes.
#' }

#' @details \strong{Input data:}
#' A minimum of 4 columns are required (the rest are considered metata info):
#' \enumerate{
#' \item \code{MARKERS}
#' \item \code{POP_ID}
#' \item \code{INDIVIDUALS}
#' \item \code{GT} and/or \code{GT_VCF_NUC} and/or \code{GT_VCF}
#' }
#'
#' \emph{How to get a tidy data frame ?}
#' \pkg{radiator} \code{\link{tidy_genomic_data}}

#' @export
#' @rdname calibrate_alleles
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

calibrate_alleles <- function(
  data,
  biallelic = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = FALSE,
  ...
) {

  # note to myself, we could add genotypes meta in args... to leave only 1 GT column
  # integer prefereably...

  # # test
  # biallelic = NULL
  # parallel.core = parallel::detectCores() - 1
  # verbose = TRUE

  message("Calibrating REF/ALT alleles...")
  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("Input file missing")

  # dotslist -------------------------------------------------------------------
  dotslist <- list(...)# Not used for the moment

  # Detecting biallelic markers and removing monomorphic markers ---------------
  if (is.null(biallelic)) {
    biallelic <- radiator::detect_biallelic_markers(data)
    if (biallelic) {
      if (verbose) message("    Data is biallellic")
    } else {
      if (verbose) message("    Data is multiallellic")
    }
  }

  # Check if REF info is present -----------------------------------------------
  ref.info <- rlang::has_name(data, "REF")

  # Check if nucleotide info is available --------------------------------------
  nuc.info <- rlang::has_name(data, "GT_VCF_NUC")

  if (ref.info) {
    letter.coding <- unique(stringi::stri_detect_regex(
      str = unique(data$REF),
      pattern = "[A-Za-z]"))

    if (letter.coding && !nuc.info) {
      data <- generate_vcf_nuc(x = data, parallel.core = parallel.core)
      nuc.info <- TRUE
    }
    if (!nuc.info && letter.coding) ref.info <- FALSE
  }

  # strip the data -------------------------------------------------------------
  strata.bk <- markers.meta.bk <- genotypes.meta.bk <- NULL
  env.arg <- rlang::current_env()
  data %<>%
    radiator::strip_rad(
      x = .,
      m = c("VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS", "COL"),
      env.arg = env.arg,
      keep.strata = FALSE,
      verbose = FALSE
    )


  # Generating REF/ALT dictionary  ---------------------------------------------
  if (verbose) message("    generating REF/ALT dictionary")
  if (ref.info) {
    old.ref <- dplyr::distinct(.data = data , M_SEQ, REF, ALT) %>%
      dplyr::select(M_SEQ, REF_OLD = REF)
  } else {
    old.ref <- NULL
  }

  conversion.df <- ref_dictionary(x = data, parallel.core = parallel.core)


  # if monomorphic markers, ALT column will have NA: check and tag
  if (anyNA(conversion.df)) {
    # fill ALT with REF
    conversion.df <- dplyr::bind_rows(
      dplyr::filter(.data = conversion.df, is.na(ALT)) %>%
        dplyr::mutate(
          ALT = REF,
          POLYMORPHIC = rep(FALSE, n())),
      dplyr::filter(.data = conversion.df, !is.na(ALT)) %>%
        dplyr::mutate(POLYMORPHIC = rep(TRUE, n()))
    ) %>%
      dplyr::arrange(M_SEQ, INTEGERS)
  }
  new.ref <- conversion.df %>%
    dplyr::distinct(M_SEQ, REF, ALT, .keep_all = TRUE) %>%
    dplyr::select(-c(ALLELES, INTEGERS))
  # want <- c("M_SEQ", "ALLELES", "INTEGERS", "POLYMORPHIC")
  want <- c("M_SEQ", "ALLELES", "INTEGERS")
  conversion.df %<>% dplyr::select(tidyselect::any_of(want))

  # Inversion: keep track of change in REF/ALT -------------------------------
  if (!is.null(old.ref)) {
    change.ref <- dplyr::full_join(
      dplyr::select(new.ref, M_SEQ, REF),
      old.ref, by = "M_SEQ") %>%
      dplyr::filter(REF != REF_OLD)
    if (nrow(change.ref) > 0) {
      inversion <- TRUE
    } else {
      inversion <- FALSE
    }
    old.ref <- NULL
    message("    number of REF/ALT switch = ", nrow(change.ref))
  } else {
    inversion <- FALSE
  }

  if (rlang::has_name(data, "REF")) data %<>% dplyr::select(-c(REF, ALT))
  data %<>% dplyr::left_join(new.ref, by = "M_SEQ")
  new.ref <- NULL

  if (verbose) message("    integrating genotypes codings...")
  data %<>% dplyr::select(-tidyselect::any_of("POLYMORPHIC"))

  data <- integrate_ref(
    x = data,
    conversion.df = conversion.df,
    biallelic = biallelic,
    parallel.core = parallel.core)
  conversion.df <- NULL

  # For biallelic dataset with REF and ALT not coded with letters ----------------
  if (biallelic) {
    letter.coding <- unique(stringi::stri_detect_regex(
      str = unique(data$REF[!is.na(data$REF)]),
      pattern = "[A-Za-z]"))

    if (!letter.coding) {
      data %<>%
        dplyr::mutate(
          REF = stringi::stri_replace_all_regex(
            str = REF,
            pattern = c("001", "002", "003", "004", "1", "2"),
            replacement = c("A", "C", "G", "T", "A", "C"),
            vectorize_all = FALSE
          ),
          ALT = stringi::stri_replace_all_regex(
            str = ALT,
            pattern = c("001", "002", "003", "004", "1", "2"),
            replacement = c("A", "C", "G", "T", "A", "C"),
            vectorize_all = FALSE)
          )
    }
  }

  # Integrate the bk------------------------------------------------------------
  data <- radiator::join_rad(
    x = data,
    s = strata.bk,
    m = markers.meta.bk,
    g = genotypes.meta.bk,
    env.arg = env.arg
    )

  # switch ALLELE_REF_DEPTH/ALLELE_ALT_DEPTH
  if (inversion & rlang::has_name(data, "ALLELE_REF_DEPTH")) {
    data <- data %>%
      dplyr::mutate(
        ALLELE_REF_DEPTH_NEW = dplyr::if_else(
          MARKERS %in% change.ref, ALLELE_ALT_DEPTH, ALLELE_REF_DEPTH),
        ALLELE_ALT_DEPTH_NEW = dplyr::if_else(
          MARKERS %in% change.ref, ALLELE_REF_DEPTH, ALLELE_ALT_DEPTH)
      ) %>%
      dplyr::select(-ALLELE_REF_DEPTH, -ALLELE_ALT_DEPTH) %>%
      dplyr::rename(ALLELE_REF_DEPTH = ALLELE_REF_DEPTH_NEW,
                    ALLELE_ALT_DEPTH = ALLELE_ALT_DEPTH_NEW)
  }
  change.ref <- NULL

  #Remove markers with REF = NA ------------------------------------------------
  #note to myself : why not use detect_all_missing ? longer ?
  if (anyNA(data$REF)) {
    all.missing <- dplyr::filter(data, is.na(REF)) %>% dplyr::distinct(MARKERS)
    readr::write_tsv(x = all.missing, file = "markers.missing.all.id.tsv")
    message("    number of markers missing in all individuals and removed: ", nrow(all.missing))
    data <- dplyr::filter(data, !MARKERS %in% all.missing$MARKERS)
  }

  if (rlang::has_name(data, "GT_BIN")) data$GT_BIN <- as.integer(data$GT_BIN)
  # Results --------------------------------------------------------------------
  return(list(input = data, biallelic = biallelic))
}#End calibrate_alleles

# Internal nested Function -----------------------------------------------------

#' @title ref_dictionary
#' @description Generate REF/ALT allele dictionary from data with/without nucleotides.
#' @rdname ref_dictionary
#' @keywords internal
#' @export
ref_dictionary <- function(x, parallel.core = parallel::detectCores() - 1) {
  generate_ref <- carrier::crate(function(x) {
    `%>%` <- magrittr::`%>%`
    `%<>%` <- magrittr::`%<>%`
    # `n()` <- dplyr::`n`

    if (!any(c("GT_VCF_NUC", "GT", "GT_BIN") %in% colnames(x))) {
      rlang::abort(message = "Problem with genotypes coding: contact author")
    }

    x %<>%
      dplyr::select(tidyselect::any_of(c("M_SEQ", "GT_VCF_NUC", "GT", "GT_BIN")))

    skip <- FALSE
    if (rlang::has_name(x, "GT_VCF_NUC")) {# with nuc.info
      x %<>% dplyr::filter(GT_VCF_NUC != "./.")
      x %<>%
        dplyr::bind_cols(
          stringi::stri_split_fixed(str = x$GT_VCF_NUC, pattern = "/", simplify = TRUE) %>%
            magrittr::set_colnames(x = ., value = c("A1", "A2")) %>%
            tibble::as_tibble()
        ) %>%
        dplyr::select(-GT_VCF_NUC)

      # problem with code below is that it won't work with multi allelic VCF
      # A1 = stringi::stri_sub(str = GT_VCF_NUC, from = 1, to = 1)
      # A2 = stringi::stri_sub(str = GT_VCF_NUC, from = 3, to = 3)
      skip <- TRUE
    }

    if (rlang::has_name(x, "GT") && !skip) {
      x %<>%
        dplyr::filter(GT != "000000") %>%
        dplyr::mutate(
          A1 = stringi::stri_sub(str = GT, from = 1, to = 3),
          A2 = stringi::stri_sub(str = GT, from = 4, to = 6),
          GT = NULL
        )
      skip <- TRUE
    }

    if (rlang::has_name(x, "GT_BIN") && !skip) {
      x %<>%
        dplyr::filter(!is.na(GT_BIN)) %>%
        dplyr::mutate(
          A1 = dplyr::if_else(GT_BIN == 0L, 1L, GT_BIN),
          A2 = dplyr::if_else(GT_BIN != 2L, GT_BIN + 1L, GT_BIN),
          GT_BIN = NULL
        )
    }

    x %<>%
      radiator::rad_long(
        x = .,
        cols = "M_SEQ",
        names_to = "ALLELES_GROUP",
        values_to = "ALLELES",
        variable_factor = FALSE
      ) %>%
      dplyr::select(-ALLELES_GROUP) %>%
      dplyr::group_by(M_SEQ, ALLELES) %>%
      dplyr::tally(.) %>%
      dplyr::arrange(-n) %>%
      dplyr::mutate(INTEGERS = seq(0, dplyr::n() - 1)) %>%
      dplyr::select(-n) %>%
      dplyr::arrange(M_SEQ, INTEGERS) %>%
      dplyr::ungroup(.)

    ref.alt <- x %>%
      dplyr::filter(INTEGERS == 0) %>%
      dplyr::mutate(REF = as.character(ALLELES)) %>%
      dplyr::distinct(M_SEQ, REF) %>%
      dplyr::bind_cols(
        x %>%
          dplyr::filter(INTEGERS != 0) %>%
          dplyr::group_by(M_SEQ) %>%
          dplyr::mutate(ALT = stringi::stri_join(ALLELES, collapse = ",")) %>%
          dplyr::ungroup(.) %>%
          dplyr::distinct(M_SEQ, ALT) %>%
          dplyr::select(-M_SEQ)
      )

    x %<>%
      dplyr::left_join(ref.alt, by = "M_SEQ") %>%
      dplyr::arrange(M_SEQ, INTEGERS)
    return(x)
  })#End generate_ref

  n.markers <- length(unique(x$M_SEQ))

  if (n.markers <= 100000) {
    x <- generate_ref(x)
  } else {

    if (n.markers > 100000) split.chunks <- 2L
    if (n.markers > 200000) split.chunks <- 5L
    if (n.markers > 300000) split.chunks <- 10L

    x <- radiator_future(
      .x = x,
      .f = generate_ref,
      flat.future = "dfr",
      split.vec = TRUE,
      split.with = "M_SEQ",
      split.chunks = split.chunks,
      parallel.core = parallel.core
    )
  }
  return(x)
}

#' @title integrate_ref
#' @description Integrate ref allele in original dataset
#' @rdname integrate_ref
#' @keywords internal
#' @export
integrate_ref <- function(
  x,
  conversion.df = NULL,
  biallelic = TRUE,
  parallel.core = parallel::detectCores() - 1
) {
  # function needed
  new_gt <- carrier::crate(function(
    x,
    conversion.df = NULL,
    biallelic = TRUE,
    parallel.core = parallel::detectCores() - 1
  ) {
    `%>%` <- magrittr::`%>%`
    `%<>%` <- magrittr::`%<>%`
    # x <- new.gt #test
    nuc.info <- rlang::has_name(x, "GT_VCF_NUC")
    skip <- FALSE

    if (nuc.info) {
      x %<>% dplyr::select(M_SEQ, GT_VCF_NUC)
      x %<>%
        dplyr::bind_cols(
          stringi::stri_split_fixed(str = x$GT_VCF_NUC, pattern = "/", simplify = TRUE) %>%
            magrittr::set_colnames(x = ., value = c("A1", "A2")) %>%
            tibble::as_tibble()
        ) %>%
        dplyr::left_join(dplyr::rename(conversion.df, A1 = ALLELES), by = c("M_SEQ", "A1")) %>%
        dplyr::rename(A1_NUC = INTEGERS)
          # A1 = stringi::stri_sub(str = GT_VCF_NUC, from = 1, to = 1),
          # A2 = stringi::stri_sub(str = GT_VCF_NUC, from = 3, to = 3)
      # if (rlang::has_name(x, "POLYMORPHIC")) x %<>% dplyr::select(-POLYMORPHIC)

      x %<>%
        dplyr::left_join(dplyr::rename(conversion.df, A2 = ALLELES), by = c("M_SEQ", "A2")) %>%
        dplyr::rename(A2_NUC = INTEGERS) %>%
        dplyr::mutate(
          GT_VCF = stringi::stri_join(A1_NUC, A2_NUC, sep = "/"),
          GT_VCF = stringi::stri_replace_na(str = GT_VCF, replacement = "./."))

      if (biallelic) {
        x %<>%
          dplyr::select(-A1_NUC, -A2_NUC) %>%
          dplyr::mutate(
            A1 = stringi::stri_replace_all_regex(
              str = A1,
              pattern = c("^A$", "^C$", "^G$", "^T$", "^.$"),
              replacement = c("001", "002", "003", "004", "000"),
              vectorize_all = FALSE
            ),
            A2 = stringi::stri_replace_all_regex(
              str = A2,
              pattern = c("^A$", "^C$", "^G$", "^T$", "^.$"),
              replacement = c("001", "002", "003", "004", "000"),
              vectorize_all = FALSE
            ),
            GT_BIN = as.numeric(stringi::stri_replace_all_fixed(
              str = GT_VCF, pattern = c("0/0", "1/1", "0/1", "1/0", "./."),
              replacement = c("0", "2", "1", "1", NA), vectorize_all = FALSE))
          ) %>%
          tidyr::unite(data = ., col = GT, A1, A2, sep = "")

      } else {
        x %<>%
          dplyr::select(-A1, -A2) %>%
          dplyr::mutate(
            A1_NUC = as.character(A1_NUC + 1),
            A2_NUC = as.character(A2_NUC + 1),
            A1_NUC = stringi::stri_replace_na(str = A1_NUC, replacement = "0"),
            A2_NUC = stringi::stri_replace_na(str = A2_NUC, replacement = "0"),
            A1_NUC = stringi::stri_pad_left(str = A1_NUC, pad = "0", width = 3),
            A2_NUC = stringi::stri_pad_left(str = A2_NUC, pad = "0", width = 3)
          ) %>%
          tidyr::unite(data = ., col = GT, A1_NUC, A2_NUC, sep = "")
      }

      skip <- TRUE
    }
    if (rlang::has_name(x, "GT") && !skip) {
      x %<>%
        dplyr::select(M_SEQ, GT) %>%
        dplyr::mutate(
          A1 = stringi::stri_sub(str = GT, from = 1, to = 3),
          A2 = stringi::stri_sub(str = GT, from = 4, to = 6)
        ) %>%
        dplyr::select(-GT)

      # if (rlang::has_name(x, "POLYMORPHIC")) x %<>% dplyr::select(-POLYMORPHIC)

      x %<>%
        dplyr::left_join(dplyr::rename(conversion.df, A1 = ALLELES), by = c("M_SEQ", "A1")) %>%
        # dplyr::select(-dplyr::one_of("POLYMORPHIC")) %>%
        dplyr::rename(A1_NUC = INTEGERS) %>%
        dplyr::left_join(dplyr::rename(conversion.df, A2 = ALLELES), by = c("M_SEQ", "A2")) %>%
        # dplyr::select(-dplyr::one_of("POLYMORPHIC")) %>%
        dplyr::rename(A2_NUC = INTEGERS) %>%
        dplyr::mutate(
          GT_VCF = stringi::stri_join(A1_NUC, A2_NUC, sep = "/"),
          GT_VCF = stringi::stri_replace_na(str = GT_VCF, replacement = "./."),
          A1_NUC = as.character(A1_NUC + 1),
          A2_NUC = as.character(A2_NUC + 1),
          A1_NUC = stringi::stri_replace_na(str = A1_NUC, replacement = "0"),
          A2_NUC = stringi::stri_replace_na(str = A2_NUC, replacement = "0"),
          A1_NUC = stringi::stri_pad_left(str = A1_NUC, pad = "0", width = 3),
          A2_NUC = stringi::stri_pad_left(str = A2_NUC, pad = "0", width = 3)
        ) %>%
        tidyr::unite(data = ., col = GT, A1_NUC, A2_NUC, sep = "") %>%
        tidyr::unite(data = ., col = ORIG_GT, A1, A2, sep = "")

      if (biallelic) {
        x %<>%
          dplyr::mutate(
            GT_BIN = as.numeric(stringi::stri_replace_all_fixed(
              str = GT_VCF, pattern = c("0/0", "1/1", "0/1", "1/0", "./."),
              replacement = c("0", "2", "1", "1", NA), vectorize_all = FALSE)))
      }
      skip <- TRUE
    }
    if (rlang::has_name(x, "GT_BIN") && !skip) {

      # note to myself: must be a way to recycle that first part to speed up time
      x %<>%
        dplyr::select(M_SEQ, GT_BIN) %>%
        dplyr::mutate(
          A1 = dplyr::if_else(GT_BIN == 0L, 1L, GT_BIN),
          A2 = dplyr::if_else(GT_BIN != 2L, GT_BIN + 1L, GT_BIN)
          # A1 = stringi::stri_pad_left(str = A1, width = 3, pad = "0"),
          # A2 = stringi::stri_pad_left(str = A2, width = 3, pad = "0"),
        ) %>%
        dplyr::rename(ORIG_GT_BIN = GT_BIN)


      x %<>%
        dplyr::left_join(conversion.df, by = c("M_SEQ", "A1" = "ALLELES")) %>%
        dplyr::rename(A1_NUC = INTEGERS) %>%
        dplyr::left_join(conversion.df, by = c("M_SEQ", "A2" = "ALLELES")) %>%
        dplyr::rename(A2_NUC = INTEGERS) %>%
        dplyr::mutate(
          A1 = NULL, A2 = NULL,
          GT_VCF = stringi::stri_join(A1_NUC, A2_NUC, sep = "/"),
          GT_VCF = stringi::stri_replace_na(str = GT_VCF, replacement = "./."),
          A1_NUC = as.character(A1_NUC + 1),
          A2_NUC = as.character(A2_NUC + 1),
          A1_NUC = stringi::stri_replace_na(str = A1_NUC, replacement = "0"),
          A2_NUC = stringi::stri_replace_na(str = A2_NUC, replacement = "0"),
          A1_NUC = stringi::stri_pad_left(str = A1_NUC, pad = "0", width = 3),
          A2_NUC = stringi::stri_pad_left(str = A2_NUC, pad = "0", width = 3)
        ) %>%
        tidyr::unite(data = ., col = GT, A1_NUC, A2_NUC, sep = "")

      if (biallelic) {
        x %<>%
          dplyr::mutate(
            GT_BIN = as.numeric(stringi::stri_replace_all_fixed(
              str = GT_VCF, pattern = c("0/0", "1/1", "0/1", "1/0", "./."),
              replacement = c("0", "2", "1", "1", NA), vectorize_all = FALSE)))
      }
    }

    return(x)
  })#End new_gt

  nuc.info <- rlang::has_name(x, "GT_VCF_NUC")
  if (nuc.info) new.gt <- dplyr::distinct(x, M_SEQ, GT_VCF_NUC)
  if (rlang::has_name(x, "GT")) new.gt <- dplyr::distinct(x, M_SEQ, GT)
  if (rlang::has_name(x, "GT_BIN")) new.gt <- dplyr::distinct(x, M_SEQ, GT_BIN)


  n.markers <- length(unique(new.gt$M_SEQ))

  if (n.markers <= 100000) {
    new.gt %<>% new_gt(x = ., conversion.df = conversion.df, biallelic = biallelic)
  } else {

    if (n.markers > 100000) split.chunks <- 2L
    if (n.markers > 200000) split.chunks <- 5L
    if (n.markers > 300000) split.chunks <- 10L

    new.gt <- radiator_future(
      .x = new.gt,
      .f = new_gt,
      flat.future = "dfr",
      split.vec = TRUE,
      split.with = NULL,
      split.chunks = split.chunks,
      parallel.core = parallel.core,
      conversion.df = conversion.df,
      biallelic = biallelic
    )
  }

  want <- c("GT", "GT_VCF_NUC", "GT_VCF", "GT_BIN")
  gt.format <- purrr::keep(.x = want, .p = want %in% colnames(x))

  if (gt.format == "GT_VCF_NUC") {
    notwanted <- c("GT_VCF", "GT_BIN", "GT")
    x %<>%
      dplyr::select(-tidyselect::any_of(notwanted)) %>%
      dplyr::left_join(new.gt, by = c("M_SEQ", "GT_VCF_NUC"))
  }

  if (gt.format == "GT") {
    notwanted <- c("GT_VCF", "GT_BIN")
    x %<>%
      dplyr::select(-tidyselect::any_of(notwanted)) %>%
      dplyr::rename(ORIG_GT = GT) %>%
      dplyr::left_join(new.gt, by = c("M_SEQ", "GT" = "ORIG_GT"))
  }

  if (gt.format == "GT_BIN") {
    notwanted <- c("GT_VCF", "GT")
    x %<>%
      dplyr::select(-tidyselect::any_of(notwanted)) %>%
      dplyr::rename(ORIG_GT_BIN = GT_BIN) %>%
      dplyr::left_join(new.gt, by = c("M_SEQ", "ORIG_GT_BIN")) %>%
      dplyr::select(-GT_BIN) %>%
      dplyr::rename(GT_BIN = ORIG_GT_BIN)
}


  new.gt <- NULL
  return(x)
} #End integrate_ref



#' @title generate_vcf_nuc
#' @description Generate GT_VCF_NUC field from REF/ALT and GT_BIN columns
#' @rdname generate_vcf_nuc
#' @keywords internal
#' @export

generate_vcf_nuc <- function(x, parallel.core = parallel::detectCores() - 1) {
  vcf_nuc <- carrier::crate(function(x) {
    `%>%` <- magrittr::`%>%`
    `%<>%` <- magrittr::`%<>%`

    x %<>%
      dplyr::mutate(
        GT_VCF_NUC = dplyr::case_when(
          GT_BIN == "0" ~ stringi::stri_join(REF, REF, sep = "/"),
          GT_BIN == "1" ~ stringi::stri_join(REF, ALT, sep = "/"),
          GT_BIN == "2" ~ stringi::stri_join(ALT, ALT, sep = "/"),
          is.na(GT_BIN) ~ "./."
        ),
        GT_BIN = NULL
      )
  })# End vcf_nuc

  n.markers <- length(unique(x$MARKERS))

  if (n.markers <= 100000) {
    x <- vcf_nuc(x)
  } else {

    if (n.markers > 100000) split.chunks <- 2L
    if (n.markers > 200000) split.chunks <- 5L
    if (n.markers > 300000) split.chunks <- 10L

    x <- radiator_future(
      .x = x,
      .f = vcf_nuc,
      flat.future = "dfr",
      split.vec = TRUE,
      split.with = NULL,
      split.chunks = split.chunks,
      parallel.core = parallel.core
    )
  }


  return(x)
}#End generate_vcf_nuc
