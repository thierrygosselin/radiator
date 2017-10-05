# Compute the ref and alt alleles of a tidy dataset

#' @name change_alleles

#' @title Change REF and ALT alleles based on count

#' @description Change REF and ALT alleles based on count.
#' The function will generate a REF and ALT columns.
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
#' @rdname change_alleles
#' @importFrom dplyr select mutate group_by ungroup rename tally filter if_else arrange summarise top_n distinct coalesce if_else full_join
#' @importFrom stringi stri_replace_all_fixed stri_join stri_sub
#' @importFrom data.table fread
#' @importFrom tibble has_name
#' @importFrom tidyr spread gather
#' @importFrom purrr flatten_chr

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

change_alleles <- function(
  data,
  biallelic = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = FALSE) {

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file missing")

  # check genotype column naming (to work with old code)
  if (tibble::has_name(data, "GENOTYPE")) {
    colnames(data) <- stringi::stri_replace_all_fixed(
      str = colnames(data),
      pattern = "GENOTYPE",
      replacement = "GT",
      vectorize_all = FALSE)
  }

  # necessary steps to make sure we work with unique markers and not duplicated LOCUS
  if (tibble::has_name(data, "LOCUS") && !tibble::has_name(data, "MARKERS")) {
    data <- dplyr::rename(.data = data, MARKERS = LOCUS)
  }

  # get number of markers
  n.catalog.locus <- dplyr::n_distinct(data$MARKERS)

  # Check if REF info is present -----------------------------------------------
  ref.info <- tibble::has_name(data, "REF")

  # Check if nucleotide info is available --------------------------------------
  nuc.info <- tibble::has_name(data, "GT_VCF_NUC")
  if (ref.info) {
    letter.coding <- unique(stringi::stri_detect_regex(
      str = unique(data$REF),
      pattern = "[A-Za-z]"))
    if (!nuc.info && letter.coding) ref.info <- FALSE
  }
  # Check for markers meta info ------------------------------------------------
  want <- c("MARKERS", "CHROM", "LOCUS", "POS")
  markers.meta <- suppressWarnings(
    dplyr::select(data, dplyr::one_of(want)) %>%
      dplyr::distinct(MARKERS, .keep_all = TRUE))
  if (ncol(markers.meta) == 1) markers.meta <- NULL

  # Detecting biallelic markers and removing monomorphic markers ---------------
  if (is.null(biallelic)) {
    biallelic <- radiator::detect_biallelic_markers(data)
    if (biallelic) {
      if (verbose) message("    Data is biallellic")
    } else {
      if (verbose) message("    Data is multiallellic")
    }
  }

  # Generating REF/ALT dictionary  ---------------------------------------------
  if (verbose) message("    Generating REF/ALT dictionary")
  if (ref.info) {
    old.ref <- data %>%
      dplyr::distinct(MARKERS, REF, ALT) %>%
      dplyr::select(MARKERS, REF_OLD = REF)
  } else {
    old.ref <- NULL
  }

  conversion.df <- ref_dictionary(
    x = data, parallel.core = parallel.core)

  # if monomorphic markers, ALT column will have NA: check and tag
  if (anyNA(conversion.df)) {
    # fill ALT with REF
    conversion.df <- dplyr::bind_rows(
      dplyr::filter(conversion.df, is.na(ALT)) %>%
        dplyr::mutate(
          ALT = REF,
          POLYMORPHIC = rep(FALSE, n())),
      dplyr::filter(conversion.df, !is.na(ALT)) %>%
        dplyr::mutate(POLYMORPHIC = rep(TRUE, n()))) %>%
      dplyr::arrange(MARKERS, INTEGERS)
  }
  new.ref <- conversion.df %>%
    dplyr::distinct(MARKERS, REF, ALT, .keep_all = TRUE) %>%
    dplyr::select(-c(ALLELES, INTEGERS))
  want <- c("MARKERS", "ALLELES", "INTEGERS", "POLYMORPHIC")
  conversion.df <- suppressWarnings(dplyr::select(conversion.df, dplyr::one_of(want)))

  # Inversion: keep track of change in REF/ALT -------------------------------
  if (!is.null(old.ref)) {
    change.ref <- dplyr::full_join(
      dplyr::select(new.ref, MARKERS, REF),
      old.ref, by = "MARKERS") %>%
      dplyr::filter(REF != REF_OLD)
    inversion <- TRUE # not yet used
    old.ref <- NULL
    message("    Number of markers with REF/ALT change(s) = ", nrow(change.ref))
  } else {
    inversion <- FALSE # not yet used
    change.ref <- NULL
  }

  if (tibble::has_name(data, "REF")) {
    # if (tibble::has_name(data, "REF") && inversion) {
    data <- dplyr::select(data, -c(REF, ALT))
  }
  data <- dplyr::left_join(data, new.ref, by = "MARKERS")
  new.ref <- NULL

  if (verbose) message("    Integrating new genotype codings...")
  data <- integrate_ref(
    x = data,
    conversion.data = conversion.df,
    biallelic = biallelic,
    parallel.core = parallel.core)
  conversion.df <- NULL

  # switch ALLELE_REF_DEPTH/ALLELE_ALT_DEPTH
  if (!is.null(change.ref) & tibble::has_name(data, "ALLELE_REF_DEPTH")) {
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

  if (!is.null(markers.meta)) {
    no.need.markers.meta <- unique(colnames(markers.meta) %in% colnames(data))
    if (length(no.need.markers.meta) == 1 && !no.need.markers.meta) {
      data <- dplyr::left_join(data, markers.meta, by = "MARKERS") %>%
        dplyr::select(MARKERS, CHROM, LOCUS, POS, POP_ID, INDIVIDUALS, dplyr::everything())
    }
    markers.meta <- NULL
  } else {
    data <- dplyr::select(
      data,
      MARKERS, POP_ID, INDIVIDUALS, dplyr::everything())
  }
  res <- list(input = data, biallelic = biallelic)
  return(res)
}#End change_alleles

# Internal nested Function -----------------------------------------------------

#' @title ref_dictionary
#' @description Generate REF/ALT allele dictionary from data with/without nucleotides.
#' @rdname ref_dictionary
#' @keywords internal
#' @export
ref_dictionary <- function(x, parallel.core = parallel::detectCores() - 1) {
  generate_ref <- function(x) {
    nuc.info <- tibble::has_name(x, "GT_VCF_NUC")
    if (nuc.info) {
      res <- dplyr::select(x, MARKERS, GT_VCF_NUC) %>%
        dplyr::filter(GT_VCF_NUC != "./.") %>%
        tidyr::separate(col = GT_VCF_NUC, into = c("A1", "A2"), sep = "/") %>%
        tidyr::gather(data = ., key = ALLELES_GROUP, value = ALLELES, -MARKERS) %>%
        dplyr::select(-ALLELES_GROUP) %>%
        dplyr::group_by(MARKERS, ALLELES) %>%
        dplyr::tally(.) %>%
        dplyr::arrange(-n) %>%
        dplyr::mutate(INTEGERS = seq(0, n() - 1)) %>%
        dplyr::select(-n) %>%
        dplyr::arrange(MARKERS, INTEGERS) %>%
        dplyr::ungroup(.)
    } else {
      res <- suppressWarnings(
        dplyr::select(x, MARKERS, GT) %>%
          dplyr::filter(GT != "000000") %>%
          dplyr::mutate(
            A1 = stringi::stri_sub(str = GT, from = 1, to = 3),
            A2 = stringi::stri_sub(str = GT, from = 4, to = 6)
          ) %>%
          dplyr::select(-GT) %>%
          tidyr::gather(
            data = .,
            key = ALLELES_GROUP, value = ALLELES, -dplyr::one_of(names(x))) %>%
          dplyr::select(MARKERS, ALLELES) %>%
          dplyr::filter(ALLELES != "000") %>%
          dplyr::group_by(MARKERS, ALLELES) %>%
          dplyr::tally(.) %>%
          dplyr::arrange(-n) %>%
          dplyr::mutate(INTEGERS = seq(0, n() - 1)) %>%
          dplyr::select(-n) %>%
          dplyr::arrange(MARKERS, INTEGERS) %>%
          dplyr::ungroup(.))
    }
    ref.alleles <- res %>%
      dplyr::filter(INTEGERS == 0) %>%
      dplyr::mutate(REF = ALLELES) %>%
      dplyr::distinct(MARKERS, REF)

    alt.alleles <- res %>%
      dplyr::filter(INTEGERS != 0) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(ALT = stringi::stri_join(ALLELES, collapse = ",")) %>%
      dplyr::ungroup(.) %>%
      dplyr::distinct(MARKERS, ALT)

    res <- dplyr::left_join(
      res,
      dplyr::left_join(ref.alleles, alt.alleles, by = "MARKERS")
      , by = "MARKERS"
    ) %>%
      dplyr::arrange(MARKERS, INTEGERS)
    return(res)
  }#End generate_ref

  res <- x %>%
    dplyr::left_join(
      dplyr::distinct(x, MARKERS) %>%
        dplyr::mutate(
          SPLIT_VEC = dplyr::ntile(x = 1:nrow(.), n = parallel.core * 3))
      , by = "MARKERS") %>%
    split(x = ., f = .$SPLIT_VEC) %>%
    .radiator_parallel(
      X = .,
      FUN = generate_ref,
      mc.cores = parallel.core
    ) %>%
    dplyr::bind_rows(.)
  return(res)
}

#' @title integrate_ref
#' @description Integrate ref allele in original dataset
#' @rdname integrate_ref
#' @keywords internal
#' @export
integrate_ref <- function(
  x,
  conversion.data = NULL,
  biallelic = TRUE,
  parallel.core = parallel::detectCores() - 1
) {
  # function needed
  new_gt <- function(
    x,
    conversion.data = NULL,
    biallelic = TRUE,
    parallel.core = parallel::detectCores() - 1
  ) {
    nuc.info <- tibble::has_name(x, "GT_VCF_NUC")

    if (nuc.info) {
      res <- dplyr::select(x, MARKERS, GT_VCF_NUC) %>%
        tidyr::separate(data = ., col = GT_VCF_NUC, into = c("A1", "A2"), sep = "/", remove = FALSE) %>%
        dplyr::left_join(dplyr::rename(conversion.data, A1 = ALLELES), by = c("MARKERS", "A1")) %>%
        dplyr::rename(A1_NUC = INTEGERS) %>%
        dplyr::left_join(dplyr::rename(conversion.data, A2 = ALLELES), by = c("MARKERS", "A2")) %>%
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
        dplyr::select(-A1, -A2)

      if (biallelic) {
        res <- res %>%
          dplyr::mutate(
            GT_BIN = as.numeric(stringi::stri_replace_all_fixed(
              str = GT_VCF, pattern = c("0/0", "1/1", "0/1", "1/0", "./."),
              replacement = c("0", "2", "1", "1", NA), vectorize_all = FALSE)))
      }
    } else {
      res <- dplyr::select(x, MARKERS, GT) %>%
        dplyr::mutate(
          A1 = stringi::stri_sub(str = GT, from = 1, to = 3),
          A2 = stringi::stri_sub(str = GT, from = 4, to = 6)
        ) %>%
        dplyr::select(-GT) %>%
        dplyr::left_join(dplyr::rename(conversion.data, A1 = ALLELES), by = c("MARKERS", "A1")) %>%
        dplyr::rename(A1_NUC = INTEGERS) %>%
        dplyr::left_join(dplyr::rename(conversion.data, A2 = ALLELES), by = c("MARKERS", "A2")) %>%
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
        res <- res %>%
          dplyr::mutate(
            GT_BIN = as.numeric(stringi::stri_replace_all_fixed(
              str = GT_VCF, pattern = c("0/0", "1/1", "0/1", "1/0", "./."),
              replacement = c("0", "2", "1", "1", NA), vectorize_all = FALSE)))
      }
    }
    return(res)
  }#End new_gt

  nuc.info <- tibble::has_name(x, "GT_VCF_NUC")
  if (nuc.info) {
    new.gt <- dplyr::distinct(x, MARKERS, GT_VCF_NUC)
  } else {
    new.gt <- dplyr::distinct(x, MARKERS, GT)
  }
  new.gt <- new.gt %>%
    dplyr::mutate(SPLIT_VEC = dplyr::ntile(x = 1:nrow(.), n = parallel.core * 3)) %>%
    split(x = ., f = .$SPLIT_VEC) %>%
    .radiator_parallel(
      X = .,
      FUN = new_gt,
      mc.cores = parallel.core,
      conversion.data = conversion.data,
      biallelic = biallelic
    ) %>%
    dplyr::bind_rows(.)

  if (nuc.info) {
    x <- dplyr::left_join(x, new.gt, by = c("MARKERS", "GT_VCF_NUC")) %>%
      dplyr::select(-ORIG_GT)
  } else {
    if (tibble::has_name(x, "GT_VCF")) x <- dplyr::select(x, -GT_VCF)
    x <- dplyr::left_join(rename(x, ORIG_GT = GT), new.gt, by = c("MARKERS", "ORIG_GT")) %>%
      dplyr::select(-GT) %>%
      dplyr::rename(GT = ORIG_GT)
  }
  return(x)
} #End integrate_ref

