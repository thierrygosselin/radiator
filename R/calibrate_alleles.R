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
#' Default: \code{biallelic = NULL}.

#' @param parallel.core (optional) The number of core used for parallel
#' execution. This is no longer used. The code is as fast as it can. Using
#' more cores will reduce the speed.
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

  # # test
  # biallelic = NULL
  # parallel.core = parallel::detectCores() - 1
  # verbose = TRUE
  # gt <- FALSE
  # gt.bin <- FALSE
  # gt.vcf <- FALSE
  # gt.vcf.nuc <- FALSE


  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("Input file missing")

  # dotslist -------------------------------------------------------------------
  rad.dots <- rlang::list2(...)
  gt <- rad.dots[["gt"]]
  gt.vcf.nuc <- rad.dots[["gt.vcf.nuc"]]
  gt.vcf <- rad.dots[["gt.vcf"]]
  gt.bin <- rad.dots[["gt.bin"]]
  if (is.null(gt.vcf.nuc)) gt.vcf.nuc <- FALSE
  if (is.null(gt.vcf)) gt.vcf <- FALSE
  if (is.null(gt)) gt <- FALSE
  if (is.null(gt.bin)) gt.bin <- FALSE



  # Detecting biallelic markers and removing monomorphic markers ---------------
  if (is.null(biallelic)) biallelic <- radiator::detect_biallelic_markers(data, verbose = FALSE)

  if (!biallelic) {
    message("Not implemented for multi-allelic data")
    return(list(input = data, biallelic = biallelic))
  }

  # check if the REF/ALT info is available
  ref.info <- rlang::has_name(data, "REF")
  # not possible to have that format if no REF/ALT info
  if (!ref.info) gt.vcf.nuc <- FALSE

  # Remove markers with REF = NA
  # note to myself : why not use detect_all_missing ? longer ?
  if (ref.info && anyNA(data$REF)) {
    all.missing <- dplyr::filter(data, is.na(REF)) %>% dplyr::distinct(MARKERS)
    readr::write_tsv(x = all.missing, file = "markers.missing.all.id.tsv")
    message("number of markers missing in all individuals and removed: ", nrow(all.missing))
    data %<>% dplyr::filter(!MARKERS %in% all.missing$MARKERS)
  }

  # original genotype format we have
  original.gt <- radiator::detect_gt(x = data, keep.one = FALSE)

  # No need to go through all the code if just GT...
  if (all(length(original.gt) == 1, original.gt == "GT")) {

    if (all(!gt.bin, !gt.vcf, !gt.vcf.nuc)) {
      return(list(input = data, biallelic = biallelic))
    } else {
      if (verbose) message("Generating new genotypes coding with calibrated alleles")
      data %<>%
        radiator::gt_recoding(
          x = .,
          gt = TRUE,
          gt.bin = gt.bin,
          gt.vcf = gt.vcf,
          gt.vcf.nuc = gt.vcf.nuc,
          arrange = FALSE
        )
      return(list(input = data, biallelic = biallelic))
    }
  }

  message("Calibrating REF/ALT alleles...")


  # what we need to make the switch
  detect.gt <- radiator::detect_gt(x = data, keep.one = TRUE, favorite = "GT_BIN")

  # format we want at the end is the original + the ones chosen in ...
  if ("GT" %in% original.gt) gt <- TRUE
  if ("GT_BIN" %in% original.gt) gt.bin <- TRUE
  if ("GT_VCF" %in% original.gt) gt.vcf <- TRUE
  if ("GT_VCF_NUC" %in% original.gt) gt.vcf.nuc <- TRUE


  if ("GT" %in% detect.gt) gt <- FALSE
  if ("GT_BIN" %in% detect.gt) gt.bin <- FALSE
  if ("GT_VCF" %in% detect.gt) gt.vcf <- FALSE
  if ("GT_VCF_NUC" %in% detect.gt) gt.vcf.nuc <- FALSE

  recoding.todo <- sum(c(gt, gt.bin, gt.vcf, gt.vcf.nuc), na.rm = TRUE)


  # check that REF is coded with letters
  # note to myself: not sure this is relevant now...
  # if (ref.info) {
  #   letter.coding <- unique(stringi::stri_detect_regex(
  #     str = unique(data$REF),
  #     pattern = "[A-Za-z]"))
  #
  #   if (letter.coding && !gt.vcf.nuc) {
  #     # data <- generate_vcf_nuc(x = data, parallel.core = parallel.core)
  #     data <- gt_recoding(x = data, gt = FALSE, gt.vcf = FALSE, arrange = FALSE) # faster...
  #     gt.vcf.nuc <- TRUE
  #   }
  #   if (!gt.vcf.nuc && letter.coding) ref.info <- FALSE
  # }



  # strip the data -------------------------------------------------------------
  # Not necessary here, faster without
  # strata.bk <- markers.meta.bk <- genotypes.meta.bk <- NULL
  # env.arg <- rlang::current_env()
  # data %<>%
  #   radiator::strip_rad(
  #     x = .,
  #     m = c("VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS", "COL"),
  #     env.arg = env.arg,
  #     keep.strata = FALSE,
  #     verbose = FALSE
  #   )

  if (!rlang::has_name(data, "MARKERS")) rlang::abort("Contact developer")

  # Generating REF/ALT dictionary  ---------------------------------------------
  n.markers <- length(unique(data$MARKERS))

  # GT_VCF_NUC
  if (detect.gt == "GT_VCF_NUC") {
    old.ref <- NULL
    if (ref.info) {
      old.ref <- dplyr::distinct(.data = data , MARKERS, REF, ALT) %>%
        dplyr::select(MARKERS, REF_OLD = REF)
    }

    switch.ref <- data %>%
      dplyr::select(tidyselect::any_of(c("MARKERS", "GT_VCF_NUC"))) %>%
      dplyr::filter(GT_VCF_NUC != "./.") %>%
      dplyr::count(MARKERS, GT_VCF_NUC) %>%
      radiator::separate_gt(
        x = .,
        gt = detect.gt,
        exclude = c("MARKERS", "n"),
        filter.missing = TRUE,
        split.chunks = 1L
      ) %>%
      dplyr::select(MARKERS, REF = ALLELES, n) %>%
      dplyr::group_by(MARKERS, REF) %>%
      dplyr::summarise(n = sum(n, na.rm = TRUE), .groups = "drop_last")

    if (nrow(switch.ref) / 2 != n.markers) {
      rlang::abort("Remove monomorphic markers in dataset")
    }
    # to make the code as precise as GT_BIN I need to account for
    # markers with equal counts for each alleles...

    switch.ref %<>%
      dplyr::filter(n == max(n, na.rm = TRUE)) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(EQUAL_COUNTS = dplyr::n()) %>%
      dplyr::ungroup(.) %>%
      dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
      dplyr::mutate(REF_OLD = old.ref$REF_OLD) %>%
      dplyr::filter(EQUAL_COUNTS == 1L) %>%
      dplyr::filter(REF != REF_OLD) %>%
      dplyr::distinct(MARKERS) %>%
      dplyr::pull()
    old.ref <- NULL
  }

  # GT_BIN
  if (detect.gt == "GT_BIN") {
    switch.ref <- data %>%
      dplyr::select(tidyselect::any_of(c("MARKERS", "GT_BIN"))) %>%
      dplyr::filter(!is.na(GT_BIN)) %>%
      dplyr::count(MARKERS, GT_BIN) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::summarise(
        REF = sum((2 * n[GT_BIN == 0L]), n[GT_BIN == 1L], na.rm = TRUE),
        ALT = sum(2 * n[GT_BIN == 2L], n[GT_BIN == 1L], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::filter(REF < ALT) %>%
      dplyr::distinct(MARKERS) %>%
      dplyr::pull()
  }

  # GT_VCF
  if (detect.gt == "GT_VCF") {
    switch.ref <- data %>%
      dplyr::select(tidyselect::any_of(c("MARKERS", "GT_VCF"))) %>%
      dplyr::filter(GT_VCF != "./.") %>%
      dplyr::count(MARKERS, GT_VCF) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::summarise(
        HET = sum(n[!GT_VCF %in% c("0/0", "1/1")], na.rm = TRUE),
        REF = sum(2 * n[GT_VCF == "0/0"], na.rm = TRUE) + HET,
        ALT = sum(2 * n[GT_VCF == "1/1"], na.rm = TRUE) + HET,
        .groups = "drop"
      ) %>%
      dplyr::filter(REF < ALT) %>%
      dplyr::distinct(MARKERS) %>%
      dplyr::pull()
  }

  if (is.factor(switch.ref)) switch.ref %<>% as.character(.)


  # Strategies -----------------------------------------------------------------
  # 1. GT and GT_VCF_NUC don't change...
  # 2. REF/ALT changes for those...
  # het from the other format don't change
  n.switch <- length(switch.ref)
  if (n.switch > 0) {
    if (verbose) message("Number of REF/ALT switch = ", n.switch)
    col.order <- colnames(data)
    remove.gt <- setdiff(c("GT", "GT_BIN", "GT_VCF", "GT_VCF_NUC"), detect.gt)
    ref.depth <- FALSE
    if (rlang::has_name(data, "ALLELE_REF_DEPTH")) ref.depth <- TRUE


    data <- dplyr::bind_rows(
      data %>%
        dplyr::select(-tidyselect::any_of(remove.gt)) %>%
        dplyr::filter(!MARKERS %in% switch.ref),
      data %>%
        dplyr::select(-tidyselect::any_of(remove.gt)) %>%
        dplyr::filter(MARKERS %in% switch.ref) %>%
        dplyr::rename(REF = ALT, ALT = REF) %>%
        {if (ref.depth) dplyr::rename(.data = ., ALLELE_REF_DEPTH = ALLELE_ALT_DEPTH, ALLELE_ALT_DEPTH = ALLELE_REF_DEPTH) else .} %>%
        dplyr::relocate(tidyselect::one_of(col.order)) %>%
        {if (detect.gt == "GT_BIN") dplyr::mutate(.data = .,GT_BIN = dplyr::recode(GT_BIN, '2' = 0L, '0' = 2L, .missing = NULL)) else .} %>%
        {if (detect.gt == "GT_VCF") dplyr::mutate(.data = .,GT_VCF = dplyr::recode(GT_VCF, "1/1" = "0/0", "0/0" = "1/1", .missing = NULL)) else .}
    )
  }
  if (recoding.todo > 0L) {
    if (verbose) message("Generating new genotypes coding")
    data %<>%
      radiator::gt_recoding(
        x = .,
        gt = gt,
        gt.bin = gt.bin,
        gt.vcf = gt.vcf,
        gt.vcf.nuc = gt.vcf.nuc,
        arrange = FALSE
      )
  }
  # Results --------------------------------------------------------------------
  return(list(input = data, biallelic = biallelic))
}#End calibrate_alleles

