# Compute the ref and alt alleles of a tidy dataset

#' @name change_alleles

#' @title Change REF and ALT alleles based on count

#' @description Change REF and ALT alleles based on count, for biallelic data
#' with REF or ALT info, the function will generate a REF and ALT columns.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.
#'
#' @param data A biallelic genomic data set in the working directory or
#' object in the global environment in wide or long (tidy) formats.
#' See details for more info.

#' @param monomorphic.out (optional) Should the monomorphic
#' markers present in the dataset be filtered out ?
#' Default: \code{monomorphic.out = TRUE}.

#' @param biallelic (optional) If \code{biallelic = TRUE/FALSE} will be use
#' during multiallelic REF/ALT decision. Used internally in
#' \href{https://github.com/thierrygosselin/radiator}{radiator}.
#' Default: \code{biallelic = NULL}.

#' @param parallel.core (optional) The number of core used for parallel
#' execution.
#' Default: \code{parallel::detectCores() - 1}.


#' @param verbose (optional, logical) \code{verbose = TRUE} to be chatty
#' during execution.
#' Default: \code{verbose = FALSE}.

#' @return A tidy data frame with 6 columns:
#' \code{MARKERS, INDIVIDUALS, REF, ALT, GT_VCF, GT_BIN}.
#' \code{GT_VCF}: the genotype in VCF format
#' \code{GT_BIN}: coding used internally to easily convert to genlight,
#' the coding \code{0, 1, 2, NA} stands for the number of ALT allele in the
#' genotype and \code{NA} for missing genotype.

#' @details \strong{Input data:}
#'
#' To discriminate the long from the wide format,
#' the function \pkg{radiator} \code{\link[radiator]{tidy_wide}} searches
#' for \code{MARKERS or LOCUS} in column names (TRUE = long format).
#' The data frame is tab delimitted.

#' \strong{Wide format:}
#' The wide format cannot store metadata info.
#' The wide format starts with these 2 id columns:
#' \code{INDIVIDUALS}, \code{POP_ID} (that refers to any grouping of individuals),
#' the remaining columns are the markers in separate columns storing genotypes.
#'
#' \strong{Long/Tidy format:}
#' The long format is considered to be a tidy data frame and can store metadata info.
#' (e.g. from a VCF see \pkg{radiator} \code{\link{tidy_genomic_data}}). A minimum of 4 columns
#' are required in the long format: \code{INDIVIDUALS}, \code{POP_ID},
#' \code{MARKERS or LOCUS} and \code{GENOTYPE or GT}. The rest are considered metata info.
#'
#' \strong{2 genotypes formats are available:}
#' 6 characters no separator: e.g. \code{001002 of 111333} (for heterozygote individual).
#' 6 characters WITH separator: e.g. \code{001/002 of 111/333} (for heterozygote individual).
#' The separator can be any of these: \code{"/", ":", "_", "-", "."}.
#'
#' \emph{How to get a tidy data frame ?}
#' \pkg{radiator} \code{\link{tidy_genomic_data}} can transform 6 genomic data formats
#' in a tidy data frame.
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
  data, monomorphic.out = TRUE, biallelic = NULL,
  parallel.core = parallel::detectCores() - 1, verbose = FALSE) {

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file missing")

  input <- data

  # check genotype column naming
  if (tibble::has_name(input, "GENOTYPE")) {
    colnames(input) <- stringi::stri_replace_all_fixed(
      str = colnames(input),
      pattern = "GENOTYPE",
      replacement = "GT",
      vectorize_all = FALSE)
  }

  # necessary steps to make sure we work with unique markers and not duplicated LOCUS
  if (tibble::has_name(input, "LOCUS") && !tibble::has_name(input, "MARKERS")) {
    input <- dplyr::rename(.data = input, MARKERS = LOCUS)
  }

  if (!tibble::has_name(input, "GT_VCF_NUC")) {
    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "INDIVIDUALS", "POP_ID",
              "REF", "ALT", "GT", "GT_HAPLO")

    input <- suppressWarnings(
      dplyr::select(input, dplyr::one_of(want)))
    need.gt.vcf <- TRUE
  } else {
    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "INDIVIDUALS", "POP_ID",
              "REF", "ALT", "GT_VCF_NUC")

    input <- suppressWarnings(
      dplyr::select(input, dplyr::one_of(want)))

    need.gt.vcf <- FALSE
    biallelic <- TRUE
  }

  # Detecting biallelic markers and removing monomorphic markers ---------------
  if (verbose) message("    Scanning for number of alleles per marker...")
  if (tibble::has_name(input, "GT")) {
    input.genotyped.split <- dplyr::select(.data = input, MARKERS, POP_ID, INDIVIDUALS, GT) %>%
      dplyr::filter(GT != "000000") %>%
      dplyr::mutate(
        A1 = stringi::stri_sub(str = GT, from = 1, to = 3),
        A2 = stringi::stri_sub(str = GT, from = 4, to = 6)
      ) %>%
      dplyr::select(-GT) %>%
      tidyr::gather(data = ., key = ALLELES_GROUP, value = ALLELES, -c(MARKERS, INDIVIDUALS, POP_ID))

    alleles.old <- dplyr::distinct(.data = input.genotyped.split, MARKERS, ALLELES) %>%
      dplyr::arrange(MARKERS, ALLELES)

    marker.type <- dplyr::distinct(.data = input.genotyped.split, MARKERS, ALLELES) %>%
      dplyr::count(x = ., MARKERS)

    # monomorphic
    if (monomorphic.out) {
      mono.markers <-  dplyr::filter(.data = marker.type, n == 1) %>%
        dplyr::select(MARKERS)
    }

    # Biallelic marker detection -------------------------------------------------
    biallelic <- unique(marker.type$n)
    # if (length(biallelic) > 4) stop("Mix of bi- and multi-allelic markers is not supported")
    biallelic <- purrr::flatten_dbl(.x = dplyr::summarise(.data = marker.type, BIALLELIC = max(n, na.rm = TRUE)))

    if (biallelic > 3) {
      biallelic <- FALSE
      if (verbose) message("    Data is multiallellic")
      if (!tibble::has_name(input, "GT_HAPLO")) {
        input <- dplyr::rename(input, GT_HAPLO = GT)
      }
    } else {
      biallelic <- TRUE
      if (verbose) message("    Data is biallellic")
    }

    marker.type <- NULL
  } # end GT prep

  if (tibble::has_name(input, "GT_HAPLO")) biallelic <- FALSE

  # Function to calculate REF\ALT --------------------------------------------
  ref_compute <- function(data, new.ref) {
    input <- data %>%
      dplyr::mutate(
        A1 = stringi::stri_sub(str = GT, from = 1, to = 3),
        A2 = stringi::stri_sub(str = GT, from = 4, to = 6)
      ) %>%
      dplyr::full_join(new.ref, by = "MARKERS") %>%
      dplyr::mutate(
        A1 = replace(A1, which(A1 == "000"), NA),
        A2 = replace(A2, which(A2 == "000"), NA),
        GT_VCF_A1 = dplyr::if_else(A1 == REF, "0", "1", missing = "."),
        GT_VCF_A2 = dplyr::if_else(A2 == REF, "0", "1", missing = "."),
        GT_VCF = stringi::stri_join(GT_VCF_A1, GT_VCF_A2, sep = "/"),
        GT_BIN = as.numeric(stringi::stri_replace_all_fixed(
          str = GT_VCF,
          pattern = c("0/0", "1/1", "0/1", "1/0", "./."),
          replacement = c("0", "2", "1", "1", NA),
          vectorize_all = FALSE
        )),
        REF = stringi::stri_replace_all_fixed(
          str = REF,
          pattern = c("001", "002", "003", "004"),
          replacement = c("A", "C", "G", "T"),
          vectorize_all = FALSE),
        ALT = stringi::stri_replace_all_fixed(
          str = ALT,
          pattern = c("001", "002", "003", "004"),
          replacement = c("A", "C", "G", "T"),
          vectorize_all = FALSE)
      ) %>%
      dplyr::select(-c(A1, A2, GT_VCF_A1, GT_VCF_A2))
  }#End ref_compute


  # Detection and change -------------------------------------------------------
  if (tibble::has_name(input, "MARKERS") && tibble::has_name(input, "CHROM") && tibble::has_name(input, "LOCUS") && tibble::has_name(input, "POS")) {
    markers.meta <- dplyr::distinct(input, MARKERS, CHROM, LOCUS, POS)
  } else {
    markers.meta <- NULL
  }

  if (tibble::has_name(input, "GT")) {
    if (verbose) message("    Generating vcf-style coding")
    alleles.new.ref <- dplyr::select(.data = input.genotyped.split, MARKERS, ALLELES) %>%
      dplyr::count(x = ., MARKERS, ALLELES) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::top_n(1, n) %>%
      dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
      dplyr::mutate(REF = rep("REF", n())) %>%
      dplyr::select(-n) %>%
      dplyr::full_join(alleles.old, by = c("MARKERS", "ALLELES")) %>%
      dplyr::mutate(REF = dplyr::coalesce(REF, "ALT")) %>% # faster than stri_replace_na
      dplyr::group_by(MARKERS) %>%
      tidyr::spread(data = ., key = REF, value = ALLELES) %>%
      dplyr::mutate(ALT = dplyr::if_else(is.na(ALT), REF, ALT))

    alleles.old <- NULL

    if (tibble::has_name(input, "REF")) {
      change.ref <- dplyr::distinct(.data = input, MARKERS, REF, ALT) %>%
        dplyr::full_join(
          alleles.new.ref %>%
            dplyr::select(MARKERS, REF_NEW = REF)
          , by = "MARKERS") %>%
        dplyr::mutate(
          REF_NEW = stringi::stri_replace_all_fixed(
            str = REF_NEW,
            pattern = c("001", "002", "003", "004"),
            replacement = c("A", "C", "G", "T"),
            vectorize_all = FALSE),
          CHANGE = dplyr::if_else(REF == REF_NEW, "identical", "different")
        ) %>%
        dplyr::filter(CHANGE == "different") %>%
        dplyr::select(MARKERS) %>%
        purrr::flatten_chr(.)

      message("    Number of markers with REF/ALT change = ", length(change.ref))

      # switch ALLELE_REF_DEPTH/ALLELE_ALT_DEPTH
      if (length(change.ref) > 0 & tibble::has_name(input, "ALLELE_REF_DEPTH")) {
        input <- input %>%
          dplyr::mutate(
            ALLELE_REF_DEPTH_NEW = dplyr::if_else(MARKERS %in% change.ref, ALLELE_ALT_DEPTH, ALLELE_REF_DEPTH),
            ALLELE_ALT_DEPTH_NEW = dplyr::if_else(MARKERS %in% change.ref, ALLELE_REF_DEPTH, ALLELE_ALT_DEPTH)
          ) %>%
          dplyr::select(-ALLELE_REF_DEPTH, -ALLELE_ALT_DEPTH) %>%
          dplyr::rename(ALLELE_REF_DEPTH = ALLELE_REF_DEPTH_NEW, ALLELE_ALT_DEPTH = ALLELE_ALT_DEPTH_NEW)
      }

      # switch REF/ALT in the dataset
      if (length(change.ref) > 0) {
        input <- dplyr::select(input, -c(REF, ALT))
        input <- ref_compute(data = input, new.ref = alleles.new.ref)
      }

    } else {
      input <- ref_compute(data = input, new.ref = alleles.new.ref)
    }

    # monomorphic filter
    if (monomorphic.out) {
      if (dplyr::n_distinct(mono.markers$MARKERS) > 0) {
        input <- dplyr::filter(input, !MARKERS %in% mono.markers$MARKERS)
      }
    }
  } else {# for vcf haplotypes and multiallelic data
    n.catalog.locus <- dplyr::n_distinct(input$MARKERS)
    if (tibble::has_name(input, "REF")) {
      old.ref <- input %>%
        dplyr::distinct(MARKERS, REF, ALT)
    } else {
      old.ref <- NULL
    }


    if (tibble::has_name(input, "GT_HAPLO")) {
      input <- dplyr::select(input, MARKERS, INDIVIDUALS, GT_HAPLO, POP_ID)

      if (n.catalog.locus > 200000) {
        input <- input %>%
          dplyr::mutate(
            SPLIT_VEC = dplyr::ntile(x = 1:nrow(.), n = parallel.core * 3)) %>%
          split(x = ., f = .$SPLIT_VEC) %>%
          .radiator_parallel(
            X = .,
            FUN = gt_haplo2gt_vcf_nuc,
            mc.cores = parallel.core
          ) %>%
          dplyr::bind_rows(.) %>%
          dplyr::select(-SPLIT_VEC)
      } else {
        input <- input %>%
          dplyr::mutate(
            GT_VCF_NUC = dplyr::if_else(
              stringi::stri_detect_fixed(
                str = GT_HAPLO, pattern = "/"),
              GT_HAPLO,
              stringi::stri_join(GT_HAPLO, GT_HAPLO, sep = "/")),
            GT_VCF_NUC = stringi::stri_replace_na(str = GT_VCF_NUC, replacement = "./.")
          ) %>%
          dplyr::select(-GT_HAPLO)
      }
    }

    conversion.df <- dplyr::select(input, MARKERS, GT_VCF_NUC) %>%
      dplyr::left_join(
        dplyr::distinct(input, MARKERS) %>%
          dplyr::mutate(
            SPLIT_VEC = dplyr::ntile(x = 1:nrow(.), n = parallel.core * 3))
        , by = "MARKERS") %>%
      split(x = ., f = .$SPLIT_VEC) %>%
      .radiator_parallel(
        X = .,
        FUN = nuc2integers,
        mc.cores = parallel.core
      ) %>%
      dplyr::bind_rows(.)

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
    ref.alt.mono <- conversion.df %>%
      dplyr::distinct(MARKERS, REF, ALT, .keep_all = TRUE) %>%
      dplyr::select(-c(ALLELES, INTEGERS))
    conversion.df <- dplyr::select(conversion.df, MARKERS, ALLELES, INTEGERS)

    if (tibble::has_name(input, "REF")) {
      input <- dplyr::select(input, -c(REF, ALT))
    }
    input <- dplyr::left_join(input, ref.alt.mono, by = "MARKERS")
    ref.alt.mono <- NULL

    # Inversion: change in REF/ALT
    if (!is.null(old.ref)) {
      new.ref <- input %>%
        dplyr::distinct(MARKERS, REF, ALT) %>%
        dplyr::select(MARKERS, REF_NEW = REF)

      change.ref <- dplyr::full_join(new.ref, old.ref, by = "MARKERS") %>%
        dplyr::mutate(
          CHANGE = dplyr::if_else(REF == REF_NEW, "identical", "different")
        ) %>%
        dplyr::filter(CHANGE == "different") %>%
        dplyr::select(MARKERS) %>%
        purrr::flatten_chr(.)

      message("    Number of markers with REF/ALT change(s) = ", length(change.ref))
      change.ref <- new.ref <- old.ref <- NULL
    }
    if (verbose) message("Integrating new genotype codings...")
    new.gt <- dplyr::distinct(input, MARKERS, GT_VCF_NUC) %>%
      dplyr::mutate(SPLIT_VEC = dplyr::ntile(x = 1:nrow(.), n = parallel.core * 3)) %>%
      split(x = ., f = .$SPLIT_VEC) %>%
      .radiator_parallel(
        # parallel::mclapply(
        X = .,
        FUN = nuc2gt,
        mc.cores = parallel.core,
        conversion.data = conversion.df,
        biallelic = biallelic
      ) %>%
      dplyr::bind_rows(.)

    input <- dplyr::left_join(input, new.gt, by = c("MARKERS", "GT_VCF_NUC"))
    new.gt <- conversion.df <- NULL


    if (!is.null(markers.meta)) {
      no.need.markers.meta <- unique(colnames(markers.meta) %in% colnames(input))
      if (length(no.need.markers.meta) == 1 && !no.need.markers.meta) {
        input <- dplyr::left_join(input, markers.meta, by = "MARKERS") %>%
          dplyr::select(MARKERS, CHROM, LOCUS, POS, POP_ID, INDIVIDUALS, dplyr::everything())
      }
      markers.meta <- NULL
    }
  }

  res <- list(input = input, biallelic = biallelic)
  return(res)
}#End change_alleles
