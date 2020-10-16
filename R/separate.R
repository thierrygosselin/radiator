# separate_markers -------------------------------------------------------------

# Separate a column (markers) into CHROM LOCUS and POS
# generate markers meta

#' @name separate_markers
#' @title Separate markers column into chrom, locus and pos

#' @description Radiator uses unique marker names by combining
#' \code{CHROM}, \code{LOCUS}, \code{POS} columns, with double underscore
#' separators, into \code{MARKERS = CHROM__LOCUS__POS}.
#'
#'
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users who need to get back to the original metadata the
#' function provides an easy way to do it.

#' @param data An object with a column named \code{MARKERS}.
#' If \code{CHROM}, \code{LOCUS}, \code{POS} are already present, the function
#' returns the dataset untouched.
#' The data can be whitelists and blacklists of markers or tidy datasets or
#' radiator GDS object.

#' @param sep (optional, character) Separator used to identify the different
#' field in the \code{MARKERS} column.
#'
#' When the \code{MARKERS} column doesn't have separator and the function is used
#' to generate markers metadata column:
#' \code{"CHROM", "LOCUS", "POS", "REF", "ALT"}, use \code{sep = NULL}.
#' Default: \code{sep = "__"}.
#'
#' @param markers.meta.lists.only (logical, optional)
#' Allows to keep only the markers metadata:
#' \code{"VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS"}, useful for whitelist
#' or blacklist.
#' Default: \code{markers.meta.lists.only = FALSE}

#' @param markers.meta.all.only (logica, optionall)
#' Allows to keep all available markers metadata:
#' \code{"VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS", "COL", "REF", "ALT"},
#' useful inside radiator.
#' Default: \code{markers.meta.all.only = FALSE}

#' @param generate.markers.metadata (logical, optional)
#' Generate missing markers metadata when missing.
#' \code{"CHROM", "LOCUS", "POS"}.
#' Default: \code{generate.markers.metadata = TRUE}
#'
#' @param generate.ref.alt (logical, optional) Generate missing REF/ALT alleles
#' with: REF = A and ALT = C (for biallelic datasets, only).
#' It is turned off automatically
#' when argument \code{markers.meta.lists.only = TRUE} and
#' on automatically when argument \code{markers.meta.all.only = TRUE}
#' Default: \code{generate.ref.alt = FALSE}

#' @param biallelic (logical) Speed up the function execution by entering
#' if the dataset is biallelic or not. Used internally for verification, before
#' generating REF/ALT info.
#' By default, the function calls \code{\link{detect_biallelic_markers}}.
#' The argument is required if \code{data} is a tidy dataset and not just
#' a whitelist/blacklist.
#' Default: \code{biallelic = NULL}
#' @inheritParams tidy_genomic_data

#' @return The same data in the global environment, with 3 new columns:
#' \code{CHROM}, \code{LOCUS}, \code{POS}. Additionnal columns may be genrated,
#' see arguments documentation.
#' @rdname separate_markers

#' @examples
#' \dontrun{
#' whitelist <- radiator::separate_markers(data = whitelist.markers)
#' tidy.data <- radiator::separate_markers(data = bluefintuna.data)
#' }
#' @export

#' @seealso \code{\link{detect_biallelic_markers}} and \code{\link{generate_markers_metadata}}

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

separate_markers <- function(
  data,
  sep = "__",
  markers.meta.lists.only = FALSE,
  markers.meta.all.only = FALSE,
  generate.markers.metadata = TRUE,
  generate.ref.alt = FALSE,
  biallelic = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE
) {
  # data.bk <- data
  # sep <- "__"

  # check if markers column is present
  if (!tibble::has_name(data, "MARKERS")) {
    rlang::abort("The data require a column named MARKERS")
  }

  n.markers <- length(unique(data$MARKERS))
  unique.markers <- nrow(data) == n.markers

  if (unique.markers && generate.ref.alt && is.null(biallelic)) {
    rlang::abort("biallelic TRUE/FALSE required")
  }

  if (markers.meta.lists.only) {
    want <- c("VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS")
    suppressWarnings(data %<>% dplyr::select(dplyr::one_of(want)) %>%
                       dplyr::distinct(MARKERS, .keep_all = TRUE))
    generate.ref.alt <- FALSE
    generate.markers.metadata <- FALSE
  }

  if (markers.meta.all.only) {
    notwanted <- c("GT_BIN", "GT", "GT_VCF", "GT_VCF_NUC", "DP", "AD", "GL",
                   "PL", "GQ", "HQ", "GOF", "NR", "NV", "CATG")
    # want <- c("FILTERS", "VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS",
    #           "COL", "REF", "ALT", "CALL_RATE", "AVG_COUNT_REF",
    #           "AVG_COUNT_SNP", "REP_AVG", "ONE_RATIO_REF", "ONE_RATIO_SNP",
    #           "SEQUENCE")
    suppressWarnings(data %<>% dplyr::select(-dplyr::one_of(notwanted)) %>%
                       dplyr::distinct(MARKERS, .keep_all = TRUE))
    generate.markers.metadata <- generate.ref.alt <- TRUE
  }

  if (!is.null(sep)) {
    rad.sep <- unique(
      stringi::stri_detect_fixed(
        str = sample(x = unique(data$MARKERS), size = min(200, length(data$MARKERS))),
        pattern = sep))

    if (length(rad.sep) != 1) rlang::abort("More than 1 separator was detected")
    if (!rad.sep) {
      message("The separator specified is not valid")
    } else {
      if (FALSE %in% unique(c("CHROM", "LOCUS", "POS") %in% colnames(data))) {
        rad.sep <- TRUE
      } else {
        rad.sep <- FALSE
      }

      if (rad.sep) {
        # Note to myself: this section could be parallelized when whole dataset are required
        want <- c("CHROM", "LOCUS", "POS")

        if (unique.markers) {
          data <- suppressWarnings(
            data %>%
              dplyr::select(-dplyr::one_of(want)) %>%
              tidyr::separate(data = ., col = "MARKERS",
                              into = want, sep = sep, remove = FALSE))
        } else {# for datasets
          temp <- suppressWarnings(tidyr::separate(
            data = dplyr::distinct(data, MARKERS),
            col = "MARKERS",
            into = want,
            sep = sep,
            remove = FALSE))

          suppressWarnings(data %<>% dplyr::select(-dplyr::one_of(want)))
          data %<>% dplyr::left_join(temp, by = intersect(colnames(data),
                                                          colnames(temp)))
          temp <- NULL
        }
      }
    }
  }# End of splitting markers column

  # Generate missing markers meta
  if (generate.markers.metadata) {
    data <- generate_markers_metadata(
      data = data,
      generate.markers.metadata = generate.markers.metadata,
      generate.ref.alt = generate.ref.alt,
      biallelic = biallelic,
      parallel.core = parallel.core, verbose = verbose)
  }
  return(data)
}#End separate_markers

# generate_markers_metadata ----------------------------------------------------
#' @name generate_markers_metadata
#' @title Generate markers metadata

#' @description Generate markers metadata: \code{CHROM, LOCUS, POS, REF, ALT}
#' when missing from tidy datasets.

#' @inheritParams separate_markers
#' @inheritParams tidy_genomic_data

#' @return Depending on argument's value, the same data is returned
#' in the global environment, with potential these additional columns:
#' \code{CHROM, LOCUS, POS, REF, ALT}.
#' @rdname generate_markers_metadata

#' @examples
#' \dontrun{
#' tidy.data <- radiator::generate_markers_metadata(data = bluefintuna.data)
#' }
#' @export
#' @seealso \code{\link{detect_biallelic_markers}} and \code{\link{separate_markers}}


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
generate_markers_metadata <- function(
  data,
  generate.markers.metadata = TRUE,
  generate.ref.alt = FALSE,
  biallelic = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE
) {
  if (!generate.markers.metadata && generate.ref.alt) {
    if (verbose) message("generate.markers.metadata: switched to TRUE automatically")
    generate.markers.metadata <- TRUE
  }

  if (generate.markers.metadata) {
    n.markers <- length(unique(data$MARKERS))

    unique.markers <- nrow(data) == n.markers


    if (unique.markers && generate.ref.alt && is.null(biallelic)) {
      rlang::abort("biallelic TRUE/FALSE required")
    }

    want <- c("FILTERS", "VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS", "COL", "REF",
              "ALT", "CALL_RATE", "AVG_COUNT_REF", "AVG_COUNT_SNP", "REP_AVG",
              "ONE_RATIO_REF", "ONE_RATIO_SNP", "SEQUENCE")
    notwanted <- c("GT_BIN", "GT", "GT_VCF", "GT_VCF_NUC", "DP", "AD", "GL",
                                "PL", "GQ", "HQ", "GOF", "NR", "NV", "CATG")
    if (!unique.markers) {
      suppressWarnings(
        markers.meta <- data %>%
          # dplyr::select(dplyr::one_of(want)) %>%
          dplyr::select(-dplyr::one_of(notwanted)) %>%
          dplyr::distinct(MARKERS, .keep_all = TRUE))
    } else {
      markers.meta <- suppressWarnings(dplyr::select(data, -dplyr::one_of(notwanted)))
      data <- NULL
    }

    if (!rlang::has_name(markers.meta, "VARIANT_ID")) {#nrow(markers.meta) == n.markers
      markers.meta %<>% dplyr::mutate(VARIANT_ID = as.integer(factor(MARKERS)))
    }

    if (!tibble::has_name(markers.meta, "CHROM")) {
      markers.meta %<>% dplyr::mutate(CHROM = rep("CHROM_1", n.markers))
      if (verbose) message("CHROM info missing: 'CHROM_1' integer was added to dataset")
    }

    # Generate LOCUS info if not present
    if (!tibble::has_name(markers.meta, "LOCUS")) {
      markers.meta %<>% dplyr::mutate(LOCUS = seq(1, n.markers, by = 1))
      if (verbose) message("LOCUS info missing: unique integers were added to dataset")
    }

    if (!tibble::has_name(markers.meta, "POS")) {
      markers.meta %<>% dplyr::mutate(CHROM = rep(1L, n.markers))
      if (verbose) message("POS info missing: dataset filled with MARKERS column")
    }

    # Generate REF/ALT allele if not in dataset
    if (generate.ref.alt) {
      if (tibble::has_name(markers.meta, "REF")) generate.ref.alt <- FALSE
      if (is.null(biallelic)) {
        biallelic <- radiator::detect_biallelic_markers(data = data, verbose = FALSE, parallel.core = parallel.core)
      }
      if (!biallelic) generate.ref.alt <- FALSE
    }

    if (generate.ref.alt) {
      markers.meta %<>% dplyr::mutate(REF = "A", ALT = "C")
      if (verbose) message("REF and ALT allele info missing: setting REF = A and ALT = C")
    }

    if (!unique.markers) {
      data %<>% dplyr::left_join(markers.meta, by = intersect(colnames(data),
                                                              colnames(markers.meta)))
    } else {
      data <- markers.meta
    }
    markers.meta <- NULL
  }
  return(data)
}#End generate_markers_metadata



# separate_gt ------------------------------------------------------------------
#' @title separate_gt
#' @description Separate genotype field
#' @rdname separate_gt
#' @keywords internal
#' @export
separate_gt <- function(
  x,
  sep = "/",
  gt = "GT_VCF_NUC",
  gather = TRUE,
  haplotypes = TRUE,
  exclude = c("LOCUS", "INDIVIDUALS", "POP_ID"),
  alleles.naming = c("ALLELE1", "ALLELE2"),
  cpu.rounds = 10,
  parallel.core = parallel::detectCores() - 1
) {
  # sep <-  "/"
  # x <- input2
  separate_genotype <- function(x, sep, gt, alleles.naming, gather, haplotypes, exclude){
    res <- tidyr::separate(
      data = x,
      col = gt,
      into = alleles.naming,
      sep = sep,
      extra = "drop", remove = TRUE
    )

    if (gather) {
      # res <- tidyr::gather(
      #   data = res,
      #   key = ALLELE_GROUP,
      #   value = HAPLOTYPES,
      #   -dplyr::one_of(exclude)
      # )

      res <- data.table::melt.data.table(
        data = data.table::as.data.table(res),
        id.vars = exclude,
        variable.name = "ALLELE_GROUP",
        value.name = "HAPLOTYPES",
        variable.factor = FALSE) %>%
        tibble::as_tibble(.)

      if (!haplotypes) {
        res  %<>% dplyr::rename(ALLELES = HAPLOTYPES)
      }
    }

    return(res)
  }

  if (parallel.core > 1) {
    n.row <- nrow(x)
    split.vec <- as.integer(floor((parallel.core * cpu.rounds * (1:n.row - 1) / n.row) + 1))
    res <- split(x = x, f = split.vec) %>%
      radiator_parallel_mc(
        X = .,
        FUN = separate_genotype,
        mc.cores = parallel.core,
        sep = sep,
        gt = gt,
        alleles.naming = alleles.naming,
        gather = gather,
        haplotypes = haplotypes,
        exclude = exclude
      ) %>%
      dplyr::bind_rows(.)
  } else {
    res <- separate_genotype(
      x = x,
      sep = sep,
      gt = gt,
      alleles.naming = alleles.naming,
      gather = gather,
      haplotypes = haplotypes,
      exclude = exclude)
  }
  return(res)
}#End separate_gt

# radiator_split_tibble ------------------------------------------------------------------
#' @title radiator_split_tibble
#' @description radiator alternative to data.table::tstrsplit
#' @rdname radiator_split_tibble
#' @keywords internal
#' @export
# radiator alternative to data.table::tstrsplit
# required function that is identical to data.table::tstrsplit
radiator_split_tibble <- function(x, parallel.core = parallel::detectCores() - 1) {

  split_gt <- function(x) {
    x %<>%
      as.character(x = .) %>%
      stringi::stri_split_fixed(str = ., pattern = "/") %>% # a bit faster than strsplit with large dataset
      # strsplit(x = ., split = "/", fixed = TRUE) %>%
      stats::setNames(object = ., nm = paste0("V", seq_along(.))) %>%
      tibble::as_tibble(x = .) %>%
      t(x = .) %>%
      tibble::as_tibble(x = .)
  }#End split_gt

  x %<>%
    tibble::as_tibble(x = .) %>%
    radiator_parallel_mc(
      X = .,
      FUN = split_gt,
      mc.cores = parallel.core
    ) %>%
    dplyr::bind_cols(.)
  return(x)
}#End radiator_split_tibble
# Note to myself:
# what was tried but was not adopted because to slow or else:

# # fast for small dataset, doesn't scale well

# split_bind_tibble <- function(x) {
#   stringi::stri_split_fixed(str = x, pattern = "/") %>%
#     # purrr::map(.x = ., .f = rbind) %>% # doesnt work because no names
#     do.call(what = rbind) %>%
#     tibble::as_tibble(.)
# }

# data.table alternative
# # fast for small dataset, doesn't scale well
# test <- purrr::map_dfc(.x = gt2, .f = data.table::tstrsplit, "/")

# parallel version is fast, small cost with < 2000 markers
# tictoc::tic()
# test <- radiator_parallel_mc(
#   X = gt2,
#   FUN = data.table::tstrsplit,
#   "/",
#   mc.cores = 12
# ) %>%
#   dplyr::bind_cols(.)
# tictoc::toc()

# melting and casting seems counter productive but it's actually pretty much the same as above
# tictoc::tic()
# test <- magrittr::set_colnames(x = gt, value = markers) %>%
#   tibble::as_tibble(x = ., rownames = "INDIVIDUALS") %>%
#   data.table::as.data.table(.) %>%
#   data.table::melt.data.table(
#     data = .,
#     id.vars = "INDIVIDUALS",
#     variable.name = "MARKERS",
#     value.name = "GT_VCF_NUC",
#     variable.factor = FALSE) %>%
#   tibble::as_tibble(.) %>%
#   # radiator::separate_gt(
#     separate_gt(
#       x = .,
#     sep = "/",
#     gather = TRUE,
#     alleles.naming = c("A1", "A2"),
#     haplotypes = FALSE,
#     exclude = c("MARKERS", "INDIVIDUALS"),
#     gt = "GT_VCF_NUC",
#     parallel.core = parallel.core) %>%
#   dplyr::mutate(
#     # ALLELE_GROUP = stringi::stri_replace_all_fixed(
#     #   str = ALLELE_GROUP,
#     #   pattern = "ALLELE",
#     #   replacement = "A",
#     #   vectorize_all = FALSE),
#     ALLELES = replace(x = ALLELES, which(ALLELES == "."), NA)
#   ) %>%
#   tidyr::unite(col = MARKERS_ALLELES, MARKERS , ALLELE_GROUP, sep = ".") %>%
#   dplyr::arrange(INDIVIDUALS, MARKERS_ALLELES) %>%
#   data.table::as.data.table(.) %>%
#   data.table::dcast.data.table(
#     data = .,
#     formula = INDIVIDUALS ~ MARKERS_ALLELES,
#     value.var = "ALLELES"
#   ) %>%
#   tibble::as_tibble(.) %>%
#   dplyr::ungroup(.)
# tictoc::toc()
