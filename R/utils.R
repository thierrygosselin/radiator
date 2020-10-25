# replace_by_na ----------------------------------------------------------------
#' @title replace_by_na
#' @description Fast removal of NA
#' @rdname replace_by_na
#' @keywords internal
#' @export
replace_by_na <- function(data, what = ".") {
  replace(data, which(data == what), NA)
}#End replace_by_na

# distance2tibble---------------------------------------------------------------
#' @title distance2tibble
#' @description melt the distance matrix into a tibble
#' @rdname distance2tibble
#' @export
#' @keywords internal
distance2tibble <- function(
  x,
  remove.diag = TRUE,
  na.diag = FALSE,
  remove.lower = TRUE,
  relative = TRUE,
  pop.levels = NULL,
  distance.class.double = TRUE
) {
  # x <- dist.computation
  x <- as.matrix(x)
  if (remove.diag || na.diag) diag(x) <- NA
  if (remove.lower) x[lower.tri(x)] <- NA
  x <- dplyr::bind_cols(tibble::tibble(ID1 = rownames(x)),
                        tibble::as_tibble(x)) %>%
    data.table::as.data.table(.) %>%
    data.table::melt.data.table(
      data = ., id.vars = "ID1", variable.name = "ID2", value.name = "DISTANCE",
      variable.factor = FALSE) %>%
    tibble::as_tibble(.)

  if (na.diag || remove.diag) x  %<>% dplyr::filter(!is.na(DISTANCE))

  if (distance.class.double) {
    x %<>% dplyr::mutate(DISTANCE = as.double(as.character(DISTANCE)))
  }

  x %<>% dplyr::arrange(DISTANCE)

  if (relative) {
    x  %<>% dplyr::mutate(DISTANCE_RELATIVE = DISTANCE / max(DISTANCE))
  }
  x %<>% dplyr::mutate(dplyr::across(.cols = c("ID1", "ID2"), .fns = as.character))

  if (!is.null(pop.levels)) {
    x  %<>% dplyr::mutate(
      ID1 = factor(x = ID1, levels = pop.levels, ordered = TRUE),
      ID2 = factor(x = ID2, levels = pop.levels, ordered = TRUE)
    )
  } else {
    x  %<>% dplyr::mutate(
      ID1 = factor(x = ID1),
      ID2 = factor(x = ID2)
    )
  }

  return(x)
}#End distance2tibble


# split_vec ----------------------------------------------------------------
#' @title split_vec
#' @description Split input into chunk for parallel processing
#' @rdname split_vec
#' @keywords internal
#' @export
split_vec <- function(x, chunks) {
  if (any(class(x) %in% c("tbl_df","tbl","data.frame"))) x <- nrow(x)
  if (length(x) > 1L) x <- length(x)
  stopifnot(is.integer(x))
  split.vec <- as.integer(floor((chunks * (seq_len(x) - 1) / x) + 1))
  # split.vec <- as.integer(floor((parallel.core * cpu.rounds * (seq_len(x) - 1) / x) + 1))
  stopifnot(length(split.vec) == x)
  return(split.vec)
}#End split_vec




#' @title split_tibble_rows
#' @description Split rows of tibble for parallel processing
#' @rdname split_tibble_rows
#' @keywords internal
#' @export
split_tibble_rows <- function(
  x,
  lines.cpu = 1000, #lines per CPU rounds
  parallel.core = parallel::detectCores() - 1,
  group.split = TRUE # does dplyr: group_by and group_split
) {
  n.row <- nrow(x)
  n.cores <- parallel::detectCores()
  if (parallel.core > n.cores) parallel.core <- n.cores
  if (n.row < parallel.core) return(x)
  if (lines.cpu > n.row) lines.cpu <- n.row
  lines.rounds <- parallel.core * lines.cpu
  x$SPLIT_VEC <- sort(rep_len(x = 1:floor(n.row / lines.rounds), length.out = n.row))
  if (group.split) {
    x %<>%
      dplyr::group_by(SPLIT_VEC) %>%
      dplyr::group_split(.tbl = ., .keep = FALSE)
  }
  return(x)
}#End split_tibble_rows
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
  separate_genotype <-  carrier::crate(function(x, sep, gt, alleles.naming, gather, haplotypes, exclude){
    `%>%` <- magrittr::`%>%`
    `%<>%` <- magrittr::`%<>%`

    res <- tidyr::separate(
      data = x,
      col = gt,
      into = alleles.naming,
      sep = sep,
      extra = "drop", remove = TRUE
    )

    if (gather) {
      res %<>%
        radiator::rad_long(
          x = .,
          cols = exclude,
          names_to = "ALLELE_GROUP",
          values_to = "ALLELES",
          variable_factor = FALSE
        )
      if (haplotypes) res %<>% dplyr::rename(HAPLOTYPES = ALLELES)
    }

    return(res)
  })#End separate_genotype

  if (parallel.core > 1) {
    res <- radiator_future(
      .x = x,
      .f = separate_genotype,
      flat.future = "dfr",
      split.vec = TRUE,
      split.with = NULL,
      split.chunks = 10L,
      parallel.core = parallel.core,
      sep = sep,
      gt = gt,
      alleles.naming = alleles.naming,
      gather = gather,
      haplotypes = haplotypes,
      exclude = exclude
    )
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

  split_gt <- carrier::crate(function(x) {
    `%>%` <- magrittr::`%>%`
    `%<>%` <- magrittr::`%<>%`

    x %<>%
      tibble::as_tibble(x = .) %>%
      as.character(x = .) %>%
      stringi::stri_split_fixed(str = ., pattern = "/") %>% # a bit faster than strsplit with large dataset
      stats::setNames(object = ., nm = paste0("V", seq_along(.))) %>%
      tibble::as_tibble(x = .) %>%
      t(x = .) %>%
      tibble::as_tibble(x = .)
  })#End split_gt

  x  <- radiator_future(
    .x = x,
    .f = split_gt,
    flat.future = "dfc",
    parallel.core = parallel.core
  )
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


# parallel_core_opt ------------------------------------------------------------
#' @title parallel_core_opt
#' @description Optimization of parallel core argument for radiator
#' @keywords internal
#' @export
parallel_core_opt <- function(parallel.core = NULL, max.core = NULL) {
  # strategy:
  # minimum of 1 core and a maximum of all the core available -2
  # even number of core
  # test
  # parallel.core <- 1
  # parallel.core <- 2
  # parallel.core <- 3
  # parallel.core <- 11
  # parallel.core <- 12
  # parallel.core <- 16
  # max.core <- 5
  # max.core <- 50
  # max.core <- NULL

  # Add-ons options
  # to control the max and min number to use...

  if (is.null(parallel.core)) {
    parallel.core <- parallel::detectCores() - 2
  } else {
    parallel.core <- floor(parallel.core / 2) * 2
    parallel.core <- max(1, min(parallel.core, parallel::detectCores() - 2))
  }

  if (is.null(max.core)) {
    parallel.core.opt <- parallel.core
  } else {
    parallel.core.opt <- min(parallel.core, floor(max.core / 2) * 2)
  }
  return(parallel.core.opt)
}#End parallel_core_opt


# radiator_future --------------------------------------------------------------
#' @name radiator_future
#' @title radiator parallel function
#' @description Updating radiator to use future
# @inheritParams future::plan
# @inheritParams future::availableCores
#' @inheritParams future.apply::future_apply
#' @rdname radiator_future
#' @keywords internal
radiator_future <- function(
  .x,
  .f,
  flat.future = c("int", "chr", "dfr", "dfc", "walk", "drop"),
  split.vec = FALSE,
  split.with = NULL,
  split.chunks = 4L,
  parallel.core = parallel::detectCores() - 1,
  ...
) {
  opt.change <- getOption("width")
  options(width = 70)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(if (parallel.core > 1L) future::plan(strategy = "sequential"), add = TRUE)

  # argument for flattening the results
  flat.future <- match.arg(
    arg = flat.future,
    choices = c("int", "chr", "dfr", "dfc", "walk", "drop"),
    several.ok = FALSE
  )

  # splitting into chunks-------------------------------------------------------
  if (split.vec && is.null(split.with)) {
    # d: data, data length, data size
    # sv: split vector
    d <- .x
    df <- FALSE
    if (any(class(d) %in% c("tbl_df","tbl","data.frame"))) {
      d <- nrow(d)
      df <- TRUE
    }
    if (length(d) > 1L) d <- length(d)
    stopifnot(is.integer(d))
    sv <- as.integer(floor((split.chunks * (seq_len(d) - 1) / d) + 1))
    # sv <- as.integer(floor((parallel.core * cpu.rounds * (seq_len(d) - 1) / d) + 1))
    stopifnot(length(sv) == d)

    # split
    if (df) {
      .x$SPLIT_VEC <- sv
      .x %<>% dplyr::group_split(.tbl = ., "SPLIT_VEC", .keep = FALSE)
    } else {
      .x %<>% split(x = ., f = sv)
    }
  }
  if (!is.null(split.with)) {
    # check
    if (length(split.with) != 1 || !is.character(split.with)) {
      rlang::abort(message = "Contact author: problem with parallel computation")
    }
    .data <- NULL
    stopifnot(rlang::has_name(.x, split.with))
    if (split.vec) {
      sv <- dplyr::distinct(.x, .data[[split.with]])
      d <- nrow(sv)
      sv$SPLIT_VEC <- as.integer(floor((split.chunks * (seq_len(d) - 1) / d) + 1))
      .x %<>%
        dplyr::left_join(sv, by = split.with) %>%
        dplyr::group_split(.tbl = ., "SPLIT_VEC", .keep = FALSE)
    } else {
      .x %<>% dplyr::group_split(.tbl = ., .data[[split.with]], .keep = TRUE)
    }
  }



  if (parallel.core == 1L) {
    future::plan(strategy = "sequential")
  } else {
    parallel.core <- parallel_core_opt(parallel.core = parallel.core)
    lx <- length(.x)
    if (lx < parallel.core) {
      future::plan(strategy = "multisession", workers = lx)
    } else {
      future::plan(strategy = "multisession", workers = parallel.core)
    }
  }

  # .x <- future.apply::future_apply(X = .x, FUN = .f, ...)
  # capture dots
  # d <- rlang::dots_list(..., .ignore_empty = "all", .preserve_empty = TRUE, .homonyms = "first")
  # if (bind.rows) .x %<>% dplyr::bind_rows(.)



  # Run the function in parallel and account for dots-dots-dots argument
  rad_map <- switch(flat.future,
                    int = {furrr::future_map_int},
                    chr = {furrr::future_map_chr},
                    dfr = {furrr::future_map_dfr},
                    dfc = {furrr::future_map_dfc},
                    walk = {furrr::future_walk},
                    drop = {furrr::future_map}
  )

  opts <- furrr::furrr_options(globals = FALSE, seed = TRUE)
  if (length(list(...)) == 0) {
    .x %<>% rad_map(.x = ., .f = .f, .options = opts)
  } else {
    .x %<>% rad_map(.x = ., .f = .f, ..., .options = opts)
  }
  return(.x)
}#End radiator_future


# PIVOT-GATHER-CAST ------------------------------------------------------------
# rationale for doing this is that i'm tired of using tidyverse or data.table semantics
# tidyr changed from gather/spread to pivot_ functions but their are still very slow compared
# to 1. the original gather/spread and data.table equivalent...

#' @title rad_long
#' @description Gather, melt and pivot_longer
#' @rdname rad_long
#' @keywords internal
#' @export

rad_long <- function(
  x,
  cols = NULL,
  measure_vars = NULL,
  names_to = NULL,
  values_to = NULL,
  variable_factor = TRUE,
  keep_rownames = FALSE,
  tidy = FALSE
){


  # tidyr
  if (tidy) {
    x %>%
      tidyr::pivot_longer(
        data = .,
        cols = -cols,
        names_to = names_to,
        values_to = values_to
      )
  } else {# data.table
    x %>%
      data.table::as.data.table(., keep.rownames = keep_rownames) %>%
      data.table::melt.data.table(
        data = .,
        id.vars = cols,
        measure.vars = measure_vars,
        variable.name = names_to,
        value.name = values_to,
        variable.factor = variable_factor
      ) %>%
      tibble::as_tibble(.)
  }
}#rad_long
#' @rdname rad_wide
#' @title rad_wide
#' @description Spread, dcast and pivot_wider
#' @keywords internal
#' @export
rad_wide <- function(
  x ,
  formula = NULL,
  names_from = NULL,
  values_from = NULL,
  values_fill = NULL,
  sep = "_",
  tidy = FALSE

){
  # tidyr
  if (tidy) {
    x %<>%
      tidyr::pivot_wider(
        data = .,
        names_from = names_from,
        values_from = values_from,
        values_fill = values_fill
      )
  } else {# data.table
    x  %>%
      data.table::as.data.table(.) %>%
      data.table::dcast.data.table(
        data = .,
        formula =  formula,
        value.var = values_from,
        sep = sep,
        fill = values_fill
      ) %>%
      tibble::as_tibble(.)
  }
}#rad_wide


# strip_rad --------------------------------------------------------------------
#' @rdname strip_rad
#' @title strip_rad
#' @description Strip a tidy data set of it's strata and markers meta. Used internally.
#' @param x The data
#' @param m (character, string) The variables part of the markers metadata.
#' @param env.arg You want to redirect \code{rlang::current_env()} to this argument.
#' @param keep.strata (logical) Keep the strata in the dataset or remove the info
#' and keep only the sample ids.
#' @param verbose (logical) The function will chat more when allowed.
# @noRd
# @keywords internal
#' @export
strip_rad <- function(
  x,
  m = c("VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS", "COL", "REF", "ALT"),
  env.arg = NULL,
  keep.strata = TRUE,
  verbose = TRUE
) {
  objs <- utils::object.size(x)

  # STRATA ----------
  strata.n <- intersect(colnames(x), c("STRATA", "POP_ID"))

  if (rlang::has_name(x, "POP_ID")) {
    strata <- radiator::generate_strata(data = x, pop.id = TRUE) %>%
      dplyr::mutate(
        ID_SEQ = seq_len(length.out = dplyr::n()),
        STRATA_SEQ = as.integer(factor(x = POP_ID, levels = unique(POP_ID)))
      )
  } else {#STRATA
    strata <- radiator::generate_strata(data = x, pop.id = FALSE) %>%
      dplyr::mutate(
        ID_SEQ = seq_len(length.out = dplyr::n()),
        STRATA_SEQ = as.integer(factor(x = STRATA, levels = unique(STRATA)))
      )
  }


  cm <- intersect(colnames(strata), colnames(x))
  x %<>%
    dplyr::left_join(strata, by = cm) %>%
    dplyr::select(-tidyselect::any_of(cm))

  if (!keep.strata) x %<>% dplyr::select(-STRATA_SEQ)

  # Note to myself: maybe write using vroom and read back when necessary ?
  assign(
    x = "strata.bk",
    value = strata,
    pos = env.arg,
    envir = env.arg
  )
  cm <- keep.strata <- pop.id <- strata.n <- strata <- NULL

  # MARKERS ---------
  x %<>%
    dplyr::mutate(
      M_SEQ = as.integer(factor(x = MARKERS, levels = unique(MARKERS)))
    )
  m <- c("M_SEQ", m)
  markers.meta <- x %>%
    dplyr::select(tidyselect::any_of(m)) %>%
    dplyr::distinct(M_SEQ, .keep_all = TRUE)

  cm <- intersect(colnames(markers.meta), colnames(x)) %>%
    purrr::discard(.x = ., .p = . %in% "M_SEQ")

  # remove markers meta
  x %<>%
    dplyr::select(-tidyselect::any_of(cm))

  # Note to myself: maybe write using vroom and read back when necessary ?
  assign(
    x = "markers.meta.bk",
    value = markers.meta,
    pos = env.arg,
    envir = env.arg
  )
  markers.meta <- cm <- m <- NULL

  # GENOTYPE META ---------
  # g <- c("READ_DEPTH", "ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH",
  # "GL", "CATG", "PL", "HQ", "GQ", "GOF","NR", "NV")
  want <- c("GT", "GT_VCF_NUC", "GT_VCF", "GT_BIN", "REF", "ALT", "MARKERS",
            "CHROM", "LOCUS", "POS", "INDIVIDUALS", "POP_ID"
  )
  g <- purrr::keep(.x = colnames(x), .p = !colnames(x) %in% want) %>%
    purrr::discard(.x = ., .p = . %in% c("STRATA_SEQ"))
  # purrr::discard(.x = ., .p = . %in% c("ID_SEQ", "STRATA_SEQ", "M_SEQ"))
  # g
  genotypes.meta <- NULL # default

  if (length(g) != 0L) {
    genotypes.meta <- x %>% dplyr::select(tidyselect::any_of(g))
    want <- c("ID_SEQ", "STRATA_SEQ", "M_SEQ", "GT", "GT_VCF_NUC", "GT_VCF", "GT_BIN", "REF", "ALT")
    x %<>% dplyr::select(tidyselect::any_of(want))
  }

  assign(
    x = "genotypes.meta.bk",
    value = genotypes.meta,
    pos = env.arg,
    envir = env.arg
  )
  g <- want <- genotypes.meta <- NULL

  if (verbose) message("Proportion of size reduction: ", round(1 - (utils::object.size(x) / objs), 2))
  return(x)
  # return(list(data = x, markers.meta = markers.meta, genotypes.meta = genotypes.meta))
}# End strip_rad


# join_rad ---------------------------------------------------------------------
#' @rdname join_rad
#' @title join_rad
#' @description Join  back the parts stripped in strip_rad. Used internally.
#' @param x The data
#' @param s The strata metadata.
#' @param m The markers metadata.
#' @param g The genotypes metadata.
#' @param env.arg You want to redirect \code{rlang::current_env()} to this argument.
# @noRd
# @keywords internal
#' @export
join_rad <- function(x, s, m, g, env.arg = NULL) {
  if (!is.null(g)) {
    x %<>% dplyr::left_join(g, by = c("ID_SEQ", "M_SEQ"))
    env.arg$genotypes.meta.bk <- NULL
  }

  want <- c("VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS",
            "COL", "REF", "ALT", "POP_ID", "STRATA","INDIVIDUALS")
  x %<>%
    dplyr::left_join(s, by = "ID_SEQ") %>%
    dplyr::select(-tidyselect::any_of(c("ID_SEQ", "STRATA_SEQ"))) %>%
    dplyr::left_join(m, by = "M_SEQ") %>%
    dplyr::select(-M_SEQ) %>%
    dplyr::select(tidyselect::any_of(want), tidyselect::everything())

  env.arg$strata.bk <- env.arg$markers.meta.bk <- NULL
  return(x)
}#End join_rad
