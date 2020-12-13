# replace_by_na ----------------------------------------------------------------
#' @title replace_by_na
#' @description Fast replacement of values by NA
#' @rdname replace_by_na
#' @keywords internal
#' @export
replace_by_na <- function(data, what = ".") {
  data  %<>% dplyr::na_if(x = ., y = what)
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
    data %<>%
      dplyr::select(tidyselect::any_of(want)) %>%
      dplyr::distinct(MARKERS, .keep_all = TRUE)
    generate.ref.alt <- FALSE
    generate.markers.metadata <- FALSE
  }

  if (markers.meta.all.only) {
    notwanted <- c("GT_BIN", "GT", "GT_VCF", "GT_VCF_NUC", "DP", "AD", "GL",
                   "PL", "GQ", "HQ", "GOF", "NR", "NV", "CATG")
    data %<>%
      dplyr::select(-tidyselect::any_of(notwanted)) %>%
      dplyr::distinct(MARKERS, .keep_all = TRUE)
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
          data %<>%
              dplyr::select(-tidyselect::any_of(want)) %>%
              tidyr::separate(data = ., col = "MARKERS", into = want, sep = sep, remove = FALSE)
        } else {# for datasets
          temp <- tidyr::separate(
            data = dplyr::distinct(data, MARKERS),
            col = "MARKERS",
            into = want,
            sep = sep,
            remove = FALSE)

          data %<>%
            dplyr::select(-tidyselect::any_of(want)) %>%
            dplyr::left_join(temp, by = intersect(colnames(data), colnames(temp)))
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
      markers.meta <- data %>%
        dplyr::select(-tidyselect::any_of(notwanted)) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE)
    } else {
      markers.meta <- dplyr::select(data, -tidyselect::any_of(notwanted))
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
  gt = "GT_VCF_NUC",
  gather = TRUE,
  haplotypes = FALSE,
  exclude = c("LOCUS", "INDIVIDUALS", "POP_ID"),
  alleles.naming = c("A1", "A2"),
  remove = TRUE,
  filter.missing = FALSE,
  split.chunks = 3,
  parallel.core = parallel::detectCores() - 1
) {

  ## TEST
  # gather = TRUE
  # haplotypes = FALSE
  # exclude = c("LOCUS", "INDIVIDUALS", "POP_ID")
  # alleles.naming = c("A1", "A2")
  # remove = TRUE
  # filter.missing = FALSE
  # split.chunks = 3


  separate_genotype <-  carrier::crate(function(x, gt, alleles.naming, remove, filter.missing, gather, haplotypes, exclude){
    `%>%` <- magrittr::`%>%`
    `%<>%` <- magrittr::`%<>%`

    # discard the other gt format
    gt.format <- c("GT", "GT_BIN", "GT_VCF", "GT_VCF_NUC")
    not.wanted <- setdiff(gt.format, gt)
    x %<>% dplyr::select(-tidyselect::any_of(not.wanted))

    if (gt == "GT_VCF_NUC") {
      if (filter.missing) x  %<>% dplyr::filter(GT_VCF_NUC != "./.")
      x %<>%
        dplyr::bind_cols(
          stringi::stri_split_fixed(str = x$GT_VCF_NUC, pattern = "/", simplify = TRUE) %>%
            magrittr::set_colnames(x = ., value = alleles.naming) %>%
            tibble::as_tibble()
        )
    }
    if (gt == "GT_VCF") {
      if (filter.missing) x  %<>% dplyr::filter(GT_VCF != "./.")
      x %<>%
        dplyr::mutate(
          A1 = stringi::stri_sub(str = GT_VCF, from = 1, to = 1),
          A2 = stringi::stri_sub(str = GT_VCF, from = 3, to = 3)
        )
    }
    if (gt == "GT_BIN") {
      if (filter.missing) x  %<>% dplyr::filter(!is.na(GT_BIN))
      x %<>%
        dplyr::mutate(
          A1 = dplyr::if_else(GT_BIN == 0L, 1L, GT_BIN),
          A2 = dplyr::if_else(GT_BIN != 2L, GT_BIN + 1L, GT_BIN)
        )
    }
    if (gt == "GT") {
      if (filter.missing) x  %<>% dplyr::filter(GT != "000000")
      x %<>%
        dplyr::mutate(
          A1 = stringi::stri_sub(str = GT, from = 1, to = 3),
          A2 = stringi::stri_sub(str = GT, from = 4, to = 6)
        )
    }

    if (remove) x %<>% dplyr::select(-tidyselect::any_of(gt.format))

    if (gather) {
      x %<>%
        radiator::rad_long(
          x = .,
          cols = exclude,
          names_to = "ALLELES_GROUP",
          values_to = "ALLELES",
          variable_factor = FALSE
        )
      if (haplotypes) x %<>% dplyr::rename(HAPLOTYPES = ALLELES)
    }

    return(x)
  })#End separate_genotype


  if (split.chunks > 1) {
    x %<>%
      radiator_future(
        .x = .,
        .f = separate_genotype,
        flat.future = "dfr",
        split.vec = TRUE,
        split.with = NULL,
        split.chunks = split.chunks,
        parallel.core = split.chunks,
        forking = TRUE,
        gt = gt,
        alleles.naming = alleles.naming,
        remove = remove,
        filter.missing = filter.missing,
        gather = gather,
        haplotypes = haplotypes,
        exclude = exclude
      )
  } else {
    x %<>%
      separate_genotype(
        x = .,
        gt = gt,
        alleles.naming = alleles.naming,
        remove = remove,
        filter.missing = filter.missing,
        gather = gather,
        haplotypes = haplotypes,
        exclude = exclude)
  }
  return(x)
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
  forking = FALSE,
  ...
) {
  os <- Sys.info()[['sysname']]
  if (os == "Windows") forking <- FALSE

  opt.change <- getOption("width")
  options(width = 70)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(if (parallel.core > 1L && !forking) future::plan(strategy = "sequential"), add = TRUE)

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
      .x %<>% dplyr::ungroup(.) %>% dplyr::group_split(.tbl = ., "SPLIT_VEC", .keep = FALSE)
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
        dplyr::ungroup(.) %>%
        dplyr::group_split(.tbl = ., "SPLIT_VEC", .keep = FALSE)
    } else {
      .x %<>% dplyr::ungroup(.) %>% dplyr::group_split(.tbl = ., .data[[split.with]], .keep = TRUE)
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
      if (!forking) future::plan(strategy = "multisession", workers = parallel.core)
    }
  }

  # Run the function in parallel and account for dots-dots-dots argument

  if (forking) {
    if (length(list(...)) == 0) {
      rad_map <- switch(
        flat.future,
        int = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, mc.cores = parallel.core) %>%
            purrr::flatten_int(.)
        },
        chr = {.x %<>%
            parallel::mclapply(X = ., FUN = .f, mc.cores = parallel.core) %>%
            purrr::flatten_chr(.)
        },
        dfr = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, mc.cores = parallel.core) %>%
            purrr::flatten_dfr(.)
        },
        dfc = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, mc.cores = parallel.core) %>%
            purrr::flatten_dfc(.)
        },
        walk = {furrr::future_walk},
        drop = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, mc.cores = parallel.core)
        }
      )
    } else {
      rad_map <- switch(
        flat.future,
        int = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, ..., mc.cores = parallel.core) %>%
            purrr::flatten_int(.)
        },
        chr = {.x %<>%
            parallel::mclapply(X = ., FUN = .f, ..., mc.cores = parallel.core) %>%
            purrr::flatten_chr(.)
        },
        dfr = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, ..., mc.cores = parallel.core) %>%
            purrr::flatten_dfr(.)
        },
        dfc = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, ..., mc.cores = parallel.core) %>%
            purrr::flatten_dfc(.)
        },
        walk = {furrr::future_walk},
        drop = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, ..., mc.cores = parallel.core)
        }
      )
    }
  } else {
    rad_map <- switch(flat.future,
                      int = {furrr::future_map_int},
                      chr = {furrr::future_map_chr},
                      dfr = {furrr::future_map_dfr},
                      dfc = {furrr::future_map_dfc},
                      walk = {furrr::future_walk},
                      drop = {furrr::future_map}
    )
    p <- NULL
    p <- progressr::progressor(along = .x)
    opts <- furrr::furrr_options(globals = FALSE, seed = TRUE)
    if (length(list(...)) == 0) {
      .x %<>% rad_map(.x = ., .f = .f, .options = opts)
    } else {
      .x %<>% rad_map(.x = ., .f = .f, ..., .options = opts)
    }
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
    x %<>%
      tidyr::pivot_longer(
        data = .,
        cols = -cols,
        names_to = names_to,
        values_to = values_to
      )
  } else {# data.table
    x %<>%
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

  if (names_to == "M_SEQ") x$M_SEQ %<>% as.integer
  if (names_to == "ID_SEQ") x$ID_SEQ %<>% as.integer
  return(x)
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
  x %<>% dplyr::select(-tidyselect::any_of(cm))

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
#' @description Join back the parts stripped in strip_rad. Used internally.
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
    g.by <- intersect(colnames(x), colnames(g))
    x %<>% dplyr::left_join(g, by = g.by)
    env.arg$genotypes.meta.bk <- NULL
  }

  want <- c("VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS",
            "COL", "REF", "ALT", "INDIVIDUALS", "POP_ID", "STRATA",
            "GT_VCF", "GT_VCF_NUC", "GT", "GT_BIN")
  s.by <- intersect(colnames(x), colnames(s))
  m.by <- intersect(colnames(x), colnames(m))
  x %<>%
    dplyr::left_join(s, by = s.by) %>%
    dplyr::select(-tidyselect::any_of(c("ID_SEQ", "STRATA_SEQ"))) %>%
    dplyr::left_join(m, by = m.by) %>%
    dplyr::select(-M_SEQ) %>%
    dplyr::select(tidyselect::any_of(want), tidyselect::everything())

  env.arg$strata.bk <- env.arg$markers.meta.bk <- NULL
  return(x)
}#End join_rad


# detect_gt ---------------------------------------------------------------------
#' @rdname detect_gt
#' @title detect_gt
#' @description Detect the genotype format used in the data set. Sample one if multiple are present.
#' @param x The data
#' @param gt.format (character)
#' Default: \code{gt.format = c("GT", "GT_BIN", "GT_VCF", "GT_VCF_NUC")}.
#' @param keep.one (logical) Will return only one format if \code{keep.one = TRUE}.
#' Default: \code{keep.one = TRUE}.
#' @param favorite If more than one format is present and \code{keep.one = TRUE},
#' the favorite will be returned, if present. Sample one if not present.
#' Default: \code{favorite = "GT_BIN"}.
#' @keywords internal
#' @export
detect_gt <- function(x, gt.format = c("GT", "GT_BIN", "GT_VCF", "GT_VCF_NUC"), keep.one = TRUE, favorite = "GT_BIN") {

  detect.gt <- purrr::keep(.x = colnames(x), .p = colnames(x) %in% gt.format)

  if (keep.one) {
    if (length(detect.gt) > 1) {
      if (favorite %in% detect.gt) {
        detect.gt <- favorite
      } else {
        detect.gt %<>% sample(x = ., size = 1)
      }
    }
  }
  return(detect.gt)
}#End detect_gt



# gt_recoding ---------------------------------------------------------------------
#' @rdname gt_recoding
#' @title gt_recoding
#' @description Detect the genotype format used in the data set. Sample one if than one is present.
#' @keywords internal
#' @export
gt_recoding <- function(x, gt = TRUE, gt.bin = TRUE, gt.vcf = TRUE, gt.vcf.nuc = TRUE, arrange = TRUE) {

  # conditions and checks
  if (arrange) x %<>% dplyr::mutate(IDTEMP = seq_len(dplyr::n()))
  if (gt.vcf.nuc && !rlang::has_name(x, "REF")) gt.vcf.nuc <- FALSE


  # what genotype format we have
  detect.gt <- detect_gt(x) #utils
  remove.extra <- FALSE
  if (gt) {
    if (!rlang::has_name(x, "A1")) {
      if (rlang::has_name(x, "REF")) {
        x  %<>%
          dplyr::mutate(
            A1 = dplyr::recode(REF, "A" = "001", "C" = "002", "G" = "003", "T" = "004"),
            A2 = dplyr::recode(ALT, "A" = "001", "C" = "002", "G" = "003", "T" = "004")
          )
        remove.extra <- TRUE
      } else {# if no REF
        gt <- FALSE
      }
    }
  }



  gt_map <- function(
    x,
    gt.format = c("GT", "GT_BIN", "GT_VCF", "GT_VCF_NUC"),
    gt = TRUE,
    gt.bin = TRUE,
    gt.vcf = TRUE,
    gt.vcf.nuc = TRUE
  ) {

    gt.format <- match.arg(
      arg = gt.format,
      choices = c("GT", "GT_BIN", "GT_VCF", "GT_VCF_NUC"),
      several.ok = FALSE
    )

    #start with missing genotypes
    if (gt.format == "GT_BIN") {
      gt.bin <- unique(x$GT_BIN)
      if (is.na(gt.bin)) {
        x %<>%
          {if (gt.vcf) dplyr::mutate(.data = ., GT_VCF = "./.") else .} %>%
          {if (gt.vcf.nuc) dplyr::mutate(.data = ., GT_VCF_NUC = "./.") else .} %>%
          {if (gt) dplyr::mutate(.data = ., GT = "000000") else .}
      } else {
        if (gt.bin == 0L) {
          x %<>%
            {if (gt.vcf) dplyr::mutate(.data = ., GT_VCF = "0/0") else .} %>%
            {if (gt.vcf.nuc) dplyr::mutate(.data = ., GT_VCF_NUC = stringi::stri_join(REF, REF, sep = "/")) else .} %>%
            {if (gt) dplyr::mutate(.data = ., GT = stringi::stri_join(A1, A1)) else .}
        }
        if (gt.bin == 1L) {
          x %<>%
            {if (gt.vcf) dplyr::mutate(.data = ., GT_VCF = "0/1") else .} %>%
            {if (gt.vcf.nuc) dplyr::mutate(.data = ., GT_VCF_NUC = stringi::stri_join(REF, ALT, sep = "/")) else .} %>%
            {if (gt) dplyr::mutate(.data = ., GT = stringi::stri_join(A1, A2)) else .}
        }
        if (gt.bin == 2L) {
          x %<>%
            {if (gt.vcf) dplyr::mutate(.data = ., GT_VCF = "1/1") else .} %>%
            {if (gt.vcf.nuc) dplyr::mutate(.data = ., GT_VCF_NUC = stringi::stri_join(ALT, ALT, sep = "/")) else .} %>%
            {if (gt) dplyr::mutate(.data = ., GT = stringi::stri_join(A2, A2)) else .}
        }
      }
    }#End GT_BIN


    if (gt.format == "GT_VCF") {
      gt.vcf <- unique(x$GT_VCF)
      if (gt.vcf == "./.") {
        x %<>%
          {if (gt.bin) dplyr::mutate(.data = ., GT_BIN = NA_integer_) else .} %>%
          {if (gt.vcf.nuc) dplyr::mutate(.data = ., GT_VCF_NUC = "./.") else .} %>%
          {if (gt) dplyr::mutate(.data = ., GT = "000000") else .}
      } else {
        if (gt.vcf == "0/0") {
          x %<>%
            {if (gt.bin) dplyr::mutate(.data = ., GT_BIN = 0L) else .} %>%
            {if (gt.vcf.nuc) dplyr::mutate(.data = ., GT_VCF_NUC = stringi::stri_join(REF, REF, sep = "/")) else .} %>%
            {if (gt) dplyr::mutate(.data = ., GT = stringi::stri_join(A1, A1)) else .}
        }
        if (gt.vcf == "1/0") {
          x %<>%
            {if (gt.bin) dplyr::mutate(.data = ., GT_BIN = 1L) else .} %>%
            {if (gt.vcf.nuc) dplyr::mutate(.data = ., GT_VCF_NUC = stringi::stri_join(REF, ALT, sep = "/")) else .} %>%
            {if (gt) dplyr::mutate(.data = ., GT = stringi::stri_join(A1, A2)) else .}
        }
        if (gt.vcf == "0/1") {
          x %<>%
            {if (gt.bin) dplyr::mutate(.data = ., GT_BIN = 1L) else .} %>%
            {if (gt.vcf.nuc) dplyr::mutate(.data = ., GT_VCF_NUC = stringi::stri_join(REF, ALT, sep = "/")) else .} %>%
            {if (gt) dplyr::mutate(.data = ., GT = stringi::stri_join(A1, A2)) else .}
        }
        if (gt.vcf == "1/1") {
          x %<>%
            {if (gt.bin) dplyr::mutate(.data = ., GT_BIN = 2L) else .} %>%
            {if (gt.vcf.nuc) dplyr::mutate(.data = ., GT_VCF_NUC = stringi::stri_join(ALT, ALT, sep = "/")) else .} %>%
            {if (gt) dplyr::mutate(.data = ., GT = stringi::stri_join(A2, A2)) else .}
        }
      }
    }#End GT_VCF

    if (gt.format %in% c("GT", "GT_VCF_NUC")) {
      message("Not implemented yet...")
    }
    return(x)
  }#End gt_map


  # split the data
  x %<>%
    dplyr::group_split(.data[[detect.gt]]) %>%
    purrr::map_dfr(.x = ., .f = gt_map, gt.format = detect.gt, gt = gt, gt.bin = gt.bin, gt.vcf = gt.vcf, gt.vcf.nuc = gt.vcf.nuc)

  if (remove.extra) x  %<>% dplyr::select(-c(A1, A2))
  if (arrange) x  %<>% dplyr::arrange(IDTEMP) %>% dplyr::select(-IDTEMP)
  return(x)
}#End gt_recoding
