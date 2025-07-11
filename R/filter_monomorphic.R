# Filter monomorphic markers

#' @name filter_monomorphic

#' @title Filter monomorphic markers

#' @description Filter monomorphic markers.
#' This filter will remove from the dataset
#' markers with just \emph{one genotype phenotype}:
#' \itemize{
#' \item genotypes are ALL homozygotes REF/REF (pp)
#' \item genotypes are ALL heterozygotes REF/ALT, ALT/REF (pq or qp)
#' \item genotypes are ALL homozygotes ALT/ALT (qq)
#' }
#'
#' \strong{Filter targets}: SNPs
#'
#' \strong{Statistics}: the number of genotype phenotypes
#'
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users who wants to keep only polymorphic markers in
#' their dataset.

#' @inheritParams radiator_common_arguments
#' @param filter.monomorphic (optional, logical)
#' Default: \code{filter.monomorphic = TRUE}.
#' @return A list with the filtered input, whitelist and blacklist of markers..

#' @export
#' @rdname filter_monomorphic
#' @seealso
#' \code{\link{filter_rad}},
#' \code{\link{tidy_genomic_data}}, \code{\link{read_vcf}},
#' \code{\link{tidy_vcf}}.

#' @examples
#' \dontrun{
#' require(SeqArray) # when using gds
#' mono <- radiator::filter_monomorphic(data = "my.radiator.gds.rad", verbose = TRUE)
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

filter_monomorphic <- function(
  data,
  filter.monomorphic = TRUE,
  parallel.core = parallel::detectCores() - 1,
  verbose = FALSE,
  ...) {

  ## Test
  # parallel.core = parallel::detectCores() - 1
  # verbose = FALSE
  # path.folder = NULL
  # parameters = NULL
  # internal = FALSE
  # obj.keeper <- c(ls(envir = globalenv()), "data")

  if (filter.monomorphic) {
    radiator_function_header(f.name = "filter_monomorphic", verbose = verbose)

    # Cleanup---------------------------------------------------------------------
    file.date <- format(Sys.time(), "%Y%m%d@%H%M")
    if (verbose) message("Execution date@time: ", file.date)
    old.dir <- getwd()
    opt.change <- getOption("width")
    options(width = 70)
    timing <- radiator_tic()
    res <- list()
    #back to the original directory and options
    on.exit(setwd(old.dir), add = TRUE)
    on.exit(options(width = opt.change), add = TRUE)
    on.exit(radiator_toc(timing, verbose = verbose), add = TRUE)
    on.exit(radiator_function_header(f.name = "filter_monomorphic", start = FALSE, verbose = verbose), add = TRUE)
    # on.exit(rm(list = setdiff(ls(envir = sys.frame(-1L)), obj.keeper), envir = sys.frame(-1L)))


    # if (verbose) message("\nScanning for monomorphic markers...")
    # message("\nScanning for monomorphic markers...")

    # Checking for missing and/or default arguments ------------------------------
    if (missing(data)) rlang::abort("Input file missing")

    # Function call and dotslist -------------------------------------------------
    rad.dots <- radiator_dots(
      func.name = as.list(sys.call())[[1]],
      fd = rlang::fn_fmls_names(),
      args.list = as.list(environment()),
      dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE),
      keepers = c("path.folder", "parameters", "internal"),
      verbose = FALSE
    )

    # Folders---------------------------------------------------------------------
    path.folder <- generate_folder(
      rad.folder = "filter_monomorphic",
      path.folder = path.folder,
      internal = internal,
      file.date = file.date,
      verbose = verbose)

    # write the dots file
    write_radiator_tsv(
      data = rad.dots,
      path.folder = path.folder,
      filename = "radiator_filter_monomorphic_args",
      date = TRUE,
      internal = internal,
      write.message = "Function call and arguments stored in: ",
      verbose = verbose
    )


    # Detect format --------------------------------------------------------------
    data.type <- radiator::detect_genomic_format(data)
    if (!data.type %in% c("tbl_df", "fst.file", "SeqVarGDSClass", "gds.file")) {
      rlang::abort("Input not supported for this function: read function documentation")
    }

    # GDS
    if (data.type %in% c("SeqVarGDSClass", "gds.file")) {
      radiator_packages_dep(package = "SeqArray", cran = FALSE, bioc = TRUE)

      if (data.type == "gds.file") {
        data <- radiator::read_rad(data, verbose = verbose)
        data.type <- "SeqVarGDSClass"
      }

      # Filter parameter file: generate and initiate -----------------------------
      filters.parameters <- radiator_parameters(
        generate = TRUE,
        initiate = TRUE,
        update = FALSE,
        parameter.obj = parameters,
        data = data,
        path.folder = path.folder,
        file.date = file.date,
        internal = internal,
        verbose = verbose)

      # Scanning for monomorphic markers------------------------------------------
      n.markers.before <- filters.parameters$info$n.snp
      bl <- count_monomorphic(x = data, parallel.core = parallel.core)
      n.markers.removed <- length(bl)
      want <- c("VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS")
      markers.meta <- extract_markers_metadata(gds = data, whitelist = FALSE)

      if (n.markers.removed > 0) {
        n.markers.after <- n.markers.before - n.markers.removed
        markers.meta %<>%
          dplyr::mutate(
            FILTERS = dplyr::if_else(VARIANT_ID %in% bl, "filter.monomorphic", FILTERS)
          )
        write_radiator_tsv(
          data = markers.meta %>% dplyr::filter(FILTERS == "filter.monomorphic"),
          path.folder = path.folder,
          filename = "blacklist.monomorphic.markers",
          date = TRUE,
          internal = internal,
          write.message = "standard",
          verbose = verbose
        )


        # wl %<>% dplyr::filter(!MARKERS %in% bl$MARKERS)


        # Update GDS
        update_radiator_gds(
          gds = data,
          node.name = "markers.meta",
          value = markers.meta,
          sync = TRUE,
          verbose = verbose
        )
      }

      # write the whitelist even if no blacklist...
      write_radiator_tsv(
        data = markers.meta %>% dplyr::filter(FILTERS == "whitelist"),
        path.folder = path.folder,
        filename = "whitelist.polymorphic.markers",
        date = TRUE,
        internal = internal,
        write.message = "standard",
        verbose = verbose
      )
    } else {# tidy data
      # Import data ---------------------------------------------------------------
      if (is.vector(data)) data <- radiator::tidy_wide(data = data, import.metadata = TRUE)
      data.type <- "tbl_df"

      # Keep whitelist and blacklist (same = same space used)
      wl <- radiator::separate_markers(
        data = data,
        sep = NULL,
        markers.meta.lists.only = TRUE,
        generate.markers.metadata = FALSE
      )
      data %<>% dplyr::left_join(wl, by = intersect(colnames(data), colnames(wl)))

      # Filter parameter file: generate and initiate ------------------------------------------
      filters.parameters <- radiator_parameters(
        generate = TRUE,
        initiate = TRUE,
        update = FALSE,
        parameter.obj = parameters,
        data = data,
        path.folder = path.folder,
        file.date = file.date,
        internal = internal,
        verbose = verbose)

      # Scanning for monomorphic markers------------------------------------------
      if (verbose) message("Scanning for monomorphic markers...")
      # n.markers.before <- nrow(wl)

      if (tibble::has_name(data, "GT_BIN")) {
        bl <- dplyr::select(.data = data, MARKERS, GT_BIN) %>%
          dplyr::filter(!is.na(GT_BIN)) %>%
          dplyr::distinct(MARKERS, GT_BIN) %>%
          dplyr::count(x = ., MARKERS) %>%
          dplyr::filter(n == 1) %>%
          dplyr::distinct(MARKERS)
      } else {
        bl <- dplyr::select(.data = data, MARKERS, GT) %>%
          dplyr::filter(GT != "000000") %>%
          dplyr::distinct(MARKERS, GT) %>%
          dplyr::mutate(
            A1 = stringi::stri_sub(GT, 1, 3),
            A2 = stringi::stri_sub(GT, 4,6)
          ) %>%
          dplyr::select(-GT) %>%
          tidyr::pivot_longer(
            data = .,
            cols = -MARKERS,
            names_to = "ALLELES_GROUP",
            values_to = "ALLELES"
          ) %>%
          dplyr::distinct(MARKERS, ALLELES) %>%
          dplyr::count(x = ., MARKERS) %>%
          dplyr::filter(n == 1) %>%
          dplyr::distinct(MARKERS)
      }
      # Remove the markers from the dataset
      n.markers.removed <- nrow(bl)

      if (n.markers.removed > 0) {
        data <- dplyr::filter(data, !MARKERS %in% bl$MARKERS)
        bl <- wl %>% dplyr::filter(MARKERS %in% bl$MARKERS)
        write_radiator_tsv(
          data = bl,
          path.folder = path.folder,
          filename = "blacklist.monomorphic.markers",
          date = TRUE,
          internal = internal,
          write.message = "standard",
          verbose = verbose
        )

      } else {
        bl <- wl[0,]
      }
      # write the whitelist even if no blacklist...
      wl %<>% dplyr::filter(!MARKERS %in% bl$MARKERS)
      write_radiator_tsv(
        data = wl,
        path.folder = path.folder,
        filename = "whitelist.polymorphic.markers",
        date = TRUE,
        internal = internal,
        write.message = "standard",
        verbose = verbose
      )

    }#Tidy data

    # Filter parameter file: update --------------------------------------------
    filters.parameters <- radiator_parameters(
      generate = FALSE,
      initiate = FALSE,
      update = TRUE,
      parameter.obj = filters.parameters,
      data = data,
      filter.name = "Filter monomorphic markers",
      param.name = "filter.monomorphic",
      values = "",
      path.folder = path.folder,
      file.date = file.date,
      internal = internal,
      verbose = verbose)

    # Return -------------------------------------------------------------------
    radiator_results_message(
      rad.message = "\nFilter monomorphic markers",
      filters.parameters,
      internal,
      verbose
    )
  }
  return(data)
}#End filter_monomorphic

#' @title count_monomorphic
#' @description count monomorphe in gds
#' @name count_monomorphic
#' @rdname count_monomorphic
#' @keywords internal
#' @export
count_monomorphic <- function(x, parallel.core = parallel::detectCores() - 1) {
  if (Sys.info()[['sysname']] == "Windows") parallel.core <- 1
  variants <- SeqArray::seqGetData(gdsfile = x, var.name = "variant.id")
  # Function proposed by Xiuwen doesnt work below line 329 it missed some markers
  # in some datasets...
  mono <- SeqArray::seqApply(
    gdsfile = x,
    var.name = "$dosage_alt",
    # FUN = function(g) all(g == 1L, na.rm = TRUE),
    FUN = function(g) length(unique(g[!is.na(g)])) == 1,
    margin = "by.variant", as.is = "logical",
    parallel = parallel.core
  )
  bl <- variants[mono]
  return(bl)
}
