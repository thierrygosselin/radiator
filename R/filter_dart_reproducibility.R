# SNP number per haplotype
#' @name filter_dart_reproducibility
#' @title Filter data based on DArT reproducibility statistics
#' @description This filter removes markers below a certain threshold.
#' Based on the repoducibility column found in DArT files.
#'
#' \strong{Filter targets}: Markers
#'
#' \strong{Statistics}: Reproducibility (established by DArT)


# Most arguments are inherited from tidy_genomic_data
#' @inheritParams radiator_common_arguments

#' @param filter.reproducibility (double, character) This is best decided after viewing the figures.
#' Usually values higher than 0.95 are not uncommon.
#' The value can also be character: \code{filter.reproducibility = "outliers"}.
#' Using this, will remove outlier markers using the lower outlier statistics.
#' Default: \code{filter.reproducibility = NULL}.


#' @rdname filter_dart_reproducibility
#' @export

#' @details
#' \strong{Interactive version}
#'
#' There are 2 steps in the interactive version to visualize and filter
#' the data based on the reproducibility value:
#'
#' Step 1. Visualization using a box plot
#'
#' Step 2. Choose the filtering threshold
#'
#'
#' @return A list in the global environment with 6 objects:
#' \enumerate{
#' \item $whitelist.markers
#' \item $blacklist.markers
#' \item $filters.parameters
#' }
#'
#' The object can be isolated in separate object outside the list by
#' following the example below.

#' @examples
#' \dontrun{
#' spotted.cod <- radiator::tidy_dart(
#'     data = "Combined_1514and1614_SNP_80Callrate.csv",
#'     strata = "strata.dart.spotted.cod.tsv"
#' )
#' turtle.filtered <- radiator::filter_dart_reproducibility(
#' data = spotted.cod
#' filter.reproducibility = 0.97
#' )
#' }

filter_dart_reproducibility <- function(
  interactive.filter = TRUE,
  data,
  filter.reproducibility = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
  ...
) {

  # interactive.filter <- TRUE
  # data <- gds
  # path.folder <- "testing_reproducility"
  # force.stats <- TRUE
  # parameters <- NULL
  # filename <- NULL
  # filter.reproducibility <- NULL
  # parallel.core <- parallel::detectCores() - 1
  # verbose = TRUE
  if (interactive.filter || !is.null(filter.reproducibility)) {
    if (interactive.filter) verbose <- TRUE
    if (verbose) {
      cat("################################################################################\n")
      cat("##################### radiator::filter_dart_reproducibility ####################\n")
      cat("################################################################################\n")
    }
    # Cleanup-------------------------------------------------------------------
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
    on.exit(if (verbose) cat("#################### completed filter_dart_reproducibility #####################\n"), add = TRUE)

    # Function call and dotslist -------------------------------------------------
    rad.dots <- radiator_dots(
      func.name = as.list(sys.call())[[1]],
      fd = rlang::fn_fmls_names(),
      args.list = as.list(environment()),
      dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE),
      keepers = c("path.folder", "parameters", "internal"),
      verbose = verbose
    )

    # Checking for missing and/or default arguments ------------------------------
    if (missing(data)) rlang::abort("data is missing")

    # Folders---------------------------------------------------------------------
    path.folder <- generate_folder(
      f = path.folder,
      rad.folder = "filter_dart_reproducibility",
      internal = internal,
      file.date = file.date,
      verbose = verbose)

    # write the dots file
    write_rad(
      data = rad.dots,
      path = path.folder,
      filename = stringi::stri_join("radiator_filter_dart_reproducibility_args_", file.date, ".tsv"),
      tsv = TRUE,
      internal = internal,
      verbose = verbose
    )

    # interactive.filter steps -------------------------------------------------
    if (interactive.filter) {
      message("\nInteractive mode: on")
      message("2 steps to visualize and filter the data based on reproducibility:")
      message("Step 1. Visualization")
      message("Step 2. Choose the filtering threshold\n\n")
    }

    # File type detection----------------------------------------------------------
    data.type <- radiator::detect_genomic_format(data)
    if (!data.type %in% c("tbl_df", "fst.file", "SeqVarGDSClass", "gds.file")) {
      rlang::abort("Input not supported for this function: read function documentation")
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
    } else {
      if (is.vector(data)) data <- radiator::tidy_wide(data = data, import.metadata = TRUE)
      data.type <- "tbl_df"
    }

    # Filter parameter file: initiate ------------------------------------------
    filters.parameters <- radiator_parameters(
      generate = TRUE,
      initiate = TRUE,
      update = FALSE,
      parameter.obj = parameters,
      data = data,
      path.folder = path.folder,
      file.date = file.date,
      verbose = verbose)

    # Whitelist and blacklist --------------------------------------------------
    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "COL", "REP_AVG")
    if (data.type == "SeqVarGDSClass") {
      wl <- bl <- extract_markers_metadata(gds = data)
    } else {
      wl <- bl <- suppressWarnings(dplyr::select(data, dplyr::one_of(want)))
    }
    # Check that required info is present in data:
    check <- tibble::has_name(wl, "REP_AVG")
    if (!check) {
      message("This filter requires REP_AVG info, skipping filtering...")
      bl %<>% dplyr::setdiff(wl)
      return(res = list(input = data,
                        whitelist.markers.snp.number = wl,
                        blacklist.markers.snp.number = bl))
    }

    # Generating statistics ----------------------------------------------------
    rep.stats <- tibble_stats(x = wl$REP_AVG, group = "reproducibility")
    readr::write_tsv(x = rep.stats, path = file.path(path.folder, "dart_reproducibility_stats.tsv"))
    if (verbose) message("File written: dart_reproducibility_stats.tsv")

    # Generate box plot ---------------------------------------------------------
    outlier.rep <- ceiling(rep.stats[[8]] * 1000) / 1000
    bp.filename <- stringi::stri_join("dart_reproducibility_boxplot_", file.date, ".pdf")
    rep.fig <- boxplot_stats(
      data = rep.stats,
      title =  "DArT marker's reproducibility",
      subtitle = stringi::stri_join("\nLower outlier: ", outlier.rep),
      x.axis.title = NULL,
      y.axis.title = "Reproducibility (proportion)",
      bp.filename = bp.filename,
      path.folder = path.folder)
    if (verbose) message("File written: ", bp.filename)

    # Helper table -------------------------------------------------------------
    if (verbose) message("Generating helper table...")
    how_many_markers <- function(threshold, x) {
      nrow(dplyr::filter(x, REP_AVG >= threshold))
    }#End how_many_markers

    snp.range <- seq(0.9, 1, 0.005)
    n.markers <- nrow(wl)
    helper.table <- tibble::tibble(REPRODUCIBILITY = snp.range) %>%
      dplyr::mutate(
        WHITELISTED_MARKERS = purrr::map_int(.x = snp.range, .f = how_many_markers, x = wl),
        BLACKLISTED_MARKERS = n.markers - WHITELISTED_MARKERS
      ) %>%
      readr::write_tsv(
        x = .,
        path = file.path(path.folder, "dart.reproducibility.helper.table.tsv"))

    # figures
    markers.plot <- ggplot2::ggplot(
      data = helper.table  %<>% tidyr::gather(
        data = .,
        key = LIST, value = MARKERS, -REPRODUCIBILITY),
      ggplot2::aes(x = REPRODUCIBILITY, y = MARKERS)) +
      ggplot2::geom_line() +
      ggplot2::geom_point(size = 2, shape = 21, fill = "white") +
      ggplot2::geom_vline(ggplot2::aes(xintercept = as.numeric(outlier.rep)), color = "yellow") +
      ggplot2::scale_x_continuous(name = "Reproducibility (proportion)", breaks = snp.range) +
      ggplot2::scale_y_continuous(name = "Number of markers") +
      ggplot2::theme_bw()+
      ggplot2::theme(
        axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica")#, angle = 90, hjust = 1, vjust = 0.5)
      ) +
      ggplot2::facet_grid(LIST ~. , scales = "free", space = "free")
    print(markers.plot)

    # save
    ggplot2::ggsave(
      filename = file.path(path.folder, "dart.reproducibility.helper.plot.pdf"),
      plot = markers.plot,
      width = 20,
      height = 15,
      dpi = 300,
      units = "cm",
      useDingbats = FALSE)
    helper.table <- markers.plot <- NULL
    if (verbose) message("Files written: helper tables and plots")


    # Step 2. Thresholds selection ---------------------------------------------
    if (interactive.filter) {
      message("\nStep 2. Filtering markers based on markers reproducibility\n")

      filter.reproducibility <- interactive_question(
        x = "Do you still want to blacklist markers? (y/n):",
        answer.opt = c("y", "n"))

      if (filter.reproducibility == "y") {
        outlier.stats <- interactive_question(
          x = "Do you want to remove markers based on the outlier statistics or not (y/n) ?
(n: next question will be to enter your own threshold)", answer.opt = c("y", "n"))
        if (outlier.stats == "y") {
          filter.reproducibility <- "outliers"
        } else {
          filter.reproducibility <- interactive_question(
            x = "Enter the proportion threshold (0-1)
The minimum reproducibility tolerated:", minmax = c(0, 1))
        }
        outlier.stats <- NULL
      } else {
        filter.reproducibility <- NULL
      }
    }


    # Filtering ----------------------------------------------------------------
    if (!is.null(filter.reproducibility)) {
      if (!purrr::is_double(filter.reproducibility)) {
        message("\nRemoving outlier markers based on reproducibility statistic: ", outlier.rep)
        filter.reproducibility <- outlier.rep
      } else {
        message("\nRemoving markers based on reproducibility statistic: ", filter.reproducibility)
      }

      # Whitelist and Blacklist of markers
      wl %<>% dplyr::filter(REP_AVG >= filter.reproducibility) %>%
        readr::write_tsv(
          x = .,
          path = file.path(path.folder, "whitelist.dart.reproducibility.tsv"),
          append = FALSE, col_names = TRUE)
      bl %<>% dplyr::setdiff(wl) %>%
        readr::write_tsv(
          x = .,
          path = file.path(path.folder, "blacklist.dart.reproducibility.tsv"),
          append = FALSE, col_names = TRUE)
      # saving whitelist and blacklist
      if (verbose) message("File written: whitelist.dart.reproducibility.tsv")
      if (verbose) message("File written: blacklist.dart.reproducibility.tsv")

      if (data.type == "SeqVarGDSClass") {
        # Update GDS
        update_radiator_gds(
          gds = data,
          node.name = "markers.meta",
          value = wl,
          sync = TRUE
        )
        # sync_gds(gds = data, markers = wl$VARIANT_ID)
        # radiator.gds <- gdsfmt::index.gdsn(
        #   node = data, path = "radiator", silent = TRUE)
        #
        # # Update metadata
        # gdsfmt::add.gdsn(
        #   node = radiator.gds,
        #   name = "markers.meta",
        #   val = wl,
        #   replace = TRUE,
        #   compress = "ZIP_RA",
        #   closezip = TRUE)

        # update blacklist.markers
        if (nrow(bl) > 0) {
          bl %<>% dplyr::select(MARKERS) %>%
            dplyr::mutate(FILTER = "filter.dart.reproducibility")
          bl.gds <- update_bl_markers(gds = data, update = bl)
        }
      } else {
        # Apply the filter to the tidy data
        data  %<>% dplyr::filter(MARKERS %in% wl$MARKERS)
      }

      # Update parameters --------------------------------------------------------
      filters.parameters <- radiator_parameters(
        generate = FALSE,
        initiate = FALSE,
        update = TRUE,
        parameter.obj = filters.parameters,
        data = data,
        filter.name = "Filter DArT reproducibility",
        param.name = "filter.reproducibility",
        values = filter.reproducibility,
        path.folder = path.folder,
        file.date = file.date,
        verbose = verbose)

      # Return -----------------------------------------------------------------------
      if (verbose) {
        cat("################################### RESULTS ####################################\n")
        message("Filter DArT reproducibility: ", filter.reproducibility)
        message("Number of individuals / strata / chrom / locus / SNP:")
        message("    Before: ", filters.parameters$filters.parameters$BEFORE)
        message("    Blacklisted: ", filters.parameters$filters.parameters$BLACKLIST)
        message("    After: ", filters.parameters$filters.parameters$AFTER)
      }
    }
    return(data)
  } else {
    return(data)
  }
} #End filter_dart_reproducibility
