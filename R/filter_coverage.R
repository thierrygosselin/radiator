# Filter markers coverage
#' @name filter_coverage
#' @title Filter markers mean coverage
#' @description This function is designed to remove/blacklist markers
#' based on mean coverage information.
#'
#' \strong{Filter targets}: SNPs
#'
#' \strong{Statistics}: mean coverage
#' ( The read depth of individual genotype is averaged across markers).

#' @param filter.coverage (optional, string) 2 options:
#' \itemize{
#' \item character string \code{filter.coverage = "outliers"} will use as
#' thresholds the lower and higher outlier values in the box plot.
#' \item integers string \code{filter.coverage = c(10, 200)}. For the
#' marker's mean coverage lower and upper bound.
#' }
#' Default: \code{filter.coverage = NULL}.

#' @param filename (optional, character)
#' Default: \code{filename = NULL}.

#' @inheritParams radiator_common_arguments

#' @section Advance mode:
#'
#' \emph{dots-dots-dots ...} allows to pass several arguments for fine-tuning the function:
#' \enumerate{
#' \item \code{filter.common.markers} (optional, logical).
#' Default: \code{filter.common.markers = FALSE},
#' Documented in \code{\link{filter_common_markers}}.
#' \item \code{filter.monomorphic} (logical, optional) Should the monomorphic
#' markers present in the dataset be filtered out ?
#' Default: \code{filter.monomorphic = TRUE}.
#' Documented in \code{\link{filter_monomorphic}}.
#' \item \code{path.folder}: to write ouput in a specific path
#' (used internally in radiator).
#' Default: \code{path.folder = getwd()}.
#' If the supplied directory doesn't exist, it's created.
#' }


#' @section Interactive version:
#'
#' To help choose a threshold for the local and global MAF
#' use the interactive version.
#'
#' 2 steps in the interactive version:
#'
#' Step 1. Visualization and helper table.
#'
#' Step 2. Filtering markers based on mean coverage


#' @rdname filter_coverage
#' @export

#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot stat_summary labs theme element_blank element_text geom_jitter scale_colour_manual scale_y_reverse theme_light geom_bar facet_grid geom_histogram aes_string scale_fill_manual theme_bw stat_smooth geom_boxplot ggsave scale_size_area
#' @importFrom stringi stri_join stri_replace_all_fixed stri_sub
#' @importFrom dplyr select distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join funs summarise_at bind_rows
#' @importFrom readr write_tsv
#' @importFrom tibble data_frame has_name
#' @importFrom tidyr complete gather unite spread nesting
#' @importFrom fst read.fst


#' @return With \code{interactive.filter = FALSE}, a list in the global environment,
#' with 7 objects:
#' \enumerate{
#' \item $tidy.filtered.mac
#' \item $whitelist.markers
#' \item $blacklist.markers
#' \item $mac.data
#' \item $filters.parameters
#' }
#'
#' With \code{interactive.filter = TRUE}, a list with 4 additionnal objects are generated.
#' \enumerate{
#' \item $distribution.mac.global
#' \item $distribution.mac.local
#' \item $mac.global.summary
#' \item $mac.helper.table
#' }

#' @examples
#' \dontrun{

#' # The minumum
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

filter_coverage <- function(
  interactive.filter = TRUE,
  data,
  filter.coverage = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
  ...
) {

  if (!is.null(filter.coverage) || interactive.filter) {

    # interactive.filter = TRUE
    # data = gds
    # filter.coverage = c(5, 150)
    # filter.coverage = NULL
    # filename = NULL
    # parallel.core = parallel::detectCores() - 1
    # verbose = TRUE
    # path.folder = "coverage_test2"
    # parameters <- NULL
    # force.stats <- NULL


    if (interactive.filter) verbose <- TRUE
    if (verbose) {
      cat("################################################################################\n")
      cat("########################### radiator::filter_coverage ##########################\n")
      cat("################################################################################\n")
    }
    if (!verbose) message("filter_coverage...")

    # Cleanup-------------------------------------------------------------------
    file.date <- format(Sys.time(), "%Y%m%d@%H%M")
    if (verbose) message("Execution date/time: ", file.date)
    old.dir <- getwd()
    opt.change <- getOption("width")
    options(width = 70)
    timing <- proc.time()# for timing
    # res <- list()
    #back to the original directory and options
    on.exit(setwd(old.dir), add = TRUE)
    on.exit(options(width = opt.change), add = TRUE)
    on.exit(timing <- proc.time() - timing, add = TRUE)
    on.exit(if (verbose) message("\nComputation time, overall: ", round(timing[[3]]), " sec"), add = TRUE)
    on.exit(if (verbose) cat("########################## completed filter_coverage ###########################\n"), add = TRUE)

    # Function call and dotslist -----------------------------------------------
    rad.dots <- radiator_dots(
      func.name = as.list(sys.call())[[1]],
      fd = rlang::fn_fmls_names(),
      args.list = as.list(environment()),
      dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE),
      keepers = c("path.folder", "parameters", "force.stats", "internal"),
      verbose = verbose
    )

    # Checking for missing and/or default arguments ----------------------------
    if (missing(data)) rlang::abort("data is missing")

    # Folders---------------------------------------------------------------------
    path.folder <- generate_folder(
      f = path.folder,
      rad.folder = "filter_coverage",
      internal = internal,
      file.date = file.date,
      verbose = verbose)

    # write the dots file
    write_rad(
      data = rad.dots,
      path = path.folder,
      filename = stringi::stri_join("radiator_filter_coverage_args_", file.date, ".tsv"),
      tsv = TRUE,
      internal = internal,
      verbose = verbose
    )

    # Message about steps taken during the process -----------------------------
    if (interactive.filter) {
      message("Interactive mode: on\n")
      message("Step 1. Visualization and helper table")
      message("Step 2. Filtering markers based on total coverage\n\n")
    }

    # Detect format --------------------------------------------------------------
    data.type <- radiator::detect_genomic_format(data)
    if (!data.type %in% c("SeqVarGDSClass", "gds.file")) {
      rlang::abort("Input not supported for this function: read function documentation")
    }

    # Import data ---------------------------------------------------------------
    if (verbose) message("Importing data ...")
    if (!"SeqVarTools" %in% utils::installed.packages()[,"Package"]) {
      rlang::abort('Please install SeqVarTools for this option:\n
             install.packages("BiocManager")
             BiocManager::install("SeqVarTools")')
    }

    if (data.type == "gds.file") {
      data <- radiator::read_rad(data, verbose = verbose)
      data.type <- "SeqVarGDSClass"
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

    # Step 1. Visuals ----------------------------------------------------------
    if (interactive.filter) message("\nStep 1. Coverage visualization and helper table\n")

    # Generate coverage stats---------------------------------------------------
    if (verbose) message("Generating coverage statistics")
    info <- generate_markers_stats(
      gds = data,
      missing = FALSE,
      heterozygosity = FALSE,
      snp.per.locus = FALSE,
      snp.position.read = FALSE,
      force.stats = force.stats,
      path.folder = path.folder,
      file.date = file.date,
      parallel.core = parallel.core
    )
    stats <- info$stats
    info <- info$info

    # Whitelist and blacklist --------------------------------------------------
    # want <- c("VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS")
    wl <- bl <- extract_markers_metadata(gds = data)

    # identify outliers: low and high -----------------------------------------
    out.low <- floor(stats$OUTLIERS_LOW[stats$GROUP == "mean coverage"]*1000)/1000
    out.high <- floor(stats$OUTLIERS_HIGH[stats$GROUP == "mean coverage"]*1000)/1000
    if (verbose) message("Generating coverage statistics: without outliers")
    variant.wo.out <- dplyr::filter(info, COVERAGE_MEAN >= out.low &
                                      COVERAGE_MEAN <= out.high) %$%
      VARIANT_ID %>%
      as.integer

    # Helper table -------------------------------------------------------------
    # filter.coverage <- "outliers"
    # filter.coverage <- c(10, 100)
    if (!is.null(filter.coverage)) {
      if (length(filter.coverage) > 1) {
        combined.info <- TRUE
        cov.low <- filter.coverage[1]
        cov.high <- filter.coverage[2]
        # if (verbose) message("\nRemoving markers based on mean coverage statistics: ", cov.low, " / ", cov.high)
      } else {
        if (is.character(filter.coverage)) {
          combined.info <- FALSE
          cov.low <- out.low
          cov.high <- out.high
          filter.coverage <- c(cov.low, cov.high)
        } else {
          rlang::abort("Unknown mean coverage thresholds used")
        }
      }
    } else {
      combined.info <- FALSE
    }

    # combined.info <- TRUE
    # cov.low <- 10
    # cov.high <- 100
    min.c <- floor(stats$MIN[stats$GROUP == "mean coverage"]*1000)/1000
    max.c <- floor(stats$MAX[stats$GROUP == "mean coverage"]*1000)/1000

    if (combined.info) {
      # include outliers
      max.out <- floor(stats$OUTLIERS_HIGH[stats$GROUP == "mean coverage"]*1000)/1000

      if (out.low != cov.low) {
        cl.range <- min(out.low,cov.low):max(out.low,cov.low)
      } else {
        cl.range <- out.low
      }

      if (out.high != cov.high) {
        ch.range <- min(out.high,cov.high):max(out.high,cov.high)
      } else {
        ch.range <- out.high
      }

    } else {
      # use min and max values with outliers
      #MIN
      if (min.c != out.low) {
        cl.range <- out.low:min.c
      } else {
        cl.range <- min.c
      }

      #MAX
      if (max.c != out.high) {
        ch.range <- out.high:max.c
      } else {
        ch.range <- max.c
      }
    }

    # reduce size of points...
    if (length(cl.range) > 100) cl.range <- seq(min(cl.range), max(cl.range), by = 10)
    if (length(ch.range) > 100) ch.range <- seq(min(ch.range), max(ch.range), by = 10)


    if (verbose) message("Generating mean coverage helper table...")
    how_many_markers <- function(threshold, x, low = FALSE) {
      if (low) {
        nrow(dplyr::filter(x, COVERAGE_MEAN >= threshold))
      } else {
        nrow(dplyr::filter(x, COVERAGE_MEAN <= threshold))
      }
    }#End how_many_markers

    n.markers <- nrow(info)
    helper.table.low <- tibble::tibble(COVERAGE_LOW = cl.range) %>%
      dplyr::mutate(
        WHITELISTED_MARKERS = purrr::map_int(.x = cl.range, .f = how_many_markers, x = info, low = TRUE),
        BLACKLISTED_MARKERS = n.markers - WHITELISTED_MARKERS
      ) %>%
      readr::write_tsv(
        x = .,
        path = file.path(path.folder, "coverage.low.helper.table.tsv"))

    if (nrow(helper.table.low) > 1) {
      markers.plot.low <- ggplot2::ggplot(
        data = tidyr::gather(
          data = helper.table.low,
          key = LIST, value = MARKERS, -COVERAGE_LOW),
        ggplot2::aes(x = COVERAGE_LOW, y = MARKERS)) +
        ggplot2::geom_line() +
        ggplot2::geom_point(size = 2, shape = 21, fill = "white") +
        ggplot2::scale_x_continuous(name = "Minimum mean coverage allowed", breaks = cl.range) +
        ggplot2::scale_y_continuous(name = "Number of markers")+
        ggplot2::theme(
          axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
          axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
          # axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica") #angle = 90, hjust = 1, vjust = 0.5),
          axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5)
        ) +
        ggplot2::theme_bw()+
        ggplot2::facet_grid(LIST ~ ., scales = "free", space = "free")

      print(markers.plot.low)
      # save
      ggplot2::ggsave(
        filename = file.path(path.folder, "coverage.low.helper.plot.pdf"),
        plot = markers.plot.low, width = 20, height = 15, dpi = 300, units = "cm", useDingbats = FALSE)
    }

    helper.table.high <- tibble::tibble(COVERAGE_HIGH = ch.range) %>%
      dplyr::mutate(
        WHITELISTED_MARKERS = purrr::map_int(.x = ch.range, .f = how_many_markers, x = info),
        BLACKLISTED_MARKERS = n.markers - WHITELISTED_MARKERS
      ) %>%
      readr::write_tsv(
        x = .,
        path = file.path(path.folder, "coverage.high.helper.table.tsv"))

    if (nrow(helper.table.high) > 1) {
      markers.plot.high <- ggplot2::ggplot(
        data = tidyr::gather(
          data = helper.table.high,
          key = LIST, value = MARKERS, -COVERAGE_HIGH),
        ggplot2::aes(x = COVERAGE_HIGH, y = MARKERS)) +
        ggplot2::geom_line() +
        ggplot2::geom_point(size = 2, shape = 21, fill = "white") +
        ggplot2::scale_x_continuous(name = "Maximum mean coverage allowed", breaks = ch.range) +
        ggplot2::scale_y_continuous(name = "Number of markers")+
        ggplot2::theme(
          axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
          axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
          axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica") #angle = 90, hjust = 1, vjust = 0.5),
          # strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
        ) +
        ggplot2::theme_bw()+
        ggplot2::facet_grid(LIST ~ ., scales = "free", space = "free")

      print(markers.plot.high)
      # save
      ggplot2::ggsave(
        filename = file.path(path.folder, "coverage.high.helper.plot.pdf"),
        plot = markers.plot.high, width = 20, height = 15, dpi = 300, units = "cm", useDingbats = FALSE)
    }
    if (verbose) message("Files written: helper tables and plots")


    # Step 2. Thresholds selection ---------------------------------------------
    if (interactive.filter) {
      filter.coverage <- c(min.c, max.c)
      message("\nStep 2. Filtering markers based on mean coverage\n")
      filter.coverage[1] <- interactive_question(
        x = "Choose the min mean coverage threshold: ", minmax = c(1, max.c))
    }
    if (interactive.filter) {
      filter.coverage[2] <- interactive_question(
        x = "Choose the max mean coverage threshold: ", minmax = c(1, max.c))
    }

    # Whitelist and Blacklist of markers
    wl %<>% dplyr::filter(COVERAGE_MEAN >= filter.coverage[1] &
                            COVERAGE_MEAN <= filter.coverage[2]) %>%
      readr::write_tsv(
        x = .,
        path = file.path(path.folder, "whitelist.markers.coverage.tsv"),
        append = FALSE, col_names = TRUE)
    bl %<>% dplyr::setdiff(wl) %>%
      readr::write_tsv(
        x = .,
        path = file.path(path.folder, "blacklist.markers.coverage.tsv"),
        append = FALSE, col_names = TRUE)
    # saving whitelist and blacklist
    if (verbose) message("File written: whitelist.markers.coverage.tsv")
    if (verbose) message("File written: blacklist.markers.coverage.tsv")

    # Filtering ----------------------------------------------------------------
    # Update GDS
    update_radiator_gds(
      gds = data,
      node.name = "markers.meta",
      value = wl,
      sync = TRUE
    )

    # update blacklist.markers
    if (nrow(bl) > 0) {
      bl %<>% dplyr::select(MARKERS) %>%
        dplyr::mutate(FILTER = "filter.mean.coverage")
      bl.gds <- update_bl_markers(gds = data, update = bl)
    }


    # Update parameters --------------------------------------------------------
    filters.parameters <- radiator_parameters(
      generate = FALSE,
      initiate = FALSE,
      update = TRUE,
      parameter.obj = filters.parameters,
      data = data,
      filter.name = "Filter coverage min / max",
      param.name = "filter.coverage",
      values = paste(filter.coverage, collapse = " / "),
      path.folder = path.folder,
      file.date = file.date,
      verbose = verbose)

    # results --------------------------------------------------------------------
    if (verbose) cat("################################### RESULTS ####################################\n")
    message("Filter mean coverage thresholds: ", paste(filter.coverage, collapse = " / "))
    message("Number of individuals / strata / chrom / locus / SNP:")
    if (verbose) message("    Before: ", filters.parameters$filters.parameters$BEFORE)
    message("    Blacklisted: ", filters.parameters$filters.parameters$BLACKLIST)
    if (verbose) message("    After: ", filters.parameters$filters.parameters$AFTER)
  }
  return(data)
}#End filter_coverage
