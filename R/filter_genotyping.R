# Filter markers genotyping
#' @name filter_genotyping
#' @title Filter markers based on genotyping/missing rate
#' @description This function is designed to remove/blacklist markers
#' based on genotyping/call rate.
#'
#' \strong{Filter targets}: SNPs
#'
#' \strong{Statistics}: mean genotyping/call rate (missingness information).
#'
#' @param interactive.filter (optional, logical) Do you want the filtering session to
#' be interactive. With default: \code{interactive.filter == TRUE}, figures and
#' tables are shown before making decisions for filtering.

#' @param filter.genotyping (optional, string) 2 options:
#' \itemize{
#' \item character string \code{filter.genotyping = "outliers"} will use as
#' thresholds the higher outlier values in the box plot.
#' \item double \code{filter.genotyping = 0.2}. Will allow up to 0.2 missing genotypes.
#' }
#' Default: \code{filter.genotyping = NULL}.

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
#' To help choose a threshold use the interactive version.
#'
#' 2 steps in the interactive version:
#'
#' Step 1. Visualization and helper table.
#'
#' Step 2. Filtering markers based on mean genotyping/missing rate


#' @rdname filter_genotyping
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

filter_genotyping <- function(
  interactive.filter = TRUE,
  data,
  filter.genotyping = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
  ...
) {

  if (interactive.filter || !is.null(filter.genotyping)) {

    # interactive.filter = TRUE
    # data = gds
    # filter.genotyping = 0.2
    # filename = NULL
    # parallel.core = parallel::detectCores() - 1
    # verbose = TRUE
    # path.folder = "genotyping_test"
    # parameters <- NULL
    # force.stats <- NULL


    if (interactive.filter) verbose <- TRUE
    if (verbose) {
      cat("################################################################################\n")
      cat("######################### radiator::filter_genotyping ##########################\n")
      cat("################################################################################\n")
    }
    # Cleanup---------------------------------------------------------------------
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
    on.exit(if (verbose) cat("######################## completed filter_genotyping ###########################\n"), add = TRUE)

    # Function call and dotslist -------------------------------------------------
    rad.dots <- radiator_dots(
      func.name = as.list(sys.call())[[1]],
      fd = rlang::fn_fmls_names(),
      args.list = as.list(environment()),
      dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE),
      keepers = c("path.folder", "parameters", "force.stats", "internal"),
      verbose = verbose
    )

    # Checking for missing and/or default arguments ------------------------------
    if (missing(data)) rlang::abort("data is missing")

    # Folders---------------------------------------------------------------------
    path.folder <- generate_folder(
      f = path.folder,
      rad.folder = "filter_genotyping",
      internal = internal,
      file.date = file.date,
      verbose = verbose)

    # write the dots file
    write_rad(
      data = rad.dots,
      path = path.folder,
      filename = stringi::stri_join("radiator_filter_genotyping_args_", file.date, ".tsv"),
      tsv = TRUE,
      internal = internal,
      verbose = verbose
    )
    # Message about steps taken during the process ---------------------------------
    if (interactive.filter) {
      message("Interactive mode: on\n")
      message("Step 1. Visualization and helper table")
      message("Step 2. Filtering markers based on maximum missing proportion allowed\n\n")
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
    if (interactive.filter) message("\nStep 1. Missing visualization and helper table\n")

    # Whitelist and blacklist --------------------------------------------------
    # want <- c("VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS")
    wl <- bl <- extract_markers_metadata(gds = data)

    # Generate coverage stats---------------------------------------------------
    if (verbose) message("Generating statistics")
    info <- generate_markers_stats(
      gds = data,
      missing = TRUE,
      coverage = FALSE, allele.coverage = FALSE, mac = FALSE,
      heterozygosity = FALSE, snp.per.locus = FALSE, snp.position.read = FALSE,
      force.stats = force.stats,
      path.folder = path.folder,
      file.date = file.date,
      parallel.core = parallel.core
    )
    stats <- info$stats
    info <- info$info

    # Helper table -------------------------------------------------------------
    if (verbose) message("Generating missingness/genotyping helper table...")
    mis_many_markers <- function(threshold, x) {
      nrow(dplyr::filter(x, MISSING_PROP <= threshold))
    }#End how_many_markers

    n.markers <- nrow(info)
    helper.table <- tibble::tibble(MISSING_PROP = seq(0, 1, by = 0.1)) %>%
      dplyr::mutate(
        WHITELISTED_MARKERS = purrr::map_int(
          .x = seq(0, 1, by = 0.1),
          .f = mis_many_markers,
          x = info),
        BLACKLISTED_MARKERS = n.markers - WHITELISTED_MARKERS
      ) %>%
      readr::write_tsv(
        x = .,
        path = file.path(path.folder, "genotyping.helper.table.tsv"))

    # checking if strata present
    strata <- extract_individuals(gds = data, ind.field.select = c("INDIVIDUALS", "STRATA"))
    if (!is.null(strata$STRATA)) {
      m.strata <- missing_per_pop(
        gds = data, strata = strata, parallel.core = parallel.core)

      round_mean <- function(x) as.integer(round(mean(x, na.rm = TRUE), 0))
      mean.pop <- m.strata %>% dplyr::group_by(MISSING_PROP) %>%
        dplyr::summarise_if(.tbl = ., .predicate = is.integer, .funs = round_mean) %>%
        dplyr::mutate(STRATA = "MEAN_POP")
      if (is.factor(strata$STRATA)) {
        strata.pop <- levels(strata$STRATA)
      } else {
        strata.pop <- unique(strata$STRATA)
      }

      strata.levels <- c(strata.pop, "MEAN_POP", "OVERALL")
      n.pop <- length(strata.levels)

      suppressWarnings(
        helper.table %<>%
          dplyr::mutate(STRATA = "OVERALL") %>%
          dplyr::bind_rows(mean.pop, m.strata) %>%
          dplyr::mutate(STRATA = factor(STRATA, levels = strata.levels, ordered = TRUE)) %>%
          dplyr::arrange(STRATA) %>%
          readr::write_tsv(
            x = .,
            path = file.path(path.folder, "markers.pop.missing.helper.table.tsv")))
      m.strata <- round_mean <- mean.pop <- NULL

      if (verbose) message("File written: markers.pop.missing.helper.table.tsv")
      helper.table  %<>% tidyr::gather(
        data = .,
        key = LIST,
        value = MARKERS,
        -c(MISSING_PROP, STRATA)
        ) %>%
        dplyr::mutate(STRATA = factor(STRATA, levels = strata.levels, ordered = TRUE)) %>%
        dplyr::arrange(STRATA)

      strata <- TRUE
      strata.levels <- NULL
    } else {
      helper.table  %<>% tidyr::gather(
        data = .,
        key = LIST, value = MARKERS, -MISSING_PROP)
      n.pop <- 1L
      strata <- FALSE
    }

    # figures
    markers.plot <- ggplot2::ggplot(
      data = helper.table,
      ggplot2::aes(x = MISSING_PROP, y = MARKERS)) +
      ggplot2::geom_line() +
      ggplot2::geom_point(size = 2, shape = 21, fill = "white") +
      ggplot2::scale_x_continuous(name = "Maximum missing proportion allowed", breaks = seq(0, 1, by = 0.1)) +
      ggplot2::scale_y_continuous(name = "Number of markers")+
      ggplot2::theme_bw()+
      ggplot2::theme(
        axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5)
      )

    if (strata) {
      markers.plot <- markers.plot + ggplot2::facet_grid(LIST ~ STRATA, scales = "free", space = "free")
    } else {
      markers.plot <- markers.plot + ggplot2::facet_grid(LIST ~. , scales = "free", space = "free")
    }
    print(markers.plot)

    # save
    ggplot2::ggsave(
      filename = file.path(path.folder, "markers.genotyping.helper.plot.pdf"),
      plot = markers.plot,
      width = 20 + 3 * n.pop,
      height = 15,
      dpi = 300,
      units = "cm",
      useDingbats = FALSE)
    helper.table <- markers.plot <- NULL
    if (verbose) message("Files written: helper tables and plots")


    # Step 2. Thresholds selection ---------------------------------------------
    if (interactive.filter) {
      filter.genotyping <- 1L
      if (verbose) message("\nStep 2. Filtering markers based on maximum missing proportion\n")
      filter.genotyping <- interactive_question(
        x = "Choose the maximum missing proportion allowed: ", minmax = c(0, 1))
    }

    # identify outliers: low and high -----------------------------------------
    if (!purrr::is_double(filter.genotyping)) {
      out.high <- floor(stats$OUTLIERS_HIGH[stats$GROUP == "missing genotypes"]*1000)/1000
      message("\nRemoving outliers markers based on genotyping statistic: ", out.high)
      filter.genotyping <- out.high
    } else {
      message("\nRemoving markers based on genotyping statistic: ", filter.genotyping)
    }

    # Whitelist and Blacklist of markers
    wl %<>% dplyr::filter(MISSING_PROP <= filter.genotyping) %>%
      readr::write_tsv(
        x = .,
        path = file.path(path.folder, "whitelist.markers.genotyping.tsv"),
        append = FALSE, col_names = TRUE)
    bl %<>% dplyr::setdiff(wl) %>%
      readr::write_tsv(
        x = .,
        path = file.path(path.folder, "blacklist.markers.genotyping.tsv"),
        append = FALSE, col_names = TRUE)
    # saving whitelist and blacklist
    if (verbose) message("File written: whitelist.markers.genotyping.tsv")
    if (verbose) message("File written: blacklist.markers.genotyping.tsv")

    # Filtering ----------------------------------------------------------------
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
        dplyr::mutate(FILTER = "filter.genotyping")
      bl.gds <- update_bl_markers(gds = data, update = bl)
    }


    # Update parameters --------------------------------------------------------
    filters.parameters <- radiator_parameters(
      generate = FALSE,
      initiate = FALSE,
      update = TRUE,
      parameter.obj = filters.parameters,
      data = data,
      filter.name = "Filter genotyping",
      param.name = "filter.genotyping",
      values = filter.genotyping,
      path.folder = path.folder,
      file.date = file.date,
      verbose = verbose)

    # if (filters.parameters$filters.parameters$BLACKLIST == 0) {
    #   file.remove(path.folder)
    #   if (verbose) message("Folder removed: ", folder_short(path.folder))
    # }

    # results --------------------------------------------------------------------
    if (verbose) cat("################################### RESULTS ####################################\n")
    if (verbose) message("Filter genotyping threshold: ", filter.genotyping)
    message("Number of individuals / strata / chrom / locus / SNP:")
    message("    Before: ", filters.parameters$filters.parameters$BEFORE)
    message("    Blacklisted: ", filters.parameters$filters.parameters$BLACKLIST)
    message("    After: ", filters.parameters$filters.parameters$AFTER)
  }
  return(data)
}#End filter_genotyping
