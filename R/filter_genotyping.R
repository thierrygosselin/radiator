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
  data,
  interactive.filter = TRUE,
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
    # path.folder = NULL
    # parameters <- NULL
    # force.stats <- TRUE


    if (interactive.filter) verbose <- TRUE
    if (verbose) {
      cat("################################################################################\n")
      cat("######################### radiator::filter_genotyping ##########################\n")
      cat("################################################################################\n")
    }
    # Cleanup---------------------------------------------------------------------
    file.date <- format(Sys.time(), "%Y%m%d@%H%M")
    if (verbose) message("Execution date@time: ", file.date)
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
      verbose = FALSE
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
      write.message = "Function call and arguments stored in: ",
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
      internal = internal,
      verbose = verbose)

    # Step 1. Visuals ----------------------------------------------------------
    if (interactive.filter) message("\nStep 1. Missing visualization and helper table\n")

    # Generate coverage stats---------------------------------------------------
    if (verbose) message("Generating statistics")
    # return(data)
    # summary_gds(data)
    info <- generate_markers_stats(
      gds = data,
      missing = TRUE,
      coverage = FALSE,
      allele.coverage = FALSE,
      mac = FALSE,
      heterozygosity = FALSE,
      snp.per.locus = FALSE,
      snp.position.read = FALSE,
      force.stats = force.stats,
      path.folder = path.folder,
      file.date = file.date,
      parallel.core = parallel.core
    )
    # summary_gds(data)

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
      )

    readr::write_tsv(
      x = helper.table,
      path = file.path(path.folder, "genotyping.helper.table.tsv")
    )

    # checking if strata present
    strata <- extract_individuals_metadata(
      gds = data,
      ind.field.select = c("INDIVIDUALS", "STRATA"),
      whitelist = TRUE)
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
      helper.table  %<>%
        tidyr::gather(
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
      helper.table  %<>%
        tidyr::gather(
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
      filter.genotyping <- radiator_question(
        x = "Choose the maximum missing proportion allowed: ", minmax = c(0, 1))
    }

    # identify outliers: low and high -----------------------------------------
    if (!purrr::is_double(filter.genotyping)) {
      out.high <- floor(stats$OUTLIERS_HIGH[stats$GROUP == "missing genotypes"]*1000)/1000
      if (verbose) message("\nRemoving outliers markers based on genotyping statistic: ", out.high)
      filter.genotyping <- out.high
    } else {
      if (verbose) message("\nRemoving markers based on genotyping statistic: ", filter.genotyping)
    }

    # Whitelist and blacklist --------------------------------------------------
    bl <- info %>%
      dplyr::filter(MISSING_PROP > filter.genotyping) %$%
      VARIANT_ID
    markers.meta <- extract_markers_metadata(gds = data, whitelist = FALSE) %>%
      dplyr::mutate(
        FILTERS = dplyr::if_else(VARIANT_ID %in% bl, "filter.genotyping", FILTERS)
      )
    # Update GDS
    update_radiator_gds(
      gds = data,
      node.name = "markers.meta",
      value = markers.meta,
      sync = TRUE
    )
    write_rad(
      data = markers.meta %>% dplyr::filter(FILTERS == "filter.genotyping"),
      path = path.folder,
      filename = stringi::stri_join("blacklist.markers.genotyping_", file.date, ".tsv"),
      tsv = TRUE, internal = internal, verbose = verbose)

    write_rad(
      data = markers.meta %>% dplyr::filter(FILTERS == "whitelist"),
      path = path.folder,
      filename = stringi::stri_join("whitelist.markers.genotyping_", file.date, ".tsv"),
      tsv = TRUE, internal = internal, verbose = verbose)

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
      internal = internal,
      verbose = verbose)


    # results --------------------------------------------------------------------
    radiator_results_message(
      rad.message = stringi::stri_join("\nFilter genotyping threshold: ", filter.genotyping),
      filters.parameters,
      internal,
      verbose
    )

    }
  return(data)
}#End filter_genotyping
