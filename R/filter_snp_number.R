# SNP number per haplotype
#' @name filter_snp_number
#' @title Filter SNP number per locus/read
#' @description This filter removes outlier markers
#' with too many SNP number per locus/read.
#' The data requires snp and locus information (e.g. from a VCF file).
#' Having a higher than "normal" SNP number is usually the results of
#' assembly artifacts or bad assembly parameters.
#' This filter is population-agnostic, but still requires a strata
#' file if a vcf file is used as input.
#'
#' \strong{Filter targets}: Markers
#'
#' \strong{Statistics}: The number of SNPs per locus.


# Most arguments are inherited from tidy_genomic_data
#' @inheritParams read_strata
#' @inheritParams radiator_common_arguments

#' @param filter.snp.number (integer) This is best decided after viewing the figures.
#' If the argument is set to 2, locus with 3 and more SNPs will be blacklisted.
#' Default: \code{filter.snp.number = NULL}.

#' @param filename (optional) Name of the filtered tidy data frame file
#' written to the working directory (ending with \code{.tsv})
#' Default: \code{filename = NULL}.


#' @rdname filter_snp_number
#' @export
#' @details
#' \strong{Interactive version}
#'
#' There are 2 steps in the interactive version to visualize and filter
#' the data based on the number of SNP on the read/locus:
#'
#' Step 1. SNP number per read/locus visualization
#'
#' Step 2. Choose the filtering thresholds
#'
#'
#' @return A list in the global environment with 6 objects:
#' \enumerate{
#' \item $snp.number.markers
#' \item $number.snp.reads.plot
#' \item $whitelist.markers
#' \item $tidy.filtered.snp.number
#' \item $blacklist.markers
#' \item $filters.parameters
#' }
#'
#' The object can be isolated in separate object outside the list by
#' following the example below.

#' @examples
#' \dontrun{
#' turtle.outlier.snp.number <- radiator::filter_snp_number(
#' data = "turtle.vcf",
#' strata = "turtle.strata.tsv",
#' max.snp.number = 4,
#' filename = "tidy.data.turtle.tsv"
#' )

#'
#' tidy.data <- turtle.outlier.snp.number$tidy.filtered.snp.number
#'
#' #Inside the same list, to isolate the markers blacklisted:
#' blacklist <- turtle.outlier.snp.number$blacklist.markers
#'
#' }

filter_snp_number <- function(
  data,
  strata = NULL,
  interactive.filter = TRUE,
  filter.snp.number = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
  ...
) {

  # interactive.filter <- TRUE
  # data <- gds
  # path.folder <- "testing_snp_number"
  # force.stats <- TRUE
  # parameters <- NULL
  # filename <- NULL
  # filter.snp.number <- NULL
  # parallel.core <- parallel::detectCores() - 1
  # verbose = TRUE

  if (!is.null(filter.snp.number) || interactive.filter) {
    if (interactive.filter) verbose <- TRUE
    if (verbose) {
      cat("################################################################################\n")
      cat("############################ radiator::filter_snp_number #######################\n")
      cat("################################################################################\n")
    }
    # Cleanup-------------------------------------------------------------------
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
    on.exit(if (verbose) cat("######################### completed filter_snp_number ##########################\n"), add = TRUE)

    # Function call and dotslist -------------------------------------------------
    rad.dots <- radiator_dots(
      func.name = as.list(sys.call())[[1]],
      fd = rlang::fn_fmls_names(),
      args.list = as.list(environment()),
      dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE),
      keepers = c("path.folder", "parameters", "internal", "force.stats"),
      verbose = verbose
    )

    # Checking for missing and/or default arguments ------------------------------
    if (missing(data)) rlang::abort("data is missing")

    # Folders---------------------------------------------------------------------
    path.folder <- generate_folder(
      f = path.folder,
      rad.folder = "filter_snp_number",
      internal = internal,
      file.date = file.date,
      verbose = verbose)

    # write the dots file
    write_rad(
      data = rad.dots,
      path = path.folder,
      filename = stringi::stri_join("radiator_filter_snp_number_args_", file.date, ".tsv"),
      tsv = TRUE,
      internal = internal,
      verbose = verbose
    )

    # Message about steps taken during the process ---------------------------------
    if (interactive.filter) {
      message("Interactive mode: on")
      message("2 steps to visualize and filter the data based on the number of SNP on the read/locus:")
      message("Step 1. Impact of SNP number per read/locus (on individual genotypes and locus/snp number potentially filtered)")
      message("Step 2. Choose the filtering thresholds")
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
      internal = internal,
      verbose = verbose)

    # Whitelist and blacklist --------------------------------------------------
    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "COL")
    if (data.type == "SeqVarGDSClass") {
      wl <- extract_markers_metadata(gds = data, whitelist = TRUE) # not optimal, currently used just to get locus info
    } else {
      wl <- bl <- dplyr::select(data, dplyr::one_of(want))
    }
    # Check that required info is present in data: snp and locus
    if (!tibble::has_name(wl, "LOCUS") || !tibble::has_name(wl, "POS")) {
      problem.data <- "This filter requires dataset with SNP (POS) and LOCUS information (columns)"
      message("\n\n", problem.data)
      readr::write_lines(
        x = problem.data,
        path = file.path(path.folder, "README"))
      return(data)
    }

    # Generate snp per locus stats----------------------------------------------
    if (verbose) message("Generating statistics")
    if (data.type == "SeqVarGDSClass") {
      wl <- generate_markers_stats(
        gds = data,
        snp.per.locus = TRUE,
        missing = FALSE,
        coverage = FALSE,
        allele.coverage = FALSE,
        mac = FALSE,
        heterozygosity = FALSE,
        snp.position.read = FALSE,
        force.stats = force.stats,
        path.folder = path.folder,
        file.date = file.date,
        plot = FALSE,
        parallel.core = parallel.core
      )
      stats <- wl$stats
      wl <- bl <- wl$info
    } else {
      bl <- wl
      wl %<>%
        dplyr::group_by(LOCUS) %>%
        dplyr::mutate(SNP_PER_LOCUS = n()) %>%
        dplyr::ungroup(.)
      stats <- tibble_stats(
        x = wl$SNP_PER_LOCUS,
        group = "SNPs per locus")
    }

    if (tibble::has_name(wl, "COL")) {
      read.length <- max(wl$COL)
      if (verbose) message("\nWith max read length taken from data: ", read.length)

      if (verbose) message("    The max number of SNP per locus correspond to:")
      if (verbose) message("    1 SNP per ", round(read.length / stats[[6]]), " bp\n")
    } else {
      read.length <- NULL
    }


    # Generate box plot ---------------------------------------------------------
    snp.per.locus.fig <- boxplot_stats(
      data = stats,
      title =  "Number of SNPs per locus",
      subtitle = if (!is.null(read.length)) {
        stringi::stri_join("Read length (max): ", read.length, " bp", "\nOutlier: ", ceiling(stats[[9]]))
      } else {
        stringi::stri_join("\nOutlier: ", ceiling(stats[[9]]))
      },
      x.axis.title = NULL,
      y.axis.title = "Number of SNPs per locus",
      bp.filename = stringi::stri_join("snp_per_locus_", file.date, ".pdf"),
      path.folder = path.folder)


    # Distribution -------------------------------------------------------------
    d.plot <- wl %>%
      dplyr::distinct(LOCUS, SNP_PER_LOCUS) %>%
      ggplot2::ggplot(data = ., ggplot2::aes(factor(SNP_PER_LOCUS))) +
      ggplot2::geom_bar() +
      ggplot2::labs(x = "Number of SNPs per locus") +
      ggplot2::labs(y = "Distribution (number of locus)") +
      ggplot2::geom_vline(ggplot2::aes(xintercept = as.numeric(ceiling(stats[[9]]))), color = "yellow") +
      ggplot2::theme_bw()+
      ggplot2::theme(
        axis.title.x = ggplot2::element_text(size = 12, face = "bold"),
        axis.title.y = ggplot2::element_text(size = 12, face = "bold"),
        legend.title = ggplot2::element_text(size = 12, face = "bold"),
        legend.text = ggplot2::element_text(size = 12, face = "bold"),
        strip.text.x = ggplot2::element_text(size = 12, face = "bold"),
        axis.text.x = ggplot2::element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5)
      )
    print(d.plot)

    # save
    d.plot.filename <- stringi::stri_join("snp_per_locus_distribution_", file.date, ".pdf")

    ggplot2::ggsave(
      filename = file.path(path.folder, d.plot.filename),
      plot = d.plot,
      width = 20, height = 10, dpi = 300, units = "cm", useDingbats = FALSE,
      limitsize = FALSE)


    # Helper table -------------------------------------------------------------
    if (verbose) message("Generating helper table...")
    how_many_markers <- function(threshold, x) {
      nrow(dplyr::filter(x, SNP_PER_LOCUS > threshold))
    }#End how_many_markers

    snp.range <- stats[[2]]:stats[[6]]
    n.markers <- nrow(wl)

    helper.table <- tibble::tibble(SNP_PER_LOCUS = snp.range) %>%
      dplyr::mutate(
        BLACKLISTED_MARKERS = purrr::map_int(
          .x = snp.range, .f = how_many_markers, x = wl),
        WHITELISTED_MARKERS = n.markers - BLACKLISTED_MARKERS
      ) %>%
      readr::write_tsv(
        x = .,
        path = file.path(path.folder, "snp.per.locus.helper.table.tsv"))

    # figures
    markers.plot <- ggplot2::ggplot(
      data = helper.table  %<>% tidyr::gather(
        data = .,
        key = LIST, value = MARKERS, -SNP_PER_LOCUS),
      ggplot2::aes(x = SNP_PER_LOCUS, y = MARKERS)) +
      ggplot2::geom_line() +
      ggplot2::geom_point(size = 2, shape = 21, fill = "white") +
      ggplot2::geom_vline(ggplot2::aes(xintercept = as.numeric(ceiling(stats[[9]]))), color = "yellow") +
      ggplot2::scale_x_continuous(name = "Number of SNPs per locus allowed", breaks = snp.range) +
      ggplot2::scale_y_continuous(name = "Number of markers") +
      ggplot2::theme_bw()+
      ggplot2::theme(
        axis.title.x = ggplot2::element_text(size = 10, face = "bold"),
        axis.title.y = ggplot2::element_text(size = 10, face = "bold"),
        axis.text.x = ggplot2::element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5)
      ) +
      ggplot2::facet_grid(LIST ~. , scales = "free", space = "free")
    print(markers.plot)

    # save
    ggplot2::ggsave(
      filename = file.path(path.folder, "snp.per.locus.helper.plot.pdf"),
      plot = markers.plot,
      width = 20 + (stats[[6]] / 10),
      height = 15,
      dpi = 300,
      units = "cm",
      useDingbats = FALSE,
      limitsize = FALSE)
    helper.table <- markers.plot <- NULL
    if (verbose) message("Files written: helper tables and plots")

    # Step 2. Thresholds selection ---------------------------------------------
    if (interactive.filter) {
      max.allowed <- stats[[6]]
      if (verbose) message("\nStep 2. Filtering markers based on the maximum of SNPs per locus\n")

      filter.snp.number <- radiator_question(
        x = "Do you still want to blacklist markers? (y/n):",
        answer.opt = c("y", "n"))

      if (filter.snp.number == "y") {
        outlier.stats <- radiator_question(
          x = "Do you want to remove markers based on the outlier statistics or not (y/n) ?
          (n: next question will be to enter your own threshold)", answer.opt = c("y", "n"))
        if (outlier.stats == "y") {
          filter.snp.number <- "outliers"
        } else {
          filter.snp.number <- radiator_question(
            x = "Enter the maximum number of SNP per locus allowed:", minmax = c(1, max.allowed))
        }
        outlier.stats <- NULL
      } else {
        filter.snp.number <- NULL
      }
    }

    # Filtering ----------------------------------------------------------------
    if (!is.null(filter.snp.number)) {
      if (!purrr::is_double(filter.snp.number)) {
        out.high <- round(stats$OUTLIERS_HIGH[stats$GROUP == "SNPs per locus"])
        if (verbose) message("\nRemoving outliers markers based on the number of SNPs per locus statistic: ", out.high)
        filter.snp.number <- out.high
      } else {
        if (verbose) message("\nRemoving markers based on the number of SNPs per locus statistic: ", filter.snp.number)
      }
    } else {
      filter.snp.number <- 1000000000000
    }

    # Whitelist and Blacklist of markers
    if (!is.null(filter.snp.number)) {
      wl %<>% dplyr::filter(SNP_PER_LOCUS <= filter.snp.number)
    }
    readr::write_tsv(
      x = wl,
      path = file.path(path.folder, "whitelist.snp.per.locus.tsv"),
      append = FALSE, col_names = TRUE)
    bl %<>% dplyr::setdiff(wl) %>% dplyr::mutate(FILTERS = "filter.snp.number")

    readr::write_tsv(
      x = bl,
      path = file.path(path.folder, "blacklist.snp.per.locus.tsv"),
      append = FALSE, col_names = TRUE)
    # saving whitelist and blacklist
    if (verbose) message("File written: whitelist.markers.genotyping.tsv")
    if (verbose) message("File written: blacklist.markers.genotyping.tsv")

    if (data.type == "SeqVarGDSClass") {

      markers.meta <- extract_markers_metadata(gds = data) %>%
        dplyr::mutate(
          FILTERS = dplyr::if_else(
            VARIANT_ID %in% bl$VARIANT_ID, "filter.snp.number", FILTERS
          )
        )

      # Update GDS
      update_radiator_gds(
        gds = data,
        node.name = "markers.meta",
        value = markers.meta,
        sync = TRUE
      )



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
      filter.name = "Filter markers snp number",
      param.name = "filter.snp.number",
      values = filter.snp.number,
      path.folder = path.folder,
      file.date = file.date,
      internal = internal,
      verbose = verbose)

    # results ------------------------------------------------------------------
    radiator_results_message(
      rad.message = stringi::stri_join("\nFilter SNPs per locus threshold: ",
                                       filter.snp.number),
      filters.parameters,
      internal,
      verbose
    )
  }
  return(data)
} #End filter_snp_number
