# Discard monomorphic markers --------------------------------------------------

#' @name discard_monomorphic_markers

#' @title Discard monomorphic markers

#' @description Discard monomorphic markers.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.


#' @param verbose (optional, logical) \code{verbose = TRUE} to be chatty
#' during execution.
#' Default: \code{verbose = FALSE}.

#' @return A list with the filtered input file and the blacklist of markers removed.

#' @export
#' @rdname discard_monomorphic_markers
#' @importFrom dplyr select mutate group_by ungroup rename tally filter semi_join n_distinct
#' @importFrom stringi stri_replace_all_fixed stri_join
#' @importFrom tibble has_name
#' @keywords internal

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

discard_monomorphic_markers <- function(data, verbose = FALSE) {
  message("\n\nDeprecated function, update your code to use: filter_monomorphic function or filter.monomorphic argument\n\n")

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("Input file missing")

  # Detect format --------------------------------------------------------------
  data.type <- radiator::detect_genomic_format(data)

  if (data.type == "SeqVarGDSClass") {
    if (verbose) message("Scanning for monomorphic markers...")
    # radiator::summary_gds(gds = data)
    # radiator::sync_gds(gds = data, samples = NULL, markers = NULL)
    markers <- SeqArray::seqGetData(data, "variant.id")
    n.markers.before <- length(markers)
    mono.markers <- markers[SeqArray::seqNumAllele(gdsfile = data) == 1]
    # mono.markers <- c(6L,7L)
    n.markers.removed <- length(mono.markers)

    if (n.markers.removed > 0) {
      whitelist.polymorphic.markers <- markers[-mono.markers]
      n.markers.after <- n.markers.before - n.markers.removed
      radiator::sync_gds(gds = data, markers = whitelist.polymorphic.markers)
    } else {
      whitelist.polymorphic.markers <- markers
      n.markers.after <- n.markers.before
    }

    # check for markers meta node
    # markers.meta <- radiator_gds_markers_meta(gds = data)
    #
    # if (is.null(markers.meta)) {
    #   suppressWarnings(gdsfmt::add.gdsn(
    #     node = radiator.gds,
    #     name = "markers.meta",
    #     val = res$markers.meta,
    #     replace = TRUE,
    #     compress = "ZIP_RA",
    #     closezip = TRUE))
    # }
    res <- list(input = data,
                blacklist.monomorphic.markers = mono.markers,
                whitelist.polymorphic.markers = whitelist.polymorphic.markers
    )
    if (verbose) {
      n.markers <- stringi::stri_join(n.markers.before, n.markers.removed, n.markers.after, sep = "/")
      message("    Number of markers before/blacklisted/after: ", n.markers)
    }



  } else {
    # Import data ---------------------------------------------------------------
    data <- radiator::tidy_wide(data = data, import.metadata = TRUE)


    if (tibble::has_name(data, "CHROM")) {
      markers.df <- dplyr::distinct(.data = data, MARKERS, CHROM, LOCUS, POS)
    }
    if (verbose) message("Scanning for monomorphic markers...")
    n.markers.before <- dplyr::n_distinct(data$MARKERS)

    if (tibble::has_name(data, "GT_BIN")) {
      mono.markers <- dplyr::select(.data = data, MARKERS, GT_BIN) %>%
        dplyr::filter(!is.na(GT_BIN)) %>%
        dplyr::distinct(MARKERS, GT_BIN) %>%
        dplyr::count(x = ., MARKERS) %>%
        dplyr::filter(n == 1) %>%
        dplyr::distinct(MARKERS)
    } else {
      mono.markers <- dplyr::select(.data = data, MARKERS, GT) %>%
        dplyr::filter(GT != "000000") %>%
        dplyr::distinct(MARKERS, GT) %>%
        dplyr::mutate(
          A1 = stringi::stri_sub(GT, 1, 3),
          A2 = stringi::stri_sub(GT, 4,6)
        ) %>%
        dplyr::select(-GT) %>%
        tidyr::gather(data = ., key = ALLELES_GROUP, value = ALLELES, -MARKERS) %>%
        dplyr::distinct(MARKERS, ALLELES) %>%
        dplyr::count(x = ., MARKERS) %>%
        dplyr::filter(n == 1) %>%
        dplyr::distinct(MARKERS)
    }
    # Remove the markers from the dataset
    n.markers.removed <- nrow(mono.markers)

    if (length(mono.markers$MARKERS) > 0) {
      data <- dplyr::anti_join(data, mono.markers, by = "MARKERS")
      if (tibble::has_name(data, "CHROM")) {
        mono.markers <- dplyr::left_join(mono.markers, markers.df, by = "MARKERS")
      }
    } else {
      mono.markers <- tibble::data_frame(MARKERS = character(0))
    }
    n.markers.after <- dplyr::n_distinct(data$MARKERS)

    want <- c("MARKERS", "CHROM", "LOCUS", "POS")
    whitelist.polymorphic.markers <- suppressWarnings(
      dplyr::select(data, dplyr::one_of(want)) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE))
    res <- list(input = data,
                blacklist.monomorphic.markers = mono.markers,
                whitelist.polymorphic.markers = whitelist.polymorphic.markers
    )
    if (verbose) {
      n.markers <- stringi::stri_join(n.markers.before, n.markers.removed, n.markers.after, sep = "/")
      message("    Number of markers before/blacklisted/after: ", n.markers)
    }
  }

  return(res)
} # end discard mono markers

# Keep markers in common between all populations

#' @name keep_common_markers

#' @title Keep markers in common between all populations

#' @description Keep markers in common between all populations.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.


#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.

#' @param plot (optional, logical) \code{plot = TRUE} will produce an
#' \href{https://github.com/hms-dbmi/UpSetR}{UpSet plot} to visualize the number
#' of markers between populations. The package is required for this to work...
#' Default: \code{plot = FALSE}.

#' @param verbose (optional, logical) \code{verbose = TRUE} to be chatty
#' during execution.
#' Default: \code{verbose = FALSE}.

#' @return A list with the filtered input data set,
#' the whitelist of markers in common and the blacklist of markers discarded.


#' @export
#' @keywords internal

#' @rdname keep_common_markers
#' @importFrom dplyr select mutate group_by ungroup rename tally filter semi_join n_distinct
#' @importFrom stringi stri_replace_all_fixed stri_join
#' @importFrom tibble has_name
# @importFrom UpSetR upset

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @section Life cycle:
#'
#' As of radiator v.0.0.22, keep.common.markers function is deprecated and was replaced by
#' \code{\link{filter_common_markers}} in an effort to have more meaningful functions names.


# keep_common_markers ----------------------------------------------------------
keep_common_markers <- function(data, plot = FALSE, verbose = FALSE) {

  message("\n\nDeprecated function, update your code to use: filter_common_markers\n\n")


  if (plot) {
    if (!requireNamespace("UpSetR", quietly = TRUE)) {
      rlang::abort("UpSetR needed for this function to work
                   Install with install.packages('UpSetR')")
    }
    }

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("Input file missing")

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    data <- radiator::tidy_wide(data = data, import.metadata = TRUE)
  }

  # check genotype column naming
  if (tibble::has_name(data, "GENOTYPE")) {
    colnames(data) <- stringi::stri_replace_all_fixed(
      str = colnames(data),
      pattern = "GENOTYPE",
      replacement = "GT",
      vectorize_all = FALSE)
  }

  # necessary steps to make sure we work with unique markers and not duplicated LOCUS
  if (tibble::has_name(data, "LOCUS") && !tibble::has_name(data, "MARKERS")) {
    data <- dplyr::rename(.data = data, MARKERS = LOCUS)
  }

  # markers.meta
  want <- c("MARKERS", "CHROM", "LOCUS", "POS")
  markers.meta <- suppressWarnings(dplyr::select(data, dplyr::one_of(want)) %>%
                                     dplyr::distinct(MARKERS, .keep_all = TRUE))

  if (verbose) message("Using markers common in all populations:")

  if (tibble::has_name(data, "GT_BIN")) {
    blacklist <- dplyr::select(.data = data, MARKERS, POP_ID, GT_BIN) %>%
      dplyr::filter(!is.na(GT_BIN))
  } else {
    blacklist <- dplyr::select(.data = data, MARKERS, POP_ID, GT) %>%
      dplyr::filter(GT != "000000")
  }

  blacklist <- dplyr::distinct(blacklist, MARKERS, POP_ID) %>%
    dplyr::count(x = ., MARKERS) %>%
    dplyr::filter(n != dplyr::n_distinct(data$POP_ID)) %>%
    dplyr::distinct(MARKERS) %>%
    dplyr::arrange(MARKERS)

  markers.data <- dplyr::n_distinct(data$MARKERS)
  blacklist.markers <- nrow(blacklist)
  markers.in.common <- markers.data - blacklist.markers

  if (verbose) {
    n.markers <- stringi::stri_join(markers.data, blacklist.markers, markers.in.common, sep = "/")
    message("    Number of markers before/blacklisted/after:", n.markers)
  }

  if (plot) {
    pops <- unique(data$POP_ID)

    if (length(pops) > 1) {
      if (tibble::has_name(data, "GT_BIN")) {
        plot.data <- dplyr::filter(data, !is.na(GT_BIN))
      } else {
        plot.data <- dplyr::filter(data, GT != "000000")
      }
      plot.data <- dplyr::distinct(plot.data, MARKERS, POP_ID) %>%
        dplyr::mutate(
          n = rep(1, n()),
          POP_ID = stringi::stri_join("POP_", POP_ID)
        ) %>%
        tidyr::spread(data = ., key = POP_ID, value = n, fill = 0) %>%
        data.frame(.)

      UpSetR::upset(plot.data, nsets = length(pops), order.by = "freq", empty.intersections = "on")
      message("    UpSet plot generated to visualize common markers")
    }
  }

  if (blacklist.markers > 0) {
    # system.time(data2 <- dplyr::anti_join(data, blacklist, by = "MARKERS"))
    data <- dplyr::filter(data, !MARKERS %in% blacklist$MARKERS)

    if (ncol(markers.meta) > 1) {
      blacklist <- dplyr::left_join(blacklist, markers.meta, by = "MARKERS") %>%
        readr::write_tsv(x = ., path = "blacklist.not.in.common.markers.tsv")
    }

  } else {
    blacklist <- tibble::data_frame(INDIVIDUALS = character(0))
  }
  want <- c("MARKERS", "CHROM", "LOCUS", "POS")
  whitelist <- suppressWarnings(dplyr::select(data, dplyr::one_of(want)) %>%
                                  dplyr::distinct(MARKERS, .keep_all = TRUE))


  if (blacklist.markers > 0) {
    readr::write_tsv(x = whitelist, path = "whitelist.common.markers.tsv")
  }

  return(res = list(input = data,
                    whitelist.common.markers = whitelist,
                    blacklist.not.in.common.markers = blacklist))
  }

# Filter dart function ---------------------------------------------------------
#' @name filter_dart

#' @title filter DArT file

#' @description Use \code{\link{filter_rad}}.
#' @param ... Read life cycle section below
#' @section Life cycle:
#'
#' As of radiator v.0.0.22, \code{filter_dart} is now replaced by
#' \code{\link{filter_rad}} in an effort to have more meaningful functions names.
#' This is THE function to rule them all.

#' @export
#' @rdname filter_dart

filter_dart <- function(...) {
  rlang::abort("As of radiator v.0.0.22, filter_dart is now replaced by filter_rad.
               \nRead documentation, for latest changes, and modify your codes!\n")
}


# filter_individual-------------------------------------------------------------
#' @name filter_individual
#' @title Individual filter
#' @description Filter individuals genotyped at a marker from a tidy data
#' set (long format) of any of these file format:
#' vcf, plink (tped/tfam), stacks haplotype file, genind,
#' genepop, data frame in wide format. The function uses
#' \code{\link[radiator]{tidy_genomic_data}} and
#' \code{\link[radiator]{tidy_wide}} to load the file. For filtering
#' you can consider the overall number of individuals (no concept of populations here),
#' or the number of genotyped individuals per pop. The threshold can be a fixed
#' number of individuals, a proportion or a percentage.

# Most arguments are inherited from tidy_genomic_data
#' @inheritParams tidy_genomic_data
#' @inheritParams read_strata
#' @inheritParams radiator_common_arguments

#' @param interactive.filter (optional, logical) Do you want the filtering session to
#' be interactive. With default: \code{interactive.filter == TRUE}, the user is
#' asked to see figures of distribution before making decisions for filtering.

#' @param ind.approach (character).
#' The approach to filter a marker, is it based on the overall number of
#' genotyped individuals (no concept of populations here),
#' \code{ind.approach = "overall"}, or
#' the number of genotyped individuals per population \code{ind.approach = "pop"}.
#' Default: \code{ind.approach = "pop"}.

#' @param ind.threshold The individual threshold, proportion, percentage or
#' number e.g. 0.70, 70, 15.
#' Default: \code{ind.threshold = 0.5}.

#' @param percent Is the threshold a percentage? TRUE or FALSE.
#' This argument is necessary to distinguish percentage from integer individual
#' threshold (e.g. 70 percent or 70 individuals).
#' Default: \code{percent = FALSE}.

#' @param prob.pop.threshold (integer, optional) Useful to incorporate problematic
#' populations dragging down polymorphism discovery, but still wanted for analysis.
#' If individuals with missing genotypes are not managed upstream,
#' use this threshold to allow variance in the number of populations passing
#' the \code{ind.threshold} argument. e.g. with \code{prob.pop.threshold = 2},
#' you tolerate a maximum of 2 populations failing to pass the
#' \code{ind.threshold}. Manage after the problematic populations
#' downstream with blacklist of individuals and/or missing data imputations. This
#' argument is not necessary with \code{ind.approach = "overall"}.
#' Default: \code{prob.pop.threshold = 0}. See details for more info.

#' @param filename (optional) The function uses \code{\link[fst]{write.fst}},
#' to write the tidy data frame in
#' the folder created in the working directory. The file extension appended to
#' the \code{filename} provided is \code{.rad}.
#' Default: \code{filename = NULL}.


#' @rdname filter_individual
#' @export
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot stat_summary labs theme element_blank element_text geom_jitter scale_colour_manual scale_y_reverse theme_light geom_bar facet_grid geom_histogram geom_bar aes_string scale_fill_manual theme_bw stat_smooth geom_boxplot ggsave
#' @importFrom stringi stri_join stri_replace_all_fixed stri_sub stri_detect_fixed
#' @importFrom dplyr select distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join funs
#' @importFrom readr write_tsv
#' @importFrom tidyr spread gather
#' @importFrom tibble data_frame has_name

#' @details
#' \strong{Interactive version}
#'
#' There is 2 steps in the interactive version to visualize and filter
#' the data based on the number of genotyped individuals:
#'
#' Step 1. Impact of individual threshold on marker discovery
#'
#' Step 2. Choose the filtering approach and thresholds
#'
#'
#'
#' \strong{prob.pop.threshold}
#'
#' If your a regular radiator user, you've seen the \code{pop.num.threshold}.
#' \code{prob.pop.threshold}, is a bit different and requires more thinking,
#' because the number of populations genotyped potentially vary across markers,
#' which make difficult, sometimes, the use of \code{pop.num.threshold}.
#' e.g. If only 2 populations out of 10 are genotyped for a marker and you use
#' \code{pop.num.threshold = 0.6} this could lead to inconsistensis. Thus,
#' it's easier to use the number of population you are willing to accept
#' failling to conform to the threshold.


#' @return With \code{interactive.filter = FALSE}, a list in the global environment,
#' with 7 objects:
#' \enumerate{
#' \item $tidy.filtered.ind
#' \item $whitelist.markers
#' \item $blacklist.markers
#' \item $strata
#' \item $filters.parameters
#' }
#'
#' With \code{interactive.filter = TRUE}, a list with 2 additionnal objects is created.
#' \enumerate{
#' \item $plot.ind.threshold
#' \item $ind.threshold.helper.table
#' }
#'
#' The object can be isolated in separate object outside the list by
#' following the example below.

#' @examples
#' \dontrun{
#' turtle.ind <- filter_individual(
#' data = "turtle.vcf",
#' strata = "turtle.strata.tsv",
#' ind.approach = "pop",
#' ind.threshold = 50,
#' percent = TRUE,
#' prob.pop.threshold = 2,
#' filename = "tidy.data.turtle"
#' )

#'
#' #If interactive.filter = TRUE, a list is created and to view the filtered tidy data:
#' tidy.data <- turtle.ind$tidy.filtered.ind
#'
#' #Inside the same list, to isolate the blacklist.genotypes:
#' bg <- turtle.ind$blacklist.genotypes
#'
#' # The remaining argument are used in tidy_genomic_data during import and allow
#' # the user to apply some filtering or selection before doing the filtering.
#' }

filter_individual <- function(
  interactive.filter = TRUE,
  data,
  strata = NULL,
  ind.approach = "pop",
  ind.threshold = 0.5,
  prob.pop.threshold = 0,
  percent = FALSE,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  ...
) {
  cat("#######################################################################\n")
  cat("##################### radiator::filter_individual #######################\n")
  cat("#######################################################################\n")
  opt.change <- getOption("width")
  options(width = 70)
  timing <- proc.time()
  # manage missing arguments -----------------------------------------------------
  if (missing(data)) rlang::abort("Input file missing")
  if (!is.null(pop.levels) & is.null(pop.labels)) {
    pop.levels <- stringi::stri_replace_all_fixed(pop.levels, pattern = " ", replacement = "_", vectorize_all = FALSE)
    pop.labels <- pop.levels
  }

  if (!is.null(pop.labels) & is.null(pop.levels)) rlang::abort("pop.levels is required if you use pop.labels")

  if (!is.null(pop.labels)) {
    if (length(pop.labels) != length(pop.levels)) rlang::abort("pop.labels and pop.levels must have the same length (number of groups)")
    pop.labels <- stringi::stri_replace_all_fixed(pop.labels, pattern = " ", replacement = "_", vectorize_all = FALSE)
  }

  if (!is.null(pop.select)) {
    pop.select <- stringi::stri_replace_all_fixed(pop.select, pattern = " ", replacement = "_", vectorize_all = FALSE)
  }
  # Message about steps taken during the process ---------------------------------
  if (interactive.filter) {
    message("Interactive mode: on")
    message("2 steps to visualize and filter the data based on the number of genotyped individuals:")
    message("Step 1. Impact of individual threshold on marker discovery")
    message("Step 2. Choose the filtering approach and thresholds")
  }

  # Folder -------------------------------------------------------------------
  # Date and time
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  folder.message <- stringi::stri_join("filter_individual_", file.date, sep = "")
  path.folder <- stringi::stri_join(getwd(),"/", folder.message, sep = "")
  dir.create(file.path(path.folder))

  message("\nFolder created: \n", folder.message)
  file.date <- NULL #unused object

  # Filter parameter file ------------------------------------------------------
  message("Parameters used in this run are stored in a file")
  filters.parameters <- list.files(path = getwd(), pattern = "filters_parameters.tsv", full.names = TRUE)
  if (length(filters.parameters) == 0) {
    filters.parameters <- tibble::data_frame(FILTERS = as.character(), PARAMETERS = as.character(), VALUES = as.integer(), BEFORE = as.character(), AFTER = as.character(), BLACKLIST = as.integer(), UNITS = as.character(), COMMENTS = as.character())
    readr::write_tsv(x = filters.parameters, path = "filters_parameters.tsv", append = FALSE, col_names = TRUE)
    message("    Created a parameter file: filters_parameters.tsv")
  } else {
    message("    Using the filters parameters file: filters_parameters.tsv")
  }

  # File type detection----------------------------------------------------------
  data.type <- radiator::detect_genomic_format(data)

  # import data ----------------------------------------------------------------
  message("Importing data ...")
  input <- radiator::tidy_genomic_data(
    data = data,
    vcf.metadata = vcf.metadata,
    blacklist.id = blacklist.id,
    blacklist.genotype = blacklist.genotype,
    whitelist.markers = whitelist.markers,
    monomorphic.out = monomorphic.out,
    max.marker = max.marker,
    snp.ld = snp.ld,
    common.markers = common.markers,
    strata = strata,
    pop.levels = pop.levels,
    pop.labels = pop.labels,
    pop.select = pop.select,
    filename = NULL,
    parallel.core = parallel.core
  )

  # create a strata.df
  strata.df <- input %>%
    dplyr::select(INDIVIDUALS, POP_ID) %>%
    dplyr::distinct(INDIVIDUALS, .keep_all = TRUE)
  strata <- strata.df
  pop.levels <- levels(input$POP_ID)
  pop.labels <- pop.levels

  # prepare filter, table and figure----------------------------------------------
  ind.total <- dplyr::n_distinct(input$INDIVIDUALS) # total number of individuals
  pop.number <- dplyr::n_distinct(input$POP_ID) # number of pop

  # individuals per pop
  ind.pop <- input %>%
    dplyr::distinct(POP_ID, INDIVIDUALS) %>%
    dplyr::group_by(POP_ID) %>%
    dplyr::tally(.) %>%
    dplyr::rename(N_IND = n)

  # input genotyped
  input.genotyped <- input %>%
    dplyr::filter(GT != "000000")

  # overall genotyped individuals
  overall <- input.genotyped %>%
    dplyr::select(MARKERS, INDIVIDUALS) %>%
    dplyr::group_by(MARKERS) %>%
    dplyr::tally(.) %>%
    dplyr::rename(GENOTYPED = n) %>%
    dplyr::mutate(PERCENT = ceiling(GENOTYPED/ind.total*100))

  # number of pop. genotyped per marker
  pop.genotyped.marker <- input.genotyped %>%
    dplyr::distinct(MARKERS, POP_ID) %>%
    dplyr::group_by(MARKERS) %>%
    dplyr::tally(.) %>%
    dplyr::rename(POP_GENOTYPED = n)

  # genotyped individuals per pop
  pop <- input.genotyped %>%
    dplyr::select(MARKERS, INDIVIDUALS, POP_ID) %>%
    dplyr::group_by(MARKERS, POP_ID) %>%
    dplyr::tally(.) %>%
    dplyr::rename(GENOTYPED = n) %>%
    dplyr::inner_join(ind.pop, by = "POP_ID") %>%
    dplyr::mutate(PERCENT = ceiling(GENOTYPED/N_IND*100))

  input.genotyped <- NULL # unused object

  # Step 1. Impact of individual threshold on marker discovery------------------
  if (interactive.filter) {
    message("Step 1. Impact of individual threshold on marker discovery")
  }
  # line.graph <- as.character(readLines(n = 1))
  # if (line.graph == "y") {
  threshold.helper.overall <- overall %>%
    dplyr::ungroup(.) %>%
    dplyr::summarise(
      `10` = length(PERCENT[PERCENT >= 10]),
      `20` = length(PERCENT[PERCENT >= 20]),
      `30` = length(PERCENT[PERCENT >= 30]),
      `40` = length(PERCENT[PERCENT >= 40]),
      `50` = length(PERCENT[PERCENT >= 50]),
      `60` = length(PERCENT[PERCENT >= 60]),
      `70` = length(PERCENT[PERCENT >= 70]),
      `80` = length(PERCENT[PERCENT >= 80]),
      `90` = length(PERCENT[PERCENT >= 90]),
      `100` = length(PERCENT[PERCENT == 100])
    ) %>%
    tidyr::gather(data = ., key = IND_THRESHOLD, value = MARKER_NUMBER) %>%
    dplyr::mutate(POP_ID = rep("OVERALL", n()))

  threshold.helper.pop <- pop %>%
    dplyr::group_by(POP_ID) %>%
    dplyr::summarise(
      `10` = length(PERCENT[PERCENT >= 10]),
      `20` = length(PERCENT[PERCENT >= 20]),
      `30` = length(PERCENT[PERCENT >= 30]),
      `40` = length(PERCENT[PERCENT >= 40]),
      `50` = length(PERCENT[PERCENT >= 50]),
      `60` = length(PERCENT[PERCENT >= 60]),
      `70` = length(PERCENT[PERCENT >= 70]),
      `80` = length(PERCENT[PERCENT >= 80]),
      `90` = length(PERCENT[PERCENT >= 90]),
      `100` = length(PERCENT[PERCENT == 100])
    ) %>%
    tidyr::gather(data = ., key = IND_THRESHOLD, value = MARKER_NUMBER, -POP_ID)

  mean.pop <- threshold.helper.pop %>%
    dplyr::group_by(IND_THRESHOLD) %>%
    dplyr::summarise(
      MARKER_NUMBER = round(mean(MARKER_NUMBER), 0)
    ) %>%
    dplyr::mutate(POP_ID = rep("MEAN_POP", n()))

  threshold.helper <- suppressWarnings(
    dplyr::bind_rows(threshold.helper.pop, mean.pop, threshold.helper.overall) %>%
      dplyr::mutate(
        IND_THRESHOLD = as.numeric(IND_THRESHOLD),
        POP_ID = factor(POP_ID, levels = c(levels(input$POP_ID), "MEAN_POP", "OVERALL"), ordered = TRUE)
      ))

  # Set the breaks for the figure
  max.markers <- dplyr::n_distinct(input$MARKERS)

  #Function to replace packageplyr round_any
  rounder <- function(x, accuracy, f = round) {
    f(x / accuracy) * accuracy
  }

  # max.markers <- 658
  if (max.markers >= 1000) {
    y.breaks.by <- rounder(max.markers/10, 100, ceiling)
    y.breaks.max <- rounder(max.markers, 1000, ceiling)
    y.breaks <- seq(0, y.breaks.max, by = y.breaks.by)
  } else {
    y.breaks.by <- rounder(max.markers/10, 10, ceiling)
    y.breaks.max <- rounder(max.markers, 100, ceiling)
    y.breaks <- seq(0, y.breaks.max, by = y.breaks.by)
  }

  plot.ind.threshold <- ggplot2::ggplot(threshold.helper, ggplot2::aes(x = IND_THRESHOLD, y = MARKER_NUMBER)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2, shape = 21, fill = "white") +
    ggplot2::scale_x_continuous(name = "Individual threshold (percent)", breaks = seq(10, 100, by = 10)) +
    ggplot2::scale_y_continuous(name = "Number of markers", breaks = y.breaks, limits = c(0, y.breaks.max)) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5),
      strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
    ) +
    ggplot2::facet_grid(~POP_ID)
  if (interactive.filter) {
    print(plot.ind.threshold)
    message("    Look at the plot and inspect the change in the number of markers
            in relation to the individual percentage thresholds\n")
  }

  # save
  ggplot2::ggsave(stringi::stri_join(path.folder, "/plot.ind.threshold.pdf"), width = pop.number*4, height = 10, dpi = 600, units = "cm", useDingbats = FALSE)
  ggplot2::ggsave(stringi::stri_join(path.folder, "/plot.ind.threshold.png"), width = pop.number*4, height = 10, dpi = 300, units = "cm")
  message("2 versions (pdf and png) of the line graph of individual threshold on marker discovery were written in the folder")

  # Helper table for individual thresholds
  ind.threshold.helper.table <- threshold.helper %>%
    dplyr::group_by(POP_ID) %>%
    tidyr::spread(data = ., key = IND_THRESHOLD, MARKER_NUMBER)

  if (interactive.filter) {
    message("\n    Inspect the table (ind.threshold.helper.table) to help you view
            the relation between individual thresholds and marker discovery:")
    print(ind.threshold.helper.table)
    message("    First column: POP_ID")
    message("    Remaining columns: the individual thresholds in percent
            with the value = the number of markers discovered")
    message("    The last 2 rows: MEAN_POP is the mean across your populations
            and OVERALL is if you had 1 large population")
  }
  readr::write_tsv(
    x = ind.threshold.helper.table,
    path = stringi::stri_join(path.folder, "/", "ind.threshold.helper.table.tsv"))
  message("    ind.threshold.helper.table was written in the folder")

  # Step 2. Choose the filtering approach and thresholds------------------------
  # 2 approach: filtering with the overall n. ind. ("overall") or by pop ("pop")
  if (interactive.filter) {
    message("Step 2. Choose the filtering approach and thresholds")
    message("    The approach to filter a marker: do you want it based on the overall
            number of genotyped individuals or
            on the number of genotyped individuals per pop ? (overall or pop):")
    ind.approach <- as.character(readLines(n = 1))

    message("    Enter the individual threshold (number, proportion or percentage).
            e.g. enter 10 (for 10 individuals), 0.8 for proportion and 80 for percentage.")
    ind.threshold <- as.numeric(readLines(n = 1))

    message("    The value you just entered for the individual threshold, is it a percentage? (TRUE/FALSE)")
    percent <- as.character(readLines(n = 1))
    if (ind.approach == "pop") {
      message("Tolerance for deviation: look at the plot produced ealier and if you see some populations dragging down
              the number of markers for certain percentage thresholds, you have 3 options:\n
              1. remove the population (use pop.select argument to keep the desired populations)
              2. remove individuals with too many missing genotypes (use blacklist.id argument)
              3. use the next threshold (below) to allow variance and then
              manage the missing values with blacklist of individuals and/or
              missing data imputations.\n
              Enter the number of problematic population that you allow to deviate from the threshold:")
      prob.pop.threshold <- as.numeric(readLines(n = 1))
    }
  }

  # Filtering ------------------------------------------------------------------
  message("    Filtering...")

  # some discrepencies need to be highllighted here. If you have entire pop not genotyped for a markers
  # this will compute them when doing the filtering:
  # summarise(GENOTYPED = length(INDIVIDUALS[GT != "000000"]))
  # if we first remove the individuals not genotyped with :
  # filter(GT != "000000")
  # the pop not genotyped are not accounted for. And this is what we want here.
  # filter_populations take care the ungenotyped pop and common.markers make sure that for certain analysis
  # you can have common markers or not.
  # so here we focus on when pop got a marker is it at 50% 60% 70%  ... genotyped?
  if (ind.approach == "overall") {
    if (stringi::stri_detect_fixed(percent, "F") & ind.threshold > 1) {
      # ind.threshold <- 100
      threshold.id <- "(Using a fixed threshold)"
      prob.pop.threshold <- "NA"
      ind.approach <- "overall individuals (no pop)"
      filter <- overall %>%
        dplyr::filter(GENOTYPED >= ind.threshold) %>%
        dplyr::distinct(MARKERS)
    } else if (stringi::stri_detect_fixed(percent, "F") & stringi::stri_detect_fixed(ind.threshold, ".") & ind.threshold < 1) {
      # ind.threshold <- 0.7
      threshold.id <- "(proportion)"
      prob.pop.threshold <- "NA"
      ind.approach <- "overall individuals (no pop)"
      filter <- overall %>%
        dplyr::filter(PERCENT >= ind.threshold*100) %>%
        dplyr::distinct(MARKERS)
    } else {# percent
      # ind.threshold <- 70
      threshold.id <- "(percent)"
      prob.pop.threshold <- "NA"
      ind.approach <- "overall individuals (no pop)"
      filter <- overall %>%
        dplyr::filter(PERCENT >= ind.threshold) %>%
        dplyr::distinct(MARKERS)
    }
  } else {# approach by pop
    if (stringi::stri_detect_fixed(percent, "F") & ind.threshold > 1) {
      message("Using a fixed threshold")
      threshold.id <- "(ind.)"
      ind.approach <- "individuals by pop"
      # ind.threshold <- 15
      # prob.pop.threshold <- 3
      filter <- pop %>%
        dplyr::ungroup(.) %>%
        dplyr::filter(GENOTYPED >= ind.threshold) %>%
        dplyr::group_by(MARKERS) %>%
        dplyr::tally(.) %>%
        dplyr::inner_join(pop.genotyped.marker, by = "MARKERS") %>%
        dplyr::mutate(PROB_POP = POP_GENOTYPED - n) %>%
        dplyr::filter(PROB_POP <= prob.pop.threshold) %>%
        dplyr::distinct(MARKERS)
    } else if (stringi::stri_detect_fixed(percent, "F") & stringi::stri_detect_fixed(ind.threshold, ".") & ind.threshold < 1) {
      message("Using a proportion threshold...")
      threshold.id <- "(proportion)"
      ind.approach <- "individuals by pop"
      # ind.threshold <- 0.6
      # prob.pop.threshold <- 3
      filter <- pop %>%
        dplyr::ungroup(.) %>%
        dplyr::filter(PERCENT >= ind.threshold*100) %>%
        dplyr::group_by(MARKERS) %>%
        dplyr::tally(.) %>%
        dplyr::inner_join(pop.genotyped.marker, by = "MARKERS") %>%
        dplyr::mutate(PROB_POP = POP_GENOTYPED - n) %>%
        dplyr::filter(PROB_POP <= prob.pop.threshold) %>%
        dplyr::distinct(MARKERS)
    } else {
      message("Using a percentage threshold...")
      threshold.id <- "(percent)"
      ind.approach <- "individuals by pop"
      # ind.threshold <- 0.6
      # prob.pop.threshold <- 3
      filter <- pop %>%
        dplyr::ungroup(.) %>%
        dplyr::filter(PERCENT >= ind.threshold) %>%
        dplyr::group_by(MARKERS) %>%
        dplyr::tally(.) %>%
        dplyr::inner_join(pop.genotyped.marker, by = "MARKERS") %>%
        dplyr::mutate(PROB_POP = POP_GENOTYPED - n) %>%
        dplyr::filter(PROB_POP <= prob.pop.threshold) %>%
        dplyr::distinct(MARKERS)
    }
  }

  # Apply the filter to the tidy data
  filter <- dplyr::left_join(x = filter, input, by = "MARKERS")

  # Update filters.parameters SNP ----------------------------------------------

  if (data.type == "haplo.file") {
    snp.before <- as.character("NA")
    snp.after <- as.character("NA")
    snp.blacklist <- as.character("NA")
  } else {
    snp.before <- as.integer(dplyr::n_distinct(input$MARKERS))
    snp.after <- as.integer(dplyr::n_distinct(filter$MARKERS))
    snp.blacklist <- as.integer(dplyr::n_distinct(input$MARKERS) - dplyr::n_distinct(filter$MARKERS))
  }

  if (tibble::has_name(input, "LOCUS")) {
    locus.before <- as.integer(dplyr::n_distinct(input$LOCUS))
    locus.after <- as.integer(dplyr::n_distinct(filter$LOCUS))
    locus.blacklist <- as.integer(dplyr::n_distinct(input$LOCUS) - dplyr::n_distinct(filter$LOCUS))
  } else {
    locus.before <- as.character("NA")
    locus.after <- as.character("NA")
    locus.blacklist <- as.character("NA")
  }


  markers.before <- stringi::stri_join(snp.before, locus.before, sep = "/")
  markers.after <- stringi::stri_join(snp.after, locus.after, sep = "/")
  markers.blacklist <- stringi::stri_join(snp.blacklist, locus.blacklist, sep = "/")


  filters.parameters <- tibble::data_frame(
    FILTERS = c("Individuals", rep(as.character(""), 2)),
    PARAMETERS = c("ind.approach", "ind.threshold", "prob.pop.threshold"),
    VALUES = c(ind.approach, paste(">=", ind.threshold, " ", threshold.id), prob.pop.threshold),
    BEFORE = c("", "", markers.before),
    AFTER = c("", "", markers.after),
    BLACKLIST = c("", "", markers.blacklist),
    UNITS = c("", "genotyped", "SNP/LOCUS"),
    COMMENTS = c("", "", "")
  )
  readr::write_tsv(x = filters.parameters, path = "filters_parameters.tsv", append = TRUE, col_names = FALSE)

  # saving tidy data
  if (!is.null(filename)) {
    tidy.name <- stringi::stri_join(filename, ".rad")
    message("Writing the filtered tidy data set: ", tidy.name)
    write_rad(data = filter, path = file.path(path.folder, tidy.name))
  }

  # saving whitelist
  message("Writing the whitelist of markers: whitelist.markers.ind.tsv")

  if (tibble::has_name(input, "CHROM")) {
    whitelist.markers <- dplyr::ungroup(filter) %>%
      dplyr::distinct(CHROM, LOCUS, POS)
  } else {
    whitelist.markers <- dplyr::ungroup(filter) %>%
      dplyr::distinct(MARKERS)
  }
  readr::write_tsv(whitelist.markers, stringi::stri_join(path.folder, "/", "whitelist.markers.ind.tsv"), append = FALSE, col_names = TRUE)


  # saving blacklist
  message("Writing the blacklist of markers: blacklist.markers.ind.tsv")
  if (tibble::has_name(input, "CHROM")) {
    blacklist.markers <- dplyr::ungroup(input) %>%
      dplyr::distinct(CHROM, LOCUS, POS) %>%
      dplyr::anti_join(whitelist.markers, by = c("CHROM", "LOCUS", "POS"))
  } else {
    blacklist.markers <- dplyr::ungroup(input) %>%
      dplyr::distinct(MARKERS) %>%
      dplyr::anti_join(whitelist.markers, by = "MARKERS")
  }
  readr::write_tsv(blacklist.markers, stringi::stri_join(path.folder, "/", "blacklist.markers.ind.tsv"), append = FALSE, col_names = TRUE)

  # results --------------------------------------------------------------------
  cat("############################### RESULTS ###############################\n")
  message("ind.approach: ", ind.approach)
  message("ind.threshold: ", ">= ", ind.threshold, " ", threshold.id)
  message("prob.pop.threshold: ", prob.pop.threshold)
  message("The number of markers (SNP/LOCUS) removed by the Individual filter:\n", markers.blacklist)
  message("The number of markers (SNP/LOCUS) before -> after the Individual filter:\n", markers.before, " -> ", markers.after)
  timing <- proc.time() - timing
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("############################## completed ##############################\n")
  res <- list()
  res$tidy.filtered.ind <- filter
  res$whitelist.markers <- whitelist.markers
  res$blacklist.markers <- blacklist.markers
  res$strata <- strata
  res$filters.parameters <- filters.parameters
  res$plot.ind.threshold <- plot.ind.threshold
  res$ind.threshold.helper.table <- ind.threshold.helper.table
  return(res)
  }




# filter_population-------------------------------------------------------------

#' @name filter_population
#' @title Population filter
#' @description Filter markers based on populations genotyped. Use a tidy data
#' set (long format) of any of these file format:
#' vcf, plink (tped/tfam), stacks haplotype file, genind,
#' genepop, data frame in wide format. The function uses
#' \code{\link[radiator]{tidy_genomic_data}} and
#' \code{\link[radiator]{tidy_wide}} to load the file. For filtering
#' The threshold can be a fixed number of population, a proportion or a percentage.

# Most arguments are inherited from tidy_genomic_data
#' @inheritParams tidy_genomic_data
#' @inheritParams read_strata
#' @inheritParams radiator_common_arguments


#' @param interactive.filter (optional, logical) Do you want the filtering session to
#' be interactive. With default: \code{interactive.filter == TRUE}, the user is
#' asked to see figures of distribution before making decisions for filtering.

#' @param pop.threshold The population threshold, proportion, percentage or
#' number e.g. 0.70, 70, 15.
#' Default: \code{pop.threshold = 100}.

#' @param percent Is the threshold a percentage? TRUE or FALSE.
#' This argument is necessary to distinguish percentage from integer population
#' threshold (e.g. 50 percent or 50 populations).
#' Default: \code{percent = TRUE}.

#' @param filename (optional) The function uses \code{\link[fst]{write.fst}},
#' to write the tidy data frame in
#' the folder created in the working directory. The file extension appended to
#' the \code{filename} provided is \code{.rad}.
#' Default: \code{filename = NULL}.

#' @rdname filter_population
#' @export
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot stat_summary labs theme element_blank element_text geom_jitter scale_colour_manual scale_y_reverse theme_light geom_bar facet_grid geom_histogram geom_bar aes_string scale_fill_manual theme_bw stat_smooth geom_boxplot ggsave
#' @importFrom stringi stri_join stri_replace_all_fixed stri_sub stri_detect_fixed
#' @importFrom dplyr select distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join funs
#' @importFrom readr write_tsv
#' @importFrom tibble data_frame has_name

#' @details
#' \strong{Interactive version}
#'
#' There is 2 steps in the interactive version to visualize and filter
#' the data based on the population representation:
#'
#' Step 1. Impact of population threshold on marker discovery
#'
#' Step 2. Choose the filtering threshold


#' @return With \code{interactive.filter = FALSE}, a list in the global environment,
#' with 7 objects:
#' \enumerate{
#' \item $tidy.filtered.ind
#' \item $whitelist.markers
#' \item $blacklist.markers
#' \item $strata
#' \item $filters.parameters
#' }
#'
#' With \code{interactive.filter = TRUE}, a list with 2 additionnal objects is created.
#' \enumerate{
#' \item $plot.pop.threshold
#' \item $pop.threshold.helper.table
#' }
#'
#' The object can be isolated in separate object outside the list by
#' following the example below.

#' @examples
#' \dontrun{
#' turtle.pop <- filter_population(
#' data = "turtle.vcf",
#' strata = "turtle.strata.tsv",
#' pop.thresholds = 100,
#' percent = TRUE,
#' filename = "tidy.data.turtle.tsv"
#' )

#'
#' #If interactive.filter = TRUE, a list is created and to view the filtered tidy data:
#' tidy.data <- turtle.ind$tidy.filtered.ind
#'
#' #Inside the same list, to isolate the blacklist.genotypes:
#' bg <- turtle.ind$blacklist.genotypes
#'
#' # The remaining argument are used in tidy_genomic_data during import and allow
#' # the user to apply some filtering or selection before doing the filtering.
#' }


filter_population <- function(
  interactive.filter = TRUE,
  data,
  pop.threshold = 100,
  percent = TRUE,
  filename = NULL,
  strata = NULL,
  parallel.core = parallel::detectCores() - 1,
  ...
) {
  cat("#######################################################################\n")
  cat("##################### radiator::filter_population #######################\n")
  cat("#######################################################################\n")
  opt.change <- getOption("width")
  options(width = 70)
  timing <- proc.time()
  # manage missing arguments -----------------------------------------------------
  if (missing(data)) rlang::abort("Input file missing")
  if (!is.null(pop.levels) & is.null(pop.labels)) {
    pop.levels <- stringi::stri_replace_all_fixed(pop.levels, pattern = " ", replacement = "_", vectorize_all = FALSE)
    pop.labels <- pop.levels
  }

  if (!is.null(pop.labels) & is.null(pop.levels)) rlang::abort("pop.levels is required if you use pop.labels")

  if (!is.null(pop.labels)) {
    if (length(pop.labels) != length(pop.levels)) rlang::abort("pop.labels and pop.levels must have the same length (number of groups)")
    pop.labels <- stringi::stri_replace_all_fixed(pop.labels, pattern = " ", replacement = "_", vectorize_all = FALSE)
  }

  if (!is.null(pop.select)) {
    pop.select <- stringi::stri_replace_all_fixed(pop.select, pattern = " ", replacement = "_", vectorize_all = FALSE)
  }

  # Message about steps taken during the process ---------------------------------
  if (interactive.filter) {
    message("Interactive mode: on")
    message("2 steps to visualize and filter the data based on the number of genotyped individuals:")
    message("Step 1. Impact of population threshold on marker discovery")
    message("Step 2. Choose the filtering threshold")
  }
  # Folder -------------------------------------------------------------------
  # Date and time
  file.date <- stringi::stri_replace_all_fixed(
    Sys.time(),
    pattern = " EDT", replacement = "") %>%
    stringi::stri_replace_all_fixed(
      str = .,
      pattern = c("-", " ", ":"), replacement = c("", "@", ""),
      vectorize_all = FALSE) %>%
    stringi::stri_sub(str = ., from = 1, to = 13)
  folder.extension <- stringi::stri_join("filter_population_", file.date, sep = "")
  path.folder <- stringi::stri_join(getwd(),"/", folder.extension, sep = "")
  dir.create(file.path(path.folder))

  message("\nFolder created: \n", folder.extension)
  file.date <- NULL #unused object

  # Filter parameter file ------------------------------------------------------
  message("Parameters used in this run are stored in a file")
  filters.parameters <- list.files(path = getwd(), pattern = "filters_parameters.tsv", full.names = TRUE)
  if (length(filters.parameters) == 0) {
    filters.parameters <- tibble::data_frame(FILTERS = as.character(), PARAMETERS = as.character(), VALUES = as.integer(), BEFORE = as.character(), AFTER = as.character(), BLACKLIST = as.integer(), UNITS = as.character(), COMMENTS = as.character())
    readr::write_tsv(x = filters.parameters, path = "filters_parameters.tsv", append = FALSE, col_names = TRUE)
    message("    Created a parameter file: filters_parameters.tsv")
  } else {
    message("    Using the filters parameters file: filters_parameters.tsv")
  }
  # File type detection----------------------------------------------------------
  data.type <- radiator::detect_genomic_format(data)

  # import data ----------------------------------------------------------------
  message("Importing data ...")
  input <- radiator::tidy_genomic_data(
    data = data,
    vcf.metadata = vcf.metadata,
    blacklist.id = blacklist.id,
    blacklist.genotype = blacklist.genotype,
    whitelist.markers = whitelist.markers,
    monomorphic.out = monomorphic.out,
    max.marker = max.marker,
    snp.ld = snp.ld,
    common.markers = common.markers,
    strata = strata,
    pop.levels = pop.levels,
    pop.labels = pop.labels,
    pop.select = pop.select,
    filename = NULL,
    parallel.core = parallel.core
  )

  # create a strata.df
  strata.df <- input %>%
    dplyr::select(INDIVIDUALS, POP_ID) %>%
    dplyr::distinct(INDIVIDUALS, .keep_all = TRUE)
  strata <- strata.df
  pop.levels <- levels(input$POP_ID)
  pop.labels <- pop.levels

  # prepare filter, table and figure--------------------------------------------
  pop.number <- dplyr::n_distinct(input$POP_ID) # number of pop

  filter.prep <- input %>%
    dplyr::filter(GT != "000000") %>%
    dplyr::distinct(MARKERS, POP_ID) %>%
    dplyr::group_by(MARKERS) %>%
    dplyr::tally(.) %>%
    dplyr::rename(POP_GENOTYPED = n)


  get_pop_number <- function(pop.indice, data) {
    marker.num <- dplyr::filter(
      .data = data,
      POP_GENOTYPED >= pop.indice) %>%
      dplyr::summarise(MARKER_NUMBER = sum(MARKER_NUMBER))
    res <- tibble::data_frame(
      POP_GENOTYPED = pop.indice,
      MARKER_NUMBER = marker.num$MARKER_NUMBER)
    return(res)
  }

  pop.threshold.helper.table <- filter.prep %>%
    dplyr::group_by(POP_GENOTYPED) %>%
    dplyr::summarise(MARKER_NUMBER = length(MARKERS))

  pop.threshold.helper.table <- purrr::map_df(.x = 1:pop.number, .f = get_pop_number, data = pop.threshold.helper.table)

  # Step 1. Impact of population threshold on marker discovery------------------
  if (interactive.filter) {
    message("Step 1. Impact of population threshold on marker discovery")
    message("    Inspect the plot and table produced showing the relationship
            between the number of markers and the number of population genotyped for
            the marker")
  }
  # line.graph <- as.character(readLines(n = 1))
  # if (line.graph == "y") {
  # Set the breaks for the figure
  max.markers <- dplyr::n_distinct(input$MARKERS)

  #Function to replace packageplyr round_any
  rounder <- function(x, accuracy, f = round) {
    f(x / accuracy) * accuracy
  }

  # max.markers <- 658
  if (max.markers >= 1000) {
    y.breaks.by <- rounder(max.markers/10, 100, ceiling)
    y.breaks.max <- rounder(max.markers, 1000, ceiling)
    y.breaks <- seq(0, y.breaks.max, by = y.breaks.by)
  } else {
    y.breaks.by <- rounder(max.markers/10, 10, ceiling)
    y.breaks.max <- rounder(max.markers, 100, ceiling)
    y.breaks <- seq(0, y.breaks.max, by = y.breaks.by)
  }

  plot.pop.threshold <- ggplot2::ggplot(
    pop.threshold.helper.table,
    ggplot2::aes(x = POP_GENOTYPED, y = MARKER_NUMBER)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2, shape = 21, fill = "white") +
    ggplot2::scale_x_continuous(name = "Number of populations genotyped") +
    ggplot2::scale_y_continuous(
      name = "Number of markers", breaks = y.breaks, limits = c(0, y.breaks.max)) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica", hjust = 1, vjust = 0.5),
      strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
    )
  # plot.pop.threshold
  if (interactive.filter) print(plot.pop.threshold)
  # save
  ggplot2::ggsave(stringi::stri_join(path.folder, "/plot.pop.threshold.pdf"), width = length(pop.threshold.helper.table$POP_GENOTYPED)*2, height = 10, dpi = 600, units = "cm", useDingbats = FALSE)
  ggplot2::ggsave(stringi::stri_join(path.folder, "/plot.pop.threshold.png"), width = length(pop.threshold.helper.table$POP_GENOTYPED)*2, height = 10, dpi = 300, units = "cm")
  message("\n    2 versions (pdf and png) of the line graph of populations threshold
          on marker discovery were writted in the folder")


  # Helper table for individual thresholds
  if (interactive.filter) {
    message("\n    A table (pop.threshold.helper.table) to help you view
            the relation between population threshold and marker discovery, with your dataset:")
    print(pop.threshold.helper.table)
  }
  readr::write_tsv(x = pop.threshold.helper.table, path = stringi::stri_join(path.folder, "/", "pop.threshold.helper.table.tsv"))
  message("    pop.threshold.helper.table also written in the folder")

  # Step 2. Choose the filtering threshold -------------------------------------
  if (interactive.filter) {
    message("Step 2. Choose the filtering threshold.")
    message("    Enter the population threshold (number, proportion or percentage).
            e.g. enter 10 (for 10 populations), 0.8 for proportion and 80 for percentage.")
    pop.threshold <- as.numeric(readLines(n = 1))

    message("    The value you just enter for the pop threshold, is it a percentage? (TRUE/FALSE)")
    percent <- as.logical(readLines(n = 1))
  }

  # Filtering ------------------------------------------------------------------
  # pop.number <- dplyr::n_distinct(input$POP_ID) # number of pop

  filter.prep <- filter.prep %>%
    dplyr::mutate(PERCENT = ceiling(POP_GENOTYPED/pop.number*100))

  if (!percent & pop.threshold >= 1) {
    message("Using a fixed threshold")
    threshold.id <- "(pop.)"
    # pop.threshold <- 12
    filter <- filter.prep %>%
      dplyr::ungroup(.) %>%
      dplyr::filter(POP_GENOTYPED >= pop.threshold) %>%
      dplyr::distinct(MARKERS)
  } else if (!percent & stringi::stri_detect_fixed(pop.threshold, ".") & pop.threshold < 1) {
    message("Using a proportion threshold...")
    threshold.id <- "(proportion)"
    # pop.threshold <- 0.6
    filter <- filter.prep %>%
      dplyr::ungroup(.) %>%
      dplyr::filter(PERCENT >= pop.threshold*100) %>%
      dplyr::distinct(MARKERS)
  } else {
    message("Using a percentage threshold...")
    threshold.id <- "(percent)"
    # pop.threshold <- 100
    filter <- filter.prep %>%
      dplyr::ungroup(.) %>%
      dplyr::filter(PERCENT >= pop.threshold) %>%
      dplyr::distinct(MARKERS)
  }

  # Apply the filter to the tidy data
  filter <- dplyr::left_join(x = filter, input, by = "MARKERS")


  # Update filters.parameters SNP ----------------------------------------------
  if (data.type == "haplo.file") {
    snp.before <- as.character("NA")
    snp.after <- as.character("NA")
    snp.blacklist <- as.character("NA")
  } else {
    snp.before <- as.integer(dplyr::n_distinct(input$MARKERS))
    snp.after <- as.integer(dplyr::n_distinct(filter$MARKERS))
    snp.blacklist <- as.integer(dplyr::n_distinct(input$MARKERS) - dplyr::n_distinct(filter$MARKERS))
  }

  if (tibble::has_name(input, "LOCUS")) {
    locus.before <- as.integer(dplyr::n_distinct(input$LOCUS))
    locus.after <- as.integer(dplyr::n_distinct(filter$LOCUS))
    locus.blacklist <- as.integer(dplyr::n_distinct(input$LOCUS) - dplyr::n_distinct(filter$LOCUS))
  } else {
    locus.before <- as.character("NA")
    locus.after <- as.character("NA")
    locus.blacklist <- as.character("NA")
  }


  markers.before <- stringi::stri_join(snp.before, locus.before, sep = "/")
  markers.after <- stringi::stri_join(snp.after, locus.after, sep = "/")
  markers.blacklist <- stringi::stri_join(snp.blacklist, locus.blacklist, sep = "/")

  filters.parameters <- tibble::data_frame(
    FILTERS = c("Populations"),
    PARAMETERS = c("pop.threshold"),
    VALUES = c(stringi::stri_join(pop.threshold, " ", threshold.id)),
    BEFORE = c(markers.before),
    AFTER = c(markers.after),
    BLACKLIST = c(markers.blacklist),
    UNITS = c("SNP/LOCUS"),
    COMMENTS = c("")
  )
  readr::write_tsv(x = filters.parameters, path = "filters_parameters.tsv", append = TRUE, col_names = FALSE)

  # saving tidy data
  if (!is.null(filename)) {
    tidy.name <- stringi::stri_join(filename, ".rad")
    message("Writing the filtered tidy data set: ", tidy.name)
    write_rad(data = filter, path = file.path(path.folder, tidy.name))
  }

  # saving whitelist
  message("Writing the whitelist of markers: whitelist.markers.pop.tsv")

  if (tibble::has_name(input, "CHROM")) {
    whitelist.markers <- dplyr::ungroup(filter) %>%
      dplyr::distinct(CHROM, LOCUS, POS)
  } else {
    whitelist.markers <- dplyr::ungroup(filter) %>%
      dplyr::distinct(MARKERS)
  }
  readr::write_tsv(whitelist.markers, stringi::stri_join(path.folder, "/", "whitelist.markers.pop.tsv"), append = FALSE, col_names = TRUE)

  # saving blacklist
  message("Writing the blacklist of markers: blacklist.markers.pop.tsv")
  if (tibble::has_name(input, "CHROM")) {
    blacklist.markers <- dplyr::ungroup(input) %>%
      dplyr::distinct(CHROM, LOCUS, POS) %>%
      dplyr::anti_join(whitelist.markers, by = c("CHROM", "LOCUS", "POS"))
  } else {
    blacklist.markers <- dplyr::ungroup(input) %>%
      dplyr::distinct(MARKERS) %>%
      dplyr::anti_join(whitelist.markers, by = "MARKERS")
  }
  readr::write_tsv(blacklist.markers, stringi::stri_join(path.folder, "/", "blacklist.markers.pop.tsv"), append = FALSE, col_names = TRUE)

  # results --------------------------------------------------------------------
  cat("############################### RESULTS ###############################\n")
  message("pop.threshold: ", ">= ", pop.threshold)
  message("The number of markers (SNP/LOCUS) removed by the Population filter:\n", markers.blacklist)
  message("The number of markers (SNP/LOCUS) before -> after the Population filter:\n", markers.before, " -> ", markers.after)
  timing <- proc.time() - timing
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("############################## completed ##############################\n")
  res <- list()
  res$tidy.filtered.pop <- filter
  res$whitelist.markers <- whitelist.markers
  res$blacklist.markers <- blacklist.markers
  res$strata <- strata
  res$filters.parameters <- filters.parameters
  res$plot.pop.threshold <- plot.pop.threshold
  res$pop.threshold.helper.table <- pop.threshold.helper.table
  return(res)
  }
