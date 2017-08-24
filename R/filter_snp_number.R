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

# Most arguments are inherited from tidy_genomic_data
#' @inheritParams tidy_genomic_data

#' @param interactive.filter (optional, logical) Do you want the filtering session to
#' be interactive. With default: \code{interactive.filter == TRUE}, the user is
#' asked to see figures of distribution before making filtering decisions.

#' @param max.snp.number (integer) This is best decide after viewing the figures.
#' If the argument is set to 2, locus with 3 and more SNPs will be blacklisted.

#' @param filename (optional) Name of the filtered tidy data frame file
#' written to the working directory (ending with \code{.tsv})
#' Default: \code{filename = NULL}.


#' @rdname filter_snp_number
#' @export
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot stat_summary labs theme element_blank element_text geom_jitter scale_colour_manual scale_y_reverse theme_light geom_bar facet_grid geom_histogram geom_bar aes_string scale_fill_manual theme_bw stat_smooth geom_boxplot ggsave
#' @importFrom stringi stri_join stri_replace_all_fixed stri_sub stri_detect_fixed
#' @importFrom dplyr select distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join funs
#' @importFrom readr write_tsv
#' @importFrom tibble data_frame has_name

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
  vcf.metadata = FALSE,
  interactive.filter = TRUE,
  max.snp.number,
  filename = NULL,
  blacklist.id = NULL,
  blacklist.genotype = NULL,
  whitelist.markers = NULL,
  monomorphic.out = TRUE,
  max.marker = NULL,
  snp.ld = NULL,
  common.markers = FALSE,
  strata = NULL,
  pop.levels = NULL,
  pop.labels = NULL,
  pop.select = NULL,
  parallel.core = parallel::detectCores() - 1
) {
  cat("#######################################################################\n")
  cat("#################### radiator::filter_snp_number ######################\n")
  cat("#######################################################################\n")
  opt.change <- getOption("width")
  options(width = 70)
  timing <- proc.time()

  # list to store results
  res <- list()

  # manage missing arguments -----------------------------------------------------
  if (missing(data)) stop("Input file missing")

  if (!is.null(pop.levels) & is.null(pop.labels)) {
    pop.levels <- stringi::stri_replace_all_fixed(pop.levels, pattern = " ", replacement = "_", vectorize_all = FALSE)
    pop.labels <- pop.levels
  }
  if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")

  if (!is.null(pop.labels)) {
    if (length(pop.labels) != length(pop.levels)) stop("pop.labels and pop.levels must have the same length (number of groups)")
    pop.labels <- stringi::stri_replace_all_fixed(pop.labels, pattern = " ", replacement = "_", vectorize_all = FALSE)
  }

  if (!is.null(pop.select)) {
    pop.select <- stringi::stri_replace_all_fixed(pop.select, pattern = " ", replacement = "_", vectorize_all = FALSE)
  }
  # Message about steps taken during the process ---------------------------------
  if (interactive.filter) {
    message("Interactive mode: on")
    message("2 steps to visualize and filter the data based on the number of SNP on the read/locus:")
    message("Step 1. Impact of SNP number per read/locus (on individual genotypes and locus/snp number potentially filtered)")
    message("Step 2. Choose the filtering thresholds")

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
    folder.extension <- stringi::stri_join("filter_snp_number_", file.date, sep = "")
    path.folder <- stringi::stri_join(getwd(),"/", folder.extension, sep = "")
    dir.create(file.path(path.folder))

    message("Folder created: \n", folder.extension)
    file.date <- NULL #unused object
  } else {
    path.folder <- getwd()
  }

  # Filter parameter file ------------------------------------------------------
  message("Parameters used in this run will be store in a file")
  filters.parameters <- list.files(path = getwd(), pattern = "filters_parameters.tsv", full.names = TRUE)
  if (length(filters.parameters) == 0) {
    filters.parameters <- tibble::data_frame(FILTERS = as.character(), PARAMETERS = as.character(), VALUES = as.integer(), BEFORE = as.character(), AFTER = as.character(), BLACKLIST = as.integer(), UNITS = as.character(), COMMENTS = as.character())
    readr::write_tsv(x = filters.parameters, path = "filters_parameters.tsv", append = FALSE, col_names = TRUE)
    message("Created a file to store filters parameters: filters_parameters.tsv")
  } else {
    message("Using the filters parameters file found in the directory: \nfilters_parameters.tsv")
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
    parallel.core = parallel.core,
    filename = NULL
  )

  # Check that required info is present in data: snp and locus

  locus.check <- tibble::has_name(input, "LOCUS")
  snp.check <- tibble::has_name(input, "POS")
  if (!locus.check || !snp.check) stop("This filter requires SNP (POS) and LOCUS information (columns) in the dataset (e.g. found in VCF files)")

  snp.before <- as.integer(dplyr::n_distinct(input$MARKERS))
  locus.before <- as.integer(dplyr::n_distinct(input$LOCUS))

  # create a strata.df
  strata.df <- input %>%
    dplyr::select(INDIVIDUALS, POP_ID) %>%
    dplyr::distinct(INDIVIDUALS, .keep_all = TRUE)
  pop.levels <- levels(input$POP_ID)
  pop.labels <- pop.levels

  # some stats
  ind.total <- dplyr::n_distinct(input$INDIVIDUALS) # total number of individuals
  pop.number <- dplyr::n_distinct(input$POP_ID) # number of pop

  # prepare the info
  res$snp.num.markers <- dplyr::distinct(.data = input, LOCUS, POS) %>%
    dplyr::group_by(LOCUS) %>%
    dplyr::tally(.) %>%
    dplyr::rename(SNP_NUMBER = n)

  res$number.snp.reads.plot <- ggplot2::ggplot(res$snp.num.markers,
                                           ggplot2::aes(factor(SNP_NUMBER))) +
    ggplot2::geom_bar() +
    ggplot2::labs(x = "Number of SNP per locus/reads") +
    ggplot2::labs(y = "Distribution (number of locus)") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
      legend.title = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
      legend.text = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
      strip.text.x = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"))

  if (interactive.filter) {
    print(res$number.snp.reads.plot)
    message("Based on the plot, choose the threshold in maximum number of SNP per locus allowed (an integer): ")
    max.snp.number <- as.integer(readLines(n = 1))
  }

  # save
  ggplot2::ggsave(stringi::stri_join(path.folder, "/number.snp.locus.plot.pdf"), width = 20, height = 10, dpi = 600, units = "cm", useDingbats = FALSE)
  ggplot2::ggsave(stringi::stri_join(path.folder, "/number.snp.locus.plot.png"), width = 20, height = 10, dpi = 300, units = "cm")
  message("2 versions (pdf and png) of the histogram of the number of SNP per locus were saved in this directory: \n", path.folder)


  # Filtering ------------------------------------------------------------------
  message("Filtering data...")
  res$whitelist.markers <- res$snp.num.markers %>%
    dplyr::filter(SNP_NUMBER <= max.snp.number) %>%
    dplyr::select(LOCUS) %>%
    dplyr::arrange(LOCUS)

  # Apply the filter to the tidy data
  res$tidy.filtered.snp.number <- suppressWarnings(
    dplyr::semi_join(input, res$whitelist.markers, by = "LOCUS"))

  # saving whitelist
  message("Writing the whitelist of markers in your working directory\nwhitelist.markers.snp.number.tsv")

  if (tibble::has_name(res$tidy.filtered.snp.number, "CHROM")) {
    res$whitelist.markers <- dplyr::ungroup(res$tidy.filtered.snp.number) %>%
      dplyr::select(MARKERS, CHROM, LOCUS, POS) %>%
      dplyr::distinct(CHROM, LOCUS, POS, .keep_all = TRUE)
  } else {
    res$whitelist.markers <- dplyr::ungroup(res$tidy.filtered.snp.number) %>%
      dplyr::distinct(MARKERS)
  }
  readr::write_tsv(res$whitelist.markers, stringi::stri_join(path.folder, "/whitelist.markers.snp.number.tsv"), append = FALSE, col_names = TRUE)

  # saving blacklist
  if (tibble::has_name(res$tidy.filtered.snp.number, "CHROM")) {
    res$blacklist.markers <- dplyr::ungroup(input) %>%
      dplyr::select(MARKERS, CHROM, LOCUS, POS) %>%
      dplyr::distinct(CHROM, LOCUS, POS, .keep_all = TRUE) %>%
      dplyr::anti_join(res$whitelist.markers, by = c("CHROM", "LOCUS", "POS"))
  } else {
    res$blacklist.markers <- dplyr::ungroup(input) %>%
      dplyr::distinct(MARKERS) %>%
      dplyr::anti_join(res$whitelist.markers, by = "MARKERS")
  }
  if (length(res$blacklist.markers$MARKERS) > 0) {
    message("Writing the blacklist of markers in your working directory\nblacklist.markers.snp.number.tsv")
    readr::write_tsv(res$blacklist.markers, stringi::stri_join(path.folder, "/blacklist.markers.snp.number.tsv"), append = FALSE, col_names = TRUE)
  }

  input <- NULL

  # Update filters.parameters SNP ----------------------------------------------
  snp.after <- as.integer(dplyr::n_distinct(res$tidy.filtered.snp.number$MARKERS))
  snp.blacklist <- as.integer(snp.before - snp.after)
  locus.after <- as.integer(dplyr::n_distinct(res$tidy.filtered.snp.number$LOCUS))
  locus.blacklist <- as.integer(locus.before - locus.after)

  markers.before <- stringi::stri_join(snp.before, locus.before, sep = "/")
  markers.after <- stringi::stri_join(snp.after, locus.after, sep = "/")
  markers.blacklist <- stringi::stri_join(snp.blacklist, locus.blacklist, sep = "/")

  res$filters.parameters <- tibble::data_frame(
    FILTERS = c("SNP number per locus"),
    PARAMETERS = c("max.snp.number"),
    VALUES = c(stringi::stri_join("<=", max.snp.number)),
    BEFORE = c(markers.before),
    AFTER = c(markers.after),
    BLACKLIST = c(markers.blacklist),
    UNITS = c("SNP/LOCUS"),
    COMMENTS = c("")
  )
  readr::write_tsv(x = res$filters.parameters, path = "filters_parameters.tsv", append = TRUE, col_names = FALSE)

  # saving tidy data
  if (!is.null(filename)) {
    message("Writing the filtered tidy data set in your working directory...")
    # if (!is.null(save.feather)) {
    # feather::write_feather(filter, stri_replace_all_fixed(filename, pattern = ".tsv", replacement = "_feather.tsv", vectorize_all = TRUE))
    # } else {
    readr::write_tsv(res$tidy.filtered.snp.number,
                     filename, append = FALSE, col_names = TRUE)
    # }
  }


  # results --------------------------------------------------------------------
  cat("############################### RESULTS ###############################\n")
    message("The number of markers blacklisted (SNP/LOCUS): ", markers.blacklist)
    message("The number of markers before -> after filter_snp_number (SNP/LOCUS)")
    message(markers.before, " -> ", markers.after)
    timing <- proc.time() - timing
    message("\nComputation time: ", round(timing[[3]]), " sec")
    cat("############################## completed ##############################\n")
  return(res)
}
