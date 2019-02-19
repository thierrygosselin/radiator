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
