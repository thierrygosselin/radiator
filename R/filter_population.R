#' @name filter_population
#' @title Population filter
#' @description Filter markes based on populations genotyped. Use a tidy data 
#' set (long format) of any of these file format: 
#' vcf, plink (tped/tfam), stacks haplotype file, genind, 
#' genepop, data frame in wide format. The function uses 
#' \code{\link[radiator]{tidy_genomic_data}} and 
#' \code{\link[radiator]{tidy_wide}} to load the file. For filtering
#' The threshold can be a fixed number of population, a proportion or a percentage.

# Most arguments are inherited from tidy_genomic_data
#' @inheritParams tidy_genomic_data

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

#' @param filename (optional) Name of the filtered tidy data frame file 
#' written to the working directory (ending with \code{.tsv})
#' Default: \code{filename = NULL}.

#' @rdname filter_population
#' @export
#' @import ggplot2
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
  data,
  vcf.metadata = FALSE,
  interactive.filter = TRUE,
  pop.threshold = 100,
  percent = TRUE,
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
  pop.select = NULL
) {
  
  
  cat("#######################################################################\n")
  cat("##################### radiator::filter_population #######################\n")
  cat("#######################################################################\n")
  
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
    message("2 steps to visualize and filter the data based on the number of genotyped individuals:")
    message("Step 1. Impact of population threshold on marker discovery")
    message("Step 2. Choose the filtering threshold")
    
    # Folder -------------------------------------------------------------------
    # Get date and time to have unique filenaming
    file.date <- stringi::stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "")
    file.date <- stringi::stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
    file.date <- stringi::stri_sub(file.date, from = 1, to = 13)
    
    path.folder <- stringi::stri_join(getwd(),"/", "filter_population_", file.date, sep = "")
    dir.create(file.path(path.folder))
    
    message(stringi::stri_join("Folder created: \n", path.folder))
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
    filename = NULL
  )
  
  # create a strata.df
  strata.df <- input %>% 
    dplyr::select(INDIVIDUALS, POP_ID) %>% 
    dplyr::distinct(INDIVIDUALS, .keep_all = TRUE)
  strata <- strata.df
  pop.levels <- levels(input$POP_ID)
  pop.labels <- pop.levels
  
  # prepare filter, table and figure----------------------------------------------
  filter.prep <- input %>% 
    dplyr::filter(GT != "000000") %>% 
    dplyr::distinct(MARKERS, POP_ID) %>% 
    dplyr::group_by(MARKERS) %>% 
    dplyr::tally(.) %>% 
    dplyr::rename(POP_GENOTYPED = n)
  
  pop.threshold.helper.table <- filter.prep %>% 
    dplyr::group_by(POP_GENOTYPED) %>% 
    dplyr::summarise(MARKER_NUMBER = length(MARKERS))
  
  
  
  # Step 1. Impact of population threshold on marker discovery------------------
  if (interactive.filter) {
    message("Step 1. Impact of population threshold on marker discovery")
    message("Show the line graph to inspect the change in the number of markers in relation to the population percentage thresholds (y/n)): ")
    line.graph <- as.character(readLines(n = 1))
    if (line.graph == "y") {
      # Set the breaks for the figure
      max.markers <- dplyr::n_distinct(input$MARKERS)
      
      #Function to replace plyr::round_any
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
      
      plot.pop.threshold <- ggplot(pop.threshold.helper.table, aes(x = POP_GENOTYPED, y = MARKER_NUMBER)) +
        geom_line() +
        geom_point(size = 2, shape = 21, fill = "white") +
        scale_x_continuous(name = "Number of populations genotyped") +
        scale_y_continuous(name = "Number of markers", breaks = y.breaks, limits = c(0, y.breaks.max)) +
        theme(
          axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
          axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
          axis.text.x = element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5),
          strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold")
        )
      # plot.pop.threshold
      print(plot.pop.threshold)
      # save
      ggsave(stringi::stri_join(path.folder, "/plot.pop.threshold.pdf"), width = length(pop.threshold.helper.table$POP_GENOTYPED)*2, height = 10, dpi = 600, units = "cm", useDingbats = F)
      ggsave(stringi::stri_join(path.folder, "/plot.pop.threshold.png"), width = length(pop.threshold.helper.table$POP_GENOTYPED)*2, height = 10, dpi = 300, units = "cm")
      message(stringi::stri_join("2 versions (pdf and png) of the line graph of populations threshold on marker discovery were saved in this directory: \n", path.folder))
      
      
      # Helper table for individual thresholds
      message("A table (pop.threshold.helper.table) to help you view the relation between population threshold and marker discovery, with your dataset:")
      print(pop.threshold.helper.table)
      readr::write_tsv(x = pop.threshold.helper.table, path = stringi::stri_join(path.folder, "/", "pop.threshold.helper.table.tsv"))
      message(stringi::stri_join("pop.threshold.helper.table was written in this directory: \n", path.folder))
    }
  }
  
  # Step 2. Choose the filtering threshold -------------------------------------
  if (interactive.filter) {
    message("Step 2. Choose the filtering threshold.")
    message("Enter the population threshold (number, proportion or percentage).
e.g. enter 10 (for 10 populations), 0.8 for proportion and 80 for percentage.")
    pop.threshold <- as.numeric(readLines(n = 1))
    
    message("The value you just enter for the pop threshold, is it a percentage? (TRUE/FALSE)")
    percent <- as.character(readLines(n = 1))
  }
  
  # Filtering ------------------------------------------------------------------
  pop.number <- dplyr::n_distinct(input$POP_ID) # number of pop
  
  filter.prep <- filter.prep %>% 
    mutate(PERCENT = ceiling(POP_GENOTYPED/pop.number*100))
  
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
  if (data.type == "vcf.file") {
    # markers.chrom.locus.pos <- input %>% distinct(MARKERS, LOCUS)
    # filter <- left_join(filter, markers.chrom.locus.pos, by = "MARKERS")
    snp.before <- as.integer(dplyr::n_distinct(input$MARKERS))
    snp.after <- as.integer(dplyr::n_distinct(filter$MARKERS))
    snp.blacklist <- as.integer(dplyr::n_distinct(input$MARKERS) - dplyr::n_distinct(filter$MARKERS))
    locus.before <- as.integer(dplyr::n_distinct(input$LOCUS))
    locus.after <- as.integer(dplyr::n_distinct(filter$LOCUS))
    locus.blacklist <- as.integer(dplyr::n_distinct(input$LOCUS) - dplyr::n_distinct(filter$LOCUS))
  } else if (data.type == "haplo.file") {
    snp.before <- as.character("NA")
    snp.after <- as.character("NA")
    snp.blacklist <- as.character("NA")
    locus.before <- as.integer(dplyr::n_distinct(input$MARKERS))
    locus.after <- as.integer(dplyr::n_distinct(filter$MARKERS))
    locus.blacklist <- as.integer(dplyr::n_distinct(input$MARKERS) - dplyr::n_distinct(filter$MARKERS))
  } else {
    snp.before <- as.integer(dplyr::n_distinct(input$MARKERS))
    snp.after <- as.integer(dplyr::n_distinct(filter$MARKERS))
    snp.blacklist <- as.integer(dplyr::n_distinct(input$MARKERS) - dplyr::n_distinct(filter$MARKERS))
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
    message("Writing the filtered tidy data set in your working directory...")
    # if (!is.null(save.feather)) {
    # feather::write_feather(filter, stringi::stri_replace_all_fixed(filename, pattern = ".tsv", replacement = "_feather.tsv", vectorize_all = TRUE))
    # } else {
    readr::write_tsv(filter, filename, append = FALSE, col_names = TRUE)
    # }
  }
  
  # saving whitelist
  message("Writing the whitelist of markers in your working directory\nwhitelist.markers.pop.tsv")
  
  if (tibble::has_name(input, "CHROM")) {
    whitelist.markers <- dplyr::ungroup(filter) %>%
      dplyr::distinct(CHROM, LOCUS, POS)
  } else {
    whitelist.markers <- dplyr::ungroup(filter) %>%
      dplyr::distinct(MARKERS)
  }
  readr::write_tsv(whitelist.markers, "whitelist.markers.pop.tsv", append = FALSE, col_names = TRUE)
  
  
  # saving blacklist
  message("Writing the blacklist of markers in your working directory\nblacklist.markers.pop.tsv")
  if (tibble::has_name(input, "CHROM")) {
    blacklist.markers <- dplyr::ungroup(input) %>%
      dplyr::distinct(CHROM, LOCUS, POS) %>% 
      dplyr::anti_join(whitelist.markers, by = c("CHROM", "LOCUS", "POS"))
  } else {
    blacklist.markers <- dplyr::ungroup(input) %>%
      dplyr::distinct(MARKERS) %>% 
      dplyr::anti_join(whitelist.markers, by = "MARKERS")
  }
  readr::write_tsv(blacklist.markers, "blacklist.markers.pop.tsv", append = FALSE, col_names = TRUE)
  
  # results --------------------------------------------------------------------
  cat("############################### RESULTS ###############################\n")
  message(stringi::stri_join("pop.threshold: ", ">= ", pop.threshold))
  if (data.type != "vcf.file" & data.type != "haplo.file") {
    message(stringi::stri_join("The number of markers removed by the Population filter: ", dplyr::n_distinct(input$MARKERS) - dplyr::n_distinct(filter$MARKERS)))
    message("The number of markers before -> after the Population filter")
    message(stringi::stri_join("SNP: ", as.integer(dplyr::n_distinct(input$MARKERS)), " -> ", as.integer(dplyr::n_distinct(filter$MARKERS))))
  } else if (data.type == "vcf.file") {
    message(stringi::stri_join("The number of markers removed by the Population filter:\nSNP: ", dplyr::n_distinct(input$POS) - dplyr::n_distinct(filter$POS), "\nLOCUS: ", dplyr::n_distinct(input$LOCUS) - dplyr::n_distinct(filter$LOCUS)))
    message("The number of markers before -> after the Population filter")
    message(stringi::stri_join("SNP: ", as.integer(dplyr::n_distinct(input$POS)), " -> ", as.integer(dplyr::n_distinct(filter$POS))))
    message(stringi::stri_join("LOCUS: ", as.integer(dplyr::n_distinct(input$LOCUS)), " -> ", as.integer(dplyr::n_distinct(filter$LOCUS))))
  } else {# for haplotype file
    message(stringi::stri_join("The number of markers/locus removed by the Population filter: ", dplyr::n_distinct(input$MARKERS) - dplyr::n_distinct(filter$MARKERS)))
    message("The number of markers before -> after the Population filter")
    message(stringi::stri_join("MARKERS/LOCUS: ", as.integer(dplyr::n_distinct(input$MARKERS)), " -> ", as.integer(dplyr::n_distinct(filter$MARKERS))))
  }
  cat("############################## completed ##############################\n")
  res <- list()
  res$tidy.filtered.pop <- filter
  res$whitelist.markers <- whitelist.markers
  res$blacklist.markers <- blacklist.markers
  res$strata <- strata
  res$filters.parameters <- filters.parameters
  if (interactive.filter) {
    res$plot.pop.threshold <- plot.pop.threshold
    res$pop.threshold.helper.table <- pop.threshold.helper.table
  }
  return(res)
}
