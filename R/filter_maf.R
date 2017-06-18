# Minor Allele Frequency
#' @name filter_maf
#' @title MAF filter
#' @description Minor Allele Frequency filter from a tidy data set (long format)
#' of any of these file format:
#' vcf, plink (tped/tfam), stacks haplotype file, genind,
#' genepop, data frame in wide format. The function uses
#' \code{\link[radiator]{tidy_genomic_data}} and
#' \code{\link[radiator]{tidy_wide}} to load the file.

# Most arguments are inherited from tidy_genomic_data
#' @inheritParams tidy_genomic_data

#' @param interactive.filter (optional, logical) Do you want the filtering session to
#' be interactive. With default: \code{interactive.filter == TRUE}, the user is
#' asked to see figures of distribution before making decisions for filtering.

#' @param filename (optional) Name of the filtered tidy data frame file
#' written to the working directory (ending with \code{.tsv})
#' Default: \code{filename = NULL}.

# @param save.feather (optional) Use the package
# \href{https://github.com/wesm/feather}{feather} to save the data frame (very fast).
# Default: \code{save.feather = NULL}.

#' @rdname filter_maf
#' @export
# @importFrom feather write_feather

#' @import ggplot2
#' @importFrom stringi stri_join stri_replace_all_fixed stri_sub
#' @importFrom dplyr select distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join funs summarise_at bind_rows
#' @importFrom readr write_tsv
#' @importFrom tibble data_frame has_name
#' @importFrom tidyr complete gather unite spread nesting

#' @details To help choose a threshold for the local and global MAF, look
#' at the function \link{diagnostic_maf}, or use the interactive version.
#'
#' There is 3 steps in the interactive version:
#'
#' Step 1. Global MAF: Inspecting the MAF globally
#'
#' Step 2. Local MAF: Inspecting the MAF at the populations level
#'
#' Step 3. Filtering markers based on the different MAF arguments

#' @return With \code{interactive.filter = FALSE}, a list in the global environment,
#' with 7 objects:
#' \enumerate{
#' \item $tidy.filtered.maf
#' \item $whitelist.markers
#' \item $blacklist.markers
#' \item $blacklist.markers
#' \item $maf.data.thresholds
#' \item $strata
#' \item $filters.parameters
#' }
#'
#' With \code{interactive.filter = TRUE}, a list with 5 additionnal objects is created.
#' \enumerate{
#' \item $violinplot.maf.global <- violinplot.maf.global
#' \item $violinplot.maf.local
#' \item $plot.distribution.maf.locall
#' \item $maf.global.summary
#' \item $maf.helper.table
#' }
#'
#' The object can be isolated in separate object outside the list by
#' following the example below.

#' @examples
#' \dontrun{
#' # The minumum
#' turtle.maf <- filter_maf(
#' data = "turtle.vcf",
#' strata = "turtle.strata.tsv",
#' maf.thresholds = c(0.04, 0.02)
#' )
#' #This will use the default: interactive version,
#' maf.approach = "SNP",
#' maf.operator = "OR",
#' maf.pop.num.threshold = 1
#'
#'
#' #If interactive.filter = TRUE, a list is created and to view the filtered tidy data:
#' tidy.data <- turtle.maf$tidy.filtered.maf
#'
#' #Inside the same list, to isolate the blacklist.genotypes:
#' bg <- turtle.maf$blacklist.genotypes
#'
#' # The remaining argument are used in tidy_genomic_data during import and allow
#' # the user to apply some filtering or selection before doing the MAF filtering.
#' }

filter_maf <- function(
  data,
  vcf.metadata = FALSE,
  interactive.filter = TRUE,
  maf.approach = "SNP",
  maf.thresholds = NULL,
  maf.operator = "OR",
  maf.pop.num.threshold = 1,
  filename = NULL,
  #save.feather = NULL,
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
  cat("######################### radiator::filter_maf ##########################\n")
  cat("#######################################################################\n")

  timing <- proc.time()
  # manage missing arguments -----------------------------------------------------
  if (!interactive.filter & missing(maf.thresholds)) stop("The required maf.thresholds argument is missing")
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
    message("Interactive mode: on\n")
    message("3 steps to visualize and filter the data based on MAF:")
    message("Step 1. Global MAF: Inspecting the MAF globally")
    message("Step 2. Local MAF: Inspecting the MAF at the populations level")
    message("Step 3. Filtering markers based on the different MAF arguments\n\n")
  }
    # Folder -------------------------------------------------------------------
    # Get date and time to have unique filenaming
    file.date <- stringi::stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "")
    file.date <- stringi::stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
    file.date <- stringi::stri_sub(file.date, from = 1, to = 13)

    path.folder <- stringi::stri_join(getwd(),"/", "filter_maf_", file.date, sep = "")
    dir.create(file.path(path.folder))

    message(stringi::stri_join("Folder created: \n", path.folder))
    file.date <- NULL #unused object

  # Filter parameter file ------------------------------------------------------
  message("\nParameters used in this run will be store in a file")
  filters.parameters <- list.files(path = getwd(), pattern = "filters_parameters.tsv", full.names = TRUE)
  if (length(filters.parameters) == 0) {
    filters.parameters <- tibble::data_frame(FILTERS = as.character(), PARAMETERS = as.character(), VALUES = as.integer(), BEFORE = as.character(), AFTER = as.character(), BLACKLIST = as.integer(), UNITS = as.character(), COMMENTS = as.character())
    readr::write_tsv(x = filters.parameters, path = "filters_parameters.tsv", append = FALSE, col_names = TRUE)
    message("Created a file to store filters parameters: filters_parameters.tsv\n")
  } else {
    message("Using the filters parameters file found in the directory: \nfilters_parameters.tsv\n")
  }

  # File type detection----------------------------------------------------------
  data.type <- radiator::detect_genomic_format(data)

  if (data.type == "haplo.file") {
    message("With stacks haplotype file the maf.approach is automatically set to: haplotype")
    maf.approach <- "SNP"
    # confusing, but because the haplotpe file doesn't have snp info, only locus info
    # it's treated as markers/snp info and filtered the same way as the approach by SNP.
    # but it's really by haplotype
  }

  if (maf.approach == "haplotype") {
    if (data.type != "vcf.file" & data.type != "haplo.file") {
      stop("The haplotype approach during MAF filtering is for VCF and
         stacks haplotypes file, only. Use the snp approach for the other file types")
    }
  }

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

  pop.levels <- levels(input$POP_ID)
  pop.labels <- pop.levels

  # Detection if data is summarized for MAF or not -----------------------------
  if ("MAF_LOCAL" %in% colnames(input)) {
    summarize.data <- TRUE
  } else {
    message("\nSummarizing the data by populations and globally")
    summarize.data <- FALSE
  }

  # Summarize data for MAF if required -----------------------------------------
  if (!summarize.data) {
    message("Calculating global and local MAF on large data set may take some time...")

    if (tibble::has_name(input, "GT_VCF")) {
      maf.local <- input %>%
        dplyr::filter(GT_VCF != "./.") %>%
        dplyr::group_by(MARKERS, POP_ID, REF, ALT) %>%
        dplyr::summarise(
          N = as.numeric(n()),
          PQ = as.numeric(length(GT_VCF[GT_VCF == "1/0" | GT_VCF == "0/1"])),
          QQ = as.numeric(length(GT_VCF[GT_VCF == "1/1"]))
        ) %>%
        dplyr::mutate(MAF_LOCAL = ((QQ * 2) + PQ) / (2 * N))

      maf.global <- maf.local %>%
        dplyr::group_by(MARKERS) %>%
        dplyr::summarise_at(.tbl = ., .vars = c("N", "PQ", "QQ"), .funs = sum) %>%
        dplyr::mutate(MAF_GLOBAL = ((QQ * 2) + PQ) / (2 * N)) %>%
        dplyr::select(MARKERS, MAF_GLOBAL)


      maf.data <- maf.global %>%
        dplyr::left_join(maf.local, by = c("MARKERS")) %>%
        dplyr::select(MARKERS, POP_ID, MAF_LOCAL, MAF_GLOBAL)

      maf.local <- maf.global <- NULL
    } else {# not vcf file
      # We split the alleles here to prep for MAF
      maf.data <- input %>%
        dplyr::filter(GT != "000000") %>%
        dplyr::select(MARKERS,POP_ID, INDIVIDUALS, GT) %>%
        dplyr::mutate(
          A1 = stringi::stri_sub(GT, 1, 3),
          A2 = stringi::stri_sub(GT, 4,6)
        ) %>%
        dplyr::select(MARKERS, POP_ID, INDIVIDUALS, A1, A2) %>%
        tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS, POP_ID)) %>%
        dplyr::group_by(MARKERS, GT, POP_ID) %>%
        dplyr::tally(.) %>%
        dplyr::ungroup() %>%
        tidyr::complete(data = ., POP_ID, tidyr::nesting(MARKERS, GT), fill = list(n = 0)) %>%
        dplyr::rename(n.al.pop = n) %>%
        dplyr::arrange(MARKERS, GT) %>%
        dplyr::group_by(MARKERS, GT) %>%
        dplyr::mutate(n.al.tot = sum(n.al.pop)) %>%
        dplyr::group_by(MARKERS) %>%
        dplyr::mutate(MAF_GLOBAL = min(n.al.tot)/sum(n.al.pop)) %>%
        dplyr::group_by(MARKERS, POP_ID) %>%
        dplyr::mutate(MAF_LOCAL = n.al.pop/sum(n.al.pop)) %>%
        dplyr::arrange(MARKERS, POP_ID, GT) %>%
        dplyr::group_by(MARKERS, POP_ID) %>%
        dplyr::filter(n.al.pop == min(n.al.pop)) %>%
        dplyr::distinct(MARKERS, POP_ID, .keep_all = TRUE) %>%
        dplyr::select(MARKERS, POP_ID, MAF_LOCAL, MAF_GLOBAL)
    }# end maf calculations
  } # end summarize data

  # Step 1. Global MAF: Inspecting the MAF globally----------------------------
  if (interactive.filter) {
    message("\nStep 1. Global MAF: Inspecting the MAF globally\n")
    # plot_1: Violin plot global MAF individuals and pop
    message("Show the violin plot for the global MAF (y/n)): ")
    violinplot <- as.character(readLines(n = 1))
    if (violinplot == "y") {

      message("Generating violin plot may take some time...")
      # plot
      OVERALL <- NULL
      global.data <- dplyr::ungroup(maf.data) %>%
        dplyr::distinct(MARKERS, MAF_GLOBAL) %>%
        dplyr::mutate(OVERALL = rep("overall", n()))

      maf.global.summary <- dplyr::ungroup(global.data) %>%
        dplyr::summarise(
          MEAN = mean(MAF_GLOBAL, na.rm = TRUE),
          MEDIAN = stats::median(MAF_GLOBAL, na.rm = TRUE),
          RANGE = stringi::stri_join(round(min(MAF_GLOBAL, na.rm = TRUE), 4), " - ", round(max(MAF_GLOBAL, na.rm = TRUE), 4))
        )

      violinplot.maf.global <- ggplot(data = global.data, aes(x = OVERALL, y = MAF_GLOBAL, na.rm = TRUE)) +
        geom_violin(trim = TRUE) +
        geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA) +
        stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white") +
        # scale_y_continuous(name = "Global MAF", breaks = c(0, 0.01, 0.02, 0.))
        labs(x = "Overall") +
        labs(y = "Global MAF") +
        theme(
          legend.position = "none",
          axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
          axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
          axis.text.x = element_blank(),
          # axis.text.x = element_text(size = 8, family = "Helvetica"),
          legend.title = element_text(size = 10, family = "Helvetica", face = "bold"),
          legend.text = element_text(size = 10, family = "Helvetica", face = "bold"),
          strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold")
        )
      print(violinplot.maf.global)
      # save
      ggsave(stringi::stri_join(path.folder, "/maf.global.pdf"), width = 10, height = 10, dpi = 600, units = "cm", useDingbats = F)
      ggsave(stringi::stri_join(path.folder, "/maf.global.png"), width = 10, height = 10, dpi = 300, units = "cm")
      message(stringi::stri_join("2 versions (pdf and png) of the violin plot for the global MAF were saved in this directory: \n", path.folder))

      global.data <- NULL # unused object

      readr::write_tsv(x = maf.global.summary, path = paste0(path.folder, "/maf.global.summary.tsv"))
      message(stringi::stri_join("The global MAF mean: ", round(maf.global.summary$MEAN, 4)))
      message(stringi::stri_join("The global MAF median: ", round(maf.global.summary$MEDIAN, 4)))
      message(stringi::stri_join("The global MAF range: ", maf.global.summary$RANGE))
      message(stringi::stri_join("maf.global.summary.tsv was saved in this directory: \n", path.folder))
    }
    violinplot.maf.global <- "not selected"
    maf.global.summary <- "not selected"
  } # end global maf

  # Step 2. Local MAF: Inspecting the MAF at the populations level--------------
  # plot_2: Violin plot local MAF
  if (interactive.filter) {
    message("\nStep 2. Local MAF: Inspecting the MAF at the population level\n")
    # plot_2: Violin plot local MAF
    message("Show the violin plot for the local MAF (y/n)): ")
    violinplot <- as.character(readLines(n = 1))
    if (violinplot == "y") {
      pop.number <- length(levels(maf.data$POP_ID))

      message("Generating violin plot may take some time...")
      # plot
      # maf.local.summary <- maf.data %>%
      #   ungroup %>%
      #   summarise(
      #     MEAN = mean(MAF_LOCAL, na.rm = TRUE),
      #     MEDIAN = stats::median(MAF_LOCAL, na.rm = TRUE),
      #     RANGE = stringi::stri_join(round(min(MAF_LOCAL, na.rm = TRUE), 4), " - ", round(max(MAF_GLOBAL, na.rm = TRUE), 4))
      #   )

      violinplot.maf.local <- ggplot(data = maf.data, aes(x = POP_ID, y = MAF_LOCAL, na.rm = TRUE)) +
        geom_violin(trim = TRUE) +
        geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA) +
        stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white") +
        labs(x = "Populations/Groupings") +
        labs(y = "Local/populations MAF") +
        theme(
          legend.position = "none",
          axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
          axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
          # axis.text.x = element_blank(),
          axis.text.x = element_text(size = 8, family = "Helvetica"),
          legend.title = element_text(size = 10, family = "Helvetica", face = "bold"),
          legend.text = element_text(size = 10, family = "Helvetica", face = "bold"),
          strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold")
        )
      print(violinplot.maf.local)
      # save
      ggsave(stringi::stri_join(path.folder, "/maf.local.violinplot.pdf"), width = pop.number, height = 10, dpi = 600, units = "cm", useDingbats = F)
      ggsave(stringi::stri_join(path.folder, "/maf.local.violinplot.png"), width = pop.number, height = 10, dpi = 300, units = "cm")
      message(stringi::stri_join("2 versions (pdf and png) of the violin plot for the global MAF were saved in this directory: \n", path.folder))
    }
    violinplot.maf.local <- "not selected"
  }

  # plot_3: Distribution of local MAF
  if (interactive.filter) {
    message("Show site frequency spectrum (y/n): ")
    spectrum <- as.character(readLines(n = 1))
    if (spectrum == "y") {
      pop.number <- length(levels(maf.data$POP_ID))
      plot.distribution.maf.local <- ggplot(data = maf.data, aes(x = MAF_LOCAL, na.rm = FALSE)) +
        geom_line(aes(y = ..scaled.., color = POP_ID), stat = "density", adjust = 1) + # pop colored
        #   scale_colour_manual(name ="Sampling sites", values = colour_palette_sites.pink) +
        scale_x_continuous(breaks = seq(0,1, by = 0.1)) +
        # labels = c("0", "0.02", "0.05", "0.10", "0.20", "0.50", "1.00"),
        # limits = c(0,1)) +
        labs(x = "Minor Allele Frequency (MAF)") +
        labs(y = "Density of MARKERS (scaled)") +
        expand_limits(y = 0) +
        theme(
          axis.title.x = element_text(size = 12, family = "Helvetica", face = "bold"),
          axis.title.y = element_text(size = 12, family = "Helvetica", face = "bold"),
          legend.position = "none",
          axis.text.x = element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5),
          # legend.title = element_text(size = 12, family = "Helvetica", face = "bold"),
          # legend.text = element_text(size = 12, family = "Helvetica", face = "bold"),
          strip.text.y = element_text(angle = 0, size = 12, family = "Helvetica", face = "bold"),
          strip.text.x = element_text(size = 12, family = "Helvetica", face = "bold")
        ) +
        facet_grid(~POP_ID)
      print(plot.distribution.maf.local)
      ggsave(stringi::stri_join(path.folder, "/maf.local.spectrum.pdf"), width = pop.number * 5, height = 15, dpi = 600, units = "cm", useDingbats = F)
      ggsave(stringi::stri_join(path.folder, "/maf.local.spectrum.png"), width = pop.number * 5, height = 10, dpi = 300, units = "cm")
      message(stringi::stri_join("2 versions (pdf and png) of the local maf spectrum plot were saved in this directory: \n", path.folder))
    }
    plot.distribution.maf.local <- "not selected"
  }

  # Helper table for global and local MAF --------------------------------------
  # number of individuals / pop
  maf.helper.table <- input %>%
    dplyr::group_by(POP_ID, INDIVIDUALS) %>%
    dplyr::distinct(POP_ID, INDIVIDUALS, .keep_all = TRUE) %>%
    dplyr::ungroup(.) %>%
    dplyr::group_by(POP_ID) %>%
    dplyr::tally(.)

  # Global
  n.ind <- dplyr::n_distinct(input$INDIVIDUALS)
  TOTAL <- tibble::data_frame(POP_ID = "TOTAL/GLOBAL", n = n.ind)

  # bind pop and global + MAF local and global for different ALT allele
  maf.helper.table <- suppressWarnings(
    dplyr::bind_rows(maf.helper.table, TOTAL) %>%
      dplyr::mutate(
        ALT_1 = 1/(2*n), ALT_2 = 2/(2*n), ALT_3 = 3/(2*n), ALT_4 = 4/(2*n),
        ALT_5 = 5/(2*n), ALT_6 = 6/(2*n), ALT_7 = 7/(2*n), ALT_8 = 8/(2*n),
        ALT_9 = 9/(2*n), ALT_10 = 10/(2*n), ALT_11 = 11/(2*n), ALT_12 = 12/(2*n),
        ALT_13 = 13/(2*n), ALT_14 = 14/(2*n), ALT_15 = 15/(2*n), ALT_16 = 16/(2*n),
        ALT_17 = 17/(2*n), ALT_18 = 18/(2*n), ALT_19 = 19/(2*n), ALT_20 = 20/(2*n)
      )
  )
    # remove unused objects
  TOTAL <- n.ind <- NULL

  readr::write_tsv(x = maf.helper.table, path = paste0(path.folder, "/maf.helper.table.tsv"))

  message(stringi::stri_join("\nTo visualize the local and global MAF, a table
(maf.helper.table.tsv) was written in this directory: \n", path.folder))
  if (interactive.filter) {
    message("\nFirst and second variable columns: POP_ID and sample size.
The last row is the TOTAL/GLOBAL observation.
The remaining columns are the variable corresponding to the number of ALT allele (range 1 to 20).
The observations in the ALT allele variable columns are the local (for the pop)
and global (last row) MAF of your dataset.\n
e.g. ALT_3 could be 3 heterozygote individuals with the ALT allele or
1 homozygote individuals for the ALT allele and 1 heterozygote individual. And so on...\n")
  }
  # Step 3. Thresholds selection -----------------------------------------------

  # maf.approach
  if (interactive.filter) {
    message("\nStep 3. Filtering markers based on the different MAF arguments\n")
    message("The maf.approach:\n
maf.approach = \"haplotype\" : looks at the minimum MAF found on the read/haplotype.
This will discard all the markers/snp on that read based on the thresholds chosen.
This method is only available for VCF and haplotype files or tidy data frame from
those file type.\n
maf.approach = \"SNP\" : treats all the SNP on the same haplotype/read as independent.")
    message("Choose the maf.approach (SNP/haplotype):")
    maf.approach <- as.character(readLines(n = 1))
    if (!maf.approach %in% c("SNP", "haplotype")) stop("maf.approach: SNP or haplotype")
  }

  # maf.thresholds
  if (interactive.filter) {
    message("Choose the maf local threshold (usually a value between 0 and 0.3):")
    maf.local.threshold <- as.character(readLines(n = 1))
    message("Choose the maf global threshold (usually a value between 0 and 0.3):")
    maf.global.threshold <- as.character(readLines(n = 1))
  }

  # maf.operator
  if (interactive.filter) {
    message("The maf.operator:
Option 1: AND = the local \"AND\" the global MAF threshold are required to pass.
Option 2: OR = EITHER the local \"OR\" the global MAF threshold are required to pass.")
    message("Choose the maf.operator (AND/OR):")
    maf.operator <- as.character(readLines(n = 1))
    if (!maf.operator %in% c("OR", "AND")) stop("maf.operator: either OR/AND")
  }

  # maf.pop.num.threshold
  if (interactive.filter) {
    message("The maf.pop.num.threshold:\n
How many populations are required to pass the thresholds to keep the locus?\n
Example: if you have 10 populations and choose maf.pop.num.threshold = 3,
3 populations out of 10 are required to have LOCAL and/or GLOBAL MAF higher than the
thresholds entered.")
    message("Choose the value for the maf.pop.num.threshold:")
    maf.pop.num.threshold <- as.numeric(readLines(n = 1))
  }

  if (!interactive.filter) {
    maf.local.threshold <- maf.thresholds[1]
    maf.global.threshold <- maf.thresholds[2]
  }


  # Update the maf.data with pass or not filter based on threshold -------------
  maf.data.thresholds <- maf.data %>%
    dplyr::mutate(
      OR = ifelse((MAF_LOCAL >= maf.local.threshold | MAF_GLOBAL >= maf.global.threshold), "pass", "pruned"),
      AND = ifelse((MAF_LOCAL >= maf.local.threshold & MAF_GLOBAL >= maf.global.threshold), "pass", "pruned")
    ) %>%
    dplyr::group_by(MARKERS) %>%
    dplyr::mutate(
  OR_POP_THRESHOLD = ifelse(length(POP_ID[OR == "pass"]) >= maf.pop.num.threshold, "pass", "pruned"),
  AND_POP_THRESHOLD  = ifelse(length(POP_ID[AND == "pass"]) >= maf.pop.num.threshold, "pass", "pruned")
  )

readr::write_tsv(
  x = maf.data.thresholds,
  path = paste0(path.folder, "/maf.data.tsv"),
  col_names = TRUE,
  append = FALSE
)
message(stringi::stri_join("\nThe MAF summary statistics (maf.data.tsv), written in this directory: \n", path.folder))

# Filtering ------------------------------------------------------------------
if (maf.approach == "haplotype") {
  filter <- tidyr::separate(
    data = maf.data,
    col = MARKERS,
    into = c("CHROM", "LOCUS", "POS"),
    sep = "__",
    remove = FALSE,
    extra = "warn"
  )

  if (maf.operator == "OR") {
    filter <- filter %>%
      dplyr::group_by(LOCUS, POP_ID) %>%
      dplyr::summarise(
        MAF_GLOBAL = mean(MAF_GLOBAL, na.rm = TRUE),
        MAF_LOCAL = min(MAF_LOCAL, na.rm = TRUE)
      ) %>%
      dplyr::filter(MAF_LOCAL >= maf.local.threshold | MAF_GLOBAL >= maf.global.threshold) %>%
      dplyr::group_by(LOCUS) %>%
      dplyr::tally(.) %>%
      dplyr::filter(n >= maf.pop.num.threshold) %>%
      dplyr::select(LOCUS) %>%
      dplyr::left_join(input, by = "LOCUS") %>%
      dplyr::arrange(LOCUS, POP_ID)
  } else {# AND operator between local and global maf
    filter <- filter %>%
      dplyr::group_by(LOCUS, POP_ID) %>%
      dplyr::summarise(
        MAF_GLOBAL = mean(MAF_GLOBAL, na.rm = TRUE),
        MAF_LOCAL = min(MAF_LOCAL, na.rm = TRUE)
      ) %>%
      dplyr::filter(MAF_LOCAL >= maf.local.threshold & MAF_GLOBAL >= maf.global.threshold) %>%
      dplyr::group_by(LOCUS) %>%
      dplyr::tally(.) %>%
      dplyr::filter(n >= maf.pop.num.threshold) %>%
      dplyr::select(LOCUS) %>%
      dplyr::left_join(input, by = "LOCUS") %>%
      dplyr::arrange(LOCUS, POP_ID)
  }
  # filter <- filter %>% select(-c(CHROM, LOCUS, POS))
} # end maf haplotype approach

if (maf.approach == "SNP") { # SNP approach
  if (maf.operator == "OR") {
    filter <- maf.data %>%
      dplyr::group_by(MARKERS, POP_ID) %>%
      dplyr::summarise(
        MAF_GLOBAL = mean(MAF_GLOBAL, na.rm = TRUE),
        MAF_LOCAL = min(MAF_LOCAL, na.rm = TRUE)
      ) %>%
      dplyr::filter(MAF_LOCAL >= maf.local.threshold | MAF_GLOBAL >= maf.global.threshold) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::tally(.) %>%
      dplyr::filter(n >= maf.pop.num.threshold) %>%
      dplyr::select(MARKERS) %>%
      dplyr::left_join(input, by = "MARKERS") %>%
      dplyr::arrange(MARKERS, POP_ID)
  } else {# AND operator between local and global maf
    filter <- maf.data %>%
      dplyr::group_by(MARKERS, POP_ID) %>%
      dplyr::summarise(
        MAF_GLOBAL = mean(MAF_GLOBAL, na.rm = TRUE),
        MAF_LOCAL = min(MAF_LOCAL, na.rm = TRUE)
      ) %>%
      dplyr::filter(MAF_LOCAL >= maf.local.threshold & MAF_GLOBAL >= maf.global.threshold) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::tally(.) %>%
      dplyr::filter(n >= maf.pop.num.threshold) %>%
      dplyr::select(MARKERS) %>%
      dplyr::left_join(input, by = "MARKERS") %>%
      dplyr::arrange(MARKERS, POP_ID)
  }
} # end maf snp approach

# unused object
maf.data <- NULL

# Update filters.parameters SNP ----------------------------------------------
if (tibble::has_name(input, "LOCUS") & maf.approach == "haplotype") {
  snp.before <- as.integer(dplyr::n_distinct(input$MARKERS))
  snp.after <- as.integer(dplyr::n_distinct(filter$MARKERS))
  snp.blacklist <- as.integer(dplyr::n_distinct(input$MARKERS) - dplyr::n_distinct(filter$MARKERS))
  locus.before <- as.integer(dplyr::n_distinct(input$LOCUS))
  locus.after <- as.integer(dplyr::n_distinct(filter$LOCUS))
  locus.blacklist <- as.integer(locus.before - locus.after)
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
  FILTERS = c("Minor Allele Frequency", rep(as.character(""), 4)),
  PARAMETERS = c("maf.approach", "maf.local.threshold", "maf.global.threshold", "maf.operator", "maf.pop.num.threshold"),
  VALUES = c(maf.approach, paste(">=", maf.local.threshold), paste(">=", maf.global.threshold), maf.operator, paste(">=", maf.pop.num.threshold)),
  BEFORE = c("", "", "", "", markers.before),
  AFTER = c("", "", "", "", markers.after),
  BLACKLIST = c("", "", "", "", markers.blacklist),
  UNITS = c("", "", "", "", "SNP/LOCUS"),
  COMMENTS = c("", "", "", "", "")
)
readr::write_tsv(x = filters.parameters, path = "filters_parameters.tsv", append = TRUE, col_names = FALSE)

# saving tidy data
if (!is.null(filename)) {
  message("\nWriting the filtered tidy data set in your working directory...")
  # if (!is.null(save.feather)) {
  # feather::write_feather(filter, stringi::stri_replace_all_fixed(filename, pattern = ".tsv", replacement = "_feather.tsv", vectorize_all = TRUE))
  # } else {
  readr::write_tsv(filter, filename, append = FALSE, col_names = TRUE)
  # }
}

# saving whitelist
message("\nWriting the whitelist of markers in your working directory\nwhitelist.markers.maf.tsv")

if (tibble::has_name(input, "CHROM")) {
  whitelist.markers <- dplyr::ungroup(filter) %>%
    dplyr::distinct(CHROM, LOCUS, POS)
} else {
  whitelist.markers <- dplyr::ungroup(filter) %>%
    dplyr::distinct(MARKERS)
}
readr::write_tsv(whitelist.markers, paste0(path.folder, "/whitelist.markers.maf.tsv"), append = FALSE, col_names = TRUE)


# saving blacklist
message("\nWriting the blacklist of markers in your working directory\nblacklist.markers.maf.tsv")
if (tibble::has_name(input, "CHROM")) {
  blacklist.markers <- dplyr::ungroup(input) %>%
    dplyr::distinct(CHROM, LOCUS, POS) %>%
    dplyr::anti_join(whitelist.markers, by = c("CHROM", "LOCUS", "POS"))
} else {
  blacklist.markers <- dplyr::ungroup(input) %>%
    dplyr::distinct(MARKERS) %>%
    dplyr::anti_join(whitelist.markers, by = "MARKERS")
}
readr::write_tsv(blacklist.markers, paste0(path.folder, "/blacklist.markers.maf.tsv"), append = FALSE, col_names = TRUE)

# results --------------------------------------------------------------------
cat("############################### RESULTS ###############################\n")
message(stringi::stri_join("maf.approach: ", maf.approach))
message(stringi::stri_join("maf.thresholds: ", "local = ", maf.local.threshold, ", global = ", maf.global.threshold))
message(stringi::stri_join("maf.operator: ", maf.operator))
message(stringi::stri_join("maf.pop.num.threshold: ", maf.pop.num.threshold))
if (tibble::has_name(input, "LOCUS") & maf.approach == "haplotype") {
  message(stringi::stri_join("The number of markers removed by the MAF filter:\nSNP: ", snp.blacklist, "\nLOCUS: ", locus.blacklist))
  message("The number of markers before -> after the MAF filter")
  message(stringi::stri_join("SNP: ", snp.before, " -> ", snp.after))
  message(stringi::stri_join("LOCUS: ", locus.before, " -> ", locus.after))
} else {
  message(stringi::stri_join("The number of markers removed by the MAF filter: ", snp.blacklist))
  message("The number of markers before -> after the MAF filter")
  message(stringi::stri_join("SNP: ", snp.before, " -> ", snp.after))
}
if (!interactive.filter) {
  timing <- proc.time() - timing
  message(stringi::stri_join("Computation time: ", round(timing[[3]]), " sec"))
}
cat("############################## completed ##############################\n")
res <- list()
res$tidy.filtered.maf <- filter
res$whitelist.markers <- whitelist.markers
res$blacklist.markers <- blacklist.markers
res$maf.data.thresholds <- maf.data.thresholds
res$maf.helper.table <- maf.helper.table
res$strata <- strata
res$filters.parameters <- filters.parameters
if (interactive.filter) {
  res$violinplot.maf.global <- violinplot.maf.global
  res$violinplot.maf.local <- violinplot.maf.local
  res$plot.distribution.maf.local <- plot.distribution.maf.local
  res$maf.global.summary <- maf.global.summary
}
return(res)
}
