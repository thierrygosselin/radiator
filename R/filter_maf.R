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
  cat("######################## radiator::filter_maf #########################\n")
  cat("#######################################################################\n")
  opt.change <- getOption("width")
  options(width = 70)
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
  # Date and time
  file.date <- stringi::stri_replace_all_fixed(
    Sys.time(),
    pattern = " EDT", replacement = "") %>%
    stringi::stri_replace_all_fixed(
      str = .,
      pattern = c("-", " ", ":"), replacement = c("", "@", ""),
      vectorize_all = FALSE) %>%
    stringi::stri_sub(str = ., from = 1, to = 13)
  folder.extension <- stringi::stri_join("filter_maf_", file.date, sep = "")
  path.folder <- stringi::stri_join(getwd(),"/", folder.extension, sep = "")
  dir.create(file.path(path.folder))
  message("\nFolder created: ", folder.extension)
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

  if (data.type == "haplo.file") {
    message("With stacks haplotype file the maf.approach is automatically set to: haplotype")
    maf.approach <- "SNP"
    # confusing, but because the haplotpe file doesn't have snp info, only locus info
    # it's treated as markers/snp info and filtered the same way as the approach by SNP.
    # but it's really by haplotype
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
    filename = NULL,
    parallel.core = parallel.core
  )

  # create a strata.df
  strata.df <- input %>%
    dplyr::select(INDIVIDUALS, POP_ID) %>%
    dplyr::distinct(INDIVIDUALS, .keep_all = TRUE)

  pop.levels <- levels(input$POP_ID)
  pop.labels <- pop.levels

  if (maf.approach == "haplotype") {
    # if (data.type != "vcf.file" & data.type != "haplo.file") {
    if (!tibble::has_name(input, "LOCUS") || !tibble::has_name(input, "POS")) {
      stop("The haplotype approach during MAF filtering is for VCF and
         stacks haplotypes file, only. Use the snp approach for the other file types")
    }
  }

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

    biallelic <- radiator::detect_biallelic_markers(input, verbose = TRUE)

    if (tibble::has_name(input, "GT_VCF") && biallelic) {
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
    } else  {# not vcf file
      if (!tibble::has_name(input, "GT_VCF_NUC")) {
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
          dplyr::ungroup(.) %>%
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
      } else {
        # Experimental, for haplotypes file
        ref.allele <- dplyr::distinct(input, MARKERS, REF)
        maf.prep <- input %>%
          dplyr::filter(GT_VCF_NUC != "./.") %>%
          dplyr::select(MARKERS, POP_ID, INDIVIDUALS, GT_VCF_NUC) %>%
          separate_gt(
            x = ., sep = "/",
            gt = "GT_VCF_NUC",
            exclude = c("MARKERS", "POP_ID", "INDIVIDUALS"),
            cpu.rounds = 10) %>%
          dplyr::select(-ALLELE_GROUP) %>%
          dplyr::rename(HAPLOTYPES = ALLELES) %>%
          dplyr::group_by(MARKERS, HAPLOTYPES, POP_ID) %>%
          dplyr::tally(.) %>%
          dplyr::ungroup(.) %>%
          tidyr::complete(data = ., POP_ID, tidyr::nesting(MARKERS, HAPLOTYPES), fill = list(n = 0)) %>%
          dplyr::rename(n.al.pop = n) %>%
          dplyr::arrange(MARKERS, HAPLOTYPES) %>%
          dplyr::group_by(MARKERS, HAPLOTYPES) %>%
          dplyr::mutate(n.al.tot = sum(n.al.pop)) %>%
          dplyr::left_join(ref.allele, by = "MARKERS") %>%
          dplyr::mutate(REF = dplyr::if_else(HAPLOTYPES == REF, "REF", "ALT")) %>%
          dplyr::group_by(MARKERS) %>%
          dplyr::mutate(N = sum(n.al.pop)) %>%
          dplyr::group_by(MARKERS, POP_ID) %>%
          dplyr::mutate(N_POP = sum(n.al.pop)) %>%
          dplyr::filter(REF == "ALT") %>%
          dplyr::group_by(MARKERS, HAPLOTYPES) %>%
          dplyr::mutate(MAF_GLOBAL = n.al.tot/N) %>%
          dplyr::group_by(MARKERS, HAPLOTYPES, POP_ID) %>%
          dplyr::mutate(MAF_LOCAL = n.al.pop/N_POP) %>%
          dplyr::arrange(MARKERS, POP_ID, HAPLOTYPES) %>%
          dplyr::ungroup(.)

        maf.data <- dplyr::select(
          maf.prep,
          MARKERS, POP_ID, HAPLOTYPES, MAF_LOCAL, MAF_GLOBAL)
      }
    }# end maf calculations
  } # end summarize data

  # Step 1. Global MAF: Inspecting the MAF globally----------------------------
  if (interactive.filter) {
    message("\nStep 1. Global MAF visualization\n")
    # plot_1: Violin plot global MAF individuals and pop
  }
  # violinplot <- as.character(readLines(n = 1))
  # if (violinplot == "y") {

  # message("Generating violin plot may take some time...")
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

  violinplot.maf.global <- ggplot2::ggplot(data = global.data, ggplot2::aes(x = OVERALL, y = MAF_GLOBAL, na.rm = TRUE)) +
    ggplot2::geom_violin(trim = TRUE) +
    ggplot2::geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA) +
    ggplot2::stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white") +
    # scale_y_continuous(name = "Global MAF", breaks = c(0, 0.01, 0.02, 0.))
    ggplot2::labs(x = "Overall") +
    ggplot2::labs(y = "Global MAF") +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.text.x = ggplot2::element_blank(),
      # axis.text.x = element_text(size = 8, family = "Helvetica"),
      legend.title = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.text = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
    )
  if (interactive.filter) {
    message("    Showing the violin plot for the global MAF")
    print(violinplot.maf.global)
  }
  # save
  ggplot2::ggsave(filename = stringi::stri_join(path.folder, "/maf.global.pdf"), plot = violinplot.maf.global, width = 10, height = 10, dpi = 600, units = "cm", useDingbats = FALSE)
  ggplot2::ggsave(filename = stringi::stri_join(path.folder, "/maf.global.png"), plot = violinplot.maf.global, width = 10, height = 10, dpi = 300, units = "cm")
  message("2 versions (pdf and png) of the violin plot for the global MAF were written in the folder")

  global.data <- NULL # unused object

  readr::write_tsv(x = maf.global.summary, path = paste0(path.folder, "/maf.global.summary.tsv"))
  message("maf.global.summary.tsv was saved in the folder")

  if (interactive.filter) {
    message("The global MAF mean: ", round(maf.global.summary$MEAN, 4))
    message("The global MAF median: ", round(maf.global.summary$MEDIAN, 4))
    message("The global MAF range: ", maf.global.summary$RANGE)
  }

  # Step 2. Local MAF: Inspecting the MAF at the populations level--------------
  # plot_2: Violin plot local MAF
  if (interactive.filter) {
    message("\nStep 2. Local MAF visualization (population level)\n")
  }
  # plot_2: Violin plot local MAF
  # violinplot <- as.character(readLines(n = 1))
  # if (violinplot == "y") {
  pop.number <- length(levels(maf.data$POP_ID))

  # message("Generating violin plot may take some time...")
  # plot
  # maf.local.summary <- maf.data %>%
  #   ungroup %>%
  #   summarise(
  #     MEAN = mean(MAF_LOCAL, na.rm = TRUE),
  #     MEDIAN = stats::median(MAF_LOCAL, na.rm = TRUE),
  #     RANGE = stringi::stri_join(round(min(MAF_LOCAL, na.rm = TRUE), 4), " - ", round(max(MAF_GLOBAL, na.rm = TRUE), 4))
  #   )

  violinplot.maf.local <- ggplot2::ggplot(data = maf.data, ggplot2::aes(x = POP_ID, y = MAF_LOCAL, na.rm = TRUE)) +
    ggplot2::geom_violin(trim = TRUE) +
    ggplot2::geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA) +
    ggplot2::stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white") +
    ggplot2::labs(x = "Populations/Groupings") +
    ggplot2::labs(y = "Local/populations MAF") +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      # axis.text.x = element_blank(),
      axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica"),
      legend.title = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.text = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
    )
  if (interactive.filter) {
    message("Showing the violin plot for the local MAF")
    print(violinplot.maf.local)
  }

  # save
  ggplot2::ggsave(filename = stringi::stri_join(path.folder, "/maf.local.violinplot.pdf"), plot = violinplot.maf.local, width = pop.number, height = 10, dpi = 600, units = "cm", useDingbats = FALSE)
  ggplot2::ggsave(filename = stringi::stri_join(path.folder, "/maf.local.violinplot.png"), plot = violinplot.maf.local, width = pop.number, height = 10, dpi = 300, units = "cm")
  message("2 versions (pdf and png) of the violin plot for the global MAF were saved in the folder")

  # plot_3: Distribution of local MAF

  # spectrum <- as.character(readLines(n = 1))
  # if (spectrum == "y") {
  pop.number <- length(levels(maf.data$POP_ID))
  plot.distribution.maf.local <- ggplot2::ggplot(data = maf.data, ggplot2::aes(x = MAF_LOCAL, na.rm = FALSE)) +
    ggplot2::geom_line(ggplot2::aes(y = ..scaled.., color = POP_ID), stat = "density", adjust = 1) + # pop colored
    #   scale_colour_manual(name ="Sampling sites", values = colour_palette_sites.pink) +
    ggplot2::scale_x_continuous(breaks = seq(0,1, by = 0.1)) +
    # labels = c("0", "0.02", "0.05", "0.10", "0.20", "0.50", "1.00"),
    # limits = c(0,1)) +
    ggplot2::labs(x = "Minor Allele Frequency (MAF)") +
    ggplot2::labs(y = "Density of MARKERS (scaled)") +
    ggplot2::expand_limits(y = 0) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
      legend.position = "none",
      axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5),
      # legend.title = element_text(size = 12, family = "Helvetica", face = "bold"),
      # legend.text = element_text(size = 12, family = "Helvetica", face = "bold"),
      strip.text.y = ggplot2::element_text(angle = 0, size = 12, family = "Helvetica", face = "bold"),
      strip.text.x = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold")
    ) +
    ggplot2::facet_grid(~POP_ID)
  if (interactive.filter) {
    message("Showing site frequency spectrum")
    print(plot.distribution.maf.local)
  }
  ggplot2::ggsave(filename = stringi::stri_join(path.folder, "/maf.local.spectrum.pdf"), plot = plot.distribution.maf.local, width = pop.number * 5, height = 15, dpi = 600, units = "cm", useDingbats = FALSE)
  ggplot2::ggsave(filename = stringi::stri_join(path.folder, "/maf.local.spectrum.png"), plot.distribution.maf.local, width = pop.number * 5, height = 10, dpi = 300, units = "cm")
  message("2 versions (pdf and png) of the local maf spectrum plot were saved in the folder")

  # Helper table for global and local MAF --------------------------------------
  # number of individuals / pop
  message("Generating maf.helper.table...")
  maf.helper.table <- input %>%
    dplyr::group_by(POP_ID, INDIVIDUALS) %>%
    dplyr::distinct(POP_ID, INDIVIDUALS, .keep_all = TRUE) %>%
    dplyr::ungroup(.) %>%
    dplyr::group_by(POP_ID) %>%
    dplyr::tally(.) %>%
    tibble::add_row(.data = ., POP_ID = "TOTAL/GLOBAL",
                    n = dplyr::n_distinct(input$INDIVIDUALS)) %>%
    dplyr::mutate(
      POP_ID = as.character(POP_ID),
      ALT_1 = 1/(2*n), ALT_2 = 2/(2*n), ALT_3 = 3/(2*n), ALT_4 = 4/(2*n),
      ALT_5 = 5/(2*n), ALT_6 = 6/(2*n), ALT_7 = 7/(2*n), ALT_8 = 8/(2*n),
      ALT_9 = 9/(2*n), ALT_10 = 10/(2*n), ALT_11 = 11/(2*n), ALT_12 = 12/(2*n),
      ALT_13 = 13/(2*n), ALT_14 = 14/(2*n), ALT_15 = 15/(2*n), ALT_16 = 16/(2*n),
      ALT_17 = 17/(2*n), ALT_18 = 18/(2*n), ALT_19 = 19/(2*n), ALT_20 = 20/(2*n)
  )

  new <- dplyr::select(maf.helper.table, -n) %>%
    dplyr::filter(POP_ID != "TOTAL/GLOBAL") %>%
    dplyr::summarise_if(.tbl = ., .predicate = is.numeric, .funs = min ) %>%
    dplyr::mutate(POP_ID = "MIN_LOCAL", n = NA) %>%
    dplyr::select(POP_ID, n, dplyr::everything())

  sort.pop <- c(levels(input$POP_ID), "MIN_LOCAL", "TOTAL/GLOBAL")

  maf.helper.table <- dplyr::bind_rows(maf.helper.table, new) %>%
    dplyr::mutate(POP_ID = factor(POP_ID, levels = sort.pop, ordered = TRUE)) %>%
    dplyr::arrange(POP_ID)


  # remove unused objects
  new <- NULL

  readr::write_tsv(
    x = maf.helper.table,
    path = stringi::stri_join(path.folder, "/maf.helper.table.tsv"))

  message("\nA table of local and global MAF (maf.helper.table.tsv) was written in the directory")
  if (interactive.filter) {
    message("\nFirst and second variable columns represents POP_ID and sample size (n).
\nThe last 2 rows are the local MAF (suggested based on the lowest pop value) and\nthe TOTAL/GLOBAL observations.
\nColumns starting with ALT are the variable corresponding\nto the number of alternative (ALT) allele (ranging from 1 to 20).
\nThe observations in the ALT allele variable columns are the local (for the pop)
and global (last row) MAF of your dataset.\n
e.g. ALT_3 can potentially represent 3 heterozygote individuals with the ALT allele or
1 homozygote individuals for the ALT allele and 1 heterozygote individual. And so on...\n")
  }

  # Step 3. Thresholds selection -----------------------------------------------

  # maf.approach
  if (interactive.filter) {
    message("\nStep 3. Filtering markers based on the different MAF arguments\n")
    message("The maf.approach:\n
maf.approach = \"haplotype\" : looks at the frequency of all ALTernative allele on the locus.\n
maf.approach = \"SNP\" : SNPs on the same haplotype/read are considered independent.")
    message("Choose the maf.approach (SNP/haplotype):")
    maf.approach <- as.character(readLines(n = 1))
    if (!maf.approach %in% c("SNP", "haplotype")) stop("maf.approach: SNP or haplotype")
  }

  # maf.thresholds
  if (interactive.filter) {
    message("Choose the maf local threshold (usually a value between 0 and 0.3):")
    maf.local.threshold <- as.character(readLines(n = 1))
  }
  if (interactive.filter) {
    message("Choose the maf global threshold (usually a value between 0 and 0.3):")
    maf.global.threshold <- as.character(readLines(n = 1))
  }

  # maf.operator
  if (interactive.filter) {
    message("The maf.operator:
Option 1: AND = to consider both the local \"AND\" the global MAF threshold (more severe).
Option 2: OR = to consider EITHER the local \"OR\" the global MAF threshold (more tolerant).")
    message("Choose the maf.operator (AND/OR):")
    maf.operator <- as.character(readLines(n = 1))
    if (!maf.operator %in% c("OR", "AND")) stop("maf.operator: either OR/AND")
  }

  # maf.pop.num.threshold
  if (interactive.filter) {
    message("Last threshold... maf.pop.num.threshold:\n
How many populations are required to pass all the thresholds and keep the locus/SNP?\n
Example: if you have 10 populations and choose maf.pop.num.threshold = 3,
3 populations out of 10 are required to pass previous thresholds")
    message("Choose the value for the maf.pop.num.threshold:")
    maf.pop.num.threshold <- as.numeric(readLines(n = 1))
  }

  if (!interactive.filter) {
    maf.local.threshold <- maf.thresholds[1]
    maf.global.threshold <- maf.thresholds[2]
  }


  # Update the maf.data with pass or not filter based on threshold -------------
  if (tibble::has_name(maf.data, "HAPLOTYPES")) {
    maf.data.thresholds <- maf.data %>%
      dplyr::mutate(
        OR = dplyr::if_else((MAF_LOCAL >= maf.local.threshold | MAF_GLOBAL >= maf.global.threshold), "pass", "pruned"),
        AND = dplyr::if_else((MAF_LOCAL >= maf.local.threshold & MAF_GLOBAL >= maf.global.threshold), "pass", "pruned")
      ) %>%
      dplyr::group_by(MARKERS, HAPLOTYPES) %>%
      dplyr::mutate(
        OR_POP_THRESHOLD = dplyr::if_else(length(POP_ID[OR == "pass"]) >= maf.pop.num.threshold, "pass", "pruned"),
        AND_POP_THRESHOLD  = dplyr::if_else(length(POP_ID[AND == "pass"]) >= maf.pop.num.threshold, "pass", "pruned")
      )
  } else {
    maf.data.thresholds <- maf.data %>%
      dplyr::mutate(
        OR = dplyr::if_else((MAF_LOCAL >= maf.local.threshold | MAF_GLOBAL >= maf.global.threshold), "pass", "pruned"),
        AND = dplyr::if_else((MAF_LOCAL >= maf.local.threshold & MAF_GLOBAL >= maf.global.threshold), "pass", "pruned")
      ) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(
        OR_POP_THRESHOLD = dplyr::if_else(length(POP_ID[OR == "pass"]) >= maf.pop.num.threshold, "pass", "pruned"),
        AND_POP_THRESHOLD  = dplyr::if_else(length(POP_ID[AND == "pass"]) >= maf.pop.num.threshold, "pass", "pruned")
      )
  }

  readr::write_tsv(
    x = maf.data.thresholds,
    path = stringi::stri_join(path.folder, "/maf.data.tsv"),
    col_names = TRUE,
    append = FALSE
  )
  message("\nThe MAF summary statistics (maf.data.tsv) were written in the folder\n")

  # Filtering ------------------------------------------------------------------
  if (maf.approach == "haplotype") {
    need.split <- dplyr::distinct(maf.data, MARKERS) %>%
      dplyr::sample_frac(tbl = ., size = 0.3) %>%
      stringi::stri_detect_fixed(str = ., pattern = "__") %>%
      unique

    if (need.split) {
      maf.data <- tidyr::separate(
        data = maf.data,
        col = MARKERS,
        into = c("CHROM", "LOCUS", "POS"),
        sep = "__",
        remove = FALSE,
        extra = "warn"
      )

      maf.data.thresholds <- tidyr::separate(
        data = maf.data.thresholds,
        col = MARKERS,
        into = c("CHROM", "LOCUS", "POS"),
        sep = "__",
        remove = FALSE,
        extra = "warn"
      )
    } else {
      if (tibble::has_name(maf.data, "MARKERS") && !tibble::has_name(maf.data, "LOCUS")) {
        maf.data <- dplyr::rename(maf.data, LOCUS = MARKERS)
        maf.data.thresholds <- dplyr::rename(maf.data.thresholds, LOCUS = MARKERS)
      }
    }
    need.split <- NULL
    if (maf.operator == "OR") {
      if (tibble::has_name(maf.data, "HAPLOTYPES")) {
        ref.allele <- dplyr::rename(ref.allele, LOCUS = MARKERS, HAPLOTYPES = REF)

        test <- maf.data.thresholds %>%
          dplyr::select(LOCUS, POP_ID, HAPLOTYPES, OR_POP_THRESHOLD) %>%
          dplyr::filter(OR_POP_THRESHOLD == "pass")

        test <- test %>%
          dplyr::distinct(LOCUS, HAPLOTYPES) %>%
          dplyr::bind_rows(ref.allele)

        polymorphism.haplo <- test %>%
          dplyr::group_by(LOCUS) %>%
          dplyr::summarise(HAPLOTYPES = stringi::stri_join(unique(HAPLOTYPES), collapse = ",")) %>%
          dplyr::mutate(POLYMORPHISM = stringi::stri_count_fixed(HAPLOTYPES, ","))

        blacklist.monomorphic.markers <- polymorphism.haplo %>%
          dplyr::filter(POLYMORPHISM == 0) %>%
          dplyr::distinct(LOCUS)


        haplo.reconstruction <- test %>% dplyr::arrange(LOCUS, HAPLOTYPES) %>%
          dplyr::anti_join(blacklist.monomorphic.markers, by = "LOCUS") %>%
          dplyr::mutate(SNP_N = stringi::stri_length(str = HAPLOTYPES))


        no.reconstruction <- dplyr::filter(haplo.reconstruction, SNP_N == 1)
        # dplyr::n_distinct(no.reconstruction$LOCUS)

        haplo.reconstruction2 <- dplyr::filter(haplo.reconstruction, SNP_N > 1)
        # dplyr::n_distinct(haplo.reconstruction2$LOCUS)






      } else{
        maf.data <- maf.data %>%
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
      }
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
  if (tibble::has_name(input, "LOCUS") && maf.approach == "haplotype") {
    snp.before <- as.integer(dplyr::n_distinct(input$MARKERS))
    snp.after <- as.integer(dplyr::n_distinct(filter$MARKERS))
    snp.blacklist <- snp.before - snp.after
    locus.before <- as.integer(dplyr::n_distinct(input$LOCUS))
    locus.after <- as.integer(dplyr::n_distinct(filter$LOCUS))
    locus.blacklist <- locus.before - locus.after
  } else {
    snp.before <- as.integer(dplyr::n_distinct(input$MARKERS))
    snp.after <- as.integer(dplyr::n_distinct(filter$MARKERS))
    snp.blacklist <- snp.before - snp.after
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

  # saving filtered tidy data --------------------------------------------------
  if (!is.null(filename)) {
    tidy.name <- stringi::stri_join(filename, ".rad")
    message("Writing the filtered tidy data set: ", tidy.name)
    fst::write.fst(x = filter, path = stringi::stri_join(path.folder, "/", tidy.name), compress = 85)
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
    message("Computation time: ", round((proc.time() - timing)[[3]]), " sec")
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
  res$violinplot.maf.global <- violinplot.maf.global
  res$violinplot.maf.local <- violinplot.maf.local
  res$plot.distribution.maf.local <- plot.distribution.maf.local
  res$maf.global.summary <- maf.global.summary
  return(res)
}


