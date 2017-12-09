# Minor Allele Frequency
#' @name filter_maf
#' @title MAF filter
#' @description Minor Allele Frequency filter.

# Most arguments are inherited from tidy_genomic_data
#' @inheritParams tidy_genomic_data

#' @param interactive.filter (optional, logical) Do you want the filtering session to
#' be interactive. With default: \code{interactive.filter == TRUE}, the user is
#' shown figures and helper tables before making decisions for filtering.

#' @rdname filter_maf
#' @export

#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot stat_summary labs theme element_blank element_text geom_jitter scale_colour_manual scale_y_reverse theme_light geom_bar facet_grid geom_histogram aes_string scale_fill_manual theme_bw stat_smooth geom_boxplot ggsave scale_size_area
#' @importFrom stringi stri_join stri_replace_all_fixed stri_sub
#' @importFrom dplyr select distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join funs summarise_at bind_rows
#' @importFrom readr write_tsv
#' @importFrom tibble data_frame has_name
#' @importFrom tidyr complete gather unite spread nesting

#' @details To help choose a threshold for the local and global MAF
#' use the interactive version.
#'
#' 2 steps in the interactive version:
#'
#' Step 1. Global and Local MAF visualization and helper table.
#'
#' Step 2. Filtering markers based on the different MAF arguments


#' @return With \code{interactive.filter = FALSE}, a list in the global environment,
#' with 7 objects:
#' \enumerate{
#' \item $tidy.filtered.maf
#' \item $whitelist.markers
#' \item $blacklist.markers
#' \item $maf.data.thresholds
#' \item $filters.parameters
#' }
#'
#' With \code{interactive.filter = TRUE}, a list with 4 additionnal objects are generated.
#' \enumerate{
#' \item $distribution.maf.global
#' \item $distribution.maf.local
#' \item $maf.global.summary
#' \item $maf.helper.table
#' }
#'
#' \strong{maf.helper.table:}
#'
#' First and second variables (in columns) represents POP_ID and sample size (n).
#'
#' The last 2 rows are the local MAF (suggested based on the lowest pop value)
#' and the TOTAL/GLOBAL observations.
#'
#' Columns starting with ALT are the variable corresponding
#' to the number of alternative (ALT) allele (ranging from 1 to 20).
#' The observations in the ALT allele variable columns are the local (for the pop)
#' and global (last row) MAF of your dataset.
#' e.g. ALT_3 can potentially represent 3 heterozygote individuals with
#' the ALT allele or 1 homozygote individuals for the ALT allele and
#' 1 heterozygote individual. And so on...


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

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

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


  # Check pop.levels, pop.labels and pop.select---------------------------------
  check.levels <- check_pop_levels(pop.levels = pop.levels, pop.labels = pop.labels, pop.select = pop.select)
  pop.levels <- check.levels$pop.levels
  pop.labels <- check.levels$pop.labels
  pop.select <- check.levels$pop.select
  check.levels <- NULL
  # Message about steps taken during the process ---------------------------------
  if (interactive.filter) {
    message("Interactive mode: on\n")
    message("Step 1. Global and Local MAF visualization and helper table")
    message("Step 2. Filtering markers based on the different MAF arguments\n\n")
  }
  # Folder -------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")# Date and time
  folder.extension <- stringi::stri_join("filter_maf_", file.date)
  path.folder <- file.path(getwd(), folder.extension)
  dir.create(path.folder)
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
    maf.approach <- "SNP"
    # confusing, but because the haplotpe file doesn't have snp/locus info,
    # the data is filtered the same way as the approach by SNP.
  }

  if (!interactive.filter && maf.approach == "haplotype") {
    if (!tibble::has_name(input, "LOCUS") || !tibble::has_name(input, "POS")) {
      stop("The haplotype approach during MAF filtering is for datasets with LOCUS
           and POS (SNP) info. e.g. VCF, haplotype VCF and stacks haplotypes file.
           Use the snp approach for the other file types")
    }
  }

  # import data ----------------------------------------------------------------
  message("Importing data ...")
  if (data.type %in% c("tbl_df", "fst.file")) {
    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "INDIVIDUALS", "POP_ID", "REF", "ALT", "GT", "GT_VCF", "GT_VCF_NUC", "GT_BIN")
    if (data.type == "tbl_df") {
      input <- suppressWarnings(dplyr::select(data, dplyr::one_of(want)))
      data <- NULL
    }
    if (data.type == "fst.file") {
      import.col <- colnames(fst::read.fst(path = data, from = 1, to = 1))
      import.col <- purrr::discard(.x = import.col, .p = !import.col %in% want)
      input <- fst::read.fst(path = data, columns = import.col)
      import.col <- want <- NULL
    }
  } else {
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
      parallel.core = parallel.core,
      verbose = FALSE
    )
  }

  # keeping markers meta -------------------------------------------------------
  markers.meta <- suppressWarnings(
    dplyr::select(input, dplyr::one_of(c("MARKERS", "CHROM", "LOCUS", "POS")))
  ) %>%
    dplyr::distinct(MARKERS, .keep_all = TRUE)

  # strata ---------------------------------------------------------------------
  strata.df <- input %>%
    dplyr::select(INDIVIDUALS, POP_ID) %>%
    dplyr::distinct(INDIVIDUALS, .keep_all = TRUE)

  pop.levels <- levels(input$POP_ID)
  pop.labels <- pop.levels

  # MAF calculation ------------------------------------------------------------
  message("Calculating global and local MAF on large data set may take some time...")
  if (data.type == "haplo.file") {
    biallelic <- FALSE
  } else {
    biallelic <- radiator::detect_biallelic_markers(input, verbose = TRUE)
  }
  markers.df <- dplyr::distinct(input, MARKERS)
  n.markers <- nrow(markers.df)
  maf.data <- dplyr::filter(input, GT != "000000")

  if (n.markers > 10000) {
    split.vec <- markers.df %>%
      dplyr::mutate(SPLIT_VEC = split_vec_row(
        markers.df,
        cpu.rounds = ceiling(n.markers/10000),
        parallel.core = parallel.core))

    maf.data <- dplyr::left_join(maf.data, split.vec, by = "MARKERS") %>%
      split(x = ., f = .$SPLIT_VEC) %>%
      .radiator_parallel_mc(
        X = .,
        FUN = compute_maf,
        mc.cores = parallel.core,
        biallelic = biallelic
      ) %>%
      dplyr::bind_rows(.)
    markers.df <- split.vec <- NULL
  } else {
    maf.data <- compute_maf(x = maf.data, biallelic = biallelic)
  }

  # Step 1. Global and Local MAF------------------------------------------------
  if (interactive.filter) {
    message("\nStep 1. Global and Local MAF visualization and helper table\n")
  }
  # plot
  OVERALL <- NULL

  if (tibble::has_name(maf.data, "HAPLOTYPES")) {
    global.data <- dplyr::ungroup(maf.data) %>%
      dplyr::filter(REF == "MAF") %>%
      dplyr::select(MARKERS, MAF_GLOBAL) %>%
      dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
      dplyr::mutate(OVERALL = rep("overall", n()))
  } else {
    global.data <- dplyr::ungroup(maf.data) %>%
      dplyr::distinct(MARKERS, MAF_GLOBAL) %>%
      dplyr::mutate(OVERALL = rep("overall", n()))
  }

  maf.global.summary <- dplyr::ungroup(global.data) %>%
    dplyr::summarise(
      MEAN = mean(MAF_GLOBAL, na.rm = TRUE),
      MEDIAN = stats::median(MAF_GLOBAL, na.rm = TRUE),
      RANGE = stringi::stri_join(round(min(MAF_GLOBAL, na.rm = TRUE), 4), " - ", round(max(MAF_GLOBAL, na.rm = TRUE), 4))
    )

  histo.maf.global <- ggplot2::ggplot(global.data, ggplot2::aes(x = MAF_GLOBAL)) +
    ggplot2::geom_histogram(bins = 30) +
    ggplot2::labs(y = "Number of markers", x = "Global MAF") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
    )
  print(histo.maf.global)
  # save
  ggplot2::ggsave(
    filename = file.path(path.folder, "distribution.maf.global.pdf"),
    plot = histo.maf.global, width = 10, height = 10, dpi = 600, units = "cm", useDingbats = FALSE)
  message("Distribution of global MAF written in the folder")

  global.data <- NULL # unused object

  readr::write_tsv(x = maf.global.summary, path = file.path(path.folder, "maf.global.summary.tsv"))
  message("maf.global.summary.tsv table saved in the folder")

  if (interactive.filter) {
    message("The global MAF mean: ", round(maf.global.summary$MEAN, 4))
    message("The global MAF median: ", round(maf.global.summary$MEDIAN, 4))
    message("The global MAF range: ", maf.global.summary$RANGE)
  }

  # Local MAF
  pop.number <- length(levels(maf.data$POP_ID))

  if (tibble::has_name(maf.data, "HAPLOTYPES")) {
    histo.maf.local <- maf.data %>%
      dplyr::filter(REF == "ALT") %>%
      dplyr::select(MARKERS, POP_ID, MAF_LOCAL) %>%
      dplyr::distinct(MARKERS, POP_ID, .keep_all = TRUE)
  } else {
    histo.maf.local <- maf.data
  }

  histo.maf.local <- ggplot2::ggplot(data = histo.maf.local, ggplot2::aes(x = MAF_LOCAL, na.rm = FALSE)) +
    # ggplot2::geom_line(ggplot2::aes(y = ..scaled.., color = POP_ID), stat = "density", adjust = 1) + # pop colored
    ggplot2::geom_histogram(bins = 30) +
    ggplot2::scale_x_continuous(breaks = seq(0,1, by = 0.1)) +
    ggplot2::labs(x = "Minor Allele Frequency (MAF)", y = "Number of markers") +
    # ggplot2::labs(y = "Density of MARKERS (scaled)") +
    ggplot2::expand_limits(y = 0) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
      legend.position = "none",
      axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5),
      strip.text.y = ggplot2::element_text(angle = 0, size = 12, family = "Helvetica", face = "bold"),
      strip.text.x = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold")
    ) +
    ggplot2::facet_grid(~POP_ID)
  print(histo.maf.local)
  ggplot2::ggsave(
    filename = file.path(path.folder, "maf.local.spectrums.pdf"),
    plot = histo.maf.local, width = pop.number * 5, height = 15,
    dpi = 600, units = "cm", useDingbats = FALSE, limitsize = FALSE)
  message("Local maf spectrums written in the folder")

  # Helper table for global and local MAF --------------------------------------
  message("Generating MAF helper table...")
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

  new <- NULL# remove unused objects

  readr::write_tsv(
    x = maf.helper.table,
    path = file.path(path.folder, "maf.helper.table.tsv"))

  message("\nWritten in the directory: maf.helper.table.tsv")

  # Step 2. Thresholds selection -----------------------------------------------
  # maf.approach
  if (interactive.filter) {
    message("\nStep 2. Filtering markers based on the different MAF arguments\n")
    message("The maf.approach:\n
maf.approach = \"haplotype\" : looks at the frequency of all ALTernative allele on the locus.
maf.approach = \"SNP\" : SNPs on the same haplotype/read are considered independent.")
    message("Choose the maf.approach (SNP/haplotype):")
    maf.approach <- as.character(readLines(n = 1))
    if (!maf.approach %in% c("SNP", "haplotype")) stop("maf.approach: SNP or haplotype")
  }

  if (maf.approach == "haplotype") {
    if (!tibble::has_name(input, "LOCUS") || !tibble::has_name(input, "POS")) {
      maf.approach <- "SNP"
      message("maf.approach = 'haplotype', switched to maf.approach = 'SNP',
because LOCUS and POS (SNP) info is not available")
      message("\nFor haplotype data, this is not important: it's all about locus/haplotypes")
    }
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
Option 1: AND: local \"AND\" global MAF thresholds are required to pass (more severe).
Option 2: OR: local \"OR\" global MAF thresholds are required to pass (more tolerant).")
    message("Choose the maf.operator (AND/OR):")
    maf.operator <- as.character(readLines(n = 1))
    if (!maf.operator %in% c("OR", "AND")) stop("maf.operator: either OR/AND")
  }

  # maf.pop.num.threshold
  if (interactive.filter) {
    message("Last threshold... \n\nmaf.pop.num.threshold:
How many populations are required to pass all the thresholds to keep the marker?\n
Example: if you have 10 populations and choose maf.pop.num.threshold = 3,
3 populations out of 10 are required to pass previous thresholds")
    message("Choose the value (integer) for the maf.pop.num.threshold:")
    maf.pop.num.threshold <- as.numeric(readLines(n = 1))
  }

  if (!interactive.filter) {
    maf.local.threshold <- maf.thresholds[1]
    maf.global.threshold <- maf.thresholds[2]
  }

  # Update the maf.data with pass or not filter based on threshold -------------
  message("Compiling all the statistics and thresholds")
  if (tibble::has_name(maf.data, "HAPLOTYPES")) {
    maf.data.thresholds <- maf.data %>%
      dplyr::filter(REF == "MAF" ) %>%
      dplyr::mutate(
        OR = dplyr::if_else((MAF_LOCAL >= maf.local.threshold | MAF_GLOBAL >= maf.global.threshold), "pass", "pruned"),
        AND = dplyr::if_else((MAF_LOCAL >= maf.local.threshold & MAF_GLOBAL >= maf.global.threshold), "pass", "pruned")
      ) %>%
      dplyr::group_by(MARKERS, HAPLOTYPES) %>%
      dplyr::mutate(
        OR = dplyr::if_else(length(POP_ID[OR == "pass"]) >= maf.pop.num.threshold, "pass", "pruned"),
        AND  = dplyr::if_else(length(POP_ID[AND == "pass"]) >= maf.pop.num.threshold, "pass", "pruned")
      ) %>%
      dplyr::select(POP_ID, MARKERS, HAPLOTYPES, OR, AND) %>%
      dplyr::right_join(maf.data, by = c("MARKERS", "POP_ID", "HAPLOTYPES")) %>%
      dplyr::mutate(
        OR = stringi::stri_replace_na(OR, "pass"),
        AND = stringi::stri_replace_na(AND, "pass")
      ) %>%
      dplyr::select(MARKERS, POP_ID, HAPLOTYPES, REF, MAF_LOCAL, OR, AND, MAF_GLOBAL)

  } else {
    maf.data.thresholds <- maf.data %>%
      dplyr::mutate(
        OR = dplyr::if_else((MAF_LOCAL >= maf.local.threshold | MAF_GLOBAL >= maf.global.threshold), "pass", "pruned"),
        AND = dplyr::if_else((MAF_LOCAL >= maf.local.threshold & MAF_GLOBAL >= maf.global.threshold), "pass", "pruned")
      ) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(
        OR = dplyr::if_else(length(POP_ID[OR == "pass"]) >= maf.pop.num.threshold, "pass", "pruned"),
        AND = dplyr::if_else(length(POP_ID[AND == "pass"]) >= maf.pop.num.threshold, "pass", "pruned")
      ) %>%
      dplyr::select(MARKERS, POP_ID, MAF_LOCAL, OR, AND, MAF_GLOBAL)
  }

  if (ncol(markers.meta) > 1) {
    markers.meta.col <- colnames(markers.meta)
    maf.data.thresholds <- suppressWarnings(
      maf.data.thresholds %>%
        dplyr::left_join(markers.meta, by = "MARKERS") %>%
        dplyr::select(dplyr::one_of(markers.meta.col), dplyr::everything(.)))

    maf.data <- suppressWarnings(
      maf.data %>%
        dplyr::left_join(markers.meta, by = "MARKERS") %>%
        dplyr::select(dplyr::one_of(markers.meta.col), dplyr::everything(.)))
  }

  readr::write_tsv(
    x = maf.data.thresholds,
    path = file.path(path.folder, "maf.data.tsv"),
    col_names = TRUE,
    append = FALSE
  )
  message("\nWritten in the folder: maf.data.tsv\n")

  # Filtering ------------------------------------------------------------------
  if (maf.approach == "haplotype") {

    # Needed switched for haplotype file (e.g. stacks haplotype)
    if (tibble::has_name(maf.data, "MARKERS") && !tibble::has_name(maf.data, "LOCUS")) {
      maf.data <- dplyr::rename(maf.data, LOCUS = MARKERS)
      maf.data.thresholds <- dplyr::rename(maf.data.thresholds, LOCUS = MARKERS)
    }

    filter <- maf.data %>%
      dplyr::group_by(LOCUS, POP_ID) %>%
      dplyr::summarise(
        MAF_GLOBAL = mean(MAF_GLOBAL, na.rm = TRUE),
        MAF_LOCAL = min(MAF_LOCAL, na.rm = TRUE)
      )

    if (maf.operator == "OR") {
      filter <- dplyr::filter(
        filter,
        MAF_LOCAL >= maf.local.threshold | MAF_GLOBAL >= maf.global.threshold)
    } else {# AND operator between local and global maf
      filter <- dplyr::filter(
        filter,
        MAF_LOCAL >= maf.local.threshold & MAF_GLOBAL >= maf.global.threshold)
    }
    filter <- filter %>%
      dplyr::group_by(LOCUS) %>%
      dplyr::tally(.) %>%
      dplyr::filter(n >= maf.pop.num.threshold) %>%
      dplyr::select(LOCUS) %>%
      dplyr::left_join(input, by = "LOCUS") %>%
      dplyr::arrange(LOCUS, POP_ID)
  } # end maf haplotype approach

  if (maf.approach == "SNP") { # SNP approach

    filter <- maf.data %>%
      dplyr::group_by(MARKERS, POP_ID) %>%
      dplyr::summarise(
        MAF_GLOBAL = mean(MAF_GLOBAL, na.rm = TRUE),
        MAF_LOCAL = min(MAF_LOCAL, na.rm = TRUE)
      )
    if (maf.operator == "OR") {
      filter <- dplyr::filter(
        filter,
        MAF_LOCAL >= maf.local.threshold | MAF_GLOBAL >= maf.global.threshold)
    } else {# AND operator between local and global maf
      filter <- dplyr::filter(
        filter,
        MAF_LOCAL >= maf.local.threshold & MAF_GLOBAL >= maf.global.threshold)
    }
    filter <- dplyr::group_by(filter, MARKERS) %>%
      dplyr::tally(.) %>%
      dplyr::filter(n >= maf.pop.num.threshold) %>%
      dplyr::select(MARKERS) %>%
      dplyr::left_join(input, by = "MARKERS") %>%
      dplyr::arrange(MARKERS, POP_ID)

    # For haplotypes file requiring reconstruction------------------------------
    if (tibble::has_name(maf.data, "HAPLOTYPES")) {
      message("Haplotype reconstruction analysis")
      haplo.reconstruction <- dplyr::distinct(filter, MARKERS) %>%
        dplyr::left_join(maf.data, by = "MARKERS") %>%
        dplyr::distinct(MARKERS, HAPLOTYPES) %>%
        dplyr::arrange(MARKERS, HAPLOTYPES) %>%
        dplyr::mutate(SNP_N = stringi::stri_length(str = HAPLOTYPES)) %>%
        dplyr::filter(SNP_N > 1) %>%
        haplotype_reconstruction(data = ., parallel.core = parallel.core) %>%
        dplyr::mutate(DIFF = dplyr::if_else(HAPLOTYPES == HAPLOTYPES_NEW, "same", "different")) %>%
        dplyr::filter(DIFF == "different") %>%
        dplyr::select(-DIFF)

      n.reconstruction <- nrow(haplo.reconstruction)

      if (n.reconstruction > 0) {
        message("Number of markers requiring reconstruction: ", n.reconstruction)
        # no.reconstruction (written to disk to save mem)
        filter %>%
          dplyr::filter(!MARKERS %in% haplo.reconstruction$MARKERS) %>%
          fst::write.fst(x = ., path = file.path(path.folder, "temp.rad"), compress = 85)
        filter <- dplyr::select(filter, MARKERS, POP_ID, INDIVIDUALS, GT_VCF_NUC) %>%
          dplyr::filter(MARKERS %in% haplo.reconstruction$MARKERS) %>%
          separate_gt(
            x = .,
            sep = "/", gt = "GT_VCF_NUC",
            exclude = c("MARKERS", "INDIVIDUALS", "POP_ID"),
            cpu.rounds = 10, parallel.core = parallel.core) %>%
          dplyr::left_join(haplo.reconstruction, by = c("MARKERS", "HAPLOTYPES"))
        haplo.reconstruction <- NULL
        filter <- filter %>%
          dplyr::select(-HAPLOTYPES, HAPLOTYPES = HAPLOTYPES_NEW) %>%
          dplyr::group_by(MARKERS, INDIVIDUALS, POP_ID) %>%
          dplyr::summarise(GT_VCF_NUC = stringi::stri_join(HAPLOTYPES, collapse = "/")) %>%
          dplyr::mutate(GT_VCF_NUC = stringi::stri_replace_na(str = GT_VCF_NUC, replacement = "./.")) %>%
          dplyr::ungroup(.) %>%
          change_alleles(data = ., biallelic = FALSE, parallel.core = parallel.core)

        filter <- filter$input %>%
          dplyr::select(MARKERS, POP_ID, INDIVIDUALS, GT_VCF_NUC, REF, ALT, GT, GT_VCF) %>%
          dplyr::bind_rows(fst::read.fst(path = file.path(path.folder, "temp.rad")))
        file.remove(file.path(path.folder, "temp.rad"))
      } else {
        haplo.reconstruction <- NULL
      }
      maf.data <- n.reconstruction <- NULL
      # Private haplotypes before filter -------------------------------------------
      old.dir <- getwd()
      setwd(path.folder)
      private.haplo.before <- private_haplotypes(data = filter, verbose = FALSE)
      setwd(old.dir)
    }
  } # end maf snp approach

  # Update filters.parameters SNP ----------------------------------------------
  snp.before <- as.integer(dplyr::n_distinct(input$MARKERS))
  snp.after <- as.integer(dplyr::n_distinct(filter$MARKERS))
  snp.blacklist <- snp.before - snp.after

  if (tibble::has_name(input, "LOCUS") && maf.approach == "haplotype") {
    locus.before <- as.integer(dplyr::n_distinct(input$LOCUS))
    locus.after <- as.integer(dplyr::n_distinct(filter$LOCUS))
    locus.blacklist <- locus.before - locus.after
  } else {
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
    fst::write.fst(x = filter, path = file.path(path.folder, tidy.name), compress = 85)
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
  readr::write_tsv(whitelist.markers, file.path(path.folder, "whitelist.markers.maf.tsv"), append = FALSE, col_names = TRUE)


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
  readr::write_tsv(blacklist.markers, file.path(path.folder, "blacklist.markers.maf.tsv"), append = FALSE, col_names = TRUE)

  # results --------------------------------------------------------------------
  cat("############################### RESULTS ###############################\n")
  message("maf.approach: ", maf.approach)
  message("maf.thresholds: ", "local = ", maf.local.threshold, ", global = ", maf.global.threshold)
  message("maf.operator: ", maf.operator)
  message("maf.pop.num.threshold: ", maf.pop.num.threshold)
  if (tibble::has_name(input, "LOCUS") & maf.approach == "haplotype") {
    message("The number of markers removed by the MAF filter:\nSNP: ", snp.blacklist, "\nLOCUS: ", locus.blacklist)
    message("The number of markers before -> after the MAF filter")
    message("SNP: ", snp.before, " -> ", snp.after)
    message("LOCUS: ", locus.before, " -> ", locus.after)
  } else {
    message("The number of markers removed by the MAF filter: ", snp.blacklist)
    message("The number of markers before -> after the MAF filter")
    message("SNP: ", snp.before, " -> ", snp.after)
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
  res$filters.parameters <- filters.parameters
  res$distribution.maf.global <- histo.maf.global
  res$distribution.maf.local <- histo.maf.local
  res$maf.global.summary <- maf.global.summary
  return(res)
}

# Internal nested functon ------------------------------------------------------
#' @title compute_maf
#' @description Compute MAF
#' @rdname compute_maf
#' @keywords internal
#' @export
compute_maf <- function(x, biallelic) {
  if (tibble::has_name(x, "GT_BIN") && biallelic) {
    maf.data <- x %>%
      dplyr::group_by(MARKERS, POP_ID) %>%
      dplyr::summarise(
        NN = as.numeric(2 * n()),
        PP = as.numeric(2 * length(GT_BIN[GT_BIN == 0])),
        PQ = as.numeric(length(GT_BIN[GT_BIN == 1])),
        QQ = as.numeric(2 * length(GT_BIN[GT_BIN == 2]))
      ) %>%
      # need this step because seen cases where the minor allele is not minor
      dplyr::mutate(
        PP = PP + PQ,
        QQ = QQ + PQ,
        PQ = NULL) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(
        NN_G = sum(NN),
        PP_G = sum(PP),
        QQ_G = sum(QQ)) %>%
      dplyr::group_by(MARKERS, POP_ID) %>%
      dplyr::mutate(
        ALT = dplyr::if_else(PP_G < QQ_G, PP, QQ),
        MAF_LOCAL = (ALT / NN),
        PP = NULL,
        QQ = NULL,
        NN = NULL) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(
        ALT = dplyr::if_else(PP_G < QQ_G, PP_G, QQ_G),
        MAF_GLOBAL = (ALT / NN_G),
        ALT = NULL,
        PP_G = NULL,
        QQ_G = NULL,
        NN_G = NULL) %>%
      dplyr::ungroup(.)
  } else {
    if (!tibble::has_name(x, "GT_VCF_NUC")) {
      maf.data <- x %>%
        dplyr::select(MARKERS,POP_ID, INDIVIDUALS, GT) %>%
        dplyr::mutate(
          A1 = stringi::stri_sub(GT, 1, 3),
          A2 = stringi::stri_sub(GT, 4,6)
        ) %>%
        dplyr::select(MARKERS, POP_ID, INDIVIDUALS, A1, A2) %>%
        tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS, POP_ID))

      maf.local <- maf.data %>%
        dplyr::group_by(MARKERS, POP_ID, GT) %>%
        dplyr::tally(.) %>%
        dplyr::ungroup(.) %>%
        tidyr::complete(data = ., POP_ID, tidyr::nesting(MARKERS, GT), fill = list(n = 0)) %>%
        dplyr::group_by(MARKERS, POP_ID) %>%
        dplyr::mutate(n.al.tot = sum(n)) %>%
        dplyr::filter(n == min(n)) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
        dplyr::summarise(MAF_LOCAL = n / n.al.tot) %>%
        dplyr::ungroup(.) %>%
        dplyr::select(MARKERS, POP_ID, MAF_LOCAL)

      maf.data <- maf.data %>%
        dplyr::group_by(MARKERS, GT) %>%
        dplyr::tally(.) %>%
        dplyr::group_by(MARKERS) %>%
        dplyr::mutate(n.al.tot = sum(n)) %>%
        dplyr::filter(n == min(n)) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
        dplyr::summarise(MAF_GLOBAL = n / n.al.tot) %>%
        dplyr::ungroup(.) %>%
        dplyr::select(MARKERS, MAF_GLOBAL) %>%
        dplyr::left_join(maf.local, by = c("MARKERS")) %>%
        dplyr::select(MARKERS, POP_ID, MAF_LOCAL, MAF_GLOBAL)
    } else {
      maf.data <- x %>%
        dplyr::select(MARKERS, POP_ID, INDIVIDUALS, GT_VCF_NUC) %>%
        tidyr::separate(
          data = .,
          col = GT_VCF_NUC, into = c("A1", "A2"),
          sep = "/",
          extra = "drop", remove = TRUE
        ) %>%
        tidyr::gather(
          data = ., key = ALLELE_GROUP, value = HAPLOTYPES,
          -dplyr::one_of(c("MARKERS", "INDIVIDUALS", "POP_ID"))) %>%
        dplyr::select(-ALLELE_GROUP) %>%
        dplyr::group_by(MARKERS, HAPLOTYPES, POP_ID) %>%
        dplyr::tally(.) %>%
        dplyr::ungroup(.) %>%
        tidyr::complete(data = ., POP_ID, tidyr::nesting(MARKERS, HAPLOTYPES), fill = list(n = 0)) %>%
        dplyr::group_by(MARKERS, POP_ID) %>%
        dplyr::mutate(N_LOCAL = sum(n)) %>%
        dplyr::group_by(MARKERS) %>%
        dplyr::mutate(N_GLOBAL = sum(n)) %>%
        dplyr::arrange(MARKERS, POP_ID) %>%
        dplyr::group_by(MARKERS, POP_ID, HAPLOTYPES) %>%
        dplyr::mutate(MAF_LOCAL = n / N_LOCAL) %>%
        dplyr::group_by(MARKERS, HAPLOTYPES) %>%
        dplyr::mutate(
          MAF_GLOBAL = sum(n) / N_GLOBAL,
          n = NULL,
          N_LOCAL = NULL,
          N_GLOBAL = NULL
        ) %>%
        dplyr::ungroup(.)

      ref.info <- dplyr::distinct(maf.data, MARKERS, HAPLOTYPES, MAF_GLOBAL) %>%
        dplyr::group_by(MARKERS) %>%
        dplyr::filter(MAF_GLOBAL == max(MAF_GLOBAL)) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
        dplyr::ungroup(.) %>%
        dplyr::mutate(REF = rep("REF", n()), MAF_GLOBAL = NULL) %>%
        dplyr::bind_rows(
          dplyr::distinct(maf.data, MARKERS, HAPLOTYPES, MAF_GLOBAL) %>%
            dplyr::group_by(MARKERS) %>%
            dplyr::filter(MAF_GLOBAL == min(MAF_GLOBAL)) %>%
            dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
            dplyr::ungroup(.) %>%
            dplyr::mutate(REF = rep("MAF", n()), MAF_GLOBAL = NULL)
        )

      maf.data <- dplyr::left_join(maf.data, ref.info, by = c("MARKERS", "HAPLOTYPES")) %>%
        dplyr::mutate(REF = stringi::stri_replace_na(REF, replacement = "ALT"))
      ref.info <- NULL

      # maf.local <- maf.data %>%
      #   dplyr::group_by(MARKERS, HAPLOTYPES, POP_ID) %>%
      #   dplyr::tally(.) %>%
      #   dplyr::ungroup(.) %>%
      #   tidyr::complete(data = ., POP_ID, tidyr::nesting(MARKERS, HAPLOTYPES), fill = list(n = 0)) %>%
      #   dplyr::group_by(MARKERS, POP_ID) %>%
      #   dplyr::mutate(n.al.tot = sum(n)) %>%
      #   dplyr::filter(n == min(n)) %>%
      #   dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
      #   dplyr::summarise(MAF_LOCAL = n / n.al.tot) %>%
      #   dplyr::ungroup(.) %>%
      #   dplyr::select(MARKERS, POP_ID, MAF_LOCAL)

      # maf.data <- maf.data %>%
      #   dplyr::group_by(MARKERS, HAPLOTYPES) %>%
      #   dplyr::tally(.) %>%
      #   dplyr::mutate(n.al.tot = sum(n)) %>%
      #   dplyr::ungroup(.)
      #
      # ref.allele <- maf.data %>%
      #   dplyr::group_by(MARKERS) %>%
      #   dplyr::filter(n == max(n)) %>%
      #   dplyr::ungroup(.) %>%
      #   dplyr::select(MARKERS, REF = HAPLOTYPES) %>%
      #   dplyr::distinct(MARKERS, .keep_all = TRUE)
      #
      # maf.data <- maf.data %>%
      #   dplyr::group_by(MARKERS) %>%
      #   dplyr::filter(n == min(n)) %>%
      #   dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
      #   dplyr::summarise(MAF_GLOBAL = n / n.al.tot) %>%
      #   dplyr::ungroup(.) %>%
      #   dplyr::select(MARKERS, MAF_GLOBAL) %>%
      #   dplyr::left_join(maf.local, by = c("MARKERS")) %>%
      #   dplyr::left_join(ref.allele, by = c("MARKERS")) %>%
      #   dplyr::select(MARKERS, POP_ID, REF, MAF_LOCAL, MAF_GLOBAL)
    }
  }
  maf.local <- NULL
  return(maf.data)
}#End maf_local
