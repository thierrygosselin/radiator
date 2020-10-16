# Minor Allele Frequency
#' @name filter_maf
#' @title MAF/MAC filter
#' @description Minor Allele Filter.
#' Remove markers based on Minor/Alternate Allele Frequency (MAF) or Count (MAC).
#'
#' Most arguments are inherited from
#' \code{\link[radiator]{tidy_genomic_data}}.
#' @inheritParams tidy_genomic_data

#' @param interactive.filter (optional, logical) Do you want the filtering session to
#' be interactive. With default: \code{interactive.filter == TRUE}, figures and
#' tables are shown before making decisions for filtering.

#' @param maf.thresholds (optional, string) 5 values required in the string.
#'
#' Example using maf.thresholds with SNP and frequency approach:
#' \code{maf.thresholds = c("SNP", 0.001, "OR", 0.001, 1)}.
#'
#' Example using maf.thresholds with locus and count approach:
#' \code{maf.thresholds = c("locus", 3, "OR", 6, 1)}.
#'
#'
#' \enumerate{
#' \item \code{maf approach} (character: "SNP"/"locus"):
#' MAF filtering is conducted by SNPs or locus.
#' \code{"SNP"}: will consider all the SNPs on the same locus/read as independent
#' and will be filtered independently of their locus id.
#' \code{"locus"}: looks at the minimum MAF found on the
#' read/locus. Using this option will discard all the markers/snp on
#' that read based on the thresholds chosen. For the locus approach to work, your dataset
#' requires SNP and Locus info (e.g. from a VCF file).
#'
#' \item \code{local} threshold (double or integer): For a frequency threshold use
#' a double (e.g. 0.05). For a count threshold, use an integer
#' (e.g. 3 for the number of alternate allele required in a population). This
#' threshold is applied by population.
#' Not sure about the threshold to use, choose the interactive mode argument.
#'
#' \item \code{operator} (character: "OR" / "AND"):
#' To consider both the local \code{"AND"} the global thresholds, use \code{"AND"}.
#' To consider the local \code{"OR"} the global thresholds, use \code{"OR"}.
#'
#' \item \code{global} threshold (double or integer): For a frequency threshold
#' use a double (e.g. 0.02). For a count threshold, use an integer
#' (e.g. 6 for the number of alternate allele required). This threshold is
#' applied at the dataset level (no population).
#' Not sure about the threshold to use, choose the interactive mode argument.
#'
#'\item \code{maf pop num threshold} (integer)
#'The number of pop required to pass the thresholds
#' to keep the marker. Usually, I always use \code{1}, for 1 pop is required to
#' pass thresholds and keep the marker.
#' }


#' @param filename (optional, character) Write to folder the MAF filtered
#' tidy dataset. The name will be appended \code{.rad}. With default, the filtered
#' data is only in the global environment.
#' Default: \code{filename = NULL}.

#' @param ... (optional) To pass further arguments for fine-tuning the function
#' and legacy arguments (see Life cycle section below).

#' @section Life cycle:
#'
#' The \strong{maf.thresholds} arguments used in the past in several of
#' radiator functions
#' is still available inside \code{filter_maf} arguments.
#' But it is in evaluation mode and the function will probaly be deprecated in
#' future releases. The arguments are complicated to use and
#' most RADseq datesets encountered so far are better filtered
#' using \code{\link{filter_mac}}.

#' @section MAC or MAF ?:
#' Using count or frequency to remove a SNPs ? The preferred choice in radiator
#' as changed from frequency to count, because we think the filtering should not
#' alter the spectrum and this is only achieved if the same criteria is applied
#' for each SNP.
#'
#'
#' Even small differences in missing data between RADseq markers
#' generates differences in MAF frequency thresholds applied.
#'
#' Example with a datset consisting of N = 36 individuals and 3 SNPs
#' with varying level of missing genotypes:
#'
#' \itemize{
#' \item \code{SNP number : number samples genotypes : REF/ALT counts}
#' \item SNP1 : 36 : 69/3
#' \item SNP2 : 30 : 65/3
#' \item SNP3 : 24 : 45/3
#' }
#'
#' Each SNPs have the same alternate allele count, corresponding to
#' 2 individuals with the polymorphism: 1 homozygote + 1 heterozygote.
#' Applying a MAF threshold of
#' 0.05 would mean that SNP3 would be blacklisted
#' (\code{24 * 2 * 0.05 = 2.4 alt alleles required to pass}).
#'
#'
#' \strong{Using count instead of frequency allows each RADseq markers,
#' with varying missing data, to be treated equally.}


#' @section Interactive version:
#'
#' To help choose a threshold for the local and global MAF
#' use the interactive version.
#'
#' 2 steps in the interactive version:
#'
#' Step 1. Global and Local MAF visualization and helper table.
#'
#' Step 2. Filtering markers based on the different MAF arguments


#' @rdname filter_maf
#' @export

#' @return With \code{interactive.filter = FALSE}, a list in the global environment,
#' with 7 objects:
#' \enumerate{
#' \item $tidy.filtered.maf
#' \item $whitelist.markers
#' \item $blacklist.markers
#' \item $maf.data
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
#' turtle.maf <- radiator::filter_maf(
#'         data = "turtle.vcf",
#'         strata = "turtle.strata.tsv")
#' # This will use the default: interactive version,
#' # a list is created and to view the filtered tidy data:
#' tidy.data <- turtle.maf$tidy.filtered.maf
#'
#' # The remaining argument are used in tidy_genomic_data during import and allow
#' # the user to apply some filtering or selection before doing the MAF filtering.
#'
#' # If I want to filter Alternate alleles based on count.
#' turtle.maf <- radiator::filter_maf(
#'         data = "turtle.vcf",
#'         strata = "turtle.strata.tsv",
#'         maf.thresholds = c("SNP", 3, "OR", "5", 1),
#'         filename = "turtle.maf")
#'
#' # This will remove monomorphic markers (by default filter.monomorphic = TRUE)
#' # it will use the SNP approache (SNPs are considered independent and
#' # filtered independent)
#' # I keep markers if: they have a local count of 3 alternate allele OR global
#' # count of 5 alternate allele.
#' # I keep the marker if at leat 1 pop pass that.
#' # finally the filtered data will be written in the directory under the name:
#' # turtle.maf.rad
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

filter_maf <- function(
  interactive.filter = TRUE,
  data,
  strata = NULL,
  maf.thresholds = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
  ...
) {

  # testing
  # strata = NULL
  # interactive.filter = TRUE
  # maf.thresholds = NULL
  # filter.monomorphic = TRUE
  # common.markers = FALSE
  # snp.ld = NULL
  # pop.select = NULL
  # filename = NULL
  # parallel.core = parallel::detectCores() - 1
  # verbose = TRUE


  if (!is.null(maf.thresholds)) {
    if (interactive.filter) verbose <- TRUE
    if (verbose) cat("#######################################################################\n")
    if (verbose) cat("######################## radiator::filter_maf #########################\n")
    if (verbose) cat("#######################################################################\n")
    opt.change <- getOption("width")
    options(width = 70)
    timing <- proc.time()

    # manage missing argument -----------------------------------------------------
    if (missing(data)) rlang::abort("Input file missing")
    if (!interactive.filter) {
      if (length(maf.thresholds) < 5) rlang::abort("Please read MAF related arguments in the function documentation")
    }
    # dot dot dot ----------------------------------------------------------------
    dotslist <- list(...)
    unknowned_param <- setdiff(
      names(dotslist),
      c("whitelist.markers", "blacklist.id", "blacklist.genotype", "pop.levels",
        "pop.labels", "vcf.metadata", "filter.short.ld", "filter.monomorphic",
        "filter.common.markers"))

    if (length(unknowned_param) > 0) {
      rlang::abort("Unknowned \"...\" parameters to filter_maf: ",
           stringi::stri_join(unknowned_param, collapse = " "))
    }
    radiator.dots <- dotslist[
      names(dotslist) %in%
        c("whitelist.markers", "blacklist.id", "blacklist.genotype", "pop.levels",
          "pop.labels", "vcf.metadata", "filter.short.ld", "filter.monomorphic",
          "filter.common.markers")
      ]

    if (!is.null(radiator.dots[["whitelist.markers"]])) {
      whitelist.markers <- radiator.dots[["whitelist.markers"]]
    } else {
      whitelist.markers = NULL
    }
    if (!is.null(radiator.dots[["blacklist.id"]])) {
      blacklist.id <- radiator.dots[["blacklist.id"]]
    } else {
      blacklist.id = NULL
    }

    if (!is.null(radiator.dots[["blacklist.genotype"]])) {
      blacklist.genotype <- radiator.dots[["blacklist.genotype"]]
    } else {
      blacklist.genotype = NULL
    }

    if (!is.null(radiator.dots[["pop.levels"]])) {
      pop.levels <- radiator.dots[["pop.levels"]]
    } else {
      pop.levels = NULL
    }
    if (!is.null(radiator.dots[["pop.labels"]])) {
      pop.labels <- radiator.dots[["pop.labels"]]
    } else {
      pop.labels = NULL
    }

    if (!is.null(radiator.dots[["vcf.metadata"]])) {
      vcf.metadata <- radiator.dots[["vcf.metadata"]]
    } else {
      vcf.metadata = FALSE
    }
    filter.short.ld <- radiator.dots[["filter.short.ld"]]
    filter.monomorphic <- radiator.dots[["filter.monomorphic"]]
    filter.common.markers <- radiator.dots[["filter.common.markers"]]

    # whitelist.markers <- blacklist.id <- blacklist.genotype <- pop.levels <- pop.labels <- vcf.metadata <- NULL

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
    if (interactive.filter) {
      file.date <- format(Sys.time(), "%Y%m%d@%H%M")# Date and time
      folder.extension <- stringi::stri_join("filter_maf_", file.date)
      path.folder <- file.path(getwd(), folder.extension)
      dir.create(path.folder)
      if (verbose) message("\nFolder created: ", folder.extension)
      file.date <- NULL #unused object
    } else {
      path.folder <- getwd()
    }
    # Filter parameter file ------------------------------------------------------
    if (verbose) message("Parameters used in this run are stored in a file")
    filters.parameters <- list.files(path = getwd(), pattern = "filters_parameters.tsv", full.names = TRUE)
    if (length(filters.parameters) == 0) {
      filters.parameters <- tibble::data_frame(FILTERS = as.character(), PARAMETERS = as.character(), VALUES = as.integer(), BEFORE = as.character(), AFTER = as.character(), BLACKLIST = as.integer(), UNITS = as.character(), COMMENTS = as.character())
      readr::write_tsv(x = filters.parameters, path = "filters_parameters.tsv", append = FALSE, col_names = TRUE)
      if (verbose) message("    Created a parameter file: filters_parameters.tsv")
    } else {
      if (verbose) message("    Using the filters parameters file: filters_parameters.tsv")
    }
    # File type detection-------------------------------------------------------
    data.type <- radiator::detect_genomic_format(data)

    if (data.type == "haplo.file") {
      maf.approach <- "SNP"
      # confusing, but because the haplotpe file doesn't have snp/locus info,
      # the data is filtered the same way as the approach by SNP.
    }

    # Interactive mode details for MAF arguments -------------------------------
    if (!interactive.filter) {
      maf.approach <- as.character(maf.thresholds[1])
      maf.local.threshold <- maf.thresholds[2]
      if (maf.local.threshold >= 1) {
        maf.local.threshold <- as.integer(maf.local.threshold)
        maf.count <- TRUE
      } else {
        maf.local.threshold <- as.numeric(maf.local.threshold)
        maf.count <- FALSE
      }
      maf.operator <- as.character(maf.thresholds[3])
      maf.global.threshold <- maf.thresholds[4]
      if (maf.global.threshold >= 1) {
        maf.global.threshold <- as.integer(maf.global.threshold)
      } else {
        maf.global.threshold <- as.numeric(maf.global.threshold)
      }
      maf.pop.num.threshold <- as.numeric(maf.thresholds[5])
    }

    # import data --------------------------------------------------------------
    if (verbose) message("Importing data ...")
    if (is.null(blacklist.id) && is.null(blacklist.genotype) && is.null(whitelist.markers) &&
        !filter.monomorphic && is.null(filter.short.ld) && !filter.common.markers && is.null(pop.select)) {
      if (data.type %in% c("tbl_df", "fst.file")) {
        skip.import <- TRUE
      } else {
        skip.import <- FALSE
      }
    } else {
      skip.import <- FALSE
    }

    if (skip.import) {
      want <- c("MARKERS", "CHROM", "LOCUS", "POS", "INDIVIDUALS", "POP_ID", "REF", "ALT", "GT", "GT_VCF", "GT_VCF_NUC", "GT_BIN")
      if (data.type == "tbl_df") {
        input <- suppressWarnings(dplyr::select(data, dplyr::one_of(want)))
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
        filter.monomorphic = filter.monomorphic,
        filter.short.ld = filter.short.ld,
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
    data <- NULL

    # keeping markers meta -------------------------------------------------------
    markers.meta <- suppressWarnings(
      dplyr::select(input, dplyr::one_of(c("MARKERS", "CHROM", "LOCUS", "POS"))) %>%
      dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
      dplyr::mutate_all(.tbl = ., .funs = as.character))

    # strata ---------------------------------------------------------------------
    strata.df <- input %>%
      dplyr::select(INDIVIDUALS, POP_ID) %>%
      dplyr::distinct(INDIVIDUALS, .keep_all = TRUE)

    pop.levels <- levels(input$POP_ID)
    pop.labels <- pop.levels

    # biallelic data detection ---------------------------------------------------
    if (data.type == "haplo.file") {
      biallelic <- FALSE
    } else {
      biallelic <- radiator::detect_biallelic_markers(input, verbose = verbose)
    }

    # MAF calculation ------------------------------------------------------------
    if (verbose) message("Calculating global and local MAF")
    # Prepare for parallel computations
    markers.df <- dplyr::distinct(input, MARKERS)
    n.markers <- nrow(markers.df)

    if (tibble::has_name(input, "GT_BIN")) {
      maf.data <- dplyr::filter(input, !is.na(GT_BIN))
    } else {
      maf.data <- dplyr::filter(input, GT != "000000")
    }
    write_rad(data = input, path = "maf.temp.rad")
    # readr::write_tsv(x = input, path = "maf.temp.rad")
    input <- NULL

    if (n.markers > 10000) {
      split.vec <- markers.df %>%
        dplyr::mutate(SPLIT_VEC = split_vec_row(
          markers.df,
          cpu.rounds = ceiling(n.markers/10000),
          parallel.core = parallel.core))

      maf.data <- dplyr::left_join(maf.data, split.vec, by = "MARKERS") %>%
        split(x = ., f = .$SPLIT_VEC) %>%
        radiator_parallel_mc(
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

      if (tibble::has_name(maf.data, "HAPLOTYPES")) {
        global.data <- dplyr::ungroup(maf.data) %>%
          dplyr::filter(REF == "MAF") %>%
          dplyr::select(MARKERS, MAF_GLOBAL, ALT_GLOBAL) %>%
          dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
          dplyr::mutate(OVERALL = rep("overall", n()))
      } else {
        global.data <- dplyr::ungroup(maf.data) %>%
          dplyr::distinct(MARKERS, MAF_GLOBAL, ALT_GLOBAL) %>%
          dplyr::mutate(OVERALL = rep("overall", n()))
      }

      maf.global.summary <- dplyr::ungroup(global.data) %>%
        dplyr::summarise(
          MEAN = mean(MAF_GLOBAL, na.rm = TRUE),
          MEDIAN = stats::median(MAF_GLOBAL, na.rm = TRUE),
          RANGE = stringi::stri_join(round(min(MAF_GLOBAL, na.rm = TRUE), 4), " - ", round(max(MAF_GLOBAL, na.rm = TRUE), 4))
        )

      global.data <- tidyr::pivot_longer(
        data = global.data,
        cols = -c("MARKERS", "OVERALL"),
        names_to = "GROUP",
        values_to = "GLOBAL"
      ) %>%
        dplyr::mutate(GROUP = stringi::stri_replace_all_fixed(str = GROUP, pattern = c("ALT_GLOBAL", "MAF_GLOBAL"), replacement = c("ALT count", "ALT frequency (MAF)"), vectorize_all = FALSE))

      histo.maf.global <- ggplot2::ggplot(global.data, ggplot2::aes(x = GLOBAL)) +
        ggplot2::geom_histogram(bins = 30) +
        ggplot2::labs(y = "Number of markers", x = "Global Minor Allele (count and frequency)") +
        ggplot2::theme(
          legend.position = "none",
          axis.title.x = ggplot2::element_text(size = 10, face = "bold"),
          axis.title.y = ggplot2::element_text(size = 10, face = "bold")
        ) +
        ggplot2::theme_minimal() +
        ggplot2::facet_grid(~ GROUP, scales = "free")
      print(histo.maf.global)

      # save
      ggplot2::ggsave(
        filename = file.path(path.folder, "distribution.maf.global.pdf"),
        plot = histo.maf.global, width = 10, height = 10, dpi = 600, units = "cm", useDingbats = FALSE)
      if (verbose) message("Distribution of global MAF written in the folder")

      global.data <- NULL # unused object

      readr::write_tsv(x = maf.global.summary, path = file.path(path.folder, "maf.global.summary.tsv"))
      if (verbose) message("maf.global.summary.tsv table written in the folder")

      message("\nThe global MAF mean: ", round(maf.global.summary$MEAN, 4))
      message("The global MAF median: ", round(maf.global.summary$MEDIAN, 4))
      message("The global MAF range: ", maf.global.summary$RANGE)

      # Local MAF
      pop.number <- length(levels(maf.data$POP_ID))

      if (tibble::has_name(maf.data, "HAPLOTYPES")) {
        histo.maf.local <- maf.data %>%
          dplyr::filter(REF == "ALT") %>%
          dplyr::select(MARKERS, POP_ID, MAF_LOCAL, ALT_LOCAL) %>%
          dplyr::distinct(MARKERS, POP_ID, .keep_all = TRUE)
      } else {
        histo.maf.local <- maf.data
      }

      histo.maf.local <- histo.maf.local %>%
        dplyr::select(MARKERS, POP_ID, `ALT count` = ALT_LOCAL, `ALT frequency` = MAF_LOCAL) %>%
        tidyr::pivot_longer(
          data = .,
          cols = -c("POP_ID", "MARKERS"),
          names_to = "GROUP",
          values_to = "LOCAL"
        ) %>%
        ggplot2::ggplot(
        data = ., ggplot2::aes(x = LOCAL, na.rm = FALSE)) +
        ggplot2::geom_histogram(bins = 30) +
        ggplot2::labs(x = "Minor Allele Frequency (MAF)", y = "Number of markers") +
        ggplot2::expand_limits(y = 0) +
        ggplot2::theme(
          axis.title.x = ggplot2::element_text(size = 12, face = "bold"),
          axis.title.y = ggplot2::element_text(size = 12, face = "bold"),
          legend.position = "none",
          axis.text.x = ggplot2::element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5),
          strip.text.y = ggplot2::element_text(angle = 0, size = 12, face = "bold"),
          strip.text.x = ggplot2::element_text(size = 12, face = "bold")
        ) +
        ggplot2::facet_grid(POP_ID ~ GROUP, scales = "free_x")
      print(histo.maf.local)
      ggplot2::ggsave(
        filename = file.path(path.folder, "maf.local.spectrums.pdf"),
        plot = histo.maf.local, width = 15, height = pop.number * 5,
        dpi = 600, units = "cm", useDingbats = FALSE, limitsize = FALSE)
      message("Local maf spectrums written in the folder")
    }

    # import back data ---------------------------------------------------------
    input <- read_rad("maf.temp.rad")
    file.remove("maf.temp.rad")

    # Helper table for global and local MAF -------------------------------------
    if (interactive.filter) {
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
    }
    # Step 2. Thresholds selection -----------------------------------------------
    # maf.approach
    if (interactive.filter) {
      message("\nStep 2. Filtering markers based on the different MAF arguments\n")
      message("The maf.approach:\n
maf.approach = \"locus\" : looks at the frequency of all ALTernative allele on the locus.
maf.approach = \"SNP\" : SNPs on the same locus/read are considered independent.")
      # message("Choose the maf.approach (SNP/locus):")
      # maf.approach <- as.character(readLines(n = 1))

      maf.approach <- radiator_question(
        x = "\nChoose the maf.approach (SNP/locus): ", answer.opt = c("SNP", "locus"))


      # if (!maf.approach %in% c("SNP", "locus")) rlang::abort("maf.approach: SNP or locus")
    }

    if (maf.approach == "locus") {
      if (!tibble::has_name(input, "LOCUS") || !tibble::has_name(input, "POS")) {
        maf.approach <- "SNP"
        if (verbose) message("maf.approach = 'locus', switched to maf.approach = 'SNP',
because LOCUS and POS (SNP) info is not available")
        if (verbose) message("\nFor haplotype data, this is not important: it's all about locus/haplotypes")
      }
    }

    # maf.thresholds
    if (interactive.filter) {
      message("\nChoose the local minor allele threshold")
      message("    Using a FREQUENCY, choose a value between 0 and 0.9")
      message("    Using COUNT of alternate allele, choose an integer >= 1")
      # maf.local.threshold <- as.character(readLines(n = 1))
      maf.local.threshold <- radiator_question(
        x = "    Enter value: ", minmax = c(0, 1000))

      if (maf.local.threshold >= 1) {
        maf.local.threshold <- as.integer(maf.local.threshold)
        maf.count <- TRUE
      } else {
        maf.local.threshold <- as.numeric(maf.local.threshold)
        maf.count <- FALSE
      }
    }

    if (interactive.filter) {
      message("\nChoose the global minor allele threshold")
      message("    Using a FREQUENCY, choose a value between 0 and 0.9")
      message("    Using COUNT of alternate allele, choose an integer >= 1")
      message("    Note: please use the same method count/freqency as the local threshold")
      # maf.global.threshold <- as.character(readLines(n = 1))
      maf.global.threshold <- radiator_question(
        x = "    Enter value: ", minmax = c(0, 1000))
      if (maf.global.threshold >= 1) {
        maf.global.threshold <- as.integer(maf.global.threshold)
      } else {
        maf.global.threshold <- as.numeric(maf.global.threshold)
      }
    }

    # maf.operator
    if (interactive.filter) {
      message("\nThe maf.operator (AND/OR):
    Option 1: local \"AND\" global MAF thresholds are required to pass (more severe).
    Option 2: local \"OR\" global MAF thresholds are required to pass (more tolerant).")
      message("No idea what to choose? Use OR")
      # message("Choose the maf.operator (AND/OR):")
      # maf.operator <- as.character(readLines(n = 1))
      maf.operator <- radiator_question(
        x = "    Choose the maf.operator (AND/OR): ", answer.opt = c("AND", "OR"))

      # if (!maf.operator %in% c("OR", "AND")) rlang::abort("maf.operator: either OR/AND")
    }

    # maf.pop.num.threshold
    if (interactive.filter) {
      message("\nLast threshold... \n\nmaf.pop.num.threshold:
    How many populations are required to pass all the thresholds to keep the marker?\n
    Example: if you have 10 populations and choose maf.pop.num.threshold = 3,
    3 populations out of 10 are required to pass previous thresholds")
      message("    Note: not sure? use 1")
      # message("\n    Choose the maf.pop.num.threshold:")
      # maf.pop.num.threshold <- as.integer(readLines(n = 1))
      maf.pop.num.threshold <- radiator_question(
        x = "    Enter threshold (integer): ", minmax = c(1, 1000))
    }


    # Update the maf.data with pass or not filter based on threshold -------------
    if (interactive.filter) {
      message("Compiling all the statistics and thresholds")
      if (maf.count) {
        if (tibble::has_name(maf.data, "HAPLOTYPES")) {
          maf.data <- maf.data %>%
            dplyr::filter(REF == "MAF" ) %>%
            dplyr::mutate(
              OR = dplyr::if_else((ALT_LOCAL >= maf.local.threshold | ALT_GLOBAL >= maf.global.threshold), "pass", "pruned"),
              AND = dplyr::if_else((ALT_LOCAL >= maf.local.threshold & ALT_GLOBAL >= maf.global.threshold), "pass", "pruned")
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
            )

        } else {
          maf.data <- maf.data %>%
            dplyr::mutate(
              OR = dplyr::if_else((ALT_LOCAL >= maf.local.threshold | ALT_GLOBAL >= maf.global.threshold), "pass", "pruned"),
              AND = dplyr::if_else((ALT_LOCAL >= maf.local.threshold & ALT_GLOBAL >= maf.global.threshold), "pass", "pruned")
            ) %>%
            dplyr::group_by(MARKERS) %>%
            dplyr::mutate(
              OR = dplyr::if_else(length(POP_ID[OR == "pass"]) >= maf.pop.num.threshold, "pass", "pruned"),
              AND = dplyr::if_else(length(POP_ID[AND == "pass"]) >= maf.pop.num.threshold, "pass", "pruned")
            )
        }
      } else {
        if (tibble::has_name(maf.data, "HAPLOTYPES")) {
          maf.data <- maf.data %>%
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
            )

        } else {
          maf.data <- maf.data %>%
            dplyr::mutate(
              OR = dplyr::if_else((MAF_LOCAL >= maf.local.threshold | MAF_GLOBAL >= maf.global.threshold), "pass", "pruned"),
              AND = dplyr::if_else((MAF_LOCAL >= maf.local.threshold & MAF_GLOBAL >= maf.global.threshold), "pass", "pruned")
            ) %>%
            dplyr::group_by(MARKERS) %>%
            dplyr::mutate(
              OR = dplyr::if_else(length(POP_ID[OR == "pass"]) >= maf.pop.num.threshold, "pass", "pruned"),
              AND = dplyr::if_else(length(POP_ID[AND == "pass"]) >= maf.pop.num.threshold, "pass", "pruned")
            )
        }
      }

      if (ncol(markers.meta) > 1) {
        markers.meta.col <- colnames(markers.meta)
        # maf.data.thresholds <- suppressWarnings(
        #   maf.data.thresholds %>%
        #     dplyr::left_join(markers.meta, by = "MARKERS") %>%
        #     dplyr::select(dplyr::one_of(markers.meta.col), dplyr::everything(.)))

        maf.data <- suppressWarnings(
          maf.data %>%
            dplyr::left_join(markers.meta, by = "MARKERS") %>%
            dplyr::select(dplyr::one_of(markers.meta.col), dplyr::everything(.)))
      }

      readr::write_tsv(
        x = maf.data,
        path = file.path(path.folder, "maf.data.tsv"),
        col_names = TRUE,
        append = FALSE
      )
      message("\nWritten in the folder: maf.data.tsv\n")
    }

    # Filtering ------------------------------------------------------------------
    # To avoid changing or duplicating code depending on user choosing count or frequency:
    if (maf.count) {
      maf.data <- dplyr::select(maf.data, -MAF_LOCAL, -MAF_GLOBAL) %>%
        dplyr::rename(MAF_LOCAL = ALT_LOCAL, MAF_GLOBAL = ALT_GLOBAL)
    } else {
      maf.data <- dplyr::select(maf.data, -ALT_LOCAL, -ALT_GLOBAL)
    }


    if (maf.approach == "locus") {

      # Needed switched for haplotype file (e.g. stacks haplotype)
      if (tibble::has_name(maf.data, "MARKERS") && !tibble::has_name(maf.data, "LOCUS")) {
        maf.data <- dplyr::rename(maf.data, LOCUS = MARKERS)
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
      input <- NULL
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
      input <- NULL
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
            write_rad(data = ., path = file.path(path.folder, "temp.rad"))
            # readr::write_tsv(x = ., path = file.path(path.folder, "temp.rad"))

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
            calibrate_alleles(data = ., biallelic = FALSE, parallel.core = parallel.core)

          filter <- filter$input %>%
            dplyr::select(MARKERS, POP_ID, INDIVIDUALS, GT_VCF_NUC, REF, ALT, GT, GT_VCF) %>%
            dplyr::bind_rows(read_rad(data = file.path(path.folder, "temp.rad")))
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
    snp.before <- as.integer(dplyr::n_distinct(markers.meta$MARKERS))
    snp.after <- as.integer(dplyr::n_distinct(filter$MARKERS))
    snp.blacklist <- snp.before - snp.after

    if (tibble::has_name(markers.meta, "LOCUS")) {
      locus.before <- as.integer(dplyr::n_distinct(markers.meta$LOCUS))
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
      if (verbose) message("Writing the MAF filtered tidy data set: ", tidy.name)
      write_rad(data = filter, path = file.path(path.folder, tidy.name))
      # readr::write_tsv(x = filter, path = tidy.name)
    }

    # saving whitelist
    if (verbose) message("Writing the whitelist of markers: whitelist.markers.maf.tsv")

    if (tibble::has_name(filter, "CHROM")) {
      whitelist.markers <- dplyr::ungroup(filter) %>%
        dplyr::distinct(CHROM, LOCUS, POS) %>%
        dplyr::mutate_all(.tbl = ., .funs = as.character)
    } else {
      whitelist.markers <- dplyr::ungroup(filter) %>%
        dplyr::distinct(MARKERS) %>%
        dplyr::mutate_all(.tbl = ., .funs = as.character)
    }
    readr::write_tsv(whitelist.markers, file.path(path.folder, "whitelist.markers.maf.tsv"), append = FALSE, col_names = TRUE)


    # saving blacklist
    if (verbose) message("\nWriting the blacklist of markers: blacklist.markers.maf.tsv")
    if (tibble::has_name(filter, "CHROM")) {
      blacklist.markers <- dplyr::setdiff(
        dplyr::select(markers.meta, CHROM, LOCUS, POS), whitelist.markers)
      # dplyr::select(markers.meta, CHROM, LOCUS, POS), dplyr::mutate_all(.tbl = whitelist.markers, .funs = as.character))

    } else {
      blacklist.markers <- dplyr::setdiff(
        dplyr::select(markers.meta, MARKERS), whitelist.markers)
    }
    readr::write_tsv(blacklist.markers, file.path(path.folder, "blacklist.markers.maf.tsv"), append = FALSE, col_names = TRUE)

    # results --------------------------------------------------------------------
    if (verbose) cat("############################### RESULTS ###############################\n")
    if (verbose) message("maf.approach: ", maf.approach)
    if (verbose) message("maf.thresholds: ", "local = ", maf.local.threshold, ", global = ", maf.global.threshold)
    if (verbose) message("maf.operator: ", maf.operator)
    if (verbose) message("maf.pop.num.threshold: ", maf.pop.num.threshold)
    n.markers <- stringi::stri_join(markers.before, markers.blacklist, markers.after, sep = "; ")
    message("Number of markers/locus, before; blacklisted; after: ", n.markers)
    if (verbose) message("\nComputation time: ", round((proc.time() - timing)[[3]]), " sec")
    if (verbose) cat("############################## completed ##############################\n")
    res <- list()
    res$tidy.filtered.maf <- filter
    res$whitelist.markers <- whitelist.markers
    res$blacklist.markers <- blacklist.markers
    if (interactive.filter) {
      res$maf.data <- maf.data
      res$maf.helper.table <- maf.helper.table
      res$distribution.maf.global <- histo.maf.global
      res$distribution.maf.local <- histo.maf.local
      res$maf.global.summary <- maf.global.summary
    }
    res$filters.parameters <- filters.parameters
    return(res)
  }
}#End maf_local
