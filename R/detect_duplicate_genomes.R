# Detect duplicate genomes

#' @name detect_duplicate_genomes
#' @title Compute pairwise genome similarity or distance between individuals
#' to highligh potential duplicate individuals
#' @description The function can compute two methods
#' to highligh potential duplicate individuals.
#' \enumerate{
#' \item distance between individuals and/or
#' \item pairwise genome similarity
#' }

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.

#' @param subsample.markers (optional, integer) To speed up computation and rapidly
#' test the function's arguments (e.g. using 200 markers).
#' Default: \code{subsample.markers = NULL}.

#' @param random.seed (integer, optional) For reproducibility, set an integer
#' for randomness when argument \code{subsample.markers} is used.
#' By default, a random number is generated and printed.
#' Default: \code{random.seed = NULL}.

#' @param distance.method (character) The distance measure used inside \code{stats::dist}
#' (<= 30000 markers) or \code{amap::Dist} (> 30000 markers).
#' This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary".
#' Using \code{distance.method = NULL} will not run this method.
#' Default: \code{distance.method = "manhattan"}. This is very fast
#' compared to the genome similarity method. It uses allele counts and the codes
#' are tailored for biallelic and multiallelic markers.

#' @param genome (logical) Computes pairwise genome similarity in parallel.
#' The proportion of the shared genotypes is averaged across shared markers between
#' each pairwise comparison. This method makes filtering easier because the
#' threshold is more intuitive with the plots produced, but it's much longer
#' to run, even in parallel, so better to run overnight.
#' Default: \code{genome = FALSE}.

#' @param parallel.core (optional) The number of core for parallel computation.
#' Default: \code{parallel.core = parallel::detectCores() - 1}.

#' @return A list with potentially 8 objects:
#' \code{$distance }: results of the distance method
#' \code{$distance.stats}: Summary statistics of the distance method
#' \code{$pairwise.genome.similarity}: results of the genome method
#' \code{$genome.stats}: Summary statistics of the genome method
#' \code{$violin.plot.distance}: violin plot showing the distribution of pairwise distances
#' \code{$manhattan.plot.distance}: same info different visual with manhattan plot
#' \code{$violin.plot.genome}: violin plot showing the distribution of pairwise genome similarities
#' \code{$manhattan.plot.genome}: same info different visual with manhattan plot
#' \code{$blacklist.id.similar}: blacklisted duplicates
#'
#' Saved in the working directory:
#' individuals.pairwise.dist.tsv, individuals.pairwise.distance.stats.tsv,
#' individuals.pairwise.genome.similarity.tsv, individuals.pairwise.genome.stats.tsv,
#' blackliste.id.similar.tsv

#' @details
#' Strategically, run the default first (\code{distance.method},
#' no \code{genome})
#'
#' \strong{\code{distance.method} argument is fast, but...}
#'
#' you don't know if the observed comparison (close or distant)
#' is influenced by missing values/the number of markers in common
#' between the pair compared. This is something that needs to be considered.
#' Be suspicious of a \emph{distant outlier} from the same pop pairwise comparison,
#' and similarly, be suspicious of a \emph{close outlier} from a different pop
#' pairwise comparisons.
#'
#' If there is no outlier, don't bother running the function again with
#' (\code{genome = TRUE}).
#'
#'
#' \strong{\code{genome = TRUE}}
#'
#' The function will run slower, but...
#' If you see outliers with the first run, take the time to run the function
#' with \code{genome = TRUE}. Because this option is much better at detecting
#' duplicated individuals and it also shows the impact of \strong{missingness}
#' or the number of \strong{shared markers} between comparisons.
#'
#' \emph{Your outlier duo could well be the result of one of the individual having
#' an extremely low number genotypes...}


#' @export
#' @rdname detect_duplicate_genomes
#' @importFrom stringi stri_paste stri_replace_all_fixed
#' @importFrom dplyr arrange rename select group_by filter mutate rename_ filter_ bind_cols bind_rows summarise n_distinct intersect desc
#' @importFrom utils combn
#' @importFrom stats na.omit var median quantile dist
#' @importFrom amap Dist
#' @importFrom readr write_tsv
#' @importFrom parallel detectCores
#' @importFrom purrr flatten map
#' @importFrom tibble as_data_frame has_name remove_rownames column_to_rownames
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot stat_summary labs theme element_blank element_text geom_jitter scale_colour_manual scale_y_reverse theme_light ggsave


#' @examples
#' \dontrun{
#' # First run and simplest way (if you have the tidy df):
#' dup <- radiator::detect_duplicate_genomes(data = "wombat_tidy.tsv")
#'
#' #If you need a tidy df:
#' dup <- radiator::tidy_genomic_data(
#' data = "wombat_tidy.tsv",
#' strata = "wombat.strata.tsv",
#' vcf.metadata = FALSE
#' ) %>%
#' radiator::detect_duplicate_genomes(data = .)
#'
#' # This will use by defaul:
#' distance.method = "manhattan"
#' genome = FALSE
#' #parallel.core = all my CPUs - 1
#'
#' # To view the manhattan plot:
#' dup$manhattan.plot.distance
#'
#' # to view the data stats
#' dup.data.stats <- dup$distance.stats
#'
#' # to view the data
#' dup.data <- dup$distance
#'
#' # Based on the look of the distribution using both manhattan and boxplot,
#' # I can filter the dataset to highlight potential duplicates.
#'
#' # To run the distance (with euclidean distance instead of the default manhattan,
#' # and also carry the second analysis (with the genome method):
#' dup <- radiator::tidy_genomic_data(
#' data = "wombat_tidy.tsv",
#' strata = "wombat.strata.tsv",
#' vcf.metadata = FALSE
#' ) %>%
#' radiator::detect_duplicate_genomes(data = ., distance.method = "euclidean", genome = TRUE)
#'
#' # to view the data of the genome data
#' dup.data <- dup$pairwise.genome.similarity
#'
#' # Based on the look of the distribution using both manhattan and boxplot,
#' # I can filter the dataset based on 98% of identical genotype proportion,
#' # to highlight potential duplicates:
#' dup.filtered <- dplyr::filter(.data = dup.data, PROP_IDENTICAL > 0.98)
#'
#' # Get the list of duplicates id
#' dup.list.names <- data.frame(INDIVIDUALS = unique(c(dup.filtered$ID1, dup.filtered$ID2)))
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

detect_duplicate_genomes <- function(
  data,
  subsample.markers = NULL,
  random.seed = NULL,
  distance.method = "manhattan",
  genome = FALSE,
  parallel.core = parallel::detectCores() - 1
) {
  cat("\n")
  cat("###############################################################################\n")
  cat("##################### radiator::detect_duplicate_genomes ######################\n")
  cat("###############################################################################\n")
  timing <- proc.time()
  opt.change <- getOption("width")
  options(width = 70)
  # Manage missing arguments ---------------------------------------------------
  if (missing(data)) stop("missing data argument")

  # folder ---------------------------------------------------------------------
  # Get date and time to have unique filenaming
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  folder.extension <- stringi::stri_join("detect_duplicate_genomes_", file.date, sep = "")
  path.folder <- stringi::stri_join(getwd(),"/", folder.extension, sep = "")
  dir.create(file.path(path.folder))
  message("Folder created: \n", folder.extension)
  file.date <- NULL #unused object

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    data <- radiator::tidy_wide(data = data, import.metadata = TRUE)
  }
  want <- c("MARKERS", "CHROM", "LOCUS", "POS", "POP_ID", "INDIVIDUALS", "GT", "GT_BIN", "REF", "ALT")
  data <- suppressWarnings(dplyr::select(data, dplyr::one_of(want)))

  # necessary steps to make sure we work with unique markers and not duplicated LOCUS
  if (tibble::has_name(data, "LOCUS") && !tibble::has_name(data, "MARKERS")) {
    data <- dplyr::rename(.data = data, MARKERS = LOCUS)
  }

  # Subsampling ----------------------------------------------------------------
  if (!is.null(subsample.markers)) {
    # Set seed for random sampling
    if (is.null(random.seed)) {
      random.seed <- sample(x = 1:1000000, size = 1)
      set.seed(random.seed)
      message("Random seed used: ", random.seed)
    } else {
      set.seed(random.seed)
    }
    sample.markers <- dplyr::distinct(data, MARKERS) %>%
      dplyr::sample_n(tbl = ., size = subsample.markers) %>%
      readr::write_tsv(x = ., path = file.path(path.folder, stringi::stri_join("subsampled.markers_random.seed_", random.seed, ".tsv"))) %>%
      purrr::flatten_chr(.)
    data <- dplyr::filter(data, MARKERS %in% sample.markers)
    sample.markers <- NULL
  }

  # strata
  strata <- dplyr::ungroup(data) %>%
    dplyr::distinct(POP_ID, INDIVIDUALS)

  # New list to prepare for results
  res <- list()

  # Preparing data for comparisons ---------------------------------------------
  message("Preparing data for analysis")

  #Genotyped stats -------------------------------------------------------------
  n.markers <- dplyr::n_distinct(data$MARKERS)
  geno.stats <- data %>%
    dplyr::filter(GT != "000000") %>%
    dplyr::group_by(INDIVIDUALS) %>%
    dplyr::summarise(GENOTYPED_PROP = length(GT) / n.markers) %>%
    readr::write_tsv(
      x = .,
      path = file.path(path.folder, "genotyped.statistics.tsv"))

  # GT_BIN available
  if (!is.null(distance.method) & tibble::has_name(data, "GT_BIN")) {
    input.prep <- dplyr::ungroup(data) %>%
      dplyr::select(MARKERS, INDIVIDUALS, ALT = GT_BIN) %>%
      dplyr::mutate(REF = 2 - ALT) %>%
      tidyr::gather(data = ., key = ALLELES, value = n, -c(MARKERS,INDIVIDUALS)) %>%
      dplyr::mutate(MARKERS_ALLELES = stringi::stri_join(MARKERS, ALLELES, sep = ".")) %>%
      dplyr::select(-ALLELES) %>%
      dplyr::arrange(MARKERS_ALLELES,INDIVIDUALS)
  }

  # GT_BIN NOT available
  if (!tibble::has_name(data, "GT_BIN")) {
    # Allele count
    missing.geno <- dplyr::select(.data = data, MARKERS, INDIVIDUALS, GT) %>%
      dplyr::filter(GT == "000000") %>%
      dplyr::select(-GT) %>%
      dplyr::mutate(MISSING = rep("blacklist", n()))

    message("Preparing data: calculating allele count")
    input.prep <- dplyr::select(data, MARKERS, INDIVIDUALS, GT) %>%
      dplyr::left_join(
        dplyr::distinct(data, MARKERS) %>%
          dplyr::mutate(
            SPLIT_VEC = dplyr::ntile(x = 1:nrow(.), n = parallel.core * 3))
        , by = "MARKERS") %>%
      split(x = ., f = .$SPLIT_VEC) %>%
      .radiator_parallel(
        X = .,
        FUN = allele_count,
        mc.cores = parallel.core
      ) %>%
      dplyr::bind_rows(.)

    if (nrow(missing.geno) > 0) {
      input.prep <- dplyr::anti_join(input.prep, missing.geno, by = c("MARKERS", "INDIVIDUALS"))
    }

    input.prep <- input.prep %>%
      dplyr::mutate(MARKERS_ALLELES = stringi::stri_join(MARKERS, GT, sep = ".")) %>%
      dplyr::select(-GT) %>%
      dplyr::arrange(MARKERS_ALLELES, INDIVIDUALS)

    missing.geno <- NULL # unused object
  }#end preparing data

  # Computing distance ---------------------------------------------------------
  if (!is.null(distance.method)) {
    message("Calculating ", distance.method, " distances between individuals...")

    res$distance <- distance_individuals(
      x = dplyr::select(input.prep, -MARKERS),
      strata = strata,
      distance.method = distance.method,
      parallel.core = parallel.core
    ) %>%
      dplyr::mutate(PAIRS = seq(from = 1, to = n(), by = 1)) %>%
      dplyr::arrange(PAIRS)


    geno <- dplyr::select(res$distance, ID1, ID2, PAIRS) %>%
      dplyr::left_join(dplyr::rename(geno.stats, ID1 = INDIVIDUALS, ID1_G = GENOTYPED_PROP), by = "ID1") %>%
      dplyr::left_join(dplyr::rename(geno.stats, ID2 = INDIVIDUALS, ID2_G = GENOTYPED_PROP), by = "ID2") %>%
      tidyr::gather(data = ., key = GENOTYPED_MAX, value = GENOTYPED_PROP, -c(ID1, ID2, PAIRS)) %>%
      dplyr::group_by(PAIRS) %>%
      dplyr::filter(GENOTYPED_PROP == max(GENOTYPED_PROP)) %>%
      dplyr::ungroup(.) %>%
      dplyr::distinct(PAIRS, .keep_all = TRUE) %>%
      dplyr::select(-GENOTYPED_PROP) %>%
      dplyr::mutate(GENOTYPED_MAX = dplyr::if_else(GENOTYPED_MAX == "ID1_G", ID1, ID2)) %>%
      dplyr::arrange(PAIRS) %>%
      dplyr::select(-ID1, -ID2, -PAIRS)

    res$distance <- dplyr::bind_cols(res$distance, geno)
    geno <- NULL
    # test <- res$distance

    readr::write_tsv(
      x = res$distance,
      path = file.path(path.folder, "individuals.pairwise.dist.tsv"),
      col_names = TRUE
    )

    # Stats
    message("Generating summary statistics")
    res$distance.stats <- res$distance %>%
      dplyr::summarise(
        MEAN = mean(DISTANCE_RELATIVE, na.rm = TRUE),
        MEDIAN = stats::median(DISTANCE_RELATIVE, na.rm = TRUE),
        SE = round(sqrt(stats::var(DISTANCE_RELATIVE, na.rm = TRUE)/length(stats::na.omit(DISTANCE_RELATIVE))), 2),
        MIN = round(min(DISTANCE_RELATIVE, na.rm = TRUE), 2),
        MAX = round(max(DISTANCE_RELATIVE, na.rm = TRUE), 2),
        QUANTILE25 = stats::quantile(DISTANCE_RELATIVE, 0.25), # quantile25
        QUANTILE75 = stats::quantile(DISTANCE_RELATIVE, 0.75)#, # quantile75
        # OUTLIERS_LOW = QUANTILE25 - (1.5 * (QUANTILE75 - QUANTILE25)), # outliers : below the outlier boxplot
        # OUTLIERS_HIGH = QUANTILE75 + (1.5 * (QUANTILE75 - QUANTILE25)) # outliers : higher the outlier boxplot
      ) %>%
      readr::write_tsv(
        x = .,
        path = stringi::stri_join(path.folder, "/individuals.pairwise.distance.stats.tsv"),
        col_names = TRUE
      )

    message("Generating plots")
    # violin plot
    res$violin.plot.distance <- ggplot2::ggplot(
      data = res$distance,
      ggplot2::aes(x = PAIRWISE, y = DISTANCE_RELATIVE, na.rm = TRUE)
    ) +
      ggplot2::geom_violin(trim = TRUE) +
      ggplot2::geom_boxplot(width = 0.1, fill = "black", outlier.colour = "black") +
      ggplot2::stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white") +
      ggplot2::labs(y = "Distance (relative)\n <- distant      close->") +
      ggplot2::labs(x = "Pairwise comparisons") +
      ggplot2::scale_y_reverse() +
      ggplot2::theme(
        # legend.position = "none",
        panel.grid.minor.x = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        # panel.grid.major.y = element_blank(),
        axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.y = ggplot2::element_text(size = 8, family = "Helvetica")
      )
    ggplot2::ggsave(
      filename = file.path(path.folder, "violin.plot.distance.pdf"),
      plot = res$violin.plot.distance,
      width = 20, height = 15, dpi = 600, units = "cm", useDingbats = FALSE)

    # Manhattan plot
    res$manhattan.plot.distance <- ggplot2::ggplot(
      data = res$distance,
      ggplot2::aes(x = PAIRWISE, y = DISTANCE_RELATIVE, colour = POP_COMP)
    ) +
      ggplot2::geom_jitter(alpha = 0.3) +
      ggplot2::labs(y = "Distance (relative)\n <- distant      close->") +
      ggplot2::labs(x = "Pairwise comparisons") +
      ggplot2::labs(colour = "Population comparisons") +
      ggplot2::scale_colour_manual(values = c("#0571b0", "black")) +
      ggplot2::scale_y_reverse() +
      ggplot2::theme_light() +
      ggplot2::theme(
        # legend.position = "none",
        panel.grid.minor.x = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        # panel.grid.major.y = element_blank(),
        axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.y = ggplot2::element_text(size = 8, family = "Helvetica")
      )
    ggplot2::ggsave(
      filename = file.path(path.folder, "manhattan.plot.distance.pdf"),
      plot = res$manhattan.plot.distance,
      width = 20, height = 15, dpi = 600, units = "cm", useDingbats = FALSE)
  } # end distance method

  # Compute genome similarity -------------------------------------------------
  if (genome) {

    # If GT_BIN available, we need a new input.prep (not the same as dist method)
    if (tibble::has_name(data, "GT_BIN")) {
      input.prep <- dplyr::filter(.data = data, !is.na(GT_BIN))
    }

    data <- NULL

    # all combination of individual pair
    id.pairwise <- utils::combn(unique(input.prep$INDIVIDUALS), 2, simplify = FALSE)

    # get the number of pairwise comp.
    number.pairwise <- length(id.pairwise)

    # Optimizing cpu usage
    if (number.pairwise <= 50) {
      round.cpu <- floor(number.pairwise / parallel.core)
    } else {
      round.cpu <- floor(number.pairwise / (50 * parallel.core))
    }
    # as.integer is usually twice as light as numeric vector...
    split.vec <- as.integer(floor((parallel.core * round.cpu * (1:number.pairwise - 1) / number.pairwise) + 1))
    id.pairwise <- split(x = id.pairwise, f = split.vec)

    message("Starting scan for duplicate genomes")
    message("Pairwise comparisons: ", number.pairwise)
    if (number.pairwise > 5000) message("    Time for coffee...")
    round.cpu <- split.vec <- number.pairwise <- NULL

    pairwise.genome.similarity <- list()
    pairwise.genome.similarity <- .radiator_parallel(
      X = id.pairwise,
      FUN = genome_similarity,
      mc.preschedule = FALSE,
      mc.silent = FALSE,
      mc.cleanup = TRUE,
      mc.cores = parallel.core,
      input.prep = input.prep
    ) %>%
      dplyr::bind_rows(.)

    input.prep <- id.pairwise <- NULL # no longer needed

    # Include population info with strata
    ID1.pop <- suppressWarnings(
      pairwise.genome.similarity %>%
        dplyr::select(INDIVIDUALS = ID1) %>%
        dplyr::inner_join(strata, by = "INDIVIDUALS") %>%
        dplyr::select(ID1_POP = POP_ID))

    ID2.pop <- suppressWarnings(
      pairwise.genome.similarity %>%
        dplyr::select(INDIVIDUALS = ID2) %>%
        dplyr::inner_join(strata, by = "INDIVIDUALS") %>%
        dplyr::select(ID2_POP = POP_ID))

    pairwise.genome.similarity <- dplyr::bind_cols(
      pairwise.genome.similarity, ID1.pop, ID2.pop
    ) %>%
      dplyr::mutate(
        POP_COMP = ifelse(ID1.pop == ID2.pop, "same pop", "different pop"),
        POP_COMP = factor(POP_COMP, levels = c("same pop", "different pop"), ordered = TRUE),
        PAIRWISE = rep("pairwise comparison", n()),
        METHOD = rep("genome similarity", n())
      )
    ID1.pop <- ID2.pop <- NULL
    res$pairwise.genome.similarity <- pairwise.genome.similarity %>%
      dplyr::arrange(dplyr::desc(PROP_IDENTICAL)) %>%
      dplyr::mutate(PAIRS = seq(from = 1, to = n(), by = 1)) %>%
      dplyr::arrange(PAIRS)
    pairwise.genome.similarity <- NULL

    geno <- dplyr::select(res$pairwise.genome.similarity, ID1, ID2, PAIRS) %>%
      dplyr::left_join(dplyr::rename(geno.stats, ID1 = INDIVIDUALS, ID1_G = GENOTYPED_PROP), by = "ID1") %>%
      dplyr::left_join(dplyr::rename(geno.stats, ID2 = INDIVIDUALS, ID2_G = GENOTYPED_PROP), by = "ID2") %>%
      tidyr::gather(data = ., key = GENOTYPED_MAX, value = GENOTYPED_PROP, -c(ID1, ID2, PAIRS)) %>%
      dplyr::group_by(PAIRS) %>%
      dplyr::filter(GENOTYPED_PROP == max(GENOTYPED_PROP)) %>%
      dplyr::ungroup(.) %>%
      dplyr::distinct(PAIRS, .keep_all = TRUE) %>%
      dplyr::select(-GENOTYPED_PROP) %>%
      dplyr::mutate(GENOTYPED_MAX = dplyr::if_else(GENOTYPED_MAX == "ID1_G", ID1, ID2)) %>%
      dplyr::arrange(PAIRS) %>%
      dplyr::select(-ID1, -ID2, -PAIRS)

    res$pairwise.genome.similarity <- dplyr::bind_cols(res$pairwise.genome.similarity, geno)
    geno <- NULL

    readr::write_tsv(
      x = res$pairwise.genome.similarity,
      path = file.path(path.folder, "individuals.pairwise.genome.similarity.tsv"),
      col_names = TRUE
    )

    # Stats
    message("Generating summary statistics")
    res$genome.stats <- res$pairwise.genome.similarity %>%
      dplyr::summarise(
        MEAN = mean(PROP_IDENTICAL, na.rm = TRUE),
        MEDIAN = stats::median(PROP_IDENTICAL, na.rm = TRUE),
        SE = round(sqrt(stats::var(PROP_IDENTICAL, na.rm = TRUE)/length(stats::na.omit(PROP_IDENTICAL))), 2),
        MIN = round(min(PROP_IDENTICAL, na.rm = TRUE), 2),
        MAX = round(max(PROP_IDENTICAL, na.rm = TRUE), 2),
        QUANTILE25 = stats::quantile(PROP_IDENTICAL, 0.25), # quantile25
        QUANTILE75 = stats::quantile(PROP_IDENTICAL, 0.75)#, # quantile75
        # OUTLIERS_LOW = QUANTILE25 - (1.5 * (QUANTILE75 - QUANTILE25)), # outliers : below the outlier boxplot
        # OUTLIERS_HIGH = QUANTILE75 + (1.5 * (QUANTILE75 - QUANTILE25)) # outliers : higher the outlier boxplot
      )
    readr::write_tsv(
      x = res$genome.stats,
      path = file.path(path.folder, "individuals.pairwise.genome.stats.tsv"),
      col_names = TRUE
    )

    # Visualization ------------------------------------------------------------
    message("Generating the plots")

    # violin plot
    res$violin.plot.genome <- ggplot2::ggplot(
      data = res$pairwise.genome.similarity,
      ggplot2::aes(x = PAIRWISE, y = PROP_IDENTICAL, na.rm = TRUE)
    ) +
      ggplot2::geom_violin(trim = TRUE) +
      ggplot2::geom_boxplot(width = 0.1, fill = "black", outlier.colour = "black") +
      ggplot2::stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white") +
      ggplot2::labs(y = "Genome similarity (proportion)") +
      ggplot2::labs(x = "Pairwise comparison") +
      ggplot2::theme(
        # legend.position = "none",
        panel.grid.minor.x = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        # panel.grid.major.y = element_blank(),
        axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.y = ggplot2::element_text(size = 8, family = "Helvetica")
      )
    ggplot2::ggsave(
      filename = file.path(path.folder, "violin.plot.genome.pdf"),
      plot = res$violin.plot.genome,
      width = 20, height = 15, dpi = 600, units = "cm", useDingbats = FALSE)

    # Manhattan plot
    res$manhattan.plot.genome <- ggplot2::ggplot(
      data = res$pairwise.genome.similarity,
      ggplot2::aes(x = PAIRWISE, y = PROP_IDENTICAL, colour = POP_COMP, size = MARKERS_COMMON)
    ) +
      ggplot2::geom_jitter(alpha = 0.3) +
      ggplot2::labs(y = "Genome similarity (proportion)") +
      ggplot2::labs(x = "Pairwise comparisons") +
      ggplot2::labs(colour = "Population comparisons") +
      ggplot2::scale_colour_manual(values = c("#0571b0", "black")) +
      ggplot2::scale_size_area(name = "Markers in common", max_size = 5) +
      ggplot2::theme_light() +
      ggplot2::theme(
        # legend.position = "none",
        panel.grid.minor.x = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        # panel.grid.major.y = element_blank(),
        axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.y = ggplot2::element_text(size = 8, family = "Helvetica")
      )

    ggplot2::ggsave(
      filename = file.path(path.folder, "manhattan.plot.genome.pdf"),
      plot = res$manhattan.plot.genome,
      width = 20, height = 15, dpi = 600, units = "cm", useDingbats = FALSE)
  } # end genome method

  # Removing duplicates
  message("\nInspect tables and figures to decide if some individual(s) need to be blacklisted")
  message("blacklist individual(s) (y/n): ")
  remove.id <- as.character(readLines(n = 1))
  if (remove.id == "y") {
    message("2 options to blacklist individuals: manually or with a threshold")
    message("    manually: the function generate a blacklist that you populate manually")
    message("    threshold: more powerful to fully remove duplicates")
    message("\n    remove (manually/threshold): ")
    remove.dup <- as.character(readLines(n = 1))
    if (remove.dup == "manually") {
      readr::write_tsv(
        x = tibble::data_frame(INDIVIDUALS = as.character()),
        path = file.path(path.folder, "blacklist.id.similar.tsv"),
        append = FALSE, col_names = TRUE)
      message("    An empty blacklist file was generated: blacklist.id.similar.tsv")
      message("    Keep column name, just add the individual(s) to blacklist(s)")
      res$blacklist.id.similar <- "check blacklist.id.similar.tsv file"
    } else {
      message("Use the distance or genome analysis to blacklist duplicates ? (distance/genome): ")
      analysis <- as.character(readLines(n = 1))
      if (analysis == "distance") {
        data <-  "individuals.pairwise.dist.tsv"
      } else {
        data <-  "individuals.pairwise.genome.similarity.tsv"
      }
      message("\n    Enter threshold to remove duplicates: ")
      dup.threshold <- as.numeric(readLines(n = 1))

      message("\n    Remove all duplicates involved in pairs from different pop/group?")
      message("\n    y: remove both samples in the pair")
      message("\n    n: removes one sample in the pair, with more missing genotypes: (y/n)")
      diff.pop.remove <- as.character(readLines(n = 1))
      if (diff.pop.remove == "n") {
        diff.pop.remove <- FALSE
      } else {
        diff.pop.remove <- TRUE
      }
      old.dir <- getwd()
      setwd(path.folder)
      duplicates <- remove_duplicates(
        data = data,
        stats = "genotyped.statistics.tsv",
        dup.threshold = dup.threshold,
        diff.pop.remove = diff.pop.remove)
      setwd(old.dir)
      res$blacklist.id.similar <- duplicates$blacklist.id
    }
  }



  # RESULTS --------------------------------------------------------------------
  cat("################################### RESULTS ###################################\n")
  message("Object in the list (if all arguments are selected):\n
$distance                         # Distance method results
$distance.stats                   # Summary statistics of the distance method
$pairwise.genome.similarity       # Genome method results
$genome.stats                     # Summary statistics of the genome method
$blacklist.id.similar             # Blacklisted duplicates\n\n
Visualization:
    $violin.plot.distance
    $manhattan.plot.distance
    $violin.plot.genome
    $manhattan.plot.genome\n
Saved in the working directory:
    blacklist.id.similar.tsv
    individuals.pairwise.dist.tsv
    individuals.pairwise.distance.stats.tsv
    individuals.pairwise.genome.similarity.tsv
    individuals.pairwise.genome.stats.tsv
    violin.plot.distance.pdf
    manhattan.plot.distance.pdf
    violin.plot.genome.pdf
    manhattan.plot.genome.pdf
")
  message("More details in: ", folder.extension)
  message("Computation time: ", round((proc.time() - timing)[[3]]), " sec")
  cat("############################## completed ##############################\n")
  options(width = opt.change)
  return(res)
} # end function detect_duplicate_genomes

# Internal nested functions: ---------------------------------------------------

# distance method --------------------------------------------------------------
#' @title Distance individuals
#' @description distance method
#' @rdname distance_individuals
#' @export
#' @keywords internal
distance_individuals <- function(
  x,
  strata = NULL,
  distance.method = "manhattan",
  parallel.core = parallel::detectCores() - 1
) {
  # Prep data
  dist.computation <- suppressWarnings(
    dplyr::ungroup(x) %>%
      dplyr::group_by(INDIVIDUALS) %>%
      tidyr::spread(data = ., key = MARKERS_ALLELES, value = n) %>%
      dplyr::ungroup(.) %>%
      tibble::as_data_frame(.) %>%
      tibble::remove_rownames(.) %>%
      tibble::column_to_rownames(df = ., var = "INDIVIDUALS")
  )
  # rownames(dist.computation) <- dist.computation[["INDIVIDUALS"]]
  # dist.computation[["INDIVIDUALS"]] <- NULL

  # compute distance
  # gain in speed between the 2 is very small on small data set
  if (dplyr::n_distinct(x$MARKERS_ALLELES) > 60000) {
    dist.computation <- suppressWarnings(
      amap::Dist(
        x = dist.computation,
        method = distance.method,
        nbproc = parallel.core
      )
    )

  } else {
    dist.computation <- stats::dist(
      x = dist.computation,
      method = distance.method
    )
  }

  # melt the dist matrice into a data frame
  dist.computation <- distance2df(dist.computation)

  # Include population info with strata
  ID1.pop <- suppressWarnings(
    dplyr::select(.data = dist.computation, INDIVIDUALS = ID1) %>%
      dplyr::inner_join(strata, by = "INDIVIDUALS") %>%
      dplyr::select(ID1_POP = POP_ID))

  ID2.pop <- suppressWarnings(
    dplyr::select(.data = dist.computation, INDIVIDUALS = ID2) %>%
      dplyr::inner_join(strata, by = "INDIVIDUALS") %>%
      dplyr::select(ID2_POP = POP_ID))

  dist.computation <- dplyr::bind_cols(dist.computation, ID1.pop, ID2.pop) %>%
    dplyr::mutate(
      POP_COMP = ifelse(ID1.pop == ID2.pop, "same pop", "different pop"),
      POP_COMP = factor(POP_COMP, levels = c("same pop", "different pop"), ordered = TRUE),
      PAIRWISE = rep("pairwise", n()),
      METHOD = rep(distance.method, n())
    )
  x <- ID1.pop <- ID2.pop <- NULL
  return(dist.computation)
}#End distance_individuals


# pairwise genome similarity method---------------------------------------------
#' @title Pgenome_similarity_map
#' @description for purrr::map
#' @rdname genome_similarity_map
#' @export
#' @keywords internal
genome_similarity_map <- function(list.pair = NULL, input.prep = NULL) {

  id.select <- list.pair
  id1 <- id.select[1]
  id2 <- id.select[2]

  if (tibble::has_name(input.prep, "GT_BIN")) {
    # filtered dataset for the 2 ind.
    input.select <- dplyr::filter(
      .data = input.prep,
      !is.na(GT_BIN) & INDIVIDUALS %in% id.select
    ) %>%
      dplyr::select(MARKERS, INDIVIDUALS, GT_BIN)

    # genotypes & markers
    id1.data <- dplyr::filter(.data = input.select, INDIVIDUALS %in% id1) %>%
      dplyr::select(-INDIVIDUALS)

    id2.data <- dplyr::filter(.data = input.select, INDIVIDUALS %in% id2) %>%
      dplyr::select(-INDIVIDUALS)

    # markers
    id1.markers <- dplyr::select(id1.data, MARKERS)

    id2.markers <- dplyr::select(id2.data, MARKERS)

    # output comparison
    genome.comparison <- tibble::data_frame(
      ID1 = id1,
      ID2 = id2,
      MARKERS_COMMON = nrow(dplyr::intersect(x = id1.markers, y = id2.markers)),
      IDENTICAL = nrow(dplyr::intersect(x = id1.data, y = id2.data)),
      DIFFERENT = MARKERS_COMMON - IDENTICAL,
      PROP_IDENTICAL = IDENTICAL / MARKERS_COMMON
    )
  } else {
    input.select <- dplyr::filter(
      .data = input.prep,
      INDIVIDUALS %in% id.select & n != 0
    )

    # genotypes & markers
    id1.data <- dplyr::filter(.data = input.select, INDIVIDUALS %in% id1) %>%
      dplyr::select(-INDIVIDUALS)

    id2.data <- dplyr::filter(.data = input.select, INDIVIDUALS %in% id2) %>%
      dplyr::select(-INDIVIDUALS)

    # markers
    id1.markers <- dplyr::distinct(id1.data, MARKERS)
    id2.markers <- dplyr::distinct(id2.data, MARKERS)

    # output comparison
    genome.comparison <- tibble::data_frame(
      ID1 = id1,
      ID2 = id2,
      MARKERS_COMMON = nrow(dplyr::intersect(x = id1.markers, y = id2.markers)),
      IDENTICAL = nrow(dplyr::intersect(x = id1.data, y = id2.data) %>% dplyr::distinct(MARKERS)),
      DIFFERENT = MARKERS_COMMON - IDENTICAL,
      PROP_IDENTICAL = IDENTICAL / MARKERS_COMMON
    )
  }
  #unused objets:
  id.select <- id1 <- id2 <- input.select <- id1.data <- id2.data <- id1.markers <- id2.markers <- NULL
  return(genome.comparison)
}#End genome_similarity_map


#' @title genome_similarity
#' @description for the parallel part
#' @rdname genome_similarity
#' @export
#' @keywords internal
genome_similarity <- function(list.pair, input.prep = NULL, ...) {
  # small.list.pair <- purrr::flatten(id.pairwise.split[1]) #test


  genome.comparison <- purrr::map(.x = list.pair, .f = genome_similarity_map, input.prep = input.prep) %>%
    dplyr::bind_rows(.)

  return(genome.comparison)
} #End genome_similarity

# calculate allele count in parallel -------------------------------------------
#' @title allele_count
#' @description to calculate allele count in parallel
#' @rdname allele_count
#' @export
#' @keywords internal
allele_count <- function(x) {
  res <- dplyr::ungroup(x) %>%
    dplyr::select(MARKERS, INDIVIDUALS, GT) %>%
    dplyr::filter(GT != "000000") %>%
    dplyr::mutate(
      A1 = stringi::stri_sub(str = GT, from = 1, to = 3),
      A2 = stringi::stri_sub(str = GT, from = 4, to = 6)
    ) %>%
    dplyr::select(-GT) %>%
    tidyr::gather(
      data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS)) %>%
    dplyr::arrange(MARKERS, INDIVIDUALS, GT) %>%
    dplyr::count(x = ., INDIVIDUALS, MARKERS, GT) %>%
    dplyr::ungroup(.) %>%
    tidyr::complete(
      data = ., INDIVIDUALS, tidyr::nesting(MARKERS, GT), fill = list(n = 0))
  return(res)
}#End allele_count

# melt the dist matrice into a data frame --------------------------------------
#' @title distance2df
#' @description melt the dist matrice into a data frame
#' @rdname distance2df
#' @export
#' @keywords internal
distance2df <- function(x) {
  res <- as.matrix(x)
  diag(res) <- NA
  res[lower.tri(res)] <- NA
  res <- dplyr::bind_cols(tibble::data_frame(ID1 = rownames(res)),
                          tibble::as_data_frame(res)) %>%
    tidyr::gather(data = ., key = ID2, value = DISTANCE, -ID1) %>%
    dplyr::filter(!is.na(DISTANCE)) %>%
    dplyr::mutate(DISTANCE_RELATIVE = DISTANCE/max(DISTANCE)) %>%
    dplyr::arrange(DISTANCE)
  return(res)
}#End distance2df


# remove_duplicates  -----------------------------------------------------------

#' @name remove_duplicates
#' @title Read tidy genomic data file ending .rad
#' @description Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' To remove duplicate individuals based on threshold established from the
#' visualization figures.

#' @param data (path) The individual's pairwise data.
#' Default: \code{data = "individuals.pairwise.dist.tsv"}.
#' @param stats (path) The genotype statistics
#' Default: \code{stats = "genotyped.statistics.tsv"}.
#' @param dup.threshold (double) The threshold to filter out duplicates
#' Default: \code{dup.threshold = 0.25}.
#' @param diff.pop.remove Remove all individuals in pairs from different pop.
#' Both samples are potentially problems. With defautl, the function will not keep
#' one sample in the duplicate pair.
#' Default: \code{diff.pop.remove = TRUE}.

#' @return A list with blacklisted duplicates. Write the blacklist in the working
#' directory.
#' @export
#' @keywords internal
#' @rdname remove_duplicates
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

remove_duplicates <- function(
  data = "individuals.pairwise.dist.tsv",
  stats = "genotyped.statistics.tsv",
  dup.threshold = 0.25,
  diff.pop.remove = TRUE
) {
  dup.filtered <- suppressWarnings(suppressMessages(readr::read_tsv(data)))

  if (tibble::has_name(dup.filtered, "DISTANCE_RELATIVE")) {
    dup.filtered <- dup.filtered %>%
      dplyr::filter(DISTANCE_RELATIVE < dup.threshold)
  } else {
    dup.filtered <- dup.filtered %>%
      dplyr::filter(PROP_IDENTICAL > dup.threshold)
  }


  if (nrow(dup.filtered) > 0) {
    dup.list.names <- tibble::data_frame(INDIVIDUALS = c(dup.filtered$ID1, dup.filtered$ID2)) %>%
      dplyr::group_by(INDIVIDUALS) %>%
      dplyr::tally(.) %>%
      dplyr::ungroup(.) %>%
      dplyr::arrange(dplyr::desc(n)) %>%
      dplyr::distinct(INDIVIDUALS) %>%
      purrr::flatten_chr(.)

    geno.stats <- readr::read_tsv(stats, col_types = "cd")

    res <- list(blacklist.id = tibble::tibble(INDIVIDUALS = character(0)),
                whitelist.id = tibble::tibble(INDIVIDUALS = character(0)))

    for (i in dup.list.names) {
      # i <- dup.list.names[1]
      dups <- dplyr::filter(dup.filtered, ID1 %in% i | ID2 %in% i)
      dups <- sort(unique(c(dups$ID1, dups$ID2)))

      # find all duplicates associated with the network
      new.dups <- 0L
      while(length(new.dups) > 0) {
        new.dups <- dplyr::filter(dup.filtered, ID1 %in% dups | ID2 %in% dups)
        new.dups <- sort(unique(c(new.dups$ID1, new.dups$ID2)))
        new.dups <- purrr::keep(.x = new.dups, .p = !new.dups %in% dups)
        if (length(new.dups) > 0) {
          dups <- c(dups, new.dups)
        }
      }
      dups <- tibble::data_frame(INDIVIDUALS = dups)

      if (nrow(res$blacklist.id) > 0) {
        dups <- dplyr::filter(dups, !INDIVIDUALS %in% res$blacklist.id$INDIVIDUALS)
      }

      if (nrow(dups) > 0) {
        if (diff.pop.remove) {
          blacklist.diff.pop <- dup.filtered %>%
            dplyr::filter(ID1 %in% dups$INDIVIDUALS | ID2 %in% dups$INDIVIDUALS) %>%
            dplyr::distinct(POP_COMP) %>%
            dplyr::filter(POP_COMP == "different pop")

          if (nrow(blacklist.diff.pop) > 0) {
            res$blacklist.id <- dplyr::bind_rows(res$blacklist.id, dups)
          }
          blacklist.diff.pop <- NULL
        } else {
          whitelist.id <- dups %>%
            dplyr::left_join(geno.stats, by = "INDIVIDUALS") %>%
            dplyr::filter(GENOTYPED_PROP == max(GENOTYPED_PROP)) %>%
            dplyr::sample_n(tbl = ., size = 1) %>% # make sure only 1 is selected
            dplyr::select(INDIVIDUALS)

          if (nrow(whitelist.id) > 0) res$whitelist.id <- dplyr::bind_rows(res$whitelist.id, whitelist.id)

          blacklist.id <- dplyr::filter(dups, !INDIVIDUALS %in% whitelist.id$INDIVIDUALS) %>%
            dplyr::select(INDIVIDUALS)

          if (nrow(blacklist.id) > 0) res$blacklist.id <- dplyr::bind_rows(res$blacklist.id, blacklist.id)
        }
      }
    }
    dups <- blacklist.id <- whitelist.id <- i <- new.dups <- NULL

    res$blacklist.id <- dplyr::distinct(res$blacklist.id, INDIVIDUALS)
    res$whitelist.id <- dplyr::distinct(res$whitelist.id, INDIVIDUALS)
    message("With threshold selected, ", nrow(res$blacklist.id) ," individual(s) blacklisted")
    readr::write_tsv(x = res$blacklist.id, path = "blacklist.id.similar.tsv")
    message("Written in the directory: blacklist.id.similar.tsv")

  } else {
    message("With threshold selected, the blacklist of duplicate individuals is empty")
    res <- list(blacklist.id = tibble::tibble(INDIVIDUALS = character(0)),
                whitelist.id = tibble::tibble(INDIVIDUALS = character(0)))
  }
  return(res)
} # End remove_duplicates
