# tidy_plink
#' @name tidy_plink
#' @title Tidy plink \code{.tped} and {.tfam} files in a tidy data frame
#' @description Transform genomic data set produced by massive parallel
#' sequencing pipeline (e.g.GBS/RADseq,
#' SNP chip, DArT, etc) into a tidy format.
#' The use of blacklist and whitelist along
#' several filtering options are available to prune the dataset.
#' Several arguments are available to make your data population-wise and easily
#' rename the pop id.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @param data The plink file ending with \code{.tped}, the corresponding
#' \code{.tfam} must be in the same directory.

#' @inheritParams read_strata
#' @inheritParams radiator_common_arguments

# @param max.marker (integer, optional) For large PLINK files,
# useful to subsample marker number. e.g. if the data set
# contains 200 000 markers and \code{max.marker = 10000}, 10000 markers are
# subsampled randomly from the 200000 markers. If you need specific markers,
# use \code{whitelist.markers} argument.
# Default: \code{max.marker = NULL}.


#' @references
#' PLINK: a tool set for whole-genome association and population-based linkage
#' analyses.
#' American Journal of Human Genetics. 2007: 81: 559–575. doi:10.1086/519795

#' @export
#' @rdname tidy_plink
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

tidy_plink <- function(data, strata = NULL, verbose = FALSE, ...) {

  # NOTE TO MYSELF: OLD CODE THAT COULD BE COMPLETELY RE_WRITTEN IF NECESSARY...
  # BUT WHO USE PLINK file as input file these days...?



  # dotslist -------------------------------------------------------------------
  dotslist <- rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE)
  want <- c(
    "whitelist.markers", "blacklist.id",
    "blacklist.genotype",
    "pop.select", "pop.levels", "pop.labels",
    "filter.snp.read.position",
    "filter.mac", "maf.thresholds",
    "filter.coverage.outliers",
    "filter.markers.missing",
    "filter.short.ld", "filter.long.ld", "long.ld.missing", "ld.method", "snp.ld",
    "filter.individuals.missing",
    "filter_common_markers", "common.markers",
    "filter_monomorphic",
    "filter.strands",
    "path.folder",
    "ref.calibration", "gt.vcf.nuc", "gt.vcf", "gt", "gt.bin",
    "filename", "keep.gds",
    "vcf.metadata", "vcf.stats"
  )
  unknowned_param <- setdiff(names(dotslist), want)

  if (length(unknowned_param) > 0) {
    rlang::abort("Unknowned \"...\" parameters ",
         stringi::stri_join(unknowned_param, collapse = " "))
  }

  radiator.dots <- dotslist[names(dotslist) %in% want]


  # argument <- radiator.dots[["argument"]]

  whitelist.markers <- radiator.dots[["whitelist.markers"]]
  blacklist.id <- radiator.dots[["blacklist.id"]]
  blacklist.genotype <- radiator.dots[["blacklist.genotype"]]

  pop.levels <- radiator.dots[["pop.levels"]]
  pop.labels <- radiator.dots[["pop.labels"]]
  pop.select <- radiator.dots[["pop.select"]]

  filter.snp.read.position <- radiator.dots[["filter.snp.read.position"]]

  filter.mac <- radiator.dots[["filter.mac"]]
  maf.thresholds <- radiator.dots[["maf.thresholds"]]
  if (!is.null(maf.thresholds)) message("Deprecated argument, see details")

  filter.coverage.outliers <- radiator.dots[["filter.coverage.outliers"]]
  filter.markers.missing <- radiator.dots[["filter.markers.missing"]]

  filter.short.ld <- radiator.dots[["filter.short.ld"]]
  filter.long.ld <- radiator.dots[["filter.long.ld"]]
  long.ld.missing <- radiator.dots[["long.ld.missing"]]
  if (is.null(long.ld.missing)) long.ld.missing <- FALSE
  ld.method <- radiator.dots[["ld.method"]]
  if (is.null(ld.method)) ld.method <- "r2"
  snp.ld <- radiator.dots[["filter.short.ld"]]
  if (!is.null(snp.ld)) message("Deprecated argument, see details")




  filter.individuals.missing <- radiator.dots[["ld.method"]]

  filter_common_markers <- radiator.dots[["filter_common_markers"]]
  if (is.null(filter_common_markers)) filter_common_markers <- TRUE
  common.markers <- radiator.dots[["filter_common_markers"]]
  if (!is.null(common.markers)) message("Deprecated argument, see details")

  filter_monomorphic <- radiator.dots[["filter_monomorphic"]]
  filter.strands <- radiator.dots[["filter.strands"]]

  path.folder <- radiator.dots[["path.folder"]]
  if (is.null(path.folder)) {
    path.folder <- getwd()
  } else {
    if (!dir.exists(path.folder)) dir.create(path.folder)
  }

  ref.calibration <- radiator.dots[["ref.calibration"]]
  if (is.null(ref.calibration)) ref.calibration <- FALSE

  gt.vcf.nuc <- radiator.dots[["gt.vcf.nuc"]]
  if (is.null(gt.vcf.nuc)) gt.vcf.nuc <- TRUE

  gt.vcf <- radiator.dots[["gt.vcf"]]
  if (is.null(gt.vcf)) gt.vcf <- TRUE

  gt <- radiator.dots[["gt"]]
  if (is.null(gt)) gt <- TRUE

  gt.bin <- radiator.dots[["gt.bin"]]
  if (is.null(gt.bin)) gt.bin <- TRUE

  filename <- radiator.dots[["filename"]]

  keep.gds <- radiator.dots[["keep.gds"]]
  if (is.null(keep.gds)) keep.gds <- TRUE

  vcf.stats <- radiator.dots[["vcf.stats"]]
  if (is.null(vcf.stats)) vcf.stats <- TRUE

  vcf.metadata <- radiator.dots[["vcf.metadata"]]
  if (is.null(vcf.metadata)) vcf.metadata <- TRUE


  # Import whitelist of markers-------------------------------------------------

  # Strata ---------------------------------------------------------------------
  strata.df <- read_strata(
    strata = strata,
    pop.levels = pop.levels,
    pop.labels = pop.labels,
    pop.select = pop.select,
    blacklist.id = blacklist.id,
    verbose = verbose)

  # Importing file -------------------------------------------------------------
  if (verbose) message("Importing the PLINK files...")
  tfam <- data.table::fread(
    input = stringi::stri_replace_all_fixed(
      str = data,
      pattern = ".tped", replacement = ".tfam", vectorize_all = FALSE
    ),
    sep = " ",
    header = FALSE,
    stringsAsFactors = FALSE,
    verbose = FALSE,
    select = c(1,2),
    colClasses = list(character = c(1,2)),
    col.names = c("POP_ID", "INDIVIDUALS"),
    showProgress = TRUE,
    data.table = FALSE) %>%
    tibble::as_data_frame(.) %>%
    dplyr::mutate_at(.tbl = ., .vars = "INDIVIDUALS",
                     .funs = clean_ind_names) %>%
    dplyr::mutate_at(.tbl = ., .vars = "POP_ID",
                     .funs = clean_pop_names)

  # if no strata tfam = strata.df
  if (is.null(strata)) {
    strata.df <- strata.df <- read_strata(
      strata = tfam,
      pop.levels = pop.levels,
      pop.labels = pop.labels,
      pop.select = pop.select,
      blacklist.id = blacklist.id,
      verbose = FALSE)
  }

  tped.header.prep <- tfam %>%
    dplyr::select(INDIVIDUALS) %>%
    dplyr::mutate(
      NUMBER = seq(1, n()),
      ALLELE1 = rep("A1", n()), ALLELE2 = rep("A2", n())
    ) %>%
    tidyr::gather(ALLELES_GROUP, ALLELES, -c(INDIVIDUALS, NUMBER)) %>%
    dplyr::arrange(NUMBER) %>%
    dplyr::select(-ALLELES_GROUP) %>%
    tidyr::unite(INDIVIDUALS_ALLELES, c(INDIVIDUALS, ALLELES), sep = "_", remove = FALSE) %>%
    dplyr::arrange(NUMBER) %>%
    dplyr::mutate(NUMBER = seq(from = (1 + 4), to = n() + 4)) %>%
    dplyr::select(-ALLELES)

  tped.header.names <- c("LOCUS", tped.header.prep$INDIVIDUALS_ALLELES)
  tped.header.integer <- c(2, tped.header.prep$NUMBER)

  # import PLINK
  input <- data.table::fread(
    input = data,
    sep = " ",
    header = FALSE,
    stringsAsFactors = FALSE,
    verbose = FALSE,
    select = tped.header.integer,
    col.names = tped.header.names,
    showProgress = TRUE,
    data.table = FALSE) %>%
    tibble::as_data_frame(.) %>%
    dplyr::mutate(LOCUS = as.character(LOCUS))

  # Filter with whitelist of markers
  if (!is.null(whitelist.markers)) {
    if (verbose) message("Filtering with whitelist of markers")
    input  %<>% suppressWarnings(
      dplyr::semi_join(
        whitelist.markers,
        by = intersect(colnames(whitelist.markers), colnames(input)))
    )
  }

  # To reduce the size of the dataset we subsample the markers with max.marker
  # if (!is.null(max.marker)) {
  #   if (verbose) message("Using the max.marker to reduce the size of the dataset")
  #   input <- dplyr::sample_n(tbl = input, size = max(as.numeric(max.marker)), replace = FALSE)
  #
  #   max.marker.subsample.select <- input %>%
  #     dplyr::distinct(LOCUS, .keep_all = TRUE) %>%
  #     dplyr::arrange(LOCUS)
  #
  #   readr::write_tsv(# save results
  #     x = max.marker.subsample.select,
  #     path = "max.marker.subsample.select.tsv")
  # }

  # Make tidy
  if (verbose) message("Tidying the PLINK file ...")
  # Filling GT and new separating INDIVIDUALS from ALLELES
  # combining alleles
  input <- tidyr::gather(data = input, key = INDIVIDUALS_ALLELES, value = GT, -LOCUS)

  # detect GT coding
  if (verbose) message("Scanning for PLINK tped genotype coding")
  detect.gt.coding <- unique(sample(x = input$GT, size = 100, replace = FALSE))
  gt.letters <- c("A", "C", "G", "T")

  if (TRUE %in% unique(gt.letters %in% detect.gt.coding)) {
    if (verbose) message("    genotypes coded with letters")
    gt.letters.df <- tibble::data_frame(GT = c("A", "C", "G", "T", "0"), NEW_GT = c("001", "002", "003", "004", "000"))
    input  %<>% dplyr::left_join(
      gt.letters.df, by = "GT") %>%
      dplyr::select(-GT) %>%
      dplyr::rename(GT = NEW_GT)
    gt.letters.df <- NULL
  } else {
    if (verbose) message("    genotypes coded with integers")
    input  %<>%
      dplyr::mutate(GT = stringi::stri_pad_left(str = GT, width = 3, pad = "0"))
  }
  detect.gt.coding <- gt.letters <- NULL


  input  %<>%
    tidyr::separate(
      data = .,
      col = INDIVIDUALS_ALLELES,
      into = c("INDIVIDUALS", "ALLELES"),
      sep = "_") %>%
    dplyr::group_by(LOCUS, INDIVIDUALS) %>%
    tidyr::spread(data = ., key = ALLELES, value = GT) %>%
    dplyr::ungroup(.) %>%
    tidyr::unite(data = ., col = GT, A1, A2, sep = "") %>%
    dplyr::select(LOCUS, INDIVIDUALS, GT)

  # population levels and strata
  if (verbose) message("Integrating the tfam/strata file...")

  input %<>% dplyr::left_join(strata.df, by = "INDIVIDUALS")

  # removing untyped markers across all-pop
  remove.missing.gt <- input %>%
    dplyr::select(LOCUS, GT) %>%
    dplyr::filter(GT != "000000")

  untyped.markers <- dplyr::n_distinct(input$LOCUS) - dplyr::n_distinct(remove.missing.gt$LOCUS)
  if (untyped.markers > 0) {
    if (verbose) message("Number of marker with 100 % missing genotypes: ", untyped.markers)
    input <- suppressWarnings(
      dplyr::semi_join(input,
                       remove.missing.gt %>%
                         dplyr::distinct(LOCUS, .keep_all = TRUE),
                       by = "LOCUS")
    )
  }

  # Unused objects
  tped.header.prep <- tped.header.integer <- tped.header.names <- remove.missing.gt <- NULL

  # detect if biallelic give vcf style genotypes
  # biallelic <- radiator::detect_biallelic_markers(input)
  input.temp <- radiator::calibrate_alleles(data = input,
                                            verbose = verbose)
  return(res = list(input = input.temp$input, biallelic = input.temp$biallelic))
} # End import PLINK
