#' @name tidy_vcf

#' @title Tidy a vcf file (bi/multi-allelic).

#' @description Used internally in
#' \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.
#' Highly recommended to use \code{\link[radiator]{tidy_genomic_data}} or
#' \code{\link[radiator]{genomic_converter}}. Those two functions allows
#' to manipulate and prune the dataset with blacklists and whitelists along
#' several other filtering options.
#'
#' @param data (character string) The VCF SNPs are biallelic or haplotypes.
#'
#'
#' To make the VCF population-ready, you have to use \code{strata} argument.
#' \strong{GATK VCF files:} Some VCF have an \code{ID} column filled with \code{.},
#' the LOCUS information is all contained along the linkage group in the
#' \code{CHROM} column. To make it work with
#' \href{https://github.com/thierrygosselin/radiator}{radiator},
#' the \code{ID} column is filled with the \code{POS} column info.
#' \strong{platypus VCF files:} Some VCF files don't have an ID filed with values,
#' here the same thing as GATK VCF files above is done.
#' \strong{stacks VCF files:} with \emph{de novo} approaches, the CHROM column is
#' filled with "1", the LOCUS column correspond to the CHROM section in stacks VCF and
#' the COL column is POS -1. With a reference genome, the ID column in stacks VCF is
#' separated into "LOCUS", "COL", "STRANDS".


#' @inheritParams tidy_genomic_data

#' @param ... (optional) To pass further argument for fine-tuning the tidying
#' (details below).

#' @export
#' @rdname tidy_vcf
# @importFrom vcfR read.vcfR extract.gt vcf_field_names
#' @importFrom rlang UQ
#' @importFrom data.table melt.data.table as.data.table

#' @return The output in your global environment is a tidy data frame.

#' @details
#' \strong{Advance mode, using \emph{dots-dots-dots ...}}
#' \enumerate{
#' \item \code{whitelist.markers: }documented in \code{\link[radiator]{tidy_genomic_data}}
#' \item \code{blacklist.id: }documented in \code{\link[radiator]{tidy_genomic_data}}
#' \item \code{pop.select: }documented in \code{\link[radiator]{tidy_genomic_data}}
#' \item \code{pop.levels: }documented in \code{\link[radiator]{tidy_genomic_data}}
#' \item \code{pop.labels: }documented in \code{\link[radiator]{tidy_genomic_data}}
#' \item \code{ref.calibration: } (optional, logical)
#' Default: \code{ref.calibration = FALSE}.
#' REF/ALT alleles are designated based on count.
#' The allele with the higher overall count is called REF.
#' When the count for REF/ALT alleles is equal, the REF allele is chosen randomly.
#' REF/ALT designation can change when individuals are blacklisted.
#' If individuals are blacklisted by
#' internal filters inside this function the argument is turned on automatically.
#' If you removed individuals in another software that doesn't check this, force
#' the calibration of REF/ALT alleles: \code{ref.calibration = TRUE}.
#' \item \code{vcf.stats: } (optional, logical) Generates individuals and markers important statistics
#' helpful for filtering. These are very fast to generate and because computational
#' cost is minimal, even for huge VCFs, the default is \code{vcf.stats = TRUE}.
#' \item \code{vcf.metadata: } (optional, logical or character string)
#' With \code{vcf.metadata = FALSE}, only the genotypes are kept (GT field).
#' With \code{vcf.metadata = TRUE},
#' all the metadata contained in the \code{FORMAT} field will be kept in
#' the tidy data file. radiator is currently keeping/cleaning only these metadata:
#' \code{"DP", "AD", "GL", "PL", "GQ", "HQ", "GOF", "NR", "NV"}.
#' e.g. you only wnat AD and PL, \code{vcf.metadata = c("AD", "PL")}.
#' If yours is not in the list, submit a request.
#' Default: \code{vcf.metadata = FALSE}.
#' \item \code{common.markers}: (logical) Default: code{common.markers = TRUE}.
#' Detailed in \code{\link[radiator]{tidy_genomic_data}}.
#' \item \code{filter.mac}: (integer) To blacklist markers below a specific MAC
#' (calculated overall).
#' \item \code{filter.coverage.outliers}: (logical) To blacklist markers with low and
#' high coverage based on outlier statistics.
#' \item \code{filter.markers.missing}: (integer) To blacklist markers with too
#' many missing data. e.g. \code{filter.markers.missing = 10}, will only keep
#' markers with missing rate <= to 10%.
#' \item \code{filter.snp.read.position}: 3 options available, \code{"outliers", "q75", "iqr"}.
#' This argument will blacklist markers based on it's position on the read.
#' \code{filter.snp.read.position = "outliers"}, will remove markers based
#' on outlier statistics of the position on the reads. e.g. if majority of SNPs
#' are found between 10 and 90 pb, and very few above and below, those markers are
#' discarded. Use this function argument with dataset with problematic assembly,
#' or \emph{de novo} assembly with undocumented or poorly selected
#' mismatch threshold.
#' \item \code{filter.short.ld}: this is the \code{snp.ld} argument
#' in \code{\link[radiator]{snp_ld}}
#' \item \code{filter.long.ld}: this is the \code{ld.threshold} argument
#' in \code{\link[radiator]{snp_ld}}
#' \item \code{long.ld.missing}: this is the \code{long.ld.missing} argument
#' in \code{\link[radiator]{snp_ld}}
#' \item \code{filter.individuals.missing}: (double) Use this argument to
#' blacklist individuals with too many missing data.
#' e.g. \code{filter.individuals.missing = 0.7}, will remove individuals with >
#' 0.7 or 70% missing genotypes. This can help discover more polymorphic markers
#' with some dataset.
#' \item \code{markers.info}: (character) With default: \code{markers.info = NULL},
#' all the variable in the vcf INFO field are imported.
#' To import only DP (the SNP total read depth) and AF (the SNP allele frequency),
#' use \code{markers.info = c("DP", "AF")}.
#' Using, \code{markers.info = character(0)} will not import INFO variables.
#' \item \code{path.folder}: to write ouput in a specific path
#' (used internally in radiator). Default: \code{path.folder = getwd()}.
#' If the supplied directory doesn't exist, it's created.
#' }
#' @examples
#' \dontrun{
#' # very basic with built-in defaults (not recommended):
#'  prep.data <- radiator::tidy_vcf(data = "populations.snps.vcf")
#'
#' # Using more arguments and filters (recommended):
#' prep.data <- radiator::tidy_vcf(
#' data = "populations.snps.vcf",
#' strata = "strata_salamander.tsv",
#' vcf.stats = TRUE,
#' filter.individuals.missing = "outlier",
#' common.markers = TRUE,
#' keep.both.strands = FALSE,
#' filter.mac = 4,
#' filter.markers.missing = 50,
#' filter.snp.read.position = "outliers",
#' filter.short.ld = "maf",
#' filter.long.ld = NULL,
#' vcf.metadata = TRUE,
#' path.folder = "salamander/prep_data",
#' verbose = TRUE)
#'
#' # To view the objects generated by the function:
#' names(prep.data)
#'
#' # if you want to check the blacklisted markers:
#' bl <- prep.data$blacklist.markers
#' # or get the numbers associated with the filters applied:
#'  bl.numbers <- dplyr::count(prep.data$blacklist.markers, FILTER)
#' }

#' @seealso \code{\link[radiator]{write_seqarray}}


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

tidy_vcf <- function(
  data,
  strata = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = FALSE,
  ...) {

  # # test
  # data = "populations.snps.vcf"
  # strata = "spis-popmap-448samples.tsv"
  # vcf.stats = TRUE
  # vcf.metadata = TRUE
  # parallel.core = 8
  # verbose = TRUE
  # blacklist.id = NULL
  # whitelist.markers = NULL
  # pop.select = NULL
  # pop.levels = NULL
  # pop.labels = NULL
  # filename = NULL

  # filter.individuals.missing = "outlier"
  # common.markers = TRUE
  # keep.both.strands = FALSE
  # filter.mac = 4
  # filter.coverage.outliers = TRUE
  # filter.markers.missing = 10
  # filter.snp.read.position = "outliers"
  # filter.short.ld = "maf"
  # filter.long.ld = 0.8
  # gt.vcf.nuc = TRUE
  # gt.vcf = TRUE
  # gt = TRUE
  # gt.bin = TRUE
  # wide = FALSE
  # ref.calibration = FALSE
  # keep.gds = TRUE

  opt.change <- getOption("width")
  options(width = 70)
  timing <- proc.time()

  # Note to myself: have to integrate this...
  wide <- FALSE


  # required packages ----------------------------------------------------------
  # Check that SeqArray is installed
  if (!"SeqArray" %in% utils::installed.packages()[,"Package"]) {
    stop('Please install SeqArray for this option:\n
         devtools::install_github("zhengxwen/SeqArray")
         or the bioconductor version:
         source("https://bioconductor.org/biocLite.R")
         biocLite("SeqArray")')
  }

  if (!"SeqVarTools" %in% utils::installed.packages()[,"Package"]) {
    stop('Please install SeqVarTools for this option:\n
         source("https://bioconductor.org/biocLite.R")
         biocLite("SeqVarTools")')
  }

  if (!"gdsfmt" %in% utils::installed.packages()[,"Package"]) {
    stop('Please install gdsfmt for this option:\n
         source("https://bioconductor.org/biocLite.R")
         biocLite("gdsfmt")')
  }

  # dotslist -------------------------------------------------------------------
  dotslist <- list(...)
  want <- c("whitelist.markers", "blacklist.id", "pop.select",
            "pop.levels", "pop.labels",
            "filter.snp.read.position", "filter.mac",
            "filter.coverage.outliers", "filter.markers.missing", "filter.short.ld",
            "filter.long.ld", "filter.individuals.missing", "common.markers",
            "keep.both.strands", "path.folder",
            "ref.calibration", "gt.vcf.nuc", "gt.vcf", "gt", "gt.bin", "vcf.stats",
            "filename", "keep.gds", "vcf.metadata")
  unknowned_param <- setdiff(names(dotslist), want)

  if (length(unknowned_param) > 0) {
    stop("Unknowned \"...\" parameters ",
         stringi::stri_join(unknowned_param, collapse = " "))
  }

  radiator.dots <- dotslist[names(dotslist) %in% want]
  whitelist.markers <- radiator.dots[["whitelist.markers"]]
  blacklist.id <- radiator.dots[["blacklist.id"]]
  pop.select <- radiator.dots[["pop.select"]]
  pop.levels <- radiator.dots[["pop.levels"]]
  pop.labels <- radiator.dots[["pop.labels"]]

  filter.snp.read.position <- radiator.dots[["filter.snp.read.position"]]
  filter.mac <- radiator.dots[["filter.mac"]]
  filter.coverage.outliers <- radiator.dots[["filter.coverage.outliers"]]
  filter.markers.missing <- radiator.dots[["filter.markers.missing"]]
  filter.short.ld <- radiator.dots[["filter.short.ld"]]
  filter.long.ld <- radiator.dots[["filter.long.ld"]]
  # markers.info <- radiator.dots[["markers.info"]]
  filter.individuals.missing <- radiator.dots[["filter.individuals.missing"]]
  common.markers <- radiator.dots[["common.markers"]]
  keep.both.strands <- radiator.dots[["keep.both.strands"]]
  path.folder <- radiator.dots[["path.folder"]]

  ref.calibration <- radiator.dots[["ref.calibration"]]
  gt.vcf.nuc <- radiator.dots[["gt.vcf.nuc"]]
  gt.vcf <- radiator.dots[["gt.vcf"]]
  gt <- radiator.dots[["gt"]]
  gt.bin <- radiator.dots[["gt.bin"]]
  vcf.stats <- radiator.dots[["vcf.stats"]]
  vcf.metadata <- radiator.dots[["vcf.metadata"]]
  filename <- radiator.dots[["filename"]]
  keep.gds <- radiator.dots[["keep.gds"]]

  if (is.null(keep.gds)) keep.gds <- TRUE
  if (is.null(vcf.stats)) vcf.stats <- TRUE
  if (is.null(ref.calibration)) ref.calibration <- FALSE
  if (is.null(gt.vcf.nuc)) gt.vcf.nuc <- TRUE
  if (is.null(gt.vcf)) gt.vcf <- TRUE
  if (is.null(gt)) gt <- TRUE
  if (is.null(gt.bin)) gt.bin <- TRUE
  if (is.null(path.folder)) path.folder <- getwd()
  if (is.null(keep.both.strands)) keep.both.strands <- FALSE
  if (is.null(filter.coverage.outliers)) filter.coverage.outliers <- FALSE
  if (is.null(common.markers)) common.markers <- TRUE

  if (!gt.vcf.nuc && !gt) {
    stop("At least one of gt.vcf.nuc or gt must be TRUE")
  }

  if (!is.null(filter.snp.read.position)) {
    filter.snp.read.position <- match.arg(
      arg = filter.snp.read.position,
      choices = c("outliers", "iqr", "q75"),
      several.ok = TRUE)
  }

  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (is.null(filename)) {
    filename.rad <- stringi::stri_join("radiator_", file.date, ".rad")
  } else {
    filename.problem <- file.exists(filename)
    if (filename.problem) {
      filename.rad <- stringi::stri_join(filename, "_", file.date, ".rad")
    } else {
      filename.rad <- stringi::stri_join(filename, ".rad")
    }
  }

  # import VCF -----------------------------------------------------------------
  data <- radiator::write_seqarray(
    data = data,
    strata = strata,
    filename = filename,
    vcf.stats = vcf.stats,
    blacklist.id = blacklist.id,
    whitelist.markers = whitelist.markers,
    verbose = TRUE,
    parallel.core = parallel.core,
    keep.gds = keep.gds,
    filter.individuals.missing = filter.individuals.missing,
    common.markers = common.markers,
    keep.both.strands = keep.both.strands,
    filter.mac = filter.mac,
    filter.coverage.outliers = filter.coverage.outliers,
    filter.markers.missing = filter.markers.missing,
    filter.snp.read.position = filter.snp.read.position,
    filter.short.ld = filter.short.ld,
    filter.long.ld = filter.long.ld,
    markers.info = character(0),
    vcf.metadata = vcf.metadata,
    path.folder = path.folder
  )

  # Number of markers and tidy approach ----------------------------------------
  if (vcf.stats) {
    id.string <- data$individuals$INDIVIDUALS[data$individuals$FILTER_INDIVIDUALS_MISSING]
  } else {
    id.string <- data$individuals$INDIVIDUALS
  }

  n.markers <- length(SeqArray::seqGetData(data$vcf.connection, "variant.id"))

  if (n.markers > 100000) {
    cat("\n\n############################# IMPORTANT ###############################\n")
    message("Tidying vcf with ", n.markers, " SNPs is not optimal")
    message("and requires a computer with > 16GB RAM")
    message("Apply basic filters detailed in the function doc")
    message("Reduce to\n< 100 000 SNPs or ideally ~ 10 000 unlinked SNPs\n\n")
    tidy.vcf <- interactive_question(
      x = "\nContinue tidying the VCF (Yes/No) ?",
      answer.opt = c("Y", "N", "Yes", "No", "YES", "NO", "yes", "no"))
    if (any(c("Y", "Yes", "YES", "yes") %in% tidy.vcf)) {
      tidy.vcf <- TRUE
    } else {
      message("\nKeeping the SeqArray GDS object and file")
      tidy.vcf <- FALSE
    }
  } else {
    tidy.vcf <- TRUE
  }

  # Tidying TRUE  --------------------------------------------------------------
  if (tidy.vcf) {
    # re-calibration of ref/alt alleles ------------------------------------------
    if (ref.calibration) {
      gt.vcf.nuc.bk <- gt.vcf.nuc
      gt.bk <- gt
      gt.vcf.bk <- gt.vcf
      gt.bin.bk <- gt.bin
    }

    if (!is.null(blacklist.id)) {
      ref.calibration <- TRUE
      if (verbose) message("\nRe-calibration of REF/ALT alleles is required...")
      #overide...
      gt.vcf.nuc.bk <- gt.vcf.nuc
      gt.bk <- gt
      gt.vcf.bk <- gt.vcf
      gt.bin.bk <- gt.bin

      gt.vcf.nuc <- TRUE
      gt <- FALSE
      gt.vcf <- FALSE
      gt.bin <- FALSE
    }

    # bi- or multi-alllelic VCF --------------------------------------------------
    biallelic <- data$biallelic

    # import genotypes -----------------------------------------------------------
    # nucleotides info required to generate genepop format 001, 002, 003, 004
    data$tidy.data$GT_VCF_NUC <- SeqVarTools::getGenotypeAlleles(
      gdsobj = data$vcf.connection, use.names = TRUE)
    if (gt) {
      if (!wide) {
        if (!gt.vcf.nuc) {
          data$tidy.data$GT <- as.vector(data$tidy.data$GT_VCF_NUC)
          data$tidy.data$GT_VCF_NUC <- NULL
        } else {
          data$tidy.data$GT <- data$tidy.data$GT_VCF_NUC <- as.vector(data$tidy.data$GT_VCF_NUC)
        }
      } else {
        data$tidy.data$GT <- as.vector(data$tidy.data$GT_VCF_NUC)
        if (!gt.vcf.nuc) data$tidy.data$GT_VCF_NUC <- NULL
      }
      data$tidy.data$GT <- stringi::stri_replace_all_fixed(
        str = data$tidy.data$GT,
        pattern = c("A", "C", "G", "T", "/"),
        replacement = c("001", "002", "003", "004", ""),
        vectorize_all = FALSE) %>%
        stringi::stri_replace_na(str = ., replacement = "000000")

      if (wide) {
        data$tidy.data$GT <- tibble::as_tibble(
          matrix(data = data$tidy.data$GT,
                 nrow = nrow(id.string),
                 ncol = nrow(data$markers.meta))
        ) %>%
          magrittr::set_colnames(x = ., value = data$markers.meta$MARKERS) %>%
          tibble::add_column(
            .data = .,
            "INDIVIDUALS" = id.string,
            .before = 1
          )
      }
    }

    # more work for gt.vcf.nuc
    if (gt.vcf.nuc) {
      if (wide) {
        data$tidy.data$GT_VCF_NUC %<>%
          magrittr::set_colnames(x = ., value = data$markers.meta$MARKERS) %>%
          tibble::as_tibble(x = ., rownames = "INDIVIDUALS")
      } else {
        data$tidy.data$GT_VCF_NUC <- as.vector(data$tidy.data$GT_VCF_NUC)
      }
      data$tidy.data$GT_VCF_NUC[is.na(data$tidy.data$GT_VCF_NUC)] <- "./."
    } else {
      data$tidy.data$GT_VCF_NUC <- NULL
    }

    # gt.vcf
    if (gt.vcf) {
      data$tidy.data$GT_VCF <- SeqVarTools::getGenotype(gdsobj = data$vcf.connection, use.names = TRUE)
      if (wide) {
        data$tidy.data$GT_VCF %<>% magrittr::set_colnames(x = ., value = data$markers.meta$MARKERS) %>%
          tibble::as_tibble(x = ., rownames = "INDIVIDUALS")
      } else {
        data$tidy.data$GT_VCF <- as.vector(data$tidy.data$GT_VCF)
      }

      data$tidy.data$GT_VCF[is.na(data$tidy.data$GT_VCF)] <- "./."
      # replace(data, which(data == what), NA)
    }

    # gt.bin (plink/alt dosage format)
    if (gt.bin) {
      data$tidy.data$GT_BIN <- SeqArray::seqGetData(
        gdsfile = data$vcf.connection, var.name = "$dosage_alt")
      if (wide) {
        data$tidy.data$GT_BIN %<>% tibble::as_tibble(x = .) %>%
          magrittr::set_colnames(x = ., value = data$markers.meta$MARKERS) %>%
          tibble::add_column(.data = .,
                             "INDIVIDUALS" = id.string,
                             .before = 1)
      } else {
        data$tidy.data$GT_BIN <- as.vector(data$tidy.data$GT_BIN)
      }
    }

    # genotypes metadata ---------------------------------------------------------
    # vcf.metadata
    if (is.logical(vcf.metadata)) {
      if (vcf.metadata) {
        overwrite.metadata <- NULL
      } else {
        overwrite.metadata <- "GT"
      }
    } else {#NULL or character
      if (is.null(vcf.metadata)) {
        overwrite.metadata <- NULL
        vcf.metadata <- TRUE
      } else {
        overwrite.metadata <- vcf.metadata
        if (!"GT" %in% overwrite.metadata) {
          message("GT field always included in vcf.metadata")
          overwrite.metadata <- c("GT", overwrite.metadata)
        }
        vcf.metadata <- TRUE
      }
    }

    if (vcf.metadata) {
      if (verbose) message("Keeping vcf genotypes metadata: yes")
      # detect FORMAT fields available
      have <-  SeqArray::seqSummary(
        gdsfile = data$vcf.connection,
        varname = "annotation/format",
        check = "none", verbose = FALSE)$ID
      # current version doesn't deal well with PL with 3 fields separated with ","
      # want <- c("DP", "AD", "GL", "PL", "HQ", "GQ", "GOF", "NR", "NV", "RO", "QR", "AO", "QA")

      if (length(have) > 0) {
        want <- c("DP", "AD", "GL", "PL", "HQ", "GQ", "GOF", "NR", "NV")

        if (!is.null(overwrite.metadata)) want <- overwrite.metadata
        if (verbose) message("    genotypes metadata: ", stringi::stri_join(want, collapse = ", "))

        parse.format.list <- purrr::keep(.x = have, .p = have %in% want)
        # work on parallelization of this part
        data$tidy.data.metadata <- purrr::map(
          .x = parse.format.list, .f = parse_gds_metadata, data = data,
          verbose = verbose, parallel.core = parallel.core) %>%
          purrr::flatten(.) %>%
          purrr::flatten_df(.)
      } else {
        if (verbose) message("    genotypes metadata: none found")
        vcf.metadata <- FALSE
        data$tidy.data.metadata <- NULL
      }
    }

    # Remove or not gds connection and file --------------------------------------
    if (!keep.gds) {
      data$vcf.connection <- NULL
      if (file.exists(data$filename)) file.remove(data$filename)
    }

    # Haplotypes or biallelic VCF-------------------------------------------------
    # # recoding genotype
    # if (biallelic) {# biallelic VCF
    #   if (verbose) message("Recoding bi-allelic VCF...")
    # } else {#multi-allelic vcf
    #   if (verbose) message("Recoding VCF haplotype...")
    # }
    # input <- dplyr::rename(input, GT_VCF_NUC = GT)

    # Generating the tidy --------------------------------------------------------

    ## Note to myself: check timig with 1M SNPs to see
    ## if this is more efficient than data.table melt...

    want <- intersect(c("MARKERS", "CHROM", "LOCUS", "POS", "COL", "REF", "ALT"),
                      names(data$markers.meta))
    data$tidy.data <- suppressWarnings(
      dplyr::select(data$markers.meta, dplyr::one_of(want))) %>%
      dplyr::bind_cols(
        tibble::as_tibble(
          matrix(
            data = NA,
            nrow = data$n.markers, ncol = data$n.individuals)) %>%
          magrittr::set_colnames(x = ., value = id.string)) %>%
      data.table::as.data.table(.) %>%
      data.table::melt.data.table(
        data = .,
        id.vars = want,
        variable.name = "INDIVIDUALS",
        value.name = "GT",
        variable.factor = FALSE) %>%
      tibble::as_data_frame(.) %>%
      dplyr::select(-GT) %>%
      dplyr::mutate(
        MARKERS = factor(x = MARKERS,
                         levels = data$markers.meta$MARKERS, ordered = TRUE),
        INDIVIDUALS = factor(x = INDIVIDUALS,
                             levels = id.string,
                             ordered = TRUE)) %>%
      dplyr::arrange(MARKERS, INDIVIDUALS) %>%
      dplyr::bind_cols(data$tidy.data)

    # Check stacks AD problem --------------------------------------------------
    # Some genotypes with missing AD...
    if (vcf.metadata) {
      stacks.ad.problem <- data$tidy.data.metadata %>%
        dplyr::select(READ_DEPTH, ALLELE_REF_DEPTH,ALLELE_ALT_DEPTH) %>%
        dplyr::filter(!is.na(READ_DEPTH)) %>%
        dplyr::filter(is.na(ALLELE_REF_DEPTH) & is.na(ALLELE_ALT_DEPTH))
      stacks.ad.problem.n <- nrow(stacks.ad.problem)

      if (stacks.ad.problem.n > 0) {
        message("\n\nStacks problem detected")
        message("    missing allele depth info")
        message("    number of genotypes with problem: ", stacks.ad.problem.n)
        message("    correcting problem by adding the read depth info into AD fields...\n\n")

        stacks.ad.problem <- dplyr::select(data$tidy.data, GT_BIN) %>%
          dplyr::bind_cols(dplyr::select(data$tidy.data.metadata, READ_DEPTH, ALLELE_REF_DEPTH, ALLELE_ALT_DEPTH)) %>%
          dplyr::mutate(COL_ID = seq(1, n())) %>%
          dplyr::filter(!is.na(READ_DEPTH)) %>%
          dplyr::filter(is.na(ALLELE_REF_DEPTH) & is.na(ALLELE_ALT_DEPTH)) %>%
          dplyr::mutate(
            ALLELE_REF_DEPTH = dplyr::if_else(GT_BIN == 0, READ_DEPTH, ALLELE_REF_DEPTH),
            ALLELE_ALT_DEPTH = dplyr::if_else(GT_BIN == 2, READ_DEPTH, ALLELE_ALT_DEPTH)
          )

        data$tidy.data.metadata$ALLELE_REF_DEPTH[stacks.ad.problem$COL_ID] <- stacks.ad.problem$ALLELE_REF_DEPTH
        data$tidy.data.metadata$ALLELE_ALT_DEPTH[stacks.ad.problem$COL_ID] <- stacks.ad.problem$ALLELE_ALT_DEPTH
        # check <- data$tidy.data.metadata[stacks.ad.problem$COL_ID, ]
      }
      stacks.ad.problem.n <- stacks.ad.problem.n <- NULL
      data$tidy.data %<>% dplyr::bind_cols(data$tidy.data.metadata)
      data$tidy.data.metadata <- NULL
    }

    # re-calibration of ref/alt alleles ------------------------------------------

    if (ref.calibration) {
      if (verbose) message("\nCalculating REF/ALT alleles...")
      data$tidy.data <- radiator::change_alleles(
        data = data$tidy.data,
        biallelic = biallelic,
        parallel.core = parallel.core,
        verbose = verbose,
        gt.vcf.nuc = gt.vcf.nuc.bk,
        gt = gt.bk,
        gt.vcf = gt.vcf.bk,
        gt.bin = gt.bin.bk
      )$input
    }

    # include strata
    if (!is.null(strata)) {
      data$tidy.data <- suppressWarnings(
        dplyr::left_join(
          data$tidy.data,
          dplyr::filter(data$individuals, INDIVIDUALS %in% id.string) %>%
            dplyr::select(INDIVIDUALS, POP_ID = STRATA),
          by = "INDIVIDUALS")
      )
    }


    # Re ordering columns
    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "ID", "COL", "INDIVIDUALS",
              "STRATA", "POP_ID",
              "REF", "ALT", "GT_VCF", "GT_VCF_NUC", "GT", "GT_BIN",
              "POLYMORPHIC")

    data$tidy.data <- suppressWarnings(
      dplyr::select(data$tidy.data, dplyr::one_of(want), dplyr::everything()))

    # Sort id
    if (!is.null(strata)) {
      data$tidy.data <- dplyr::arrange(data$tidy.data, POP_ID, INDIVIDUALS)
    } else {
      data$tidy.data <- dplyr::arrange(data$tidy.data, INDIVIDUALS)
    }


    #write tidy
    radiator::write_rad(data = data$tidy.data, path = file.path(path.folder, filename.rad))
  } #tidy.vcf

  # Results --------------------------------------------------------------------
  timing <- proc.time() - timing
  if (verbose) {
    if (vcf.stats && tidy.vcf) {
      message("\nTiming for tidying and generating stats of vcf: ", round(timing[[3]]), " sec")
    } else {
      if (tidy.vcf) {
        message("\nTiming for tidying vcf: ", round(timing[[3]]), " sec")
      } else {
        message("\nTiming for generating GDS: ", round(timing[[3]]), " sec")
      }
    }
  }
  options(width = opt.change)
  return(data)
}#End tidy_vcf


# Internal nested Function -----------------------------------------------------
#' @title parse_gds_metadata
#' @description function to parse the format field and tidy the results of VCF
#' @rdname parse_gds_metadata
#' @keywords internal
#' @export
parse_gds_metadata <- function(
  x, data = NULL, verbose = TRUE, parallel.core = parallel::detectCores() - 1) {

  res <- list()
  format.name <- x
  # format.name <- x <- "DP"
  # format.name <- x <- "AD"
  # format.name <- x <- "GL"
  # format.name <- x <- "PL"
  # format.name <- x <- "HQ"
  # format.name <- x <- "GQ"
  # format.name <- x <- "GOF"
  # format.name <- x <- "NR"
  # format.name <- x <- "NV"
  # format.name <- x <- "RO" # not yet implementer
  # format.name <- x <- "QR" # not yet implementer
  # format.name <- x <- "AO" # not yet implementer
  # format.name <- x <- "QA" # not yet implementer

  if (verbose) message("\nParsing and tidying: ", format.name)

  if (format.name == "AD") {
    if (verbose) message("AD column: splitting into ALLELE_REF_DEPTH and ALLELE_ALT_DEPTH")
    res$AD <- SeqArray::seqGetData(gdsfile = data$vcf.connection,
                                   var.name = "annotation/format/AD")$data %>%
      tibble::as_tibble(.)
    column.vec <- seq_along(res$AD)
    res$AD <- tibble::tibble(ALLELE_REF_DEPTH = res$AD[, column.vec %% 2 == 1] %>%
                               as.matrix(.) %>%
                               as.vector(.),
                             ALLELE_ALT_DEPTH = res$AD[, column.vec %% 2 == 0] %>%
                               as.matrix(.) %>%
                               as.vector(.))
    split.vec <- split_vec_row(x = res$AD, 3, parallel.core = parallel.core)
    res$AD$ALLELE_REF_DEPTH <- clean_ad(x = res$AD$ALLELE_REF_DEPTH, split.vec = split.vec,
                                        parallel.core = parallel.core)

    res$AD$ALLELE_ALT_DEPTH <- clean_ad(x = res$AD$ALLELE_ALT_DEPTH, split.vec = split.vec,
                                        parallel.core = parallel.core)
    split.vec <- NULL
    # test1 <- res$AD
  } # End AD


  # Read depth
  if (format.name == "DP") {
    if (verbose) message("DP column: cleaning and renaming to READ_DEPTH")
    res$DP <- tibble::tibble(READ_DEPTH = SeqArray::seqGetData(
      gdsfile = data$vcf.connection,
      var.name = "annotation/format/DP")$data %>% as.vector(.))
    # test <- res$DP
  } # End DP

  # Cleaning HQ: Haplotype quality as phred score
  if (format.name == "HQ") {
    res$HQ <- tibble::tibble(HQ = SeqArray::seqGetData(
      gdsfile = data$vcf.connection,
      var.name = "annotation/format/HQ")$data %>% as.vector(.))
    # test <- res$HQ

    # check HQ and new stacks version with no HQ
    all.missing <- nrow(res$HQ)
    if (all.missing != 0) {
      if (verbose) message("HQ column: Haplotype Quality")
    } else {
      message("HQ values are all missing: removing column")
      res$HQ <- NULL
    }
  } # End HQ

  # Cleaning GQ: Genotype quality as phred score
  if (format.name == "GQ") {
    if (verbose) message("GQ column: Genotype Quality")
    res$GQ <- tibble::tibble(GQ = SeqArray::seqGetData(
      gdsfile = data$vcf.connection,
      var.name = "annotation/format/GQ")$data %>% as.vector(.))
    # test <- res$GQ
  } # End GQ

  # GL cleaning
  if (format.name == "GL") {
    if (verbose) message("GL column: cleaning Genotype Likelihood column")
    res$GL <- SeqArray::seqGetData(gdsfile = data$vcf.connection,
                                   var.name = "annotation/format/GL")$data %>%
      tibble::as_tibble(.)
    column.vec <- seq_along(res$GL)
    res$GL <- tibble::tibble(GL_HOM_REF = res$GL[, column.vec %% 3 == 1] %>%
                               as.matrix(.) %>%
                               as.vector(.),
                             GL_HET = res$GL[, column.vec %% 3 == 2] %>%
                               as.matrix(.) %>%
                               as.vector(.),
                             GL_HOM_ALT = res$GL[, column.vec %% 3 == 0] %>%
                               as.matrix(.) %>%
                               as.vector(.))
    res$GL[res$GL == "NaN"] <- NA
  } # End GL

  # Cleaning GOF: Goodness of fit value
  if (format.name == "GOF") {
    if (verbose) message("GOF column: Goodness of fit value")
    res$GOF <- tibble::tibble(GOF = SeqArray::seqGetData(
      gdsfile = data$vcf.connection,
      var.name = "annotation/format/GOF")$data %>% as.vector(.))
    # test <- res$GOF
  } # End GOF

  # Cleaning NR: Number of reads covering variant location in this sample
  if (format.name == "NR") {
    if (verbose) message("NR column: splitting column into the number of variant")
    res$NR <- tibble::tibble(NR = SeqArray::seqGetData(
      gdsfile = data$vcf.connection,
      var.name = "annotation/format/NR")$data %>% as.vector(.))
    # test <- res$NR
  }#End cleaning NR column

  # Cleaning NV: Number of reads containing variant in this sample
  if (format.name == "NV") {
    if (verbose) message("NV column: splitting column into the number of variant")
    res$NR <- tibble::tibble(NV = SeqArray::seqGetData(
      gdsfile = data$vcf.connection,
      var.name = "annotation/format/NV")$data %>% as.vector(.))
    # test <- res$NV
  }#End cleaning NV column

  # RO
  # Cleaning RO: Reference allele observation count
  # if (format.name == "RO") {
  #   if (verbose) message("RO column: splitting column into the number of variant")
  #   res$RO <- tibble::tibble(RO = SeqArray::seqGetData(
  #     gdsfile = data$vcf.connection,
  #     var.name = "annotation/format/RO")$data %>% as.vector(.))
  #   # test <- res$RO
  # }#End cleaning RO column

  # # Cleaning QR: Sum of quality of the reference observations
  # if (format.name == "QR") {
  #   if (verbose) message("QR column: splitting column into the number of variant")
  #   res$QR <- tibble::tibble(QR = SeqArray::seqGetData(
  #     gdsfile = data$vcf.connection,
  #     var.name = "annotation/format/QR")$data %>% as.vector(.))
  #   # test <- res$QR
  # }#End cleaning QR column


  # # Cleaning AO: Alternate allele observation count
  # if (format.name == "AO") {
  #   if (verbose) message("AO column: splitting column into the number of variant")
  #   res$AO <- tibble::tibble(AO = SeqArray::seqGetData(
  #     gdsfile = data$vcf.connection,
  #     var.name = "annotation/format/AO")$data %>% as.vector(.))
  #   # test <- res$AO
  # }#End cleaning AO column

  # # Cleaning QA: Sum of quality of the alternate observations
  # if (format.name == "QA") {
  #   if (verbose) message("QA column: splitting column into the number of variant")
  #   res$QA <- tibble::tibble(QA = SeqArray::seqGetData(
  #     gdsfile = data$vcf.connection,
  #     var.name = "annotation/format/QA")$data %>% as.vector(.))
  #   # test <- res$QA
  # }#End cleaning QA column



  # test <- res$AD %>%
  #   dplyr::bind_cols(READ_DEPTH = SeqArray::seqGetData(
  #   gdsfile = data$vcf.connection,
  #   var.name = "annotation/format/DP")$data %>% as.vector(.))

  # if (gather.data) {
  #   x <- dplyr::mutate(x, ID = seq(1, n()))
  #   x <- data.table::melt.data.table(
  #     data = data.table::as.data.table(x),
  #     id.vars = "ID",
  #     variable.name = "INDIVIDUALS",
  #     value.name = format.name,
  #     variable.factor = FALSE) %>%
  #     tibble::as_data_frame(.) %>%
  #     dplyr::select(-ID, -INDIVIDUALS)
  # }
  return(res)
}#End parse_gds_metadata

#' @title split_vcf_id
#' @description split VCF ID in parallel
#' @rdname split_vcf_id
#' @keywords internal
#' @export

split_vcf_id <- function(x) {
  res <- dplyr::rename(x, ID = LOCUS) %>%
    tidyr::separate(data = ., col = ID, into = c("LOCUS", "COL"),
                    sep = "_", extra = "drop", remove = FALSE) %>%
    dplyr::mutate_at(.tbl = ., .vars = c("CHROM", "POS", "LOCUS"), .funs = as.character) %>%
    dplyr::mutate_at(.tbl = ., .vars = c("CHROM", "POS", "LOCUS"), .funs = clean_markers_names) %>%
    tidyr::unite(MARKERS, c(CHROM, LOCUS, POS), sep = "__", remove = FALSE)
  return(res)
}#End split_vcf_id



#' @title clean_ad
#' @description Clean allele depth field in VCF.
#' AD column: splitting coverage info into ALLELE_REF_DEPTH and ALLELE_ALT_DEPTH
#' @rdname clean_ad
#' @keywords internal
#' @export
clean_ad <- function(x, split.vec, parallel.core = parallel::detectCores() - 1) {
  clean <- function(x) {
    x <- as.integer(replace_by_na(data = x, what = 0))
  }

  x <- split(x = x, f = split.vec) %>%
    .radiator_parallel(
      X = ., FUN = clean, mc.cores = parallel.core) %>%
    purrr::flatten_int(.)
  return(x)
}#End clean_ad


#' @title clean_pl
#' @description Clean PL column.
#' PL column (normalized, phred-scaled likelihoods for genotypes):
#' separating into PROB_HOM_REF, PROB_HET and PROB_HOM_ALT")
#' Value 1: probability that the site is homozgyous REF
#' Value 2: probability that the sample is heterzygous
#' Value 2: probability that it is homozygous ALT
#' @rdname clean_pl
#' @keywords internal
#' @export
clean_pl <- function(x, split.vec, parallel.core = parallel::detectCores() - 1) {
  clean <- function(x) {
    res <- x %>%
      dplyr::mutate(
        PL = dplyr::if_else(GT_VCF == "./.", NA_character_, PL)) %>%
      tidyr::separate(
        data = ., PL, c("PROB_HOM_REF", "PROB_HET", "PROB_HOM_ALT"),
        sep = ",", extra = "drop", remove = FALSE) %>%
      dplyr::mutate_at(
        .tbl = ., .vars = c("PROB_HOM_REF", "PROB_HET", "PROB_HOM_ALT"),
        .funs = as.numeric) %>%
      dplyr::select(-GT_VCF)
    return(res)
  }#End clean
  x <- dplyr::bind_cols(
    dplyr::select(x, -PL),
    dplyr::ungroup(x) %>%
      dplyr::select(GT_VCF, PL) %>%
      split(x = ., f = split.vec) %>%
      .radiator_parallel(
        X = ., FUN = clean, mc.cores = parallel.core) %>%
      dplyr::bind_rows(.))
  return(x)
}#End clean_pl

#' @title clean_gl
#' @description Clean GL column.
#' GL column: cleaning Genotype Likelihood column
#' @rdname clean_gl
#' @keywords internal
#' @export

clean_gl <- function(x, split.vec, parallel.core = parallel::detectCores() - 1) {
  x <- x %>%
    dplyr::mutate(
      GL = dplyr::if_else(GT_VCF == "./.", NA_character_, GL),
      GL = suppressWarnings(
        stringi::stri_replace_all_fixed(
          GL, c(".,.,.", ".,", ",."), c("NA", "", ""), vectorize_all = FALSE))
    )

  # check GL and new stacks version with no GL
  all.missing <- all(is.na(x$GL))

  if (!all.missing) {
    gl.clean <- max(
      unique(stringi::stri_count_fixed(
        str = unique(sample(x = x$GL, size = 100, replace = FALSE)),
        pattern = ",")
      ), na.rm = TRUE
    )

    if (gl.clean == 2) {
      message("GL column: separating into PROB_HOM_REF, PROB_HET and PROB_HOM_ALT")
      # Value 1: probability that the site is homozgyous REF
      # Value 2: probability that the sample is heterzygous
      # Value 2: probability that it is homozygous ALT
      # system.time(input2 <- input %>%
      #   tidyr::separate(data = ., GL, c("PROB_HOM_REF", "PROB_HET", "PROB_HOM_ALT"), sep = ",", extra = "drop", remove = FALSE) %>%
      #   dplyr::mutate_at(.tbl = ., .vars = c("PROB_HOM_REF", "PROB_HET", "PROB_HOM_ALT"), .funs = as.numeric)
      # )
      clean <- function(x) {
        res <- x %>%
          tidyr::separate(
            data = ., GL, c("PROB_HOM_REF", "PROB_HET", "PROB_HOM_ALT"),
            sep = ",", extra = "drop", remove = FALSE) %>%
          dplyr::mutate_at(
            .tbl = ., .vars = c("PROB_HOM_REF", "PROB_HET", "PROB_HOM_ALT"),
            .funs = as.numeric)
        return(res)
      }

      x <- dplyr::bind_cols(
        dplyr::select(x, -GL),
        dplyr::ungroup(x) %>%
          dplyr::select(GL) %>%
          split(x = ., f = split.vec) %>%
          .radiator_parallel(
            X = ., FUN = clean, mc.cores = parallel.core) %>%
          dplyr::bind_rows(.))

    } else {
      x$GL <- suppressWarnings(as.numeric(x$GL))
    }
  } else {
    message("    GL values are all missing: removing column")
    x <- dplyr::select(x, -GL)
  }
  return(x)
}#End clean_gl

#' @title clean_nr
#' @description Cleaning NR: Number of reads covering variant location in
#' this sample
#' @rdname clean_nr
#' @keywords internal
#' @export
clean_nr <- function(x, split.vec, parallel.core = parallel::detectCores() - 1){
  nr.col <- max(unique(stringi::stri_count_fixed(str = unique(x$NR), pattern = ","))) + 1
  nr.col.names <- stringi::stri_join(rep("NR_", nr.col), seq(1:nr.col))
  nr.col <- NULL

  clean <- function(x, nr.col.names = NULL) {
    res <- tidyr::separate(data = x, col = NR, into = nr.col.names,
                           sep = ",", extra = "drop", remove = FALSE)
    return(res)
  }#End clean

  x <- dplyr::bind_cols(
    dplyr::select(x, -NR),
    dplyr::ungroup(x) %>%
      dplyr::mutate(NR = dplyr::if_else(GT_VCF == "./.", NA_character_, NR)) %>%
      dplyr::select(NR) %>%
      split(x = ., f = split.vec) %>%
      .radiator_parallel_mc(
        X = ., FUN = clean, mc.cores = parallel.core,
        nr.col.names = nr.col.names) %>%
      dplyr::bind_rows(.))
  return(x)
}#End clean_nr

#' @title clean_nv
#' @description Cleaning NV: Number of reads containing variant in this sample
#' @rdname clean_nv
#' @keywords internal
#' @export
clean_nv <- function(x, split.vec, parallel.core = parallel::detectCores() - 1) {

  nv.col <- max(unique(stringi::stri_count_fixed(str = unique(x$NV), pattern = ","))) + 1
  nv.col.names <- stringi::stri_join(rep("NV_", nv.col), seq(1:nv.col))
  nv.col <- NULL

  clean <- function(x, nv.col.names = NULL) {
    res <- tidyr::separate(
      data = x, col = NV, into = nv.col.names,
      sep = ",", extra = "drop", remove = FALSE)
    return(res)
  }

  x <- dplyr::bind_cols(
    dplyr::select(x, -NV),
    dplyr::ungroup(x) %>%
      dplyr::mutate(NV = dplyr::if_else(GT_VCF == "./.", NA_character_, NV)) %>%
      dplyr::select(NV) %>%
      split(x = ., f = split.vec) %>%
      .radiator_parallel_mc(
        X = ., FUN = clean, mc.cores = parallel.core,
        nv.col.names = nv.col.names) %>%
      dplyr::bind_rows(.))
  return(x)
}#End clean_nv

#' @title strata_vcf
#' @description Manage strata
#' @rdname strata_vcf
#' @keywords internal
#' @export
strata_vcf <- function(strata, input, blacklist.id) {

  if (is.null(strata)) {
    message("No strata file provided")
    message("    generating a strata with 1 grouping")
    strata.df <- dplyr::distinct(input, INDIVIDUALS) %>%
      dplyr::mutate(STRATA = rep("pop1", n()))
  } else {
    if (is.vector(strata)) {
      strata.df <- suppressMessages(readr::read_tsv(
        file = strata, col_names = TRUE,
        col_types = readr::cols(.default = readr::col_character())
      ))
    } else {
      strata.df <- strata
      strata.df <- dplyr::mutate_all(.tbl = strata.df, .funs = as.character)
    }
  }

  colnames(strata.df) <- stringi::stri_replace_all_fixed(
    str = colnames(strata.df),
    pattern = "STRATA",
    replacement = "POP_ID",
    vectorize_all = FALSE
  )

  # Remove potential whitespace in pop_id
  strata.df$POP_ID <- clean_pop_names(strata.df$POP_ID)
  colnames.strata <- colnames(strata.df)

  # clean ids
  strata.df$INDIVIDUALS <- clean_ind_names(strata.df$INDIVIDUALS)

  strata.df <- dplyr::distinct(strata.df, POP_ID, INDIVIDUALS, .keep_all = TRUE)

  if (!is.null(strata)) {
    id.vcf <- dplyr::distinct(input, INDIVIDUALS) %>%
      dplyr::mutate(INDIVIDUALS = clean_ind_names(INDIVIDUALS)) %>%
      purrr::flatten_chr(.)

    strata.df <- dplyr::filter(strata.df, INDIVIDUALS %in% id.vcf)
  }


  # filtering the strata if blacklist id available
  if (!is.null(blacklist.id)) {
    strata.df <- dplyr::anti_join(x = strata.df, y = blacklist.id, by = "INDIVIDUALS")
  }

  return(strata.df)
}#End strata_vcf
