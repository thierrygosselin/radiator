#' @name tidy_vcf

#' @title Tidy a vcf file (bi and multi-allelic).

#' @description Used internally in
#' \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.
#' Highly recommended to use \code{\link[radiator]{filter_rad}} to reduce
#' the number of markers. Advance options below are also available to
#' to manipulate and prune the dataset with blacklists and whitelists along
#' several other filtering options.
#'
#' @param data (VCF file, character string) The VCF SNPs are biallelic or haplotypes.
#' To make the VCF population-ready, the argument \code{strata} is required.
#'
#'
#' @inheritParams radiator_common_arguments
#' @inheritParams tidy_genomic_data

#' @param ... (optional) To pass further argument for fine-tuning the tidying
#' (read below).

#' @export
#' @rdname tidy_vcf

#' @return The output in your global environment is a tidy data frame, the GDS file
#' generated is in the working directory under the name given during function execution.


#' @section VCF file format:
#'
#' \strong{GATK:} Some VCF have an \code{ID} column filled with \code{.},
#' the LOCUS information is all contained along the linkage group in the
#' \code{CHROM} column. To make it work with
#' \href{https://github.com/thierrygosselin/radiator}{radiator},
#' the \code{ID} column is filled with the \code{POS} column info.
#'
#' \strong{platypus:} Some VCF files don't have an ID filed with values,
#' here the same thing as GATK VCF files above is done.
#'
#' \strong{stacks:} with \emph{de novo} approaches, the CHROM column is
#' filled with "1", the LOCUS column correspond to the CHROM section in stacks VCF and
#' the COL column is POS -1. With a reference genome, the ID column in stacks VCF is
#' separated into "LOCUS", "COL", "STRANDS".
#'
#' \emph{stacks problem: } current version as some intrinsic problem with
#' missing allele depth info, during the tidying process a message will
#' highlight the number of genotypes impacted by the problem. When possible, the
#' problem is corrected by adding the read depth info into the allele depth field.


#' @section Advance mode, using \emph{dots-dots-dots ...}:
#'
#' The arguments below are not available using code completion (e.g. with TAB),
#' consequently any misspelling will generate an error or be ignored.
#'
#' \emph{dots-dots-dots ...} arguments names and values are reported and written
#' in the working directory when \code{internal = FALSE} and \code{verbose = TRUE}.
#'
#' \strong{General arguments: }
#' \enumerate{
#' \item \code{path.folder}: to write ouput in a specific path
#' (used internally in radiator). Default: \code{path.folder = getwd()}.
#' If the supplied directory doesn't exist, it's created.
#' \item \code{internal: } (optional, character)
#' Default (\code{internal = FALSE}). A folder is generated to write the files.
#' \item \code{random.seed}: (integer, optional) For reproducibility, set an integer
#' that will be used inside codes that uses randomness. With default,
#' a random number is generated, printed and written in the directory.
#' Default: \code{random.seed = NULL}.
#' \item \code{parameters} It's a parameter file where radiator output results of
#' filtering. Used internally.
#' Default: \code{parameters = NULL}.
#' }

#' \strong{tidying arguments/behavior:}
#' \enumerate{
#' \item \code{tidy.vcf:} (optional, logical)
#' Default: \code{tidy.vcf = TRUE}. But you can always stop the process after
#' the creation of the GDS file (equivalent of running \code{\link{write_seqarray}}).
#' \item \code{tidy.check:} (optional, logical)
#' Default: \code{tidy.check = TRUE}. By default, the number of markers just before
#' tidying is checked. Tidying a VCF file with more than 20000 markers is
#' sub-optimal:
#' \itemize{
#' \item a computer with lots of RAM is required
#' \item it's very slow to generate
#' \item it's very slow to run codes after
#' \item for most non model species this number of markers is not realistic...
#' }
#' Consequently, the function execution is suspended and user are asked if they
#' still want to continue with the tidying or stop and keep the GDS file/object.
#'
#' This behavior can be annoying, \emph{if the user knows what he's doing}, to turn off
#' use: \code{tidy.check = FALSE}.
#' \item \code{ref.calibration: } (optional, logical)
#' Default: \code{ref.calibration = FALSE}.
#' Documented in \code{\link[radiator]{calibrate_alleles}}.
#' \item \code{vcf.stats: } (optional, logical) Generates individuals and
#' markers statistics helpful for filtering.
#' These are very fast to generate and because computational
#' cost is minimal, even for huge VCFs, the default is \code{vcf.stats = TRUE}.
#' \item \code{vcf.metadata: } (optional, logical or character string)
#' With \code{vcf.metadata = FALSE}, only the genotypes are kept (GT field)
#' in the tidy dataset.
#' With \code{vcf.metadata = TRUE},
#' all the metadata contained in the \code{FORMAT} field will be kept in
#' the tidy data file. radiator is currently keeping and cleaning these metadata:
#' \code{"DP", "AD", "GL", "PL", "GQ", "HQ", "GOF", "NR", "NV"}.
#' e.g. you only want AD and PL, \code{vcf.metadata = c("AD", "PL")}.
#' Need another metadata ? Submit a request on github...
#' Default: \code{vcf.metadata = TRUE}.
#' }
#'
#' \strong{Filtering arguments:}
#' \enumerate{
#' \item \code{blacklist.id: } (optional, character)
#' Default (\code{blacklist.id = NULL}).
#' Documented in \code{\link[radiator]{tidy_genomic_data}}.
#' \item \code{filter.strands}: (optional, character)
#'  Default (\code{filter.strands = "blacklist"}).
#'  documented in \code{\link[radiator]{write_seqarray}}.
#' \item \code{whitelist.markers: }(optional, path)
#' Default: \code{whitelist.markers = NULL}.
#' Documented in \code{\link[radiator]{filter_whitelist}}.
#' \item \code{filter.individuals.missing}: (double)
#' Default: \code{filter.individuals.missing = NULL}.
#' Documented in \code{\link[radiator]{filter_individuals}}.
#' \item \code{filter.monomorphic}: (logical)
#' Default: \code{filter.monomorphic = TRUE}.
#' Documented in \code{\link[radiator]{filter_monomorphic}}.
#' Required package: \code{UpSetR}.
#' \item \code{filter.common.markers}: (logical)
#' Default: \code{filter.common.markers = TRUE}.
#' Documented in \code{\link[radiator]{filter_common_markers}}.
#' Required package: \code{UpSetR}.
#' \item \code{filter.mac}: (integer)
#' Default: \code{filter.mac = NULL}.
#' Documented in \code{\link[radiator]{filter_mac}}.
#' \item \code{filter.coverage}: (logical)
#' Default: \code{filter.coverage = NULL}.
#' Documented in \code{\link[radiator]{filter_coverage}}.
#' \item \code{filter.genotyping}: (integer)
#' Default: \code{filter.genotyping = NULL}.
#' Documented in \code{\link[radiator]{filter_genotyping}}.
#' \item \code{filter.snp.position.read: } (optional, character, integer)
#' Default: \code{filter.snp.position.read = NULL}.
#' Documented in \code{\link[radiator]{filter_snp_position_read}}.
#' \item \code{filter.snp.number: } (optional, character, integer)
#' Default: \code{filter.snp.number = NULL}.
#' Documented in \code{\link[radiator]{filter_snp_number}}.
#' \item \code{filter.short.ld}: (optional, character)
#' Default: \code{filter.short.ld = NULL}.
#' Documented in \code{\link[radiator]{filter_ld}}.
#' \item \code{filter.long.ld}: (optional, character)
#' Default: \code{filter.long.ld = NULL}.
#' Documented in \code{\link[radiator]{filter_ld}}.
#' Required package: \code{SNPRelate}.
#' \item \code{long.ld.missing}: Documented in \code{\link[radiator]{filter_ld}}.
#' Default: \code{long.ld.missing = FALSE}.
#' \item \code{ld.method}: Documented in \code{\link[radiator]{filter_ld}}.
#' Default: \code{ld.method = "r2"}.
#' }
#'

#'
#' @examples
#' \dontrun{
#' # very basic with built-in defaults (not recommended):
#' prep.data <- radiator::tidy_vcf(data = "populations.snps.vcf")
#'
#' # Using more arguments and filters (recommended):
#' tidy.data <- radiator::tidy_vcf(
#'     data = "populations.snps.vcf",
#'     strata = "strata_salamander.tsv",
#'     filter.individuals.missing = "outlier",
#'     filter.mac = 4,
#'     filter.genotyping = 0.1,
#'     filter.snp.position.read = "outliers",
#'     filter.short.ld = "mac",
#'     path.folder = "salamander/prep_data",
#'     verbose = TRUE)
#' }

#' @seealso \code{\link[radiator]{write_seqarray}},
#' \code{\link[radiator]{tidy_genomic_data}},
#' \code{\link[radiator]{filter_rad}}

#' @references Danecek P, Auton A, Abecasis G et al. (2011)
#' The variant call format and VCFtools.
#' Bioinformatics, 27, 2156-2158.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

tidy_vcf <- function(
  data,
  strata = NULL,
  filename = NULL,
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
  # filename = NULL

  # filter.individuals.missing = "outlier"
  # filter.common.markers = TRUE
  # filter.monomorphic = TRUE
  # filter.strands = FALSE
  # filter.mac = 4
  # filter.coverage = TRUE
  # filter.genotyping = 10
  # filter.snp.position.read = "outliers"
  # filter.short.ld = "maf"
  # filter.long.ld = 0.8
  # gt.vcf.nuc = TRUE
  # gt.vcf = TRUE
  # gt = TRUE
  # gt.bin = TRUE
  # wide = FALSE
  # ref.calibration = FALSE
  # keep.gds = TRUE

  # Cleanup---------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (verbose) message("Execution date@time: ", file.date)
  old.dir <- getwd()
  opt.change <- getOption("width")
  options(width = 70)
  timing <- proc.time()# for timing
  res <- list()
  #back to the original directory and options
  on.exit(setwd(old.dir), add = TRUE)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(timing <- proc.time() - timing, add = TRUE)
  on.exit(if (verbose) message("\nComputation time, overall: ", round(timing[[3]]), " sec"), add = TRUE)
  on.exit(if (verbose) cat("############################## completed tidy_vcf ##############################\n"), add = TRUE)


  # Required package -----------------------------------------------------------
  if (!"SeqVarTools" %in% utils::installed.packages()[,"Package"]) {
    rlang::abort('Please install SeqVarTools for this option:\n
                 install.packages("BiocManager")
                 BiocManager::install("SeqVarTools")')
  }

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("vcf file missing")

  # Function call and dotslist -------------------------------------------------
  rad.dots <- radiator_dots(
    func.name = as.list(sys.call())[[1]],
    fd = rlang::fn_fmls_names(),
    args.list = as.list(environment()),
    dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE),
    keepers = c(
      "blacklist.id",
      "filter.individuals.missing",
      "filter.individuals.coverage.total",
      "filter.common.markers",
      "filter.monomorphic",
      "filter.mac",
      "filter.coverage",
      "filter.genotyping",
      "filter.snp.position.read",
      "filter.snp.number",
      "filter.short.ld", "filter.long.ld", "long.ld.missing", "ld.method",
      "filter.strands",
      "random.seed",
      "path.folder",
      "parameters",
      "gt", "gt.bin", "gt.vcf", "gt.vcf.nuc",
      "ref.calibration",
      "vcf.metadata", "vcf.stats",
      "whitelist.markers",
      "internal",
      "tidy.check", "tidy.vcf"
    ),
    verbose = verbose
  )

  if (!is.null(filter.snp.position.read) ||
      !is.null(filter.mac) ||
      !is.null(filter.coverage) ||
      !is.null(filter.genotyping) ||
      !is.null(filter.short.ld) ||
      !is.null(filter.long.ld) ||
      !is.null(long.ld.missing) ||
      !is.null(filter.individuals.missing) ||
      !is.null(filter.individuals.coverage.total) ||
      !is.null(filter.snp.number)) subsample.markers.stats <- 1

  if (!is.null(ld.method)) {
    ld.method <- match.arg(ld.method, c("composite", "r", "r2", "dprime", "corr"))
  }
  if (!is.null(filter.snp.position.read)) {
    filter.snp.position.read <- match.arg(
      arg = filter.snp.position.read,
      choices = c("outliers", "iqr", "q75"),
      several.ok = TRUE)
  }
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

  if (is.null(random.seed)) {
    random.seed <- sample(x = 1:1000000, size = 1)
    set.seed(random.seed)
  } else {
    set.seed(random.seed)
  }

  # LD
  # currently: will only filter long ld if short ld as been taken care of first...
  if (!is.null(filter.long.ld) && is.null(filter.short.ld)) {
    message("\nfilter.short.ld argument set by default to: maf")
    filter.short.ld <- "mac"
  }

  if (is.null(tidy.vcf)) tidy.vcf <- TRUE
  if (is.null(tidy.check)) tidy.check <- TRUE

  # Folders---------------------------------------------------------------------
  wf <- path.folder <- generate_folder(
    f = path.folder,
    rad.folder = "tidy_vcf",
    prefix_int = FALSE,
    internal = internal,
    file.date = file.date,
    verbose = verbose)

  # import VCF -----------------------------------------------------------------
  data <- radiator::write_seqarray(
    data = data,
    strata = strata,
    filename = filename,
    verbose = FALSE,
    parallel.core = parallel.core,
    internal = TRUE,
    vcf.stats = vcf.stats,
    vcf.metadata = vcf.metadata,
    path.folder = path.folder,
    random.seed = random.seed,
    parameters = parameters,
    blacklist.id = blacklist.id,
    filter.strands = filter.strands,
    filter.monomorphic = filter.monomorphic,
    filter.common.markers = filter.common.markers,
    whitelist.markers = whitelist.markers,
    filter.individuals.missing = filter.individuals.missing,
    filter.individuals.coverage.total = filter.individuals.coverage.total,
    filter.mac = filter.mac,
    filter.coverage = filter.coverage,
    filter.genotyping = filter.genotyping,
    filter.snp.position.read = filter.snp.position.read,
    filter.snp.number = filter.snp.number,
    filter.short.ld = filter.short.ld,
    filter.long.ld = filter.long.ld,
    long.ld.missing = long.ld.missing,
    ld.method = ld.method
  )

  # tidy_vcf folder ------------------------------------------------------------
  tidy.folder <- generate_folder(
    f = path.folder,
    rad.folder = "tidy_vcf",
    prefix_int = TRUE,
    internal = FALSE,
    file.date = file.date,
    verbose = verbose)

  # write the dots file: after the GDS import...
  write_rad(
    data = rad.dots,
    path = tidy.folder,
    filename = stringi::stri_join("radiator_tidy_vcf_args_", file.date, ".tsv"),
    tsv = TRUE,
    internal = internal,
    verbose = verbose
  )

  # Random seed ----------------------------------------------------------------
  readr::write_lines(x = random.seed, path = file.path(tidy.folder, "random.seed"))
  if (verbose) message("File written: random.seed (", random.seed,")")

  # Tidy the data --------------------------------------------------------------
  if (tidy.vcf) {
    # Get info markers and individuals -----------------------------------------
    wl.v <- SeqArray::seqGetData(data, "variant.id")
    markers.meta <- extract_markers_metadata(gds = data) %>%
      dplyr::filter(VARIANT_ID %in% wl.v)
    wl.s <- SeqArray::seqGetData(data, "sample.id")

    individuals <- extract_individuals(
      gds = data, ind.field.select = "INDIVIDUALS"
    ) %>%
      dplyr::filter(INDIVIDUALS %in% wl.s) %$% INDIVIDUALS
    n.markers <- length(markers.meta$VARIANT_ID)
    n.individuals <- length(individuals)

    # STRATEGY tidy vcf  ---------------------------------------------------------
    # Depending on the number of markers ...
    # All this can be overwritten in ... argument
    # gt.bin is the dosage of ALT allele: 0, 1, 2 NA
    if (is.null(gt.bin)) {
      if (n.markers < 5000) gt.bin <- TRUE
      if (n.markers >= 5000 && n.markers < 30000) gt.bin <- TRUE
      if (n.markers >= 30000) gt.bin <- FALSE
    }
    # gt.vcf is genotype coding in the VCF: 0/0, 0/1, 1/1, ./.
    if (is.null(gt.vcf)) {
      if (n.markers < 5000) gt.vcf <- TRUE
      if (n.markers >= 5000 && n.markers < 30000) gt.vcf <- TRUE
      if (n.markers >= 30000) gt.vcf <- FALSE
    }
    # gt.vcf.nuc is genotype coding in the VCF but with nucleotides: A/C, ./.
    if (is.null(gt.vcf.nuc)) {
      if (n.markers < 5000) gt.vcf.nuc <- TRUE
      if (n.markers >= 5000 && n.markers < 30000) gt.vcf.nuc <- FALSE
      if (n.markers >= 30000) gt.vcf.nuc <- FALSE
    }
    # gt is genotype coding a la genepop: 001002, 000000
    if (is.null(gt)) {
      if (n.markers < 5000) gt <- TRUE
      if (n.markers >= 5000 && n.markers < 30000) gt <- FALSE
      if (n.markers >= 30000) gt <- FALSE
    }

    # check
    # gt.bin
    # gt.vcf
    # gt.vcf.nuc
    # gt

    # Tidy check ----------  ---------------------------------------------------
    if (tidy.check && n.markers > 20000) {
      cat("\n\n################################## IMPORTANT ###################################\n")
      message("Tidying vcf with ", n.markers, " SNPs is not optimal:")
      message("    1. a computer with lots of RAM is required")
      message("    2. it's very slow to generate")
      message("    3. it's very slow to run codes after")
      message("    4. for most non model species this number of markers is not realistic...")
      message("\nRecommendations:")
      message("    1. use advance features available in this function (read doc)")
      message("    2. filter your dataset. e.g. with filter_rad")
      message("\nIdeally target a maximum of ~ 10 000 - 20 0000 unlinked SNPs\n")

      if (n.markers > 20000) tidy.vcf <- FALSE
      tidy.vcf <- interactive_question(
        x = "\nContinue tidying the VCF (y/n) ?",
        answer.opt = c("Y", "N", "Yes", "No", "YES", "NO", "yes", "no", "y", "n"))
      if (any(c("y", "Y", "Yes", "YES", "yes") %in% tidy.vcf)) {
        tidy.vcf <- TRUE
        message("Tidying the large vcf...")
      } else {
        message("\nKeeping the GDS object/file")
        tidy.vcf <- FALSE
        return(data)
      }
    }

    # Print genotypes tidying
    message("\nGenotypes formats generated with ", n.markers, " SNPs: ")
    message("    GT_BIN (the dosage of ALT allele: 0, 1, 2 NA): ", gt.bin)
    message("    GT_VCF (the genotype coding VCFs: 0/0, 0/1, 1/1, ./.): ", gt.vcf)
    message("    GT_VCF_NUC (the genotype coding in VCFs, but with nucleotides: A/C, ./.): ", gt.vcf.nuc)
    message("    GT (the genotype coding 'a la genepop': 001002, 001001, 000000): ", gt)

    # Tidying TRUE  --------------------------------------------------------------
    if (tidy.vcf) {
      tibble.size <- n.markers * n.individuals
      tidy.data <- tibble::tibble(GT_BIN = numeric(tibble.size))
      # When IDs are blacklisted...you want to recalibrate REF/ALT alleles
      if (!is.null(blacklist.id)) {
        ref.calibration <- TRUE
        if (verbose) message("\nRe-calibration of REF/ALT alleles: TRUE")
      }

      # bi- or multi-alllelic VCF ------------------------------------------------
      biallelic <- detect_biallelic_markers(data = data)

      # import genotypes ---------------------------------------------------------
      # gt.bin
      if (gt.bin) {
        tidy.data %<>%
          dplyr::mutate(
            GT_BIN = as.vector(SeqArray::seqGetData(
              gdsfile = data, var.name = "$dosage_alt")
            )
          )
        # if (wide) {
        #   tidy.data$GT_BIN <- tibble::as_tibble(x = SeqArray::seqGetData(
        #     gdsfile = data, var.name = "$dosage_alt")) %>%
        #     magrittr::set_colnames(x = ., value = data$markers.meta$MARKERS) %>%
        #     tibble::add_column(.data = .,
        #                        "INDIVIDUALS" = id.string,
        #                        .before = 1)
        # }
      } else {
        tidy.data %<>% dplyr::select(-GT_BIN)
      }

      if (gt) gt.vcf.nuc <- TRUE
      # nucleotides info required to generate genepop format 001, 002, 003, 004
      # when gt is TRUE, gt.vcf.nuc is always TRUE
      if (gt.vcf.nuc) {
        tidy.data %<>%
          dplyr::mutate(
            GT_VCF_NUC = as.vector(SeqVarTools::getGenotypeAlleles(
              gdsobj = data, use.names = TRUE))
          )
      }

      if (gt.vcf) {
        tidy.data %<>%
          dplyr::mutate(
            GT_VCF = as.vector(SeqVarTools::getGenotype(gdsobj = data, use.names = TRUE))
          )
        tidy.data$GT_VCF[is.na(tidy.data$GT_VCF)] <- "./."
      } else {
        tidy.data %<>% dplyr::select(-GT_VCF)
      }

      if (gt) {
        tidy.data %<>%
          dplyr::mutate(
            GT = stringi::stri_replace_all_fixed(
              str = GT_VCF_NUC,
              pattern = c("A", "C", "G", "T", "/"),
              replacement = c("001", "002", "003", "004", ""),
              vectorize_all = FALSE) %>%
              stringi::stri_replace_na(str = ., replacement = "000000")
          )
      }
      # if (wide) {
      #   # id.string
      #   #   id.string <- data$individuals$INDIVIDUALS
      #
      #   data$tidy.data$GT <- tibble::as_tibble(
      #     matrix(data = data$tidy.data$GT,
      #            nrow = nrow(id.string),
      #            ncol = nrow(data$markers.meta))
      #   ) %>%
      #     magrittr::set_colnames(x = ., value = data$markers.meta$MARKERS) %>%
      #     tibble::add_column(
      #       .data = .,
      #       "INDIVIDUALS" = id.string,
      #       .before = 1
      #     )
      # }

      # replace NA in gt.vcf.nuc
      if (gt.vcf.nuc) {
        tidy.data$GT_VCF_NUC[is.na(tidy.data$GT_VCF_NUC)] <- "./."
      }
      # more work for gt.vcf.nuc
      # if (gt.vcf.nuc) {
      #   if (wide) {
      #     tidy.data$GT_VCF_NUC %<>%
      #       magrittr::set_colnames(x = ., value = markers.meta$MARKERS) %>%
      #       tibble::as_tibble(x = ., rownames = "INDIVIDUALS")
      #   } else {
      #     tidy.data$GT_VCF_NUC <- as.vector(tidy.data$GT_VCF_NUC)
      #   }
      #   tidy.data$GT_VCF_NUC[is.na(tidy.data$GT_VCF_NUC)] <- "./."
      # } else {
      #   tidy.data$GT_VCF_NUC <- NULL
      # }

      # gt.vcf
      # if (gt.vcf) {
      #   data$tidy.data$GT_VCF <- SeqVarTools::getGenotype(gdsobj = gds, use.names = TRUE)
      #   if (wide) {
      #     data$tidy.data$GT_VCF %<>% magrittr::set_colnames(x = ., value = data$markers.meta$MARKERS) %>%
      #       tibble::as_tibble(x = ., rownames = "INDIVIDUALS")
      #   } else {
      #     data$tidy.data$GT_VCF <- as.vector(data$tidy.data$GT_VCF)
      #   }
      #
      #   data$tidy.data$GT_VCF[is.na(data$tidy.data$GT_VCF)] <- "./."
      #   # replace(data, which(data == what), NA)
      # }

      # check missing genotypes
      # head(tidy.data$GT[is.na(tidy.data$GT_BIN)])
      # head(tidy.data$GT_VCF_NUC[is.na(tidy.data$GT_BIN)])
      # head(tidy.data$GT_VCF[is.na(tidy.data$GT_BIN)])

      # genotypes metadata ---------------------------------------------------------
      # Check vcf.metadata
      # vcf.metadata <- TRUE
      # vcf.metadata <- FALSE
      # vcf.metadata <- NULL
      # vcf.metadata <- "GT"

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
      #Check
      # vcf.metadata
      # overwrite.metadata


      if (vcf.metadata) {
        if (verbose) message("\nKeeping vcf genotypes metadata: yes")
        # detect FORMAT fields available
        have <-  SeqArray::seqSummary(
          gdsfile = data,
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
          tidy.data %<>%
            dplyr::bind_cols(
              genotypes.meta = purrr::map(
                .x = parse.format.list, .f = parse_gds_metadata, gds = data,
                verbose = verbose, parallel.core = parallel.core) %>%
                purrr::flatten(.) %>%
                purrr::flatten_df(.)
            )
          # check
          # test <- tidy.data
        } else {
          if (verbose) message("    genotypes metadata: none found")
          vcf.metadata <- FALSE
        }
      }

      # # Remove or not gds connection and file
      # if (!keep.gds) {
      #   gds <- NULL
      #   if (file.exists(data$filename)) file.remove(data$filename)
      # }

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
                        names(markers.meta))
      tidy.data <- suppressWarnings(
        dplyr::select(markers.meta, dplyr::one_of(want))) %>%
        dplyr::bind_cols(
          tibble::as_tibble(
            matrix(
              data = NA,
              nrow = n.markers, ncol = n.individuals)) %>%
            magrittr::set_colnames(x = ., value = individuals)) %>%
        data.table::as.data.table(.) %>%
        data.table::melt.data.table(
          data = .,
          id.vars = want,
          variable.name = "INDIVIDUALS",
          value.name = "GT",
          variable.factor = FALSE) %>%
        tibble::as_tibble(.) %>%
        dplyr::select(-GT) %>%
        dplyr::mutate(
          MARKERS = factor(x = MARKERS,
                           levels = markers.meta$MARKERS, ordered = TRUE),
          INDIVIDUALS = factor(x = INDIVIDUALS,
                               levels = individuals,
                               ordered = TRUE)) %>%
        dplyr::arrange(MARKERS, INDIVIDUALS) %>%
        dplyr::bind_cols(tidy.data)
      # Check stacks AD problem --------------------------------------------------
      # Some genotypes with missing AD...
      if (vcf.metadata) {
        stacks.ad.problem <- tidy.data %>%
          dplyr::select(READ_DEPTH, ALLELE_REF_DEPTH,ALLELE_ALT_DEPTH) %>%
          dplyr::filter(!is.na(READ_DEPTH)) %>%
          dplyr::filter(is.na(ALLELE_REF_DEPTH) & is.na(ALLELE_ALT_DEPTH))
        stacks.ad.problem.n <- nrow(stacks.ad.problem)

        if (stacks.ad.problem.n > 0) {
          message("\n\nStacks problem detected")
          message("    missing allele depth info")
          message("    number of genotypes with problem: ", stacks.ad.problem.n)
          message("    correcting problem by adding the read depth info into AD fields...\n\n")

          stacks.ad.problem <- dplyr::select(tidy.data, GT_BIN) %>%
            dplyr::bind_cols(dplyr::select(tidy.data, READ_DEPTH, ALLELE_REF_DEPTH, ALLELE_ALT_DEPTH)) %>%
            dplyr::mutate(COL_ID = seq(1, n())) %>%
            dplyr::filter(!is.na(READ_DEPTH)) %>%
            dplyr::filter(is.na(ALLELE_REF_DEPTH) & is.na(ALLELE_ALT_DEPTH)) %>%
            dplyr::mutate(
              ALLELE_REF_DEPTH = dplyr::if_else(GT_BIN == 0, READ_DEPTH, ALLELE_REF_DEPTH),
              ALLELE_ALT_DEPTH = dplyr::if_else(GT_BIN == 2, READ_DEPTH, ALLELE_ALT_DEPTH)
            )

          tidy.data$ALLELE_REF_DEPTH[stacks.ad.problem$COL_ID] <- stacks.ad.problem$ALLELE_REF_DEPTH
          tidy.data$ALLELE_ALT_DEPTH[stacks.ad.problem$COL_ID] <- stacks.ad.problem$ALLELE_ALT_DEPTH
          # check <- tidy.data[stacks.ad.problem$COL_ID, ]
        }
        stacks.ad.problem.n <- stacks.ad.problem.n <- NULL
      }

      # re-calibration of ref/alt alleles ------------------------------------------

      if (ref.calibration) {
        if (verbose) message("\nCalculating REF/ALT alleles...")
        tidy.data <- radiator::calibrate_alleles(
          data = tidy.data,
          biallelic = biallelic,
          parallel.core = parallel.core,
          verbose = verbose,
          gt.vcf.nuc = gt.vcf.nuc,
          gt = gt,
          gt.vcf = gt.vcf,
          gt.bin = gt.bin
        ) %$% input
      }

      # include strata ---------------------------------------------------------
      strata <- extract_individuals(gds = data, ind.field.select = c("STRATA", "INDIVIDUALS"))

      if (!is.null(strata)) {
        tidy.data <- join_strata(
          data = tidy.data,
          strata = strata,
          pop.id = TRUE,
          verbose = TRUE
          )
      } else {
        tidy.data %<>% dplyr::mutate(POP_ID = 1L)
      }


      # Re ordering columns
      want <- c("MARKERS", "CHROM", "LOCUS", "POS", "ID", "COL", "INDIVIDUALS",
                "STRATA", "POP_ID",
                "REF", "ALT", "GT_VCF", "GT_VCF_NUC", "GT", "GT_BIN",
                "POLYMORPHIC")

      suppressWarnings(
        tidy.data %<>% dplyr::select(dplyr::one_of(want), dplyr::everything())
        )

      # Sort id
      tidy.data %<>% dplyr::arrange(POP_ID, INDIVIDUALS)

      #write tidy
      # path.folder <- generate_folder(
      #   f = wf,
      #   rad.folder = "tidy_data",
      #   internal = FALSE,
      #   file.date = file.date,
      #   verbose = verbose)

      filename.rad <- generate_filename(
        path.folder = tidy.folder,
        extension = "rad")
      write_rad(data = tidy.data, path = filename.rad$filename)

      if (verbose) message("Updating GDS with genotypes.meta values")
      update_radiator_gds(gds = data, node.name = "genotypes.meta", value = tidy.data)
      message("\nTidy data file written: ", filename.rad$filename.short)
      } #tidy.vcf
  }#tidy.vcf


  # Results --------------------------------------------------------------------
  return(tidy.data)
}#End tidy_vcf


# Internal nested Function -----------------------------------------------------
#' @title parse_gds_metadata
#' @description function to parse the format field and tidy the results of VCF
#' @rdname parse_gds_metadata
#' @keywords internal
#' @export
parse_gds_metadata <- function(
  x,
  gds = NULL,
  verbose = TRUE,
  parallel.core = parallel::detectCores() - 1) {

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
  # format.name <- x <- "RO" # not yet implemented
  # format.name <- x <- "QR" # not yet implemented
  # format.name <- x <- "AO" # not yet implemented
  # format.name <- x <- "QA" # not yet implemented

  if (verbose) message("\nParsing and tidying: ", format.name)

  if (format.name == "AD") {
    if (verbose) message("AD column: splitting into ALLELE_REF_DEPTH and ALLELE_ALT_DEPTH")
    res$AD <- SeqArray::seqGetData(
      gdsfile = gds,
      var.name = "annotation/format/AD"
    )$data %>%
      tibble::as_tibble(.)
    column.vec <- seq_along(res$AD)
    res$AD <- tibble::tibble(
      ALLELE_REF_DEPTH = res$AD[, column.vec %% 2 == 1] %>%
        as.matrix(.) %>%
        as.vector(.),
      ALLELE_ALT_DEPTH = res$AD[, column.vec %% 2 == 0] %>%
        as.matrix(.) %>%
        as.vector(.))
    split.vec <- split_vec_row(x = res$AD, 3, parallel.core = parallel.core)
    res$AD$ALLELE_REF_DEPTH <- clean_ad(
      x = res$AD$ALLELE_REF_DEPTH,
      split.vec = split.vec,
      parallel.core = parallel.core)

    res$AD$ALLELE_ALT_DEPTH <- clean_ad(
      x = res$AD$ALLELE_ALT_DEPTH,
      split.vec = split.vec,
      parallel.core = parallel.core)
    split.vec <- NULL
    # test1 <- res$AD
  } # End AD


  # Read depth
  if (format.name == "DP") {
    if (verbose) message("DP column: cleaning and renaming to READ_DEPTH")
    res$DP <- tibble::tibble(READ_DEPTH = SeqArray::seqGetData(
      gdsfile = gds,
      var.name = "annotation/format/DP")$data %>% as.vector(.))
    # test <- res$DP
  } # End DP

  # Cleaning HQ: Haplotype quality as phred score
  if (format.name == "HQ") {
    res$HQ <- tibble::tibble(HQ = SeqArray::seqGetData(
      gdsfile = gds,
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
      gdsfile = gds,
      var.name = "annotation/format/GQ")$data %>% as.vector(.))
    # test <- res$GQ
  } # End GQ

  # GL cleaning
  if (format.name == "GL") {
    if (verbose) message("GL column: cleaning Genotype Likelihood column")
    res$GL <- SeqArray::seqGetData(gdsfile = gds,
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
      gdsfile = gds,
      var.name = "annotation/format/GOF")$data %>% as.vector(.))
    # test <- res$GOF
  } # End GOF

  # Cleaning NR: Number of reads covering variant location in this sample
  if (format.name == "NR") {
    if (verbose) message("NR column: splitting column into the number of variant")
    res$NR <- tibble::tibble(NR = SeqArray::seqGetData(
      gdsfile = gds,
      var.name = "annotation/format/NR")$data %>% as.vector(.))
    # test <- res$NR
  }#End cleaning NR column

  # Cleaning NV: Number of reads containing variant in this sample
  if (format.name == "NV") {
    if (verbose) message("NV column: splitting column into the number of variant")
    res$NR <- tibble::tibble(NV = SeqArray::seqGetData(
      gdsfile = gds,
      var.name = "annotation/format/NV")$data %>% as.vector(.))
    # test <- res$NV
  }#End cleaning NV column

  # RO
  # Cleaning RO: Reference allele observation count
  # if (format.name == "RO") {
  #   if (verbose) message("RO column: splitting column into the number of variant")
  #   res$RO <- tibble::tibble(RO = SeqArray::seqGetData(
  #     gdsfile = gds,
  #     var.name = "annotation/format/RO")$data %>% as.vector(.))
  #   # test <- res$RO
  # }#End cleaning RO column

  # # Cleaning QR: Sum of quality of the reference observations
  # if (format.name == "QR") {
  #   if (verbose) message("QR column: splitting column into the number of variant")
  #   res$QR <- tibble::tibble(QR = SeqArray::seqGetData(
  #     gdsfile = gds,
  #     var.name = "annotation/format/QR")$data %>% as.vector(.))
  #   # test <- res$QR
  # }#End cleaning QR column


  # # Cleaning AO: Alternate allele observation count
  # if (format.name == "AO") {
  #   if (verbose) message("AO column: splitting column into the number of variant")
  #   res$AO <- tibble::tibble(AO = SeqArray::seqGetData(
  #     gdsfile = gds,
  #     var.name = "annotation/format/AO")$data %>% as.vector(.))
  #   # test <- res$AO
  # }#End cleaning AO column

  # # Cleaning QA: Sum of quality of the alternate observations
  # if (format.name == "QA") {
  #   if (verbose) message("QA column: splitting column into the number of variant")
  #   res$QA <- tibble::tibble(QA = SeqArray::seqGetData(
  #     gdsfile = gds,
  #     var.name = "annotation/format/QA")$data %>% as.vector(.))
  #   # test <- res$QA
  # }#End cleaning QA column



  # test <- res$AD %>%
  #   dplyr::bind_cols(READ_DEPTH = SeqArray::seqGetData(
  #   gdsfile = gds,
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
