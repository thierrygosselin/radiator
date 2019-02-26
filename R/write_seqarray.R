# write a SeqArray object from a tidy data frame

#' @name write_seqarray
#' @title Write a SeqArray GDS file from a vcf file and generate a connection object.
#' @description Write a SeqArray \href{http://zhengxwen.github.io/SeqArray/}{SeqArray}
#' object of class \code{SeqVarGDSClass} (Zheng et al. 2017)
#' from a vcf file and generate a connection object with the
#' Genomic Data Structure (GDS) file format
#' \href{http://zhengxwen.github.io/gdsfmt}{gdsfmt}
#'
#' The function as an advance mode that allows various filtering arguments to
#' be used. Handy to prune the dataset based on various QCs but also to
#' remove artifacts, duplicated information and thus lower the overall noise in
#' the VCF.
#'
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.


#' @param data (path, character) The VCF markers are biallelic SNPs or haplotypes.
#' To make the VCF population-ready, you have to use \code{strata} argument.
#'
#' \itemize{
#' \item \strong{GATK, platypus and ipyrad VCFs:}
#' Some VCFs have an \code{ID} column filled with \code{.},
#' the LOCUS information is all contained along the linkage group in the
#' \code{CHROM} column. Consequently, the short read locus information is unknown.
#' To make it work with
#' \href{https://github.com/thierrygosselin/radiator}{radiator},
#' the \code{ID} column is filled with the \code{POS} column info.
#' \item \strong{stacks VCFs:} with \emph{de novo} approaches, the CHROM column is
#' filled with "1", the LOCUS column correspond to the CHROM section in stacks VCF and
#' the COL column is POS -1. With a reference genome, the ID column in stacks VCF is
#' separated into "LOCUS", "COL", "STRANDS".
#' }


#' @param filename (optional) The file name of the Genomic Data Structure (GDS) file.
#' radiator will append \code{.gds.rad} to the filename.
#' If the filename chosen exists in the working directory,
#' the default \code{radiator_datetime.gds} is chosen.
#' Default: \code{filename = NULL}.

#' @param vcf.stats (optional, logical) Generates individuals and markers
#' statistics helpful for filtering.
#' These are very fast to generate and because computational
#' cost is minimal, even for huge VCFs, the default is \code{vcf.stats = TRUE}.
#' Individual's missing genotype proportion, averaged heterozygosity,
#' total coverage, mean genotype coverage and marker's metadata
#' along count for ref and alt alleles and mean
#' coverage is generated and written in the working directory.
#' Default: \code{vcf.stats = TRUE}.


#' @inheritParams tidy_genomic_data
#' @inheritParams radiator_common_arguments


#' @details
#' A vcf file of 35 GB with ~4 millions SNPs take about ~7 min with 8 CPU.
#' A vcf file of 21 GB with ~2 millions SNPs take about ~5 min with 7 CPU.
#'
#' After the file is generated, you can close your computer and
#' come back to it a month later and it's now a matter of sec to open a connection.
#'
#'
#' \strong{markers heterozygosity and inbreeding coefficient}:
#' The heterozygosity statistics (het obs and Fis) presented in this function
#' are calculated globally across strata, without the possibility to filter out
#' markers. We think this filtering for these statistics are best
#' addressed in \code{\link{filter_het}}, \code{\link{filter_fis}} and
#' \code{\link{filter_hwe}}, after outlier individuals and
#' markers (for other stats) are blacklisted. Why ? The reason is that an excess of
#' heterozygotes and negative Fis values are not a bad thing \emph{per se}.
#' Genomic regions under balancing selection may contain such markers and
#' statistics.
#'
#'
#' \strong{Advance mode, using \emph{dots-dots-dots ...}}
#' \enumerate{
#' \item \code{whitelist.markers}: detailed in \code{\link[radiator]{filter_whitelist}}.
#' \item \code{blacklist.id}: detailed in \code{\link[radiator]{tidy_genomic_data}}.
#' \item \code{pop.select}: detailed in \code{\link[radiator]{tidy_genomic_data}}.
#' \item \code{pop.levels}: detailed in \code{\link[radiator]{tidy_genomic_data}}.
#' \item \code{pop.labels}: detailed in \code{\link[radiator]{tidy_genomic_data}}.
#' \item \code{filter.strands}: (optional, character) Filter duplicate SNPs
#' found on different strands (+/-), options are:
#' \code{"keep.both", "best.stats", "blacklist"}. \code{"keep.both"}: does nothing
#' and duplicated markers are kept (not recommended, but here for testing purposes),
#' \code{"best.stats"}: will keep only one, based on the best statistics
#' (MAC and missingness). \code{"blacklist"}: discard all duplicated markers.
#' Default (\code{filter.strands = "blacklist"}).
#' \item \code{filter.common.markers}: (logical) Default: \code{filter.common.markers = TRUE}.
#' Documented in \code{\link{filter_common_markers}}.
#' \item \code{filter.mac}: (integer) Default: code{filter.mac = NULL}.
#' To blacklist markers below a specific Minor Allele Count (calculated overall/global).
#' Documented in \code{\link{filter_mac}}.
#'
#' \item \code{filter.coverage}: (optional, string)
#' Default: code{filter.coverage = NULL}. To blacklist markers based on mean coverage.
#' Documented in \code{\link{filter_coverage}}.
#'
#' \item \code{filter.genotyping}: (integer) To blacklist markers with too
#' many missing data. e.g. \code{filter.genotyping = 10}, will only keep
#' markers with missing rate <= to 10 percent.
#' \item \code{filter.snp.position.read}: 3 options available, \code{"outliers", "q75", "iqr"}.
#' This argument will blacklist markers based on it's position on the read.
#' \code{filter.snp.position.read = "outliers"}, will remove markers based
#' on outlier statistics of the position on the reads. e.g. if majority of SNPs
#' are found between 10 and 90 pb, and very few above and below, those markers are
#' discarded. Use this function argument with dataset with problematic assembly,
#' or \emph{de novo} assembly with undocumented or poorly selected
#' mismatch threshold.
#' \item \code{filter.short.ld}: Described in \code{\link[radiator]{filter_ld}}
#' \item \code{filter.long.ld}: Described in  \code{\link[radiator]{filter_ld}}
#' \item \code{long.ld.missing}: Described in \code{\link[radiator]{filter_ld}}
#' \item \code{ld.method}: (optional, character) The values available are
#' \code{"composite"}, for LD composite measure, \code{"r"} for R coefficient
#' (by EM algorithm assuming HWE, it could be negative), \code{"r2"} for r^2,
#' \code{"dprime"} for D',
#' \code{"corr"} for correlation coefficient. The method corr and composite are
#' equivalent when SNPs are coded based on the presence of the alternate allele
#' (\code{0, 1, 2}).
#' Default: \code{ld.method = "r2"}.

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
#' \item \code{random.seed}: (integer, optional) For reproducibility, set an integer
#' that will be used inside codes that uses randomness. With default,
#' a random number is generated, printed and written in the directory.
#' Default: \code{random.seed = NULL}.
#' \item \code{subsample.markers.stats}: By default, when no filters are
#' requested and that the number of markers is > 200K,
#' 0.2 of markers are randomly selected to generate the
#' statistics (individuals and markers). This is an all-around and
#' reliable number.
#' In doubt, overwrite this value by using 1 (all markers selected) and
#' expect a small computational cost.
#' }

#' @return
#' The function returns a list with:
#' \enumerate{
#' \item \code{vcf.connection}: the name of the GDS file.
#' To close the connection with the GDS file: use \code{SeqArray::seqClose(OBJECT_NAME$vcf.connection)}
#' \item \code{individuals}: a tibble with the names of samples in the original
#' VCF \code{INDIVIDUALS_VCF} and corrected for use in radiator \code{INDIVIDUALS}
#' \item \code{biallelic}: is the data biallelic or not.
#' \item \code{markers.meta}: a tibble with the markers metadata, including:
#' \code{VARIANT_ID, CHROM, LOCUS, POS, COL, MARKERS, REF, ALT} when available.
#' \item \code{n.ind}: the number of individuals
#' \item \code{n.markers}: the number of markers
#' \item \code{n.chromosome}: the numer of chromosome
#' \item \code{n.locus}: the number of locus
#' \item \code{filename}: the name of the file used.
#' }


#' @export
#' @rdname write_seqarray
#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join
#' @importFrom stringi stri_replace_all_fixed
#' @importFrom tibble has_name
#' @importFrom tidyr spread

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

#' @references Zheng X, Gogarten S, Lawrence M, Stilp A, Conomos M, Weir BS,
#' Laurie C, Levine D (2017). SeqArray -- A storage-efficient high-performance
#' data format for WGS variant calls.
#' Bioinformatics.
#'
#' @references Danecek P, Auton A, Abecasis G et al. (2011)
#' The variant call format and VCFtools.
#' Bioinformatics, 27, 2156-2158.

#' @seealso \code{\link{filter_ld}} and \code{\link[radiator]{tidy_genomic_data}}

#' @examples
#' \dontrun{
#' #require(SeqVarTools)
#' # with built-in defaults:
#'  prep.data <- radiator::write_seqarray(data = "populations.snps.vcf")
#'
#' # Using more arguments and filters (recommended):
#' prep.data <- radiator::write_seqarray(
#'     data = "populations.snps.vcf",
#'     strata = "strata_salamander.tsv",
#'     vcf.stats = TRUE,
#'     filter.individuals.missing = "outlier",
#'     filter.common.markers = TRUE,
#'     filter.strands = "blacklist",
#'     filter.mac = 4,
#'     filter.genotyping = 50,
#'     filter.snp.position.read = "outliers",
#'     filter.short.ld = "maf",
#'     filter.long.ld = NULL,
#'     vcf.metadata = TRUE,
#'     path.folder = "salamander/prep_data",
#'     verbose = TRUE)
#' }

write_seqarray <- function(
  data,
  strata = NULL,
  filename = NULL,
  vcf.stats = TRUE,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
  ...
) {

  # #Test
  # data = "populations.snps.vcf"
  # strata <- NULL
  # filename <- NULL
  # parallel.core <- parallel::detectCores() - 1
  # verbose <- TRUE
  # path.folder = NULL
  # random.seed <- NULL

  # vcf.stats <- TRUE
  # vcf.metadata = TRUE
  # filter.strands = "blacklist"
  # filter.coverage = NULL
  # filter.genotyping <- NULL
  # filter.snp.position.read <- NULL
  # filter.mac <- NULL
  # filter.common.markers = FALSE
  # filter.short.ld <- NULL
  # filter.long.ld <- NULL
  # long.ld.missing <- TRUE
  # ld.method <- "r2"
  # blacklist.id = NULL
  # pop.select = NULL
  # pop.levels = NULL
  # pop.labels = NULL
  # whitelist.markers = NULL
  # keep.gds <- TRUE
  # markers.info = NULL
  # subsample.markers.stats = 0.2
  # parameters <- NULL
  # filter.individuals.missing <- NULL
  # filter.individuals.coverage.total <- NULL


  ## With filters
  # filter.individuals.missing <- "outliers"
  # filter.individuals.coverage.total <- "outliers"
  # filter.short.ld <- "mac"
  # filter.long.ld <- 0.3
  # filter.common.markers = TRUE
  # filter.mac <- 4
  # filter.coverage = c(5, 150)
  # filter.genotyping <- 0.15
  # filter.snp.position.read <- "outliers"
  # filter.snp.number <- "outliers"

  # Cleanup---------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (verbose) message("Execution date/time: ", file.date)
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
  on.exit(if (verbose) cat("########################### completed write_seqarray ###########################\n"), add = TRUE)

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
    keepers = c("whitelist.markers",
                "filter.mac", "filter.snp.position.read", "filter.snp.number",
                "filter.coverage", "filter.genotyping", "filter.short.ld",
                "filter.long.ld", "long.ld.missing", "ld.method",
                "filter.individuals.missing", "filter.individuals.coverage.total",
                "filter.common.markers", "filter.monomorphic",
                "filter.strands",
                "blacklist.id", "pop.select", "pop.levels", "pop.labels", "keep.gds",
                "path.folder",
                "markers.info", "vcf.metadata",
                "subsample.markers.stats", "parameters", "random.seed", "internal"),
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

  # Folders---------------------------------------------------------------------
  path.folder <- generate_folder(
    f = path.folder,
    rad.folder = "radiator",
    internal = internal,
    file.date = file.date,
    verbose = verbose)

  # write the dots file
  write_rad(
    data = rad.dots,
    path = path.folder,
    filename = stringi::stri_join("radiator_write_seqarray_args_", file.date, ".tsv"),
    tsv = TRUE,
    internal = internal,
    verbose = verbose
  )


  if (is.null(filename)) {
    ind.file <- stringi::stri_join("vcf_individuals_info_", file.date, ".tsv")
    markers.file <- stringi::stri_join("vcf_markers_metadata_", file.date, ".tsv")
    filename <- stringi::stri_join("radiator_", file.date, ".gds")
    blacklist.markers <- stringi::stri_join("blacklist.markers_", file.date, ".tsv")
    blacklist.id.filename <- stringi::stri_join("blacklist.individuals_", file.date, ".tsv")
  } else {
    ind.file <- stringi::stri_join(filename, "_vcf_individuals_info_", file.date, ".tsv")
    markers.file <- stringi::stri_join(filename, "_vcf_markers_metadata_", file.date, ".tsv")
    blacklist.markers <- stringi::stri_join(filename, "_blacklist.markers_", file.date, ".tsv")
    blacklist.id.filename <- stringi::stri_join(filename, "_blacklist.individuals_", file.date, ".tsv")
    filename.problem <- file.exists(stringi::stri_join(filename, ".gds"))
    if (filename.problem) {
      filename <- stringi::stri_join(filename, "_", file.date, ".gds")
    } else {
      filename <- stringi::stri_join(filename, ".gds")
    }
    filename.problem <- file.exists(filename)
    if (filename.problem) {
      filename <- stringi::stri_join("duplicated_", filename)
    }
  }

  filename.short <- filename
  filename <- file.path(path.folder, filename)

  # Random seed ----------------------------------------------------------------
  readr::write_lines(x = random.seed, path = file.path(path.folder, "random.seed"))
  if (verbose) message("File written: random.seed (", random.seed,")")

  # Filter parameter file: generate --------------------------------------------
  filters.parameters <- radiator_parameters(
    generate = TRUE,
    initiate = FALSE,
    update = FALSE,
    parameter.obj = parameters,
    path.folder = path.folder,
    file.date = file.date,
    verbose = verbose)

  # Read vcf -------------------------------------------------------------------
  timing.vcf <- proc.time()

  # Get file size
  big.vcf <- file.size(data)
  if (big.vcf <= 10000000) {
    message("\nReading VCF")
  }

  if (big.vcf > 10000000 &&  big.vcf < 500000000) {
    message("\nReading VCF...you have time for a espresso...")
  }

  if (big.vcf >= 500000000) {
    message("\nReading a large VCF...you actually have time for coffee, take your time!\n")
  }
  # Check for bad header generated by stacks
  detect.source <- check_header_source_vcf(data)
  stacks.2 <- detect.source$source
  source <- stringi::stri_replace_all_regex(
    str = detect.source$check.header$header$value
    [detect.source$check.header$header$id == "source"],
    pattern = "[\"]", replacement = "", vectorize_all = FALSE)
  check.header <- detect.source$check.header
  dp <- "DP" %in% detect.source$check.header$format$ID # Check that DP is trere

  if (!is.null(detect.source$markers.info)) markers.info <- detect.source$markers.info
  if (!is.null(detect.source$overwrite.metadata)) overwrite.metadata <- detect.source$overwrite.metadata

  # for small VCF SeqArray not that good with parallel...
  parallel.temp <- parallel.core
  if (big.vcf < 10000000) parallel.temp <- 1
  big.vcf <- NULL
  gds <- SeqArray::seqVCF2GDS(
    vcf.fn = data,
    out.fn = filename,
    parallel = parallel.temp,
    storage.option = "ZIP_RA",
    verbose = FALSE,
    header = check.header,
    info.import = markers.info, # characters, the variable name(s) in the INFO field for import; or NULL for all variables
    fmt.import = overwrite.metadata# characters, the variable name(s) in the FORMAT field for import; or NULL for all variables#
  ) %>%
    SeqArray::seqOpen(gds.fn = ., readonly = FALSE)

  parallel.temp <- check.header <- detect.source <- NULL #no longer used

  # Summary --------------------------------------------------------------------
  sum <- summary_gds(gds, verbose = verbose)
  n.ind <- sum$n.ind
  if (verbose && keep.gds) message("\nFile written: ", filename.short)
  if (verbose) message("Import timing: ", round((proc.time() - timing.vcf)[[3]]), " sec\n\n")

  if (verbose)  message("Analyzing the data...")
  # radiator skeleton folder ---------------------------------------------------
  radiator.gds <- radiator_gds_skeleton(gds)

  # source
  update_radiator_gds(gds = gds, node.name = "source", value = source)

  # Blacklist of individuals: generate
  bl.i <- update_bl_individuals(gds = gds, generate = TRUE)
  # Blacklist of markers: generate
  bl.gds <- update_bl_markers(gds = gds, generate = TRUE)

  # bi- or multi-alllelic VCF --------------------------------------------------
  biallelic <- detect_biallelic_markers(data = gds, verbose = verbose)

  # clean sample id in VCF -----------------------------------------------------
  individuals.vcf <- tibble::tibble(
    INDIVIDUALS_VCF = SeqArray::seqGetData(gds, "sample.id")) %>%
    dplyr::mutate(INDIVIDUALS_CLEAN = radiator::clean_ind_names(INDIVIDUALS_VCF))

  if (!identical(individuals.vcf$INDIVIDUALS_VCF, individuals.vcf$INDIVIDUALS_CLEAN)) {
    if (verbose) message("Cleaning VCF's sample names")
    clean.id.filename <- stringi::stri_join("cleaned.vcf.id.info_", file.date, ".tsv")
    readr::write_tsv(x = individuals.vcf,
                     path = file.path(path.folder, clean.id.filename))
    gdsfmt::add.gdsn(
      node = radiator.gds,
      name = "id.clean",
      val = individuals.vcf,
      replace = TRUE,
      compress = "ZIP_RA",
      closezip = TRUE)
  }

  # replace id in VCF
  gdsfmt::add.gdsn(
    node = gds,
    name = "sample.id",
    val = individuals.vcf$INDIVIDUALS_CLEAN,
    replace = TRUE,
    compress = "ZIP_RA",
    closezip = TRUE)


  individuals <- dplyr::select(individuals.vcf, INDIVIDUALS = INDIVIDUALS_CLEAN)

  # Add a individuals node
  update_radiator_gds(gds = gds, node.name = "individuals", value = individuals)

  # Sync id in STRATA and VCF --------------------------------------------------
  # Strata ---------------------------------------------------------------------
  strata <- radiator::read_strata(
    strata = strata,
    pop.levels = pop.levels,
    pop.labels = pop.labels,
    pop.select = pop.select,
    blacklist.id = blacklist.id) %$% strata

  if (!is.null(strata)) {
    if (verbose) message("Synchronizing samples in VCF and strata...")

    strata %<>% dplyr::filter(INDIVIDUALS %in% individuals.vcf$INDIVIDUALS_CLEAN)
    bl <- individuals.vcf %>%
      dplyr::filter(!INDIVIDUALS_CLEAN %in% strata$INDIVIDUALS) %>%
      dplyr::select(INDIVIDUALS = INDIVIDUALS_CLEAN) %>%
      dplyr::distinct(INDIVIDUALS) %>%
      dplyr::mutate(FILTER = "blacklisted by strata")
    # strata %<>% dplyr::filter(INDIVIDUALS %in% individuals.vcf$INDIVIDUALS_CLEAN)
    if (nrow(strata) == 0) {
      rlang::abort("No more individuals in your data, check VCF and strata ID names...")
    }
    individuals.vcf <- NULL
    blacklist.strata <- nrow(bl)
    if (blacklist.strata != 0) {
      bl.i <- update_bl_individuals(gds = radiator.gds, update = bl, bl.i.gds = bl.i)
      if (verbose) message("    number of sample blacklisted by the stata: ", blacklist.strata)
      # the param file is updated after markers metadata below
    }
    sync_gds(gds = gds, samples = strata$INDIVIDUALS)
    individuals <- strata
    strata <- TRUE
  } else {
    strata <- FALSE
    blacklist.strata <- 0L
    individuals %<>% dplyr::mutate(STRATA = 1L)
  }

  #Update GDS node
  update_radiator_gds(gds = gds, node.name = "individuals", value = individuals)

  # Markers metadata  ----------------------------------------------------------
  markers.meta <- extract_markers_metadata(
    gds = gds, radiator.node = FALSE)

  # reference genome or de novo ------------------------------------------------
  # system.time(ref.genome <- detect_ref_genome(chromosome = markers.meta$CHROM, verbose = verbose))
  ref.genome <- detect_ref_genome(data = gds, verbose = verbose)

  # Stacks specific adjustments
  if (!ref.genome) {
    if (stacks.2) {
      markers.meta %<>% dplyr::mutate(
        LOCUS = CHROM,
        CHROM = "1",
        COL = POS - 1)
    } else {
      markers.meta %<>% dplyr::mutate(
        CHROM = "1")
    }
  }

  # GATK, platypus and freebayes specific adjustment
  # Locus with NA or . or ""
  weird.locus <- length(unique(markers.meta$LOCUS)) <= 1
  if (weird.locus && !stacks.2) {
    if (verbose) message("LOCUS field empty... adding unique id instead")
    markers.meta$LOCUS <- markers.meta$VARIANT_ID
  }

  # LOCUS cleaning and Strands detection ---------------------------------------
  if (stringi::stri_detect_regex(str = markers.meta[1,3], pattern = "[^[:alnum:]]+")) {
    locus.sep <- unique(stringi::stri_extract_all_regex(
      str = markers.meta[1,3],
      pattern = "[^a-zA-Z0-9-+-]",
      omit_no_match = TRUE)[[1]])

    markers.meta <- suppressWarnings(
      tidyr::separate(
        data = markers.meta,
        col = LOCUS, into = c("LOCUS", "COL", "STRANDS"),
        sep = locus.sep,
        # sep = "",
        extra = "drop", fill = "warn",
        remove = TRUE, convert = TRUE)
    )
    if (anyNA(markers.meta$STRANDS)) markers.meta$STRANDS <- NA_character_
    locus.sep <- NULL

    detect.strand <- any(stringi::stri_detect_fixed(str = unique(markers.meta$STRANDS), pattern = "+"))
    if (anyNA(detect.strand)) detect.strand <- FALSE

    if (detect.strand) {
      blacklist.strands <- dplyr::distinct(markers.meta, CHROM, LOCUS, POS) %>%
        dplyr::group_by(CHROM, POS) %>%
        dplyr::tally(.) %>%
        dplyr::filter(n > 1) %>%
        dplyr::ungroup(.) %>%
        dplyr::select(CHROM, POS)

      if (nrow(blacklist.strands) == 0) {
        blacklist.strands <- NULL
      } else {
        blacklist.strands %<>%
          dplyr::mutate_at(.tbl = .,
                           .vars = c("CHROM", "POS"),
                           .funs = radiator::clean_markers_names)
        if (verbose) {
          message("\nDetected ", nrow(blacklist.strands)," duplicate SNPs on different strands (+/-)")
          message("    By default radiator prune those SNPs")
          message("    To change this behavior, use argument: filter.strands\n")
        }
      }

    } else {
      markers.meta %<>% dplyr::select(-STRANDS)
      blacklist.strands <- NULL
    }
    detect.strand <- NULL
  } else {
    blacklist.strands <- NULL
  }

  # Generate MARKERS column and fix types --------------------------------------
  markers.meta %<>%
    dplyr::mutate_at(
      .tbl = .,
      .vars = c("CHROM", "LOCUS", "POS"),
      .funs = radiator::clean_markers_names) %>%
    dplyr::mutate(
      MARKERS = stringi::stri_join(CHROM, LOCUS, POS, sep = "__"),
      REF = SeqArray::seqGetData(gdsfile = gds, var.name = "$ref"),
      ALT = SeqArray::seqGetData(gdsfile = gds, var.name = "$alt")
    )

  # ADD MARKERS META to GDS
  update_radiator_gds(gds = gds, node.name = "markers.meta", value = markers.meta)
  # suppressWarnings(gdsfmt::add.gdsn(
  #   node = radiator.gds,
  #   name = "markers.meta",
  #   val = markers.meta,
  #   replace = TRUE,
  #   compress = "ZIP_RA",
  #   closezip = TRUE))

  # replace chromosome info in GDS
  # Why ? well snp ld e.g. will otherwise be performed by chromosome and with de novo data = by locus...
  gdsfmt::add.gdsn(
    node = gds,
    name = "chromosome",
    val = markers.meta$CHROM,
    replace = TRUE,
    compress = "ZIP_RA",
    closezip = TRUE)

  # radiator_parameters: initiate --------------------------------------------
  # with original VCF's values
  filters.parameters <- radiator_parameters(
    generate = FALSE,
    initiate = TRUE,
    update = TRUE,
    parameter.obj = filters.parameters,
    data = gds,
    filter.name = "vcf",
    param.name = "original values in VCF + strata",
    values = "",
    path.folder = path.folder,
    file.date = file.date,
    verbose = FALSE)

  # MONOMORPHIC MARKERS --------------------------------------------------
  gds <- filter_monomorphic(
    data = gds,
    filter.monomorphic = filter.monomorphic,
    parallel.core = parallel.core,
    verbose = verbose,
    parameters = filters.parameters,
    path.folder = path.folder,
    internal = internal)

  # PRE-FILTERING --------------------------------------------------------------
  # whitelist of markers, the filter column in the vcf and blacklisted strands
  # The order here doesnt matter

  # Whitelist------------------------------------------------------------
  gds <- filter_whitelist(
    data = gds,
    whitelist.markers = whitelist.markers,
    verbose = verbose,
    path.folder = path.folder,
    parameters = filters.parameters,
    biallelic = biallelic,
    markers.meta = markers.meta)

  # vcf FILTER column ---------------------------------------
  # This is automatically triggered
  markers.meta <- extract_markers_metadata(gds = gds)
  markers.meta$FILTER <- SeqArray::seqGetData(
    gds, "annotation/filter")
  filter.check.unique <- unique(markers.meta$FILTER)

  if (length(filter.check.unique) > 1) {
    message("Filtering markers based on VCF FILTER column")
    n.markers.before <- nrow(markers.meta)

    bl <- markers.meta %>%
      dplyr::filter(FILTER != "PASS") %>%
      dplyr::select(MARKERS) %>%
      dplyr::mutate(FILTER = "vcf filter column")

    bl.gds <- update_bl_markers(gds = radiator.gds, update = bl)

    markers.meta %<>% dplyr::filter(FILTER == "PASS")
    n.markers.after <- nrow(markers.meta)
    n.markers <- stringi::stri_join(n.markers.before, n.markers.before - n.markers.after, n.markers.after, sep = " / ")
    if (verbose) message("    Number of SNPs before / blacklisted / after: ", n.markers)
    sync_gds(gds = gds, markers = markers.meta$VARIANT_ID)

    # Filter parameter file: update
    filters.parameters <- radiator_parameters(
      generate = FALSE,
      initiate = FALSE,
      update = TRUE,
      parameter.obj = filters.parameters,
      data = gds,
      filter.name = "vcf filter column",
      param.name = "PASS or not",
      path.folder = path.folder,
      file.date = file.date,
      verbose = verbose)

  }
  filter.check.unique <- NULL
  markers.meta %<>% dplyr::select(-FILTER)

  # FILTER STRANDS -----------------------------------------------------------
  # This is automatically triggered filter.strands = "blacklist" by default
  if (!is.null(blacklist.strands)) {
    if (verbose) message("\nFiltering duplicated markers on different strands")
    blacklist.strands <- dplyr::right_join(markers.meta, blacklist.strands,
                                           by = c("CHROM", "POS")) %>%
      dplyr::select(VARIANT_ID, MARKERS, CHROM, POS)

    if (filter.strands == "keep.both") {
      message("Keeping both duplicated markers")
    }

    if (filter.strands == "best.stats") {
      sync_gds(gds = gds, markers = blacklist.strands$VARIANT_ID)
      blacklist.strands <- SeqArray::seqAlleleCount(
        gdsfile = gds,
        ref.allele = NULL,
        .progress = TRUE,
        parallel = parallel.core) %>%
        unlist(.) %>%
        matrix(
          data = .,
          nrow = nrow(blacklist.strands), ncol = 2, byrow = TRUE,
          dimnames = list(rownames = blacklist.strands$VARIANT_ID,
                          colnames = c("REF_COUNT", "ALT_COUNT"))) %>%
        tibble::as_tibble(., rownames = "VARIANT_ID") %>%
        tibble::add_column(.data = .,
                           MARKERS = blacklist.strands$MARKERS,
                           CHROM = blacklist.strands$CHROM,
                           POS = blacklist.strands$POS,
                           MISSING_PROP = SeqArray::seqMissing(
                             gdsfile = gds,
                             per.variant = TRUE, .progress = TRUE, parallel = parallel.core),
                           .after = 1) %>%
        dplyr::mutate(
          MAC = dplyr::if_else(ALT_COUNT < REF_COUNT, ALT_COUNT, REF_COUNT),
          ALT_COUNT = NULL, REF_COUNT = NULL) %>%
        dplyr::group_by(CHROM, POS) %>%
        dplyr::filter(MISSING_PROP == max(MISSING_PROP)) %>%
        dplyr::filter(MAC == min(MAC)) %>%
        dplyr::ungroup(.) %>%
        dplyr::arrange(CHROM, POS) %>%
        dplyr::distinct(CHROM, POS, .keep_all = TRUE) %>%
        dplyr::select(MARKERS)
    }# End best stats

    if (filter.strands == "blacklist") {
      blacklist.strands %<>% dplyr::distinct(MARKERS)
    }

    if (filter.strands != "keep.both") {
      message("Blacklisted markers: ", nrow(blacklist.strands))

      # update the blacklist
      blacklist.strands %<>%
        dplyr::mutate(FILTER = "filter.strands")

      bl.gds <- update_bl_markers(gds = radiator.gds, update = blacklist.strands)
      # Update markers.meta
      markers.meta %<>% dplyr::filter(!MARKERS %in% blacklist.strands$MARKERS)
      # check <- markers.meta

      # Update GDS
      update_radiator_gds(gds = gds, node.name = "markers.meta", value = markers.meta, sync = TRUE)

      # Filter parameter file: update
      filters.parameters <- radiator_parameters(
        generate = FALSE,
        initiate = FALSE,
        update = TRUE,
        parameter.obj = filters.parameters,
        data = gds,
        filter.name = "duplicated markers on different strands",
        param.name = "filter.strands",
        values = filter.strands,
        path.folder = path.folder,
        file.date = file.date,
        verbose = verbose)
    }
  }
  blacklist.strands <- NULL



  #______________________________###############################################
  #-----------------------------------  VCF STATS   ----------------------------


  if (vcf.stats) {
    # SUBSAMPLE markers --------------------------------------------------------
    markers.meta <- extract_markers_metadata(gds = gds)#update
    n.markers <- length(markers.meta$VARIANT_ID)
    if (n.markers < 200000) subsample.markers.stats <- 1

    if (subsample.markers.stats < 1) {
      markers.subsampled <- dplyr::sample_frac(tbl = markers.meta,
                                               size = subsample.markers.stats)
      variant.select <- markers.subsampled$VARIANT_ID
      subsample.filename <- stringi::stri_join("markers.subsampled_", file.date, ".tsv")
      dplyr::select(markers.subsampled, MARKERS) %>%
        dplyr::mutate(RANDOM_SEED = random.seed) %>%
        readr::write_tsv(x = .,
                         path = file.path(path.folder, subsample.filename))
      markers.subsampled <- NULL
    } else {
      variant.select <- NULL
    }

    # Individuals stats --------------------------------------------------------
    if (verbose) message("\nGenerating individual stats")
    # Note to my self: the test you ran show that timing of missing prop is
    # not really impacted by number of markers, but heterozygosity is...
    id.stats <- generate_id_stats(
      gds = gds,
      subsample = variant.select,
      depth.info = dp,
      path.folder = path.folder,
      file.date = file.date,
      parallel.core = parallel.core,
      verbose = verbose)

    # names(id.stats)
    print(id.stats$fig)
    individuals <- id.stats$info %>%
      readr::write_tsv(x = ., path = file.path(path.folder, ind.file))
    ind.stats.file <- stringi::stri_join("vcf_individuals_stats_", file.date, ".tsv")
    ind.stats <- id.stats$stats %>%
      readr::write_tsv(x = ., path = file.path(path.folder, ind.stats.file))

    # filter.individuals.missing -----------------------------------------------
    gds <- filter_individuals(
      interactive.filter = FALSE,
      data = gds,
      filter.individuals.missing = filter.individuals.missing,
      filter.individuals.heterozygosity = NULL,
      filter.individuals.coverage.total = filter.individuals.coverage.total,
      parallel.core = parallel.core,
      verbose = verbose,
      path.folder = path.folder,
      subsample = variant.select,
      parameters = filters.parameters,
      dp = TRUE,
      id.stats = id.stats
    )


    # Markers stats ------------------------------------------------------------
    if (verbose) message("Updating markers metadata and stats\n")
    temp <- generate_markers_stats(
      gds = gds,
      path.folder = path.folder,
      filename = markers.file,
      file.date = file.date,
      parallel.core = parallel.core,
      subsample = variant.select)
    stats.sub <- temp$info

    # Note to my self: might no longer be necessary if always using the meta inside the GDS
    if (is.null(variant.select)) {
      markers.meta <- temp$info
      stats.sub <- NULL
    } else {
      stats.sub <- temp$info
    }
    m.stats <- temp$stats# required below
    temp <- NULL

    ##______________________________________________________________________####
    # Filters ------------------------------------------------------------------
    # Filter common markers ----------------------------------------------------
    gds <- filter_common_markers(
      data = gds,
      filter.common.markers = filter.common.markers,
      plot = TRUE,
      parallel.core = parallel.core,
      verbose = verbose,
      parameters = filters.parameters,
      path.folder = path.folder,
      internal = internal)

    # FILTER MAC----------------------------------------------------------------
    gds <- filter_mac(
      interactive.filter = FALSE,
      data = gds,
      filter.mac = filter.mac,
      filename = NULL,
      parallel.core = parallel.core,
      verbose = verbose,
      parameters = filters.parameters,
      path.folder = path.folder)

    # Filter mean coverage------------------------------------------------------
    gds <- filter_coverage(
      interactive.filter = FALSE,
      data = gds,
      filter.coverage = filter.coverage,
      filename = NULL,
      parallel.core = parallel.core,
      verbose = verbose,
      parameters = filters.parameters,
      path.folder = path.folder)

    # Filter genotyping---------------------------------------------------------
    gds <- filter_genotyping(
      interactive.filter = FALSE,
      data = gds,
      filter.genotyping = filter.genotyping,
      filename = NULL,
      parallel.core = parallel.core,
      verbose = verbose,
      parameters = filters.parameters,
      path.folder = path.folder)

    # SNP position on the read--------------------------------------------------
    gds <- filter_snp_position_read(
      interactive.filter = FALSE,
      data = gds,
      filter.snp.position.read = filter.snp.position.read,
      filename = NULL,
      parallel.core = parallel.core,
      verbose = verbose,
      parameters = filters.parameters,
      path.folder = path.folder)

    # number of SNPs per locus -------------------------------------------------
    gds <- filter_snp_number(
      interactive.filter = FALSE,
      data = gds,
      filter.snp.number = filter.snp.number,
      filename = NULL,
      parallel.core = parallel.core,
      verbose = verbose,
      parameters = filters.parameters,
      path.folder = path.folder)

    # Filter Linkage disequilibrium --------------------------------------------
    gds <- filter_ld(
      interactive.filter = FALSE,
      data = gds,
      filter.short.ld = filter.short.ld,
      filter.long.ld = filter.long.ld,
      parallel.core = parallel.core,
      filename = NULL,
      verbose = verbose,
      long.ld.missing = long.ld.missing,
      ld.method = ld.method,
      parameters = filters.parameters,
      path.folder = path.folder)

    # Final Sync GDS -----------------------------------------------------------
    markers.meta <- extract_markers_metadata(gds)
    sync_gds(gds = gds, samples = individuals$INDIVIDUALS,
             markers = markers.meta$VARIANT_ID)


    # For summary at the end of the function:
    ind.missing <- round(mean(individuals$MISSING_PROP, na.rm = TRUE), 2)
    if (dp) ind.cov.total <- round(mean(individuals$COVERAGE_TOTAL, na.rm = TRUE), 0)
    if (dp) ind.cov.mean <- round(mean(individuals$COVERAGE_MEAN, na.rm = TRUE), 0)

    if (!is.null(stats.sub)) {
      markers.missing <- round(mean(stats.sub$MISSING_PROP, na.rm = TRUE), 2)
      if (dp) markers.cov <- round(mean(stats.sub$COVERAGE_MEAN, na.rm = TRUE), 0) # same as above because NA...
    } else {
      markers.missing <- round(mean(markers.meta$MISSING_PROP, na.rm = TRUE), 2)
      if (dp) markers.cov <- round(mean(markers.meta$COVERAGE_MEAN, na.rm = TRUE), 0) # same as above because NA...
    }
  }#End stats

  # RESULTS --------------------------------------------------------------------
  number.info <- SeqArray::seqGetFilter(gdsfile = gds)
  if (verbose) {
    cat("################################### SUMMARY ####################################\n")
    if (vcf.stats) {
      message("\n\nSummary (before filtering):\nMissing data: ")
      message("    markers: ", markers.missing)
      message("    individuals: ", ind.missing)
      if (dp) message("\n\nCoverage info:")
      if (dp) message("    individuals mean read depth: ", ind.cov.total)
      if (dp) message("    individuals mean genotype coverage: ", ind.cov.mean)
      if (dp) message("    markers mean coverage: ", markers.cov)
    }
  }
  message("\n\nNumber of chromosome/contig/scaffold: ", dplyr::n_distinct(markers.meta$CHROM))
  message("Number of locus: ", dplyr::n_distinct(markers.meta$LOCUS))
  message("Number of markers: ", length(number.info$variant.sel[number.info$variant.sel]))
  message("Number of individuals: ", length(number.info$sample.sel[number.info$sample.sel]))

  individuals <- stats <- markers.meta <- NULL

  message("radiator Genomic Data Structure (GDS) file: ", folder_short(filename))
  return(gds)
} # End write_SeqArray

