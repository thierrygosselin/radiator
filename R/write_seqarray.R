# write a SeqArray object from a tidy data frame

#' @name write_seqarray
#' @title Write a SeqArray GDS file from a vcf file and generate a connection object.
#' @description Write a SeqArray \href{http://zhengxwen.github.io/SeqArray/}{SeqArray}
#' file (Zheng et al. 2017) from a vcf file and generate a connection object.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.


#' @param data (character string) The VCF SNPs are biallelic or haplotypes.
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
#' radiator will append \code{.gds} to the filename.
#' If filename chosen is already present in the
#' working directory, the default \code{radiator_datetime.gds} is chosen.
#' Default: \code{filename = NULL}.

#' @param vcf.stats (optional, logical) Generates individuals and markers
#' important statistics helpful for filtering.
#' These are very fast to generate and because computational
#' cost is minimal, even for huge VCFs, the default is \code{vcf.stats = TRUE}.
#' Starts: individual's missing genotype proportion, averaged heterozygosity,
#' total coverage, mean genotype coverage and marker's metadata
#' along count for ref and alt alleles and mean
#' coverage is generated and written in the working directory.
#' Default: \code{vcf.stats = TRUE}.


#' @param ... (optional) Advance mode that allows to pass further arguments
#' for fine-tuning the function (see details).

#' @inheritParams tidy_genomic_data


#' @details
#' A vcf file of 35 GB with ~4 millions SNPs take about ~7 min with 8 CPU.
#' A vcf file of 21 GB with ~2 millions SNPs take about ~5 min with 7 CPU.
#'
#' After the file is generated, you can close your computer and
#' come back to it a month later and it's now a matter of sec to open a connection.
#'
#' \strong{Advance mode, using \emph{dots-dots-dots ...}}
#' \enumerate{
#' \item \code{whitelist.markers}: detailed in \code{\link[radiator]{tidy_genomic_data}}.
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
#' \item \code{common.markers}: (logical) Default: code{common.markers = TRUE}.
#' Detailed in \code{\link[radiator]{tidy_genomic_data}}.
#' \item \code{filter.mac}: (integer) To blacklist markers below a specific MAC
#' (calculated overall).
#' \item \code{filter.coverage.outliers}: (logical) To blacklist markers with low and
#' high coverage based on outlier statistics.
#' \item \code{filter.markers.missing}: (integer) To blacklist markers with too
#' many missing data. e.g. \code{filter.markers.missing = 10}, will only keep
#' markers with missing rate <= to 10 percent.
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

#' @seealso \code{\link{snp_ld}} and \code{\link[radiator]{tidy_genomic_data}}

#' @examples
#' \dontrun{
#' # with built-in defaults:
#'  prep.data <- radiator::write_seqarray(data = "populations.snps.vcf")
#'
#' # Using more arguments and filters (recommended):
#' prep.data <- radiator::write_seqarray(
#'     data = "populations.snps.vcf",
#'     strata = "strata_salamander.tsv",
#'     vcf.stats = TRUE,
#'     filter.individuals.missing = "outlier",
#'     common.markers = TRUE,
#'     filter.strands = "blacklist",
#'     filter.mac = 4,
#'     filter.markers.missing = 50,
#'     filter.snp.read.position = "outliers",
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

  ##Test
  # data = "populations.snps.vcf"
  # strata <- "spis-popmap-448samples.tsv"
  # filename <- NULL
  # vcf.stats <- TRUE
  # parallel.core <- parallel::detectCores() - 1
  # verbose <- TRUE
  # vcf.metadata = TRUE
  # filter.strands = "blacklist"
  # common.markers = TRUE
  # filter.mac <- 4
  # filter.coverage.outliers = TRUE
  # filter.markers.missing <- 10
  # filter.snp.read.position <- c("outliers", "q75", "iqr")
  # filter.short.ld <- "maf"
  # filter.long.ld <- 0.8
  # long.ld.missing <- TRUE
  # ld.method <- "r2"
  # filter.individuals.missing <- 0.7
  # blacklist.id = NULL
  # pop.select = NULL
  # pop.levels = NULL
  # pop.labels = NULL
  # whitelist.markers = NULL
  # keep.gds <- TRUE
  # markers.info = NULL
  # path.folder = NULL
  # subsample.markers.stats = 0.2

  res <- list() #store the results
  timing.import <- proc.time()
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

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("vcf file missing")

  # dotslist -------------------------------------------------------------------
  dotslist <- list(...)
  want <- c("whitelist.markers",
            "filter.snp.read.position", "filter.mac",
            "filter.coverage.outliers", "filter.markers.missing", "filter.short.ld",
            "filter.long.ld", "long.ld.missing", "ld.method",
            "filter.individuals.missing", "common.markers",
            "filter.strands",
            "blacklist.id", "pop.select", "pop.levels", "pop.labels", "keep.gds",
            "path.folder", "markers.info", "vcf.metadata",
            "subsample.markers.stats")
  unknowned_param <- setdiff(names(dotslist), want)

  if (length(unknowned_param) > 0) {
    stop("Unknowned \"...\" parameters ",
         stringi::stri_join(unknowned_param, collapse = " "))
  }

  radiator.dots <- dotslist[names(dotslist) %in% want]
  whitelist.markers <- radiator.dots[["whitelist.markers"]]

  filter.snp.read.position <- radiator.dots[["filter.snp.read.position"]]
  filter.mac <- radiator.dots[["filter.mac"]]
  filter.coverage.outliers <- radiator.dots[["filter.coverage.outliers"]]
  filter.markers.missing <- radiator.dots[["filter.markers.missing"]]
  filter.short.ld <- radiator.dots[["filter.short.ld"]]
  filter.long.ld <- radiator.dots[["filter.long.ld"]]
  long.ld.missing <- radiator.dots[["long.ld.missing"]]
  ld.method <- radiator.dots[["ld.method"]]
  filter.individuals.missing <- radiator.dots[["filter.individuals.missing"]]
  common.markers <- radiator.dots[["common.markers"]]
  filter.strands <- radiator.dots[["filter.strands"]]
  markers.info <- radiator.dots[["markers.info"]]
  vcf.metadata <- radiator.dots[["vcf.metadata"]]
  subsample.markers.stats <- radiator.dots[["subsample.markers.stats"]]

  blacklist.id <- radiator.dots[["blacklist.id"]]
  pop.select <- radiator.dots[["pop.select"]]
  pop.levels <- radiator.dots[["pop.levels"]]
  pop.labels <- radiator.dots[["pop.labels"]]
  keep.gds <- radiator.dots[["keep.gds"]]
  path.folder <- radiator.dots[["path.folder"]]

  # useful outside this function
  if (is.null(keep.gds)) keep.gds <- TRUE
  if (is.null(filter.strands)) filter.strands <- "blacklist"
  if (is.null(long.ld.missing)) long.ld.missing <- FALSE

  if (is.null(ld.method)) {
    ld.method <- "r2"
  } else {
    ld.method <- match.arg(ld.method, c("composite", "r", "r2", "dprime", "corr"))
  }

  if (is.null(filter.coverage.outliers)) filter.coverage.outliers <- FALSE
  if (is.null(common.markers)) common.markers <- TRUE
  if (is.null(path.folder)) {
    path.folder <- getwd()
  } else {
    if (!dir.exists(path.folder)) dir.create(path.folder)
  }

  if (!is.null(filter.snp.read.position)) {
    filter.snp.read.position <- match.arg(
      arg = filter.snp.read.position,
      choices = c("outliers", "iqr", "q75"),
      several.ok = TRUE)
  }
  if (is.null(subsample.markers.stats)) subsample.markers.stats <- 0.2

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


  # LD
  # currently: will only filter long ld if short ld as been taken care of first...
  if (!is.null(filter.long.ld) && is.null(filter.short.ld)) {
    message("\nfilter.short.ld argument set by default to: maf")
    filter.short.ld <- "maf"
  }

  # Start --------
  message("\nReading VCF...")

  # Get file size
  big.vcf <- file.size(data)

  if (verbose) {
    if (big.vcf > 500000000) message("Large vcf file may take several minutes...")
    if (big.vcf > 5000000000) message("    you actually have time for a coffee!\n")
  }

  # Filename -------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
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

  # Filter parameter file: initiate --------------------------------------------
  filters.parameters <- update_parameters(initiate = TRUE,
                                          path.folder = path.folder,
                                          file.date = file.date,
                                          verbose = verbose)

  # Read vcf -------------------------------------------------------------------
  timing.vcf <- proc.time()

  # Check for bad header generated by stacks
  detect.source <- check_header_source(data)
  stacks.2 <- detect.source$source
  check.header <- detect.source$check.header
  dp <- "DP" %in% detect.source$check.header$format$ID # Check that DP is trere

  if (!is.null(detect.source$markers.info)) markers.info <- detect.source$markers.info
  if (!is.null(detect.source$overwrite.metadata)) overwrite.metadata <- detect.source$overwrite.metadata

  res$vcf.connection <- SeqArray::seqVCF2GDS(
    vcf.fn = data,
    out.fn = filename,
    parallel = parallel.core,
    storage.option = "ZIP_RA",
    verbose = FALSE,
    header = check.header,
    info.import = markers.info, # characters, the variable name(s) in the INFO field for import; or NULL for all variables
    fmt.import = overwrite.metadata# characters, the variable name(s) in the FORMAT field for import; or NULL for all variables#
  ) %>%
    SeqArray::seqOpen(gds.fn = ., readonly = FALSE)

  check.header <- detect.source <- NULL

  # Summary --------------------------------------------------------------------
  # vcf.sum <- SeqArray::seqSummary(res$vcf.connection, verbose = FALSE)
  # this might be safer, because using seqSummary seems to forget if filters were used
  n.ind <- length(SeqArray::seqGetData(res$vcf.connection, "sample.id"))
  n.markers <- length(SeqArray::seqGetData(res$vcf.connection, "variant.id"))
  message("\nNumber of SNPs: ", n.markers)
  message("Number of samples: ", n.ind)
  # vcf.sum <- NULL

  if (verbose) message("\nconversion timing: ", round((proc.time() - timing.vcf)[[3]]), " sec")

  if (verbose && keep.gds) {
    message("\nGDS file generated: ", filename.short)
  }
  # radiator gds folder --------------------------------------------------------
  radiator.gds <- gdsfmt::addfolder.gdsn(
    node = res$vcf.connection,
    name = "radiator",
    replace = TRUE)

  # bi- or multi-alllelic VCF --------------------------------------------------
  # Haplotypes or SNPs
  ref <- SeqArray::seqGetData(gdsfile = res$vcf.connection, var.name = "$ref")
  sample.size <- min(length(unique(ref)), 100)
  res$biallelic <- max(unique(stringi::stri_count_regex(
    str = sample(
      x = unique(ref),
      size = sample.size), pattern = "[A-Z]"))) == 1
  sample.size <- ref <- NULL

  if (res$biallelic) {
    message("\nVCF: biallelic SNPs")
  } else {
    message("\nVCF: haplotypes")
  }

  # add biallelic info to GDS
  gdsfmt::add.gdsn(
    node = radiator.gds,
    name = "biallelic",
    val = res$biallelic,
    replace = TRUE,
    compress = "ZIP_RA",
    closezip = TRUE)

  # Strata ---------------------------------------------------------------------
  # import strata and filter with blacklist of id if present...
  strata <- radiator::read_strata(
    strata = strata,
    pop.levels = pop.levels,
    pop.labels = pop.labels,
    pop.select = pop.select,
    blacklist.id = blacklist.id) %$% strata

  # clean sample id in VCF -----------------------------------------------------
  individuals.vcf <- tibble::tibble(
    INDIVIDUALS_VCF = SeqArray::seqGetData(res$vcf.connection, "sample.id")) %>%
    dplyr::mutate(INDIVIDUALS_CLEAN = radiator::clean_ind_names(INDIVIDUALS_VCF))

  if (!identical(individuals.vcf$INDIVIDUALS_VCF, individuals.vcf$INDIVIDUALS_CLEAN)) {
    if (verbose) message("Cleaning VCF sample names")
    clean.id.filename <- stringi::stri_join("cleaned.vcf.id.info_", file.date, ".tsv")
    readr::write_tsv(x = individuals.vcf,
                     path = stringi::stri_join(path.folder, "/", clean.id.filename))
    gdsfmt::add.gdsn(
      node = radiator.gds,
      name = "vcf.id.clean",
      val = individuals.vcf,
      replace = TRUE,
      compress = "ZIP_RA",
      closezip = TRUE)
  }

  # replace id in VCF
  gdsfmt::add.gdsn(
    node = res$vcf.connection,
    name = "sample.id",
    val = individuals.vcf$INDIVIDUALS_CLEAN,
    replace = TRUE,
    compress = "ZIP_RA",
    closezip = TRUE)


  res$individuals <- dplyr::select(individuals.vcf, INDIVIDUALS = INDIVIDUALS_CLEAN)

  # Add a individuals node
  gdsfmt::add.gdsn(
    node = radiator.gds,
    name = "individuals",
    val = res$individuals,
    replace = TRUE,
    compress = "ZIP_RA",
    closezip = TRUE)

  # Sync id in STRATA and VCF --------------------------------------------------
  if (!is.null(strata)) {
    if (verbose) message("\nSynchronizing sample IDs in VCF and strata...")
    strata %<>% dplyr::filter(INDIVIDUALS %in% individuals.vcf$INDIVIDUALS_CLEAN)
    individuals.vcf <- NULL

    blacklist.strata <- n.ind - nrow(strata)
    if (blacklist.strata != 0) {
      if (verbose) message("    number of sample blacklisted by the stata: ", blacklist.strata)
      # the param file is updated after markers metadata below
    }

    # SeqArray::seqSetFilter(object = res$vcf.connection,
    #                        sample.id = strata$INDIVIDUALS,
    #                        verbose = FALSE)

    sync_gds(gds = res$vcf.connection, samples = strata$INDIVIDUALS)

    res$individuals <- strata
    gdsfmt::add.gdsn(
      node = radiator.gds,
      name = "individuals",
      val = res$individuals,
      replace = TRUE,
      compress = "ZIP_RA",
      closezip = TRUE)

    strata <- TRUE
    # # Add STRATA to GDS
    # gdsfmt::add.gdsn(
    #   node = radiator.gds,
    #   name = "STRATA",
    #   val = strata,
    #   # val = strata$STRATA,
    #   replace = TRUE,
    #   compress = "ZIP_RA",
    #   closezip = TRUE)
  } else {
    strata <- FALSE
    blacklist.strata <- 0L

    res$individuals %<>% dplyr::mutate(STRATA = 1L)
    gdsfmt::add.gdsn(
      node = radiator.gds,
      name = "individuals",
      val = res$individuals,
      replace = TRUE,
      compress = "ZIP_RA",
      closezip = TRUE)
  }


  # Markers metadata  ----------------------------------------------------------
  res$markers.meta <- tibble::tibble(
    VARIANT_ID = SeqArray::seqGetData(res$vcf.connection, "variant.id"),
    CHROM = SeqArray::seqGetData(res$vcf.connection, "chromosome"),
    LOCUS = SeqArray::seqGetData(res$vcf.connection, "annotation/id"),
    POS = SeqArray::seqGetData(res$vcf.connection, "position"))

  # check <- res$markers.meta

  # reference genome or de novo ------------------------------------------------
  ref.genome <- sample(x = unique(res$markers.meta$CHROM),
                       size = min(length(unique(res$markers.meta$CHROM)), 100),
                       replace = FALSE)

  # if the chrom.unique > 1 more likely not to be de novo assembly (e.g. with old stacks version)
  chrom.unique <- length(unique(ref.genome)) == 1

  # presence of underscore or other separator: more likely ref genome
  chrom.sep <- TRUE %in%
    stringi::stri_detect_regex(str = ref.genome, pattern = "[^[:alnum:]]+") %>%
    unique

  # presence of letters = more likely ref genome
  chrom.alpha <- TRUE %in%
    stringi::stri_detect_regex(str = ref.genome, pattern = "[[:alpha:]]+") %>%
    unique

  if (chrom.unique) ref.genome <- FALSE
  if (chrom.alpha || chrom.sep) ref.genome <- TRUE
  if (ref.genome) {
    if (verbose) message("Reads assembly: reference-assisted")
  } else {
    if (verbose) message("Reads assembly: de novo")
  }

  chrom.unique <- chrom.alpha <- chrom.sep <- NULL
  # Add to GDS
  gdsfmt::add.gdsn(
    node = radiator.gds,
    name = "reference.genome",
    val = ref.genome,
    replace = FALSE,
    compress = "ZIP_RA",
    closezip = TRUE)

  # Stacks specific adjustments
  if (!ref.genome) {
    if (stacks.2) {
      res$markers.meta %<>% dplyr::mutate(
        LOCUS = CHROM,
        CHROM = "1",
        COL = POS - 1)
    } else {
      res$markers.meta %<>% dplyr::mutate(
        CHROM = "1")
    }
  }

  # GATK, platypus and freebayes specific adjustment
  # Locus with NA or . or ""
  weird.locus <- length(unique(res$markers.meta$LOCUS)) <= 1
  if (weird.locus && !stacks.2) {
    if (verbose) message("LOCUS field empty... adding unique id instead")
    res$markers.meta$LOCUS <- res$markers.meta$VARIANT_ID
  }

  # LOCUS cleaning and Strands detection -------------------------------------
  if (stringi::stri_detect_regex(str = res$markers.meta[1,3], pattern = "[^[:alnum:]]+")) {
    locus.sep <- unique(stringi::stri_extract_all_regex(
      str = res$markers.meta[1,3],
      pattern = "[^a-zA-Z0-9-+-]",
      omit_no_match = TRUE)[[1]])

    res$markers.meta <- suppressWarnings(
      tidyr::separate(
        data = res$markers.meta,
        col = LOCUS, into = c("LOCUS", "COL", "STRANDS"),
        sep = locus.sep,
        # sep = "",
        extra = "drop", fill = "warn",
        remove = TRUE, convert = TRUE)
    )
    if (anyNA(res$markers.meta$STRANDS)) res$markers.meta$STRANDS <- NA_character_
    locus.sep <- NULL

    detect.strand <- any(stringi::stri_detect_fixed(str = unique(res$markers.meta$STRANDS), pattern = "+"))
    if (anyNA(detect.strand)) detect.strand <- FALSE

    if (detect.strand) {
      blacklist.strands <- dplyr::distinct(res$markers.meta, CHROM, LOCUS, POS) %>%
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
          message("\n\nDetected ", nrow(blacklist.strands)," duplicate SNPs on different strands (+/-)")
          message("    By default radiator prune those SNPs")
          message("    To change this behavior, use argument: filter.strands\n\n")
        }
      }

    } else {
      res$markers.meta %<>% dplyr::select(-STRANDS)
      blacklist.strands <- NULL
    }
    detect.strand <- NULL
  } else {
    blacklist.strands <- NULL
  }

  # Generate MARKERS column and fix types
  res$markers.meta %<>%
    dplyr::mutate_at(
      .tbl = .,
      .vars = c("CHROM", "LOCUS", "POS"),
      .funs = radiator::clean_markers_names) %>%
    dplyr::mutate(
      MARKERS = stringi::stri_join(CHROM, LOCUS, POS, sep = "__"),
      REF = SeqArray::seqGetData(gdsfile = res$vcf.connection, var.name = "$ref"),
      ALT = SeqArray::seqGetData(gdsfile = res$vcf.connection, var.name = "$alt")
    )

  # # ADD MARKERS META to GDS --------------------------------------------------
  suppressWarnings(gdsfmt::add.gdsn(
    node = radiator.gds,
    name = "markers.meta",
    val = res$markers.meta,
    replace = TRUE,
    compress = "ZIP_RA",
    closezip = TRUE))

  # replace chromosome info in GDS
  # Why ? well snp ld e.g. will otherwise be performed by chromosome and with de novo data = by locus...
  gdsfmt::add.gdsn(
    node = res$vcf.connection,
    name = "chromosome",
    val = res$markers.meta$CHROM,
    replace = TRUE,
    compress = "ZIP_RA",
    closezip = TRUE)

  # update filter param ----------------------------------------------------------
  # The original VCF's values
  filters.parameters <- update_parameters(
    parameter.obj = filters.parameters,
    gds.obj = res,
    filter.name = "vcf",
    param.name = "original values in vcf",
    values = "")

  if (blacklist.strata != 0 && strata) {
    filters.parameters <- update_parameters(
      parameter.obj = filters.parameters,
      gds.obj = res,
      filter.name = "vcf",
      param.name = "original values in vcf + strata",
      values = "")
  }

  # PRE-FILTERING --------------------------------------------------------------
  # whitelist of markers, the filter column in the vcf and blacklisted strands
  # The order here doesnt matter

  # Generate a blacklist of markers to store info ----------------------------
  res$blacklist.markers <- tibble::tibble(MARKERS = character(0), FILTER = character(0))

  # Filter with whitelist of markers--------------------------------------------
  if (!is.null(whitelist.markers)) {

    # Import whitelist of markers
    if (is.vector(whitelist.markers)) {
      whitelist.markers <- suppressMessages(
        readr::read_tsv(whitelist.markers, col_names = TRUE) %>%
          dplyr::mutate_all(.tbl = ., .funs = as.character))
    }
    columns.names.whitelist <- colnames(whitelist.markers)
    nrow.before <- nrow(whitelist.markers)
    whitelist.markers <- dplyr::distinct(whitelist.markers)
    nrow.after <- nrow(whitelist.markers)
    duplicate.whitelist.markers <- nrow.before - nrow.after
    if (duplicate.whitelist.markers > 0) {
      message("Whitelist of markers with ", duplicate.whitelist.markers, " duplicated identifiers...")
      message("    Creating unique whitelist")
      message("    Warning: downstream results might be impacted by this, check how you made your VCF file...")
    }
    nrow.before <- duplicate.whitelist.markers <- NULL

    whitelist.markers <- dplyr::mutate_all(
      .tbl = whitelist.markers, .funs = clean_markers_names)

    if (!biallelic) {
      if (ncol(whitelist.markers) >= 3) {
        message("Note: whitelist with CHROM LOCUS POS columns and VCF haplotype:
                If the whitelist was not created from this VCF,
                the filtering could result in losing all the markers.
                The POS column is different in biallelic and multiallelic file...\n")

        message("Discarding the POS column in the whitelist")
        whitelist.markers <- dplyr::select(whitelist.markers, -POS)
      }

      if (ncol(whitelist.markers) == 1 && tibble::has_name(whitelist.markers, "MARKERS")) {
        message("Note: whitelist MARKERS column and VCF haplotype:
                If the whitelist was not created from this VCF,
                the filtering could result in losing all the markers.
                The POS column used in the MARKERS column is different in biallelic and multiallelic file...\n")
      }
    }
    nrow.after <- nrow(whitelist.markers)
    message("Filtering: ", nrow.after, " markers in whitelist")
    columns.names.whitelist <- colnames(whitelist.markers)

    variant.select <- suppressWarnings(
      dplyr::semi_join(
        res$markers.meta, whitelist.markers, by = columns.names.whitelist)) %>%
      dplyr::select(VARIANT_ID) %>%
      purrr::flatten_int(.)

    SeqArray::seqSetFilter(object = res$vcf.connection,
                           variant.id = variant.select,
                           verbose = FALSE)


    if (length(SeqArray::seqGetData(res$vcf.connection, "variant.id")) == 0) {
      stop("No markers left in the dataset, check whitelist...")
    }

    # update filter param ----------------------------------------------------------
    filters.parameters <- update_parameters(
      parameter.obj = filters.parameters,
      gds.obj = res,
      filter.name = "whitelist markers",
      param.name = "whitelist.markers",
      values = nrow.after)
  }# End whitelist markers

  # Scan and filter with FILTER column ---------------------------------------
  res$markers.meta$FILTER <- SeqArray::seqGetData(
    res$vcf.connection, "annotation/filter")
  filter.check.unique <- unique(res$markers.meta$FILTER)

  if (length(filter.check.unique) > 1) {
    message("Filtering markers based on VCF FILTER column")
    n.markers.before <- nrow(res$markers.meta)

    blacklist.vcf.filter <- res$markers.meta %>%
      dplyr::filter(FILTER != "PASS") %>%
      dplyr::select(MARKERS)

    # update the blacklist
    res$blacklist.markers %<>% dplyr::bind_rows(
      dplyr::mutate(blacklist.vcf.filter, FILTER = "vcf filter column")
    )

    res$markers.meta %<>% dplyr::filter(FILTER == "PASS")
    n.markers.after <- nrow(res$markers.meta)
    n.markers <- stringi::stri_join(n.markers.before, n.markers.before - n.markers.after, n.markers.after, sep = " / ")
    if (verbose) message("    Number of SNPs before / blacklisted / after: ", n.markers)
    # SeqArray::seqSetFilter(object = res$vcf.connection,
    #                        variant.id = res$markers.meta$VARIANT_ID,
    #                        verbose = FALSE)
    sync_gds(gds = res$vcf.connection, markers = res$markers.meta$VARIANT_ID)


    # update filter param ----------------------------------------------------------
    filters.parameters <- update_parameters(
      parameter.obj = filters.parameters,
      gds.obj = res,
      filter.name = "vcf filter column",
      param.name = "PASS or not",
      values = "")
  }
  filter.check.unique <- NULL
  res$markers.meta %<>% dplyr::select(-FILTER)
  # check <- res$markers.meta

  # FILTER STRANDS -----------------------------------------------------------
  if (!is.null(blacklist.strands)) {
    blacklist.strands <- dplyr::right_join(res$markers.meta, blacklist.strands,
                                           by = c("CHROM", "POS")) %>%
      dplyr::select(VARIANT_ID, MARKERS, CHROM, POS)

    if (filter.strands == "keep.both") {
      message("Number of duplicated markers on different strands removed: 0")
    }

    if (filter.strands == "best.stats") {
      # SeqArray::seqSetFilter(object = res$vcf.connection,
      #                        variant.id = blacklist.strands$VARIANT_ID,
      #                        verbose = FALSE)
      sync_gds(gds = res$vcf.connection, markers = blacklist.strands$VARIANT_ID)


      blacklist.strands <- SeqArray::seqAlleleCount(
        gdsfile = res$vcf.connection,
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
                             gdsfile = res$vcf.connection,
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
      message("Number of duplicated markers on different strands blacklisted: ", nrow(blacklist.strands))

      # update the blacklist
      res$blacklist.markers %<>% dplyr::bind_rows(
        dplyr::mutate(blacklist.strands, FILTER = "filter.strands")
      )

      # Update markers.meta
      res$markers.meta %<>% dplyr::filter(!MARKERS %in% blacklist.strands$MARKERS)
      # check <- res$markers.meta

      # Update GDS
      gdsfmt::add.gdsn(
        node = radiator.gds,
        name = "markers.meta",
        val = res$markers.meta,
        replace = TRUE,
        compress = "ZIP_RA",
        closezip = TRUE)

      # SeqArray::seqSetFilter(object = res$vcf.connection,
      #                        variant.id = res$markers.meta$VARIANT_ID,
      #                        verbose = FALSE)
      sync_gds(gds = res$vcf.connection, markers = res$markers.meta$VARIANT_ID)


      # update filter param ----------------------------------------------------------
      filters.parameters <- update_parameters(
        parameter.obj = filters.parameters,
        gds.obj = res,
        filter.name = "duplicated markers on different strands",
        param.name = "filter.strands",
        values = filter.strands)
    }
  }

  blacklist.strands <- NULL

  #-----------------------------------  VCF STATS   ----------------------------
  if (vcf.stats) {
    # SUBSAMPLE markers when the number if high --------------------------------
    # 10-20% seems to give identical results to 1M

    whitelist.variant.id <- res$markers.meta$VARIANT_ID
    n.markers <- length(whitelist.variant.id)
    if (n.markers > 200000) {
      m.subsample <- TRUE
      variant.select <- sample(x = res$markers.meta$VARIANT_ID,
                               size = round(subsample.markers.stats * length(res$markers.meta$VARIANT_ID), 0))
      res$markers.meta %<>% dplyr::mutate(SUBSAMPLE_STATS = dplyr::if_else(VARIANT_ID %in% variant.select, TRUE, FALSE))
      # check <- res$markers.meta
    } else {
      m.subsample <- FALSE
      res$markers.meta %<>% dplyr::mutate(SUBSAMPLE_STATS = TRUE)
    }

    # variant.select.vcf.info <- sample(x = res$markers.meta$VARIANT_ID,
    #                                   size = min(10L, length(res$markers.meta$VARIANT_ID)))

    # Individuals stats --------------------------------------------------------
    if (verbose) message("Generating individual stats")
    # Note to my self: the test you ran show that timing of missing prop is
    # not really impacted by number of markers, but heterozygosity is...

    res$individuals %<>% dplyr::mutate(
      MISSING_PROP = round(SeqArray::seqMissing(
        gdsfile = res$vcf.connection, per.variant = FALSE,
        .progress = TRUE,
        parallel = parallel.core), 6))

    if (m.subsample) {
      SeqArray::seqSetFilter(object = res$vcf.connection,
                             variant.id = variant.select,
                             sample.id = res$individuals$INDIVIDUALS,
                             action = "push+set",
                             verbose = FALSE)

      # check <- SeqArray::seqGetFilter(res$vcf.connection)
      # length(check$sample.sel[check$sample.sel])
      # length(check$variant.sel[check$variant.sel])
      # n.markers
    }
    res$individuals %<>% dplyr::mutate(
      HETEROZYGOSITY = round(SeqVarTools::heterozygosity(
        gdsobj = res$vcf.connection, margin = "by.sample", use.names = FALSE), 6)
    )

    if (m.subsample) {
      SeqArray::seqSetFilter(res$vcf.connection, action = "pop", verbose = TRUE)
      # check <- SeqArray::seqGetFilter(res$vcf.connection)
      # length(check$sample.sel[check$sample.sel])
      # length(check$variant.sel[check$variant.sel])
      # n.markers
    }
    # check <- res$individuals

    if (!dp && verbose && vcf.stats) {
      message("DP genotype metadata in VCF: no")
      message("    skipping coverage summary statistics (individuals and markers")
    }


    if (stacks.2 && !res$biallelic) {
      if (verbose) message("\nHaplotype VCF generated by stacks doesn't have coverage info")
    } else {
      if (dp) {
        # if (m.subsample) {
        #   SeqArray::seqSetFilter(object = res$vcf.connection,
        #                          variant.id = variant.select,
        #                          action = "push+set",
        #                          verbose = FALSE)
        # }
        id.stats.before.filter <- extract_coverage(data = res$vcf.connection, markers = FALSE)
        # if (m.subsample) {
        #   SeqArray::seqSetFilter(res$vcf.connection, action = "pop", verbose = TRUE)
        # }
        res$individuals %<>% dplyr::mutate(
          COVERAGE_TOTAL = id.stats.before.filter$ind.cov.tot,
          COVERAGE_MEAN = id.stats.before.filter$ind.cov.mean
        )
      }
    }
    # check <- res$individuals

    readr::write_tsv(x = res$individuals, path = file.path(path.folder, ind.file))
    id.stats.before.filter <- NULL

    # test <- res$individuals

    if (stacks.2 && !res$biallelic) {
      res$stats$ind.stats <- tibble_stats(
        x = res$individuals$MISSING_PROP,
        group = "individual's missing genotypes") %>%
        dplyr::bind_rows(
          tibble_stats(
            x = res$individuals$HETEROZYGOSITY,
            group = "individual's heterozygosity")
        )
    } else {
      res$stats$ind.stats <- tibble_stats(
        x = res$individuals$MISSING_PROP,
        group = "individual's missing genotypes") %>%
        dplyr::bind_rows(
          tibble_stats(
            x = res$individuals$HETEROZYGOSITY,
            group = "individual's heterozygosity"))

      if (dp) {
        res$stats$ind.stats %<>% dplyr::bind_rows(
          tibble_stats(
            x = as.numeric(res$individuals$COVERAGE_TOTAL),
            group = "individual's total coverage"),
          tibble_stats(
            x = as.numeric(res$individuals$COVERAGE_MEAN),
            group = "individual's mean coverage")
        )
      }
    }

    # test <- res$stats$ind.stats
    if (m.subsample) {
      title.ind.stats <- stringi::stri_join("Individual's QC stats\n", "het subsample markers: ", length(variant.select))
    } else {
      title.ind.stats <- "Individual's QC stats"
    }

    res$figures$ind.stats.fig <- boxplot_stats(
      data = res$stats$ind.stats,
      title = title.ind.stats,
      x.axis.title = NULL,
      y.axis.title = "Statistics",
      facet.columns = TRUE,
      bp.filename = stringi::stri_join("vcf.individuals.qc_", file.date, ".pdf"),
      path.folder = path.folder)

    # test <- res$stats$ind.stats
    # test <- res$individuals

    if (!is.null(filter.individuals.missing)) {

      if (!purrr::is_double(filter.individuals.missing)) {
        outlier.id.missing <- floor(res$stats$ind.stats$OUTLIERS_HIGH[1]*100)/100
        message("Removing outlier individuals based on genotyping statistics: ", outlier.id.missing)
        filter.individuals.missing <- outlier.id.missing
      }

      blacklist.id <- res$individuals %>%
        dplyr::filter(MISSING_PROP > filter.individuals.missing) %>%
        dplyr::ungroup(.) %>%
        dplyr::distinct(INDIVIDUALS)
      blacklisted.id <- nrow(blacklist.id)
      if (blacklisted.id > 0) {
        if (verbose) message("    number of individuals blacklisted based on missing genotypes: ", blacklisted.id)
        res$blacklist.id <- blacklist.id %>%
          readr::write_tsv(x = ., path = file.path(path.folder, blacklist.id.filename))

        res$individuals %<>% dplyr::filter(!INDIVIDUALS %in% blacklist.id$INDIVIDUALS)
        # res$individuals %<>% dplyr::mutate(
        #   FILTER_INDIVIDUALS_MISSING = dplyr::if_else(
        #     INDIVIDUALS %in% blacklist.id$INDIVIDUALS, FALSE, TRUE)) %>%
        #   readr::write_tsv(x = ., path = file.path(path.folder, ind.file))

        # update the GDS
        gdsfmt::add.gdsn(
          node = radiator.gds,
          name = "individuals",
          val = res$individuals,
          replace = TRUE,
          compress = "ZIP_RA",
          closezip = TRUE)

        # if (!strata) {
        #   strata %<>% dplyr::filter(!INDIVIDUALS %in% blacklist.id$INDIVIDUALS)
        #   gdsfmt::add.gdsn(
        #     node = radiator.gds,
        #     name = "STRATA",
        #     val = strata,
        #     # val = strata$STRATA,
        #     replace = TRUE,
        #     compress = "ZIP_RA",
        #     closezip = TRUE)
        # }
        # SeqArray::seqSetFilter(
        #   object = res$vcf.connection,
        #   sample.id = res$individuals$INDIVIDUALS,
        #   verbose = FALSE)
        sync_gds(gds = res$vcf.connection, samples = res$individuals$INDIVIDUALS)


        # update filter param ----------------------------------------------------------
        filters.parameters <- update_parameters(
          parameter.obj = filters.parameters,
          gds.obj = res,
          filter.name = "Filter individuals based on missingness (with outlier stats or values)",
          param.name = "filter.individuals.missing",
          values = filter.individuals.missing)
      }
    }#filter.individuals.missing
    blacklist.id <- NULL

    # Markers stats ------------------------------------------------------------
    if (verbose) message("Updating markers metadata and stats")
    n.markers <- dplyr::n_distinct(res$markers.meta$VARIANT_ID)
    want <- c("VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS",
              "COL", "STRANDS", "REF", "ALT", "MISSING_PROP", "MISSING_POP",
              "REF_COUNT", "MAC", "SNP_PER_LOCUS")
    res$markers.meta <- suppressWarnings(
      dplyr::bind_cols(
        res$markers.meta,
        SeqArray::seqAlleleCount(
          gdsfile = res$vcf.connection,
          ref.allele = NULL,
          .progress = TRUE,
          parallel = parallel.core) %>%
          unlist(.) %>%
          matrix(
            data = .,
            nrow = n.markers, ncol = 2, byrow = TRUE,
            dimnames = list(rownames = res$markers.meta$VARIANT_ID,
                            colnames = c("REF_COUNT", "ALT_COUNT"))) %>%
          tibble::as_tibble(.),
        tibble::tibble(
          MISSING_PROP = SeqArray::seqMissing(
            gdsfile = res$vcf.connection,
            per.variant = TRUE, .progress = TRUE, parallel = parallel.core))) %>%
        dplyr::mutate(
          MAC = dplyr::if_else(ALT_COUNT < REF_COUNT, ALT_COUNT, REF_COUNT),
          ALT_COUNT = NULL)
    )


    if (res$biallelic) {
      res$markers.meta %<>%
        dplyr::group_by(LOCUS) %>%
        dplyr::mutate(SNP_PER_LOCUS = n()) %>%
        dplyr::ungroup(.)
    } else {
      message("With haplotype vcf, number of SNP/locus is counted based on the REF allele")
      message("    this stat is good only if nuc length is equal between REF and ALT haplotypes")
      message("    stacks haplotype vcf: ok")
      message("    dDocent/freebayes haplotype vcf: be careful")

      res$markers.meta %<>%
        dplyr::mutate(
          SNP_PER_LOCUS = stringi::stri_length(
            str = SeqArray::seqGetData(gdsfile = res$vcf.connection, var.name = "$ref"))
        )
    }

    # test <- res$markers.meta
    if (!strata) {
      suppressWarnings(res$markers.meta %<>% dplyr::select(dplyr::one_of(want)))
    } else {
      if (verbose) message("Calculating markers missingness per strata...")
      suppressWarnings(
        res$markers.meta %<>%
          dplyr::left_join(
            missing_per_pop(
              vcf.connection = res$vcf.connection,
              strata = dplyr::distinct(res$individuals, INDIVIDUALS, STRATA),
              parallel.core = parallel.core)
            , by = "MARKERS") %>%
          dplyr::select(dplyr::one_of(want))
      )
    }
    # test <- res$markers.meta

    # sync GDS with sample -----------------------------------------------------
    # SeqArray::seqSetFilter(
    #   object = res$vcf.connection,
    #   sample.id = res$individuals$INDIVIDUALS,
    #   verbose = FALSE)
    sync_gds(gds = res$vcf.connection, samples = res$individuals$INDIVIDUALS)


    # COMMON MARKERS -----------------------------------------------------------
    # Remove markers not in common if option selected
    if (common.markers && strata) {
      message("Scanning for markers in common...")
      not.common.markers <- dplyr::select(res$markers.meta, MARKERS, MISSING_POP) %>%
        tidyr::unnest(MISSING_POP) %>%
        dplyr::filter(MISSING_PROP == 1) %>%
        dplyr::distinct(MARKERS)
      not.in.common <- nrow(not.common.markers)

      if (not.in.common > 0) {
        message("    number of markers blacklisted: ", not.in.common)
        # update the blacklist
        res$blacklist.markers <- dplyr::bind_rows(
          res$blacklist.markers,
          dplyr::mutate(not.common.markers, FILTER = "common.markers"))
        res$markers.meta %<>% dplyr::filter(!MARKERS %in% not.common.markers$MARKERS)
        # Update GDS
        gdsfmt::add.gdsn(
          node = radiator.gds,
          name = "markers.meta",
          val = res$markers.meta,
          replace = TRUE,
          compress = "ZIP_RA",
          closezip = TRUE)
        # SeqArray::seqSetFilter(object = res$vcf.connection,
        #                        variant.id = res$markers.meta$VARIANT_ID,
        #                        verbose = FALSE)
        sync_gds(gds = res$vcf.connection, markers = res$markers.meta$VARIANT_ID)

        n.markers <- length(res$markers.meta$VARIANT_ID)

      } else {
        message("    all markers in common between strata")
      }
      # update filter param ----------------------------------------------------------
      filters.parameters <- update_parameters(
        parameter.obj = filters.parameters,
        gds.obj = res,
        filter.name = "Filter for common markers",
        param.name = "common.markers",
        values = common.markers)
      not.in.common <- not.common.markers <- NULL
    }

    # FILTER MAC----------------------------------------------------------------
    if (!is.null(filter.mac)) {
      blacklist.mac <- dplyr::filter(res$markers.meta, MAC < filter.mac) %>%
        dplyr::select(MARKERS)

      if (nrow(blacklist.mac) > 0) {
        message("Number of markers with low MAC blacklisted: ", nrow(blacklist.mac))
        # update the blacklist
        res$blacklist.markers <- dplyr::bind_rows(
          res$blacklist.markers,
          dplyr::mutate(blacklist.mac, FILTER = "filter.mac"))
        # check <- res$blacklist.markers

        # Update markers.meta
        res$markers.meta %<>% dplyr::filter(!MARKERS %in% blacklist.mac$MARKERS)
        # check <- res$markers.meta

        # Update GDS
        gdsfmt::add.gdsn(
          node = radiator.gds,
          name = "markers.meta",
          val = res$markers.meta,
          replace = TRUE,
          compress = "ZIP_RA",
          closezip = TRUE)

        # SeqArray::seqSetFilter(object = res$vcf.connection,
        #                        variant.id = res$markers.meta$VARIANT_ID,
        #                        verbose = FALSE)
        sync_gds(gds = res$vcf.connection, markers = res$markers.meta$VARIANT_ID)

        n.markers <- length(res$markers.meta$VARIANT_ID)
      }
      blacklist.mac <- NULL
      # update filter param ----------------------------------------------------------
      filters.parameters <- update_parameters(
        parameter.obj = filters.parameters,
        gds.obj = res,
        filter.name = "Filter for low minor allele count",
        param.name = "filter.mac",
        values = filter.mac)
    }
    # test <- res$markers.meta

    # Coverage Stats -----------------------------------------------------------
    if (dp) {
      if (verbose) message("Generating coverage stats")
      coverage.info <- extract_coverage(data = res$vcf.connection, ind = FALSE)

      res$stats$coverage.stats <- tibble_stats(
        x = coverage.info$markers.mean,
        group =  "unfiltered")


      res$stats$coverage.stats <- dplyr::bind_rows(
        res$stats$coverage.stats,
        tibble_stats(
          x = coverage.info$markers.mean[coverage.info$markers.mean >= res$stats$coverage.stats$OUTLIERS_LOW &
                                           coverage.info$markers.mean <= res$stats$coverage.stats$OUTLIERS_HIGH],
          group = "filtered outliers")
      ) %>%
        dplyr::mutate(GROUP = factor(
          x = GROUP,
          levels = c("unfiltered", "filtered outliers"))) %>%
        dplyr::arrange(GROUP)

      # Generate box plot
      res$figures$coverage.stats.fig <- boxplot_stats(
        data = res$stats$coverage.stats,
        title = "Impact of filtering outliers\non coverage statistics",
        x.axis.title = NULL,
        y.axis.title = "Markers coverage (mean read depth)",
        facet.columns = TRUE,
        bp.filename = stringi::stri_join("coverage.statistics_", file.date, ".pdf"),
        path.folder = path.folder)

      res$markers.meta$COVERAGE_MEAN <- coverage.info$markers.mean

      if (filter.coverage.outliers) {
        blacklist.coverage <- res$markers.meta %>%
          dplyr::filter(
            COVERAGE_MEAN < res$stats$coverage.stats$OUTLIERS_LOW[1] |
              COVERAGE_MEAN > res$stats$coverage.stats$OUTLIERS_HIGH[1]) %>%
          dplyr::select(MARKERS)

        if (nrow(blacklist.coverage) > 0) {
          message("    number of markers with outlier coverage blacklisted: ", nrow(blacklist.coverage))
          # update the blacklist
          res$blacklist.markers <- dplyr::bind_rows(
            res$blacklist.markers,
            dplyr::mutate(blacklist.coverage, FILTER = "filter.coverage.outliers"))
          # check <- res$blacklist.markers

          # Update markers.meta
          res$markers.meta %<>% dplyr::filter(!MARKERS %in% blacklist.coverage$MARKERS)
          # check <- res$markers.meta

          # Update GDS
          gdsfmt::add.gdsn(
            node = radiator.gds,
            name = "markers.meta",
            val = res$markers.meta,
            replace = TRUE,
            compress = "ZIP_RA",
            closezip = TRUE)

          # SeqArray::seqSetFilter(object = res$vcf.connection,
          #                        variant.id = res$markers.meta$VARIANT_ID,
          #                        verbose = FALSE)
          sync_gds(gds = res$vcf.connection, markers = res$markers.meta$VARIANT_ID)


        }
        blacklist.coverage <- NULL
      }
      # update filter param ----------------------------------------------------------
      filters.parameters <- update_parameters(
        parameter.obj = filters.parameters,
        gds.obj = res,
        filter.name = "Filter for low or high coverage markers",
        param.name = "filter.coverage.outliers: low / high",
        values = paste(res$stats$coverage.stats$OUTLIERS_LOW[1], res$stats$coverage.stats$OUTLIERS_HIGH[1], sep = " / "))
      coverage.info <- NULL
    }

    # Filter markers genotyping/missing ----------------------------------------
    if (strata) {
      res$figures$missing.markers.fig <- markers_genotyped_helper(
        x = dplyr::select(res$markers.meta, MARKERS, MISSING_POP) %>%
          tidyr::unnest(MISSING_POP) %>%
          dplyr::rename(PERCENT = MISSING_PROP, POP_ID = STRATA) %>%
          dplyr::mutate(PERCENT = PERCENT * 100),
        y = dplyr::select(res$markers.meta, MARKERS, PERCENT = MISSING_PROP) %>%
          dplyr::mutate(PERCENT = PERCENT * 100),
        overall.only = FALSE
      )
      n.pop <- dplyr::n_distinct(res$individuals$STRATA) + 2
      ggplot2::ggsave(
        plot = res$figures$missing.markers.fig,
        filename = file.path(path.folder, stringi::stri_join("plot.markers.genotyping.rate_", file.date, ".pdf")),
        width = n.pop * 10, height = 10,
        dpi = 300, units = "cm", useDingbats = FALSE, limitsize = FALSE)
    } else {#overall only
      res$figures$missing.markers.fig <- markers_genotyped_helper(
        x = NULL,
        y = dplyr::select(res$markers.meta, MARKERS, PERCENT = MISSING_PROP) %>%
          dplyr::mutate(PERCENT = PERCENT * 100),
        overall.only = TRUE
      )
      ggplot2::ggsave(
        plot = res$figures$missing.markers.fig,
        filename = file.path(path.folder, stringi::stri_join("plot.markers.genotyping.rate_", file.date, ".pdf")),
        width = 15, height = 10,
        dpi = 300, units = "cm", useDingbats = FALSE)
    }

    if (!is.null(filter.markers.missing)) {
      res$markers.meta %<>%
        dplyr::mutate(
          MARKERS_MISSING = dplyr::if_else(MISSING_PROP <= filter.markers.missing / 100, TRUE, FALSE))

      blacklist.markers.missing <- res$markers.meta %>%
        dplyr::filter(MISSING_PROP > filter.markers.missing / 100) %>%
        dplyr::select(MARKERS)

      if (nrow(blacklist.markers.missing) > 0) {
        message("    number of markers blacklisted based on missing individuals/genotypes: ", nrow(blacklist.markers.missing))
        # update the blacklist
        res$blacklist.markers <- dplyr::bind_rows(
          res$blacklist.markers,
          dplyr::mutate(blacklist.markers.missing, FILTER = "filter.markers.missing"))
        # check <- res$blacklist.markers

        # Update markers.meta
        res$markers.meta %<>% dplyr::filter(!MARKERS %in% blacklist.markers.missing$MARKERS)
        # check <- res$markers.meta

        # Update GDS
        gdsfmt::add.gdsn(
          node = radiator.gds,
          name = "markers.meta",
          val = res$markers.meta,
          replace = TRUE,
          compress = "ZIP_RA",
          closezip = TRUE)

        # SeqArray::seqSetFilter(object = res$vcf.connection,
        #                        variant.id = res$markers.meta$VARIANT_ID,
        #                        verbose = FALSE)
        sync_gds(gds = res$vcf.connection, markers = res$markers.meta$VARIANT_ID)

      }
      blacklist.markers.missing <- NULL
      # update filter param ----------------------------------------------------------
      filters.parameters <- update_parameters(
        parameter.obj = filters.parameters,
        gds.obj = res,
        filter.name = "Filter markers based on missingness",
        param.name = "filter.markers.missing",
        values = filter.markers.missing)
    }
    # test <- res$markers.meta

    # position of the SNPs on the read -----------------------------------------
    if (tibble::has_name(res$markers.meta, "COL")) {
      if (verbose) message("Generating SNP position on read stats")
      res$stats$snp.col.stats <- tibble_stats(
        x = dplyr::filter(res$markers.meta) %>% # using unfiltered dataset
          dplyr::distinct(MARKERS,COL) %$%
          COL,
        group = "snp position on read")
      snp.col.iqr.threshold <- c(res$stats$snp.col.stats$Q25, res$stats$snp.col.stats$Q75)

      # outliers
      if ("outliers" %in% filter.snp.read.position) {
        res$markers.meta %<>%
          dplyr::mutate(
            SNP_POS_READ_OUTLIERS = dplyr::if_else(
              COL < res$stats$snp.col.stats$OUTLIERS_HIGH, TRUE, FALSE))
      } else {
        res$markers.meta %<>% dplyr::mutate(SNP_POS_READ_OUTLIERS = TRUE)
      }

      # Q75
      if ("q75" %in% filter.snp.read.position) {
        res$markers.meta %<>%
          dplyr::mutate(
            SNP_POS_READ_Q75 = dplyr::if_else(
              COL <= snp.col.iqr.threshold[2], TRUE, FALSE))
      } else {
        res$markers.meta %<>% dplyr::mutate(SNP_POS_READ_Q75 = TRUE)
      }

      # IQR
      if ("iqr" %in% filter.snp.read.position) {
        res$markers.meta %<>%
          dplyr::mutate(
            SNP_POS_READ_IQR = dplyr::if_else(
              COL >= snp.col.iqr.threshold[1] & COL <= snp.col.iqr.threshold[2], TRUE, FALSE))
      } else {
        res$markers.meta %<>% dplyr::mutate(SNP_POS_READ_IQR = TRUE)
      }

      # Update with ALL filters available
      res$stats$snp.col.stats %<>% dplyr::bind_rows(
        tibble_stats(
          x = unique.col <- res$markers.meta %>%
            dplyr::filter(SNP_POS_READ_OUTLIERS) %>% # outlier snp on read filter
            dplyr::filter(SNP_POS_READ_Q75) %>% # snp below Q75
            dplyr::filter(SNP_POS_READ_IQR) %>% # focus on snp within IQR range
            dplyr::distinct(MARKERS, COL) %$% COL,
          group = "snp position on read filtered")) #%>% dplyr::select(-OUTLIERS_LOW))
      snp.col.iqr.threshold <- NULL
      # Generate box plot
      res$figures$snp.pos.read.fig <- boxplot_stats(
        data = res$stats$snp.col.stats,
        title = "Impact of filter\non the SNP position on the read",
        x.axis.title = NULL,
        y.axis.title = "SNP position (base pair)",
        facet.columns = TRUE,
        bp.filename = stringi::stri_join("vcf.snp.position.read_", file.date, ".pdf"),
        path.folder = path.folder)

      # Filtering ...
      if (!is.null(filter.snp.read.position)) {
        blacklist.snp.read.position <- res$markers.meta %>%
          dplyr::filter(!SNP_POS_READ_OUTLIERS | !SNP_POS_READ_Q75 | !SNP_POS_READ_IQR) %>%
          dplyr::select(MARKERS)

        if (nrow(blacklist.snp.read.position) > 0) {
          message("    number of markers blacklisted based on position on the read: ", nrow(blacklist.snp.read.position))
          # update the blacklist
          res$blacklist.markers <- dplyr::bind_rows(
            res$blacklist.markers,
            dplyr::mutate(blacklist.snp.read.position, FILTER = "filter.snp.read.position"))
          # check <- res$blacklist.markers

          # Update markers.meta
          res$markers.meta %<>% dplyr::filter(SNP_POS_READ_OUTLIERS) %>%
            dplyr::filter(SNP_POS_READ_Q75) %>%
            dplyr::filter(SNP_POS_READ_IQR) %>%
            dplyr::select(-c(SNP_POS_READ_OUTLIERS, SNP_POS_READ_Q75, SNP_POS_READ_IQR))
          # check <- res$markers.meta

          # Update GDS
          gdsfmt::add.gdsn(
            node = radiator.gds,
            name = "markers.meta",
            val = res$markers.meta,
            replace = TRUE,
            compress = "ZIP_RA",
            closezip = TRUE)

          # SeqArray::seqSetFilter(object = res$vcf.connection,
          #                        variant.id = res$markers.meta$VARIANT_ID,
          #                        verbose = FALSE)
          sync_gds(gds = res$vcf.connection, markers = res$markers.meta$VARIANT_ID)

        }
        blacklist.snp.read.position <- NULL
        # update filter param ----------------------------------------------------------
        filters.parameters <- update_parameters(
          parameter.obj = filters.parameters,
          gds.obj = res,
          filter.name = "Filter SNPs based on position on read",
          param.name = "filter.snp.read.position",
          values = paste(filter.snp.read.position, sep = " / "))
      }
      # test <- res$markers.meta
    }
    # number of SNPs per locus -------------------------------------------------
    # Note to myself: there is no filtering here, just an output of the stats...
    # before and after filtering
    res$stats$snp.locus.stats <- tibble_stats(
      x = dplyr::distinct(res$markers.meta, LOCUS, SNP_PER_LOCUS) %$%
        SNP_PER_LOCUS,
      group =  "filtered")

    # read.length <- max(res$markers.meta$COL)
    # Generate box plot
    res$figures$snp.per.locus.fig <- boxplot_stats(
      data = res$stats$snp.locus.stats,
      title = "Number of SNPs per locus",
      x.axis.title = NULL,
      y.axis.title = "Number of SNPs per locus",
      bp.filename = stringi::stri_join("vcf.number.snp.per.locus_", file.date, ".pdf"),
      path.folder = path.folder)

    # generate a variant.id vector ---------------------------------------------
    variant.id.select <- res$markers.meta$VARIANT_ID

    # SNP LD -------------------------------------------------------------------
    if (!is.null(filter.short.ld)) {

      # data = res$vcf.connection
      # snp.ld = filter.short.ld
      # maf.data = NULL
      # ld.threshold = filter.long.ld
      # filename = NULL
      # keep.gds <- TRUE
      # ld.figures <- TRUE
      # long.ld.missing = TRUE
      # ld.method = "r2"

      filter.ld <- snp_ld(
        data = res$vcf.connection,
        snp.ld = filter.short.ld,
        maf.data = NULL,
        ld.threshold = filter.long.ld,
        parallel.core = parallel.core,
        filename = NULL,
        verbose = verbose,
        long.ld.missing = long.ld.missing,
        ld.method = ld.method,
        path.folder = path.folder
      )
      # names(filter.ld)
      blacklist.snp.ld <- res$markers.meta %>%
        dplyr::filter(!MARKERS %in% filter.ld$whitelist.snp.ld$MARKERS) %>%
        dplyr::select(MARKERS)

      if (nrow(blacklist.snp.ld) > 0) {
        message("Total number of markers blacklisted based LD: ", nrow(blacklist.snp.ld))
        # update the blacklist
        res$blacklist.markers <- dplyr::bind_rows(
          res$blacklist.markers,
          dplyr::mutate(blacklist.snp.ld, FILTER = "filter LD (short/long)"))
        # check <- res$blacklist.markers

        # Update markers.meta
        res$markers.meta %<>%
          dplyr::filter(MARKERS %in% filter.ld$whitelist.snp.ld$MARKERS)
        # check <- res$markers.meta

        # Update GDS
        gdsfmt::add.gdsn(
          node = radiator.gds,
          name = "markers.meta",
          val = res$markers.meta,
          replace = TRUE,
          compress = "ZIP_RA",
          closezip = TRUE)

        sync_gds(gds = res$vcf.connection, markers = res$markers.meta$VARIANT_ID)
        # SeqArray::seqSetFilter(object = res$vcf.connection,
        #                        variant.id = res$markers.meta$VARIANT_ID,
        #                        verbose = FALSE)
      }
      blacklist.snp.ld <- filter.ld <- NULL

      # sync GDS with sample -----------------------------------------------------
      # SeqArray::seqSetFilter(
      #   object = res$vcf.connection,
      #   sample.id = res$individuals$INDIVIDUALS,
      #   verbose = FALSE)
      sync_gds(gds = res$vcf.connection, samples = res$individuals$INDIVIDUALS)

      # update filter param ----------------------------------------------------------
      filters.parameters <- update_parameters(
        parameter.obj = filters.parameters,
        gds.obj = res,
        filter.name = "Filter for short and long LD",
        param.name = "filter.short.ld / filter.long.ld / long.ld.missing / ld.method",
        values = paste(filter.short.ld, filter.long.ld, long.ld.missing, ld.method, sep = " / "))
    }

    # For summary at the end of the function:
    res$stats$ind.missing <- round(mean(res$individuals$MISSING_PROP, na.rm = TRUE), 2)
    if (dp) res$stats$ind.cov.total <- round(mean(res$individuals$COVERAGE_TOTAL, na.rm = TRUE), 0)
    if (dp) res$stats$ind.cov.mean <- round(mean(res$individuals$COVERAGE_MEAN, na.rm = TRUE), 0)
    res$stats$markers.missing <- round(mean(res$markers.meta$MISSING_PROP, na.rm = TRUE), 2)
    if (dp) res$stats$markers.cov <- round(mean(res$markers.meta$COVERAGE_MEAN, na.rm = TRUE), 0) # same as above because NA...
  }#End stats

  # Final Sync GDS -------------------------------------------------------------------

  # blacklist.id
  gdsfmt::add.gdsn(
    node = radiator.gds,
    name = "blacklist.id",
    val = res$blacklist.id,
    replace = TRUE,
    compress = "ZIP_RA",
    closezip = TRUE)

  # blacklist.markers
  gdsfmt::add.gdsn(
    node = radiator.gds,
    name = "blacklist.markers",
    val = res$blacklist.markers,
    replace = TRUE,
    compress = "ZIP_RA",
    closezip = TRUE)

  sync_gds(gds = res$vcf.connection, samples = res$individuals$INDIVIDUALS,
           markers = res$markers.meta$VARIANT_ID)



  # RESULTS --------------------------------------------------------------------
  number.info <- SeqArray::seqGetFilter(gdsfile = res$vcf.connection)
  res$n.individuals <- length(number.info$sample.sel[number.info$sample.sel])
  res$n.markers <- length(number.info$variant.sel[number.info$variant.sel])
  res$n.chromosome <- dplyr::n_distinct(res$markers.meta$CHROM)
  res$n.locus <- dplyr::n_distinct(res$markers.meta$LOCUS)
  res$filename <- filename

  if (verbose) {
    if (vcf.stats) {
      message("\n\nSummary:\nMissing data (averaged): ")
      message("    markers: ", res$stats$markers.missing)
      message("    individuals: ", res$stats$ind.missing)
      if (dp) message("\n\nCoverage info:")
      if (dp) message("    individuals mean read depth: ", res$stats$ind.cov.total)
      if (dp) message("    individuals mean genotype coverage: ", res$stats$ind.cov.mean)
      if (dp) message("    markers mean coverage: ", res$stats$markers.cov)
    }
    message("\n\nNumber of chromosome/contig/scaffold: ", res$n.chromosome)
    message("Number of locus: ", res$n.locus)
    message("Number of markers: ", res$n.markers)
    message("Number of individuals: ", res$n.individuals)
  }
  # prepare data for the SeqVarTools package
  # seqOptimize("tmp.gds", target="by.sample")
  timing.import <- proc.time() - timing.import
  if (verbose) message("\nWorking time: ", round(timing.import[[3]]), " sec\n")
  return(res)
} # End write_SeqArray

# Internal nested Function -----------------------------------------------------
#' @title extract_info
#' @description Extract vcf information
#' @rdname extract_info
#' @keywords internal
#' @export
extract_info <- function(vcf) {
  res <- list()
  # print(vcf, all=TRUE, attribute=TRUE)
  # tic()
  vcf.info <- SeqArray::seqVCF_Header(vcf.fn = vcf, getnum = TRUE)
  # toc()
  res$vcf.source <- vcf.info$header$value[2]
  res$n.ind <- vcf.info$num.sample
  res$n.markers <- vcf.info$num.variant
  res$sample.id <- vcf.info$sample.id
  return(res)
}#End extract_info

#' @title tibble_stats
#' @description Generate a tibble of statistics
#' @rdname tibble_stats
#' @keywords internal
#' @export

tibble_stats <- function(x, group, subsample = NULL) {
  if (is.null(subsample)) subsample <- 1L
  Q <- stats::quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  res <- tibble::tibble(
    GROUP = group,
    MIN = min(x, na.rm = TRUE),
    # Q25 = stats::quantile(x, 0.25, na.rm = TRUE),
    Q25 = Q[1],
    MEDIAN = stats::median(x, na.rm = TRUE),
    Q75 = Q[2],
    MAX = max(x, na.rm = TRUE),
    IQR = abs(diff(Q)),
    # IQR = stats::IQR(depth, na.rm = TRUE),
    OUTLIERS_LOW = Q25 - (1.5 * IQR),
    OUTLIERS_HIGH =  Q75 + (1.5 * IQR)) %>%
    dplyr::mutate_if(.tbl = ., .predicate = is.integer, .funs = as.numeric) %>%
    dplyr::mutate(
      OUTLIERS_LOW = dplyr::if_else(OUTLIERS_LOW < MIN, MIN, OUTLIERS_LOW),
      OUTLIERS_HIGH = dplyr::if_else(OUTLIERS_HIGH > MAX, MAX, OUTLIERS_HIGH),
      SUBSAMPLE = subsample
    )
  Q <- NULL
  return(res)
}#End tibble_stats

#' @title extract_coverage
#' @description Extract coverage information from a GDS file
#' @rdname extract_coverage
#' @keywords internal
#' @export
extract_coverage <- function(data = NULL, markers = TRUE, ind = TRUE) {#, variant.id.select = NULL) {
  # data <- res$vcf.connection
  coverage.info <- list()
  depth <- SeqArray::seqGetData(data, "annotation/format/DP")

  # test <- depth$data
  if (ind) {
    coverage.info$ind.cov.tot <- as.integer(round(rowSums(x = depth$data, na.rm = TRUE, dims = 1L), 0))
    coverage.info$ind.cov.mean <- as.integer(round(rowMeans(x = depth$data, na.rm = TRUE, dims = 1L), 0))
  }

  if (markers) {
    coverage.info$markers.mean <- as.integer(round(colMeans(x = depth$data, na.rm = TRUE, dims = 1L), 0))
  }
  # coverage.info$ind.cov.tot <- as.integer(round(rowSums(x = depth$data[,variant.id.select], na.rm = TRUE, dims = 1L), 0))
  # coverage.info$ind.cov.mean <- as.integer(round(rowMeans(x = depth$data[,variant.id.select], na.rm = TRUE, dims = 1L), 0))
  return(coverage.info)
}#End extract_coverage

# boxplot of stats
#' @title boxplot_stats
#' @description Generate a boxplot
#' @rdname boxplot_stats
#' @keywords internal
#' @export
boxplot_stats <- function(data,
                          title,
                          x.axis.title = NULL,
                          y.axis.title,
                          facet.columns = FALSE,
                          facet.rows = FALSE,
                          bp.filename = NULL,
                          path.folder = NULL) {
  # data <- test
  # x.axis.title = NULL
  # x.axis.title <- "SNP position on the read groupings"
  # title <- "Individual's QC stats"
  # title <- "Impact of filter on the SNP position in the read"
  # y.axis.title <- "Statistics"
  # y.axis.title <- "SNP position (base pair)"
  # bp.filename <- "vcf.snp.position.read.pdf"
  # bp.filename <- "test.pdf"
  # facet.columns = TRUE
  # facet.rows = FALSE
  # path.folder = NULL

  if (is.null(path.folder)) path.folder <- getwd()

  n.group <- dplyr::n_distinct(data$GROUP)
  element.text <- ggplot2::element_text(size = 10,
                                        family = "Helvetica", face = "bold")

  if (facet.columns) {
    data <- dplyr::mutate(data, X = "1")
    fig.boxplot <- ggplot2::ggplot(data = data, ggplot2::aes(X))
  } else {
    fig.boxplot <- ggplot2::ggplot(data = data, ggplot2::aes(GROUP))
  }


  fig.boxplot <- fig.boxplot +
    ggplot2::geom_boxplot(
      ggplot2::aes(ymin = OUTLIERS_LOW, lower = Q25, middle = MEDIAN, upper = Q75,
                   ymax = OUTLIERS_HIGH), stat = "identity") +
    ggplot2::labs(
      y = y.axis.title,
      title = title)

  # Draw upper outliers
  if (facet.columns) {
    fig.boxplot <- fig.boxplot +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = "1", xend = "1",
          y = OUTLIERS_HIGH, yend = MAX),
        linetype = "dashed")
  } else {
    fig.boxplot <- fig.boxplot +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = GROUP, xend = GROUP,
          y = OUTLIERS_HIGH, yend = MAX),
        linetype = "dashed")
  }

  # Draw lower outliers
  if (facet.columns) {
    fig.boxplot <- fig.boxplot +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = "1", xend = "1",
          y = OUTLIERS_LOW, yend = MIN),
        linetype = "dashed")
  } else {
    fig.boxplot <- fig.boxplot +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = GROUP, xend = GROUP,
          y = OUTLIERS_LOW, yend = MIN),
        linetype = "dashed")
  }

  fig.boxplot <- fig.boxplot +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, family = "Helvetica",
                                         face = "bold", hjust = 0.5),
      legend.position = "none",
      axis.title.y = element.text,
      axis.text.y = element.text
    ) +
    ggplot2::theme_bw()

  if (is.null(x.axis.title)) {
    fig.boxplot <- fig.boxplot +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  } else {
    fig.boxplot <- fig.boxplot +
      ggplot2::xlab(x.axis.title) +
      ggplot2::theme(
        axis.title.x = element.text,
        # axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  }

  if (facet.columns) {
    fig.boxplot <- fig.boxplot + ggplot2::facet_grid(GROUP ~ ., scales = "free")
    n.facet <- n.group * 2
    width <- 10
    height <- 10 + (4 * n.group)
  }

  # else {
  # width <-  10 + (5 * n.group) + 1
  # height <-  10
  # }

  if (facet.rows) {
    fig.boxplot <- fig.boxplot + ggplot2::facet_grid(FACET_ROWS ~ ., scales = "free")
    n.facet <- n.group * 2
    width <- 10
    height <- 10 + (4 * n.group)
  }

  if (!facet.rows && !facet.columns) {
    width <-  10 + (5 * n.group) + 1
    height <-  10
  }

  print(fig.boxplot)
  if (!is.null(bp.filename)) {
    suppressMessages(ggplot2::ggsave(
      filename = file.path(path.folder, bp.filename),
      plot = fig.boxplot,
      width = width,
      height = height,
      dpi = 300, units = "cm", useDingbats = FALSE))
  }
  return(fig.boxplot)
}#Endboxplot_stats

# missingness per markers per pop
#' @title missing_per_pop
#' @description Generate missingness per markers per pop
#' @rdname missing_per_pop
#' @keywords internal
#' @export
missing_per_pop <- function(
  vcf.connection,
  strata,
  parallel.core = parallel::detectCores() - 1
) {
  # vcf.connection <- res$vcf.connection

  missing_pop <- function(
    id.select, vcf.connection,
    parallel.core = parallel::detectCores() - 1
  ) {
    SeqArray::seqSetFilter(object = vcf.connection,
                           sample.id = id.select$INDIVIDUALS,
                           action = "set",
                           verbose = FALSE)
    res <- tibble::tibble(
      MARKERS = gdsfmt::read.gdsn(gdsfmt::index.gdsn(node = vcf.connection, path = "radiator/markers.meta/MARKERS")),
      MISSING_PROP = SeqArray::seqMissing(
        gdsfile = vcf.connection,
        per.variant = TRUE, .progress = TRUE, parallel = parallel.core - 1))
    SeqArray::seqResetFilter(
      object = vcf.connection, sample = TRUE, variant = FALSE, verbose = FALSE)
    return(res)
  }#End missing_pop

  res <- strata %>%
    dplyr::group_by(STRATA) %>%
    tidyr::nest(data = ., .key = id.select) %>%
    dplyr::mutate(MISSING_POP = purrr::map(
      .x = .$id.select,
      .f = missing_pop,
      vcf.connection = vcf.connection,
      parallel.core = parallel.core)) %>%
    tidyr::unnest(data = ., MISSING_POP) %>%
    dplyr::mutate(STRATA = as.character(STRATA)) %>%
    dplyr::group_by(MARKERS) %>%
    tidyr::nest(data = ., .key = MISSING_POP)
  return(res)
}#End missing_per_pop


# Sync GDS
#' @title sync_gds
#' @description Synchronize gds with samples and markers
#' @rdname sync_gds
#' @keywords internal
#' @export
sync_gds <- function(gds, samples = NULL, markers = NULL, verbose = FALSE) {
  SeqArray::seqSetFilter(
    object = gds,sample.id = samples, variant.id = markers, verbose = verbose)
}#End sync_gds


# summarize GDS
#' @title summary_gds
#' @description Summary of gds object: number of samples and markers
#' @rdname summary_gds
#' @keywords internal
#' @export
summary_gds <- function(gds) {
  message("GDS summary: ")
  check <- SeqArray::seqGetFilter(gds)
  message("    number of samples: ", length(check$sample.sel[check$sample.sel]))
  message("    number of markers: ", length(check$variant.sel[check$variant.sel]))
}


