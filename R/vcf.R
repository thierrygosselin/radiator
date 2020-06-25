# Read VCF with SeqArray--------------------------------------------------------
# write a SeqArray object from a tidy data frame
#' @name read_vcf
#' @title Read VCF files and write a GDS file
#' @description The function reads VCF files for radiator and
#' generate a connection SeqArray \href{https://github.com/zhengxwen/SeqArray}{SeqArray}
#' GDS object/file of class \code{SeqVarGDSClass} (Zheng et al. 2017)
#' The Genomic Data Structure (GDS) file format is detailed in
#' \href{https://github.com/zhengxwen/gdsfmt}{gdsfmt}.
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
#' \item \strong{DArT VCFs:} CHROM with \code{.} are replaced by \code{denovo}.
#' \code{POS} with \code{NA} are replaced by {50}.
#' \code{COL}: the position of the SNP on the read is extracted from the
#' \code{LOCUS} as any other DArT data.
#' \code{LOCUS}: are the first group of digits the rest of the info is discarded
#' and \code{LOCUS} is then joined with \code{POS} by a \code{_}.
#' \code{POS} is replaced by \code{COL}, (col info is duplicated).
#' This is necessary to make the lines in the VCF unique.
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

#' @section VCF file format:
#'
#' \strong{PLINK:} radiator fills the \code{LOCUS} column of PLINK VCFs with
#' a unique integer based on the \code{CHROM} column
#' (\code{as.integer(factor(x = CHROM))}).
#' The \code{COL} column is filled with 1L for lack of bettern info on this.
#' Not what you need ? Open an issue on GitHub for a request.
#'
#' \strong{ipyrad:} the pattern \code{locus_} in the \code{CHROM} column
#' is removed and used. The \code{COL} column is filled with the same value as
#' \code{POS}.
#'
#' \strong{GATK:} Some VCF have an \code{ID} column filled with \code{.},
#' the LOCUS information is all contained along the linkage group in the
#' \code{CHROM} column. To make it work with
#' \href{https://github.com/thierrygosselin/radiator}{radiator},
#' the \code{ID} column is filled with the \code{POS} column info.
#'
#' \strong{platypus:} Some VCF files don't have an ID filed with values,
#' here the same thing is done as GATK VCF files above.
#'
#' \strong{freebayes:} Some VCF files don't have an ID filed with values,
#' here the same thing is done as GATK VCF files above.
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

#' @section Advance mode:
#'
#' \emph{dots-dots-dots ...} allows to pass several arguments for fine-tuning the function:
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
#' \item \code{filter.mac}: (integer) Default: \code{filter.mac = NULL}.
#' To blacklist markers below a specific Minor Allele Count (calculated overall/global).
#' Documented in \code{\link{filter_mac}}.
#'
#' \item \code{filter.coverage}: (optional, string)
#' Default: \code{filter.coverage = NULL}. To blacklist markers based on mean coverage.
#' Documented in \code{\link{filter_coverage}}.
#'
#' \item \code{filter.genotyping}: (integer) To blacklist markers with too
#' many missing data. e.g. \code{filter.genotyping = 0.1}, will only keep
#' markers with missing rate <= to 10 percent.
#' Documented in \code{\link{filter_genotyping}}.
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
#' The function returns a GDS object.

#' @export
#' @rdname read_vcf

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

#' @references Zheng X, Gogarten S, Lawrence M, Stilp A, Conomos M, Weir BS,
#' Laurie C, Levine D (2017). SeqArray -- A storage-efficient high-performance
#' data format for WGS variant calls.
#' Bioinformatics.
#'
#' @references Danecek P, Auton A, Abecasis G et al. (2011)
#' The variant call format and VCFtools.
#' Bioinformatics, 27, 2156-2158.

#' @seealso
#' \code{\link{filter_ld}}
#' \code{\link[radiator]{tidy_genomic_data}}
#' \code{\link[radiator]{tidy_vcf}}

#' @examples
#' \dontrun{
#' #require(SeqVarTools)
#' # with built-in defaults:
#'  prep.data <- radiator::read_vcf(data = "populations.snps.vcf")
#'
#' # Using more arguments and filters (recommended):
#' prep.data <- radiator::read_vcf(
#'     data = "populations.snps.vcf",
#'     strata = "strata_salamander.tsv",
#'     path.folder = "salamander",
#'     filter.individuals.missing = "outliers",
#'     filter.common.markers = TRUE,
#'     filter.strands = "blacklist",
#'     filter.mac = 4,
#'     filter.genotyping = 0.3,
#'     filter.snp.position.read = "outliers",
#'     filter.short.ld = "mac",
#'     filter.long.ld = NULL,
#'     verbose = TRUE
#'     )
#' }

read_vcf <- function(
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
  # internal <- FALSE
  # random.seed <- NULL
  # markers.info = NULL
  # vcf.stats <- TRUE
  # vcf.metadata = TRUE
  # filter.strands = "blacklist"
  # filter.coverage = NULL
  # filter.genotyping <- NULL
  # filter.snp.position.read <- NULL
  # filter.mac <- NULL
  # filter.common.markers = FALSE
  # filter.monomorphic <- FALSE
  # filter.short.ld <- NULL
  # filter.long.ld <- NULL
  # long.ld.missing <- TRUE
  # ld.method <- "r2"
  # blacklist.id = NULL
  # pop.select = NULL
  # pop.levels = NULL
  # pop.labels = NULL
  # whitelist.markers = NULL
  # subsample.markers.stats = 0.2
  # parameters <- NULL
  # filter.individuals.missing <- NULL
  # filter.individuals.coverage.total <- NULL
  # filter.snp.number <- NULL

  ## With filters
  # filter.individuals.missing <- "outliers"
  # filter.individuals.coverage.total <- "outliers"
  # filter.short.ld <- "mac"
  # filter.long.ld <- 0.3
  # filter.common.markers = TRUE
  # filter.monomorphic <- TRUE
  # filter.mac <- 4
  # filter.coverage = c(5, 150)
  # filter.genotyping <- 0.15
  # filter.snp.position.read <- "outliers"
  # filter.snp.number <- "outliers"

  # Cleanup---------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (verbose) message("Execution date@time: ", file.date)
  old.dir <- getwd()
  opt.change <- getOption("width")
  options(future.globals.maxSize= Inf)
  options(width = 70)
  timing <- radiator_tic()
  res <- list()
  #back to the original directory and options
  on.exit(setwd(old.dir), add = TRUE)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(radiator_toc(timing), add = TRUE)
  on.exit(radiator_function_header(f.name = "read_vcf", start = FALSE, verbose = verbose), add = TRUE)
  # Required package -----------------------------------------------------------
  radiator_packages_dep(package = "SeqVarTools", cran = FALSE, bioc = TRUE)

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
                "blacklist.id", "pop.select", "pop.levels", "pop.labels",
                "path.folder",
                "markers.info", "vcf.metadata",
                "subsample.markers.stats", "parameters", "random.seed", "internal"),
    verbose = FALSE
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
  wf <- path.folder <- generate_folder(
    f = path.folder,
    rad.folder = "read_vcf",
    prefix_int = FALSE,
    internal = internal,
    file.date = file.date,
    verbose = verbose)

  radiator.folder <- generate_folder(
    f = path.folder,
    rad.folder = "import_gds",
    prefix_int = TRUE,
    internal = internal,
    file.date = file.date,
    verbose = verbose)

  # write the dots file
  write_rad(
    data = rad.dots,
    path = radiator.folder,
    filename = stringi::stri_join("radiator_read_vcf_args_", file.date, ".tsv"),
    tsv = TRUE,
    internal = internal,
    write.message = "Function call and arguments stored in: ",
    verbose = verbose
  )


  if (is.null(filename)) {
    ind.file <- stringi::stri_join("vcf_individuals_info_", file.date, ".tsv")
    markers.file <- stringi::stri_join("vcf_markers_metadata_", file.date, ".tsv")
    blacklist.markers <- stringi::stri_join("blacklist.markers_", file.date, ".tsv")
    blacklist.id.filename <- stringi::stri_join("blacklist.individuals_", file.date, ".tsv")
  } else {
    ind.file <- stringi::stri_join(filename, "_vcf_individuals_info_", file.date, ".tsv")
    markers.file <- stringi::stri_join(filename, "_vcf_markers_metadata_", file.date, ".tsv")
    blacklist.markers <- stringi::stri_join(filename, "_blacklist.markers_", file.date, ".tsv")
    blacklist.id.filename <- stringi::stri_join(filename, "_blacklist.individuals_", file.date, ".tsv")
  }

  filename <- generate_filename(
    name.shortcut = filename,
    path.folder = radiator.folder,
    extension = "gds")

  filename.short <- filename$filename.short
  filename <- filename$filename

  # Random seed ----------------------------------------------------------------
  readr::write_lines(x = random.seed, path = file.path(radiator.folder, "random.seed"))
  if (verbose) message("File written: random.seed (", random.seed,")")

  # VCF: Read ------------------------------------------------------------------
  timing.vcf <- proc.time()

  # Get file size
  big.vcf <- file.size(data)
  if (big.vcf <= 500000000) {
    message("\nReading VCF")
  }

  if (big.vcf > 500000000 &&  big.vcf < 900000000) {
    message("\nReading VCF...you might have time for a espresso...")
  }

  if (big.vcf >= 900000000) {
    message("\nReading a large VCF...you actually have time for coffee or tea!\n")
  }
  # Check for bad header generated by stacks
  detect.source <- check_header_source_vcf(data)
  stacks.2 <- detect.source$data.source
  data.source <- stringi::stri_replace_all_regex(
    str = detect.source$check.header$header$value
    [detect.source$check.header$header$id == "source"],
    pattern = "[\"]", replacement = "", vectorize_all = FALSE)
  if (length(data.source) == 0) data.source <- "unknown"
  check.header <- detect.source$check.header
  dp <- "DP" %in% detect.source$check.header$format$ID # Check that DP is trere

  if (!is.null(detect.source$markers.info)) markers.info <- detect.source$markers.info
  if (!is.null(detect.source$overwrite.metadata)) overwrite.metadata <- detect.source$overwrite.metadata

  # for small VCF SeqArray not that good with parallel...
  parallel.temp <- parallel.core
  if (big.vcf < 10000000) {
    parallel.temp <- FALSE
    # parallel_core_opt(parallel.core = parallel.core, max.core = parallel::detectCores())
  }
  big.vcf <- NULL
  gds <- suppressWarnings(
    SeqArray::seqVCF2GDS(
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
  )

  # VCF: Summary ----------------------------------------------------------------
  summary_gds(gds, verbose = TRUE)
  message("done! timing: ", round((proc.time() - timing.vcf)[[3]]), " sec\n")
  if (verbose) message("\nFile written: ", filename.short)
  parallel.temp <- check.header <- detect.source <- NULL #no longer used

  # Generate radiator skeleton -------------------------------------------------
  if (verbose)  message("\nAnalyzing the data...")
  radiator.gds <- radiator_gds_skeleton(gds)

  # VCF: source ----------------------------------------------------------------
  if (length(data.source) > 0) {
    update_radiator_gds(gds = gds, node.name = "data.source", value = data.source)
    if (verbose) message("VCF source: ", data.source)
  } else {
    data.source <- "unknown"
  }

  # VCF: bi- or multi-alllelic--------------------------------------------------
  biallelic <- detect_biallelic_markers(data = gds, verbose = verbose)

  # stacks haplotype VCF...
  if (!biallelic && stacks.2) dp <- FALSE


  # VCF clean sample id---------------------------------------------------------
  individuals.vcf <- tibble::tibble(
    INDIVIDUALS_VCF = SeqArray::seqGetData(gds, "sample.id")) %>%
    dplyr::mutate(INDIVIDUALS_CLEAN = radiator::clean_ind_names(INDIVIDUALS_VCF))

  if (!identical(individuals.vcf$INDIVIDUALS_VCF, individuals.vcf$INDIVIDUALS_CLEAN)) {
    if (verbose) message("Cleaning VCF's sample names")
    clean.id.filename <- stringi::stri_join("cleaned.vcf.id.info_", file.date, ".tsv")
    readr::write_tsv(x = individuals.vcf,
                     path = file.path(radiator.folder, clean.id.filename))

    update_radiator_gds(gds = gds, node.name = "id.clean", value = individuals.vcf)
  }
  # replace id in VCF
  update_radiator_gds(
    gds = gds,
    radiator.gds = FALSE,
    node.name = "sample.id",
    value = individuals.vcf$INDIVIDUALS_CLEAN,
    replace = TRUE)

  individuals <- dplyr::select(individuals.vcf, INDIVIDUALS = INDIVIDUALS_CLEAN)
  individuals.vcf <- NULL

  # VCF sync id with STRATA------------------------------------------------------
  strata <- radiator::read_strata(
    strata = strata,
    pop.levels = pop.levels,
    pop.labels = pop.labels,
    pop.select = pop.select,
    blacklist.id = blacklist.id,
    keep.two = FALSE) %$% strata

  if (!is.null(strata)) {
    id.levels <- individuals$INDIVIDUALS
    individuals %<>%
      dplyr::left_join(
        join_strata(individuals, strata, verbose = verbose) %>%
          dplyr::mutate(FILTERS = "whitelist")
        , by = "INDIVIDUALS"
      ) %>%
      dplyr::mutate(FILTERS = tidyr::replace_na(data = FILTERS, replace = "filter.stata"))

    bl <- dplyr::filter(individuals, FILTERS != "whitelist")
    if (nrow(bl) != 0) {
      if (verbose) message("    Number of sample blacklisted by the stata: ", nrow(bl))
    }

    # renaming ?
    if (rlang::has_name(individuals, "NEW_ID")) {
      if (!identical(id.levels, individuals$INDIVIDUALS)) {
        rlang::abort("Wrong id order, contact author")
      }
      # replace id in GDS
      update_radiator_gds(
        gds = gds,
        radiator.gds = FALSE,
        node.name = "sample.id",
        value = individuals$NEW_ID,
        replace = TRUE)
      individuals %<>% dplyr::rename(INDIVIDUALS = NEW_ID, OLD_ID = INDIVIDUALS)
    }
  } else {
    individuals %<>% dplyr::mutate(STRATA = "1pop", FILTERS = "whitelist")
  }

  strata <- generate_strata(data = dplyr::filter(individuals, FILTERS == "whitelist"), pop.id = FALSE)
  #Update GDS node
  update_radiator_gds(gds = gds, node.name = "individuals.meta", value = individuals, sync = TRUE)
  # summary_gds(gds, check.sync = TRUE, verbose = TRUE)

  # VCF: Markers metadata  ------------------------------------------------------
  markers.meta <- extract_markers_metadata(gds = gds)

  # VCF: reference genome or de novo -------------------------------------------
  ref.genome <- detect_ref_genome(data = gds, verbose = verbose)

  # Stacks specific adjustments
  if (!ref.genome) {
    if (stacks.2) {
      markers.meta %<>%
        dplyr::mutate(
          LOCUS = CHROM,
          CHROM = "1",
          COL = as.integer(POS) - 1
        )
    } else {
      markers.meta %<>% dplyr::mutate(CHROM = "1")
    }
  }

  # ipyrad VCF
  if (stringi::stri_detect_fixed(str = data.source, pattern = "ipyrad")) {
    markers.meta %<>% dplyr::mutate(
      LOCUS = stringi::stri_replace_all_fixed(
        str = CHROM,
        pattern = "locus_",
        replacement = "",
        vectorize_all = FALSE
      ),
      COL = as.integer(POS))
  }

  # PLINK vcf
  if (stringi::stri_detect_fixed(str = data.source, pattern = "PLINK")) {
    markers.meta %<>% dplyr::mutate(
      LOCUS = as.integer(factor(x = CHROM)),
      COL = 1L)
  }

  # GATK, platypus and freebayes specific adjustment
  # Locus with NA or . or ""
  weird.locus <- length(unique(markers.meta$LOCUS)) <= 1
  if (weird.locus && !stacks.2) {
    if (verbose) message("LOCUS field empty... adding unique id instead")
    markers.meta$LOCUS <- markers.meta$VARIANT_ID
  }

  # VCF: LOCUS cleaning and Strands detection ----------------------------------
  if (isTRUE(unique(stringi::stri_detect_fixed(
    str = markers.meta[1,3],
    pattern = c("|", "-", ":", ">")
  )))
  ) {
    data.source <- "dart.vcf"
    update_radiator_gds(gds = gds, node.name = "data.source", value = data.source)
    if (verbose) message("VCF source: ", data.source)
  }

  if (data.source != "dart.vcf" && stringi::stri_detect_regex(str = markers.meta[1,3], pattern = "[^[:alnum:]]+")) {
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
  } else {
    detect.strand <- FALSE
  }
  # VCF: DART LOCUS cleaning----------------------------------------------------
  if (data.source == "dart.vcf") {
    markers.meta %<>%
      dplyr::mutate(
        CHROM = stringi::stri_replace_all_fixed(
          str = CHROM, pattern = ".", replacement = "denovo", vectorize_all = FALSE),
        POS = stringi::stri_replace_na(str = POS, replacement = "50"),
        COL = stringi::stri_extract_first_regex(str = LOCUS, pattern = "[-][0-9]+[\\:]"),
        COL = stringi::stri_replace_all_fixed(str = COL, pattern = c("-", ":"), replacement = c("", ""), vectorize_all = FALSE),
        COL = as.integer(COL),
        LOCUS = stringi::stri_extract_first_regex(str = LOCUS, pattern = "^[0-9]+"),
        LOCUS = stringi::stri_join(LOCUS, POS, sep = "_"),
        POS = COL
      )
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

  # PLINK file with duplicate markers... sometimes tagged ISOFORMS...
  dup.markers <- length(markers.meta$MARKERS) - length(unique(markers.meta$MARKERS))
  if (dup.markers > 0) {
    message("\nNumber of duplicate MARKERS id: ", dup.markers)
    message("Adding integer to differentiate...")
    markers.meta %<>%
      dplyr::arrange(MARKERS) %>%
      dplyr::mutate(MARKERS_NEW = MARKERS) %>%
      dplyr::group_by(MARKERS_NEW) %>%
      dplyr::mutate(
        MARKERS = stringi::stri_join(MARKERS, seq(1, n(), by = 1), sep = "_")
      ) %>%
      dplyr::ungroup(.) %>%
      dplyr::select(-MARKERS_NEW)
  }

  # ADD MARKERS META to GDS
  update_radiator_gds(gds = gds, node.name = "markers.meta", value = markers.meta)

  # replace chromosome info in GDS
  # Why ? well snp ld e.g. will otherwise be performed by chromosome and with de novo data = by locus...
  update_radiator_gds(
    gds = gds,
    radiator.gds = FALSE,
    node.name = "chromosome",
    value = markers.meta$CHROM,
    replace = TRUE
  )
  # # radiator_parameters: generate --------------------------------------------
  filters.parameters <- radiator_parameters(
    generate = TRUE,
    initiate = FALSE,
    update = FALSE,
    parameter.obj = parameters,
    path.folder = radiator.folder,
    file.date = file.date,
    verbose = verbose,
    internal = internal)

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
    internal = internal,
    verbose = verbose
  )


  ##______________________________________________________________________######
  # Filters --------------------------------------------------------------------
  # Filter duplicated SNPs on different strands---------------------------------
  if (detect.strand) {
    blacklist.strands <- markers.meta %>%
      dplyr::distinct(VARIANT_ID, MARKERS, CHROM, LOCUS, POS) %>%
      dplyr::group_by(CHROM, POS) %>%
      dplyr::mutate(n = n()) %>%
      dplyr::ungroup(.) %>%
      dplyr::filter(n > 1) %>%
      dplyr::select(-n) %>%
      dplyr::distinct(CHROM, POS, .keep_all = TRUE)

    if (nrow(blacklist.strands) > 0) {
      if (verbose) {
        message("\nDetected ", nrow(blacklist.strands)," duplicate SNPs on different strands (+/-)")
        message("    By default radiator prune those SNPs")
        message("    To change this behavior, use argument: filter.strands")
      }
      # if (verbose) message("\nFiltering duplicated markers on different strands")
      if (filter.strands == "keep.both") {
        message("Keeping both duplicated markers")
      }

      if (filter.strands == "best.stats") {
        # bk variant
        variant.bk <- markers.meta$VARIANT_ID
        sync_gds(gds = gds, variant.id = blacklist.strands$VARIANT_ID)
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

        # come back to original variant
        sync_gds(gds = gds, variant.id = variant.bk, verbose = FALSE)
      }# End best stats

      if (filter.strands == "blacklist") {
        blacklist.strands %<>% dplyr::distinct(MARKERS)
      }

      if (filter.strands != "keep.both") {
        # Folders---------------------------------------------------------------------
        path.folder.strands <- generate_folder(
          f = path.folder,
          rad.folder = "filter_duplicated_markers",
          internal = FALSE,
          file.date = file.date,
          verbose = verbose)

        # update the blacklist
        blacklist.strands %<>%
          dplyr::mutate(FILTER = "filter.strands")

        write_rad(
          data = blacklist.strands,
          path = path.folder.strands,
          filename = stringi::stri_join("blacklist.duplicated.markers.strands_",
                                        file.date, ".tsv"),
          tsv = TRUE, internal = FALSE, verbose = verbose)

        # bl.gds <- update_bl_markers(gds = gds, update = blacklist.strands)
        # Update markers.meta
        markers.meta %<>%
          dplyr::mutate(
            FILTERS = dplyr::if_else(
              MARKERS %in% blacklist.strands$MARKERS, "filter.strands", FILTERS
            )
          )
        # check <- markers.meta
        write_rad(
          data = dplyr::filter(markers.meta, FILTERS == "whitelist"),
          path = path.folder.strands,
          filename = stringi::stri_join("whitelist.duplicated.markers.strands_",
                                        file.date, ".tsv"),
          tsv = TRUE, internal = FALSE, verbose = verbose)

        # Update GDS
        update_radiator_gds(
          gds = gds,
          node.name = "markers.meta",
          value = markers.meta,
          sync = TRUE)

        # Filter parameter file: update
        filters.parameters <- radiator_parameters(
          generate = FALSE,
          initiate = FALSE,
          update = TRUE,
          parameter.obj = filters.parameters,
          data = gds,
          filter.name = "filter duplicated markers on different strands",
          param.name = "filter.strands",
          values = filter.strands,
          path.folder = path.folder,
          file.date = file.date,
          internal = internal,
          verbose = verbose)
        # results ------------------------------------------------------------------
        radiator_results_message(
          rad.message = stringi::stri_join("\nFilter duplicated markers on different strands: ",
                                           filter.strands),
          filters.parameters,
          internal,
          verbose
        )
      }
    }
  }#Here
  blacklist.strands <- path.folder.strands <- NULL

  # Filter the vcf's FILTER column ---------------------------------------------
  # markers.meta$FILTER_VCF
  filter.vcf <- tibble::tibble(
    FILTER_VCF = SeqArray::seqGetData(
      gds, "annotation/filter")
  ) %>%
    dplyr::bind_cols(
      markers.meta %>%
        dplyr::filter(FILTERS == "whitelist") %>%
        dplyr::select(VARIANT_ID)
    )
  filter.check.unique <- unique(filter.vcf$FILTER_VCF)

  if (length(filter.check.unique) > 1) {
    message("Filtering markers based on VCF FILTER column")
    filter.vcf %<>% dplyr::filter(FILTER_VCF != "PASS")
    markers.meta %<>%
      dplyr::mutate(
        FILTERS = dplyr::if_else(
          VARIANT_ID %in% filter.vcf$VARIANT_ID,
          "vcf filter column",
          FILTERS
        )
      )

    update_radiator_gds(
      gds = gds,
      node.name = "markers.meta",
      value = markers.meta,
      sync = TRUE)

    # radiator_parameters: update
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

    # results ------------------------------------------------------------------
    radiator_results_message(
      rad.message = stringi::stri_join("\nFilter vcf filter column: PASS or not"),
      filters.parameters,
      internal,
      verbose
    )
  }

  filter.check.unique <- NULL
  markers.meta <- individuals <- NULL

  # Filter_monomorphic----------------------------------------------------------
  gds <- filter_monomorphic(
    data = gds,
    filter.monomorphic = filter.monomorphic,
    parallel.core = parallel.core,
    verbose = verbose,
    parameters = filters.parameters,
    path.folder = wf,
    internal = FALSE)

  # Filter common markers ------------------------------------------------------
  gds <- filter_common_markers(
    data = gds,
    filter.common.markers = filter.common.markers,
    fig = TRUE,
    parallel.core = parallel.core,
    verbose = verbose,
    parameters = filters.parameters,
    path.folder = wf,
    internal = FALSE)

  # Filter Whitelist------------------------------------------------------------
  gds <- filter_whitelist(
    data = gds,
    whitelist.markers = whitelist.markers,
    verbose = verbose,
    path.folder = wf,
    parameters = filters.parameters,
    biallelic = biallelic,
    internal = FALSE)

  # Filter_individuals----------------------------------------------------------
  gds <- filter_individuals(
    data = gds,
    interactive.filter = FALSE,
    filter.individuals.missing = filter.individuals.missing,
    filter.individuals.heterozygosity = NULL,
    filter.individuals.coverage.total = filter.individuals.coverage.total,
    parallel.core = parallel.core,
    verbose = verbose,
    path.folder = wf,
    parameters = filters.parameters,
    internal = FALSE
  )

  # Filter MAC----------------------------------------------------------------
  gds <- filter_mac(
    data = gds,
    interactive.filter = FALSE,
    filter.mac = filter.mac,
    filename = NULL,
    parallel.core = parallel.core,
    verbose = verbose,
    parameters = filters.parameters,
    path.folder = wf,
    internal = FALSE)

  # Filter mean coverage------------------------------------------------------
  gds <- filter_coverage(
    data = gds,
    interactive.filter = FALSE,
    filter.coverage = filter.coverage,
    filename = NULL,
    parallel.core = parallel.core,
    verbose = verbose,
    parameters = filters.parameters,
    path.folder = wf,
    internal = FALSE)

  # Filter genotyping---------------------------------------------------------
  gds <- filter_genotyping(
    data = gds,
    interactive.filter = FALSE,
    filter.genotyping = filter.genotyping,
    filename = NULL,
    parallel.core = parallel.core,
    verbose = verbose,
    parameters = filters.parameters,
    path.folder = wf,
    internal = FALSE)

  # Filter SNP position on the read---------------------------------------------
  gds <- filter_snp_position_read(
    data = gds,
    interactive.filter = FALSE,
    filter.snp.position.read = filter.snp.position.read,
    filename = NULL,
    parallel.core = parallel.core,
    verbose = verbose,
    parameters = filters.parameters,
    path.folder = wf,
    internal = FALSE)

  # Filter number of SNPs per locus --------------------------------------------
  gds <- filter_snp_number(
    data = gds,
    interactive.filter = FALSE,
    filter.snp.number = filter.snp.number,
    filename = NULL,
    parallel.core = parallel.core,
    verbose = verbose,
    parameters = filters.parameters,
    path.folder = wf,
    internal = FALSE)

  # Filter Linkage disequilibrium --------------------------------------------
  gds <- filter_ld(
    data = gds,
    interactive.filter = FALSE,
    filter.short.ld = filter.short.ld,
    filter.long.ld = filter.long.ld,
    parallel.core = parallel.core,
    filename = NULL,
    verbose = verbose,
    long.ld.missing = long.ld.missing,
    ld.method = ld.method,
    parameters = filters.parameters,
    path.folder = wf,
    internal = FALSE)

  # Final Sync GDS -----------------------------------------------------------
  if (verbose) message("\nPreparing output files...")
  markers.meta <- extract_markers_metadata(gds, whitelist = TRUE)
  strata <- extract_individuals_metadata(
    gds = gds,
    ind.field.select = c("INDIVIDUALS", "STRATA"),
    whitelist = TRUE)
  sync_gds(gds = gds)
  # summary_gds(gds)
  # generate a folder to put the stats...
  path.folder <- generate_folder(
    f = wf,
    rad.folder = "filtered",
    internal = FALSE,
    file.date = file.date,
    verbose = verbose)

  # Whitelist
  write_rad(data = markers.meta,
            path = path.folder,
            filename = "whitelist.markers.tsv",
            tsv = TRUE,
            write.message = "standard",
            verbose = verbose)
  # blacklist
  bl <- extract_markers_metadata(gds, blacklist = TRUE)
  if (nrow(bl) > 0) {
    write_rad(data = bl,
              path = path.folder,
              filename = "blacklist.markers.tsv",
              tsv = TRUE,
              write.message = "standard",
              verbose = verbose)
  }


  # writing the blacklist of id
  blacklist.id <- extract_individuals_metadata(gds, blacklist = TRUE)
  if (nrow(blacklist.id) > 0) {
    write_rad(data = blacklist.id,
              path = path.folder,
              filename = "blacklist.id.tsv",
              tsv = TRUE,
              write.message = "standard",
              verbose = verbose)
  }

  # Generate new strata
  write_rad(data = strata,
            path = path.folder,
            filename = "strata.filtered.tsv",
            tsv = TRUE,
            write.message = "Writing the filtered strata: strata.filtered.tsv",
            verbose = verbose)

  #______________________________###############################################
  #-----------------------------------  VCF STATS   ----------------------------
  if (vcf.stats) {
    if (verbose) message("\nGenerating statistics after filtering")
    # SUBSAMPLE markers
    n.markers <- length(markers.meta$VARIANT_ID)
    if (n.markers < 200000) subsample.markers.stats <- 1

    if (subsample.markers.stats < 1) {
      markers.subsampled <- dplyr::sample_frac(
        tbl = markers.meta,
        size = subsample.markers.stats)
      variant.select <- markers.subsampled$VARIANT_ID
      subsample.filename <- stringi::stri_join("markers.subsampled_", file.date, ".tsv")
      dplyr::select(markers.subsampled, MARKERS) %>%
        dplyr::mutate(RANDOM_SEED = random.seed) %>%
        readr::write_tsv(
          x = .,
          path = file.path(path.folder, subsample.filename))
      markers.subsampled <- NULL
    } else {
      variant.select <- NULL
    }

    # Individuals stats
    message("Generating individual stats...")
    id.stats <- generate_id_stats(
      gds = gds,
      subsample = variant.select,
      path.folder = path.folder,
      file.date = file.date,
      parallel.core = parallel.core,
      verbose = verbose) %$% info

    # Markers stats
    message("Generating markers stats...")
    markers.stats <- generate_markers_stats(
      gds = gds,
      path.folder = path.folder,
      filename = markers.file,
      file.date = file.date,
      parallel.core = parallel.core,
      verbose = verbose,
      subsample = variant.select) %$% info

    # For summary at the end of the function:
    ind.missing <- round(mean(id.stats$MISSING_PROP, na.rm = TRUE), 2)
    if (dp) ind.cov.total <- round(mean(id.stats$COVERAGE_TOTAL, na.rm = TRUE), 0)
    if (dp) ind.cov.mean <- round(mean(id.stats$COVERAGE_MEAN, na.rm = TRUE), 0)

    # if (!is.null(stats.sub)) {
    markers.missing <- round(mean(markers.stats$MISSING_PROP, na.rm = TRUE), 2)
    if (dp) markers.cov <- round(mean(markers.stats$COVERAGE_MEAN, na.rm = TRUE), 0) # same as above because NA...
    # } else {
    # markers.missing <- round(mean(markers.stats$MISSING_PROP, na.rm = TRUE), 2)
    # if (dp) markers.cov <- round(mean(markers.stats$COVERAGE_MEAN, na.rm = TRUE), 0) # same as above because NA...
    # }
  }#End stats

  #RESULTS --------------------------------------------------------------------
  number.info <- SeqArray::seqGetFilter(gdsfile = gds)
  if (verbose) {
    cat("################################### SUMMARY ####################################\n")
    if (vcf.stats) {
      message("\n\nSummary (AFTER filtering):\nMissing data: ")
      message("    markers: ", markers.missing)
      message("    individuals: ", ind.missing)
      if (dp) message("\n\nCoverage info:")
      if (dp) message("    individuals mean read depth: ", ind.cov.total)
      if (dp) message("    individuals mean genotype coverage: ", ind.cov.mean)
      if (dp) message("    markers mean coverage: ", markers.cov)
    }
  }
  message("\n\nNumber of chromosome/contig/scaffold: ", length(unique(markers.meta$CHROM)))
  message("Number of locus: ", length(unique(markers.meta$LOCUS)))
  message("Number of markers: ", length(number.info$variant.sel[number.info$variant.sel]))
  summary_strata(strata)
  id.stats <- markers.stats <- markers.meta <- NULL

  message("radiator Genomic Data Structure (GDS) file: ", folder_short(filename))
  return(gds)
} # End read_vcf

# tidy_vcf ---------------------------------------------------------------------
#' @name tidy_vcf

#' @title Tidy vcf file
#' @description The function allows to tidy a VCF file.
#'
#' Used internally in
#' \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.
#'
#' It is highly recommended to use \code{\link[radiator]{filter_rad}} to reduce
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
#' \strong{PLINK:} radiator fills the \code{LOCUS} column of PLINK VCFs with
#' a unique integer based on the \code{CHROM} column
#' (\code{as.integer(factor(x = CHROM))}).
#' The \code{COL} column is filled with 1L for lack of bettern info on this.
#' Not what you need ? Open an issue on GitHub for a request.
#'
#' \strong{ipyrad:} the pattern \code{locus_} in the \code{CHROM} column
#' is removed and used. The \code{COL} column is filled with the same value as
#' \code{POS}.
#'
#' \strong{GATK:} Some VCF have an \code{ID} column filled with \code{.},
#' the LOCUS information is all contained along the linkage group in the
#' \code{CHROM} column. To make it work with
#' \href{https://github.com/thierrygosselin/radiator}{radiator},
#' the \code{ID} column is filled with the \code{POS} column info.
#'
#' \strong{platypus:} Some VCF files don't have an ID filed with values,
#' here the same thing is done as GATK VCF files above.
#'
#' \strong{freebayes:} Some VCF files don't have an ID filed with values,
#' here the same thing is done as GATK VCF files above.
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
#' the creation of the GDS file (equivalent of running \code{\link{read_vcf}}).
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
#' \code{"DP", "AD", "GL", "PL", "GQ", "HQ", "GOF", "NR", "NV", "CATG"}.
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
#'  documented in \code{\link[radiator]{read_vcf}}.
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

#' @seealso \code{\link[radiator]{read_vcf}},
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
  # internal = FALSE
  # path.folder = NULL
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
  # random.seed = NULL
  # parameters = NULL

  # Cleanup---------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (verbose) message("Execution date@time: ", file.date)
  old.dir <- getwd()
  opt.change <- getOption("width")
  options(width = 70)
  timing <- radiator_tic()
  res <- list()
  #back to the original directory and options
  on.exit(setwd(old.dir), add = TRUE)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(radiator_toc(timing), add = TRUE)
  on.exit(radiator_function_header(f.name = "tidy_vcf", start = FALSE, verbose = verbose), add = TRUE)

  # Required package -----------------------------------------------------------
  radiator_packages_dep(package = "SeqVarTools", cran = FALSE, bioc = TRUE)

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
    verbose = FALSE
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
  data <- radiator::read_vcf(
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
    markers.meta <- extract_markers_metadata(gds = data, whitelist = TRUE)
    individuals <- extract_individuals_metadata(
      gds = data,
      ind.field.select = "INDIVIDUALS",
      whitelist = TRUE
    )
    n.markers <- nrow(markers.meta)
    n.individuals <- nrow(individuals)

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
      if (ref.calibration) {
        gt.vcf.nuc <- TRUE
      } else {
        if (n.markers < 5000) gt.vcf.nuc <- TRUE
        if (n.markers >= 5000 && n.markers < 30000) gt.vcf.nuc <- FALSE
        if (n.markers >= 30000) gt.vcf.nuc <- FALSE
      }
    }
    # gt is genotype coding a la genepop: 001002, 000000
    if (is.null(gt)) {
      if (n.markers < 5000) gt <- TRUE
      if (n.markers >= 5000 && n.markers < 30000) gt <- FALSE
      if (n.markers >= 30000) gt <- FALSE
    }

    if (gt || gt.vcf.nuc || gt.vcf) ref.calibration <- TRUE

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
      tidy.vcf <- radiator_question(
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

    # Tidying TRUE  ------------------------------------------------------------
    if (tidy.vcf) {
      if (!is.null(blacklist.id)) {
        ref.calibration <- TRUE
        if (verbose) message("\nRe-calibration of REF/ALT alleles: TRUE")
      }

      tidy.data <- gds2tidy(
        gds = data,
        markers.meta = markers.meta,
        calibrate.alleles = FALSE
      )

      # bi- or multi-alllelic VCF ------------------------------------------------
      biallelic <- detect_biallelic_markers(data = data)

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

      # stacks VCF haplotypes
      if (!biallelic) {
        data.source <- radiator::extract_data_source(gds = data)
        if (stringi::stri_detect_fixed(str = data.source, pattern = "Stacks")) {
          overwrite.metadata <- "GT"
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
          want <- c("DP", "AD", "GL", "PL", "HQ", "GQ", "GOF", "NR", "NV", "CATG")

          if (!is.null(overwrite.metadata)) want <- overwrite.metadata

          parse.format.list <- purrr::keep(.x = have, .p = have %in% want)
          if (verbose) message("    genotypes metadata: ", stringi::stri_join(parse.format.list, collapse = ", "))
          # work on parallelization of this part
          tidy.data %<>%
            dplyr::bind_cols(
              purrr::map(
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

      # want <- intersect(c("MARKERS", "CHROM", "LOCUS", "POS", "COL", "REF", "ALT"),
      #                   names(markers.meta))
      # tidy.data <- suppressWarnings(
      #   dplyr::select(markers.meta, dplyr::one_of(want))) %>%
      #   dplyr::bind_cols(
      #     tibble::as_tibble(
      #       matrix(
      #         data = NA,
      #         nrow = n.markers, ncol = n.individuals)) %>%
      #       magrittr::set_colnames(x = ., value = individuals)) %>%
      #   data.table::as.data.table(.) %>%
      #   data.table::melt.data.table(
      #     data = .,
      #     id.vars = want,
      #     variable.name = "INDIVIDUALS",
      #     value.name = "GT",
      #     variable.factor = FALSE) %>%
      #   tibble::as_tibble(.) %>%
      #   dplyr::select(-GT) %>%
      #   dplyr::mutate(
      #     MARKERS = factor(x = MARKERS,
      #                      levels = markers.meta$MARKERS, ordered = TRUE),
      #     INDIVIDUALS = factor(x = INDIVIDUALS,
      #                          levels = individuals,
      #                          ordered = TRUE)) %>%
      #   dplyr::arrange(MARKERS, INDIVIDUALS) %>%
      #   dplyr::bind_cols(tidy.data)
      # Check stacks AD problem --------------------------------------------------
      # Some genotypes with missing AD...

      # NOTE TO MYSELF: might be faster to screen stacks here in data.source...
      if (vcf.metadata) {

        catg.depth <- rlang::has_name(tidy.data, "C_DEPTH")
        if (catg.depth) {
          tidy.data %<>%
            dplyr::mutate(
              ALLELE_REF_DEPTH = dplyr::case_when(
                REF == "C" ~ C_DEPTH,
                REF == "A" ~ A_DEPTH,
                REF == "T" ~ T_DEPTH,
                REF == "G" ~ G_DEPTH
              ),
              ALLELE_ALT_DEPTH = dplyr::case_when(
                ALT == "C" ~ C_DEPTH,
                ALT == "A" ~ A_DEPTH,
                ALT == "T" ~ T_DEPTH,
                ALT == "G" ~ G_DEPTH
              )
            )
        }#catg.depth

        if (rlang::has_name(tidy.data, "ALLELE_REF_DEPTH") &&
            rlang::has_name(tidy.data, "ALLELE_ALT_DEPTH") &&
            rlang::has_name(tidy.data, "READ_DEPTH")
        ) {
          stacks.ad.problem <- tidy.data %>%
            dplyr::select(READ_DEPTH, ALLELE_REF_DEPTH, ALLELE_ALT_DEPTH) %>%
            dplyr::filter(!is.na(READ_DEPTH)) %>%
            dplyr::filter(is.na(ALLELE_REF_DEPTH) & is.na(ALLELE_ALT_DEPTH))
          stacks.ad.problem.n <- nrow(stacks.ad.problem)

          if (stacks.ad.problem.n > 0) {
            non.missing.gt <- length(tidy.data$GT_BIN[!is.na(tidy.data$GT_BIN)])
            stacks.ad.problem.prop <- round(stacks.ad.problem.n / non.missing.gt, 4)

            message("\n\nStacks problem detected")
            message("    missing allele depth info")
            message("    number of genotypes with problem: ", stacks.ad.problem.n, " (prop: ", stacks.ad.problem.prop,")")
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
      strata <- extract_individuals_metadata(
        gds = data,
        ind.field.select = c("STRATA", "INDIVIDUALS"),
        whitelist = TRUE)

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
                "REF", "ALT", "GT_VCF", "GT_VCF_NUC", "GT", "GT_BIN")

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
      return(tidy.data)
    } #tidy.vcf
  } else {
    return(data)
  }#tidy.vcf
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
  parallel.core = parallel::detectCores() - 1
) {

  # x <- parse.format.list
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
  res <- list()
  format.name <- x
  if (verbose) message("\nParsing and tidying: ", format.name)

  # Allele Depth
  if (format.name == "AD") {
    # NOTE TO MYSELF: THIS SHOULD BE DONE FOR BIALLELIC DATA ONLY...
    if (verbose) message("AD column: splitting into ALLELE_REF_DEPTH and ALLELE_ALT_DEPTH")
    # PLAN A:
    res$AD <- SeqArray::seqGetData(
      gdsfile = gds,
      var.name = "annotation/format/AD"
    ) $data %>%
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
    split.vec <- column.vec <- NULL
    # test1 <- res$AD
    #PLAN B:
    # SeqArray::seqApply(
    #   gdsfile = gds,
    #   var.name = "annotation/format/AD",
    #   FUN = function(x) {
    #     ad.sum <- colSums(x, na.rm = TRUE)
    #     ref.col <- which(ad.sum == max(ad.sum, na.rm = TRUE))
    #     ref <- mean(x[, ref.col], na.rm = TRUE)
    #     alt <- mean(x[, -ref.col], na.rm = TRUE)
    #     # return(list(ad.ref = ref, ad.alt = alt))
    #     return(tibble::tibble(ALLELE_REF_DEPTH = ref, ALLELE_ALT_DEPTH = alt))
    #   },
    #   margin = "by.variant", as.is = "list",
    #   parallel = 12) %>%
    #   dplyr::bind_rows(.)
  } # End AD


  # Read depth
  if (format.name == "DP") {
    if (verbose) message("DP column: cleaning and renaming to READ_DEPTH")
    res$DP <- tibble::tibble(READ_DEPTH = SeqArray::seqGetData(
      gdsfile = gds,
      var.name = "annotation/format/DP") %>% as.vector(.))

    # as of SeqArray version 1.28.1 the $data is no longer necessary...
    # res$DP <- tibble::tibble(READ_DEPTH = SeqArray::seqGetData(
    #   gdsfile = gds,
    #   var.name = "annotation/format/DP")$data %>% as.vector(.))

    # test <- res$DP
    # depth <- SeqArray::seqGetData(gds, "annotation/format/DP")
    #
    # if (ind) {
    #   coverage.info$ind.cov.tot <- as.integer(round(rowSums(x = depth$data, na.rm = TRUE, dims = 1L), 0))
    #   coverage.info$ind.cov.mean <- as.integer(round(rowMeans(x = depth$data, na.rm = TRUE, dims = 1L), 0))
    # }
    # if (markers) {
    #   coverage.info$markers.mean <- as.integer(round(colMeans(x = depth$data, na.rm = TRUE, dims = 1L), 0))
    #   coverage.info$markers.tot <- as.integer(round(colSums(x = depth$data, na.rm = TRUE, dims = 1L), 0))
    # }
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
      var.name = "annotation/format/GQ") %>% as.vector(.))
    # as of SeqArray version 1.28.0 breaking change
    # res$GQ <- tibble::tibble(GQ = SeqArray::seqGetData(
    #   gdsfile = gds,
    #   var.name = "annotation/format/GQ")$data %>% as.vector(.))
    # test <- res$GQ
  } # End GQ

  # GL cleaning
  if (format.name == "GL") {
    if (verbose) message("GL column: cleaning Genotype Likelihood column")
    gl <- unique(SeqArray::seqGetData(gdsfile = gds,
                                      var.name = "annotation/format/GL")$length)
    if (gl > 0) {
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
      column.vec <- NULL
    }
  } # End GL

  # Cleaning GOF: Goodness of fit value
  if (format.name == "GOF") {
    if (verbose) message("GOF column: Goodness of fit value")
    res$GOF <- tibble::tibble(GOF = SeqArray::seqGetData(
      gdsfile = gds,
      var.name = "annotation/format/GOF") %>% as.vector(.))
    # as of SeqArray version 1.28.0 breaking change
    # res$GOF <- tibble::tibble(GOF = SeqArray::seqGetData(
    #   gdsfile = gds,
    #   var.name = "annotation/format/GOF")$data %>% as.vector(.))


    # test <- res$GOF
  } # End GOF

  # Cleaning NR: Number of reads covering variant location in this sample
  if (format.name == "NR") {
    if (verbose) message("NR column: splitting column into the number of variant")
    res$NR <- tibble::tibble(NR = SeqArray::seqGetData(
      gdsfile = gds,
      var.name = "annotation/format/NR") %>% as.vector(.))
    # as of SeqArray version 1.28.0 breaking change
    # res$NR <- tibble::tibble(NR = SeqArray::seqGetData(
    #   gdsfile = gds,
    #   var.name = "annotation/format/NR")$data %>% as.vector(.))
    # test <- res$NR
  }#End cleaning NR column

  # Cleaning NV: Number of reads containing variant in this sample
  if (format.name == "NV") {
    if (verbose) message("NV column: splitting column into the number of variant")
    res$NR <- tibble::tibble(NV = SeqArray::seqGetData(
      gdsfile = gds,
      var.name = "annotation/format/NV") %>% as.vector(.))
    # as of SeqArray version 1.28.0 breaking change
    # res$NR <- tibble::tibble(NV = SeqArray::seqGetData(
    #   gdsfile = gds,
    #   var.name = "annotation/format/NV")$data %>% as.vector(.))
    # test <- res$NV
  }#End cleaning NV column


  # CATG from ipyrad...
  if (format.name == "CATG") {
    if (verbose) message("CATG columns: cleaning and renaming C_DEPTH, A_DEPTH, T_DEPTH, G_DEPTH")
    res$CATG <- tibble::as_tibble(
      SeqArray::seqGetData(gdsfile = gds,var.name = "annotation/format/CATG")$data
    )

    column.vec <- seq_along( res$CATG)

    res$CATG <- tibble::tibble(
      C_DEPTH =  res$CATG[, column.vec %% 4 == 1] %>%
        as.matrix(.) %>%
        as.integer(.),
      A_DEPTH =  res$CATG[, column.vec %% 4 == 2] %>%
        as.matrix(.) %>%
        as.integer(.),
      T_DEPTH =  res$CATG[, column.vec %% 4 == 3] %>%
        as.matrix(.) %>%
        as.integer(.),
      G_DEPTH =  res$CATG[, column.vec %% 4 == 0] %>%
        as.matrix(.) %>%
        as.integer(.)
    )
    column.vec <- NULL
  } # End CATG



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
  #     tibble::as_tibble(.) %>%
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
  # future::plan(multiprocess)
  x <- split(x = x, f = split.vec) %>%
    # furrr::future_map(.x = ., .f = clean, .progress = TRUE) %>%
    .radiator_parallel_mc(
      X = ., FUN = clean,
      mc.cores = parallel.core
    ) %>%
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
      .radiator_parallel_mc(
        X = .,
        FUN = clean,
        mc.cores = parallel.core
        # max.vector.size = 1000000000000
      ) %>%
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
          .radiator_parallel_mc(
            X = .,
            FUN = clean,
            mc.cores = parallel.core
            # max.vector.size = 1000000000000
          ) %>%
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
        X = .,
        FUN = clean,
        mc.cores = parallel.core,
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
        nv.col.names = nv.col.names
      ) %>%
      dplyr::bind_rows(.)
  )
  return(x)
}#End clean_nv


# write_vcf---------------------------------------------------------------------
# write a vcf file from a tidy data frame

#' @name write_vcf
#' @title Write a vcf file from a tidy data frame
#' @description Write a vcf file (file format version 4.3, see details below)
#' from a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.

#' @param source source of vcf
#' @param empty empty generate empty vcf
#' @param pop.info (optional, logical) Should the population information be
#' included in the FORMAT field (along the GT info for each samples ?). To make
#' the VCF population-ready use \code{pop.info = TRUE}. The population information
#' must be included in the \code{POP_ID} column of the tidy dataset.
#' Default: \code{pop.info = FALSE}. Experimental.

#' @param filename (optional) The file name prefix for the vcf file
#' written to the working directory. With default: \code{filename = NULL},
#' the date and time is appended to \code{radiator_vcf_file_}.

#' @details \strong{VCF file format version:}
#'
#' If you need a different file format version than the current one, just change
#' the version inside the newly created VCF, that should do the trick.
#' \href{https://vcftools.github.io/specs.html}{For more
#' information on Variant Call Format specifications}.


#' @export
#' @rdname write_vcf

#' @references Danecek P, Auton A, Abecasis G et al. (2011)
#' The variant call format and VCFtools.
#' Bioinformatics, 27, 2156-2158.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

write_vcf <- function(
  data,
  pop.info = FALSE,
  filename = NULL,
  source = NULL,
  empty = FALSE
) {
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")

  if (empty) {
    # output <- tibble::tibble(
    #   '#CHROM' = c("chrom_1", "chrom_1"),
    #   POS = c(50L, 100L),
    #   ID = c(1000299L, 1000299L),
    #   REF = c("A", "G"),
    #   ALT = c("T", "C"),
    #   QUAL = c(".", "."),
    #   FILTER = c("PASS", "PASS"),
    #   INFO = c("NS=3", "NS=3"),
    #   FORMAT = c("GT:DP:AD:GL", "GT:DP:AD:GL"),
    #   ID1 = c("./.:50:30,20:-0.00,-3.97,-17.8", "1/1:50:30,20:-0.00,-3.97,-17.8"),
    #   ID2 = c("0/1:50:30,20:-0.00,-3.97,-17.8", "0/1:50:30,20:-0.00,-3.97,-17.8"),
    #   ID3 = c("1/1:50:30,20:-0.00,-3.97,-17.8", "0/0:50:30,20:-0.00,-3.97,-17.8")
    # )
    output <- tibble::tibble(
      '#CHROM' = character(0),
      POS = integer(0),
      ID = integer(0),
      REF = character(0),
      ALT = character(0),
      QUAL = character(0),
      FILTER = character(0),
      INFO = character(0),
      FORMAT = character(0)
    )
  } else {
    # Import data ---------------------------------------------------------------
    if (is.vector(data)) {
      data <- radiator::tidy_wide(data = data, import.metadata = TRUE)
    }

    # REF/ALT Alleles and VCF genotype format ------------------------------------
    if (!tibble::has_name(data, "GT_VCF")) {
      data <- radiator::calibrate_alleles(data = data)$input
    }

    # Include CHROM, LOCUS, POS --------------------------------------------------
    if (!tibble::has_name(data, "CHROM")) {
      data <- dplyr::mutate(
        .data = data,
        CHROM = rep("1", n()),
        LOCUS = MARKERS,
        POS = MARKERS
      )
    }

    # Order/sort by pop and ind --------------------------------------------------
    if (tibble::has_name(data, "POP_ID")) {
      data <- dplyr::arrange(data, POP_ID, INDIVIDUALS)
    } else {
      data <- dplyr::arrange(data, INDIVIDUALS)
    }

    id.string <- unique(data$INDIVIDUALS)# keep to sort vcf columns
    # Remove the POP_ID column ---------------------------------------------------
    if (tibble::has_name(data, "POP_ID") && (!pop.info)) {
      data <- dplyr::select(.data = data, -POP_ID)
    }

    # Info field -----------------------------------------------------------------
    info.field <- suppressWarnings(
      dplyr::select(.data = data, MARKERS, GT_VCF) %>%
        dplyr::filter(GT_VCF != "./.") %>%
        dplyr::count(x = ., MARKERS) %>%
        dplyr::mutate(INFO = stringi::stri_join("NS=", n, sep = "")) %>%
        dplyr::select(-n)
    )

    # VCF body  ------------------------------------------------------------------
    GT_VCF_POP_ID <- NULL
    if (pop.info) {
      output <- suppressWarnings(
        dplyr::left_join(data, info.field, by = "MARKERS") %>%
          dplyr::select(MARKERS, CHROM, LOCUS, POS, REF, ALT, INFO, INDIVIDUALS, GT_VCF, POP_ID) %>%
          dplyr::mutate(GT_VCF_POP_ID = stringi::stri_join(GT_VCF, POP_ID, sep = ":")) %>%
          dplyr::select(-c(GT_VCF, POP_ID)) %>%
          dplyr::group_by(MARKERS, CHROM, LOCUS, POS, INFO, REF, ALT) %>%
          tidyr::spread(data = ., key = INDIVIDUALS, value = GT_VCF_POP_ID) %>%
          dplyr::ungroup(.) %>%
          dplyr::mutate(
            QUAL = rep(".", n()),
            FILTER = rep("PASS", n()),
            FORMAT = rep("GT:POP", n())
          )
      )

    } else {
      output <- suppressWarnings(
        dplyr::left_join(data, info.field, by = "MARKERS") %>%
          dplyr::select(MARKERS, CHROM, LOCUS, POS, REF, ALT, INDIVIDUALS, GT_VCF, INFO) %>%
          dplyr::group_by(MARKERS, CHROM, LOCUS, POS, INFO, REF, ALT) %>%
          tidyr::spread(data = ., key = INDIVIDUALS, value = GT_VCF) %>%
          dplyr::ungroup(.) %>%
          dplyr::mutate(
            QUAL = rep(".", n()),
            FILTER = rep("PASS", n()),
            FORMAT = rep("GT", n())
          )
      )
    }

    # Transform the REF/ALT format back to A/C/G/T if 001, 002, etc is found
    ref.change <- TRUE %in% unique(c("001", "002", "003", "004") %in% unique(output$REF))

    if (ref.change) {
      output <- output %>%
        dplyr::mutate(
          REF = stringi::stri_replace_all_fixed(
            str = REF,
            pattern = c("001", "002", "003", "004"),
            replacement = c("A", "C", "G", "T"),
            vectorize_all = FALSE),
          ALT = stringi::stri_replace_all_fixed(
            str = ALT,
            pattern = c("001", "002", "003", "004"),
            replacement = c("A", "C", "G", "T"),
            vectorize_all = FALSE)
        )
    }

    if (tibble::has_name(output, "COL")) {
      output <- output %>% dplyr::mutate(LOCUS = stringi::stri_join(LOCUS, COL, sep = "_"))
    } else {
      output <- output %>% dplyr::mutate(LOCUS = stringi::stri_join(LOCUS, as.numeric(POS) - 1, sep = "_"))
    }

    # Keep the required columns
    output <- dplyr::ungroup(output) %>%
      dplyr::arrange(CHROM, LOCUS, POS) %>%
      dplyr::select(-MARKERS) %>%
      dplyr::select('#CHROM' = CHROM, POS, ID = LOCUS, REF, ALT, QUAL, FILTER, INFO, FORMAT, id.string)
    # dplyr::select('#CHROM' = CHROM, POS, ID = LOCUS, REF, ALT, QUAL, FILTER, INFO, FORMAT, dplyr::everything())
  }

  # Filename ------------------------------------------------------------------
  if (is.null(filename)) {
    # Get date and time to have unique filenaming
    filename <- stringi::stri_join("radiator_vcf_file_", file.date, ".vcf")
  } else {
    filename <- stringi::stri_join(filename, ".vcf")
  }

  # File format ----------------------------------------------------------------
  readr::write_delim(
    x = tibble::tibble("##fileformat=VCFv4.3"),
    path = filename, delim = " ", append = FALSE, col_names = FALSE)

  # File date ------------------------------------------------------------------
  x <- paste0("##fileDate=", file.date)
  readr::write_delim(
    x = tibble::tibble(x),
    path = filename,
    delim = " ",
    append = TRUE,
    col_names = FALSE
  )

  # Source ---------------------------------------------------------------------
  if (is.null(source)) {
    source <- stringi::stri_join("##source=radiator_v.",
                                 as.character(utils::packageVersion("radiator")))
    readr::write_delim(
      x = tibble::tibble(source),
      path = filename,
      delim = " ",
      append = TRUE,
      col_names = FALSE)
  } else {
    readr::write_delim(
      x = tibble::tibble(
        stringi::stri_join('##source=', rlang::quo_name(source))
      ),
      path = filename,
      delim = " ",
      append = TRUE,
      col_names = FALSE
    )
  }


  # Info field 1 ---------------------------------------------------------------
  info1 <- as.data.frame('##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">')
  utils::write.table(x = info1, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)


  # Format field 1 -------------------------------------------------------------
  format1 <- '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
  format1 <- as.data.frame(format1)
  utils::write.table(x = format1, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)

  # Format field 2 ---------------------------------------------------------------
  if (pop.info && !empty) {
    format2 <- as.data.frame('##FORMAT=<ID=POP_ID,Number=1,Type=Character,Description="Population identification of Sample">')
    utils::write.table(x = format2, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  }
  # Format field DP ---------------------------------------------------------------
  # if (empty) {
  #   format3 <- as.data.frame('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">')
  #   utils::write.table(x = format3, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  # }

  # Format field AD ---------------------------------------------------------------
  # if (empty) {
  #   format4 <- as.data.frame('##FORMAT=<ID=AD,Number=.,Type=String,Description="Allele Depth">')
  #   utils::write.table(x = format4, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  # }

  # Format field GL ---------------------------------------------------------------
  # if (empty) {
  #   format5 <- as.data.frame('##FORMAT=<ID=GL,Number=.,Type=String,Description="Genotype Likelihood">')
  #   utils::write.table(x = format5, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  # }

  #else {
  #   # Empty VCF-----------------------------------------------------------------
  #   output <- tibble::tibble(
  #     `#CHROM` = character(0),
  #     POS = character(0),
  #     ID = character(0),
  #     REF = character(0),
  #     ALT = character(0),
  #     QUAL = character(0),
  #     FILTER = character(0),
  #     INFO = character(0),
  #     FORMAT = character(0)
  #   )
  # }

  # Write the prunned vcf to the file ------------------------------------------
  suppressWarnings(readr::write_tsv(x = output, path = filename, append = TRUE, col_names = TRUE))

}# end write_vcf


# extract_individuals_vcf-------------------------------------------------------
#' @title Extract individuals from vcf file
#' @description Function that returns the individuals present in a vcf file.
#' Useful to create a strata file or
#' to make sure you have the right individuals in your VCF.
#' @param data (character) The path to the vcf file.
#' @rdname extract_individuals_vcf
#' @export
#' @return A tibble with a column: \code{INDIVIDUALS}.
#' @seealso \pkg{radiator} \code{\link{read_strata}}
#' @author Thierry Gosselin \email{thierrygosselin@icloud.com}
#'
#' @examples
#' \dontrun{
#' # Built a strata file:
#' strata <- radiator::extract_individuals_vcf("my.vcf") %>%
#'     dplyr::mutate(STRATA = "fill this") %>%
#'     readr::write_tsv(x = ., path = "my.new.vcf.strata.tsv")
#' }

extract_individuals_vcf <- function(data) {
  temp.file <-
    suppressWarnings(suppressMessages(
      readr::read_table(file = data, n_max = 200, col_names = "HEADER")
    ))
  skip.number <- which(stringi::stri_detect_fixed(str = temp.file$HEADER,
                                                  pattern = "#CHROM")) - 1
  # for some VCF file that put all the markers (usually contigs info) in the header...
  if (length(skip.number) == 0L) {
    temp.file <-
      suppressWarnings(suppressMessages(
        readr::read_table(file = data, col_names = "HEADER")
      ))
    skip.number <- which(stringi::stri_detect_fixed(str = temp.file$HEADER,
                                                    pattern = "#CHROM")) - 1
  }
  temp.file <- NULL
  remove <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
  id <- tibble::tibble(INDIVIDUALS = colnames(readr::read_tsv(
    file = data,
    n_max = 1,
    skip = skip.number,
    col_types = readr::cols(.default = readr::col_character())) %>%
      dplyr::select(-dplyr::one_of(remove))))
  return(id)
}#End extract_individuals_vcf




# extract_info_vcf-------------------------------------------------------
#' @title extract_info_vcf
#' @description Extract vcf information
#' @rdname extract_info_vcf
#' @keywords internal
#' @export
extract_info_vcf <- function(vcf) {
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
}#End extract_info_vcf


# check_header_source_vcf-------------------------------------------------------
#' @title Check the vcf header and detect vcf source
#' @description Check the vcf header and detect vcf source
#' @rdname check_header_source_vcf
#' @keywords internal
#' @export
check_header_source_vcf <- function(vcf) {

  check.header <- SeqArray::seqVCF_Header(vcf)
  problematic.id <- c("AD", "AO", "QA", "GL", "CATG", "RO", "QR", "MIN_DP")
  problematic.id <-
    purrr::keep(
      .x = problematic.id,
      .p = problematic.id %in% check.header$format$ID
    )
  for (p in problematic.id) {
    check.header$format[check.header$format$ID == p, "Number"] <- "."
  }
  # check.header$format

  # DArT vcf problem
  probl.dart <- check.header$format[check.header$format$ID == "GT", "Type"] == "Integer"
  if (length(probl.dart) == 0) probl.dart <- FALSE
  if (probl.dart) check.header$format[check.header$format$ID == "GT", "Type"] <- "String"

  check.source <- check.header$header$value[check.header$header$id == "source"]
  if (length(check.source) == 0) {
    is.stacks <- FALSE
    check.source <- "unknown"
  } else {
    is.stacks <- stringi::stri_detect_fixed(str = check.source, pattern = "Stacks")
  }

  dirty.freebayes <- FALSE
  if (stringi::stri_detect_fixed(str = check.source, pattern = "freeBayes") &&
      stringi::stri_detect_fixed(str = check.source, pattern = "dirty")) {
    message("\n\nIMPORTANT: VCF from freeBayes Dirty version: only GT and DP fields will be kept...\n\n")
    check.header$format <- dplyr::filter(check.header$format, ID %in% c("GT", "DP"))
  }

  if (is.stacks) {
    stacks.2 <- keep.stacks.gl <- stringi::stri_detect_fixed(
      str = check.source,
      pattern = "Stacks v2")
    keep.stacks.gl <- TRUE
    if (!keep.stacks.gl) {
      check.header$format <- dplyr::filter(check.header$format, ID != "GL")
    }
    markers.info <- NULL
    overwrite.metadata <- NULL
  } else {
    stacks.2 <- FALSE
    markers.info <- NULL
    overwrite.metadata <- NULL
  }
  return(
    res = list(
      data.source = stacks.2,
      check.header = check.header,
      markers.info = markers.info,
      overwrite.metadata = overwrite.metadata
    )
  )
}#End check_header_source_vcf

## Split vcf--------------------------------------------------------------------

#' @title Split a VCF file
#' @description This function allows to split a VCF file in several VCFs,
#' based on individuals or populations.

#' @param strata A file identical to the strata file usually used in radiator,
#' with an additional column named: \code{SPLIT}.
#' This new column contains numerical values
#' (e.g. 1, 1, 1, ..., 2, 2, 2, 2, ..., 3, 3, ...),
#' that indicate for each INDIVIDUALS/STRATA, how to split.
#' The number of VCF to split to is based on the max value found in the column
#' \code{SPLIT}, above this would result in 3 VCF files created).

#' @inheritParams tidy_genomic_data
#' @inheritParams radiator_common_arguments


#' @return The function returns in the global environment a list with
#' the different tidy dataset from the split vcf. In the working directory,
#' the splitted VCF files with \code{"_1", "_2"} in the name.

#' @examples
#' \dontrun{
#' split.data <- radiator::split_vcf(
#' data = "batch_1.vcf",
#' strata = "strata.split.tsv",
#' blacklist.id = "blacklisted.id.txt",
#' whitelist.markers = "whitelist.loci.txt")
#' }


#' @rdname split_vcf
#' @export
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

split_vcf <- function(
  data,
  strata,
  parallel.core = parallel::detectCores() - 1,
  ...
) {
  cat("#######################################################################\n")
  cat("########################### radiator::split_vcf #########################\n")
  cat("#######################################################################\n")
  timing <- proc.time()

  # manage missing arguments -----------------------------------------------------
  if (missing(data)) rlang::abort("data/vcf file missing")
  if (missing(strata)) rlang::abort("strata file missing")

  # if (!is.null(pop.levels) & is.null(pop.labels)) {
  #   pop.levels <- stringi::stri_replace_all_fixed(
  #     pop.levels, pattern = " ", replacement = "_", vectorize_all = FALSE)
  #   pop.labels <- pop.levels
  # }
  # if (!is.null(pop.labels) & is.null(pop.levels)) rlang::abort("pop.levels is required if you use pop.labels")
  # if (!is.null(pop.labels)) {
  #   if (length(pop.labels) != length(pop.levels)) rlang::abort("pop.labels and pop.levels must have the same length (number of groups)")
  #   pop.labels <- stringi::stri_replace_all_fixed(pop.labels, pattern = " ", replacement = "_", vectorize_all = FALSE)
  # }

  # Filename -------------------------------------------------------------------
  # Get date and time to have unique filenaming
  file.date <- stringi::stri_replace_all_fixed(
    Sys.time(),
    pattern = " EDT",
    replacement = "",
    vectorize_all = FALSE
  )
  file.date <- stringi::stri_replace_all_fixed(
    file.date,
    pattern = c("-", " ", ":"),
    replacement = c("", "@", ""),
    vectorize_all = FALSE
  )
  file.date <- stringi::stri_sub(file.date, from = 1, to = 13)

  filename <- stringi::stri_join("radiator_split_vcf_", file.date)

  # import data ----------------------------------------------------------------
  # if (!is.null(blacklist.id)) {# With blacklist of ID
  #   if (is.vector(blacklist.id)) {
  #     suppressMessages(blacklist <- readr::read_tsv(blacklist.id, col_names = TRUE))
  #   } else {
  #     if (!tibble::has_name(blacklist.id, "INDIVIDUALS")) {
  #       rlang::abort("Blacklist of individuals should have 1 column named: INDIVIDUALS")
  #     }
  #     blacklist <- blacklist.id
  #   }
  #   blacklist$INDIVIDUALS <- stringi::stri_replace_all_fixed(
  #     str = blacklist$INDIVIDUALS,
  #     pattern = c("_", ":"),
  #     replacement = c("-", "-"),
  #     vectorize_all = FALSE
  #   )
  #
  #   # remove potential duplicate id
  #   blacklist <- dplyr::distinct(.data = blacklist, INDIVIDUALS)
  # }


  split <- suppressMessages(readr::read_tsv(file = strata))
  strata <- dplyr::select(split, -SPLIT)
  split <- dplyr::select(split, -STRATA)

  # if (!is.null(blacklist.id)) {
  #   split <- dplyr::anti_join(x = split, y = blacklist, by = "INDIVIDUALS")
  # }

  split$INDIVIDUALS <- stringi::stri_replace_all_fixed(
    str = split$INDIVIDUALS,
    pattern = c("_", ":"),
    replacement = c("-", "-"),
    vectorize_all = FALSE
  )

  # Function required
  split_vcf <- function(data, filename) {
    split.id <- unique(data$SPLIT)
    filename <- stringi::stri_join(filename, "_", split.id)
    radiator::write_vcf(data = dplyr::select(data, -SPLIT),
                        pop.info = FALSE, filename = filename)
  }


  message("Importing and tidying the vcf...")
  input <- suppressMessages(
    radiator::tidy_genomic_data(
      data = data,
      strata = strata,
      vcf.metadata = FALSE,
      whitelist.markers = whitelist.markers,
      filename = NULL,
      verbose = FALSE) %>%
      dplyr::full_join(split, by = "INDIVIDUALS") %>%
      split(x = ., f = .$SPLIT)
  )

  split <- strata <- blacklist <- NULL

  .radiator_parallel_mc(
    X = input,
    FUN = split_vcf,
    mc.cores = parallel.core,
    filename = filename#,max.vector.size = 1000000000000
  )


  # results --------------------------------------------------------------------
  message("Split VCFs were written in the working directory")
  timing <- proc.time() - timing
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("############################## completed ##############################\n")
  return(input)
}#End split_vcf


# Merge vcf---------------------------------------------------------------------

# @title Merge VCF files
# @description This function allows to merge 2 VCF files.

# @param vcf1 First VCF file.
# @param strata1 strata file for vcf1.
# @param vcf2 Second VCF file.
# @param strata2 strata file for vcf2.
# @param filename Name of the merged VCF file.
# With the default, the function gives a filename based on date and time.
# Default: \code{filename = NULL}.
# @inheritParams tidy_genomic_data

# @return The function returns in the global environment a tidy dataset with
# the merged VCF files and the merged VCF in the working directory.

# @examples
# \dontrun{
# # The simplest way to run the function:
# sum <- radiator::merge_vcf(
# vcf1 = "batch_1.vcf", strata1 = "strata1_brook_charr.tsv",
# vcf1 = "batch_2.vcf", strata2 = "strata2_brook_charr.tsv",
# pop.select = c("QC", "ON", "NE"),
# maf.thresholds = c(0.002, 0.001),
# maf.pop.num.threshold = 1,
# maf.approach = "SNP",maf.operator = "OR",
# filename = "my_new_VCF.vcf"
# }

# @rdname merge_vcf
# @export
# @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

# merge_vcf <- function(
#   vcf1, strata1,
#   vcf2, strata2,
#   whitelist.markers = NULL,
#   filename = NULL,
#   parallel.core = parallel::detectCores() - 1
# ) {
#   cat("#######################################################################\n")
#   cat("########################### radiator::merge_vcf #########################\n")
#   cat("#######################################################################\n")
#   timing <- proc.time()
#
#   # manage missing arguments
#   if (missing(vcf1)) rlang::abort("vcf1 file missing")
#   if (missing(vcf2)) rlang::abort("vcf2 file missing")
#   if (missing(strata1)) rlang::abort("strata1 file missing")
#   if (missing(strata2)) rlang::abort("strata2 file missing")
#
#
#   # if (!is.null(pop.levels) & is.null(pop.labels)) {
#   #   pop.levels <- stringi::stri_replace_all_fixed(
#   #     pop.levels, pattern = " ", replacement = "_", vectorize_all = FALSE)
#   #   pop.labels <- pop.levels
#   # }
#   # if (!is.null(pop.labels) & is.null(pop.levels)) rlang::abort("pop.levels is required if you use pop.labels")
#   # if (!is.null(pop.labels)) {
#   #   if (length(pop.labels) != length(pop.levels)) rlang::abort("pop.labels and pop.levels must have the same length (number of groups)")
#   #   pop.labels <- stringi::stri_replace_all_fixed(pop.labels, pattern = " ", replacement = "_", vectorize_all = FALSE)
#   # }
#
#   # Filename
#   # Get date and time to have unique filenaming
#   if (is.null(filename)) {
#     file.date <- format(Sys.time(), "%Y%m%d@%H%M")
#   }
#
#   # import data
#   message("Importing and tidying the vcf1...")
#   input <- suppressMessages(radiator::tidy_genomic_data(
#     data = vcf1,
#     strata = strata1,
#     vcf.metadata = FALSE,
#     whitelist.markers = whitelist.markers,
#     filename = NULL,
#     verbose = FALSE))
#
#   message("Importing and tidying the vcf2...")
#   # Also Using pop.levels and pop.labels info if present
#   input <- suppressWarnings(
#     dplyr::bind_rows(
#       input,
#       suppressMessages(
#         radiator::tidy_genomic_data(
#           data = vcf2,
#           strata = strata2,
#           vcf.metadata = FALSE,
#           whitelist.markers = whitelist.markers,
#           filename = NULL,
#           verbose = FALSE))))
#
#   message("Adjusting REF/ALT alleles...")
#   input <- radiator::calibrate_alleles(
#     data = input,
#     parallel.core = parallel.core,
#     verbose = TRUE)$input
#
#   # if (filter.monomorphic) {
#   #   input <- radiator::filter_monomorphic(data = input, verbose = TRUE)
#   # }
#   #
#   # if (filter.common.markers) {
#   #   input <- radiator::filter_common_markers(data = input, verbose = TRUE)$input
#   # }
#   #
#   # if (!is.null(filter.mac)) {
#   #   input <- filter_maf(
#   #     data = input,
#   #     interactive.filter = FALSE,
#   #     filter.mac = filter.mac,
#   #     parallel.core = parallel.core,
#   #     verbose = FALSE)$tidy.filtered.mac
#   # }
#
#   # Write VCF in the working directory
#   radiator::write_vcf(data = input, pop.info = FALSE, filename = filename)
#
#   # results
#   message("Merged VCF in the working directory: ", filename, ".vcf")
#   timing <- proc.time() - timing
#   message("\nComputation time: ", round(timing[[3]]), " sec")
#   cat("############################## completed ##############################\n")
#   return(input)
# }#End merge_vcf



# vcf_strata -------------------------------------------------------------------
#' @name vcf_strata
#' @title Join stratification metadata to a VCF (population-aware VCF)
#' @description Include stratification metadata, e.g. population-level information,
#' to the \code{FORMAT} field of a VCF file.
#' @param data A VCF file

#' @param strata (optional) A tab delimited file at least 2 columns
#' with header:
#' \code{INDIVIDUALS} and \code{STRATA}.
#' The \code{STRATA} and any other columns can be any hierarchical grouping.
#' To create a strata file see \code{\link[radiator]{individuals2strata}}.

#' @param filename (optional) The file name for the modifed VCF,
#' written to the working directory. Default: \code{filename = NULL} will make a
#' custom filename with data and time.

#' @export
#' @rdname vcf_strata

#' @return A VCF file in the working directory with new \code{FORMAT} field(s)
#' correponding to the strata column(s).

#' @seealso
#' \href{https://vcftools.github.io}{VCF web page}
#'
#' \href{VCF specification page}{https://vcftools.github.io/specs.html}

#' @references Danecek P, Auton A, Abecasis G et al. (2011)
#' The variant call format and VCFtools.
#' Bioinformatics, 27, 2156-2158.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @keywords internal


vcf_strata <- function(data, strata, filename = NULL) {
  # data <- "batch_1.vcf"
  # strata <- "strata.sturgeon.12pop.tsv"
  # filename <- NULL
  # data <- "example_vcf2dadi_ferchaud_2015.vcf"
  # strata <- "strata.stickleback.tsv"

  cat("#######################################################################\n")
  cat("######################### radiator: vcf_strata ##########################\n")
  cat("#######################################################################\n")

  # Checking for missing and/or default arguments ******************************
  if (missing(data)) rlang::abort("Input file missing")
  if (missing(strata)) rlang::abort("Strata file missing")

  # import the first 50 lines
  quick.scan <- readr::read_lines(file = data, n_max = 75)

  # Function to detect where CHROM line starts
  detect_vcf_header <- function(x) {
    stringi::stri_detect_fixed(str = x, pattern = "CHROM", negate = FALSE)
  }

  # Detect the index
  max.vcf.header <- purrr::detect_index(.x = quick.scan, .p = detect_vcf_header) - 1

  # import VCF header and add a row containing the new format field
  vcf.header <- readr::read_delim(file = data, n_max = max.vcf.header, col_names = "VCF_HEADER", delim = "\n")

  # import the vcf file, no filters etc.
  message("Importing the VCF file")
  input <- data.table::fread(
    input = data,
    sep = "\t",
    stringsAsFactors = FALSE,
    header = TRUE,
    skip = "CHROM",
    showProgress = TRUE,
    verbose = FALSE
  ) %>%
    tibble::as_tibble()

  # transform in long format
  input <- data.table::melt.data.table(
    data = data.table::as.data.table(input),
    id.vars = c("#CHROM", "POS", "ID",  "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"),
    variable.name = "INDIVIDUALS",
    variable.factor = FALSE,
    value.name = "FORMAT_ID"
  ) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      INDIVIDUALS = stringi::stri_replace_all_fixed(
        str = INDIVIDUALS,
        pattern = c("_", ":"),
        replacement = c("-", "-"),
        vectorize_all = FALSE)
    )

  # population levels and strata  ---------------------------------------------
  message("Importing the strata file")
  if (is.vector(strata)) {
    # message("strata file: yes")
    number.columns.strata <- max(utils::count.fields(strata, sep = "\t"))
    col.types <- stringi::stri_join(rep("c", number.columns.strata), collapse = "")
    strata.df <- readr::read_tsv(file = strata, col_names = TRUE, col_types = col.types) %>%
      dplyr::rename(POP_ID = STRATA)
  } else {
    # message("strata object: yes")
    colnames(strata) <- stringi::stri_replace_all_fixed(
      str = colnames(strata),
      pattern = "STRATA",
      replacement = "POP_ID",
      vectorize_all = FALSE
    )
    strata.df <- strata
  }

  strata.number <- length(strata.df) - 1
  strata.colnames <- purrr::discard(.x = colnames(strata.df), .p = colnames(strata.df) %in% "INDIVIDUALS")

  # Replace unwanted whitespace pattern in the strata
  strata.df <- strata.df %>%
    dplyr::mutate(
      INDIVIDUALS = stringi::stri_replace_all_fixed(
        str = INDIVIDUALS,
        pattern = c("_", ":"),
        replacement = c("-", "-"),
        vectorize_all = FALSE
      ),
      POP_ID = stringi::stri_replace_all_fixed(
        str = POP_ID,
        pattern = " ",
        replacement = "_",
        vectorize_all = FALSE
      )
    )


  # Join strata to input and merge strata column to FORMAT field
  message("Joining the strata to the VCF into new field format...")
  input <- input %>%
    dplyr::left_join(strata.df, by = "INDIVIDUALS") %>%
    tidyr::unite_(
      data = .,
      col = "FORMAT_ID",
      from = c("FORMAT_ID", strata.colnames),
      sep = ":",
      remove = TRUE
    ) %>%
    dplyr::group_by(`#CHROM`, POS, ID,  REF, ALT, QUAL, FILTER, INFO, FORMAT) %>%
    tidyr::spread(data = ., key = INDIVIDUALS, value = FORMAT_ID) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(
      FORMAT = stringi::stri_join(
        FORMAT,
        stringi::stri_join(strata.colnames, collapse = ":"),
        sep = ":",
        collapse = NULL
      )
    )

  # Filename ------------------------------------------------------------------
  message("Writing to the working directory...")
  if (is.null(filename)) {
    # Get date and time to have unique filenaming
    file.date <- stringi::stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "", vectorize_all = FALSE)
    file.date <- stringi::stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
    file.date <- stringi::stri_sub(file.date, from = 1, to = 13)
    filename <- stringi::stri_join("radiator_vcf_file_", file.date, ".vcf")
  } else {
    filename <- stringi::stri_join(filename, ".vcf")
  }
  # File format ----------------------------------------------------------------
  # write_delim(x = data_frame("##fileformat=VCFv4.2"), path = filename, delim = " ", append = FALSE, col_names = FALSE)
  vcf.header[1,] <- "##fileformat=VCFv4.3"

  # File date ------------------------------------------------------------------
  file.date <- stringi::stri_replace_all_fixed(Sys.Date(), pattern = "-", replacement = "")
  file.date <- stringi::stri_join("##fileDate=", file.date, sep = "")
  # write_delim(x = data_frame(file.date), path = filename, delim = " ", append = TRUE, col_names = FALSE)
  vcf.header[2,] <- file.date

  # Source ---------------------------------------------------------------------
  # write_delim(x = data_frame(stringi::stri_join("##source=radiator_v.", utils::packageVersion("radiator"))), path = filename, delim = " ", append = TRUE, col_names = FALSE)
  # vcf.header[3,] <- stringi::stri_replace_all_fixed(str = vcf.header[3,], pattern = '"', replacement = "", vectorize_all = FALSE)
  # vcf.header[3,]<- stringi::stri_join(vcf.header[3,], "and radiator v.", utils::packageVersion("radiator"))

  # New FORMAT -----------------------------------------------------------------
  for (i in strata.colnames) {
    vcf.header <- tibble::add_row(
      .data = vcf.header,
      VCF_HEADER = stringi::stri_join(
        "##FORMAT=<ID=", i, ',Number=1,Type=Character,Description="New strata",Source="radiator",Version="', utils::packageVersion("radiator"), '">')
    )
  }
  # VCF HEADER  ------------------------------------------------------------------
  utils::write.table(x = vcf.header, file = filename, sep = " ", append = FALSE, col.names = FALSE, quote = FALSE, row.names = FALSE)

  # Write the data   -------------------------------------------------------------
  suppressWarnings(readr::write_tsv(x = input, path = filename, append = TRUE, col_names = TRUE))

  cat("############################## completed ##############################\n")
}#vcf_strata
