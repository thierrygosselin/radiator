# Read VCF with SeqArray ------------------------------------------------------
#' @title Read VCF files and write a radiator GDS file
#'
#' @name read_vcf
#' @rdname read_vcf
#'
#' @description
#' Read a VCF file and convert it to a
#' \href{https://github.com/zhengxwen/SeqArray}{SeqArray} GDS object
#' (\code{SeqVarGDSClass}; Zheng et al. 2017), with a
#' \code{/radiator} node containing radiator-specific metadata.
#'
#' The function has an "advanced" mode (via \code{...}) that allows a number of
#' filters to be applied on the fly (MAC, coverage, genotyping rate, LD, etc.),
#' as well as several VCF-specific clean-ups (Stacks, ipyrad, PLINK, DArT,
#' FreeBayes, GATK, samtools, etc.).
#'
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users who want a fast and robust VCF → GDS
#' import with basic QC.
#'
#'
#' @param data (character)
#' Path to a VCF file (optionally bgzipped, \code{*.vcf} or \code{*.vcf.gz}).
#' Markers can be biallelic SNPs or haplotypes; downstream filters
#' (e.g. \code{filter.haplotype.format}) will normalise variants to the
#' desired representation.
#'
#' @param strata (optional)
#' Strata definition, passed to \code{\link[radiator]{read_strata}}.
#' Can be a path to a strata file or an object.
#' Default: \code{strata = NULL}.
#'
#' @param filename (optional, character)
#' Base name of the GDS file to generate. Radiator will append
#' \code{.gds.rad} to the filename. If the chosen filename already exists
#' in \code{path.folder}, a timestamped default name is used instead.
#' Default: \code{filename = NULL}.
#'
#' @param vcf.stats (logical, optional)
#' Generate basic statistics for individuals and markers (missingness, coverage,
#' etc.) and write them to disk.
#' Computational cost can be high for very large unfiltered VCF
#' Default: \code{vcf.stats = FALSE}.
#'
#' @param parallel.core (integer, optional)
#' Approximate number of cores to use where parallelism is supported
#' (e.g. \code{SeqArray::seqApply}, some filter helpers).
#' Default: \code{parallel.core = parallel::detectCores() - 1}.
#'
#' @param verbose (logical, optional)
#' When \code{TRUE}, the function prints progress messages and a summary.
#' Default: \code{verbose = TRUE}.
#'
#' @inheritParams tidy_genomic_data
#' @inheritParams radiator_common_arguments
#'
#' @details
#' Typical performance (rough order of magnitude):
#' \itemize{
#'   \item a 35 GB VCF with ~4M SNPs: \eqn{\sim}7 minutes with 8 CPU;
#'   \item a 21 GB VCF with ~2M SNPs: \eqn{\sim}5 minutes with 7 CPU.
#' }
#'
#' The resulting GDS file can be reopened almost instantly in a later R session
#' with \code{radiator::read_rad()}.
#'
#' \strong{Heterozygosity and inbreeding (Fis):}
#' The heterozygosity statistics generated here are global across strata and
#' primarily descriptive. For filtering on het/Fis/HWE, we recommend using
#' \code{\link{filter_het}}, \code{\link{filter_fis}} and
#' \code{\link{filter_hwe}} \emph{after} obvious outlier individuals and
#' markers have been removed.
#'
#' @section VCF file format behaviour:
#'
#' \strong{PLINK:}
#' \itemize{
#'   \item \code{LOCUS} is filled with an integer based on \code{CHROM}
#'   (\code{as.integer(factor(x = CHROM))});
#'   \item \code{COL} is set to \code{1L} (no within-read position available).
#' }
#'
#' \strong{ipyrad:}
#' \itemize{
#'   \item the pattern \code{"locus_"} is stripped from \code{CHROM};
#'   \item \code{COL} is set to \code{POS}.
#' }
#'
#' \strong{GATK / platypus / FreeBayes / samtools:}
#' \itemize{
#'   \item if the VCF \code{ID} column is \code{.}, it is replaced with
#'   the position (\code{POS});
#'   \item short-read locus identity is assumed to be encoded in
#'   \code{CHROM} + \code{POS}.
#' }
#'
#' \strong{Stacks:}
#' \itemize{
#'   \item \emph{de novo}: \code{CHROM} is typically "1";
#'     \code{LOCUS} corresponds to "CHROM" in the Stacks VCF;
#'     \code{COL} is \code{POS - 1};
#'   \item \emph{reference}: \code{ID} is split into \code{LOCUS}, \code{COL},
#'     \code{STRANDS}.
#' }
#'
#' \strong{DArT VCFs:}
#' \itemize{
#'   \item \code{CHROM == "."} is replaced by \code{"denovo"};
#'   \item missing \code{POS} (\code{NA}) are set to \code{50};
#'   \item \code{COL} is extracted from \code{LOCUS} (read position);
#'   \item \code{LOCUS} is the first group of digits, then joined with \code{POS}
#'     using \code{"_"}, and \code{POS} is replaced by \code{COL}.
#' }
#'
#' @section Advanced mode (\code{...}):
#'
#' The \code{...} lets you pass many additional arguments used by radiator’s
#' filtering framework, for example:
#' \itemize{
#'   \item \code{whitelist.markers}, \code{blacklist.id},
#'     \code{pop.select}, \code{pop.levels}, \code{pop.labels};
#'   \item \code{filter.strands} – handle duplicate SNPs on opposite strands;
#'   \item \code{filter.common.markers} – keep markers common to all strata;
#'   \item \code{filter.ma} – global MAC/MAF/DP filtering;
#'   \item \code{filter.coverage}, \code{filter.genotyping};
#'   \item \code{filter.snp.position.read}, \code{filter.snp.number};
#'   \item \code{filter.short.ld}, \code{filter.long.ld}, \code{ld.method};
#'   \item \code{filter.individuals.missing}, coverage filters by individual;
#'   \item \code{markers.info}, \code{vcf.metadata};
#'   \item \code{path.folder}, \code{random.seed}, \code{subsample.markers.stats};
#'   \item \code{filter.monomorphic}, \code{filter.haplotype.format}, etc.
#' }
#'
#' See the documentation of the individual filter functions
#' (e.g. \code{\link{filter_ma}}, \code{\link{filter_ld}},
#' \code{\link{filter_monomorphic}}) for more detail.
#'
#' @return
#' A \code{SeqVarGDSClass} object (radiator GDS) with a \code{/radiator} node
#' containing per-marker and per-individual metadata and a record of all
#' filters applied.
#'
#' @seealso
#' \code{\link{filter_ld}},
#' \code{\link[radiator]{tidy_genomic_data}},
#' \code{\link[radiator]{tidy_vcf}}
#'
#' @references
#' Zheng X, Gogarten S, Lawrence M, Stilp A, Conomos M, Weir BS, Laurie C,
#' Levine D (2017). SeqArray -- A storage-efficient high-performance data
#' format for WGS variant calls. *Bioinformatics*.
#'
#' Danecek P, Auton A, Abecasis G et al. (2011) The variant call format and
#' VCFtools. *Bioinformatics* 27:2156–2158.
#'
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#'
#' @examples
#' \dontrun{
#' # Simple import, no strata, defaults:
#' gds <- radiator::read_vcf(data = "populations.snps.vcf")
#'
#' # With strata and a few filters:
#' gds <- radiator::read_vcf(
#'   data                      = "populations.snps.vcf",
#'   strata                    = "strata_salamander.tsv",
#'   path.folder               = "salamander",
#'   filter.individuals.missing = "outliers",
#'   filter.common.markers      = TRUE,
#'   filter.strands             = "blacklist",
#'   filter.ma                  = 4,
#'   filter.genotyping          = 0.3,
#'   filter.snp.position.read   = "outliers",
#'   filter.short.ld            = "mac",
#'   filter.long.ld             = NULL,
#'   verbose                    = TRUE
#' )
#'
#' # Later, in a new R session, reopen the GDS:
#' gds <- radiator::read_rad(data = "radiator_20200911@0748.gds")
#' }
#'
#' @export
read_vcf <- function(
    data,
    strata        = NULL,
    filename      = NULL,
    vcf.stats     = FALSE,
    parallel.core = parallel::detectCores() - 1,
    verbose       = TRUE,
    ...
) {

  # Common startup -------------------------------------------------------------
  .start   <- radiator_startup(f.name = "read_vcf", verbose = verbose)
  file.date <- .start$file.date
  on.exit(radiator_teardown(.start), add = TRUE)
  res <- list()

  # Required package ----------------------------------------------------------
  radiator_packages_dep(package = "SeqArray", cran = FALSE, bioc = TRUE)

  # Checking for missing arguments --------------------------------------------
  if (missing(data)) {
    rlang::abort("vcf file missing")
  }

  # Function call and dotslist ------------------------------------------------
  rad.dots <- radiator_dots(
    func.name = as.list(sys.call())[[1]],
    fd        = rlang::fn_fmls_names(),
    args.list = as.list(environment()),
    dotslist  = rlang::dots_list(
      ...,
      .homonyms     = "error",
      .check_assign = TRUE
    ),
    keepers    = dots_keepers_read_vcf(),
    deprecated = dots_deprecated_read_vcf(),
    env        = environment(),
    assign     = TRUE,
    verbose    = FALSE
  )

  # Normalise optional ... arguments: create as NULL if missing
  dots.optional <- c(
    "whitelist.markers",
    "blacklist.id",
    "pop.select",
    "pop.levels",
    "pop.labels",
    "filter.strands",
    "filter.common.markers",
    "filter.ma",
    "filter.coverage",
    "filter.genotyping",
    "filter.snp.position.read",
    "filter.snp.number",
    "filter.short.ld", "filter.long.ld", "ld.method",
    "filter.individuals.missing",
    "filter.individuals.coverage.total",
    "filter.individuals.coverage.median",
    "filter.individuals.coverage.iqr",
    "markers.info",
    "vcf.metadata",
    "path.folder",
    "random.seed",
    "subsample.markers.stats",
    "parameters",
    "filter.monomorphic",
    "filter.haplotype.format",
    "long.ld.missing",
    "filter.individuals.heterozygosity",  # downstream use
    "internal"
  )

  ## capture the environment we want to populate
  target_env <- environment()

  missing_vars <- !vapply(
    dots.optional,
    exists,
    logical(1),
    envir    = target_env,
    inherits = FALSE
  )

  purrr::walk(
    dots.optional[missing_vars],
    ~ assign(.x, NULL, envir = target_env)
  )

  # path.folder: default to getwd() if still NULL
  if (is.null(path.folder)) {
    path.folder <- getwd()
  }

  # ld.method
  ld.allowed <- c("composite", "r", "r2", "dprime", "corr")

  if (!is.null(ld.method)) {
    ld.method <- match.arg(ld.method, ld.allowed)
  }

  # LD: if long LD filter requested but no short LD, force short LD = "mac"
  if (!is.null(filter.long.ld) && is.null(filter.short.ld)) {
    message("\nfilter.short.ld argument set by default to: mac")
    filter.short.ld <- "mac"
  }

  # filter.snp.position.read
  if (!is.null(filter.snp.position.read)) {
    filter.snp.position.read <- match.arg(
      filter.snp.position.read,
      c("outliers", "iqr", "q75")
    )
  }

  # filter.strands
  if (!is.null(filter.strands)) {
    filter.strands <- match.arg(
      filter.strands,
      c("blacklist", "whitelist", "ignore")
    )
  } else {
    filter.strands <- "blacklist"
  }

  # monomorphic and common markers ...
  if (is.null(filter.monomorphic)) {
    filter.monomorphic <- TRUE
  }

  if (is.null(filter.common.markers)) {
    filter.common.markers <- TRUE
  }


  # internal
  if (is.null(internal)) internal <- FALSE

  # If any heavy filtering is requested, we want full statistics (no subsample)
  # is_set <- function(x) !is.null(x)
  #
  # heavy.filters <- c(
  #   filter.snp.position.read,
  #   filter.ma,
  #   filter.coverage,
  #   filter.genotyping,
  #   filter.short.ld,
  #   filter.long.ld,
  #   long.ld.missing,
  #   filter.individuals.missing,
  #   filter.individuals.coverage.total,
  #   filter.snp.number
  # )

  # if (any(is_set(heavy.filters))) {
  #   message("stats by default to 1")
  #   subsample.markers.stats <- 1
  # }

  # Folders -------------------------------------------------------------------
  prefix.int <- FALSE
  if (TRUE %in% unique(stringi::stri_detect_fixed(str = path.folder, pattern = "@"))) {
    prefix.int <- TRUE
    internal   <- FALSE
  }

  wf <- path.folder <- generate_folder(
    rad.folder  = "read_vcf",
    path.folder = path.folder,
    prefix.int  = prefix.int,
    internal    = internal,
    file.date   = file.date,
    verbose     = verbose
  )

  radiator.folder <- generate_folder(
    rad.folder  = "import_gds",
    path.folder = path.folder,
    prefix.int  = TRUE,
    internal    = internal,
    file.date   = file.date,
    verbose     = verbose
  )

  # Store function call + args ------------------------------------------------
  write_radiator_tsv(
    data          = rad.dots,
    path.folder   = path.folder,
    filename      = "radiator_read_vcf_args",
    date          = TRUE,
    internal      = internal,
    write.message = "Function call and arguments stored in: ",
    verbose       = FALSE
  )

  if (is.null(filename)) {
    ind.file            <- stringi::stri_join("vcf_individuals_info_", file.date, ".tsv")
    markers.file        <- stringi::stri_join("vcf_markers_metadata_", file.date, ".tsv")
    blacklist.markers   <- stringi::stri_join("blacklist.markers_", file.date, ".tsv")
    blacklist.id.filename <- stringi::stri_join("blacklist.individuals_", file.date, ".tsv")
  } else {
    ind.file            <- stringi::stri_join(filename, "_vcf_individuals_info_", file.date, ".tsv")
    markers.file        <- stringi::stri_join(filename, "_vcf_markers_metadata_", file.date, ".tsv")
    blacklist.markers   <- stringi::stri_join(filename, "_blacklist.markers_", file.date, ".tsv")
    blacklist.id.filename <- stringi::stri_join(filename, "_blacklist.individuals_", file.date, ".tsv")
  }

  filename <- generate_filename(
    name.shortcut = filename,
    path.folder   = radiator.folder,
    extension     = "gds"
  )

  filename.short <- filename$filename.short
  filename       <- filename$filename

  # Random seed ---------------------------------------------------------------
  if (is.null(random.seed)) {
    random.seed <- sample(x = 1:1000000, size = 1)
  }
  set.seed(random.seed)

  readr::write_lines(
    x    = random.seed,
    file = file.path(radiator.folder, "random.seed")
  )
  if (verbose) {
    message("File written: random.seed (", random.seed, ")")
  }

  # Check VCF indexing ---------------------------------------------------------
  file.origin <- data

  if (detect_indexing(vcf = data)) {
    if (verbose) cli::cli_alert_warning("VCF is not ready for parallel reading: running bgzip + tabix")
    data <- indexing_vcf(vcf = data)
  }

  # VCF HEADER info ------------------------------------------------------------
  detect.source <- check_header_source_vcf(
    vcf          = data,
    markers.info = if (exists("markers.info")) markers.info else NULL,
    vcf.metadata = if (exists("vcf.metadata")) vcf.metadata else NULL
  )
  stacks.checks <- detect.source$stacks.check
  data.source <- detect.source$data.source
  check.header <- detect.source$check.header
  dp <- "DP" %in% check.header$format$ID # Check that DP is valid
  markers.info <- detect.source$markers.info
  overwrite.metadata <- detect.source$overwrite.metadata # may be NULL
  # vcf.source.raw <- detect.source$vcf.source.raw # optional, if you want it later

  if (verbose) message("VCF source: ", data.source)

  # VCF -> GDS using SeqArray --------------------------------------------------

  # VCF size
  parallel.temp <- parallel.core
  # For small VCF, SeqArray parallel overhead is not worth it
  if (file.size(data) < 10000000) {
    parallel.temp <- FALSE
  }

  cli::cli_progress_step("Reading VCF")

  vcf_read_temp <- function(
    data,
    filename,
    parallel.temp,
    check.header,
    markers.info,
    overwrite.metadata
  ) {
    SeqArray::seqVCF2GDS(
      vcf.fn         = data,
      out.fn         = filename,
      parallel       = parallel.temp,
      storage.option = "ZIP_RA",
      verbose        = FALSE,
      header         = check.header,
      info.import    = markers.info,        # INFO fields (NULL = all)
      fmt.import     = overwrite.metadata   # FORMAT fields (NULL = all)
    )
  }

  safe_vcf_read <- purrr::safely(.f = vcf_read_temp)

  data.safe <- safe_vcf_read(
    data              = data,
    filename          = filename,
    parallel.temp     = parallel.temp,
    check.header      = check.header,
    markers.info      = markers.info,
    overwrite.metadata = overwrite.metadata
  )

  if (is.null(data.safe$error)) {
    gds <- SeqArray::seqOpen(gds.fn = data.safe$result, readonly = FALSE)
  } else {

    # When an error is detected while reading the VCF...
    cli::cli_alert_warning("Parallel read failed, turning parallel OFF")

    os <- Sys.info()[["sysname"]]
    if (os == "Windows") {
      cli::cli_alert_info("Common problem with Windows machines")
    } else {
      cli::cli_alert_info("Not normal with UNIX machines: check your SeqArray installation")
    }

    gds <- SeqArray::seqVCF2GDS(
      vcf.fn         = data,
      out.fn         = filename,
      parallel       = FALSE,
      storage.option = "ZIP_RA",
      verbose        = FALSE,
      header         = check.header,
      info.import    = markers.info,
      fmt.import     = overwrite.metadata
    ) %>%
      SeqArray::seqOpen(gds.fn = ., readonly = FALSE)
  }

  # no longer needed
  vcf_read_temp <- data.safe <- safe_vcf_read <- parallel.temp <- NULL
  check.header <- detect.source <- markers.info <- overwrite.metadata <- NULL

  cli::cli_progress_done()


  if (verbose) message("Analyzing VCF")

  # Radiator skeleton ---------------------------------------------------------
  radiator.gds <- radiator_gds_skeleton(gds)

  # Track input file and source -----------------------------------------------
  update_radiator_gds(gds = gds, node.name = "input.file",  value = file.origin)
  update_radiator_gds(gds = gds, node.name = "data.source", value = data.source)


  # Sample IDs -----------------------------------------------------------------
  individuals.vcf <- tibble::tibble(
    INDIVIDUALS_VCF   = SeqArray::seqGetData(gds, "sample.id")
  ) %>%
    dplyr::mutate(
      INDIVIDUALS_CLEAN = radiator::clean_ind_names(INDIVIDUALS_VCF)
    )

  if (!identical(individuals.vcf$INDIVIDUALS_VCF,
                 individuals.vcf$INDIVIDUALS_CLEAN)) {
    if (verbose) message("Cleaning VCF's sample names")
    clean.id.filename <- stringi::stri_join("cleaned.vcf.id.info_", file.date, ".tsv")
    readr::write_tsv(
      x    = individuals.vcf,
      file = file.path(radiator.folder, clean.id.filename)
    )
    if (verbose) message("File written: ", clean.id.filename)

    update_radiator_gds(
      gds        = gds,
      node.name  = "id.clean",
      value      = individuals.vcf
    )
  }

  # Replace id in GDS ---------------------------------------------------------
  update_radiator_gds(
    gds         = gds,
    radiator.gds = FALSE,
    node.name   = "sample.id",
    value       = individuals.vcf$INDIVIDUALS_CLEAN,
    replace     = TRUE
  )

  individuals <- dplyr::select(
    individuals.vcf,
    INDIVIDUALS = INDIVIDUALS_CLEAN
  )
  individuals.vcf <- NULL

  # Strata handling -----------------------------------------------------------
  strata <- radiator::read_strata(
    strata      = strata,
    pop.levels  = pop.levels,
    pop.labels  = pop.labels,
    pop.select  = pop.select,
    blacklist.id = blacklist.id,
    keep.two    = FALSE
  ) %$% strata

  if (!is.null(strata)) {
    id.levels <- individuals$INDIVIDUALS


    # highlight samples in the strata but not in the vcf...
    check.strata <- strata %>%
      dplyr::filter(!INDIVIDUALS %in% individuals$INDIVIDUALS)

    # synchronize data and strata
    individuals %<>%
      dplyr::left_join(
        join_strata(individuals, strata, verbose = verbose) %>%
          dplyr::mutate(FILTERS = "whitelist"),
        by = "INDIVIDUALS"
      ) %>%
      dplyr::mutate(
        FILTERS = tidyr::replace_na(FILTERS, "filter.stata")
      )

    # samples blacklisted by the strata
    bl <- dplyr::filter(individuals, FILTERS != "whitelist")


    if (nrow(bl) != 0) {
      if (verbose) {
        cli::cli_alert_info("    Number of sample blacklisted by the strata: {nrow(bl)}")
      }
    }

    missing.samples <- FALSE
    if (nrow(check.strata) != 0) {
      if (verbose) {
        cli::cli_alert_info("    Number of sample in strata but not in data: {nrow(check.strata)}")
      }
      missing.samples <- TRUE
    }

    if (rlang::has_name(individuals, "NEW_ID")) {
      if (!identical(id.levels, individuals$INDIVIDUALS)) {
        rlang::abort("Wrong id order, contact author")
      }
      update_radiator_gds(
        gds         = gds,
        radiator.gds = FALSE,
        node.name   = "sample.id",
        value       = individuals$NEW_ID,
        replace     = TRUE
      )
      individuals %<>% dplyr::rename(
        INDIVIDUALS = NEW_ID,
        OLD_ID      = INDIVIDUALS
      )
    }
  } else {
    individuals %<>%
      dplyr::mutate(
        STRATA  = "1pop",
        FILTERS = "whitelist"
      )
  }

  strata <- generate_strata(
    data   = dplyr::filter(individuals, FILTERS == "whitelist"),
    pop.id = FALSE
  )

  if (!is.factor(individuals$STRATA)) {
    individuals$STRATA <- factor(individuals$STRATA)
  }

  individuals %<>%
    dplyr::mutate(
      ID_SEQ     = seq_len(dplyr::n()),
      STRATA_SEQ = as.integer(factor(STRATA, levels = levels(STRATA)))
    )

  update_radiator_gds(
    gds       = gds,
    node.name = "individuals.meta",
    value     = individuals,
    sync      = TRUE
  )

  # Reference genome or de novo -----------------------------------------------
  ref.genome <- detect_ref_genome(data = gds, verbose = verbose)

  # Markers metadata -----------------------------------------------------------
  markers.meta <- extract_markers_metadata(gds = gds)

  cleaned.vcf.meta <- clean_vcf_markers_meta(
    markers.meta  = markers.meta,
    gds           = gds,
    data.source   = data.source,
    ref.genome    = ref.genome,
    stacks.checks = stacks.checks,
    verbose       = verbose
  )

  markers.meta   <- cleaned.vcf.meta$markers.meta
  detect.strand  <- cleaned.vcf.meta$detect.strand
  cleaned.vcf.meta <- NULL

  update_radiator_gds(
    gds       = gds,
    node.name = "markers.meta",
    value     = markers.meta
  )

  # Replace chromosome in GDS with cleaned CHROM -------------------------------
  update_radiator_gds(
    gds         = gds,
    radiator.gds = FALSE,
    node.name   = "chromosome",
    value       = markers.meta$CHROM,
    replace     = TRUE
  )

  # Bi- vs multi-allelic -------------------------------------------------------
  biallelic <- detect_biallelic_markers(data = gds, verbose = verbose)

  # radiator_parameters: generate + initiate -----------------------------------
  filters.parameters <- radiator_parameters(
    generate      = TRUE,
    initiate      = FALSE,
    update        = FALSE,
    parameter.obj = parameters,
    path.folder   = radiator.folder,
    file.date     = file.date,
    verbose       = verbose,
    internal      = internal
  )

  filters.parameters <- radiator_parameters(
    generate      = FALSE,
    initiate      = TRUE,
    update        = TRUE,
    parameter.obj = filters.parameters,
    data          = gds,
    filter.name   = "vcf",
    param.name    = "original values in VCF + strata",
    values        = "",
    path.folder   = path.folder,
    file.date     = file.date,
    internal      = internal,
    verbose       = FALSE
  )


  # VCF-specific filters -------------------------------------------------------
if (verbose) message("VCF-specific filters")

  # 1) Duplicate SNPs on different strands -------------------------------------
  if (detect.strand) {
    strand.out <- filter_duplicated_markers_strands(
      gds                = gds,
      markers.meta       = markers.meta,
      filters.parameters = filters.parameters,
      path.folder        = path.folder,
      file.date          = file.date,
      filter.strands     = filter.strands,
      detect.strand      = detect.strand,
      parallel.core      = parallel.core,
      internal           = internal,
      verbose            = verbose
    )

    markers.meta       <- strand.out$markers.meta
    filters.parameters <- strand.out$filters.parameters
    strand.out         <- NULL
  }

  # 2) VCF FILTER column -------------------------------------------------------
  clean.vcf <- filter_vcf_filter_column(
    gds                = gds,
    markers.meta       = markers.meta,
    filters.parameters = filters.parameters,
    path.folder        = wf,
    file.date          = file.date,
    verbose            = verbose
  )

  filters.parameters <- clean.vcf$filters.parameters
  markers.meta       <- NULL
  clean.vcf          <- NULL

  # 3) Haplotype / non-biallelic cleanup ---------------------------------------
  haplo.out <- filter_haplotype_format(
    gds                     = gds,
    filters.parameters      = filters.parameters,
    path.folder             = path.folder,
    file.date               = file.date,
    filter.haplotype.format = filter.haplotype.format,
    internal                = internal,
    verbose                 = verbose
  )

  filters.parameters <- haplo.out$filters.parameters
  haplo.out          <- NULL

  if (isTRUE(filter.haplotype.format)) {
    biallelic <- detect_biallelic_markers(data = gds, verbose = verbose)
  }

  if (!biallelic && stacks.checks) {
    dp <- FALSE
  }


  # Optional radiator filters---------------------------------------------------
  # Monomorphic markers --------------------------------------------------------
  gds <- filter_monomorphic(
    data               = gds,
    filter.monomorphic = filter.monomorphic,
    parallel.core      = parallel.core,
    verbose            = FALSE,
    parameters         = filters.parameters,
    path.folder        = wf,
    internal           = FALSE
  )

  # Common markers -------------------------------------------------------------
  gds <- filter_common_markers(
    data                  = gds,
    filter.common.markers = filter.common.markers,
    fig                   = TRUE,
    parallel.core         = parallel.core,
    verbose               = FALSE,
    parameters            = filters.parameters,
    path.folder           = wf,
    internal              = FALSE
  )

  # Whitelist markers ----------------------------------------------------------
  gds <- filter_whitelist(
    data              = gds,
    whitelist.markers = whitelist.markers,
    verbose           = FALSE,
    path.folder       = wf,
    parameters        = filters.parameters,
    biallelic         = biallelic,
    internal          = FALSE
  )

  # Individuals filters --------------------------------------------------------
  gds <- filter_individuals(
    data                               = gds,
    interactive.filter                 = FALSE,
    filter.individuals.missing         = filter.individuals.missing,
    filter.individuals.heterozygosity  = NULL,
    filter.individuals.coverage.total  = filter.individuals.coverage.total,
    filter.individuals.coverage.median = filter.individuals.coverage.median,
    filter.individuals.coverage.iqr    = filter.individuals.coverage.iqr,
    parallel.core                      = parallel.core,
    verbose                            = FALSE,
    path.folder                        = wf,
    parameters                         = filters.parameters,
    internal                           = FALSE
  )

  # MAC filter -----------------------------------------------------------------
  gds <- filter_ma(
    data               = gds,
    interactive.filter = FALSE,
    filter.ma          = filter.ma,
    filename           = NULL,
    parallel.core      = parallel.core,
    verbose            = FALSE,
    parameters         = filters.parameters,
    path.folder        = wf,
    internal           = FALSE
  )

  # Coverage filter ------------------------------------------------------------
  gds <- filter_coverage(
    data               = gds,
    interactive.filter = FALSE,
    filter.coverage    = filter.coverage,
    filename           = NULL,
    parallel.core      = parallel.core,
    verbose            = FALSE,
    parameters         = filters.parameters,
    path.folder        = wf,
    internal           = FALSE
  )

  # Genotyping rate filter -----------------------------------------------------
  gds <- filter_genotyping(
    data               = gds,
    interactive.filter = FALSE,
    filter.genotyping  = filter.genotyping,
    filename           = NULL,
    parallel.core      = parallel.core,
    verbose            = FALSE,
    parameters         = filters.parameters,
    path.folder        = wf,
    internal           = FALSE
  )

  # SNP position on read -------------------------------------------------------
  gds <- filter_snp_position_read(
    data                 = gds,
    interactive.filter   = FALSE,
    filter.snp.position.read = filter.snp.position.read,
    filename             = NULL,
    parallel.core        = parallel.core,
    verbose              = FALSE,
    parameters           = filters.parameters,
    path.folder          = wf,
    internal             = FALSE
  )

  # Number of SNPs per locus ---------------------------------------------------
  gds <- filter_snp_number(
    data               = gds,
    interactive.filter = FALSE,
    filter.snp.number  = filter.snp.number,
    filename           = NULL,
    parallel.core      = parallel.core,
    verbose            = FALSE,
    parameters         = filters.parameters,
    path.folder        = wf,
    internal           = FALSE
  )

  # LD filters -----------------------------------------------------------------
  gds <- filter_ld(
    data               = gds,
    interactive.filter = FALSE,
    filter.short.ld    = filter.short.ld,
    filter.long.ld     = filter.long.ld,
    parallel.core      = parallel.core,
    filename           = NULL,
    verbose            = FALSE,
    long.ld.missing    = long.ld.missing,
    ld.method          = ld.method,
    parameters         = filters.parameters,
    path.folder        = wf,
    internal           = FALSE
  )


  # Final sync and outputs -----------------------------------------------------

  if (verbose) message("Preparing output files...")

  markers.meta <- extract_markers_metadata(gds, whitelist = TRUE)
  strata       <- extract_individuals_metadata(
    gds                = gds,
    ind.field.select   = c("INDIVIDUALS", "STRATA"),
    whitelist          = TRUE
  )

  sync_gds(gds = gds)

  path.folder.filtered <- generate_folder(
    rad.folder  = "filtered",
    path.folder = wf,
    internal    = FALSE,
    file.date   = file.date,
    verbose     = verbose
  )

  # Whitelist markers ----------------------------------------------------------
  write_radiator_tsv(
    data          = markers.meta,
    path.folder   = path.folder.filtered,
    filename      = "whitelist.markers",
    date          = TRUE,
    internal      = internal,
    write.message = "standard",
    verbose       = verbose
  )

  # Blacklist markers ----------------------------------------------------------
  bl <- extract_markers_metadata(gds, blacklist = TRUE)
  if (nrow(bl) > 0) {
    write_radiator_tsv(
      data          = bl,
      path.folder   = path.folder.filtered,
      filename      = "blacklist.markers",
      date          = TRUE,
      internal      = internal,
      write.message = "standard",
      verbose       = verbose
    )
  }

  # Blacklist individuals ------------------------------------------------------
  blacklist.id.out <- extract_individuals_metadata(gds, blacklist = TRUE)
  if (nrow(blacklist.id.out) > 0) {
    write_radiator_tsv(
      data          = blacklist.id.out,
      path.folder   = path.folder.filtered,
      filename      = "blacklist.id",
      date          = TRUE,
      internal      = internal,
      write.message = "standard",
      verbose       = verbose
    )
  }

  if (missing.samples) {
    write_radiator_tsv(
      data          = check.strata,
      path.folder   = path.folder.filtered,
      filename      = "check.strata.missing.id.in.data",
      date          = TRUE,
      internal      = internal,
      write.message = "standard",
      verbose       = verbose
    )
  }



  # Filtered strata ------------------------------------------------------------
  write_radiator_tsv(
    data          = strata,
    path.folder   = path.folder.filtered,
    filename      = "strata.filtered",
    date          = TRUE,
    internal      = internal,
    write.message = "standard",
    verbose       = verbose
  )


  # VCF stats (optional) -------------------------------------------------------

  if (vcf.stats) {
    message("statistics was turned on somehow...")

    # large VCFs not really the tools or time to do that
    # subsampling the number of SNPs is accurate
    # will give us an idea of the state of the VCF...
    # Even subsampling 10% gives an pretty accurate picture...

    n.markers <- length(markers.meta$VARIANT_ID)
    n.ind     <- length(strata$INDIVIDUALS)
    n.snp     <- as.numeric(n.markers)

    if (is.null(subsample.markers.stats)) {
      subsample.markers.stats <- choose_markers_subsample_prop(
        n.markers        = n.snp,
        n.ind            = n.ind,
        vector.size.limit = 2^31,
        max.mem.gb       = 10,
        bytes.per.cell   = 4,
        n.matrices       = 4,      # set this to what your vcf.stats step actually holds
        small.vcf.cutoff = 200000,
        cap.prop = 0.3
      )
    }

    # sanitize user-provided values too
    subsample.markers.stats <- as.numeric(subsample.markers.stats)
    if (!is.finite(subsample.markers.stats) || subsample.markers.stats <= 0) subsample.markers.stats <- 0.10
    if (subsample.markers.stats > 1) subsample.markers.stats <- 1

    if (subsample.markers.stats < 1) {
      if (verbose) {
        message(
          "Subsampling turned ON for fast coverage statistics (prop = ",
          subsample.markers.stats, ")"
        )
      }
      markers.subsampled <- withr::with_seed(
        seed = random.seed,
        code = dplyr::slice_sample(markers.meta, prop = subsample.markers.stats)
      )

      variant.select <- markers.subsampled$VARIANT_ID

      subsample.filename <- stringi::stri_join(
        "markers.subsampled_",
        file.date,
        ".tsv"
      )

      dplyr::tibble(
        VARIANT_ID  = markers.subsampled$VARIANT_ID,
        RANDOM_SEED = random.seed
      ) %>%
        readr::write_tsv(file = file.path(path.folder.filtered, subsample.filename))

      markers.subsampled <- subsample.filename <- NULL
    } else {
      variant.select <- NULL
    }


    # currently coverage stats will not work with multi allelic data...
    if (ref.genome) {
      snp.per.locus.check <- FALSE
    } else {
      snp.per.locus.check <- TRUE
    }

    # stats...
    i.m.stats <- generate_stats(
      gds          = gds,
      snp.per.locus = snp.per.locus.check,
      subsample    = variant.select,
      exhaustive   = TRUE,
      force.stats  = TRUE,
      path.folder  = path.folder.filtered,
      plot         = TRUE,
      digits       = 6,
      file.date    = file.date,
      parallel.core = parallel.core,
      verbose      = verbose
    )

    ind.missing <- round(
      i.m.stats$i.stats$MEAN[i.m.stats$i.stats$GROUP == "missing genotypes"],
      2
    )
    if (dp) {
      ind.cov.total <- round(
        i.m.stats$i.stats$MEAN[i.m.stats$i.stats$GROUP == "total coverage"], 0
      )
      ind.cov.mean <- round(
        i.m.stats$i.stats$MEAN[i.m.stats$i.stats$GROUP == "mean coverage"], 0
      )
    }
    markers.missing <- round(
      i.m.stats$m.stats$MEAN[i.m.stats$m.stats$GROUP == "missing genotypes"],
      2
    )
    if (dp) {
      markers.cov <- round(
        i.m.stats$m.stats$MEAN[i.m.stats$m.stats$GROUP == "mean coverage"],
        0
      )
    }
  }


  # Summary --------------------------------------------------------------------

  number.info <- SeqArray::seqGetFilter(gdsfile = gds)


  if (verbose) {
    if (vcf.stats) {
      message("statistics was turned on somehow...")

      message("\nVCF summary\nMissing data: ")
      message("    markers: ", markers.missing)
      message("    individuals: ", ind.missing)
      if (dp) {
        message("\n\nCoverage info:")
        message("    individuals mean total coverage: ", ind.cov.total)
        message("    individuals mean genotype coverage: ", ind.cov.mean)
        message("    markers mean coverage: ", markers.cov)
      }
    }
  }

  message("\nVCF info:")
  message("Number of chromosome/contig/scaffold: ",
          length(unique(markers.meta$CHROM)))
  message("Number of locus: ", length(unique(markers.meta$LOCUS)))
  message("Number of markers: ",
          length(number.info$variant.sel[number.info$variant.sel]))
  summary_strata(strata)

  message("radiator Genomic Data Structure (GDS) file: ", basename(filename))

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
#' GATK with a mix of multi- and bi-allelic dataset won't generate VCF stats.
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
#' in the working directory.
#'
#' \strong{General arguments: }
#' \enumerate{
#' \item \code{path.folder}: to write ouput in a specific path
#' (used internally in radiator). Default: \code{path.folder = getwd()}.
#' If the supplied directory doesn't exist, it's created.
#' \item \code{random.seed}: (integer, optional) For reproducibility, set an integer
#' that will be used inside codes that uses randomness. With default,
#' a random number is generated, printed and written in the appropriate directory.
#' Random seed is recycled inside the function that will import the VCF file before
#' tidying.
#' Default: \code{random.seed = NULL}.
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
#' \item \code{calibrate.alleles: } (optional, logical)
#' Default: \code{calibrate.alleles = FALSE}.
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
#' \item \code{filter.ma}: (integer)
#' Default: \code{filter.ma = NULL}.
#' Documented in \code{\link[radiator]{filter_ma}}.
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
#'     filter.ma = 4,
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
  # filter.ma = 4
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
  # calibrate.alleles = FALSE
  # random.seed = NULL
  # parameters = NULL
  # subsample.markers.stats <- 1
  # ld.method <- "r2"
  # tidy.vcf <- TRUE
  # tidy.check <- TRUE


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
  radiator_packages_dep(package = "SeqArray", cran = FALSE, bioc = TRUE)

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
      "filter.ma", "ma.stats",
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
      "calibrate.alleles",
      "vcf.metadata", "vcf.stats",
      "whitelist.markers",
      "internal",
      "tidy.check",
      "tidy.vcf"
    ),
    verbose = FALSE
  )

  if (!is.null(filter.snp.position.read) ||
      !is.null(filter.ma) ||
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
    rad.folder = "tidy_vcf",
    path.folder = path.folder,
    prefix.int = FALSE,
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
    internal = FALSE,
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
    filter.ma = filter.ma,
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
    rad.folder = "tidy_vcf",
    path.folder = path.folder,
    prefix.int = TRUE,
    internal = FALSE,
    file.date = file.date,
    verbose = verbose)

  # write the dots file: after the GDS import...
  write_radiator_tsv(
    data = rad.dots,
    path.folder = path.folder,
    filename = "radiator_tidy_vcf_args",
    date = TRUE,
    internal = internal,
    write.message = "Function call and arguments stored in: ",
    verbose = verbose
  )

  # Random seed ----------------------------------------------------------------
  readr::write_lines(x = random.seed, file = file.path(tidy.folder, "random.seed"))
  if (verbose) message("File written: random.seed (", random.seed,")")

  # Tidy the data --------------------------------------------------------------
  if (!tidy.vcf) return(data)

  # Get info markers and individuals -----------------------------------------
  markers.meta <- extract_markers_metadata(gds = data, whitelist = TRUE)
  # markers.meta <- extract_markers_metadata(gds = data, markers.meta.select = "M_SEQ", whitelist = TRUE)
  n.markers <- nrow(markers.meta)

  # STRATEGY -----------------------------------------------------------------
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
      message("Tidying the large vcf...")
    } else {
      message("\nKeeping the GDS object/file")
      return(data)
    }
  }

  # Print genotypes tidying
  gt.bin <- TRUE
  message("\nGenotypes formats generated with ", n.markers, " SNPs: ")
  message("    GT_BIN (the dosage of ALT allele: 0, 1, 2 NA): ", gt.bin)
  message("    GT_VCF (the genotype coding VCFs: 0/0, 0/1, 1/1, ./.): ", gt.vcf)
  message("    GT_VCF_NUC (the genotype coding in VCFs, but with nucleotides: A/C, ./.): ", gt.vcf.nuc)
  message("    GT (the genotype coding 'a la genepop': 001002, 001001, 000000): ", gt)


  if (!is.null(blacklist.id)) calibrate.alleles <- TRUE

  tidy.data <- gds2tidy(
    gds = data,
    markers.meta = markers.meta,
    strip.rad = TRUE,
    pop.id = FALSE,
    calibrate.alleles = FALSE # not done here
  )

  # bi- or multi-alllelic VCF ------------------------------------------------
  biallelic <- detect_biallelic_markers(data = data)

  if (is.logical(vcf.metadata)) {
    overwrite.metadata <- "GT"
    if (vcf.metadata) overwrite.metadata <- NULL
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
  # vcf.metadata ---------------------------------------------------------------
  if (vcf.metadata) {
    if (verbose) message("\nKeeping vcf genotypes metadata: yes")

    # detect FORMAT fields available
    have <-  SeqArray::seqSummary(
      gdsfile = data,
      varname = "annotation/format",
      check = "none", verbose = FALSE)$ID

    if (length(have) > 0) {
      want <- c("DP", "AD", "GL", "PL", "HQ", "GQ", "GOF", "NR", "NV", "CATG")
      if (!is.null(overwrite.metadata)) want <- overwrite.metadata
      parse.format.list <- purrr::keep(.x = have, .p = have %in% want)
      # work on parallelization of this part
      meta <- parse_gds_metadata(x = parse.format.list, gds = data, strip.rad = TRUE, verbose = verbose)
      tidy.data %<>% dplyr::left_join(meta, by = intersect(colnames(tidy.data), colnames(meta)))
      meta <- NULL
    } else {
      if (verbose) message("    genotypes metadata: none found")
      vcf.metadata <- FALSE
    }

    # Check stacks AD problem
    # Some genotypes with missing AD...

    # NOTE TO MYSELF: might be faster to screen stacks here in data.source...

    if (rlang::has_name(tidy.data, "C_DEPTH")) {
      ref <- extract_markers_metadata(gds = data, markers.meta.select = c("M_SEQ", "REF", "ALT"), whitelist = TRUE)
      tidy.data %<>% dplyr::left_join(ref, by = "M_SEQ")
      ref <- NULL

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

      # temp <- tidy.data %>%
      #   dplyr::group_by(M_SEQ) %>%
      #   dplyr::summarise(
      #     A_DEPTH = sum(A_DEPTH, na.rm = TRUE),
      #     C_DEPTH = sum(C_DEPTH, na.rm = TRUE),
      #     G_DEPTH = sum(G_DEPTH, na.rm = TRUE),
      #     T_DEPTH = sum(T_DEPTH, na.rm = TRUE)
      #   ) %>%
      #   radiator::rad_long(cols = "M_SEQ", names_to = "ACGT_DEPTH", values_to = "DEPTH", tidy = TRUE) %>%
      #   dplyr::arrange(M_SEQ, -DEPTH) %>%
      #   dplyr::group_by(M_SEQ) %>%
      #   dplyr::slice_head(n = 2) %>%
      #   dplyr::ungroup() #%>% dplyr::bind_rows()
      #
      #
      # # check that we have 2
      # n.row <- nrow(temp)
      # check <- length(unique(temp$M_SEQ))
      #
      # if (check != n.row / 2) {
      #   rlang::abort("Contact author problem with Allele depth")
      # }
      #
      # temp %<>%
      #   dplyr::mutate(
      #     NUCLEOTIDE = stringi::stri_sub(ACGT_DEPTH, from = 1, to = 1),
      #     ACGT_DEPTH = NULL,
      #     ALLELE = rep(x = c("REF", "ALT"), times = n.row / 2),
      #     ALLELE_DEPTH = rep(x = c("ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH"), times = n.row / 2)
      #   ) %>%
      #   dplyr::group_by(M_SEQ) %>%
      #   radiator::rad_wide(names_from = c("ALLELE", "ALLELE_DEPTH"), values_from = c("NUCLEOTIDE", "DEPTH"), tidy = TRUE) %>%
      #   dplyr::rename(REF = NUCLEOTIDE_REF_ALLELE_REF_DEPTH, ALT = NUCLEOTIDE_ALT_ALLELE_ALT_DEPTH, ALLELE_REF_DEPTH = DEPTH_REF_ALLELE_REF_DEPTH, ALLELE_ALT_DEPTH = DEPTH_ALT_ALLELE_ALT_DEPTH)
      #
      #
      # tidy.data %<>% dplyr::left_join(temp, by = "M_SEQ")
      # n.row <- check <- temp <- NULL
    }#catg.depth


    if (all(rlang::has_name(tidy.data, c("ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH", "READ_DEPTH")))) {
      stacks.ad <- tidy.data %>%
        dplyr::select(READ_DEPTH, ALLELE_REF_DEPTH, ALLELE_ALT_DEPTH) %>%
        dplyr::filter(!is.na(READ_DEPTH)) %>%
        dplyr::filter(is.na(ALLELE_REF_DEPTH) & is.na(ALLELE_ALT_DEPTH))
      stacks.prob <- nrow(stacks.ad)

      if (stacks.prob > 0) {
        non.missing.gt <- length(tidy.data$GT_BIN[!is.na(tidy.data$GT_BIN)])
        stacks.ad.prop <- round(stacks.prob / non.missing.gt, 4)

        message("\n\nStacks problem detected")
        message("    missing allele depth info")
        message("    number of genotypes with problem: ", stacks.prob, " (prop: ", stacks.ad.prop,")")
        message("    correcting problem by adding the read depth info into AD fields...\n\n")

        stacks.ad <- dplyr::select(tidy.data, GT_BIN) %>%
          dplyr::bind_cols(dplyr::select(tidy.data, READ_DEPTH, ALLELE_REF_DEPTH, ALLELE_ALT_DEPTH)) %>%
          dplyr::mutate(COL_ID = seq(1, n())) %>%
          dplyr::filter(!is.na(READ_DEPTH)) %>%
          dplyr::filter(is.na(ALLELE_REF_DEPTH) & is.na(ALLELE_ALT_DEPTH)) %>%
          dplyr::mutate(
            ALLELE_REF_DEPTH = dplyr::if_else(GT_BIN == 0, READ_DEPTH, ALLELE_REF_DEPTH),
            ALLELE_ALT_DEPTH = dplyr::if_else(GT_BIN == 2, READ_DEPTH, ALLELE_ALT_DEPTH)
          )

        tidy.data$ALLELE_REF_DEPTH[stacks.ad$COL_ID] <- stacks.ad$ALLELE_REF_DEPTH
        tidy.data$ALLELE_ALT_DEPTH[stacks.ad$COL_ID] <- stacks.ad$ALLELE_ALT_DEPTH
      }
      stacks.prob <- stacks.prob <- NULL
    }
  }

  # re-calibration of ref/alt alleles ------------------------------------------
  if (calibrate.alleles) {
    tidy.data %<>% radiator::calibrate_alleles(data = ., biallelic = biallelic, verbose = verbose) %$% input
  }

  # Join back the info ---------------------------------------------------------
  tidy.data %<>%
    join_rad(
      x = .,
      s = extract_individuals_metadata(gds = data, whitelist = TRUE),
      m = extract_markers_metadata(gds = data, whitelist = TRUE),
      g = NULL,
      env.arg = rlang::current_env()
    )

  filename.rad <- generate_filename(
    path.folder = tidy.folder,
    extension = "arrow.parquet")
  write_rad(data = tidy.data, filename = filename.rad$filename, verbose = verbose)

  if (verbose) message("Updating GDS with genotypes.meta values")
  update_radiator_gds(gds = data, node.name = "genotypes.meta", value = tidy.data)
  # message("\nTidy data file written: ", filename.rad$filename.short)
  if (verbose) message("Closing GDS file connection")
  SeqArray::seqClose(data) # close the connection
  return(tidy.data)
}#End tidy_vcf


# Internal nested Function -----------------------------------------------------
# Context-specific keepers/deprecated ------------------------------------------

#' @title Dot-keepers for read_vcf
#' @description
#' Return the set of `...` arguments supported by \code{read_vcf()}.
#' This builds on the global core keepers and adds/removes any
#' read_vcf-specific arguments.
#' @keywords internal
#' @export
dots_keepers_read_vcf <- function(extra = NULL, exclude = NULL) {
  base <- c(
    dots_keepers_core(),
    # explicitly used inside read_vcf():
    "whitelist.markers",
    "filter.ma", "filter.snp.position.read", "filter.snp.number",
    "filter.coverage", "filter.genotyping", "filter.short.ld",
    "filter.long.ld", "long.ld.missing", "ld.method",
    "filter.individuals.missing", "filter.individuals.coverage.total",
    "filter.individuals.coverage.median",
    "filter.individuals.coverage.iqr",
    "filter.common.markers", "filter.monomorphic",
    "filter.strands", "filter.haplotype.format",
    "blacklist.id", "pop.select", "pop.levels", "pop.labels",
    "path.folder",
    "markers.info", "vcf.metadata",
    "subsample.markers.stats",
    "parameters", "random.seed", "internal"
  )
  # read_vcf has vcf.stats as a formal arg; prevent .../parameters from overriding it
  base <- setdiff(base, "vcf.stats")

  if (!is.null(extra)) {
    base <- c(base, extra)
  }
  if (!is.null(exclude)) {
    base <- setdiff(base, exclude)
  }

  unique(base)
}#END dots_keepers_read_vcf


#' @title Deprecated dot-arguments for read_vcf
#' @description
#' Return the set of deprecated `...` arguments understood by
#' \code{read_vcf()}. Currently this is just the global deprecated set,
#' but the helper allows you to add read_vcf-specific deprecated args
#' later if needed.
#' @keywords internal
#' @export
dots_deprecated_read_vcf <- function(extra = NULL, exclude = NULL) {
  base <- dots_deprecated_core()
  if (!is.null(extra)) {
    base <- c(base, extra)
  }
  if (!is.null(exclude)) {
    base <- setdiff(base, exclude)
  }
  unique(base)
}#END dots_deprecated_read_vcf


#' @title choose markers subsample prop
#' @description
#' Choose the best subsampling prop based on different metrics
#' @keywords internal
#' @rdname choose_markers_subsample_prop
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
choose_markers_subsample_prop <- function(n.markers,
                                          n.ind,
                                          vector.size.limit = 2^31,
                                          max.mem.gb = 10,
                                          bytes.per.cell = 4,
                                          n.matrices = 1,
                                          small.vcf.cutoff = 200000,
                                          cap.prop = 0.30,
                                          allowed = c(1, 0.30, 0.20, 0.10, 0.05, 0.02, 0.01),
                                          safety.margin = 0.90) {

  n.markers <- as.numeric(n.markers)
  n.ind     <- as.numeric(n.ind)

  if (n.markers < small.vcf.cutoff) return(1)

  prop.max.vector <- (vector.size.limit / (n.markers * n.ind)) * safety.margin
  prop.max.vector <- max(min(prop.max.vector, 1), 0)

  max.mem.bytes <- max.mem.gb * 1024^3
  denom.bytes   <- n.markers * n.ind * bytes.per.cell * n.matrices
  prop.max.mem  <- (max.mem.bytes / denom.bytes) * safety.margin
  prop.max.mem  <- max(min(prop.max.mem, 1), 0)

  prop.max <- min(prop.max.vector, prop.max.mem)

  # Key change: aim for the maximum safe proportion, capped at cap.prop
  desired <- min(prop.max, cap.prop)

  pick <- max(allowed[allowed <= desired], na.rm = TRUE)
  if (!is.finite(pick)) pick <- min(allowed[allowed > 0], na.rm = TRUE)

  pick
}#END choose_markers_subsample_prop

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
#' @param empty generate an empty vcf.
#' Default: \code{empty = FALSE}.
#' @param strata (optional, logical) Should the strata information be
#' included in the FORMAT field (along the GT info for each samples ?). To make
#' the VCF population-ready use \code{strata = TRUE}. The strata information
#' must be included in the \code{STRATA} column of the tidy dataset.
#' Default: \code{strata = FALSE}. Experimental.

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
    strata = FALSE,
    filename = NULL,
    source = NULL,
    empty = FALSE
) {
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")

  if (empty) {
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
    if (is.vector(data)) data <- radiator::tidy_wide(data = data, import.metadata = TRUE)

    # REF/ALT Alleles and VCF genotype format ------------------------------------
    if (!tibble::has_name(data, "GT_VCF")) {
      data %<>% radiator::calibrate_alleles(data = ., gt.vcf = TRUE, gt.vcf.nuc = TRUE) %$% input
    }

    # Include CHROM, LOCUS, POS --------------------------------------------------
    if (!tibble::has_name(data, "CHROM")) {
      data %<>%
        dplyr::mutate(
          CHROM = rep("1", n()),
          LOCUS = MARKERS,
          POS = MARKERS
        )
    }

    # Order/sort by pop and ind --------------------------------------------------
    if (tibble::has_name(data, "POP_ID")) {
      data %<>% dplyr::arrange(POP_ID, INDIVIDUALS)
    } else {
      data %<>% dplyr::arrange(INDIVIDUALS)
    }

    id.string <- unique(data$INDIVIDUALS)# keep to sort vcf columns
    # Remove the POP_ID column ---------------------------------------------------
    if (tibble::has_name(data, "POP_ID") && (!strata)) data %<>% dplyr::select(-POP_ID)

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
    if (strata) {
      output <- suppressWarnings(
        dplyr::left_join(data, info.field, by = "MARKERS") %>%
          dplyr::select(MARKERS, CHROM, LOCUS, POS, REF, ALT, INFO, INDIVIDUALS, GT_VCF, POP_ID) %>%
          dplyr::mutate(GT_VCF_POP_ID = stringi::stri_join(GT_VCF, POP_ID, sep = ":")) %>%
          dplyr::select(-c(GT_VCF, POP_ID)) %>%
          radiator::rad_wide(
            x = .,
            formula = "MARKERS + CHROM + LOCUS + POS + INFO + REF + ALT ~ INDIVIDUALS",
            names_from = "GT_VCF_POP_ID") %>%
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
          radiator::rad_wide(x = ., formula = "MARKERS + CHROM + LOCUS + POS + INFO + REF + ALT ~ INDIVIDUALS", values_from = "GT_VCF") %>%
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
      output %<>%
        dplyr::mutate(
          REF = dplyr::recode(REF, "A" = "001", "C" = "002", "G" = "003", "T" = "004"),
          ALT = dplyr::recode(ALT, "A" = "001", "C" = "002", "G" = "003", "T" = "004")
        )
    }

    if (tibble::has_name(output, "COL")) {
      output %<>% dplyr::mutate(LOCUS = stringi::stri_join(LOCUS, COL, sep = "_"))
    } else {
      output %<>% dplyr::mutate(LOCUS = stringi::stri_join(LOCUS, as.numeric(POS) - 1, sep = "_"))
    }

    # Keep the required columns
    output %<>%
      dplyr::arrange(CHROM, LOCUS, POS) %>%
      dplyr::select(-MARKERS) %>%
      dplyr::select('#CHROM' = CHROM, POS, ID = LOCUS, REF, ALT, QUAL, FILTER, INFO, FORMAT, id.string)
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
    file = filename,
    delim = " ",
    append = FALSE,
    col_names = FALSE
  )

  # File date ------------------------------------------------------------------
  x <- paste0("##fileDate=", file.date)
  readr::write_delim(
    x = tibble::tibble(x),
    file = filename,
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
      file = filename,
      delim = " ",
      append = TRUE,
      col_names = FALSE)
  } else {
    readr::write_delim(
      x = tibble::tibble(
        stringi::stri_join('##source=', rlang::quo_name(source))
      ),
      file = filename,
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
  if (strata && !empty) {
    format2 <- as.data.frame('##FORMAT=<ID=POP_ID,Number=1,Type=Character,Description="Population identification of Sample">')
    utils::write.table(x = format2, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  }

  # Write the prunned vcf to the file ------------------------------------------
  suppressWarnings(readr::write_tsv(x = output, file = filename, append = TRUE, col_names = TRUE))

}# end write_vcf


# extract_individuals_vcf-------------------------------------------------------
#' @title Extract individuals from vcf file
#' @description Function that returns the individuals present in a vcf file.
#' Useful to create a strata file or
#' to make sure you have the right individuals in your VCF.
#' @param vcf (character, path) The path to the vcf file.
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
#'     readr::write_tsv(x = ., file = "my.new.vcf.strata.tsv")
#' }

extract_individuals_vcf <- function(vcf) {
  id <- tibble::tibble(
    INDIVIDUALS = SeqArray::seqVCF_Header(vcf.fn = vcf, getnum = TRUE)$sample.id
  )
  return(id)
}#End extract_individuals_vcf


# extract_info_vcf-------------------------------------------------------
#' @title extract_info_vcf
#' @description Extract vcf information
#' @rdname extract_info_vcf
#' @keywords internal
#' @export
extract_info_vcf <- function(vcf) {
  vcf.info <- SeqArray::seqVCF_Header(vcf.fn = vcf, getnum = TRUE)
  res <- list(
    vcf.source = vcf.info$header$value[2],
    n.ind = vcf.info$num.sample,
    n.markers = vcf.info$num.variant,
    sample.id = vcf.info$sample.id
  )
  return(res)
}#End extract_info_vcf


# check_header_source_vcf-------------------------------------------------------
#' @title Check a VCF header and detect its source (caller)
#'
#' @description
#' Inspect the VCF header (via \code{SeqArray::seqVCF_Header()}) to:
#' \itemize{
#'   \item fix known FORMAT/INFO inconsistencies (Stacks, samtools/freebayes-style),
#'   \item detect the variant caller / pipeline (Stacks, freebayes, bcftools, GATK, etc.),
#'   \item flag whether the VCF comes from Stacks,
#'   \item suggest reasonable INFO/FORMAT subsets for import,
#'   \item return a cleaned header object for safe import into GDS.
#' }
#'
#' The function preserves the original \code{##source=} line (if present)
#' and returns it as \code{vcf.source.raw}.
#' Classification into a standardised descriptor is returned as \code{data.source}.
#'
#' @param vcf (character) Path to a VCF file (optionally bgzipped).
#'
#' @param markers.info (optional, character)
#' Names of INFO fields to import. If \code{NULL}, caller-specific defaults are
#' used when available; if those are also \code{NULL}, all INFO fields are
#' imported. Any values not present in the header are silently dropped.
#'
#' @param vcf.metadata (optional)
#' Controls which FORMAT fields are imported:
#' \itemize{
#'   \item \code{NULL}: use caller-specific defaults (if any), otherwise import
#'         all FORMAT fields;
#'   \item logical: \code{TRUE} = import all FORMAT fields,
#'                  \code{FALSE} = import only \code{GT};
#'   \item character: explicit list of FORMAT field IDs to import; \code{GT} is
#'         always added if missing. Any fields not present in the header are
#'         silently dropped.
#' }
#'
#' @return A list with:
#' \itemize{
#'   \item \code{data.source} – inferred caller (e.g. "bcftools", "freebayes", …);
#'   \item \code{vcf.source.raw} – exact raw \code{##source=} string or \code{NA};
#'   \item \code{stacks.check} – whether this is a Stacks VCF (TRUE/FALSE);
#'   \item \code{check.header} – cleaned \code{SeqArray::seqVCF_Header()} object;
#'   \item \code{markers.info} – final INFO fields to import (or \code{NULL} = all);
#'   \item \code{overwrite.metadata} – final FORMAT fields to import
#'         (or \code{NULL} = all).
#' }
#'
#' @rdname check_header_source_vcf
#' @keywords internal
#' @export
check_header_source_vcf <- function(
    vcf,
    markers.info = NULL,
    vcf.metadata = NULL
) {

  # read VCF header
  check.header <- SeqArray::seqVCF_Header(vcf.fn = vcf)

  # generic FORMAT fixes (freebayes/samtools/Stacks inconsistencies)
  problematic.id <- c("AD","AO","QA","GL","CATG","RO","QR","MIN_DP","PL","DPR")
  problematic.id <- purrr::keep(problematic.id, problematic.id %in% check.header$format$ID)
  for (p in problematic.id) {
    check.header$format[check.header$format$ID == p, "Number"] <- "."
  }

  # DArT / Stacks GT type fix
  probl.dart <- check.header$format[check.header$format$ID == "GT", "Type"] == "Integer"
  if (length(probl.dart) > 0 && isTRUE(probl.dart)) {
    check.header$format[check.header$format$ID == "GT", "Type"] <- "String"
  }

  # extract raw source line
  source.lines <- check.header$header$value[check.header$header$id == "source"]
  vcf.source.raw <- if (length(source.lines) > 0L) source.lines[[1]] else NA_character_
  source.lower <- tolower(vcf.source.raw %||% "")

  # helper to search any meta-information line
  header_has <- function(pattern) {
    any(grepl(pattern, check.header$header$id,    ignore.case = TRUE)) ||
      any(grepl(pattern, check.header$header$value, ignore.case = TRUE))
  }

  # initialise
  data.source  <- "unknown"
  stacks.check <- FALSE

  # Stacks
  if (grepl("stacks", source.lower) || header_has("stacks v2")) {
    data.source  <- "stacks"
    stacks.check <- TRUE
  }

  # freebayes + dirty build
  if (grepl("freebayes", source.lower)) {
    if (grepl("dirty", source.lower)) {
      data.source <- "freebayes_dirty"
    } else {
      data.source <- "freebayes"
    }
  }

  # bcftools (bcftools rarely sets source=)
  if (data.source == "unknown" &&
      (header_has("bcftoolsVersion") ||
       header_has("bcftools_callVersion") ||
       header_has("bcftools_viewVersion") ||
       header_has("bcftoolsCommand"))) {
    data.source <- "bcftools"
  }

  # GATK
  if (data.source == "unknown" &&
      (header_has("GATKCommandLine") ||
       grepl("gatk", source.lower) ||
       grepl("haplotypecaller", source.lower))) {
    data.source <- "gatk"
  }

  # DeepVariant
  if (data.source == "unknown" &&
      (header_has("DeepVariant") || grepl("deepvariant", source.lower))) {
    data.source <- "deepvariant"
  }

  # DRAGEN
  if (data.source == "unknown" &&
      (header_has("DRAGEN") || grepl("dragen", source.lower))) {
    data.source <- "dragen"
  }

  # GLnexus
  if (data.source == "unknown" &&
      (header_has("GLnexus") || grepl("glnexus", source.lower))) {
    data.source <- "glnexus"
  }

  # vcftools
  if (data.source == "unknown" && header_has("^vcftools_")) {
    data.source <- "vcftools"
  }

  # samtools
  if (data.source == "unknown" &&
      (header_has("samtoolsVersion") || header_has("samtoolsCommand"))) {
    data.source <- "samtools"
  }

  # special handling for freebayes_dirty
  if (identical(data.source, "freebayes_dirty")) {
    message(
      "\nDetected freeBayes 'dirty' build — restricting FORMAT to GT/DP ",
      "and setting problematic INFO Number fields to '.'.\n"
    )

    # retain only GT + DP
    check.header$format <- dplyr::filter(
      check.header$format,
      .data$ID %in% c("GT", "DP")
    )

    # problematic INFO tags known from dirty builds
    prob.info <- c(
      "AN","DP","DPB","EPPR","GTI","MQMR","NS","NUMALT","ODDS","PAIREDR",
      "PQR","PRO","QR","RO","RPPR","SRF","SRP","SRR"
    )
    prob.info <- purrr::keep(prob.info, prob.info %in% check.header$info$ID)
    for (p in prob.info) {
      check.header$info[check.header$info$ID == p, "Number"] <- "."
    }
  }

  # caller-aware suggestions for INFO / FORMAT import

  info.ids <- check.header$info$ID
  fmt.ids  <- check.header$format$ID

  default_info <- NULL
  default_fmt  <- NULL

  if (identical(data.source, "stacks")) {
    default_fmt  <- c("GT", "DP", "AD", "GL", "GQ")
    default_info <- c("NS", "AF", "MAF", "HWE", "FIS", "DP")

  } else if (identical(data.source, "freebayes") || identical(data.source, "freebayes_dirty")) {
    default_fmt  <- c("GT", "DP", "AD", "GQ", "GL")
    default_info <- c(
      "DP", "AC", "AN", "NS",
      "AB", "ABP", "AO", "RO", "QR", "QA",
      "SRF", "SRR", "SRP", "SAF", "SAR", "SAP",
      "ODDS"
    )

  } else if (identical(data.source, "bcftools")) {
    default_fmt  <- c("GT", "AD", "DP", "PL", "GQ")
    default_info <- c(
      "DP", "AC", "AN", "NS",
      "MQ", "MQ0F", "DP4", "VDB", "F_MISSING"
    )

  } else if (identical(data.source, "gatk")) {
    default_fmt  <- c("GT", "AD", "DP", "GQ", "PL")
    default_info <- c(
      "DP", "AC", "AN", "MQ", "QD", "FS", "SOR",
      "MQRankSum", "ReadPosRankSum"
    )

  } else if (identical(data.source, "deepvariant")) {
    default_fmt  <- c("GT", "DP", "GQ", "AD")
    default_info <- c("DP", "AC", "AN", "VAF")

  } else if (identical(data.source, "samtools")) {
    default_fmt  <- c("GT", "DP", "GQ")
    default_info <- c("DP", "AC", "AN", "MQ", "MQ0F", "DP4")

  } else if (identical(data.source, "vcftools")) {
    default_fmt  <- NULL
    default_info <- c("DP", "AC", "AN", "NS")

  } else if (identical(data.source, "glnexus") || identical(data.source, "dragen")) {
    default_fmt  <- c("GT", "DP", "GQ", "AD")
    default_info <- c("DP", "AC", "AN", "MQ", "QD")

  } else {
    # unknown/other → import everything
    default_fmt  <- NULL
    default_info <- NULL
  }

  # INFO fields import control (markers.info)
  if (is.null(markers.info)) {
    markers.info.import <- if (is.null(default_info)) NULL else intersect(default_info, info.ids)
  } else {
    markers.info.import <- intersect(markers.info, info.ids)
  }
  if (!is.null(markers.info.import) && !length(markers.info.import)) {
    markers.info.import <- NULL
  }

  # FORMAT fields import control (vcf.metadata / overwrite.metadata)
  if (is.logical(vcf.metadata)) {
    if (vcf.metadata) {
      overwrite.metadata <- NULL       # import all FORMAT fields
    } else {
      overwrite.metadata <- "GT"       # only genotype
    }
  } else {
    if (is.null(vcf.metadata)) {
      overwrite.metadata <- default_fmt
    } else {
      overwrite.metadata <- vcf.metadata
      if (!"GT" %in% overwrite.metadata) {
        message("GT field always included in vcf.metadata")
        overwrite.metadata <- c("GT", overwrite.metadata)
      }
      vcf.metadata <- TRUE
    }
  }

  if (!is.null(overwrite.metadata)) {
    overwrite.metadata <- intersect(overwrite.metadata, fmt.ids)
    if (!length(overwrite.metadata)) {
      overwrite.metadata <- NULL
    }
  }

  list(
    data.source        = data.source,
    vcf.source.raw     = vcf.source.raw,
    stacks.check       = stacks.check,
    check.header       = check.header,
    markers.info       = markers.info.import,
    overwrite.metadata = overwrite.metadata
  )
} # END check_header_source_vcf


# indexing_vcf -------------------------------------------------------
# detect_indexing------------------------------------------------------
#' @title Check and ensure that a VCF is bgzipped and indexed
#'
#' @description
#' Internal helpers used by \code{read_vcf()} to work with VCF files that
#' are suitable for parallel reading with \pkg{SeqArray}, i.e. bgzip-
#' compressed and Tabix-indexed.
#'
#' \code{detect_indexing()} checks whether a VCF file:
#' \itemize{
#'   \item is bgzip-compressed (\code{.vcf.gz});
#'   \item has an existing Tabix index (\code{.tbi} or \code{.csi});
#'   \item and that the index can be opened without error.
#' }
#'
#' \code{indexing_vcf()} ensures that a VCF file is bgzip-compressed and
#' indexed with Tabix. If the file is not \code{.vcf.gz}, it is compressed
#' with \code{Rsamtools::bgzip()}; if no index is found, it is created via
#' \code{Rsamtools::indexTabix()}.
#'
#' The typical workflow is to first call \code{detect_indexing()} and, if it
#' returns \code{TRUE}, to run \code{indexing_vcf()} to make the VCF
#' compatible with parallel chunked reading.
#'
#' @details
#' These functions are designed for internal use and are not exported.
#' \code{indexing_vcf()} emits user-friendly progress messages via \pkg{cli}
#' when compression or indexing is required, but remains quiet when files
#' are already compliant.
#'
#' @param vcf (character) Path to a VCF file. Can be plain text (\code{.vcf})
#'   or bgzipped (\code{.vcf.gz}).
#'
#' @return
#' \itemize{
#'   \item \code{detect_indexing()}: (logical) \code{TRUE} if the VCF needs
#'     compression and/or indexing (i.e. is not ready for parallel reading),
#'     \code{FALSE} if it is already bgzipped and indexed with a readable
#'     Tabix index.
#'   \item \code{indexing_vcf()}: (character) The path to the bgzipped and
#'     Tabix-indexed VCF.
#' }
#'
#' @rdname indexing_vcf
#' @keywords internal
#' @export

#' @examples
#' \dontrun{
#' if (detect_indexing(vcf)) {
#'   vcf <- indexing_vcf(vcf)
#' }
#' }
detect_indexing <- function(vcf) {
  # Not bgzipped => needs compression + index
  if (!grepl("\\.vcf\\.gz$", vcf, ignore.case = TRUE)) {
    return(TRUE)
  }

  tbi <- paste0(vcf, ".tbi")
  csi <- sub("\\.vcf\\.gz$", ".csi", vcf, ignore.case = TRUE)

  # No index file exists => needs index
  if (!file.exists(tbi) && !file.exists(csi)) {
    return(TRUE)
  }

  # Deeper test: is index readable?
  out <- tryCatch({
    suppressWarnings({
      rs <- Rsamtools::TabixFile(vcf)
      Rsamtools::open.TabixFile(rs)
      Rsamtools::close.TabixFile(rs)
    })
    FALSE  # index is fine
  }, error = function(e) TRUE) # index exists but broken

  out
} # End detect_indexing


# indexing_vcf---------------------------------------------------------
#' @rdname indexing_vcf
#' @keywords internal
#' @export
indexing_vcf <- function(vcf, bcftools.path = "bcftools", verbose = TRUE) {

  # bcftools is needed for compression, pruning and merge/index
  radiator::bcftools_require(bcftools.path = bcftools.path)

  bcftools.path <- Sys.which(bcftools.path)
  if (!nzchar(bcftools.path)) rlang::abort("bcftools not found on PATH.")
  bcftools.path <- normalizePath(path = bcftools.path, winslash = "/", mustWork = TRUE)
  vcf <- normalizePath(vcf, winslash = "/", mustWork = TRUE)

  has.gz  <- grepl("\\.vcf\\.gz$", vcf, ignore.case = TRUE)
  has.vcf <- grepl("\\.vcf$", vcf, ignore.case = TRUE)

  if (has.gz) {
    tmp.sorted <- sub("\\.vcf\\.gz$", ".sorted.vcf.gz", vcf, ignore.case = TRUE)
  } else if (has.vcf) {
    tmp.sorted <- sub("\\.vcf$", ".sorted.vcf.gz", vcf, ignore.case = TRUE)
  } else {
    rlang::abort("VCF filename must end with .vcf or .vcf.gz")
  }

  # sort VCF
  if (verbose) message("Sorting ", basename(vcf))

  sort.log   <- paste0(vcf, ".sort.log")

  res.sort <- radiator::bcftools_exec(
    bcftools.path = bcftools.path,
    args = c("sort", "-Oz", "-o", tmp.sorted, vcf),
    log.file = sort.log,
    label = "bcftools sort",
    verbose = FALSE,
    fail_on_status = FALSE
  )
  # check sorting
  if (!identical(res.sort$status, 0L)) {
    rlang::abort(
      paste0("bcftools sort failed on: ", basename(tmp.sorted), "\nSee log: ", basename(sort.log))
    )
  }


  # Final index (optional)

  # index ---------------------------------------------------

  if (verbose) message("Indexing VCF with tabix for parallel reading")
  index.log <- paste0(vcf, ".index.log")
  res.index.final <- radiator::bcftools_exec(
    bcftools.path  = bcftools.path,
    args           = c("index", "-t", tmp.sorted),
    log.file       = index.log,
    label          = "index_final",
    verbose        = FALSE,
    fail_on_status = FALSE
  )
  if (!identical(res.index.final$status, 0L)) {
    warning("bcftools index failed on final VCF; see log: ", basename(index.log), call. = FALSE)
  } else if (verbose) {
    if (verbose) message("Final VCF indexed with tabix: ", basename(tmp.sorted), ".tbi")
  }
  tmp.sorted
} # End indexing_vcf


# clean_vcf_markers_meta--------------------------------------------------------
#' @title Clean and normalise markers metadata from a VCF import
#'
#' @description
#' Internal helper used by \code{read_vcf()} to harmonise the markers metadata
#' (\code{markers.meta}) coming from different VCF-producing pipelines
#' (Stacks, ipyrad, PLINK, DArT VCF, FreeBayes, GATK, etc.).
#'
#' The function:
#' \itemize{
#'   \item normalises \code{CHROM}, \code{LOCUS}, \code{POS}, \code{COL},
#'     and (optionally) \code{STRANDS};
#'   \item detects and adjusts special cases (Stacks, ipyrad, PLINK, DArT VCF);
#'   \item generates a unique \code{MARKERS} identifier and a sequential
#'     \code{M_SEQ} index;
#'   \item attaches reference and alternate alleles from the GDS.
#' }
#'
#' It is designed for internal use in \pkg{radiator} and is not exported.
#'
#' @param markers.meta A tibble/data.frame of markers metadata containing at
#'   least \code{VARIANT_ID}, \code{CHROM}, \code{LOCUS}, \code{POS}.
#' @param gds A SeqArray GDS object with the VCF data already imported.
#' @param data.source (character) A short tag describing the originating
#'   pipeline (e.g. \code{"Stacks"}, \code{"PLINK"}, \code{"ipyrad"},
#'   \code{"freeBayes"}). May be updated inside the function (e.g. to
#'   \code{"dart.vcf"} if the structure matches DArT VCF).
#' @param ref.genome (logical) Output of \code{detect_ref_genome()}, indicating
#'   whether the data were produced with a reference genome.
#' @param stacks.checks (logical) Whether to apply additional Stacks-specific
#'   adjustments when \code{ref.genome = FALSE}.
#' @param verbose (logical) Display messages during cleaning.
#'   Default: \code{verbose = TRUE}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{markers.meta}: cleaned tibble with harmonised fields;
#'   \item \code{data.source}: (possibly updated) data source tag;
#'   \item \code{detect.strand}: logical flag indicating whether strand info
#'     was detected from the LOCUS field.
#' }
#'
#' @keywords internal
#' @export
clean_vcf_markers_meta <- function(
    markers.meta,
    gds,
    data.source,
    ref.genome,
    stacks.checks = FALSE,
    verbose = TRUE
) {

  if (verbose) cli::cli_alert_success("Cleaning and normalizing markers metadata")

  # Stacks specific adjustments -------------------------------------------------
  if (!ref.genome) {
    if (stacks.checks) {
      markers.meta %<>%
        dplyr::mutate(
          LOCUS = CHROM,
          CHROM = "1",
          COL   = as.integer(POS) - 1L
        )
    } else {
      markers.meta %<>% dplyr::mutate(CHROM = "1")
    }
  }

  # ipyrad VCF ------------------------------------------------------------------
  if (stringi::stri_detect_fixed(str = data.source, pattern = "ipyrad")) {
    markers.meta %<>%
      dplyr::mutate(
        LOCUS = stringi::stri_replace_all_fixed(
          str           = CHROM,
          pattern       = "locus_",
          replacement   = "",
          vectorize_all = FALSE
        ),
        COL = as.integer(POS)
      )
  }

  # PLINK VCF -------------------------------------------------------------------
  if (stringi::stri_detect_fixed(str = data.source, pattern = "PLINK")) {
    markers.meta %<>%
      dplyr::mutate(
        LOCUS = as.integer(factor(x = CHROM)),
        COL   = 1L
      )
  }

  # GATK, platypus, freebayes, samtools: weird LOCUS handling ------------------
  weird.locus <- length(unique(markers.meta$LOCUS)) <= 1L
  if (weird.locus && !stacks.checks) {
    if (verbose) cli::cli_alert_info("LOCUS field empty... adding unique id instead")
    markers.meta$LOCUS <- markers.meta$VARIANT_ID
  }

  # VCF: LOCUS cleaning and strands detection ----------------------------------
  # Use LOCUS by name instead of column index to avoid fragile [1,3] usage.
  locus.first <- markers.meta$LOCUS[1]

  if (!is.na(locus.first) && isTRUE(unique(stringi::stri_detect_fixed(
    str     = locus.first,
    pattern = c("|", "-", ":", ">")
  )))) {
    data.source <- "dart.vcf"
    update_radiator_gds(gds = gds, node.name = "data.source", value = data.source)
  }

  if (data.source != "dart.vcf" &&
      !is.na(locus.first) &&
      stringi::stri_detect_regex(str = locus.first, pattern = "[^[:alnum:]]+")) {

    locus.sep <- unique(
      stringi::stri_extract_all_regex(
        str            = locus.first,
        pattern        = "[^a-zA-Z0-9-+-]",
        omit_no_match  = TRUE
      )[[1]]
    )

    markers.meta <- suppressWarnings(
      tidyr::separate(
        data   = markers.meta,
        col    = LOCUS,
        into   = c("LOCUS", "COL", "STRANDS"),
        sep    = locus.sep,
        extra  = "drop",
        fill   = "warn",
        remove = TRUE,
        convert = TRUE
      )
    )

    if (anyNA(markers.meta$STRANDS)) {
      markers.meta$STRANDS <- NA_character_
    }
    locus.sep <- NULL

    detect.strand <- any(stringi::stri_detect_fixed(
      str     = unique(markers.meta$STRANDS),
      pattern = "+"
    ))
    if (anyNA(detect.strand)) detect.strand <- FALSE

  } else {
    detect.strand <- FALSE
  }

  # VCF: DArT LOCUS cleaning ---------------------------------------------------
  if (identical(data.source, "dart.vcf")) {
    markers.meta %<>%
      dplyr::mutate(
        CHROM = stringi::stri_replace_all_fixed(
          str           = CHROM,
          pattern       = ".",
          replacement   = "denovo",
          vectorize_all = FALSE
        ),
        POS = stringi::stri_replace_na(str = POS, replacement = "50"),
        COL = stringi::stri_extract_first_regex(
          str     = LOCUS,
          pattern = "[-][0-9]+[\\:]"
        ),
        COL = stringi::stri_replace_all_fixed(
          str           = COL,
          pattern       = c("-", ":"),
          replacement   = c("", ""),
          vectorize_all = FALSE
        ),
        COL   = as.integer(COL),
        LOCUS = stringi::stri_extract_first_regex(
          str     = LOCUS,
          pattern = "^[0-9]+"
        ),
        LOCUS = stringi::stri_join(LOCUS, POS, sep = "_"),
        POS   = COL
      )
  }

  # Generate MARKERS column and fix types --------------------------------------
  markers.meta %<>%
    dplyr::mutate(
      dplyr::across(
        .cols = c(CHROM, LOCUS, POS),
        .fns  = radiator::clean_markers_names
      )
    ) %>%
    dplyr::mutate(
      MARKERS = stringi::stri_join(CHROM, LOCUS, POS, sep = "__"),
      REF     = SeqArray::seqGetData(gdsfile = gds, var.name = "$ref"),
      ALT     = SeqArray::seqGetData(gdsfile = gds, var.name = "$alt")
    )

  # PLINK file with duplicate MARKERS (ISOFORM-like tags) ----------------------
  dup.markers <- length(markers.meta$MARKERS) - length(unique(markers.meta$MARKERS))
  if (dup.markers > 0L) {
    cli::cli_alert_info("\nNumber of duplicate MARKERS id: {dup.markers}")
    cli::cli_alert_info("Adding integer to differentiate")

    markers.meta %<>%
      dplyr::arrange(MARKERS) %>%
      dplyr::mutate(MARKERS_NEW = MARKERS) %>%
      dplyr::group_by(MARKERS_NEW) %>%
      dplyr::mutate(
        MARKERS = stringi::stri_join(MARKERS, seq_len(dplyr::n()), sep = "_")
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(-MARKERS_NEW)
  }

  # strip markers meta, here we could think of using variant id but let's wait.
  markers.meta %<>%
    dplyr::mutate(
      M_SEQ = as.integer(factor(x = MARKERS, levels = unique(MARKERS)))
    )

  list(
    markers.meta   = markers.meta,
    data.source    = data.source,
    detect.strand  = detect.strand
  )
} # End clean_vcf_markers_meta

# filter_vcf_filter_column------------------------------------------------------
#' @title Filter markers based on the VCF FILTER column
#'
#' @description
#' Internal helper used by \code{read_vcf()} to synchronise radiator's
#' \code{FILTERS} field with the VCF \code{FILTER} annotation.
#'
#' If the VCF contains more than one distinct FILTER value, markers with
#' FILTER != "PASS" (e.g., low quality, failed filters, etc.) are tagged in
#' radiator as:
#' \itemize{
#'   \item \code{"vcf filter column"} in \code{FILTERS}
#' }
#'
#' The function updates both \code{markers.meta} inside the GDS and the
#' \code{filters.parameters} tracking object.
#'
#' @param gds A SeqArray GDS object containing the imported VCF.
#' @param markers.meta The current markers metadata tibble (from radiator).
#' @param filters.parameters The filters.parameters object to update.
#' @param path.folder Path to the filtering results folder.
#' @param file.date Character string used for naming result files.
#' @param verbose (logical) Verbosity.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{markers.meta}: updated tibble;
#'   \item \code{filters.parameters}: updated parameters object.
#' }
#' @name read_vcf
#' @rdname read_vcf
#' @keywords internal
#' @export
filter_vcf_filter_column <- function(
    gds,
    markers.meta,
    filters.parameters,
    path.folder,
    file.date,
    verbose = TRUE
) {

  # Extract FILTER column from VCF ---------------------------------------------
  filter.vcf <- tibble::tibble(
    FILTER_VCF = SeqArray::seqGetData(gds, "annotation/filter")
  ) %>%
    dplyr::bind_cols(
      markers.meta %>%
        dplyr::filter(FILTERS == "whitelist") %>%
        dplyr::select(VARIANT_ID)
    )

  filter.check.unique <- unique(filter.vcf$FILTER_VCF)

  # If all markers have same FILTER: nothing to do -----------------------------
  if (length(filter.check.unique) <= 1) {
    return(list(
      markers.meta = markers.meta,
      filters.parameters = filters.parameters
    ))
  }

  # Otherwise apply VCF FILTER --------------------------------------------------
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

  # Save updated markers.meta into GDS -----------------------------------------
  update_radiator_gds(
    gds       = gds,
    node.name = "markers.meta",
    value     = markers.meta,
    sync      = TRUE
  )

  # Update filters.parameters ---------------------------------------------------
  filters.parameters <- radiator_parameters(
    generate      = FALSE,
    initiate      = FALSE,
    update        = TRUE,
    parameter.obj = filters.parameters,
    data          = gds,
    filter.name   = "vcf filter column",
    param.name    = "PASS or not",
    path.folder   = path.folder,
    file.date     = file.date,
    verbose       = verbose
  )

  radiator_results_message(
    rad.message = stringi::stri_join("\nFilter vcf filter column: PASS or not"),
    filters.parameters,
    internal = TRUE,
    verbose
  )

  # Cleanup + return -----------------------------------------------------------
  return(list(
    markers.meta = markers.meta,
    filters.parameters = filters.parameters
  ))
}#END filter_vcf_filter_column

# filter_duplicated_markers_strands---------------------------------------------
#' @title Filter duplicated SNPs on different strands (+/-)
#'
#' @description
#' Internal helper used by \code{read_vcf()} to detect and optionally filter
#' duplicated SNPs that appear on different strands (e.g. + / -) at the same
#' genomic position.
#'
#' When \code{detect.strand = TRUE}, the function:
#' \itemize{
#'   \item identifies duplicated SNPs by \code{CHROM} and \code{POS};
#'   \item depending on \code{filter.strands}, either keeps both, blacklists
#'     all, or uses a "best.stats" rule based on missingness and minor allele
#'     count (MAC);
#'   \item updates \code{markers.meta\$FILTERS} for pruned markers;
#'   \item writes blacklist/whitelist TSV files;
#'   \item updates the GDS and \code{filters.parameters}.
#' }
#'
#' @param gds A SeqArray GDS object.
#' @param markers.meta Tibble of markers metadata (must include \code{VARIANT_ID},
#'   \code{MARKERS}, \code{CHROM}, \code{LOCUS}, \code{POS}, \code{FILTERS}).
#' @param filters.parameters The filters.parameters tracking object.
#' @param path.folder Path to the main filtering folder.
#' @param file.date Character string used in filenames.
#' @param filter.strands (character) Strategy to handle duplicated strands:
#'   \itemize{
#'     \item \code{"keep.both"} – detect but do not filter;
#'     \item \code{"blacklist"} – blacklist all duplicated markers;
#'     \item \code{"best.stats"} – select markers based on missingness and MAC
#'       (see details in code) and blacklist them (current behaviour).
#'   }
#' @param detect.strand (logical) Whether strand information was detected in
#'   \code{LOCUS} / \code{STRANDS}. If \code{FALSE}, the function returns
#'   without modification.
#' @param parallel.core (integer) Number of cores to use for
#'   \code{SeqArray::seqAlleleCount()} and \code{SeqArray::seqMissing()} when
#'   \code{filter.strands = "best.stats"}.
#' @param internal (logical) Passed to \code{radiator_results_message()} and
#'   \code{radiator_parameters()}.
#' @param verbose (logical) Verbosity of messages.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{markers.meta}: updated tibble;
#'   \item \code{filters.parameters}: updated parameters object.
#' }
#'
#' @keywords internal
#' @name filter_duplicated_markers_strands
#' @rdname filter_duplicated_markers_strands
#' @export
filter_duplicated_markers_strands <- function(
    gds,
    markers.meta,
    filters.parameters,
    path.folder,
    file.date,
    filter.strands = c("keep.both", "blacklist", "best.stats"),
    detect.strand = FALSE,
    parallel.core = 1L,
    internal = TRUE,
    verbose = TRUE
) {

  filter.strands <- match.arg(filter.strands)

  # If no strand information, nothing to do ------------------------------------
  if (!detect.strand) {
    return(list(
      markers.meta       = markers.meta,
      filters.parameters = filters.parameters
    ))
  }

  # Identify duplicated SNPs on different strands ------------------------------
  blacklist.strands <- markers.meta %>%
    dplyr::distinct(VARIANT_ID, MARKERS, CHROM, LOCUS, POS) %>%
    dplyr::group_by(CHROM, POS) %>%
    dplyr::mutate(n = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(n > 1L) %>%
    dplyr::select(-n) %>%
    dplyr::distinct(CHROM, POS, .keep_all = TRUE)

  if (nrow(blacklist.strands) == 0L) {
    return(list(
      markers.meta       = markers.meta,
      filters.parameters = filters.parameters
    ))
  }

  if (verbose) {
    cli::cli_alert_warning("Detected {nrow(blacklist.strands)} duplicate SNPs on different strands (+/-)")
    cli::cli_alert_info("    By default radiator prunes those SNPs")
    cli::cli_alert_info("    To change this behavior, use argument: filter.strands")
  }

  if (filter.strands == "keep.both") {
    message("Keeping both duplicated markers")
    return(list(
      markers.meta       = markers.meta,
      filters.parameters = filters.parameters
    ))
  }

  # filter.strands == "best.stats" ---------------------------------------------
  if (filter.strands == "best.stats") {

    # backup full variant vector
    variant.bk <- markers.meta$VARIANT_ID

    # restrict GDS to duplicated variants
    sync_gds(gds = gds, variant.id = blacklist.strands$VARIANT_ID)

    blacklist.strands <- SeqArray::seqAlleleCount(
      gdsfile    = gds,
      ref.allele = NULL,
      parallel   = parallel.core
    ) %>%
      unlist() %>%
      matrix(
        data    = .,
        nrow    = nrow(blacklist.strands),
        ncol    = 2,
        byrow   = TRUE,
        dimnames = list(
          rownames = blacklist.strands$VARIANT_ID,
          colnames = c("REF_COUNT", "ALT_COUNT")
        )
      ) %>%
      tibble::as_tibble(rownames = "VARIANT_ID") %>%
      tibble::add_column(
        .data  = .,
        MARKERS = blacklist.strands$MARKERS,
        CHROM   = blacklist.strands$CHROM,
        POS     = blacklist.strands$POS,
        MISSING_PROP = SeqArray::seqMissing(
          gdsfile    = gds,
          per.variant = TRUE,
          parallel   = parallel.core
        ),
        .after = 1
      ) %>%
      dplyr::mutate(
        MAC = dplyr::if_else(ALT_COUNT < REF_COUNT, ALT_COUNT, REF_COUNT),
        ALT_COUNT = NULL,
        REF_COUNT = NULL
      ) %>%
      dplyr::group_by(CHROM, POS) %>%
      dplyr::filter(MISSING_PROP == max(MISSING_PROP)) %>%  # current behaviour
      dplyr::filter(MAC == min(MAC)) %>%                    # current behaviour
      dplyr::ungroup() %>%
      dplyr::arrange(CHROM, POS) %>%
      dplyr::distinct(CHROM, POS, .keep_all = TRUE) %>%
      dplyr::select(MARKERS)

    # restore original variant selection in GDS
    sync_gds(gds = gds, variant.id = variant.bk, verbose = FALSE)
  }

  # filter.strands == "blacklist" ----------------------------------------------
  if (filter.strands == "blacklist") {
    blacklist.strands %<>% dplyr::distinct(MARKERS)
  }

  # If we’re here, we are blacklisting something -------------------------------
  if (nrow(blacklist.strands) == 0L) {
    return(list(
      markers.meta       = markers.meta,
      filters.parameters = filters.parameters
    ))
  }

  # Create folder for blacklist/whitelist --------------------------------------
  path.folder.strands <- generate_folder(
    rad.folder  = "filter_duplicated_markers",
    path.folder = path.folder,
    internal    = FALSE,
    file.date   = file.date,
    verbose     = verbose
  )

  # Save blacklist -------------------------------------------------------------
  blacklist.strands %<>%
    dplyr::mutate(FILTER = "filter.strands")

  write_radiator_tsv(
    data          = blacklist.strands,
    path.folder   = path.folder.strands,
    filename      = "blacklist.duplicated.markers.strands",
    date          = TRUE,
    internal      = FALSE,
    write.message = "standard",
    verbose       = FALSE
  )

  # Update markers.meta --------------------------------------------------------
  markers.meta %<>%
    dplyr::mutate(
      FILTERS = dplyr::if_else(
        MARKERS %in% blacklist.strands$MARKERS,
        "filter.strands",
        FILTERS
      )
    )

  write_radiator_tsv(
    data          = dplyr::filter(markers.meta, FILTERS == "whitelist"),
    path.folder   = path.folder.strands,
    filename      = "whitelist.duplicated.markers.strands",
    date          = TRUE,
    internal      = FALSE,
    write.message = "standard",
    verbose       = FALSE
  )

  # Update GDS ------------------------------------------------------------------
  update_radiator_gds(
    gds       = gds,
    node.name = "markers.meta",
    value     = markers.meta,
    sync      = TRUE
  )

  # Update filters.parameters ---------------------------------------------------
  filters.parameters <- radiator_parameters(
    generate      = FALSE,
    initiate      = FALSE,
    update        = TRUE,
    parameter.obj = filters.parameters,
    data          = gds,
    filter.name   = "filter duplicated markers on different strands",
    param.name    = "filter.strands",
    values        = filter.strands,
    path.folder   = path.folder,
    file.date     = file.date,
    internal      = internal,
    verbose       = verbose
  )

  # Results message ------------------------------------------------------------
  radiator_results_message(
    rad.message = stringi::stri_join(
      "Filter duplicated markers on different strands: ",
      filter.strands
    ),
    filters.parameters,
    internal,
    verbose
  )

  list(
    markers.meta       = markers.meta,
    filters.parameters = filters.parameters
  )
} #END filter_duplicated_markers_strands

# filter_haplotype_format ------------------------------------------------------
#' @title Filter non-biallelic / haplotype-style variants in a radiator GDS
#'
#' @description
#' Internal helper used by \code{read_vcf()} and related workflows to detect and
#' blacklist variants that are not simple biallelic SNPs.
#'
#' The function:
#' \itemize{
#'   \item re-imports \code{markers.meta} from the radiator node;
#'   \item identifies clean biallelic SNPs as:
#'     \itemize{
#'       \item \code{REF_LEN == 1},
#'       \item \code{ALT_LEN == 1},
#'       \item \code{ALT} has no comma (single ALT allele);
#'     }
#'   \item blacklists all variants that are \emph{not} biallelic SNPs
#'     (haplotypes, INDELs, multi-ALT sites, etc.);
#'   \item updates \code{markers.meta$FILTERS}, the GDS, a blacklist file, and a
#'     whitelist of remaining biallelic markers;
#'   \item updates \code{filters.parameters} with the \code{filter_haplotype_format}
#'     step.
#' }
#'
#' When \code{filter.haplotype.format = FALSE}, the function only re-imports
#' and returns \code{markers.meta} and leaves \code{filters.parameters} untouched.
#'
#' @param gds
#' A radiator GDS object (\code{SeqVarGDSClass}) with a \code{/radiator/markers.meta}
#' node containing at least \code{MARKERS}, \code{REF} and \code{ALT}.
#'
#' @param filters.parameters
#' The current filter-parameter object generated by \code{radiator_parameters()}.
#'
#' @param path.folder (character)
#' Root folder used to store filter outputs.
#'
#' @param file.date (character)
#' Timestamp string (e.g. \code{"20251203@1919"}) used in filenames.
#'
#' @param filter.haplotype.format (logical)
#' Activate or deactivate this filter.
#' Default: \code{filter.haplotype.format = FALSE}.
#'
#' @param internal (logical)
#' Indicate if the function is called internally by other radiator functions,
#' used for controlling side-effect messages and paths.
#' Default: \code{internal = TRUE}.
#'
#' @param verbose (logical)
#' When \code{TRUE}, the function prints messages about detected haplotypes and
#' written files.
#' Default: \code{verbose = TRUE}.
#'
#' @return
#' A list with:
#' \itemize{
#'   \item \code{markers.meta} – updated markers metadata tibble;
#'   \item \code{filters.parameters} – updated filter-parameter object.
#' }
#'
#' @name filter_haplotype_format
#' @rdname filter_haplotype_format
#' @export
filter_haplotype_format <- function(
    gds,
    filters.parameters,
    path.folder,
    file.date,
    filter.haplotype.format = FALSE,
    internal = TRUE,
    verbose = TRUE
) {

  # If filter is disabled, just return current markers.meta + filters.parameters
  if (!isTRUE(filter.haplotype.format)) {
    markers.meta <- radiator::extract_markers_metadata(gds = gds)
    return(list(
      markers.meta       = markers.meta,
      filters.parameters = filters.parameters
    ))
  }

  # Re-import full markers.meta from radiator node -----------------------------
  markers.meta <- radiator::extract_markers_metadata(gds = gds)

  # Identify non-biallelic / haplotype-style variants -------------------------
  # Logic for "true" biallelic SNP:
  #   - REF length == 1
  #   - ALT length == 1
  #   - ALT has no comma (single ALT allele)
  # Everything else is blacklisted.
  blacklist.multiallelic <- radiator::extract_markers_metadata(
    gds       = gds,
    whitelist = TRUE
  ) %>%
    dplyr::mutate(
      REF_LEN   = stringi::stri_length(str = REF),
      ALT_LEN   = stringi::stri_length(str = ALT),
      ALT_MULTI = stringi::stri_detect_fixed(ALT, ","), # multi-allelic ALT
      IS_SNP    = REF_LEN == 1L & ALT_LEN == 1L & !ALT_MULTI
    ) %>%
    dplyr::filter(!IS_SNP)

  n.snp.m <- nrow(blacklist.multiallelic)

  if (n.snp.m > 0L) {

    if (verbose) {
      message("\nDetected ", n.snp.m,
              " non-biallelic markers (haplotypes / multi-ALT / INDELs)")
      message("blacklisting non-biallelic variants")
    }

    # Folders ------------------------------------------------------------------
    path.folder.biallelic <- generate_folder(
      rad.folder  = "filter_haplotype_format",
      path.folder = path.folder,
      internal    = FALSE,
      file.date   = file.date,
      verbose     = verbose
    )

    # Blacklist table ----------------------------------------------------------
    blacklist.multiallelic %<>%
      dplyr::mutate(FILTER = "filter.haplotype.format")

    write_radiator_tsv(
      data          = blacklist.multiallelic,
      path.folder   = path.folder.biallelic,
      filename      = "blacklist.haplotype.markers",
      date          = TRUE,
      internal      = FALSE,
      write.message = "standard",
      verbose       = verbose
    )

    # Update markers.meta FILTERS ----------------------------------------------
    markers.meta %<>%
      dplyr::mutate(
        FILTERS = dplyr::if_else(
          MARKERS %in% blacklist.multiallelic$MARKERS,
          "filter.haplotype.format",
          FILTERS
        )
      )

    # Whitelist of remaining biallelic markers --------------------------------
    write_radiator_tsv(
      data          = dplyr::filter(markers.meta, FILTERS == "whitelist"),
      path.folder   = path.folder.biallelic,
      filename      = "whitelist.biallelic.markers",
      date          = TRUE,
      internal      = FALSE,
      write.message = "standard",
      verbose       = verbose
    )

    # Update GDS ---------------------------------------------------------------
    update_radiator_gds(
      gds       = gds,
      node.name = "markers.meta",
      value     = markers.meta,
      sync      = TRUE
    )

    # Update filters.parameters ------------------------------------------------
    filters.parameters <- radiator_parameters(
      generate      = FALSE,
      initiate      = FALSE,
      update        = TRUE,
      parameter.obj = filters.parameters,
      data          = gds,
      filter.name   = "filter_haplotype_format",
      param.name    = "filter.haplotype.format",
      values        = "TRUE",
      path.folder   = path.folder,
      file.date     = file.date,
      internal      = internal,
      verbose       = verbose
    )

    # Results message ----------------------------------------------------------
    radiator_results_message(
      rad.message = stringi::stri_join(
        "\nFilter haplotype format: ",
        filter.haplotype.format
      ),
      filters.parameters,
      internal,
      verbose
    )
  } else {
    if (verbose) {
      message("\nNo non-biallelic / haplotype-style markers detected")
    }
  }

  # Return updated markers.meta and filters.parameters -------------------------
  list(
    markers.meta       = markers.meta,
    filters.parameters = filters.parameters
  )
} # End filter_haplotype_format


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

#' @keywords internal
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
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
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
  split_vcf <- carrier::crate(function(data, filename) {
    split.id <- unique(data$SPLIT)
    filename <- stringi::stri_join(filename, "_", split.id)
    radiator::write_vcf(data = dplyr::select(data, -SPLIT),
                        strata = FALSE, filename = filename)
  })#End split_vcf


  message("Importing and tidying the vcf...")
  input <- suppressMessages(
    radiator::tidy_genomic_data(
      data = data,
      strata = strata,
      vcf.metadata = FALSE,
      whitelist.markers = whitelist.markers,
      filename = NULL,
      verbose = FALSE) %>%
      dplyr::full_join(split, by = "INDIVIDUALS")
  )

  split <- strata <- blacklist <- NULL
  radiator_future(
    .x = input,
    .f = split_vcf,
    split.vec = FALSE,
    split.with = "SPLIT",
    flat.future = "walk",
    parallel.core = parallel.core,
    filename = filename
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
#   # if (!is.null(filter.ma)) {
#   #   input <- filter_maf(
#   #     data = input,
#   #     interactive.filter = FALSE,
#   #     filter.ma = filter.ma,
#   #     parallel.core = parallel.core,
#   #     verbose = FALSE)$tidy.filtered.mac
#   # }
#
#   # Write VCF in the working directory
#   radiator::write_vcf(data = input, strata = FALSE, filename = filename)
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
  input  %<>%
    radiator::rad_long(
      x = .,
      cols = c("#CHROM", "POS", "ID",  "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"),
      names_to = "INDIVIDUALS",
      values_to = "FORMAT_ID",
      variable_factor = FALSE
    ) %>%
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
    radiator::rad_wide(
      x = .,
      formula = "`#CHROM` + POS + ID +  REF + ALT + QUAL + FILTER + INFO + FORMAT ~ INDIVIDUALS",
      values_from = "FORMAT_ID") %>%
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
  # write_delim(x = data_frame("##fileformat=VCFv4.2"), file = filename, delim = " ", append = FALSE, col_names = FALSE)
  vcf.header[1,] <- "##fileformat=VCFv4.3"

  # File date ------------------------------------------------------------------
  file.date <- stringi::stri_replace_all_fixed(Sys.Date(), pattern = "-", replacement = "")
  file.date <- stringi::stri_join("##fileDate=", file.date, sep = "")
  # write_delim(x = data_frame(file.date), file = filename, delim = " ", append = TRUE, col_names = FALSE)
  vcf.header[2,] <- file.date

  # Source ---------------------------------------------------------------------
  # write_delim(x = data_frame(stringi::stri_join("##source=radiator_v.", utils::packageVersion("radiator"))), file = filename, delim = " ", append = TRUE, col_names = FALSE)
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
  suppressWarnings(readr::write_tsv(x = input, file = filename, append = TRUE, col_names = TRUE))

  cat("############################## completed ##############################\n")
}#vcf_strata
