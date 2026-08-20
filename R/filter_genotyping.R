# Filter markers genotyping ----------------------------------------------------
#' @name filter_genotyping
#' @title Filter markers based on genotyping / missing rate
#'
#' @description
#' Filter markers based on their genotyping (call) rate, i.e. the proportion of
#' non-missing genotypes per marker. This is a convenient way to remove SNPs
#' with too much missing data before downstream analyses.
#'
#' \strong{Filter targets}: SNPs.
#'
#' \strong{Statistics}: marker-level missingness
#' (\code{MISSING_PROP} from \code{\link{generate_stats}}).
#'
#' @param data
#' A Genomic Data Structure (GDS) file or connection
#' (\code{SeqVarGDSClass} or \code{gds.file}) created by
#' \code{\link{read_vcf}} or other \pkg{radiator} import functions.
#'
#' @param interactive.filter (optional, logical).
#' Do you want the filtering session to be interactive ?
#' With \code{interactive.filter = TRUE}, helper tables and plots are shown,
#' and the user is prompted for a threshold.
#' Default: \code{interactive.filter = TRUE}.
#'
#' @param filter.genotyping (optional, character or double).
#' Two modes are available:
#' \itemize{
#'   \item character string:
#'   \code{filter.genotyping = "outliers"} uses the upper outlier value from the
#'   missingness boxplot as threshold (i.e. removes high-missingness outliers);
#'   \item double:
#'   \code{filter.genotyping = 0.2} allows up to \code{0.2} missing genotypes
#'   (i.e. keeps markers with \code{MISSING_PROP <= 0.2}).
#' }
#' Default: \code{filter.genotyping = NULL}.
#'
#' If \code{interactive.filter = FALSE} and \code{filter.genotyping = NULL},
#' the function returns \code{data} unchanged.
#'
#' @param filename (optional, character).
#' Basename for some output files (plots, tables). If \code{NULL}, a name is
#' generated automatically from the function and date.
#' Default: \code{filename = NULL}.
#'
#' @param parallel.core (optional, integer).
#' Number of CPU cores to use in helper statistics.
#' Default: \code{parallel.core = parallel::detectCores() - 1}.
#'
#' @param verbose (optional, logical).
#' Show messages and progress.
#' Default: \code{verbose = TRUE}.
#'
#' @inheritParams radiator_common_arguments
#'
#' @details
#' Internally, the function uses \code{\link{generate_stats}} to compute
#' marker-level missingness (\code{MISSING_PROP}). Markers with missing
#' proportion greater than the chosen threshold are blacklisted and have their
#' \code{FILTERS} column updated to \code{"filter.genotyping"} in
#' \code{markers.meta}.
#'
#' With \code{interactive.filter = TRUE}, a helper table and plot are generated
#' to assist in choosing a threshold:
#' \itemize{
#'   \item \code{genotyping.helper.table.tsv}:
#'   number of markers kept/removed for thresholds from 0 to 1 by 0.1;
#'   \item if strata are present, \code{markers.pop.missing.helper.table.tsv}
#'   summarises missingness per population;
#'   \item a PDF/PNG helper plot summarising these patterns.
#' }
#'
#' This function operates on the GDS representation of the data and is meant
#' for genotyping/missingness-based pruning \emph{after} VCF import
#' (e.g. after \code{\link{read_vcf}}), in contrast to VCF-level slimming using
#' tools like \code{\link{filter_monomorphic_vcf}}.
#'
#' @return
#' The filtered dataset, of the same type as the input:
#' \itemize{
#'   \item if \code{data} is a GDS connection (\code{SeqVarGDSClass}), the
#'   same connection is returned with updated \code{markers.meta};
#'   \item if \code{data} is a GDS file path, the path is returned (the file is
#'   updated on disk).
#' }
#'
#' Side-effects:
#' \itemize{
#'   \item helper tables and plots are written under
#'   \code{filter_genotyping_YYYYMMDD@HHMM};
#'   \item \code{markers.meta} inside the GDS is updated and synchronised;
#'   \item the \code{filters.parameters} file is updated to record the filter.
#' }
#'
#' @seealso
#' \code{\link{filter_common_markers}},
#' \code{\link{filter_monomorphic}},
#' \code{\link{filter_ma}},
#' \code{\link{read_vcf}},
#' \code{\link{filter_rad}}.
#'
#' @examples
#' \dontrun{
#' # Filter markers with more than 20% missing genotypes
#' gds <- radiator::read_vcf(data = "populations.snps.vcf") %>%
#'        radiator::filter_genotyping(
#'             data               = .,
#'             interactive.filter = FALSE,
#'             filter.genotyping  = 0.2
#' )
#' }
#'
#' @export
#' @rdname filter_genotyping
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
filter_genotyping <- function(
    data,
    interactive.filter = TRUE,
    filter.genotyping  = NULL,
    filename           = NULL,
    parallel.core      = parallel::detectCores() - 1,
    verbose            = TRUE,
    ...
) {

  # Early exit: nothing to do ---------------------------------------------------
  if (is.null(filter.genotyping) && !interactive.filter) {
    return(data)
  }

  if (interactive.filter) verbose <- TRUE

  # Common startup --------------------------------------------------------------
  .start    <- radiator_startup(f.name = "filter_genotyping", verbose = verbose)
  file.date <- .start$file.date
  on.exit(radiator_teardown(.start), add = TRUE)

  # Function call and dotslist --------------------------------------------------
  rad.dots <- radiator_dots(
    func.name = as.list(sys.call())[[1]],
    fd        = rlang::fn_fmls_names(),
    args.list = as.list(environment()),
    dotslist  = rlang::dots_list(
      ...,
      .homonyms     = "error",
      .check_assign = TRUE
    ),
    keepers   = c("path.folder", "parameters", "force.stats", "internal"),
    verbose   = FALSE
  )

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) {
    rlang::abort("data is missing")
  }

  # Folders ---------------------------------------------------------------------
  path.folder <- generate_folder(
    rad.folder  = "filter_genotyping",
    path.folder = path.folder,
    internal    = internal,
    file.date   = file.date,
    verbose     = verbose
  )

  # Write the dots file ---------------------------------------------------------
  write_radiator_tsv(
    data          = rad.dots,
    path.folder   = path.folder,
    filename      = "radiator_filter_genotyping_args",
    date          = TRUE,
    internal      = internal,
    write.message = "Function call and arguments stored in: ",
    verbose       = verbose
  )

  # Message about steps taken during the process --------------------------------
  if (interactive.filter) {
    message("Interactive mode: on\n")
    message("Step 1. Visualization and helper table")
    message("Step 2. Filtering markers based on maximum missing proportion allowed\n")
  }

  # Detect format ---------------------------------------------------------------
  data.type <- radiator::detect_genomic_format(data)

  if (!data.type %in% c("SeqVarGDSClass", "gds.file")) {
    rlang::abort("Input not supported for this function: read function documentation")
  }

  # Import data -----------------------------------------------------------------
  if (verbose) message("Importing data ...")
  radiator_packages_dep(package = "SeqArray", cran = FALSE, bioc = TRUE)

  if (data.type == "gds.file") {
    data      <- radiator::read_rad(data, verbose = verbose)
    data.type <- "SeqVarGDSClass"
  }

  # Filter parameter file: initiate ---------------------------------------------
  filters.parameters <- radiator_parameters(
    generate      = TRUE,
    initiate      = TRUE,
    update        = FALSE,
    parameter.obj = parameters,
    data          = data,
    path.folder   = path.folder,
    file.date     = file.date,
    internal      = internal,
    verbose       = verbose
  )

  # Step 1. Visuals -------------------------------------------------------------
  if (interactive.filter) {
    message("\nStep 1. Missing visualization and helper table\n")
  }

  # Generate missingness stats --------------------------------------------------
  info <- generate_stats(
    gds             = data,
    individuals     = FALSE,
    missing         = TRUE,
    coverage        = FALSE,
    allele.coverage = FALSE,
    mac             = FALSE,
    heterozygosity  = FALSE,
    snp.per.locus   = FALSE,
    snp.position.read = FALSE,
    force.stats     = TRUE,
    plot            = FALSE,
    path.folder     = path.folder,
    file.date       = file.date,
    parallel.core   = parallel.core,
    verbose         = verbose
  )

  stats <- info$m.stats
  info  <- info$m.info

  # Helper table ----------------------------------------------------------------
  if (verbose) message("Generating missingness/genotyping helper table...")

  mis_many_markers <- function(threshold, x) {
    nrow(dplyr::filter(x, MISSING_PROP <= threshold))
  }

  n.markers <- nrow(info)

  helper.table <- tibble::tibble(MISSING_PROP = seq(0, 1, by = 0.1)) %>%
    dplyr::mutate(
      WHITELISTED_MARKERS = purrr::map_int(
        .x = MISSING_PROP,
        .f = mis_many_markers,
        x  = info
      ),
      BLACKLISTED_MARKERS = n.markers - WHITELISTED_MARKERS
    )

  readr::write_tsv(
    x    = helper.table,
    file = file.path(path.folder, "genotyping.helper.table.tsv")
  )

  # Strata-specific helper table -----------------------------------------------
  strata.meta <- extract_individuals_metadata(
    gds               = data,
    ind.field.select  = c("INDIVIDUALS", "STRATA"),
    whitelist         = TRUE
  )

  has.strata <- !is.null(strata.meta$STRATA)

  if (has.strata) {
    m.strata <- missing_per_pop(
      gds           = data,
      strata        = strata.meta,
      parallel.core = parallel.core
    )

    round_mean <- function(x) as.integer(round(mean(x, na.rm = TRUE), 0))

    mean.pop <- m.strata %>%
      dplyr::group_by(MISSING_PROP) %>%
      dplyr::summarise_if(.predicate = is.integer, .funs = round_mean) %>%
      dplyr::mutate(STRATA = "MEAN_POP")

    if (is.factor(strata.meta$STRATA)) {
      strata.pop <- levels(strata.meta$STRATA)
    } else {
      strata.pop <- unique(strata.meta$STRATA)
    }

    strata.levels <- c(strata.pop, "MEAN_POP", "OVERALL")

    suppressWarnings(
      helper.table %>%
        dplyr::mutate(STRATA = "OVERALL") %>%
        dplyr::bind_rows(mean.pop, m.strata) %>%
        dplyr::mutate(STRATA = factor(STRATA, levels = strata.levels, ordered = TRUE)) %>%
        dplyr::arrange(STRATA) %>%
        readr::write_tsv(
          x    = .,
          file = file.path(path.folder, "markers.pop.missing.helper.table.tsv")
        )
    )

    if (verbose) {
      message("File written: markers.pop.missing.helper.table.tsv")
    }

    helper.table.long <- helper.table %>%
      dplyr::mutate(STRATA = "OVERALL") %>%
      tidyr::pivot_longer(
        data      = .,
        cols      = -c("MISSING_PROP", "STRATA"),
        names_to  = "LIST",
        values_to = "MARKERS"
      ) %>%
      dplyr::mutate(STRATA = factor(STRATA, levels = strata.levels, ordered = TRUE)) %>%
      dplyr::arrange(STRATA)

  } else {

    helper.table.long <- helper.table %>%
      tidyr::pivot_longer(
        data      = .,
        cols      = -MISSING_PROP,
        names_to  = "LIST",
        values_to = "MARKERS"
      )
  }

  # Figures ---------------------------------------------------------------------
  markers.plot <- radiator_helper_plot(
    data          = helper.table.long,
    strata        = has.strata,
    stats         = "MISSING_PROP",
    x.axis.title  = "Maximum missing proportion allowed",
    x.breaks      = seq(0, 1, by = 0.1),
    plot.filename = file.path(path.folder, "markers.genotyping.helper.plot")
  )

  print(markers.plot)
  helper.table      <- NULL
  helper.table.long <- NULL
  markers.plot      <- NULL

  if (verbose) {
    message("Files written: helper tables and plots")
  }

  # Step 2. Threshold selection -------------------------------------------------
  if (interactive.filter) {
    if (verbose) {
      message("\nStep 2. Filtering markers based on maximum missing proportion\n")
    }
    filter.genotyping <- radiator_question(
      x      = "Choose the maximum missing proportion allowed: ",
      minmax = c(0, 1)
    )
  }

  # Identify threshold: outliers vs explicit value ------------------------------
  if (!purrr::is_double(filter.genotyping)) {
    out.high <- floor(
      stats$OUTLIERS_HIGH[stats$GROUP == "missing genotypes"] * 1000
    ) / 1000
    if (verbose) {
      message("\nRemoving outlier markers based on genotyping statistic: ", out.high)
    }
    filter.genotyping <- out.high
  } else {
    if (verbose) {
      message("\nRemoving markers based on genotyping statistic: ", filter.genotyping)
    }
  }

  # Whitelist and blacklist -----------------------------------------------------
  bl <- info %>%
    dplyr::filter(MISSING_PROP > filter.genotyping) %$%
    VARIANT_ID

  markers.meta <- extract_markers_metadata(gds = data, whitelist = FALSE) %>%
    dplyr::mutate(
      FILTERS = dplyr::if_else(
        VARIANT_ID %in% bl,
        "filter.genotyping",
        FILTERS
      )
    )

  # Update GDS ------------------------------------------------------------------
  update_radiator_gds(
    gds       = data,
    node.name = "markers.meta",
    value     = markers.meta,
    sync      = TRUE
  )

  write_radiator_tsv(
    data          = dplyr::filter(markers.meta, FILTERS == "filter.genotyping"),
    path.folder   = path.folder,
    filename      = "blacklist.markers.genotyping",
    date          = TRUE,
    internal      = internal,
    write.message = "standard",
    verbose       = verbose
  )

  write_radiator_tsv(
    data          = dplyr::filter(markers.meta, FILTERS == "whitelist"),
    path.folder   = path.folder,
    filename      = "whitelist.markers.genotyping",
    date          = TRUE,
    internal      = internal,
    write.message = "standard",
    verbose       = verbose
  )

  # Update parameters -----------------------------------------------------------
  filters.parameters <- radiator_parameters(
    generate      = FALSE,
    initiate      = FALSE,
    update        = TRUE,
    parameter.obj = filters.parameters,
    data          = data,
    filter.name   = "Filter genotyping",
    param.name    = "filter.genotyping",
    values        = filter.genotyping,
    path.folder   = path.folder,
    file.date     = file.date,
    internal      = internal,
    verbose       = verbose
  )

  # Results ---------------------------------------------------------------------
  radiator_results_message(
    rad.message       = stringi::stri_join(
      "\nFilter genotyping threshold: ",
      filter.genotyping
    ),
    filters.parameters = filters.parameters,
    internal          = internal,
    verbose           = verbose
  )

  return(data)
} # End filter_genotyping
