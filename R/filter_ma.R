# Minor Allele Count filtering
#' @name filter_ma
#' @title MAC, MAF and MAD filter
#' @description The Allele Frequency is the relative frequency of an allele for a markers/SNP.
#' Calculations usually involve classifying the alleles into Major/Minor,
#' Reference (REF)/Alternate (ALT) and/or Ancestral or dirived allele.
#' radiator provides a filter based of Minor Allele Count, frequency and depth function.
#' Remove/blacklist markers based on Minor/Alternate Allele Count (MAC), Frequency (MAF)
#' or Depth of coverage (MAD).
#' Use it to remove noise, sequencing errors or low marker polymorphism.
#' Some analysis performs better with the full spectrum of allele frequency,
#' consequently it's important to understand the limite of using the
#' thresholds to avoid generating biases.
#'
#' \strong{Filter targets}: Marker's alternate allele(s)
#'
#' \strong{Statistics}: count, frequency or depth of the allele per markers.
#' This filter as NO concept of strata/populations. It is computed globally.

#' @param ma.stats (optional, character) The statistic of the alternate (minor)
#' allele to keep a SNP (see details). Options are: \code{"mac", "maf", "mad"}.
#' Default: \code{ma.stats = "mac"}.

#' @param filter.ma (optional, integer or double) The threshold to
#' blacklist Minor Allele (see details). e.g. with \code{ma.stats = "mac"} and
#' \code{filter.ma = 3} will blacklist markers based on minor allele count
#' less or equal to 3.
#' This argument as no concept of locus, local or strata or pop.
#' It's applied globally, by SNPs.
#' Default: \code{filter.ma = NULL}.


#' @param calibrate.alleles (character) When ancestral allele
#' information is not available, REF / ALT alleles can be calibrated using
#' alleles count or depth (sequencing coverage for REF and ALT alleles).
#' Removing individuals will impact REF/ALT alleles. Ideally,
#' it should be based on observed read depth, when this information is available
#' for both alleles. By selecting \code{calibrate.alleles = "depth"}, radiator will
#' check if the information is available, if not, it will use
#' \code{calibrate.alleles = "count"}.
#' Using \code{calibrate.alleles = "ancestral"} will turn off calibration and
#' use the information preset available in the data.
#' Default: \code{calibrate.alleles = "depth"}.


#' @param filename (optional, character) Write to folder the MAF filtered
#' tidy dataset. The name will be appended \code{.rad}. With default, the filtered
#' data is only in the global environment.
#' Default: \code{filename = NULL}.

#' @inheritParams radiator_common_arguments

#' @section MAC, MAF or MAD ?:
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

#' @section Advance mode:
#'
#' \emph{dots-dots-dots ...} allows to pass several arguments for fine-tuning the function:
#' \enumerate{
#' \item \code{filter.common.markers} (optional, logical).
#' Default: \code{filter.common.markers = FALSE},
#' Documented in \code{\link{filter_common_markers}}.
#' \item \code{filter.monomorphic} (logical, optional) Should the monomorphic
#' markers present in the dataset be filtered out ?
#' Default: \code{filter.monomorphic = TRUE}.
#' Documented in \code{\link{filter_monomorphic}}.
#' \item \code{path.folder}: to write ouput in a specific path
#' (used internally in radiator).
#' Default: \code{path.folder = getwd()}.
#' If the supplied directory doesn't exist, it's created.
#' }


#' @section Interactive version:
#'
#' To help choose a threshold for the local and global MAF
#' use the interactive version.
#'
#' 2 steps in the interactive version:
#'
#' Step 1. Visualization and helper table.
#'
#' Step 2. Filtering markers based on MAC, MAF or MAD


#' @rdname filter_ma
#' @export
#' @return With \code{interactive.filter = FALSE}, a list in the global environment,
#' with 7 objects:
#' \enumerate{
#' \item $tidy.filtered.mac
#' \item $whitelist.markers
#' \item $blacklist.markers
#' \item $mac.data
#' \item $filters.parameters
#' }
#'
#' With \code{interactive.filter = TRUE}, a list with 4 additionnal objects are generated.
#' \enumerate{
#' \item $distribution.mac.global
#' \item $distribution.mac.local
#' \item $mac.global.summary
#' \item $mac.helper.table
#' }
#'
#' \strong{mac.helper.table:}
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
#' mac <- radiator::filter_ma(data = turtle.tidy.data)
#' # This will use the default: interactive version,
#' # a list is created and to view the filtered tidy data:
#' mac.tidy.data <- mac$tidy.filtered.mac
#'
#' # No user interaction
#'
#' # Using filter.ma = 4, is a good practice to remove mostly sequencing errors
#' # and assembly artifacts, because it requires the markers to be genotyped in
#' # 4 heterozygote individuals or 2 homozygote individuals for the alternate
#' # allele or 2 heterozygote individuals and 1 homozygote individual for the
#' # alternate allele. Overall, never less than 2 indiduals are required.
#'
#' mac <- radiator::filter_ma(
#'         data = turtle.gds, # using gds object
#'         ma.stats = "mac",
#'         filter.ma = 4,
#'         filename = "turtle.mac")
#'
#' # This will remove monomorphic markers (by default filter.monomorphic = TRUE)
#' # The filtered data will be written in the directory under the name:
#' # turtle.maf.rad
#' }



#' @references Good reading: Linck, E., Battey, C. (2019).
#' Minor allele frequency thresholds strongly affect population structure
#' inference with genomic data sets.
#' Molecular Ecology Resources 19(3), 639-647.
#' https://dx.doi.org/10.1111/1755-0998.12995


#' @note Thanks to Charles Perrier and Jeremy Gaudin for very useful comments
#' on previous version of this function.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


filter_ma <- function(
  data,
  interactive.filter = TRUE,
  ma.stats = "mac",
  filter.ma = NULL,
  calibrate.alleles = "depth",
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
  ...
) {

  ## testing
  # data <- "radiator_20221205@1134.gds"
  # interactive.filter = TRUE
  # ma.stats = "mac"
  # filter.ma = NULL
  # calibrate.alleles = "depth"
  # filename = NULL
  # parallel.core = parallel::detectCores() - 1
  # verbose = TRUE
  # path.folder = NULL
  # parameters <- NULL
  # internal <- FALSE


  if (is.null(filter.ma) && !interactive.filter) return(data)
  if (interactive.filter) verbose <- TRUE

  # Cleanup-------------------------------------------------------------------
  radiator_function_header(f.name = "filter_ma", verbose = verbose)
  if (!verbose) message("filter_ma...")
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (verbose) message("Execution date@time: ", file.date)
  old.dir <- getwd()
  opt.change <- getOption("width")
  options(width = 70)
  timing <- radiator_tic()
  on.exit(setwd(old.dir), add = TRUE)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(radiator_toc(timing), add = TRUE)
  on.exit(radiator_function_header(f.name = "filter_ma", start = FALSE, verbose = verbose), add = TRUE)

  # Function call and dotslist -------------------------------------------------
  rad.dots <- radiator_dots(
    func.name = as.list(sys.call())[[1]],
    fd = rlang::fn_fmls_names(),
    args.list = as.list(environment()),
    dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE),
    keepers = c("path.folder", "parameters", "internal"),
    verbose = FALSE
  )

  calibrate.alleles <- match.arg(
    arg = calibrate.alleles,
    choices = c("count", "depth", "ancestral"),
    several.ok = FALSE
  )

  ma.stats <- match.arg(
    arg = ma.stats,
    choices = c("mac", "maf", "mad"),
    several.ok = FALSE
  )

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("data is missing")

  # Folders---------------------------------------------------------------------
  path.folder <- generate_folder(
    f = path.folder,
    rad.folder = "filter_ma",
    internal = internal,
    file.date = file.date,
    verbose = verbose)

  # write the dots file
  write_rad(
    data = rad.dots,
    path = path.folder,
    filename = stringi::stri_join("radiator_filter_ma_args_", file.date, ".tsv"),
    tsv = TRUE,
    internal = internal,
    write.message = "Function call and arguments stored in: ",
    verbose = verbose
  )

  # Message about steps taken during the process ---------------------------------
  if (interactive.filter) {
    message("Interactive mode: on")
    message("Step 1. Visualization and helper table")
    message("Step 2. Filtering markers based on MA stats\n\n")
  }


  # Import data ---------------------------------------------------------------
  if (verbose) cli::cli_progress_step("Importing data ")

  # Detect format
  data.type <- radiator::detect_genomic_format(data)
  if (!data.type %in% c("tbl_df", "fst.file", "SeqVarGDSClass", "gds.file")) {
    rlang::abort("Input not supported for this function: read function documentation")
  }

  if (data.type %in% c("tbl_df", "fst.file")) {
    if (is.vector(data)) data <- radiator::tidy_wide(data, import.metadata = TRUE)
    data.type <- "tbl_df"
    n.diplo.samples <- dplyr::n_distinct(data$INDIVIDUALS) * 2
  }

  if (data.type %in% c("SeqVarGDSClass", "gds.file")) {
    radiator_packages_dep(package = "SeqArray", cran = FALSE, bioc = TRUE)

    if (data.type == "gds.file") {
      data <- radiator::read_rad(data, verbose = FALSE)
      data.type <- "SeqVarGDSClass"
    }
    n.diplo.samples <- nrow(
      extract_individuals_metadata(gds = data, ind.field.select = "INDIVIDUALS", whitelist = TRUE)
    ) * 2
  }
  cli::cli_progress_done()

  # Filter parameter file: initiate ------------------------------------------
  filters.parameters <- radiator_parameters(
    generate = TRUE,
    initiate = TRUE,
    update = FALSE,
    parameter.obj = parameters,
    data = data,
    path.folder = path.folder,
    file.date = file.date,
    internal = internal,
    verbose = verbose)

  # Whitelist and blacklist --------------------------------------------------
  want <- c("VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS")
  if (data.type == "tbl_df") {
    wl <- bl <- radiator::separate_markers(
      data = data,
      sep = "__",
      biallelic = radiator::detect_biallelic_markers(data = data),
      markers.meta.all.only = TRUE,
      parallel.core = parallel.core)
  }

  # Step 1. Visuals ----------------------------------------------------------
  if (interactive.filter) message("\nStep 1. Minor Allele Statistics visualization and helper table\n")

  # MAF calculation ----------------------------------------------------------
  mac.data <- minor_allele_stats(
    data = data,
    markers.meta = NULL, # will import inside function
    calibrate.alleles = calibrate.alleles,
    path.folder = path.folder,
    parallel.core = parallel.core,
    verbose = verbose
  )

  readr::write_tsv(x = mac.data, file = file.path(path.folder, "ma.global.tsv"))
  if (verbose) message("File written: ma.global.tsv")
  n.markers <- nrow(mac.data)

  # Stats
  want <- c("MAC_GLOBAL", "MAC_GLOBAL_CORR", "MAD_GLOBAL", "MAD_GLOBAL_CORR", "MAF_GLOBAL", "MAF_GLOBAL_COUNT_CORR", "MAF_GLOBAL_DEPTH_CORR")
  have <- names(mac.data)
  have <- purrr::keep(.x = have, .p = have %in% want)
  mac.stats <- purrr::map2_dfr(.x = mac.data[have], .y = have, .f = tibble_stats) %>%
    dplyr::mutate(GROUP = factor(x = GROUP, levels = have)) %>%
    dplyr::arrange(GROUP)

  readr::write_tsv(x = mac.stats, file = file.path(path.folder, "ma.summary.stats.tsv"))
  if (verbose) message("File written: ma.summary.stats.tsv")



  # Figures --------------------------------------------------------------------
  if (verbose) cli::cli_progress_step("Generating MAC helper tables and plots")
  q25 <- mac.stats$Q25[mac.stats$GROUP == "MAC_GLOBAL_CORR"]

  # boxplot
  bp.title <- stringi::stri_join("Minor Allele Stats (Global)\nQ25 = ", q25)
  mac.fig <- boxplot_stats(
    data = mac.stats,
    title = bp.title,
    x.axis.title = NULL,
    y.axis.title = "Minor Allele Stats (Global)",
    flip = TRUE,
    scale.log = TRUE,
    scientific.format = FALSE,
    bp.filename = stringi::stri_join("ma.boxplot_", file.date, ".pdf"),
    path.folder = path.folder,
    facet.columns = TRUE
    )
  suppressWarnings(print(mac.fig))

  # Histo
  histo.mac.global <- plot_histogram(
    data = mac.data, variable = "MAC_GLOBAL_CORR",
    x.axis.title =  "Minor Allele Count corrected\n(MAC global)",
    y.axis.title = "Number of SNPs"
  )
  print(histo.mac.global)

  # save
  ggplot2::ggsave(
    filename = file.path(path.folder, "distribution.mac.global.pdf"),
    plot = histo.mac.global,
    width = 10,
    height = 10,
    dpi = 300,
    units = "cm",
    limitsize = FALSE,
    useDingbats = FALSE
  )


  if (rlang::has_name(mac.data, "MAF_GLOBAL_DEPTH_CORR")) {
    histo.maf.global <- plot_histogram(
      data = mac.data, variable = "MAF_GLOBAL_DEPTH_CORR",
      x.axis.title = "Minor Allele Frequency depth corrected\n(MAF global)",
      y.axis.title = "Number of SNPs"
      )
    print(histo.maf.global)
    ggplot2::ggsave(
      filename = file.path(path.folder, "distribution.maf.global.pdf"),
      plot = histo.maf.global,
      width = 10,
      height = 10,
      dpi = 300,
      units = "cm",
      limitsize = FALSE,
      useDingbats = FALSE
    )
  }

  mac.min <- mac.stats$MIN[mac.stats$GROUP == "MAC_GLOBAL_CORR"]
  mac.max <- mac.stats$MAX[mac.stats$GROUP == "MAC_GLOBAL_CORR"]
  maf.min <- format(round(mac.stats$MIN[mac.stats$GROUP == "MAF_GLOBAL_COUNT_CORR"], 5), scientific = FALSE)
  maf.max <- format(round(mac.stats$MAX[mac.stats$GROUP == "MAF_GLOBAL_COUNT_CORR"], 5), scientific = FALSE)

  # Helper table for global and local MAF -------------------------------------

  how_many_markers <- function(threshold, x, stats = c("MAC_GLOBAL_CORR", "MAF_GLOBAL_COUNT_CORR", "MAF_GLOBAL_DEPTH_CORR", "MAD_GLOBAL_CORR")) {
    stats <- match.arg(
      arg = stats,
      choices = c("MAC_GLOBAL_CORR", "MAF_GLOBAL_COUNT_CORR", "MAF_GLOBAL_DEPTH_CORR", "MAD_GLOBAL_CORR"),
      several.ok = FALSE
    )
  if (rlang::has_name(x, stats)) nrow(dplyr::filter(.data = x, x[[stats]] >= threshold))
  }#End how_many_markers

  end.seq <- ceiling(0.2 * n.diplo.samples)


  # MAC
  mac.helper.table <- tibble::tibble(MAC = 1:end.seq) %>%
    dplyr::mutate(
      MAF = format(round(MAC / n.diplo.samples, 4), scientific = FALSE),
      WHITELISTED_MARKERS = purrr::map_int(.x = 1:end.seq, .f = how_many_markers, x = mac.data, stats = "MAC_GLOBAL_CORR"),
      BLACKLISTED_MARKERS = n.markers - WHITELISTED_MARKERS
    ) %>%
    readr::write_tsv(
      x = .,
      file = file.path(path.folder, "mac.helper.table.tsv")
      )

  # if (verbose) message("File written: mac.helper.table.tsv")

  mac.markers.plot <- radiator_helper_plot(
    data = rad_long(
      x = mac.helper.table[1:25, -2],
      cols = "MAC",
      names_to = "LIST",
      values_to = "MARKERS"
    ),
    stats = "MAC",
    x.axis.title =  "Minor Allele Count threshold (MAC)",
    x.breaks = seq(1, 25, by = 1),
    plot.filename = file.path(path.folder, "mac.markers.plot")
    )
  print(mac.markers.plot)

  # MAF
  if (rlang::has_name(mac.data, "MAF_GLOBAL_DEPTH_CORR")) {

    min.maf <- min(mac.data$MAF_GLOBAL_DEPTH_CORR, na.rm = TRUE)
    maf.round <- signif(min.maf, digits = 1)
    maf.digits <- nchar(maf.round) - 2

    maf.seq <- unique(signif(seq(from = maf.round, to = 0.2, by = signif(((0.2 - maf.round) / 100), digits = 1)), digits = 2))

    maf.helper.table <- tibble::tibble(MAF = maf.seq) %>%
      dplyr::mutate(
        WHITELISTED_MARKERS = purrr::map_int(.x = maf.seq, .f = how_many_markers, x = mac.data, stats = "MAF_GLOBAL_DEPTH_CORR"),
        BLACKLISTED_MARKERS = n.markers - WHITELISTED_MARKERS
      ) %>%
      readr::write_tsv(
        x = .,
        file = file.path(path.folder, "maf.helper.table.tsv"))

    # if (verbose) message("File written: maf.helper.table.tsv")


    maf.markers.plot <- radiator_helper_plot(
      data = rad_long(
        x = maf.helper.table[1:25,],
        cols = "MAF",
        names_to = "LIST",
        values_to = "MARKERS"
      ),
      stats = "MAF",
      x.axis.title =  "Minor Allele Frequency threshold (MAF)",
      plot.filename = file.path(path.folder, "maf.markers.plot")
    )
    print(maf.markers.plot)
    min.maf <- maf.round <- maf.digits <- maf.seq <- NULL
  }

  # MAD
  if (rlang::has_name(mac.data, "MAD_GLOBAL_CORR")) {
    mad.helper.table <- tibble::tibble(MAD = 1:end.seq) %>%
      dplyr::mutate(
        WHITELISTED_MARKERS = purrr::map_int(.x = 1:end.seq, .f = how_many_markers, x = mac.data, stats = "MAD_GLOBAL_CORR"),
        BLACKLISTED_MARKERS = n.markers - WHITELISTED_MARKERS
      ) %>%
      readr::write_tsv(
        x = .,
        file = file.path(path.folder, "mad.helper.table.tsv"))

    # if (verbose) message("File written: mad.helper.table.tsv")


    mad.markers.plot <- radiator_helper_plot(
      data = rad_long(
        x = mad.helper.table[1:25,],
        cols = "MAD",
        names_to = "LIST",
        values_to = "MARKERS"
      ),
      stats = "MAD",
      x.axis.title =  "Minor Allele Depth threshold (MAD)",
      x.breaks = seq(1, 25, by = 1),
      plot.filename = file.path(path.folder, "mad.markers.plot")
    )
    print(mad.markers.plot)
  }
  cli::cli_progress_done()


  if (verbose) message("MAC range (calibrated with count): [", mac.min, " - ", mac.max, "]")
  if (verbose) message("MAF range (calibrated with count): [", maf.min, " - ", maf.max, "]")

  if (rlang::has_name(mac.data, "MAF_GLOBAL_DEPTH_CORR")) {
    maf.min <- format(round(mac.stats$MIN[mac.stats$GROUP == "MAF_GLOBAL_DEPTH_CORR"], 5), scientific = FALSE)
    maf.max <- format(round(mac.stats$MAX[mac.stats$GROUP == "MAF_GLOBAL_DEPTH_CORR"], 5), scientific = FALSE)
    if (verbose) message("MAF range (calibrated with depth): [", maf.min, " - ", maf.max, "]")
  }

  # Step 2. Thresholds selection ---------------------------------------------
  if (interactive.filter) {
    want <- c("MAC_GLOBAL_CORR", "MAD_GLOBAL_CORR", "MAF_GLOBAL_DEPTH_CORR", "MAF_GLOBAL_COUNT_CORR")
    have <- names(mac.data)
    available.options <- purrr::keep(.x = have, .p = have %in% want)
    available.options <- stringi::stri_replace_all_fixed(
      str = available.options,
      pattern = want,
      replacement = c("mac", "mad", "maf", "maf"),
      vectorize_all = FALSE
    )

    message("\nStep 2. Minor Allele statistic\n")
    message("Your available options: ", stringi::stri_join(available.options, collapse = ", "))
    ma.stats <- radiator_question(
      x = "Choose the statistic: ", answer.opt = available.options)
  }


  if (interactive.filter) {
    if (ma.stats == "mad") {
      min.t <- 1L
      max.t <- max(mac.data$MAD_GLOBAL_CORR, na.rm = TRUE)
    }
    if (ma.stats == "mac") {
      min.t <- 1L
      max.t <- max(mac.data$MAC_GLOBAL_CORR, na.rm = TRUE)
    }
    if (ma.stats == "maf") {
      min.t <- 0
      max.t <- 1
    }
    message("\nStep 3. Filtering markers based on ", ma.stats, "\n")
    rqm <- paste0("Choose the ", ma.stats, " threshold: ")
    filter.ma <- radiator_question(
      x = rqm, minmax = c(min.t, max.t))
  }

  if (ma.stats == "mad") mac.data %<>% dplyr::filter(MAD_GLOBAL_CORR < filter.ma)
  if (ma.stats == "mac") mac.data %<>% dplyr::filter(MAC_GLOBAL_CORR < filter.ma)
  if (ma.stats == "maf") {
    if ("MAF_GLOBAL_DEPTH_CORR" %in% have) {
      mac.data %<>% dplyr::filter(MAF_GLOBAL_DEPTH_CORR < filter.ma)
    }
    if ("MAF_GLOBAL_COUNT_CORR" %in% have) {
      mac.data %<>% dplyr::filter(MAF_GLOBAL_COUNT_CORR < filter.ma)
    }
  }

  # Filtering ----------------------------------------------------------------
  if (data.type == "tbl_df") {# Whitelist and Blacklist of markers
    bl %<>%
      dplyr::filter(MARKERS %in% mac.data$MARKERS) %>%
      dplyr::mutate(FILTERS = "filter.ma")
    wl %<>% dplyr::filter(!MARKERS %in% mac.data$MARKERS)
    data %<>% dplyr::filter(MARKERS %in% wl$MARKERS)
  }

  if (data.type == "SeqVarGDSClass") {
    markers.meta <- extract_markers_metadata(gds = data, whitelist = FALSE) %>%
      dplyr::mutate(
        FILTERS = dplyr::if_else(
          VARIANT_ID %in% mac.data$VARIANT_ID, "filter.ma", FILTERS
        )
      )

    update_radiator_gds(
      gds = data,
      node.name = "markers.meta",
      value = markers.meta,
      sync = TRUE
    )

    wl <- dplyr::filter(markers.meta, FILTERS == "whitelist")
    bl <- dplyr::filter(markers.meta, FILTERS == "filter.ma")
  }

  #blacklist

  if (nrow(bl) > 0L) {
    write_rad(
      data = bl,
      path = path.folder,
      filename = stringi::stri_join("blacklist.markers.ma_", file.date, ".tsv"),
      tsv = TRUE, internal = internal, verbose = verbose
    )
  }

  # whitelist
  if (nrow(wl) > 0L) {
    write_rad(
      data = wl,
      path = path.folder,
      filename = stringi::stri_join("whitelist.markers.ma_", file.date, ".tsv"),
      tsv = TRUE, internal = internal, verbose = verbose)
  }

  if (!is.null(filename)) {
    if (data.type == "tbl_df") {
      tidy.name <- stringi::stri_join(filename, ".rad")
      if (verbose) message("Writing the MA filtered tidy data set: ", tidy.name)
      write_rad(data = data, path = file.path(path.folder, tidy.name))
    }
  }

  # Update parameters --------------------------------------------------------
  filters.parameters <- radiator_parameters(
    generate = FALSE,
    initiate = FALSE,
    update = TRUE,
    parameter.obj = filters.parameters,
    data = data,
    filter.name = "Filter MA",
    param.name = paste0(c("filter.ma", "ma.stats"), collapse = " & "),
    values = paste0(c(filter.ma, ma.stats), collapse = ", "),
    path.folder = path.folder,
    file.date = file.date,
    internal = internal,
    verbose = verbose)


  # results --------------------------------------------------------------------
  radiator_results_message(
    rad.message = stringi::stri_join("\nFilter ma stats: ", ma.stats, " and threshold: ", filter.ma),
    filters.parameters,
    internal,
    verbose
  )
  return(data)
}#End filter_ma


#' @title compute_maf
#' @description Compute MAF
#' @rdname compute_maf
#' @keywords internal
#' @export
compute_maf <- carrier::crate(function(x, biallelic) {
  `%>%` <- magrittr::`%>%`
  `%<>%` <- magrittr::`%<>%`

  if (rlang::has_name(x, "GT_BIN") && biallelic) {
    x %<>%
      dplyr::group_by(MARKERS, POP_ID) %>%
      dplyr::summarise(
        NN = as.numeric(2 * n()),
        PP = as.numeric(2 * length(GT_BIN[GT_BIN == 0])),
        PQ = as.numeric(length(GT_BIN[GT_BIN == 1])),
        QQ = as.numeric(2 * length(GT_BIN[GT_BIN == 2]))
      ) %>%
      # need this step because seen cases where the minor allele is not minor
      dplyr::mutate(
        PP = PP + PQ,
        QQ = QQ + PQ,
        PQ = NULL) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(
        NN_G = sum(NN),
        PP_G = sum(PP),
        QQ_G = sum(QQ)) %>%
      dplyr::group_by(MARKERS, POP_ID) %>%
      dplyr::mutate(
        ALT = dplyr::if_else(PP_G < QQ_G, PP, QQ),
        MAF_LOCAL = (ALT / NN),
        PP = NULL,
        # QQ = NULL,
        NN = NULL) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(
        ALT = dplyr::if_else(PP_G < QQ_G, PP_G, QQ_G),
        MAF_GLOBAL = (ALT / NN_G),
        ALT = NULL,
        PP_G = NULL,
        # QQ_G = NULL,
        NN_G = NULL) %>%
      dplyr::ungroup(.) %>%
      dplyr::rename(ALT_LOCAL = QQ, ALT_GLOBAL = QQ_G)
  } else {
    if (!rlang::has_name(x, "GT_VCF_NUC")) {
      x %<>%
        dplyr::select(MARKERS,POP_ID, INDIVIDUALS, GT) %>%
        dplyr::mutate(
          A1 = stringi::stri_sub(GT, 1, 3),
          A2 = stringi::stri_sub(GT, 4,6),
          GT = NULL
        ) %>%
        radiator::rad_long(
          x = .,
          cols = c("POP_ID", "INDIVIDUALS", "MARKERS"),
          names_to = "ALLELES",
          values_to = "GT"
        )

      maf.local <- x %>%
        dplyr::group_by(MARKERS, POP_ID, GT) %>%
        dplyr::tally(.) %>%
        dplyr::ungroup(.) %>%
        tidyr::complete(data = ., POP_ID, tidyr::nesting(MARKERS, GT), fill = list(n = 0)) %>%
        dplyr::group_by(MARKERS, POP_ID) %>%
        dplyr::mutate(n.al.tot = sum(n)) %>%
        dplyr::filter(n == min(n)) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
        dplyr::summarise(MAF_LOCAL = n / n.al.tot, ALT_LOCAL = n, .groups = "drop") %>%
        dplyr::select(MARKERS, POP_ID, MAF_LOCAL, ALT_LOCAL)

      x %<>%
        dplyr::group_by(MARKERS, GT) %>%
        dplyr::tally(.) %>%
        dplyr::group_by(MARKERS) %>%
        dplyr::mutate(n.al.tot = sum(n)) %>%
        dplyr::filter(n == min(n)) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
        dplyr::summarise(MAF_GLOBAL = n / n.al.tot, ALT_GLOBAL = n, .groups = "drop") %>%
        dplyr::select(MARKERS, MAF_GLOBAL, ALT_GLOBAL) %>%
        dplyr::left_join(maf.local, by = c("MARKERS")) %>%
        dplyr::select(MARKERS, POP_ID, MAF_LOCAL, ALT_LOCAL, MAF_GLOBAL, ALT_GLOBAL)
      maf.local <- NULL
    } else {
      x %<>%
        dplyr::select(MARKERS, POP_ID, INDIVIDUALS, GT_VCF_NUC) %>%
        tidyr::separate(
          data = .,
          col = GT_VCF_NUC, into = c("A1", "A2"),
          sep = "/",
          extra = "drop", remove = TRUE
        ) %>%
        radiator::rad_long(
          x = .,
          cols = tidyselect::any_of(c("MARKERS", "INDIVIDUALS", "POP_ID")),
          names_to = "ALLELE_GROUP",
          values_to = "HAPLOTYPES"
        ) %>%
        dplyr::select(-ALLELE_GROUP) %>%
        dplyr::group_by(MARKERS, HAPLOTYPES, POP_ID) %>%
        dplyr::tally(.) %>%
        dplyr::ungroup(.) %>%
        tidyr::complete(data = ., POP_ID, tidyr::nesting(MARKERS, HAPLOTYPES), fill = list(n = 0)) %>%
        dplyr::group_by(MARKERS, POP_ID) %>%
        dplyr::mutate(N_LOCAL = sum(n)) %>%
        dplyr::group_by(MARKERS) %>%
        dplyr::mutate(N_GLOBAL = sum(n)) %>%
        dplyr::arrange(MARKERS, POP_ID) %>%
        dplyr::group_by(MARKERS, POP_ID, HAPLOTYPES) %>%
        dplyr::mutate(MAF_LOCAL = n / N_LOCAL) %>%
        dplyr::group_by(MARKERS, HAPLOTYPES) %>%
        dplyr::mutate(
          ALT_GLOBAL = sum(n),
          MAF_GLOBAL = ALT_GLOBAL / N_GLOBAL,
          N_LOCAL = NULL,
          N_GLOBAL = NULL
        ) %>%
        dplyr::rename(ALT_LOCAL = n) %>%
        dplyr::ungroup(.)

      ref.info <- dplyr::distinct(x, MARKERS, HAPLOTYPES, MAF_GLOBAL) %>%
        dplyr::group_by(MARKERS) %>%
        dplyr::filter(MAF_GLOBAL == max(MAF_GLOBAL)) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
        dplyr::ungroup(.) %>%
        dplyr::mutate(REF = rep("REF", n()), MAF_GLOBAL = NULL) %>%
        dplyr::bind_rows(
          dplyr::distinct(x, MARKERS, HAPLOTYPES, MAF_GLOBAL) %>%
            dplyr::group_by(MARKERS) %>%
            dplyr::filter(MAF_GLOBAL == min(MAF_GLOBAL)) %>%
            dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
            dplyr::ungroup(.) %>%
            dplyr::mutate(REF = rep("MAF", n()), MAF_GLOBAL = NULL)
        )

      x %<>% dplyr::left_join(ref.info, by = c("MARKERS", "HAPLOTYPES")) %>%
        dplyr::mutate(REF = stringi::stri_replace_na(REF, replacement = "ALT"))
      ref.info <- NULL
    }
  }
  return(x)
})#End compute_maf

#' @title minor_allele_stats
#' @description Generate Minor Allele necessary statistics
#' @rdname minor_allele_stats
#' @export
#' @keywords internal
minor_allele_stats <- function(
  data,
  calibrate.alleles = c("count", "depth", "ancestral"),
  blacklist.markers = NULL,
  markers.meta = NULL,
  path.folder = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE
) {
  if (verbose) cli::cli_progress_step("Calculating Minor Allele statistics")
  n.inversion.count <- NULL
  n.inversion.depth <- NULL

  data.type <- detect_genomic_format(data) # gds or tbl_df

  # required stuff -------------------------------------------------------------
  if (data.type == "tbl_df") {
    if (rlang::has_name(data, "GT_BIN")) {
      data %<>% dplyr::filter(!is.na(GT_BIN))
    } else {
      data %<>% dplyr::filter(GT != "000000")
    }

    markers.meta <- dplyr::distinct(data, MARKERS)
    if (!is.null(blacklist.markers)) {
      markers.meta %<>% dplyr::filter(!MARKERS %in% blacklist.markers$MARKERS)
    }
    n.markers <- nrow(markers.meta)
    markers.meta <- NULL

    if (n.markers > 30000) {
      mac.data <- radiator_future(
        .x = data,
        .f = ma_one,
        flat.future = "dfr",
        split.vec = TRUE,
        split.with = "MARKERS",
        split.chunks = 10L,
        parallel.core = parallel.core,
        forking = TRUE
      )
    } else {
      mac.data <- ma_one(x = data)
    }
    mac.data %<>% dplyr::arrange(MARKERS)
  }# End tibble

  # GDS
  if (data.type == "SeqVarGDSClass") {
    if (is.null(markers.meta)) {
      mac.data <- extract_markers_metadata(
        gds = data,
        markers.meta.select = c("VARIANT_ID", "MARKERS"),
        whitelist = TRUE)
    } else {
      mac.data <- markers.meta
      markers.meta <- NULL
      if (rlang::has_name(mac.data, "FILTERS")) {
        mac.data %<>% dplyr::filter(FILTERS == "whitelist")
      }
      mac.data %<>% dplyr::select(VARIANT_ID, MARKERS)
    }

    # blacklist markers
    if (!is.null(blacklist.markers)) {
      mac.data %<>% dplyr::filter(!MARKERS %in% blacklist.markers$MARKERS)
      sync_gds(gds = data, variant.id = mac.data$VARIANT_ID)
    }
    n.markers <- nrow(mac.data)


    # calibration of alleles ---------------------------------------------------
    ad <- NULL
    if (calibrate.alleles != "ancestral") {
      if (calibrate.alleles == "depth") {
        # check if AD is available...
        ad <- check_coverage(
          gds = data,
          stacks.haplo.check = TRUE,
          genotypes.metadata.check = TRUE,
          dart.check = TRUE
        )
        if ("AD" %in% ad) {
          # get the read depth for each alleles
          calibrate.alleles <- "depth"
          ad <- extract_coverage(
            gds = data,
            individuals = FALSE,
            markers = TRUE,
            dp = FALSE,
            ad = TRUE,
            coverage.stats = "sum",
            subsample.info = 1,
            verbose = FALSE,
            parallel.core = parallel.core
          ) %$%
            m.info
        } else {
          calibrate.alleles <- "count"
          message("depth info not available, switching to count")
        }
      }
    }

    # get the count information
    # extract the allele count from GDS
    ac <- SeqArray::seqAlleleCount(
      gdsfile = data,
      ref.allele = NULL,
      parallel = parallel.core)

    # check the number of alt alleles
    n.al.max <- max(purrr::map_int(.x = ac, .f = length), na.rm = TRUE)

    # when more than 2 alternate alleles...
    if (n.al.max > 2) {
      message("\n\nCaution: More than 2 alleles detected in the dataset")
      wanted.info <- c("VARIANT_ID", "CHROM", "LOCUS", "POS", "MARKERS", "REF", "ALT", "N_ALLELES")

      prob.markers <- mac.data %>%
        dplyr::mutate(N_ALLELES = purrr::map_int(.x = ac, .f = length)) %>%
        dplyr::filter(N_ALLELES > 2) %>%
        dplyr::select(dplyr::any_of(wanted.info)) %>%
        write_rad(
          data = .,
          path = path.folder,
          filename = "markers_number_alleles_problem.tsv",
          tsv = TRUE,
          write.message = "standard",
          verbose = TRUE
        )
    }
    n.al.max <- NULL

    # make sure ac and ad are the same length
    if (calibrate.alleles == "depth") {
      if (length(ac) != nrow(ad)) rlang::abort("Problem with the data, contact author")
    }

    # compile the info----------------------------------------------------------
    # Add count info
    ac %<>%
      unlist(.) %>%
      matrix(
        data = .,
        nrow = n.markers, ncol = 2, byrow = TRUE,
        dimnames = list(rownames = mac.data$VARIANT_ID,
                        colnames = c("REF_COUNT", "ALT_COUNT"))) %>%
      tibble::as_tibble(.)
    mac.data %<>% dplyr::bind_cols(ac)
    ac <- NULL


    # Add depth info
    if (calibrate.alleles == "depth") {
      mac.data %<>% dplyr::bind_cols(ad)
      ad <- NULL
    }

    # calibrate ----------------------------------------------------------------
    # monomorphic marker check
    mono <- dplyr::filter(.data = mac.data, REF_COUNT == 0)
    # will have to deal with that...

    mac.data %<>% dplyr::filter(REF_COUNT != 0)

    if (calibrate.alleles == "ancestral") {
      mac.data %<>%
        dplyr::mutate(
          TOTAL_COUNT = ALT_COUNT + REF_COUNT,
          MAF_GLOBAL = ALT_COUNT / TOTAL_COUNT,
          TOTAL_COUNT = NULL,
          REF_COUNT = NULL
        ) %>%
        dplyr::rename(MAC_GLOBAL = ALT_COUNT)
    }

    if (calibrate.alleles != "ancestral") {
      # Report number of inversion based on count
      n.inversion.count <- dplyr::select(mac.data, ALT_COUNT, REF_COUNT) %>%
        dplyr::filter(ALT_COUNT >= REF_COUNT) %>%
        nrow()



      # NOTE: report tibble of mac based on original mac / maf ?

      mac.data %<>%
        dplyr::mutate(
          INVERSION = NULL,
          MAC_GLOBAL_CORR = dplyr::if_else(ALT_COUNT < REF_COUNT, ALT_COUNT, REF_COUNT),
          TOTAL_COUNT = ALT_COUNT + REF_COUNT, # doesnt change
          MAF_GLOBAL = ALT_COUNT / TOTAL_COUNT,
          MAF_GLOBAL_COUNT_CORR = MAC_GLOBAL_CORR / TOTAL_COUNT,
          TOTAL_COUNT = NULL,
          REF_COUNT = NULL
        ) %>%
        dplyr::rename(MAC_GLOBAL = ALT_COUNT)

      if (calibrate.alleles == "depth") {
        n.inversion.depth <- dplyr::select(mac.data, ALT_DEPTH_TOTAL, REF_DEPTH_TOTAL) %>%
          dplyr::filter(ALT_DEPTH_TOTAL >= REF_DEPTH_TOTAL) %>%
          nrow()

        mac.data %<>%
          dplyr::mutate(
            MAD_GLOBAL_CORR = dplyr::if_else(ALT_DEPTH_TOTAL < REF_DEPTH_TOTAL, ALT_DEPTH_TOTAL, REF_DEPTH_TOTAL),
            MAF_GLOBAL_DEPTH_CORR = MAD_GLOBAL_CORR / (ALT_DEPTH_TOTAL + REF_DEPTH_TOTAL),
            REF_DEPTH_TOTAL = NULL
          ) %>%
          dplyr::rename(MAD_GLOBAL = ALT_DEPTH_TOTAL) %>%
          dplyr::relocate(MAC_GLOBAL, MAC_GLOBAL_CORR, MAD_GLOBAL, MAD_GLOBAL_CORR, MAF_GLOBAL, MAF_GLOBAL_COUNT_CORR, MAF_GLOBAL_DEPTH_CORR, .after = tidyselect::last_col())
      }
    }
  } # END GDS
  cli::cli_progress_done()

  # Report number of inversion based on depth
  if (!is.null(n.inversion.count) && n.inversion.count > 0) {
    if (verbose) message("Calibration information")
    if (verbose) message("Number of REF/ALT inversions based on count: ", n.inversion.count)
  }
  if (!is.null(n.inversion.depth) && n.inversion.depth > 0) {
    if (verbose) message("Number of REF/ALT inversions based on read depth: ", n.inversion.depth)
  }

  return(mac.data)
}#End minor_allele_stats

#' @title ma_one
#' @description mac without parallel
#' @rdname ma_one
#' @keywords internal
#' @export
ma_one <- carrier::crate(function(x) {
  `%>%` <- magrittr::`%>%`
  `%<>%` <- magrittr::`%<>%`

  skip <- FALSE

  # GT_BIN
  if (rlang::has_name(x, "GT_BIN")) {
    mac.data <- x %>%
      dplyr::select(MARKERS, INDIVIDUALS, GT_BIN) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::summarise(
        PP = as.numeric(2 * length(GT_BIN[GT_BIN == 0])),
        PQ = as.numeric(length(GT_BIN[GT_BIN == 1])),
        QQ = as.numeric(2 * length(GT_BIN[GT_BIN == 2])),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        PP = PP + PQ,
        QQ = QQ + PQ,
        PQ = NULL,
        MAC_GLOBAL = dplyr::if_else(PP < QQ, PP, QQ),
        PP = NULL, QQ = NULL)
    skip <- TRUE
  }

  # GT
  if (skip && rlang::has_name(x, "GT")) {
    mac.data <- x %>%
      dplyr::select(MARKERS, INDIVIDUALS, GT) %>%
      radiator::separate_gt(x = ., gt = "GT", gather = TRUE, exclude = c("MARKERS", "INDIVIDUALS"), split.chunks = 1L) %>%
      dplyr::count(MARKERS, ALLELES) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::filter(n == min(n)) %>%
      dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
      dplyr::summarise(MAC_GLOBAL = n, .groups = "drop") %>%
      dplyr::select(MARKERS, MAC_GLOBAL)
  }
  mac.data$MAC_GLOBAL %<>% as.integer(.)
  return(mac.data)
})#End ma_one
