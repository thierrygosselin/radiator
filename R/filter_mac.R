# Minor Allele Count filtering
#' @name filter_mac
#' @title MAC filter
#' @description The Minor Allele Count function
#' remove/blacklist markers based on Minor/Alternate Allele Count (MAC). Use it
#' to remove noise, sequencing errors or low polymorphism markers. Some
#' analysis performs better with the full spectrum of allele frequency, so careful
#' with high threshold that inevitably results in biaises.
#'
#' \strong{Filter targets}: Marker's alternate allele(s)
#'
#' \strong{Statistics}: count of the allele per markers. This filter as NO concept
#' of strata/populations. It is computed globally.

#' @param filter.mac (optional, integer) The number of alternate (minor) allele
#' to keep a SNP (see details). e.g. \code{filter.mac = 3}. This argument as no
#' concept of locus, local or strata or pop. It's applied globally, by SNPs.
#' Default: \code{filter.mac = NULL}.

#' @param filename (optional, character) Write to folder the MAF filtered
#' tidy dataset. The name will be appended \code{.rad}. With default, the filtered
#' data is only in the global environment.
#' Default: \code{filename = NULL}.

#' @inheritParams radiator_common_arguments

#' @section MAC or MAF ?:
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
#' Step 2. Filtering markers based on MAC


#' @rdname filter_mac
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
#' mac <- radiator::filter_mac(data = turtle.tidy.data)
#' # This will use the default: interactive version,
#' # a list is created and to view the filtered tidy data:
#' mac.tidy.data <- mac$tidy.filtered.mac
#'
#' # No user interaction
#'
#' # Using filter.mac = 4, is a good practice to remove mostly sequencing errors
#' # and assembly artifacts, because it requires the markers to be genotyped in
#' # 4 heterozygote individuals or 2 homozygote individuals for the alternate
#' # allele or 2 heterozygote individuals and 1 homozygote individual for the
#' # alternate allele. Overall, never less than 2 indiduals are required.
#'
#' mac <- radiator::filter_mac(
#'         data = turtle.gds, # using gds object
#'         filter.mac = 4,
#'         filename = "turtle.mac")
#'
#' # This will remove monomorphic markers (by default filter.monomorphic = TRUE)
#' # The filtered data will be written in the directory under the name:
#' # turtle.maf.rad
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

filter_mac <- function(
  interactive.filter = TRUE,
  data,
  filter.mac = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
  ...
) {

  ## testing
  # interactive.filter = TRUE
  # filter.mac = 5
  # filename = NULL
  # parallel.core = parallel::detectCores() - 1
  # verbose = TRUE
  # path.folder = "00_filter_mac"
  # parameters <- NULL
  if (!is.null(filter.mac) || interactive.filter) {
    if (interactive.filter) verbose <- TRUE
    if (verbose) {
      cat("################################################################################\n")
      cat("############################## radiator::filter_mac ############################\n")
      cat("################################################################################\n")
    }
    if (!verbose) message("filter_mac...")

    # Cleanup-------------------------------------------------------------------
    file.date <- format(Sys.time(), "%Y%m%d@%H%M")
    if (verbose) message("Execution date@time: ", file.date)
    old.dir <- getwd()
    opt.change <- getOption("width")
    options(width = 70)
    timing <- proc.time()# for timing
    # res <- list()
    #back to the original directory and options
    on.exit(setwd(old.dir), add = TRUE)
    on.exit(options(width = opt.change), add = TRUE)
    on.exit(timing <- proc.time() - timing, add = TRUE)
    on.exit(if (verbose) message("\nComputation time, overall: ", round(timing[[3]]), " sec"), add = TRUE)
    on.exit(if (verbose) cat("############################ completed filter_mac ##############################\n"), add = TRUE)

    # Function call and dotslist -------------------------------------------------
    rad.dots <- radiator_dots(
      func.name = as.list(sys.call())[[1]],
      fd = rlang::fn_fmls_names(),
      args.list = as.list(environment()),
      dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE),
      keepers = c("path.folder", "parameters", "internal"),
      verbose = verbose
    )

    # Checking for missing and/or default arguments ------------------------------
    if (missing(data)) rlang::abort("data is missing")

    # Folders---------------------------------------------------------------------
    path.folder <- generate_folder(
      f = path.folder,
      rad.folder = "filter_mac",
      internal = internal,
      file.date = file.date,
      verbose = verbose)

    # write the dots file
    write_rad(
      data = rad.dots,
      path = path.folder,
      filename = stringi::stri_join("radiator_filter_mac_args_", file.date, ".tsv"),
      tsv = TRUE,
      internal = internal,
      verbose = verbose
    )

    # Message about steps taken during the process ---------------------------------
    if (interactive.filter) {
      message("Interactive mode: on\n")
      message("Step 1. Visualization and helper table")
      message("Step 2. Filtering markers based on MAC\n\n")
    }

    # Detect format --------------------------------------------------------------
    data.type <- radiator::detect_genomic_format(data)
    if (!data.type %in% c("tbl_df", "fst.file", "SeqVarGDSClass", "gds.file")) {
      rlang::abort("Input not supported for this function: read function documentation")
    }

    # Import data ---------------------------------------------------------------
    if (verbose) message("Importing data ...")
    if (data.type %in% c("tbl_df", "fst.file")) {
      if (is.vector(data)) data <- radiator::tidy_wide(data, import.metadata = TRUE)
      data.type <- "tbl_df"
      n.diplo.samples <- dplyr::n_distinct(data$INDIVIDUALS) * 2
    }

    if (data.type %in% c("SeqVarGDSClass", "gds.file")) {
      if (!"SeqVarTools" %in% utils::installed.packages()[,"Package"]) {
        rlang::abort('Please install SeqVarTools for this option:\n
             install.packages("BiocManager")
             BiocManager::install("SeqVarTools")')
      }

      if (data.type == "gds.file") {
        data <- radiator::read_rad(data, verbose = verbose)
        data.type <- "SeqVarGDSClass"
      }
      n.diplo.samples <- nrow(
        extract_individuals_metadata(gds = data, ind.field.select = "INDIVIDUALS", whitelist = TRUE)
      ) * 2
    }

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
      wl <- bl <- radiator::separate_markers(data = data, sep = "__",
                                             markers.meta.all.only = TRUE,
                                             parallel.core = parallel.core)
      markers.meta <- NULL
    } else {
      markers.meta <- extract_markers_metadata(gds = data)
    }

    # MAF calculation ----------------------------------------------------------
    if (verbose) message("Calculating GLOBAL MAC")
    mac.data <- compute_mac(data = data,
                            markers.meta = markers.meta,
                            parallel.core = parallel.core)

    readr::write_tsv(x = mac.data, path = file.path(path.folder, "mac.global.tsv"))
    if (verbose) message("File written: maf.global.tsv")
    n.markers <- nrow(mac.data)

    # Step 1. Visuals ----------------------------------------------------------
    if (interactive.filter) message("\nStep 1. MAC visualization and helper table\n")

    # Stats
    mac.stats <- tibble_stats(x = mac.data$MAC_GLOBAL, group = "unfiltered")

    q25 <- mac.stats$Q25
    mac.stats.low <- tibble_stats(x = mac.data$MAC_GLOBAL[mac.data$MAC_GLOBAL <= q25], group = "below Q25")

    mac.stats.data <- dplyr::bind_rows(mac.stats, mac.stats.low) %>%
      dplyr::mutate(GROUP = factor(x = GROUP,
                                   levels = c("unfiltered", "below Q25"))) %>%
      dplyr::arrange(GROUP)

    readr::write_tsv(x = mac.stats.data, path = file.path(path.folder, "mac.summary.stats.tsv"))
    if (verbose) message("File written: mac.summary.stats.tsv")


    # boxplot
    bp.title <- stringi::stri_join("Minor Allele Count (Global)\nQ25 = ", q25)
    mac.fig <- boxplot_stats(
      data = mac.stats.data,
      title = bp.title,
      x.axis.title = NULL,
      y.axis.title = "Minor Allele Count (Global)",
      bp.filename = stringi::stri_join("mac.boxplot_", file.date, ".pdf"),
      path.folder = path.folder, facet.columns = TRUE)
    print(mac.fig)

    # Histo
    # hist(mac.data$MAC_GLOBAL)
    histo.mac.global <- ggplot2::ggplot(mac.data, ggplot2::aes(x = MAC_GLOBAL)) +
      ggplot2::geom_histogram(bins = 30) +
      ggplot2::labs(y = "Number of SNPs", x = "Minor Allele Count (MAC)") +
      ggplot2::theme(
        legend.position = "none",
        axis.title.x = ggplot2::element_text(size = 12, face = "bold"),
        axis.title.y = ggplot2::element_text(size = 12, face = "bold")
      ) +
      ggplot2::theme_minimal()
    print(histo.mac.global)

    # save
    ggplot2::ggsave(
      filename = file.path(path.folder, "distribution.mac.global.pdf"),
      plot = histo.mac.global, width = 10, height = 10, dpi = 300, units = "cm", useDingbats = FALSE)

    if (verbose) message("MAC range: [", mac.stats$MIN, " - ", mac.stats$MAX, "]")
    if (verbose) message("MAF range: [", format(round(mac.stats$MIN/n.diplo.samples, 4), scientific = FALSE), " - ", format(round(mac.stats$MAX/n.diplo.samples, 4), scientific = FALSE), "]")

    # Helper table for global and local MAF -------------------------------------
    if (verbose) message("Generating MAC helper table...")
    how_many_markers <- function(threshold, x) {
      nrow(dplyr::filter(x, MAC_GLOBAL >= threshold))
    }#End how_many_markers

    end.seq <- ceiling(0.2 * n.diplo.samples)

    maf.helper.table <- tibble::tibble(MAC = 1:end.seq) %>%
      dplyr::mutate(
        MAF = format(round(MAC / n.diplo.samples, 4), scientific = FALSE),
        WHITELISTED_MARKERS = purrr::map_int(.x = 1:end.seq, .f = how_many_markers, x = mac.data),
        BLACKLISTED_MARKERS = n.markers - WHITELISTED_MARKERS
      ) %>%
      readr::write_tsv(
        x = .,
        path = file.path(path.folder, "mac.helper.table.tsv"))

    if (verbose) message("File written: maf.helper.table.tsv")

    mac.markers.plot <- ggplot2::ggplot(
      data = tidyr::gather(
        data = maf.helper.table[1:25, -2],
        key = LIST, value = MARKERS, -MAC),
      ggplot2::aes(x = MAC, y = MARKERS)) +
      ggplot2::geom_line() +
      ggplot2::geom_point(size = 2, shape = 21, fill = "white") +
      ggplot2::scale_x_continuous(name = "Minor Allele Count threshold (MAC)", breaks = seq(1, 25, by = 1)) +
      ggplot2::scale_y_continuous(name = "Number of markers")+#, breaks = y.breaks, limits = c(0, y.breaks.max)) +
      ggplot2::theme(
        axis.title.x = ggplot2::element_text(size = 10, face = "bold"),
        axis.title.y = ggplot2::element_text(size = 10, face = "bold"),
        axis.text.x = ggplot2::element_text(size = 8), #angle = 90, hjust = 1, vjust = 0.5),
        strip.text.x = ggplot2::element_text(size = 10, face = "bold")
      ) +
      ggplot2::theme_bw()+
      ggplot2::facet_grid(LIST ~. , scales = "free", space = "free")
    print(mac.markers.plot)
    # save
    ggplot2::ggsave(
      filename = file.path(path.folder, "mac.markers.plot.pdf"),
      plot = mac.markers.plot, width = 20, height = 15, dpi = 300, units = "cm", useDingbats = FALSE)

    # Step 2. Thresholds selection ---------------------------------------------
    if (interactive.filter) {
      message("\nStep 2. Filtering markers based on MAC\n")
      filter.mac <- radiator_question(
        x = "Choose the filter.mac threshold: ", minmax = c(1, n.markers))
    }

    mac.data %<>% dplyr::filter(MAC_GLOBAL < filter.mac)

    # Filtering ----------------------------------------------------------------
    if (data.type == "tbl_df") {# Whitelist and Blacklist of markers
      bl %<>%
        dplyr::filter(MARKERS %in% mac.data$MARKERS) %>%
        dplyr::mutate(FILTERS = "filter.mac")
      wl %<>% dplyr::setdiff(bl)
      data %<>% dplyr::filter(MARKERS %in% wl$MARKERS)
    } else {
      # Whitelist and Blacklist of markers
      # saving whitelist and blacklist
      markers.meta %<>%
        dplyr::mutate(
          FILTERS = dplyr::if_else(
            VARIANT_ID %in% mac.data$VARIANT_ID, "filter.mac", FILTERS
          )
        )
      update_radiator_gds(
        gds = data,
        node.name = "markers.meta",
        value = markers.meta,
        sync = TRUE
      )
      wl <- dplyr::filter(markers.meta, FILTERS == "whitelist")
      bl <- dplyr::filter(markers.meta, FILTERS == "filter.mac")
    }
    readr::write_tsv(
      x = wl,
      path = file.path(path.folder, "whitelist.markers.mac.tsv"),
      append = FALSE, col_names = TRUE)
    readr::write_tsv(
      x = bl,
      path = file.path(path.folder, "blacklist.markers.mac.tsv"),
      append = FALSE, col_names = TRUE)
    if (verbose) message("File written: whitelist.markers.mac.tsv")
    if (verbose) message("File written: blacklist.markers.mac.tsv")


    if (!is.null(filename)) {
      if (data.type == "tbl_df") {
        tidy.name <- stringi::stri_join(filename, ".rad")
        if (verbose) message("Writing the MAF filtered tidy data set: ", tidy.name)
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
      filter.name = "Filter MAC",
      param.name = "filter.mac",
      values = filter.mac,
      path.folder = path.folder,
      file.date = file.date,
      internal = internal,
      verbose = verbose)

    # if (filters.parameters$filters.parameters$BLACKLIST == 0) {
    #   file.remove(path.folder)
    #   if (verbose) message("Folder removed: ", folder_short(path.folder))
    # }

    # results --------------------------------------------------------------------
    radiator_results_message(
      rad.message = stringi::stri_join("\nFilter mac threshold: ", filter.mac),
      filters.parameters,
      internal,
      verbose
    )
  }
  return(data)
}#End filter_mac


#' @title compute_maf
#' @description Compute MAF
#' @rdname compute_maf
#' @keywords internal
#' @export
compute_maf <- function(x, biallelic) {
  if (tibble::has_name(x, "GT_BIN") && biallelic) {
    x <- x %>%
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
    if (!tibble::has_name(x, "GT_VCF_NUC")) {
      x <- x %>%
        dplyr::select(MARKERS,POP_ID, INDIVIDUALS, GT) %>%
        dplyr::mutate(
          A1 = stringi::stri_sub(GT, 1, 3),
          A2 = stringi::stri_sub(GT, 4,6)
        ) %>%
        dplyr::select(MARKERS, POP_ID, INDIVIDUALS, A1, A2) %>%
        tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS, POP_ID))

      maf.local <- x %>%
        dplyr::group_by(MARKERS, POP_ID, GT) %>%
        dplyr::tally(.) %>%
        dplyr::ungroup(.) %>%
        tidyr::complete(data = ., POP_ID, tidyr::nesting(MARKERS, GT), fill = list(n = 0)) %>%
        dplyr::group_by(MARKERS, POP_ID) %>%
        dplyr::mutate(n.al.tot = sum(n)) %>%
        dplyr::filter(n == min(n)) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
        dplyr::summarise(MAF_LOCAL = n / n.al.tot, ALT_LOCAL = n) %>%
        dplyr::ungroup(.) %>%
        dplyr::select(MARKERS, POP_ID, MAF_LOCAL, ALT_LOCAL)

      x <- x %>%
        dplyr::group_by(MARKERS, GT) %>%
        dplyr::tally(.) %>%
        dplyr::group_by(MARKERS) %>%
        dplyr::mutate(n.al.tot = sum(n)) %>%
        dplyr::filter(n == min(n)) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
        dplyr::summarise(MAF_GLOBAL = n / n.al.tot, ALT_GLOBAL = n) %>%
        dplyr::ungroup(.) %>%
        dplyr::select(MARKERS, MAF_GLOBAL, ALT_GLOBAL) %>%
        dplyr::left_join(maf.local, by = c("MARKERS")) %>%
        dplyr::select(MARKERS, POP_ID, MAF_LOCAL, ALT_LOCAL, MAF_GLOBAL, ALT_GLOBAL)
      maf.local <- NULL
    } else {
      x <- x %>%
        dplyr::select(MARKERS, POP_ID, INDIVIDUALS, GT_VCF_NUC) %>%
        tidyr::separate(
          data = .,
          col = GT_VCF_NUC, into = c("A1", "A2"),
          sep = "/",
          extra = "drop", remove = TRUE
        ) %>%
        tidyr::gather(
          data = ., key = ALLELE_GROUP, value = HAPLOTYPES,
          -dplyr::one_of(c("MARKERS", "INDIVIDUALS", "POP_ID"))) %>%
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

      x <- dplyr::left_join(x, ref.info, by = c("MARKERS", "HAPLOTYPES")) %>%
        dplyr::mutate(REF = stringi::stri_replace_na(REF, replacement = "ALT"))
      ref.info <- NULL
    }
  }
  return(x)
}#End compute_maf

#' @title compute_mac
#' @description Compute global MAC (Minor Allele Count)
#' @rdname compute_mac
#' @export
#' @keywords internal
compute_mac <- function (
  data,
  blacklist.markers = NULL,
  markers.meta = NULL,
  parallel.core = parallel::detectCores() - 1
) {
  data.type <- detect_genomic_format(data) # gds or tbl_df
  if (data.type == "tbl_df") {
    if (tibble::has_name(data, "GT_BIN")) {
      data %<>% dplyr::filter(!is.na(GT_BIN))
    } else {
      data %<>% dplyr::filter(GT != "000000")
    }
    markers.meta <- dplyr::distinct(data, MARKERS)
    if (!is.null(blacklist.markers)) {
      markers.meta %<>% dplyr::filter(!MARKERS %in% blacklist.markers$MARKERS)
    }
    n.markers <- nrow(markers.meta)
    if (n.markers > 10000) {
      mac.data <- data %>%
        dplyr::left_join(
          markers.meta %>%
            dplyr::mutate(SPLIT_VEC = split_vec_row(
              markers.meta,
              cpu.rounds = ceiling(n.markers/10000),
              parallel.core = parallel.core))
          , by = "MARKERS") %>%
        split(x = ., f = .$SPLIT_VEC) %>%
        .radiator_parallel_mc(
          X = .,
          FUN = mac_one,
          mc.cores = parallel.core
        ) %>%
        dplyr::bind_rows(.)
    } else {
      mac.data <- mac_one(x = data)
    }
    mac.data %<>% dplyr::arrange(MARKERS)
    markers.meta <- NULL
  } else {# GDS
    if (is.null(markers.meta)) {
      markers.meta <- extract_markers_metadata(
        gds = data,
        markers.meta.select = c("VARIANT_ID", "MARKERS"),
        whitelist = TRUE)
    } else {
      if (rlang::has_name(markers.meta, "FILTERS")) {
        markers.meta %<>% dplyr::filter(FILTERS == "whitelist")
      }
      markers.meta %<>% dplyr::select(VARIANT_ID, MARKERS)
    }

    # blacklist markers
    if (!is.null(blacklist.markers)) {
      markers.meta %<>% dplyr::filter(!MARKERS %in% blacklist.markers$MARKERS)
      sync_gds(gds = data, markers = markers.meta$VARIANT_ID)
    }
    n.markers <- nrow(markers.meta)
    ad <- extract_coverage(gds = data, ind = FALSE)

    if (is.null(ad$ref.tot) || is.null(ad$alt.tot)) {
      ad <- NULL
    } else {
      ad <- tibble::tibble(R_DEPTH = ad$ref.tot, A_DEPTH = ad$alt.tot)
    }

    mac.data <- suppressWarnings(
      dplyr::bind_cols(
        markers.meta,
        #MAC
        SeqArray::seqAlleleCount(
          gdsfile = data,
          ref.allele = NULL,
          .progress = TRUE,
          parallel = parallel.core) %>%
          unlist(.) %>%
          matrix(
            data = .,
            nrow = n.markers, ncol = 2, byrow = TRUE,
            dimnames = list(rownames = markers.meta$VARIANT_ID,
                            colnames = c("REF_COUNT", "ALT_COUNT"))) %>%
          tibble::as_tibble(.)))

    if (!is.null(ad)) {
      suppressWarnings(
        mac.data %<>%
          dplyr::bind_cols(ad) %>%
          dplyr::mutate(
            MAC_GLOBAL = dplyr::if_else(ALT_COUNT < REF_COUNT, ALT_COUNT, REF_COUNT),
            MAF_GLOBAL = MAC_GLOBAL / (REF_COUNT + ALT_COUNT),
            ALT_CHECK = dplyr::if_else(MAC_GLOBAL == ALT_COUNT, "ALT", "REF"),
            ALT_COUNT = NULL,
            REF_DEPTH = dplyr::if_else(ALT_CHECK == "ALT", R_DEPTH, A_DEPTH),
            ALT_DEPTH = dplyr::if_else(ALT_CHECK == "ALT", A_DEPTH, R_DEPTH),
            R_DEPTH = NULL, A_DEPTH = NULL, ALT_CHECK = NULL
          ) %>%
          dplyr::arrange(MARKERS)
      )

    } else {
      suppressWarnings(
        mac.data %<>%
          dplyr::mutate(
            MAC_GLOBAL = dplyr::if_else(ALT_COUNT < REF_COUNT, ALT_COUNT, REF_COUNT),
            MAF_GLOBAL = MAC_GLOBAL / (REF_COUNT + ALT_COUNT),
            ALT_CHECK = dplyr::if_else(MAC_GLOBAL == ALT_COUNT, "ALT", "REF"),
            ALT_COUNT = NULL, ALT_CHECK = NULL
          ) %>%
          dplyr::arrange(MARKERS)
      )
    }
  }
  ad <- markers.meta <- NULL
  return(mac.data)
}#End compute_mac

#' @title mac_one
#' @description mac without parallel
#' @rdname mac_one
#' @keywords internal
#' @export
mac_one <- function(x) {
  if (tibble::has_name(x, "GT_BIN")) {
    mac.data <- x %>%
      dplyr::select(MARKERS, INDIVIDUALS, GT_BIN) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::summarise(
        PP = as.numeric(2 * length(GT_BIN[GT_BIN == 0])),
        PQ = as.numeric(length(GT_BIN[GT_BIN == 1])),
        QQ = as.numeric(2 * length(GT_BIN[GT_BIN == 2]))
      ) %>%
      dplyr::ungroup(.) %>%
      dplyr::mutate(
        PP = PP + PQ,
        QQ = QQ + PQ,
        PQ = NULL,
        MAC_GLOBAL = dplyr::if_else(PP < QQ, PP, QQ),
        PP = NULL, QQ = NULL)
  } else {
    mac.data <- x %>%
      dplyr::select(MARKERS, INDIVIDUALS, GT) %>%
      dplyr::mutate(
        A1 = stringi::stri_sub(GT, 1, 3),
        A2 = stringi::stri_sub(GT, 4,6)
      ) %>%
      dplyr::select(-GT) %>%
      tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS)) %>%
      dplyr::group_by(MARKERS, GT) %>%
      dplyr::tally(.) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::filter(n == min(n)) %>%
      dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
      dplyr::summarise(MAC_GLOBAL = n) %>%
      dplyr::ungroup(.) %>%
      dplyr::select(MARKERS, MAC_GLOBAL)
  }
  mac.data$MAC_GLOBAL %<>% as.integer(.)
  return(mac.data)
}#End mac_one

