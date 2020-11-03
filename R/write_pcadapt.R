# write a pcadapt file from a tidy data frame

#' @name write_pcadapt
#' @title Write a \href{https://github.com/bcm-uga/pcadapt}{pcadapt}
#' file from a tidy data frame

#' @description Write a
#' \href{https://github.com/bcm-uga/pcadapt}{pcadapt}
#' file from a tidy data frame. The data is biallelic.
#' Used internally in
#' \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.
#'

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.

#' @inheritParams read_strata
#' @inheritParams tidy_genomic_data

#' @param filename (optional) The file name prefix for the pcadapt file
#' written to the working directory. With default: \code{filename = NULL},
#' the date and time is appended to \code{radiator_pcadapt_}.

#' @details \emph{Integrated filters:}
#' \enumerate{
#' \item by defaults only markers found in common between populations are used
#' (See advance section).
#' \item by defaults monomorphic markers are automatically removed before
#' generating the pcadapt file.
#' }

#' @section Advance mode:
#'
#' \emph{dots-dots-dots ...} allows to pass several arguments for fine-tuning the function:
#' \enumerate{
#' \item \strong{Filtering for linkage disequilibrium}: 3 arguments
#' \code{filter.long.ld, long.ld.missing, ld.method}
#' described in \code{\link{filter_ld}} are available.
#' Reducing linkage before running genome scan is essential. At least start by
#' removing SNPs on the same RADseq locus (short linkage disequilibrium).
#'
#' \item \strong{Filtering markers with low Minor Allele Count} : use the argument
#' \code{filter.mac} to evaluate the impact of MAC/MAF on genome scans.
#' The function \code{\link{filter_mac}} is called.
#'
#' \item Turning off the filter that keeps markers in common between strata:
#' This is not recommended, but users who wants to explore the impact of such filtering
#' and know the biais it can potentially generate can use the argument
#' \code{filter.common.markers}.
#' The function \code{\link{filter_common_markers}} is called.
#' Default: \code{filter.common.markers = NULL}
#' }

#' @return A pcadapt file is written in the working directory a genotype matrix
#' object is also generated in the global environment.

#' @export
#' @rdname write_pcadapt
#' @references Luu, K., Bazin, E., & Blum, M. G. (2017).
#' pcadapt: an R package to perform genome scans for selection based on principal component analysis.
#' Molecular Ecology Resources, 17(1), 67-77.

#' @references Duforet-Frebourg, N., Luu, K., Laval, G., Bazin, E., & Blum, M. G. (2015).
#' Detecting genomic signatures of natural selection with principal component analysis: application to the 1000 Genomes data.
#' Molecular biology and evolution, msv334.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_pcadapt <- function(
  data,
  pop.select = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  ...
) {

  message("Generating pcadapt file...")
  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("Input file is missing")

  # dotslist -------------------------------------------------------------------
  dotslist <- rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE)
  want <- c("filter.common.markers", "filter.mac", "filter.short.ld",
            "filter.long.ld", "long.ld.missing", "ld.method")

  unknowned_param <- setdiff(names(dotslist), want)

  if (length(unknowned_param) > 0) {
    rlang::abort("Unknowned \"...\" parameters ",
         stringi::stri_join(unknowned_param, collapse = " "))
  }

  radiator.dots <- dotslist[names(dotslist) %in% want]


  # argument <- radiator.dots[["argument"]]
  filter.mac <- radiator.dots[["filter.mac"]]
  filter.short.ld <- radiator.dots[["filter.short.ld"]]
  filter.long.ld <- radiator.dots[["filter.long.ld"]]
  long.ld.missing <- radiator.dots[["long.ld.missing"]]
  if (is.null(long.ld.missing)) long.ld.missing <- FALSE
  ld.method <- radiator.dots[["ld.method"]]
  if (is.null(ld.method)) ld.method <- "r2"
  filter.common.markers <- radiator.dots[["filter.common.markers"]]
  if (is.null(filter.common.markers)) filter.common.markers <- TRUE
  if (!filter.common.markers) {
    message("Not recommended: I hope you really know what you're doing with pcadapt")
  }

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) data %<>% radiator::tidy_wide(data = ., import.metadata = TRUE)

  # pop.select -----------------------------------------------------------------
  if (!is.null(pop.select)) {
    message("pop.select: ")
    data %<>% dplyr::filter(POP_ID %in% pop.select)
    if (is.factor(data$POP_ID)) data$POP_ID <- droplevels(data$POP_ID)
  }

  # Keeping common markers -----------------------------------------------------
  if (filter.common.markers) {
    data <- radiator::filter_common_markers(data = data, verbose = TRUE, internal = TRUE)
  }

  # Removing monomorphic markers -----------------------------------------------
  data <- radiator::filter_monomorphic(data = data, verbose = TRUE, internal = TRUE)

  # detect biallelic markers ---------------------------------------------------
  biallelic <- radiator::detect_biallelic_markers(data = data)

  if (!biallelic) rlang::abort("\npcadapt only work with biallelic dataset")

  # MAC ------------------------------------------------------------------------
  if (!is.null(filter.mac)) { # with MAF
    data <- radiator::filter_mac(
      data = data,
      interactive.filter = FALSE,
      filter.mac = filter.mac,
      parallel.core = parallel.core,
      verbose = FALSE) %$% input
  } # End of MAC filters

  # Linkage disequilibrium -----------------------------------------------------
  if (!is.null(filter.short.ld) || !is.null(filter.long.ld)) {
    data <- filter_ld(
      data = data,
      interactive.filter = FALSE,
      filter.short.ld = filter.short.ld,
      filter.long.ld = filter.long.ld,
      long.ld.missing = long.ld.missing,
      ld.method = ld.method
    )
  }

  # Biallelic and GT_BIN -------------------------------------------------------

  n.ind <- dplyr::n_distinct(data$INDIVIDUALS)
  n.pop <- dplyr::n_distinct(data$POP_ID)
  n.markers <- dplyr::n_distinct(data$MARKERS)


  if (!rlang::has_name(data, "GT_BIN")) {
    data %<>% radiator::calibrate_alleles(
      data = ., biallelic = TRUE, gt.bin = TRUE) %$% input
  }

  data  %<>% dplyr::select(MARKERS, INDIVIDUALS, POP_ID, GT_BIN)

  pop.string <- data %>%
    dplyr::distinct(POP_ID, INDIVIDUALS) %>%
    dplyr::arrange(POP_ID, INDIVIDUALS) %>%
    dplyr::select(POP_ID)

  pop.string <- pop.string$POP_ID

  data %<>%
    dplyr::select(MARKERS, INDIVIDUALS, POP_ID, GT_BIN) %>%
    dplyr::arrange(POP_ID, INDIVIDUALS, MARKERS) %>%
    dplyr::select(-POP_ID) %>%
    dplyr::mutate(GT_BIN = replace(GT_BIN, which(is.na(GT_BIN)), 9)) %>%
    rad_wide(x = ., formula = "MARKERS ~ INDIVIDUALS", values_from = "GT_BIN") %>% # could be the other way ...
    dplyr::select(-MARKERS)

  # writing file to directory  ------------------------------------------------
  # Filename: date and time to have unique filenaming
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")

  if (is.null(filename)) {
    filename <- stringi::stri_join("radiator_pcadapt_", file.date, ".txt")
  } else {
    filename.problem <- file.exists(filename)
    if (filename.problem) {
      filename <- stringi::stri_join(filename, "_pcadapt_", file.date, ".txt")
    } else {
      filename <- stringi::stri_join(filename, "_pcadapt", ".txt")
    }
  }

  message("writing pcadapt file with:
    Number of populations: ", n.pop, "\n    Number of individuals: ", n.ind,
          "\n    Number of markers: ", n.markers)

  readr::write_delim(x = data, file = filename, col_names = FALSE,
                     append = FALSE, delim = " ")

  data <- as.matrix(data)
  res <- list(genotype.matrix = data, pop.string = pop.string)
  return(res)
}# End write_pcadapt


