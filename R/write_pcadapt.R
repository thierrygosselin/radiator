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
#' Also integrated in this function:
#' \enumerate{
#' \item \emph{pop.select: } to choose population to include in the pcadapt file.
#' \item \emph{maf: } evaluate the impact of Minor Allele Frequency on
#' genome scans with radiator maf module related arguments.
#' \item \emph{snp.ld: } SNPs on the same locus should not be considered
#' independant (short linkage disequilibrium) and the snp.ld argument
#' integrated here give users the ability to adress this quickly before
#' running a genome scan.
#' \item by defaults only markers found in common between populations are used
#' \item by defaults monomorphic markers are automatically removed before
#' generating the pcadapt file.
#' }


#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.

#' @inheritParams tidy_genomic_data

#' @param filename (optional) The file name prefix for the pcadapt file
#' written to the working directory. With default: \code{filename = NULL},
#' the date and time is appended to \code{radiator_pcadapt_}.

#' @details \emph{Integrated filters:} Only markers found in common between
#' populations are used and monomorphic markers are automatically removed
#' before generating pcadapt file.


#' @return A pcadapt file is written in the working directory a genotype matrix
#' object is also generated in the global environment.

#' @export
#' @rdname write_pcadapt
#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join
#' @importFrom stringi stri_join stri_replace_all_fixed stri_sub
#' @importFrom tidyr spread unite complete nesting separate
#' @importFrom readr write_delim write_tsv write_file
#' @importFrom parallel detectCores

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
  snp.ld = NULL,
  maf.thresholds = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1
) {

  message("Generating pcadapt file...")
  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file is missing")

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    data <- radiator::tidy_wide(data = data, import.metadata = TRUE)
  }

  # necessary steps to make sure we work with unique markers and not duplicated LOCUS
  if (tibble::has_name(data, "LOCUS") && !tibble::has_name(data, "MARKERS")) {
    data <- dplyr::rename(.data = data, MARKERS = LOCUS)
  }

  # pop.select -----------------------------------------------------------------
  if (!is.null(pop.select)) {
    message("pop.select: ")
    data <- dplyr::filter(data, POP_ID %in% pop.select)
    if (is.factor(data$POP_ID)) data$POP_ID <- droplevels(data$POP_ID)
  }

  # Keeping common markers -----------------------------------------------------
  data <- radiator::keep_common_markers(data = data, verbose = TRUE)$input

  # Removing monomorphic markers -----------------------------------------------
  data <- radiator::discard_monomorphic_markers(data = data, verbose = TRUE)$input

  # detect biallelic markers ---------------------------------------------------
  biallelic <- radiator::detect_biallelic_markers(data = data)

  if (!biallelic) {
    stop("\npcadapt only work with biallelic dataset")
  }

  # Linkage disequilibrium -----------------------------------------------------
  if (!is.null(snp.ld)) {
    if (tibble::has_name(data, "LOCUS") && tibble::has_name(data, "POS")) {
      message("Short distance linkage disequilibrium pruning")
      message("    snp.ld: ", snp.ld)
      data <- radiator::snp_ld(data = data, snp.ld = snp.ld)
    }
  }


  # MAF ------------------------------------------------------------------------
  if (!is.null(maf.thresholds)) { # with MAF
    data <- radiator::filter_maf(
      data = data,
      interactive.filter = FALSE,
      maf.thresholds = maf.thresholds,
      parallel.core = parallel.core,
      verbose = FALSE)$tidy.filtered.maf
  } # End of MAF filters

  # Biallelic and GT_BIN -------------------------------------------------------

  n.ind <- dplyr::n_distinct(data$INDIVIDUALS)
  n.pop <- dplyr::n_distinct(data$POP_ID)
  n.markers <- dplyr::n_distinct(data$MARKERS)


  if (!tibble::has_name(data, "GT_BIN")) {
    data <- radiator::change_alleles(
      data = dplyr::select(data, MARKERS, INDIVIDUALS, POP_ID, GT),
      biallelic = TRUE,
      parallel.core = parallel.core, verbose = TRUE)$input
  }


  pop.string <- data %>%
    dplyr::distinct(POP_ID, INDIVIDUALS) %>%
    dplyr::arrange(POP_ID, INDIVIDUALS) %>%
    dplyr::select(POP_ID)

  pop.string <- pop.string$POP_ID

  data <- dplyr::select(data, MARKERS, INDIVIDUALS, POP_ID, GT_BIN) %>%
    dplyr::arrange(POP_ID, INDIVIDUALS, MARKERS) %>%
    dplyr::select(-POP_ID) %>%
    dplyr::mutate(GT_BIN = replace(GT_BIN, which(is.na(GT_BIN)), 9)) %>%
    dplyr::group_by(INDIVIDUALS) %>%
    tidyr::spread(data = ., key = INDIVIDUALS, value = GT_BIN) %>%
    dplyr::ungroup(.) %>%
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

  readr::write_delim(x = data, path = filename, col_names = FALSE,
                     append = FALSE, delim = " ")

  data <- as.matrix(data)
  res <- list(genotype.matrix = data, pop.string = pop.string)
  return(res)
}# End write_pcadapt


