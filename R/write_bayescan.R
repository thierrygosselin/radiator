# write a bayescan file from a tidy data frame

#' @name write_bayescan
#' @title Write a \href{http://cmpg.unibe.ch/software/BayeScan/}{BayeScan}
#' file from a tidy data frame

#' @description Write a \href{http://cmpg.unibe.ch/software/BayeScan/}{BayeScan}
#' file from a tidy data frame. The data is bi- or multi-allelic.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.

#' @inheritParams tidy_genomic_data

#' @param filename (optional) The file name prefix for the bayescan file
#' written to the working directory. With default: \code{filename = NULL},
#' the date and time is appended to \code{radiator_bayescan_}.

#' @return A bayescan file is written in the working directory.

#' @export
#' @rdname write_bayescan
#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join
#' @importFrom stringi stri_join stri_replace_all_fixed stri_sub
#' @importFrom purrr walk
#' @importFrom tidyr spread unite complete nesting separate
#' @importFrom readr write_delim write_tsv write_file
#' @importFrom parallel detectCores
#' @importFrom utils write.table

#' @references Foll, M and OE Gaggiotti (2008) A genome scan method to identify
#' selected loci appropriate
#' for both dominant and codominant markers: A Bayesian perspective.
#' Genetics 180: 977-993

#' @references Foll M, Fischer MC, Heckel G and L Excoffier (2010)
#' Estimating population structure from
#' AFLP amplification intensity. Molecular Ecology 19: 4638-4647

#' @references Fischer MC, Foll M, Excoffier L and G Heckel (2011) Enhanced AFLP
#' genome scans detect
#' local adaptation in high-altitude populations of a small rodent (Microtus arvalis).
#' Molecular Ecology 20: 1450-1462

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_bayescan <- function(
  data,
  pop.select = NULL,
  snp.ld = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1
  ) {

  message("Generating BayeScan file...")
  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file necessary to write the bayescan file is missing")

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    input <- radiator::tidy_wide(data = data, import.metadata = TRUE)
  } else {
    input <- data
  }

  # check genotype column naming
  colnames(input) <- stringi::stri_replace_all_fixed(
    str = colnames(input),
    pattern = "GENOTYPE",
    replacement = "GT",
    vectorize_all = FALSE
  )

  # necessary steps to make sure we work with unique markers and not duplicated LOCUS
  if (tibble::has_name(input, "LOCUS") && !tibble::has_name(input, "MARKERS")) {
    input <- dplyr::rename(.data = input, MARKERS = LOCUS)
  }

  # pop.select -----------------------------------------------------------------
  if (!is.null(pop.select)) {
    message("pop.select: ")
    input <- dplyr::filter(input, POP_ID %in% pop.select)
    input$POP_ID <- droplevels(input$POP_ID)
  }

  # Keeping common markers -----------------------------------------------------
  input <- radiator::keep_common_markers(data = input, verbose = TRUE)

  # Removing monomorphic markers -----------------------------------------------
  input <- radiator::discard_monomorphic_markers(data = input, verbose = TRUE)$input

  # detect biallelic markers ---------------------------------------------------
  biallelic <- radiator::detect_biallelic_markers(data = input)

  if (!biallelic) {
    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "INDIVIDUALS", "POP_ID", "GT_VCF_NUC", "GT")
    input <- suppressWarnings(dplyr::select(input, dplyr::one_of(want)))
    if (tibble::has_name(input, "GT_VCF_NUC")) {
      want <- c("MARKERS", "CHROM", "LOCUS", "POS", "INDIVIDUALS", "POP_ID", "GT_VCF_NUC")
      input <- suppressWarnings(dplyr::select(input, dplyr::one_of(want))) %>%
        dplyr::rename(GT_HAPLO = GT_VCF_NUC)
    } else {
      want <- c("MARKERS", "CHROM", "LOCUS", "POS", "INDIVIDUALS", "POP_ID", "GT")
      input <- suppressWarnings(dplyr::select(input, dplyr::one_of(want))) %>%
        dplyr::rename(GT_HAPLO = GT)
    }

    input <- change_alleles(
        biallelic = FALSE,
        data = input, monomorphic.out = TRUE,
        parallel.core = parallel.core, verbose = TRUE)$input
  }

  # Linkage disequilibrium -----------------------------------------------------
  if (!is.null(snp.ld)) {
    if (tibble::has_name(input, "LOCUS") && tibble::has_name(input, "POS")) {
      message("Short distance linkage disequilibrium pruning")
      message("    snp.ld: ", snp.ld)
      input <- radiator::snp_ld(data = input, snp.ld = snp.ld)
    }
  }

  # Biallelic and GT_BIN -------------------------------------------------------


  if (biallelic) {
    input <- dplyr::select(input, MARKERS, INDIVIDUALS, POP_ID, GT)
    input <- radiator::change_alleles(
      data = input,
      monomorphic.out = TRUE,
      biallelic = TRUE,
      parallel.core = parallel.core, verbose = TRUE)$input
    input <- dplyr::select(input, MARKERS, INDIVIDUALS, POP_ID, GT_BIN)
  } else {
    input <- input
  }

  # prep data wide format ------------------------------------------------------
  n.ind <- dplyr::n_distinct(input$INDIVIDUALS)
  n.pop <- dplyr::n_distinct(input$POP_ID)
  n.markers <- dplyr::n_distinct(input$MARKERS)
  # ind.per.pop <- dplyr::distinct(input, POP_ID, INDIVIDUALS) %>%
  #   dplyr::group_by(POP_ID) %>%
  #   dplyr::tally(.)

  input <- dplyr::ungroup(input) %>%
    dplyr::mutate(
      BAYESCAN_POP = POP_ID,
      BAYESCAN_POP = as.integer(BAYESCAN_POP),
      BAYESCAN_MARKERS = factor(MARKERS),
      BAYESCAN_MARKERS = as.integer(BAYESCAN_MARKERS)
    )

  pop.dictionary <- dplyr::distinct(input, POP_ID, BAYESCAN_POP)
  markers.dictionary <- dplyr::distinct(input, MARKERS, BAYESCAN_MARKERS) %>%
    dplyr::arrange(BAYESCAN_MARKERS)

  input <- dplyr::select(input, -POP_ID, -MARKERS)

  # writting file to directory  ------------------------------------------------
  # Filename: date and time to have unique filenaming
  file.date <- stringi::stri_replace_all_fixed(
    Sys.time(),
    pattern = " EDT", replacement = "") %>%
    stringi::stri_replace_all_fixed(
      str = .,
      pattern = c("-", " ", ":"), replacement = c("", "@", ""),
      vectorize_all = FALSE) %>%
    stringi::stri_sub(str = ., from = 1, to = 13)

  if (is.null(filename)) {
    filename <- stringi::stri_join("radiator_bayescan_", file.date, ".txt")
  } else {
    filename.problem <- file.exists(filename)
    if (filename.problem) {
      filename <- stringi::stri_join(filename, "_", file.date, ".txt")
    } else {
      filename <- stringi::stri_join(filename, ".txt")
    }
  }

  if (biallelic) {
    markers.type <- "biallelic"
  } else {
    markers.type <- "multiallelic"
  }

  message("Writting BayeScan file with:
    Number of populations: ", n.pop, "\n    Number of individuals: ", n.ind,
          "\n    Number of ", markers.type, " markers: ", n.markers)

  # Number of markers
  readr::write_file(x = stringi::stri_join("[loci]=", n.markers, "\n\n"), path = filename, append = FALSE)

  # Number of populations
  readr::write_file(x = stringi::stri_join("[populations]=", n.pop, "\n\n"), path = filename, append = TRUE)
  pop.string <- unique(input$BAYESCAN_POP)
  generate_bayescan_biallelic <- function(pop, data) {
    # pop <- "BEA"
    data.pop <- dplyr::filter(data, BAYESCAN_POP %in% pop) %>%
      dplyr::filter(!is.na(GT_BIN)) %>%
      dplyr::group_by(BAYESCAN_MARKERS) %>%
      dplyr::summarise(
        REF = (length(GT_BIN[GT_BIN == 0]) * 2) + (length(GT_BIN[GT_BIN == 1])),
        ALT = (length(GT_BIN[GT_BIN == 2]) * 2) + (length(GT_BIN[GT_BIN == 1]))
      ) %>%
      dplyr::mutate(GENE_N = REF + ALT, ALLELE_N = rep(2, n())) %>%
      dplyr::select(BAYESCAN_MARKERS, GENE_N, ALLELE_N, REF, ALT)
    readr::write_file(x = stringi::stri_join("[pop]=", pop, "\n"), path = filename, append = TRUE)
    readr::write_delim(x = data.pop, path = filename, append = TRUE, delim = "  ")
    readr::write_file(x = stringi::stri_join("\n"), path = filename, append = TRUE)
  }
  generate_bayescan_multiallelic <- function(data) {
    pop <- unique(data$BAYESCAN_POP)
    data.pop <- dplyr::select(data, -BAYESCAN_POP)
    readr::write_file(x = stringi::stri_join("[pop]=", pop, "\n"), path = filename, append = TRUE)
    # readr::write_delim(x = data.pop, path = filename, append = TRUE, delim = "  " )
    utils::write.table(x = data.pop, file = filename, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
    readr::write_file(x = stringi::stri_join("\n"), path = filename, append = TRUE)
  }

  if (!biallelic) {
    data.prep <- input %>%
      dplyr::select(GT_VCF, BAYESCAN_MARKERS, BAYESCAN_POP) %>%
      dplyr::filter(GT_VCF != "./.") %>%
      tidyr::separate(data = ., col = GT_VCF, into = c("A1", "A2"), sep = "/") %>%
      tidyr::gather(data = ., key = ALLELES_GROUP, value = ALLELES, -c(BAYESCAN_MARKERS, BAYESCAN_POP)) %>%
      dplyr::select(-ALLELES_GROUP)

    allele.count <- data.prep %>%
      dplyr::distinct(BAYESCAN_MARKERS, ALLELES) %>%
      dplyr::group_by(BAYESCAN_MARKERS) %>%
      dplyr::tally(.) %>%
      dplyr::rename(COUNT = n)

    data.prep <- data.prep %>%
      dplyr::group_by(BAYESCAN_MARKERS, BAYESCAN_POP, ALLELES) %>%
      dplyr::tally(.) %>%
      dplyr::ungroup(.) %>%
      tidyr::complete(data = ., BAYESCAN_POP, tidyr::nesting(BAYESCAN_MARKERS, ALLELES), fill = list(n = 0)) %>%
      dplyr::ungroup(.)

    alleles.markers <- data.prep %>%
      dplyr::group_by(BAYESCAN_MARKERS, BAYESCAN_POP) %>%
      dplyr::summarise(GENE_N = sum(n)) %>%
      dplyr::ungroup(.) %>%
      dplyr::left_join(allele.count, by = "BAYESCAN_MARKERS") %>%
      dplyr::group_by(BAYESCAN_MARKERS, BAYESCAN_POP) %>%
      dplyr::summarise(GENE_N = stringi::stri_join(GENE_N, COUNT, sep = " "))

    data.prep <- data.prep %>%
      dplyr::arrange(BAYESCAN_MARKERS, BAYESCAN_POP, ALLELES) %>%
      dplyr::group_by(BAYESCAN_MARKERS, BAYESCAN_POP) %>%
      dplyr::summarise(ALLELES = stringi::stri_join(n, collapse = " ")) %>%
      dplyr::arrange(BAYESCAN_MARKERS, BAYESCAN_POP)

    input <- dplyr::left_join(
      alleles.markers, data.prep, by = c("BAYESCAN_MARKERS", "BAYESCAN_POP")) %>%
      dplyr::arrange(BAYESCAN_MARKERS, BAYESCAN_POP) %>%
      dplyr::group_by(BAYESCAN_MARKERS, BAYESCAN_POP) %>%
      dplyr::summarise(GT = stringi::stri_join(GENE_N, ALLELES, sep = " ")) %>%
      dplyr::ungroup(.) %>%
      dplyr::arrange(BAYESCAN_MARKERS, BAYESCAN_POP) %>%
      split(x = ., f = .$BAYESCAN_POP)
    data.prep <- alleles.markers <- allele.count <- NULL

    purrr::walk(.x = input, .f = generate_bayescan_multiallelic)
  } else {
    purrr::walk(.x = pop.string, .f = generate_bayescan_biallelic, data = input)
  }


  message("Writting populations dictionary")
  readr::write_tsv(
    x = pop.dictionary,
    path = stringi::stri_replace_all_fixed(
      str = filename, pattern = ".txt",
      replacement = "_pop_dictionary.tsv", vectorize_all = FALSE))
  message("Writting markers dictionary")
  readr::write_tsv(
    x = markers.dictionary,
    path = stringi::stri_replace_all_fixed(
      str = filename, pattern = ".txt",
      replacement = "_markers_dictionary.tsv", vectorize_all = FALSE))

  res <- list(pop.dictionary = pop.dictionary, markers.dictionary = markers.dictionary)
  return(res)
}# End write_bayescan


