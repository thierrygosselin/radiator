# Detect if markers are biallelic

#' @name allele_frequencies

#' @title Compute allele frequencies per markers and populations
#' @description Compute allele frequencies per markers and populations.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.

#' @param verbose (optional, logical) \code{verbose = TRUE} to be chatty
#' during execution.
#' Default: \code{verbose = TRUE}.
#' @inheritParams radiator_common_arguments

#' @return A list with allele frequencies in a data frame in long and wide format,
#' and a matrix. Local (pop) and global minor allele frequency (MAF) is also computed.

#' @export
#' @rdname allele_frequencies
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

allele_frequencies <- function(data, verbose = TRUE, parallel.core = parallel::detectCores() - 1) {

  if (verbose) {
    cat("#######################################################################\n")
    cat("#################### radiator::allele_frequencies #####################\n")
    cat("#######################################################################\n")
    timing <- proc.time()
    message("Calculating allele frequencies...")
  }
  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file missing")

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    input <- radiator::tidy_wide(data = data, import.metadata = TRUE)
  } else {
    input <- data
  }

  # check genotype column naming
  if (rlang::has_name(input, "GENOTYPE")) {
    colnames(input) <- stringi::stri_replace_all_fixed(
      str = colnames(input),
      pattern = "GENOTYPE",
      replacement = "GT",
      vectorize_all = FALSE)
  }

  # necessary steps to make sure we work with unique markers and not duplicated LOCUS
  if (rlang::has_name(input, "LOCUS") && !rlang::has_name(input, "MARKERS")) {
    input <- dplyr::rename(.data = input, MARKERS = LOCUS)
  }

  if (rlang::has_name(input, "CHROM")) {
    metadata.markers <- dplyr::distinct(input, MARKERS, CHROM, LOCUS, POS)
    metadata <- TRUE
  } else {
    metadata <- FALSE
  }

  biallelic <- radiator::detect_biallelic_markers(input, verbose = verbose)
  markers.df <- dplyr::distinct(input, MARKERS)
  n.markers <- nrow(markers.df)
  maf.data <- dplyr::filter(input, GT != "000000")

  if (n.markers > 10000) {
    split.vec <- markers.df %>%
      dplyr::mutate(SPLIT_VEC = split_vec_row(
        markers.df,
        cpu.rounds = ceiling(n.markers/10000),
        parallel.core = parallel::detectCores() - 1))

    maf.data <- dplyr::left_join(maf.data, split.vec, by = "MARKERS") %>%
      radiator_future(
        X = .,
        FUN = compute_maf,
        parallel.core = parallel.core,
        split.tibble = .$SPLIT_VEC,
        bind.rows = TRUE,
        biallelic = biallelic
      )
      #
      # split(x = ., f = .$SPLIT_VEC) %>%
      # radiator_parallel_mc(
      #   X = .,
      #   FUN = compute_maf,
      #   mc.cores = parallel::detectCores() - 1,
      #   biallelic = biallelic
      # ) %>%
      # dplyr::bind_rows(.)
    markers.df <- split.vec <- NULL
  } else {
    maf.data <- compute_maf(x = maf.data, biallelic = biallelic)
  }

  maf <- maf.data
  maf.data <- NULL
  maf.local.wide <- dplyr::select(.data = maf, MARKERS, POP_ID, MAF_LOCAL) %>%
    tidyr::pivot_wider(data = ., names_from = "MARKERS", values_from = "MAF_LOCAL")


  maf.global <- dplyr::distinct(.data = maf, MARKERS, MAF_GLOBAL)

  if (rlang::has_name(input, "GT_VCF")) {

    freq <- maf %>%
      dplyr::mutate(FREQ_REF = 1 - MAF_LOCAL) %>%
      dplyr::select(MARKERS, POP_ID, FREQ_REF, MAF_LOCAL, MAF_GLOBAL)

    freq.wide <- dplyr::ungroup(freq) %>%
      dplyr::select(MARKERS, POP_ID, REF = FREQ_REF, ALT = MAF_LOCAL) %>%
      tidyr::pivot_longer(
        data = .,
        cols = -c("POP_ID", "MARKERS"),
        names_to = "ALLELES",
        values_to = "FREQ"
      ) %>%
      dplyr::mutate(MARKERS_ALLELES = stringi::stri_join(MARKERS, ALLELES, sep = ".")) %>%
      dplyr::select(-MARKERS, -ALLELES) %>%
      dplyr::arrange(MARKERS_ALLELES, POP_ID) %>%
      dplyr::group_by(POP_ID) %>%
      tidyr::pivot_wider(data = ., names_from = "MARKERS_ALLELES", values_from = "FREQ")

    freq.mat <- suppressWarnings(
      freq.wide %>%
        tibble::remove_rownames(.data = .) %>%
        tibble::column_to_rownames(.data = ., var = "POP_ID") %>%
        as.matrix(.)
    )

    if (metadata) {
      freq <- dplyr::full_join(freq, metadata.markers, by = "MARKERS") %>%
        dplyr::select(MARKERS, CHROM, LOCUS, POS, POP_ID, FREQ_REF, MAF_LOCAL, MAF_GLOBAL)
    }

  } else {

    freq <- dplyr::select(.data = input, MARKERS, POP_ID, INDIVIDUALS, GT) %>%
      dplyr::filter(GT != "000000") %>%
      dplyr::mutate(
        A1 = stringi::stri_sub(GT, 1, 3),
        A2 = stringi::stri_sub(GT, 4,6)
      ) %>%
      dplyr::select(MARKERS, POP_ID, INDIVIDUALS, A1, A2) %>%
      tidyr::pivot_longer(
        data = .,
        cols = -c("POP_ID", "INDIVIDUALS", "MARKERS"),
        names_to = "ALLELES_GROUP",
        values_to = "ALLELES"
      ) %>%
      dplyr::group_by(MARKERS, ALLELES, POP_ID) %>%
      dplyr::tally(.) %>%
      dplyr::ungroup(.) %>%
      tidyr::complete(data = ., POP_ID, tidyr::nesting(MARKERS, ALLELES), fill = list(n = 0)) %>%
      dplyr::group_by(MARKERS, POP_ID) %>%
      dplyr::mutate(
        NAPL = sum(n),
        FREQ = n / NAPL # Frequency of alleles per pop and locus
      ) %>%
      # dplyr::group_by(MARKERS, ALLELES) %>%
      # dplyr::mutate(FREQ = sum(n) / sum(NAPL)) %>% #Frequency of alleles per locus
      dplyr::select(MARKERS, POP_ID, ALLELES, FREQ) %>%
      dplyr::arrange(MARKERS, POP_ID, ALLELES)

    freq.wide <- dplyr::ungroup(freq) %>%
      dplyr::select(MARKERS, ALLELES, POP_ID, FREQ) %>%
      dplyr::mutate(MARKERS_ALLELES = stringi::stri_join(MARKERS, ALLELES, sep = ".")) %>%
      dplyr::select(-MARKERS, -ALLELES) %>%
      dplyr::arrange(MARKERS_ALLELES, POP_ID) %>%
      dplyr::group_by(POP_ID) %>%
      tidyr::pivot_wider(data = ., names_from = "MARKERS_ALLELES", values_from = "FREQ")

    freq.mat <- suppressWarnings(
      freq.wide %>%
        tibble::remove_rownames(df = .) %>%
        tibble::column_to_rownames(df = ., var = "POP_ID") %>%
        as.matrix(.)
    )

    freq <- dplyr::full_join(freq, maf, by = c("MARKERS", "POP_ID"))

    if (metadata) {
      freq <- dplyr::full_join(freq, metadata.markers, by = "MARKERS") %>%
        dplyr::select(MARKERS, CHROM, LOCUS, POS, POP_ID, ALLELES, FREQ, MAF_LOCAL, MAF_GLOBAL)
    }

  }
  if (verbose) {
    message(stringi::stri_join("Computation time: ", round((proc.time() - timing)[[3]]), " sec"))
    cat("############################## completed ##############################\n")
  }
  res <- list(
    freq.long = freq,
    freq.wide = freq.wide,
    mat = freq.mat,
    maf.local.wide = maf.local.wide,
    maf.global = maf.global
  )
  return(res)
}
