# MAF module

#' @name radiator_maf_module

#' @title MAF computation from a tidy data frame

#' @description Compute the minor allele frequency (local and global) using
#' a genomic tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.
#'
#' @param data A biallelic genomic data set in the working directory or
#' object in the global environment in wide or long (tidy) formats.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.

#' @inheritParams tidy_genomic_data

#' @return A tidy data frame with 6 columns:
#' \code{MARKERS, INDIVIDUALS, REF, ALT, GT_VCF, GT_BIN}.
#' \code{GT_VCF}: the genotype in VCF format
#' \code{GT_BIN}: coding used internally to easily convert to genlight,
#' the coding \code{0, 1, 2, NA} stands for the number of ALT allele in the
#' genotype and \code{NA} for missing genotype.


#' @export
#' @rdname radiator_maf_module
#' @importFrom dplyr select distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join funs summarise_at
#' @importFrom stringi stri_replace_all_fixed stri_join stri_sub
#' @importFrom tibble has_name
#' @importFrom tidyr gather complete separate nesting

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

radiator_maf_module <- function(
  data,
  maf.thresholds = NULL,
  maf.pop.num.threshold = 1,
  maf.approach = "SNP",
  maf.operator = "OR",
  parallel.core = parallel::detectCores() - 1
) {

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file missing")

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    input <- radiator::tidy_wide(data = data, import.metadata = TRUE)
  } else {
    input <- data
  }

  # check genotype column naming
  if (tibble::has_name(input, "GENOTYPE")) {
    colnames(input) <- stringi::stri_replace_all_fixed(
      str = colnames(input),
      pattern = "GENOTYPE",
      replacement = "GT",
      vectorize_all = FALSE)
  }

  # necessary steps to make sure we work with unique markers and not duplicated LOCUS
  if (tibble::has_name(input, "LOCUS") && !tibble::has_name(input, "MARKERS")) {
    input <- dplyr::rename(.data = input, MARKERS = LOCUS)
  }

  # keep markers metadata
  if (tibble::has_name(input, "CHROM")) {
    markers.meta <- dplyr::distinct(input, MARKERS, CHROM, LOCUS, POS)
    markers.df <- dplyr::select(markers.meta, MARKERS)
  } else {
    markers.meta <- NULL
    markers.df <- dplyr::distinct(input, MARKERS)
  }

  # Minor Allele Frequency filter ----------------------------------------------
  # maf.thresholds <- c(0.05, 0.1) # test
  if (!is.null(maf.thresholds)) { # with MAF
    maf.local.threshold <- maf.thresholds[1]
    maf.global.threshold <- maf.thresholds[2]
    message("MAF filter: yes")
  }
  message("    Calculating Minor Allele Frequency...")
  if (tibble::has_name(input, "GT_VCF")) {
    maf.local <- input %>%
      dplyr::filter(GT_VCF != "./.") %>%
      dplyr::group_by(MARKERS, POP_ID, REF, ALT) %>%
      dplyr::summarise(
        N = as.numeric(n()),
        PQ = as.numeric(length(GT_VCF[GT_VCF == "1/0" | GT_VCF == "0/1"])),
        QQ = as.numeric(length(GT_VCF[GT_VCF == "1/1"]))
      ) %>%
      dplyr::mutate(MAF_LOCAL = ((QQ * 2) + PQ) / (2 * N))

    maf.global <- maf.local %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::summarise_at(.tbl = ., .vars = c("N", "PQ", "QQ"), dplyr::funs(sum)) %>%
      dplyr::mutate(MAF_GLOBAL = ((QQ * 2) + PQ) / (2 * N)) %>%
      dplyr::select(MARKERS, MAF_GLOBAL)

    maf.data <- maf.global %>%
      dplyr::left_join(maf.local, by = c("MARKERS")) %>%
      dplyr::select(MARKERS, POP_ID, MAF_LOCAL, MAF_GLOBAL)

    maf.local <- maf.global <- NULL
  } else {

    number.markers <- dplyr::n_distinct(input$MARKERS)
    # Optimizing cpu usage
    if (number.markers <= 2000) {
      round.cpu <- floor(number.markers / parallel.core)
    } else {
      round.cpu <- floor(number.markers / (1000 * parallel.core))
    }
    # as.integer is usually twice as light as numeric vector...
    markers.df$SPLIT_VEC <- as.integer(floor((parallel.core * round.cpu * (1:number.markers - 1) / number.markers) + 1))

    # We split the alleles here to prep for MAF
    maf_gt <- function(x) {
      res <- x %>%
        dplyr::mutate(
          A1 = stringi::stri_sub(GT, 1, 3),
          A2 = stringi::stri_sub(GT, 4,6)
        ) %>%
        dplyr::select(MARKERS, POP_ID, INDIVIDUALS, A1, A2) %>%
        tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS, POP_ID)) %>%
        dplyr::group_by(MARKERS, GT, POP_ID) %>%
        dplyr::tally(.) %>%
        dplyr::ungroup(.) %>%
        tidyr::complete(data = ., POP_ID, tidyr::nesting(MARKERS, GT), fill = list(n = 0)) %>%
        dplyr::rename(n.al.pop = n) %>%
        dplyr::arrange(MARKERS, GT) %>%
        dplyr::group_by(MARKERS, GT) %>%
        dplyr::mutate(n.al.tot = sum(n.al.pop)) %>%
        dplyr::group_by(MARKERS) %>%
        dplyr::mutate(MAF_GLOBAL = min(n.al.tot)/sum(n.al.pop)) %>%
        dplyr::group_by(MARKERS, POP_ID) %>%
        dplyr::mutate(MAF_LOCAL = n.al.pop/sum(n.al.pop)) %>%
        dplyr::arrange(MARKERS, POP_ID, GT) %>%
        dplyr::group_by(MARKERS, POP_ID) %>%
        dplyr::filter(n.al.pop == min(n.al.pop)) %>%
        dplyr::distinct(MARKERS, POP_ID, .keep_all = TRUE) %>%
        dplyr::select(MARKERS, POP_ID, MAF_LOCAL, MAF_GLOBAL)
      return(res)
    }#End maf_gt

    maf.data <- input %>%
      dplyr::filter(GT != "000000") %>%
      dplyr::select(MARKERS, POP_ID, INDIVIDUALS, GT) %>%
      dplyr::left_join(markers.df, by = "MARKERS") %>%
      split(x = ., f = .$SPLIT_VEC) %>%
      .radiator_parallel(
        X = .,
        FUN = maf_gt,
        mc.preschedule = FALSE,
        mc.silent = FALSE,
        mc.cleanup = TRUE,
        mc.cores = parallel.core
      ) %>%
      dplyr::bind_rows(.)
    markers.df <- round.cpu <- NULL
    }# end maf calculations

  if (!is.null(markers.meta)) {
    maf.data <- dplyr::left_join(maf.data, markers.meta, by = "MARKERS") %>%
      dplyr::select(MARKERS, CHROM, LOCUS, POS, POP_ID, MAF_GLOBAL, MAF_LOCAL)
  }

  readr::write_tsv(
    x = maf.data,
    path = "maf.data.tsv",
    col_names = TRUE,
    append = FALSE
  )
  message("    The MAF table was written in your folder")

  if (!is.null(maf.thresholds)) {
    if (maf.approach == "haplotype") {
      if (!tibble::has_name(input, "CHROM")) {
        stop("The haplotype approach for maf needs locus and snp info from vcf")
      }
      # vcf.maf <- tidyr::separate(
      #   data = maf.data,
      #   col = MARKERS,
      #   into = c("CHROM", "LOCUS", "POS"),
      #   sep = "__",
      #   remove = FALSE,
      #   extra = "warn"
      # )

      if (maf.operator == "OR") {
        vcf.maf <- maf.data %>%
          dplyr::group_by(LOCUS, POP_ID) %>%
          dplyr::summarise(
            MAF_GLOBAL = mean(MAF_GLOBAL, na.rm = TRUE),
            MAF_LOCAL = min(MAF_LOCAL, na.rm = TRUE)
          ) %>%
          dplyr::filter(MAF_LOCAL >= maf.local.threshold | MAF_GLOBAL >= maf.global.threshold) %>%
          dplyr::group_by(LOCUS) %>%
          dplyr::tally(.) %>%
          dplyr::filter(n >= maf.pop.num.threshold) %>%
          dplyr::select(LOCUS) %>%
          dplyr::left_join(input, by = "LOCUS") %>%
          dplyr::arrange(LOCUS, POP_ID)
      } else {# AND operator between local and global maf
        vcf.maf <- maf.data %>%
          dplyr::group_by(LOCUS, POP_ID) %>%
          dplyr::summarise(
            MAF_GLOBAL = mean(MAF_GLOBAL, na.rm = TRUE),
            MAF_LOCAL = min(MAF_LOCAL, na.rm = TRUE)
          ) %>%
          dplyr::filter(MAF_LOCAL >= maf.local.threshold & MAF_GLOBAL >= maf.global.threshold) %>%
          dplyr::group_by(LOCUS) %>%
          dplyr::tally(.) %>%
          dplyr::filter(n >= maf.pop.num.threshold) %>%
          dplyr::select(LOCUS) %>%
          dplyr::left_join(input, by = "LOCUS") %>%
          dplyr::arrange(LOCUS, POP_ID)
      }
    } # end maf haplotype approach

    if (maf.approach == "SNP") { # SNP approach
      if (maf.operator == "OR") {
        vcf.maf <- maf.data %>%
          dplyr::group_by(MARKERS, POP_ID) %>%
          dplyr::summarise(
            MAF_GLOBAL = mean(MAF_GLOBAL, na.rm = TRUE),
            MAF_LOCAL = min(MAF_LOCAL, na.rm = TRUE)
          ) %>%
          dplyr::filter(MAF_LOCAL >= maf.local.threshold | MAF_GLOBAL >= maf.global.threshold) %>%
          dplyr::group_by(MARKERS) %>%
          dplyr::tally(.) %>%
          dplyr::filter(n >= maf.pop.num.threshold) %>%
          dplyr::select(MARKERS) %>%
          dplyr::left_join(input, by = "MARKERS") %>%
          dplyr::arrange(MARKERS, POP_ID)
      } else {# AND operator between local and global maf
        vcf.maf <- maf.data %>%
          dplyr::group_by(MARKERS, POP_ID) %>%
          dplyr::summarise(
            MAF_GLOBAL = mean(MAF_GLOBAL, na.rm = TRUE),
            MAF_LOCAL = min(MAF_LOCAL, na.rm = TRUE)
          ) %>%
          dplyr::filter(MAF_LOCAL >= maf.local.threshold & MAF_GLOBAL >= maf.global.threshold) %>%
          dplyr::group_by(MARKERS) %>%
          dplyr::tally(.) %>%
          dplyr::filter(n >= maf.pop.num.threshold) %>%
          dplyr::select(MARKERS) %>%
          dplyr::left_join(input, by = "MARKERS") %>%
          dplyr::arrange(MARKERS, POP_ID)
      }
      if (!is.null(markers.meta)) {
        vcf.maf <- dplyr::left_join(vcf.maf, markers.meta, by = "MARKERS")
      }
    } # end maf snp approach

    #reorder columns
    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "REF", "ALT", "POP_ID", "INDIVIDUALS", "GT", "GT_VCF", "GT_BIN", "GT_VCF_NUC")
    vcf.maf <- suppressWarnings(dplyr::select(vcf.maf, dplyr::one_of(want), dplyr::everything()))

    markers.before <- dplyr::n_distinct(input$MARKERS)
    markers.after <- dplyr::n_distinct(vcf.maf$MARKERS)
    message("    The number of MARKERS before the MAF filters = ", markers.before)
    message("    The number of MARKERS removed by the MAF filters = ", markers.before - markers.after)
    message("    The number of MARKERS after the MAF filters = ", markers.after)
  } else {
    vcf.maf <- "filtering not selected"
  }#End maf filtering
  res <- list(input = vcf.maf, maf.data = maf.data)
  return(res)
}#End radiator_maf_module
