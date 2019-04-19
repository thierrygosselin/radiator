#' @title Erase genotypes
#' @description This function uses the information in the vcf tidy to i) make a
#' blacklist of individual genotypes to erase based on coverage and
#' genotype likelihood thresholds and ii) erase those blacklisted genotypes from
#'  a tidy vcf file and a STACKS haplotypes file.
#' @param tidy.vcf.file A data frame object or file (using ".tsv") of a tidy vcf.
#' @param haplotypes.file (optional) The 'batch_x.haplotypes.tsv'. If you want to
#' erase the genotypes that don't pass the threshold.
#' @param read.depth.threshold (integer) Threshold for the read depth.
#' @param allele.depth.threshold (integer) Threhold for the min depth
#' of REF or ALT alleles.
#' @param allele.imbalance.threshold (numeric) Threshold of ratio between ALT and REF
#' depth of coverage. See details below.
#' @param filename (optional) Name of the file written to the working directory.
#' @details Genotypes below average quality i.e. below threshold for the
#' coverage of REF and/or ALT allele and genotype likelihood
#' are zeroed from the file. The function erase SNP in the VCF file and loci
#' in the haplotypes file.
#' Also creates a blacklist of genotypes erased based on the genotype
#' likelihood threshold, the REF and ALT allele coverage threshold.
#' The ratio is calculated :
#' (read depth ALT allele - read depth REF allele)/(read depth ALT allele + read
#' depth REF allele).
#' e.g. REF = 3 and ALT = 2 the ratio = -0.2. For the function to work properly,
#' use positive values, the function will calculate the +/- imbalance.
#' @return The function returns the blacklisted individuals genotypes,
#' by loci, position (SNP, POS in stacks), populations and individuals.
#' For VCF, return the tidy vcf in the global environment only and in
#' the directory with \code{filename}. For haplotype file the original
#' filename with "_erased_geno" is appended.
#' @rdname erase_genotypes
#' @export

erase_genotypes <- function(tidy.vcf.file, haplotypes.file, read.depth.threshold, allele.depth.threshold, allele.imbalance.threshold, filename) {

  if (is.vector(tidy.vcf.file) == "TRUE") {
    tidy.vcf.file <- readr::read_tsv(tidy.vcf.file, col_names = T)
    message("Using the tidy vcf file in your directory")
  } else {
    tidy.vcf.file <- tidy.vcf.file
    message("Using the tidy vcf file from your global environment")
  }

  blacklist <- tidy.vcf.file %>%
    dplyr::filter(GT != "./." & GT != "0/0" & GT != "1/1" ) %>%
    dplyr::filter(READ_DEPTH == read.depth.threshold) %>%
    dplyr::filter(ALLELE_REF_DEPTH < allele.depth.threshold | ALLELE_ALT_DEPTH < allele.depth.threshold) %>%
    dplyr::filter(ALLELE_COVERAGE_RATIO < -allele.imbalance.threshold | ALLELE_COVERAGE_RATIO > allele.imbalance.threshold) %>%
    dplyr::select(LOCUS, POS, POP_ID, INDIVIDUALS) %>%
    dplyr::arrange(LOCUS, POS, POP_ID, INDIVIDUALS)

  readr::write_tsv(x = blacklist, path = "blacklist.genotypes.erased.txt", append = FALSE, col_names = T)

  # interesting stats.
  erased.genotype.number <- length(blacklist$INDIVIDUALS)
  total.genotype.number <- length(tidy.vcf.file$GT[tidy.vcf.file$GT != "./."])
  percent <- paste(round(((erased.genotype.number/total.genotype.number)*100), 2), "%", sep = " ")
  message.vcf.erasing.geno <- stringi::stri_join("Out of a total of ", total.genotype.number, " genotypes, ", percent, " (", erased.genotype.number, ")"," will be erased")
  message(message.vcf.erasing.geno)

  # VCF ------------------------------------------------------------------------
  message("Erasing... Erasing...")
  message("Step 1: Genotypes and Read Depth")

  tidy.vcf.file <- suppressWarnings(
    tidy.vcf.file %>%
      dplyr::mutate(
        GT = ifelse(LOCUS %in% blacklist$LOCUS & POS %in% blacklist$POS & INDIVIDUALS %in% blacklist$INDIVIDUALS, "./.", GT),
        READ_DEPTH = as.numeric(ifelse(LOCUS %in% blacklist$LOCUS & POS %in% blacklist$POS & INDIVIDUALS %in% blacklist$INDIVIDUALS, "NA", READ_DEPTH))
      )
  )

  message("Step 2: Reference and Alternate Alleles Depth")

  tidy.vcf.file <- suppressWarnings(
    tidy.vcf.file %>%
      dplyr::mutate(
        ALLELE_REF_DEPTH = as.numeric(ifelse(LOCUS %in% blacklist$LOCUS & POS %in% blacklist$POS & INDIVIDUALS %in% blacklist$INDIVIDUALS, "NA", ALLELE_REF_DEPTH)),
        ALLELE_ALT_DEPTH = as.numeric(ifelse(LOCUS %in% blacklist$LOCUS & POS %in% blacklist$POS & INDIVIDUALS %in% blacklist$INDIVIDUALS, "NA", ALLELE_ALT_DEPTH))
      )
  )

  message("Step 3: Genotype Likelihood")

  tidy.vcf.file <- suppressWarnings(
    tidy.vcf.file %>%
      dplyr::mutate(
        ALLELE_COVERAGE_RATIO = as.numeric(ifelse(LOCUS %in% blacklist$LOCUS & POS %in% blacklist$POS & INDIVIDUALS %in% blacklist$INDIVIDUALS, "NA", ALLELE_COVERAGE_RATIO)),
        GL = as.numeric(ifelse(LOCUS %in% blacklist$LOCUS & POS %in% blacklist$POS & INDIVIDUALS %in% blacklist$INDIVIDUALS, "NA", GL))
      )
  )
  if (missing(filename)) {
    message("Saving the tidy vcf not selected")
  } else {
    message("Writing the file to your working directory, this may take some time...")
    readr::write_tsv(x = tidy.vcf.file, path = filename, append = FALSE, col_names = TRUE)
  }

  # Haplotype file -------------------------------------------------------------
  if (missing(haplotypes.file)) {
    message("STACKS haplotypes file not provided")
  } else {
    message("Using the STACKS haplotypes file in your directory")
    haplo <- readr::read_tsv(haplotypes.file, col_names = TRUE) %>%
      dplyr::rename(LOCUS =`Catalog ID`)

    # haplotypes file preparation
    haplo.prep <- haplo %>%
      tidyr::gather(INDIVIDUALS, HAPLOTYPES, -c(LOCUS, Cnt)) %>%
      dplyr::mutate(INDIVIDUALS = as.character(INDIVIDUALS))

    # consensus
    consensus.pop <- haplo.prep %>%
      dplyr::mutate(CONSENSUS = stringi::stri_count_fixed(HAPLOTYPES, "consensus")) %>%
      dplyr::group_by(LOCUS) %>%
      dplyr::summarise(CONSENSUS_MAX = max(CONSENSUS)) %>%
      dplyr::filter(CONSENSUS_MAX > 0) %>%
      dplyr::select(LOCUS)

    haplo.number <- haplo.prep %>%
      dplyr::filter(HAPLOTYPES != "-") %>%
      dplyr::select(HAPLOTYPES)

    total.genotype.number.haplo <- length(haplo.number$HAPLOTYPES)
    percent.haplo <- paste(round(((erased.genotype.number/total.genotype.number.haplo)*100), 2), "%", sep = " ")
    message.haplo.erasing.geno <- stringi::stri_join("Out of a total of ", total.genotype.number.haplo, " genotypes, ", percent.haplo, " (", erased.genotype.number, ")"," will be erased")
    message(message.haplo.erasing.geno)
    message("Erasing... Erasing...")

    # Erasing genotype with the blacklist
    erase <- blacklist %>%
      dplyr::select(LOCUS, INDIVIDUALS) %>%
      dplyr::left_join(haplo.prep, by = c("LOCUS", "INDIVIDUALS")) %>%
      dplyr::mutate(HAPLOTYPES = rep("-", n()))

    keep <- haplo.prep %>%
      dplyr::anti_join(consensus.pop, by = "LOCUS") %>%
      dplyr::anti_join(
        blacklist %>%
          dplyr::select(LOCUS, INDIVIDUALS),
        by = c("LOCUS", "INDIVIDUALS")
      )

    haplo.erased <- dplyr::bind_rows(erase, keep) %>%
      dplyr::arrange(LOCUS, INDIVIDUALS) %>%
      dplyr::rename(`Catalog ID` = LOCUS) %>%
      tidyr::spread(INDIVIDUALS, HAPLOTYPES)

    haplo.erase.filename <- stringi::stri_replace_all_fixed(str = haplotypes.file, pattern = ".tsv", replacement = "_erased_geno.tsv", vectorize_all = FALSE)
    message("Saving the modified haplotypes file in your working directory")
    readr::write_tsv(haplo.erased, haplo.erase.filename, append = FALSE, col_names = TRUE)
  }

  res <- list()
  res$blacklist.genotypes <- blacklist
  res$tidy.vcf.erased.geno <- tidy.vcf.file

  invisible(cat(sprintf(
    "In your global environment, a list with the blacklisted genotypes and the tidy vcf"
  )))
  return(res)
}
