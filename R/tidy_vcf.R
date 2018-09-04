#' @name tidy_vcf

#' @title Tidy a vcf file (bi/multi-allelic).

#' @description Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.
#' Highly recommended to use \code{\link[radiator]{tidy_genomic_data}} or
#' \code{\link[radiator]{genomic_converter}}. Those two functions allows
#' to manipulate and prune the dataset with blacklists and whitelists along
#' several other filtering options.

#' @inheritParams tidy_genomic_data

#' @param ... (optional) To pass further argument for fine-tuning the tidying
#' (details below).

#' @export
#' @rdname tidy_vcf
# @importFrom vcfR read.vcfR extract.gt vcf_field_names
#' @importFrom rlang UQ
#' @importFrom data.table melt.data.table as.data.table

#' @return The output in your global environment is a tidy data frame.

#' @details
#' \strong{... :dot dot dot arguments}
#' 5 arguments are available:
#' \enumerate{
#' \item whitelist.markers
#' \item blacklist.id
#' \item pop.select
#' \item pop.levels
#' \item pop.labels
#' \item ref.calibration
#' \item vcf.stats
#' }
#' Documentation for these arguments is detailed
#' in \code{\link[radiator]{tidy_genomic_data}}, the only difference here, the
#' arguments are objects are in the global environment, not as file in the directory.
#' If files are essential in your pipeline, please use \code{\link[radiator]{tidy_genomic_data}}.


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

tidy_vcf <- function(
  data,
  strata = NULL,
  vcf.metadata = FALSE,
  parallel.core = parallel::detectCores() - 1,
  verbose = FALSE,
  ...) {

  # # test
  # data = "test.vcf" # data = "populations.snps.vcf"
  # strata = "StLaw_popmap_thierry.tsv"
  # vcf.metadata = TRUE
  # parallel.core = 8
  # verbose = TRUE
  # blacklist.id = "blacklist.id.missing.50.tsv"
  # whitelist.markers = NULL
  # pop.select = NULL
  # pop.levels = NULL
  # pop.labels = NULL
  # filename = NULL
  # vcf.stats = TRUE
  # snp.read.position.filter = c("outliers", "q75", "iqr")
  # mac.threshold = 4
  # gt.vcf.nuc = TRUE
  # gt.vcf = TRUE
  # gt = TRUE
  # gt.bin = TRUE
  # wide = FALSE
  # ref.calibration = FALSE
  # keep.gds = TRUE

  opt.change <- getOption("width")
  options(width = 70)
  timing <- proc.time()

  # Note to myself: have to integrate this...
  wide <- FALSE


  # required packages ----------------------------------------------------------
  # Check that SeqArray is installed
  if (!"SeqArray" %in% utils::installed.packages()[,"Package"]) {
    stop('Please install SeqArray for this option:\n
         devtools::install_github("zhengxwen/SeqArray")')
  }

  if (!"SeqVarTools" %in% utils::installed.packages()[,"Package"]) {
    stop('Please install SeqVarTools for this option:\n
         source("https://bioconductor.org/biocLite.R")
         biocLite("SeqVarTools")')
  }

  if (!"gdsfmt" %in% utils::installed.packages()[,"Package"]) {
    stop('Please install gdsfmt for this option:\n
         source("https://bioconductor.org/biocLite.R")
         biocLite("gdsfmt")')
  }

  # dotslist -------------------------------------------------------------------
  dotslist <- list(...)
  want <- c("whitelist.markers", "blacklist.id", "pop.select",
            "pop.levels", "pop.labels", "snp.read.position.filter", "mac.threshold",
            "ref.calibration", "gt.vcf.nuc", "gt.vcf", "gt", "gt.bin", "vcf.stats",
            "filename", "keep.gds")
  unknowned_param <- setdiff(names(dotslist), want)

  if (length(unknowned_param) > 0) {
    stop("Unknowned \"...\" parameters ",
         stringi::stri_join(unknowned_param, collapse = " "))
  }

  radiator.dots <- dotslist[names(dotslist) %in% want]
  whitelist.markers <- radiator.dots[["whitelist.markers"]]
  blacklist.id <- radiator.dots[["blacklist.id"]]
  pop.select <- radiator.dots[["pop.select"]]
  pop.levels <- radiator.dots[["pop.levels"]]
  pop.labels <- radiator.dots[["pop.labels"]]
  snp.read.position.filter <- radiator.dots[["snp.read.position.filter"]]
  mac.threshold <- radiator.dots[["mac.threshold"]]
  ref.calibration <- radiator.dots[["ref.calibration"]]
  gt.vcf.nuc <- radiator.dots[["gt.vcf.nuc"]]
  gt.vcf <- radiator.dots[["gt.vcf"]]
  gt <- radiator.dots[["gt"]]
  gt.bin <- radiator.dots[["gt.bin"]]
  vcf.stats <- radiator.dots[["vcf.stats"]]
  filename <- radiator.dots[["filename"]]
  keep.gds <- radiator.dots[["keep.gds"]]

  if (is.null(keep.gds)) keep.gds <- TRUE
  if (is.null(vcf.stats)) vcf.stats <- TRUE
  if (is.null(ref.calibration)) ref.calibration <- FALSE
  if (is.null(gt.vcf.nuc)) gt.vcf.nuc <- TRUE
  if (is.null(gt.vcf)) gt.vcf <- TRUE
  if (is.null(gt)) gt <- TRUE
  if (is.null(gt.bin)) gt.bin <- TRUE


  if (!gt.vcf.nuc && !gt) {
    stop("At least one of gt.vcf.nuc or gt must be TRUE")
  }

  if (!is.null(snp.read.position.filter)) {
    snp.read.position.filter <- match.arg(
      arg = snp.read.position.filter,
      choices = c("outliers", "iqr", "q75"),
      several.ok = TRUE)
  }

  # import VCF -----------------------------------------------------------------
  data <- radiator::write_seqarray(
    vcf = data,
    strata = strata,
    filename = filename,
    vcf.stats = vcf.stats,
    blacklist.id = blacklist.id,
    snp.read.position.filter = snp.read.position.filter,
    mac.threshold = mac.threshold,
    whitelist.markers = whitelist.markers,
    verbose = TRUE,
    parallel.core = parallel.core,
    keep.gds = keep.gds)

  # Number of markers and tidy approach ----------------------------------------
  if (data$n.markers > 50000) {
    cat("\n\n############################# IMPORTANT ###############################\n")
    message("Tidying vcf with ", data$n.markers, " SNPs is not optimal")
    message("use radiator::filter_rad to reduce to ~ 10 000 unlinked SNPs\n\n")
  }

  # re-calibration of ref/alt alleles ------------------------------------------
  if (!is.null(blacklist.id)) {
    ref.calibration <- TRUE
    if (verbose) message("\nRe-calibration of REF/ALT alleles is required...")
    #overide...
    gt.vcf.nuc.bk <- gt.vcf.nuc
    gt.bk <- gt
    gt.vcf.bk <- gt.vcf
    gt.bin.bk <- gt.bin

    gt.vcf.nuc <- TRUE
    gt <- FALSE
    gt.vcf <- FALSE
    gt.bin <- FALSE
  }

  # bi- or multi-alllelic VCF --------------------------------------------------
  biallelic <- data$biallelic

  # import genotypes -----------------------------------------------------------
  # nucleotides info required to generate genepop format 001, 002, 003, 004
  data$genotypes$GT_VCF_NUC <- SeqVarTools::getGenotypeAlleles(
    gdsobj = data$vcf.connection, use.names = TRUE)
  if (gt) {
    if (!wide) {
      if (!gt.vcf.nuc) {
        data$genotypes$GT <- as.vector(data$genotypes$GT_VCF_NUC)
        data$genotypes$GT_VCF_NUC <- NULL
      } else {
        data$genotypes$GT <- data$genotypes$GT_VCF_NUC <- as.vector(data$genotypes$GT_VCF_NUC)
      }
    } else {
      data$genotypes$GT <- as.vector(data$genotypes$GT_VCF_NUC)
      if (!gt.vcf.nuc) data$genotypes$GT_VCF_NUC <- NULL
    }
    data$genotypes$GT <- stringi::stri_replace_all_fixed(
      str = data$genotypes$GT,
      pattern = c("A", "C", "G", "T", "/"),
      replacement = c("001", "002", "003", "004", ""),
      vectorize_all = FALSE) %>%
      stringi::stri_replace_na(str = ., replacement = "000000")

    if (wide) {
      data$genotypes$GT <- tibble::as_tibble(
        matrix(data = data$genotypes$GT,
               nrow = nrow(data$individuals),
               ncol = nrow(data$markers.meta))
      ) %>%
        magrittr::set_colnames(x = ., value = data$markers.meta$MARKERS) %>%
        tibble::add_column(
          .data = .,
          "INDIVIDUALS" = data$individuals$INDIVIDUALS,
          .before = 1
        )
    }
  }

  # more work for gt.vcf.nuc
  if (gt.vcf.nuc) {
    if (wide) {
      data$genotypes$GT_VCF_NUC %<>%
        magrittr::set_colnames(x = ., value = data$markers.meta$MARKERS) %>%
        tibble::as_tibble(x = ., rownames = "INDIVIDUALS")
    } else {
      data$genotypes$GT_VCF_NUC <- as.vector(data$genotypes$GT_VCF_NUC)
    }
    data$genotypes$GT_VCF_NUC[is.na(data$genotypes$GT_VCF_NUC)] <- "./."
  } else {
    data$genotypes$GT_VCF_NUC <- NULL
  }

  # gt.vcf
  if (gt.vcf) {
    data$genotypes$GT_VCF <- SeqVarTools::getGenotype(gdsobj = data$vcf.connection, use.names = TRUE)
    if (wide) {
      data$genotypes$GT_VCF %<>% magrittr::set_colnames(x = ., value = data$markers.meta$MARKERS) %>%
        tibble::as_tibble(x = ., rownames = "INDIVIDUALS")
    } else {
      data$genotypes$GT_VCF <- as.vector(data$genotypes$GT_VCF)
    }

    data$genotypes$GT_VCF[is.na(data$genotypes$GT_VCF)] <- "./."
    # replace(data, which(data == what), NA)
  }

  # gt.bin (plink/alt dosage format)
  if (gt.bin) {
    data$genotypes$GT_BIN <- SeqArray::seqGetData(
      gdsfile = data$vcf.connection, var.name = "$dosage_alt")
    if (wide) {
      data$genotypes$GT_BIN %<>% tibble::as_tibble(x = .) %>%
        magrittr::set_colnames(x = ., value = data$markers.meta$MARKERS) %>%
        tibble::add_column(.data = .,
                           "INDIVIDUALS" = data$individuals$INDIVIDUALS,
                           .before = 1)
    } else {
      data$genotypes$GT_BIN <- as.vector(data$genotypes$GT_BIN)
    }
  }

  # genotypes metadata ---------------------------------------------------------
  if (is.logical(vcf.metadata)) {
    overwrite.metadata <- NULL
  } else {
    overwrite.metadata <- vcf.metadata
    vcf.metadata <- TRUE
  }

  if (vcf.metadata) {
    if (verbose) message("Keeping vcf genotypes metadata: yes")
    # detect FORMAT fields available
    have <-  SeqArray::seqSummary(
      gdsfile = data$vcf.connection,
      varname = "annotation/format",
      check = "none", verbose = FALSE)$ID
    # current version doesn't deal well with PL with 3 fields separated with ","
    want <- c("DP", "AD", "GL", "PL", "HQ", "GQ", "GOF", "NR", "NV")
    if (!is.null(overwrite.metadata)) want <- overwrite.metadata
    parse.format.list <- purrr::keep(.x = have, .p = have %in% want)

    # work on parallelization of this part
    data$genotypes.metadata <- purrr::map(
      .x = parse.format.list, .f = parse_gds_metadata, data = data,
      verbose = verbose, parallel.core = parallel.core) %>%
      purrr::flatten(.) %>%
      purrr::flatten_df(.)
  }

  # Remove or not gds connection and file --------------------------------------
  if (!keep.gds) {
    data$vcf.connection <- NULL
    if (file.exists(data$filename)) file.remove(data$filename)
  }

  # Haplotypes or biallelic VCF-------------------------------------------------
  # # recoding genotype
  # if (biallelic) {# biallelic VCF
  #   if (verbose) message("Recoding bi-allelic VCF...")
  # } else {#multi-allelic vcf
  #   if (verbose) message("Recoding VCF haplotype...")
  # }
  # input <- dplyr::rename(input, GT_VCF_NUC = GT)

  # Generating the tidy --------------------------------------------------------

  ## Note to myself: check timig with 1M SNPs to see
  ## if this is more efficient than data.table melt...

  want <- c("MARKERS", "CHROM", "LOCUS", "POS", "COL", "REF", "ALT")
  data$genotypes <- suppressWarnings(
    dplyr::select(data$markers.meta, dplyr::one_of(want))) %>%
    dplyr::bind_cols(
      tibble::as_tibble(
        matrix(
          data = NA,
          nrow = data$n.markers, ncol = data$n.individuals)) %>%
        magrittr::set_colnames(x = ., value = data$individuals$INDIVIDUALS)) %>%
    data.table::as.data.table(.) %>%
    data.table::melt.data.table(
      data = .,
      id.vars = want,
      variable.name = "INDIVIDUALS",
      value.name = "GT",
      variable.factor = FALSE) %>%
    tibble::as_data_frame(.) %>%
    dplyr::select(-GT) %>%
    dplyr::mutate(
      MARKERS = factor(x = MARKERS,
                       levels = data$markers.meta$MARKERS, ordered = TRUE),
      INDIVIDUALS = factor(x = INDIVIDUALS,
                           levels = data$individuals$INDIVIDUALS,
                           ordered = TRUE)) %>%
    dplyr::arrange(MARKERS, INDIVIDUALS) %>%
    dplyr::bind_cols(data$genotypes)

  # data$input <- suppressWarnings(
  #   dplyr::select(data$markers.meta, dplyr::one_of(want)))
  #
  # data$input <- dplyr::bind_rows(replicate(n = length(data$individuals$INDIVIDUALS),
  #                                     expr = data$input, simplify = FALSE)) %>%
  #   dplyr::bind_cols(data$genotypes)

  if (vcf.metadata) {
    data$genotypes %<>% dplyr::bind_cols(data$genotypes.metadata)
    data$genotypes.metadata <- NULL
  }

  # re-calibration of ref/alt alleles ------------------------------------------

  if (ref.calibration) {
    if (verbose) message("\nCalculating REF/ALT alleles...")
    data$genotypes <- radiator::change_alleles(
      data = data$genotypes,
      biallelic = biallelic,
      parallel.core = parallel.core,
      verbose = verbose,
      gt.vcf.nuc = gt.vcf.nuc.bk,
      gt = gt.bk,
      gt.vcf = gt.vcf.bk,
      gt.bin = gt.bin.bk
    )$input
  }

  # include strata
  data$genotypes <- suppressWarnings(
    dplyr::left_join(
      data$genotypes,
      dplyr::select(data$individuals, INDIVIDUALS, POP_ID = STRATA),
      by = "INDIVIDUALS")
  )

  # Re ordering columns
  want <- c("MARKERS", "CHROM", "LOCUS", "POS", "ID", "COL", "INDIVIDUALS",
            "STRATA", "POP_ID",
            "REF", "ALT", "GT_VCF", "GT_VCF_NUC", "GT", "GT_BIN",
            "POLYMORPHIC")

  data$genotypes <- suppressWarnings(
    dplyr::select(data$genotypes, dplyr::one_of(want), dplyr::everything()))

  # Sort id
  data$genotypes <- dplyr::arrange(data$genotypes, POP_ID, INDIVIDUALS)

  timing <- proc.time() - timing
  if (verbose) message("\nTidying vcf time: ", round(timing[[3]]), " sec")
  options(width = opt.change)
  return(data)
}#End tidy_vcf


# Internal nested Function -----------------------------------------------------
#' @title parse_gds_metadata
#' @description function to parse the format field and tidy the results of VCF
#' @rdname parse_gds_metadata
#' @keywords internal
#' @export
parse_gds_metadata <- function(
  x, data = NULL, verbose = TRUE, parallel.core = parallel::detectCores() - 1) {

  res <- list()
  format.name <- x
  # format.name <- x <- "DP"
  # format.name <- x <- "AD"
  # format.name <- x <- "GL"
  # format.name <- x <- "PL"
  # format.name <- x <- "HQ"
  # format.name <- x <- "GQ"
  # format.name <- x <- "GOF"
  # format.name <- x <- "NR"
  # format.name <- x <- "NV"

  if (verbose) message("\nParsing and tidying: ", format.name)

  if (format.name == "AD") {
    if (verbose) message("AD column: splitting into ALLELE_REF_DEPTH and ALLELE_ALT_DEPTH")
    res$AD <- SeqArray::seqGetData(gdsfile = data$vcf.connection,
                                   var.name = "annotation/format/AD")$data %>%
      tibble::as_tibble(.)
    column.vec <- seq_along(res$AD)
    res$AD <- tibble::tibble(ALLELE_REF_DEPTH = res$AD[, column.vec %% 2 == 1] %>%
                               as.matrix(.) %>%
                               as.vector(.),
                             ALLELE_ALT_DEPTH = res$AD[, column.vec %% 2 == 0] %>%
                               as.matrix(.) %>%
                               as.vector(.))
    split.vec <- split_vec_row(x = res$AD, 3, parallel.core = parallel.core)
    res$AD$ALLELE_REF_DEPTH <- clean_ad(x = res$AD$ALLELE_REF_DEPTH, split.vec = split.vec,
                                        parallel.core = parallel.core)

    res$AD$ALLELE_ALT_DEPTH <- clean_ad(x = res$AD$ALLELE_ALT_DEPTH, split.vec = split.vec,
                                        parallel.core = parallel.core)
    split.vec <- NULL
    # test1 <- res$AD
  } # End AD


  # Read depth
  if (format.name == "DP") {
    if (verbose) message("DP column: cleaning and renaming to READ_DEPTH")
    res$DP <- tibble::tibble(READ_DEPTH = SeqArray::seqGetData(
      gdsfile = data$vcf.connection,
      var.name = "annotation/format/DP")$data %>% as.vector(.))
    # test <- res$DP
  } # End DP



  # Haplotype Quality
  # Cleaning HQ: Haplotype quality as phred score
  if (format.name == "HQ") {
    res$HQ <- tibble::tibble(HQ = SeqArray::seqGetData(
      gdsfile = data$vcf.connection,
      var.name = "annotation/format/HQ")$data %>% as.vector(.))
    # test <- res$HQ

    # check HQ and new stacks version with no HQ
    all.missing <- nrow(res$HQ)
    if (all.missing != 0) {
      if (verbose) message("HQ column: Haplotype Quality")
    } else {
      message("HQ values are all missing: removing column")
      res$HQ <- NULL
    }
  } # End HQ

  # Genotypes Quality
  # Cleaning GQ: Genotype quality as phred score
  if (format.name == "GQ") {
    if (verbose) message("GQ column: Genotype Quality")
    res$GQ <- tibble::tibble(GQ = SeqArray::seqGetData(
      gdsfile = data$vcf.connection,
      var.name = "annotation/format/GQ")$data %>% as.vector(.))
    # test <- res$GQ
  } # End GQ

  # GL cleaning
  if (format.name == "GL") {
    if (verbose) message("GL column: cleaning Genotype Likelihood column")
    res$GL <- SeqArray::seqGetData(gdsfile = data$vcf.connection,
                                   var.name = "annotation/format/GL")$data %>%
      tibble::as_tibble(.)
    column.vec <- seq_along(res$GL)
    res$GL <- tibble::tibble(GL_HOM_REF = res$GL[, column.vec %% 3 == 1] %>%
                               as.matrix(.) %>%
                               as.vector(.),
                             GL_HET = res$GL[, column.vec %% 3 == 2] %>%
                               as.matrix(.) %>%
                               as.vector(.),
                             GL_HOM_ALT = res$GL[, column.vec %% 3 == 0] %>%
                               as.matrix(.) %>%
                               as.vector(.))
    res$GL[res$GL == "NaN"] <- NA
  } # End GL

  # Cleaning GOF: Goodness of fit value
  if (format.name == "GOF") {
    if (verbose) message("GOF column: Goodness of fit value")
    res$GOF <- tibble::tibble(GOF = SeqArray::seqGetData(
      gdsfile = data$vcf.connection,
      var.name = "annotation/format/GOF")$data %>% as.vector(.))
    # test <- res$GOF
  } # End GOF

  # Cleaning NR: Number of reads covering variant location in this sample
  if (format.name == "NR") {
    if (verbose) message("NR column: splitting column into the number of variant")
    res$NR <- tibble::tibble(NR = SeqArray::seqGetData(
      gdsfile = data$vcf.connection,
      var.name = "annotation/format/NR")$data %>% as.vector(.))
    # test <- res$NR
  }#End cleaning NR column

  #
  #     input <- clean_nr(x = input,
  #                       split.vec = split.vec,
  #                       parallel.core = parallel.core)

  # NV
  # Cleaning NV: Number of reads containing variant in this sample
  if (format.name == "NV") {
    if (verbose) message("NV column: splitting column into the number of variant")
    res$NR <- tibble::tibble(NV = SeqArray::seqGetData(
      gdsfile = data$vcf.connection,
      var.name = "annotation/format/NV")$data %>% as.vector(.))
    # test <- res$NV
  }#End cleaning NV column


  # input <- clean_nv(x = input,
  #                   split.vec = split.vec,
  #                   parallel.core = parallel.core)




  # test <- res$AD %>%
  #   dplyr::bind_cols(READ_DEPTH = SeqArray::seqGetData(
  #   gdsfile = data$vcf.connection,
  #   var.name = "annotation/format/DP")$data %>% as.vector(.))

  # if (gather.data) {
  #   x <- dplyr::mutate(x, ID = seq(1, n()))
  #   x <- data.table::melt.data.table(
  #     data = data.table::as.data.table(x),
  #     id.vars = "ID",
  #     variable.name = "INDIVIDUALS",
  #     value.name = format.name,
  #     variable.factor = FALSE) %>%
  #     tibble::as_data_frame(.) %>%
  #     dplyr::select(-ID, -INDIVIDUALS)
  # }
  return(res)
}#End parse_gds_metadata

#' @title split_vcf_id
#' @description split VCF ID in parallel
#' @rdname split_vcf_id
#' @keywords internal
#' @export

split_vcf_id <- function(x) {
  res <- dplyr::rename(x, ID = LOCUS) %>%
    tidyr::separate(data = ., col = ID, into = c("LOCUS", "COL"),
                    sep = "_", extra = "drop", remove = FALSE) %>%
    dplyr::mutate_at(.tbl = ., .vars = c("CHROM", "POS", "LOCUS"), .funs = as.character) %>%
    dplyr::mutate_at(.tbl = ., .vars = c("CHROM", "POS", "LOCUS"), .funs = clean_markers_names) %>%
    tidyr::unite(MARKERS, c(CHROM, LOCUS, POS), sep = "__", remove = FALSE)
  return(res)
}#End split_vcf_id



#' @title clean_ad
#' @description Clean allele depth field in VCF.
#' AD column: splitting coverage info into ALLELE_REF_DEPTH and ALLELE_ALT_DEPTH
#' @rdname clean_ad
#' @keywords internal
#' @export
clean_ad <- function(x, split.vec, parallel.core = parallel::detectCores() - 1) {
  clean <- function(x) {
    x <- as.integer(replace_by_na(data = x, what = 0))
  }

  x <- split(x = x, f = split.vec) %>%
    .radiator_parallel(
      X = ., FUN = clean, mc.cores = parallel.core) %>%
    purrr::flatten_int(.)
  return(x)
}#End clean_ad


#' @title clean_pl
#' @description Clean PL column.
#' PL column (normalized, phred-scaled likelihoods for genotypes):
#' separating into PROB_HOM_REF, PROB_HET and PROB_HOM_ALT")
#' Value 1: probability that the site is homozgyous REF
#' Value 2: probability that the sample is heterzygous
#' Value 2: probability that it is homozygous ALT
#' @rdname clean_pl
#' @keywords internal
#' @export
clean_pl <- function(x, split.vec, parallel.core = parallel::detectCores() - 1) {
  clean <- function(x) {
    res <- x %>%
      dplyr::mutate(
        PL = dplyr::if_else(GT_VCF == "./.", NA_character_, PL)) %>%
      tidyr::separate(
        data = ., PL, c("PROB_HOM_REF", "PROB_HET", "PROB_HOM_ALT"),
        sep = ",", extra = "drop", remove = FALSE) %>%
      dplyr::mutate_at(
        .tbl = ., .vars = c("PROB_HOM_REF", "PROB_HET", "PROB_HOM_ALT"),
        .funs = as.numeric) %>%
      dplyr::select(-GT_VCF)
    return(res)
  }#End clean
  x <- dplyr::bind_cols(
    dplyr::select(x, -PL),
    dplyr::ungroup(x) %>%
      dplyr::select(GT_VCF, PL) %>%
      split(x = ., f = split.vec) %>%
      .radiator_parallel(
        X = ., FUN = clean, mc.cores = parallel.core) %>%
      dplyr::bind_rows(.))
  return(x)
}#End clean_pl

#' @title clean_gl
#' @description Clean GL column.
#' GL column: cleaning Genotype Likelihood column
#' @rdname clean_gl
#' @keywords internal
#' @export

clean_gl <- function(x, split.vec, parallel.core = parallel::detectCores() - 1) {
  x <- x %>%
    dplyr::mutate(
      GL = dplyr::if_else(GT_VCF == "./.", NA_character_, GL),
      GL = suppressWarnings(
        stringi::stri_replace_all_fixed(
          GL, c(".,.,.", ".,", ",."), c("NA", "", ""), vectorize_all = FALSE))
    )

  # check GL and new stacks version with no GL
  all.missing <- all(is.na(x$GL))

  if (!all.missing) {
    gl.clean <- max(
      unique(stringi::stri_count_fixed(
        str = unique(sample(x = x$GL, size = 100, replace = FALSE)),
        pattern = ",")
      ), na.rm = TRUE
    )

    if (gl.clean == 2) {
      message("GL column: separating into PROB_HOM_REF, PROB_HET and PROB_HOM_ALT")
      # Value 1: probability that the site is homozgyous REF
      # Value 2: probability that the sample is heterzygous
      # Value 2: probability that it is homozygous ALT
      # system.time(input2 <- input %>%
      #   tidyr::separate(data = ., GL, c("PROB_HOM_REF", "PROB_HET", "PROB_HOM_ALT"), sep = ",", extra = "drop", remove = FALSE) %>%
      #   dplyr::mutate_at(.tbl = ., .vars = c("PROB_HOM_REF", "PROB_HET", "PROB_HOM_ALT"), .funs = as.numeric)
      # )
      clean <- function(x) {
        res <- x %>%
          tidyr::separate(
            data = ., GL, c("PROB_HOM_REF", "PROB_HET", "PROB_HOM_ALT"),
            sep = ",", extra = "drop", remove = FALSE) %>%
          dplyr::mutate_at(
            .tbl = ., .vars = c("PROB_HOM_REF", "PROB_HET", "PROB_HOM_ALT"),
            .funs = as.numeric)
        return(res)
      }

      x <- dplyr::bind_cols(
        dplyr::select(x, -GL),
        dplyr::ungroup(x) %>%
          dplyr::select(GL) %>%
          split(x = ., f = split.vec) %>%
          .radiator_parallel(
            X = ., FUN = clean, mc.cores = parallel.core) %>%
          dplyr::bind_rows(.))

    } else {
      x$GL <- suppressWarnings(as.numeric(x$GL))
    }
  } else {
    message("    GL values are all missing: removing column")
    x <- dplyr::select(x, -GL)
  }
  return(x)
}#End clean_gl

#' @title clean_nr
#' @description Cleaning NR: Number of reads covering variant location in
#' this sample
#' @rdname clean_nr
#' @keywords internal
#' @export
clean_nr <- function(x, split.vec, parallel.core = parallel::detectCores() - 1){
  nr.col <- max(unique(stringi::stri_count_fixed(str = unique(x$NR), pattern = ","))) + 1
  nr.col.names <- stringi::stri_join(rep("NR_", nr.col), seq(1:nr.col))
  nr.col <- NULL

  clean <- function(x, nr.col.names = NULL) {
    res <- tidyr::separate(data = x, col = NR, into = nr.col.names,
                           sep = ",", extra = "drop", remove = FALSE)
    return(res)
  }#End clean

  x <- dplyr::bind_cols(
    dplyr::select(x, -NR),
    dplyr::ungroup(x) %>%
      dplyr::mutate(NR = dplyr::if_else(GT_VCF == "./.", NA_character_, NR)) %>%
      dplyr::select(NR) %>%
      split(x = ., f = split.vec) %>%
      .radiator_parallel_mc(
        X = ., FUN = clean, mc.cores = parallel.core,
        nr.col.names = nr.col.names) %>%
      dplyr::bind_rows(.))
  return(x)
}#End clean_nr

#' @title clean_nv
#' @description Cleaning NV: Number of reads containing variant in this sample
#' @rdname clean_nv
#' @keywords internal
#' @export
clean_nv <- function(x, split.vec, parallel.core = parallel::detectCores() - 1) {

  nv.col <- max(unique(stringi::stri_count_fixed(str = unique(x$NV), pattern = ","))) + 1
  nv.col.names <- stringi::stri_join(rep("NV_", nv.col), seq(1:nv.col))
  nv.col <- NULL

  clean <- function(x, nv.col.names = NULL) {
    res <- tidyr::separate(
      data = x, col = NV, into = nv.col.names,
      sep = ",", extra = "drop", remove = FALSE)
    return(res)
  }

  x <- dplyr::bind_cols(
    dplyr::select(x, -NV),
    dplyr::ungroup(x) %>%
      dplyr::mutate(NV = dplyr::if_else(GT_VCF == "./.", NA_character_, NV)) %>%
      dplyr::select(NV) %>%
      split(x = ., f = split.vec) %>%
      .radiator_parallel_mc(
        X = ., FUN = clean, mc.cores = parallel.core,
        nv.col.names = nv.col.names) %>%
      dplyr::bind_rows(.))
  return(x)
}#End clean_nv

#' @title strata_vcf
#' @description Manage strata
#' @rdname strata_vcf
#' @keywords internal
#' @export
strata_vcf <- function(strata, input, blacklist.id) {

  if (is.null(strata)) {
    message("No strata file provided")
    message("    generating a strata with 1 grouping")
    strata.df <- dplyr::distinct(input, INDIVIDUALS) %>%
      dplyr::mutate(STRATA = rep("pop1", n()))
  } else {
    if (is.vector(strata)) {
      strata.df <- suppressMessages(readr::read_tsv(
        file = strata, col_names = TRUE,
        col_types = readr::cols(.default = readr::col_character())
      ))
    } else {
      strata.df <- strata
      strata.df <- dplyr::mutate_all(.tbl = strata.df, .funs = as.character)
    }
  }

  colnames(strata.df) <- stringi::stri_replace_all_fixed(
    str = colnames(strata.df),
    pattern = "STRATA",
    replacement = "POP_ID",
    vectorize_all = FALSE
  )

  # Remove potential whitespace in pop_id
  strata.df$POP_ID <- clean_pop_names(strata.df$POP_ID)
  colnames.strata <- colnames(strata.df)

  # clean ids
  strata.df$INDIVIDUALS <- clean_ind_names(strata.df$INDIVIDUALS)

  strata.df <- dplyr::distinct(strata.df, POP_ID, INDIVIDUALS, .keep_all = TRUE)

  if (!is.null(strata)) {
    id.vcf <- dplyr::distinct(input, INDIVIDUALS) %>%
      dplyr::mutate(INDIVIDUALS = clean_ind_names(INDIVIDUALS)) %>%
      purrr::flatten_chr(.)

    strata.df <- dplyr::filter(strata.df, INDIVIDUALS %in% id.vcf)
  }


  # filtering the strata if blacklist id available
  if (!is.null(blacklist.id)) {
    strata.df <- dplyr::anti_join(x = strata.df, y = blacklist.id, by = "INDIVIDUALS")
  }

  return(strata.df)
}#End strata_vcf
