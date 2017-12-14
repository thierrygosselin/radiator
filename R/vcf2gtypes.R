# Write a gtypes object from a VCF file

#' @name vcf2gtypes
#' @title VCF to \code{gtypes} with filters and data imputation

#' @description For full details of the function, please use
#' \pkg{radiator} \code{\link[radiator]{genomic_converter}}. This function is a shorcut
#' to output only gtypes object.
#' @inheritParams genomic_converter
#' @inheritParams tidy_genomic_data
#' @inheritParams write_genepop
#' @inheritParams write_genind
#' @inheritParams write_genlight
#' @inheritParams write_structure
#' @inheritParams write_plink
#' @inheritParams write_vcf
#' @inheritParams write_gtypes
#' @inheritParams write_hierfstat
#' @inheritParams radiator_imputations_module

#' @export
#' @rdname vcf2gtypes
#' @import dplyr

#' @seealso \code{\link[radiator]{genomic_converter}}

#' @references Eric Archer, Paula Adams and Brita Schneiders (2016).
#' strataG: Summaries and Population Structure Analyses of
#' Genetic Data. R package version 1.0.5. https://CRAN.R-project.org/package=strataG

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

vcf2gtypes <- function(
  data,
  output,
  filename = NULL,
  blacklist.id = NULL,
  blacklist.genotype = NULL,
  whitelist.markers = NULL,
  monomorphic.out = TRUE,
  snp.ld = NULL,
  common.markers = TRUE,
  maf.thresholds = NULL,
  max.marker = NULL,
  strata = NULL,
  pop.levels = NULL,
  pop.labels = NULL,
  pop.select = NULL,
  imputation.method = NULL,
  hierarchical.levels = "populations",
  num.tree = 50,
  pred.mean.matching = 0,
  random.seed = NULL,
  verbose = FALSE,
  parallel.core = detectCores() - 1
) {

  res <- genomic_converter(
    data,
    output = "gtypes",
    filename = filename,
    blacklist.id = blacklist.id,
    blacklist.genotype = blacklist.genotype,
    whitelist.markers = whitelist.markers,
    monomorphic.out = monomorphic.out,
    snp.ld = snp.ld,
    common.markers = common.markers,
    maf.thresholds = maf.thresholds,
    max.marker = max.marker,
    strata = strata,
    pop.levels = pop.levels,
    pop.labels = pop.labels,
    pop.select = pop.select,
    imputation.method = imputation.method,
    hierarchical.levels = hierarchical.levels,
    num.tree = num.tree,
    pred.mean.matching = pred.mean.matching,
    random.seed = random.seed,
    verbose = verbose,
    parallel.core = parallel.core
  )
  return(res)
}
