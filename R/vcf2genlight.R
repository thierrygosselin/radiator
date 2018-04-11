# Write a adegenet genlight object from VCF file

#' @name vcf2genlight
#' @title VCF to \code{adegenet genlight} object with filters and data imputation

#' @description For full details of the function, please use
#' \pkg{radiator} \code{\link[radiator]{genomic_converter}}. This function is a shorcut
#' to output only genlight object.
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
#' @rdname vcf2genlight

#' @references Jombart T (2008) adegenet: a R package for the multivariate
#' analysis of genetic markers. Bioinformatics, 24, 1403-1405.
#' @references Jombart T, Ahmed I (2011) adegenet 1.3-1:
#' new tools for the analysis of genome-wide SNP data.
#' Bioinformatics, 27, 3070-3071.

#' @seealso \code{adegenet} is available on CRAN \url{http://cran.r-project.org/web/packages/adegenet/} and github \url{https://github.com/thibautjombart/}
#' \code{\link[radiator]{genomic_converter}}

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

vcf2genlight <- function(
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
    output = "genlight",
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

