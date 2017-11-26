# Convert genomic dataset to other useful genomic format with filter and imputation

#' @name genomic_converter

#' @title Conversion tool among several genomic formats

#' @description The arguments in the \code{genomic_converter} function were tailored for the
#' reality of GBS/RADseq data while maintaining a reproducible workflow.
#'
#' \itemize{
#'   \item \strong{Input file:} 12 file formats are supported (see \code{data} argument below)
#'   \item \strong{Filters:} genotypes, markers, individuals and populations can be
#'   filtered and/or selected in several ways using blacklist,
#'   whitelist and other arguments
#'   \item \strong{Imputations:} Map-independent imputation of missing genotype/alleles
#'   using Random Forest or the most frequent category.
#'   \item \strong{Parallel:} Some parts of the function are designed to be conduncted on multiple CPUs
#'   \item \strong{Output:} 18 output file formats are supported (see \code{output} argument below)
#' }

#' @param output 19 genomic data formats can be exported: tidy (by default),
#' genind, genlight, vcf (for file format version, see details below), plink, genepop,
#' structure, arlequin, hierfstat, gtypes (strataG), bayescan, betadiv, pcadapt,
#' hzar, snprelate (for SNPRelate's GDS, see details).
#' Use a character string,
#' e.g. \code{output = c("genind", "genepop", "structure")}, to have preferred
#' output formats generated. The tidy format is generated automatically.
#' Default: \code{output = NULL}.

#' @param filename (optional) The filename prefix for the object in the global environment
#' or the working directory. Default: \code{filename = NULL}. A default name will be used,
#' customized with the output file(s) selected.

#' @inheritParams tidy_genomic_data
#' @inheritParams write_genepop
#' @inheritParams write_genind
#' @inheritParams write_genlight
#' @inheritParams write_structure
#' @inheritParams write_arlequin
#' @inheritParams write_plink
#' @inheritParams write_vcf
#' @inheritParams write_gtypes
#' @inheritParams write_hierfstat
#' @inheritParams write_bayescan
#' @inheritParams write_pcadapt
#' @inheritParams write_hzar
#' @inheritParams radiator_imputations_module
# @inheritParams write_snprelate


#' @details
#' \strong{Input files:} Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.
#'
#'
#' \strong{Imputations details:}
#'
#' The imputations using Random Forest requires more time to compute and can take several
#' minutes and hours depending on the size of the dataset and polymorphism of
#' the species used. e.g. with a low polymorphic taxa, and a data set
#' containing 30\% missing data, 5 000 haplotypes loci and 500 individuals
#' will require 15 min.
#'
#'
#' \strong{VCF file format version:}
#'
#' If you need a different file format version than the current one, just change
#' the version inside the newly created VCF, that should do the trick.
#' \href{https://vcftools.github.io/specs.html}{For more
#' information on Variant Call Format specifications}.
#'
#'
#' \strong{SNPRelate:}
#'
#' I received too many emails with students having problem with SNPRelate, so
#' I'm removing the option to output in SNPRelate from genomic_coverter.
#' If you really want to do it,
#' \href{https://github.com/zhengxwen/SNPRelate}{install SNPRelate}
#' separately following the author's instruction. From a tidy data set created
#' in genomic_converter you can download my code
#' \href{https://www.dropbox.com/s/7xujizkvpi0ddac/write_snprelate.R?dl=0}{here}
#' and source it to get a SNPRelate object. However, I'm no longer offering
#' support for this.


#' @return The function returns an object (list). The content of the object
#' can be listed with \code{names(object)} and use \code{$} to isolate specific
#' object (see examples). Some output format will write the output file in the
#' working directory. The tidy genomic data frame is generated automatically.

#' @export
#' @rdname genomic_converter
# @importFrom adegenet df2genind
#' @importFrom dplyr n_distinct summarise group_by ungroup mutate select tally distinct summarise
#' @importFrom stringi stri_join stri_replace_all_fixed stri_sub
#' @importFrom purrr flatten_chr
#' @importFrom tidyr gather
#' @importFrom parallel detectCores
#' @importFrom utils installed.packages

#' @examples
#' \dontrun{
#' # The simplest form of the function:
#' snowcrab <- genomic_converter(
#'     data = "batch_1.vcf",
#'     output = c("genlight", "genepop"),
#'     strata = "snowcrab.strata.tsv"
#'     )
#' # With imputations using random forest:
#' snowcrab <- genomic_converter(
#'     data = "batch_1.vcf",
#'     output = c("genlight", "genepop"),
#'     strata = "snowcrab.strata.tsv",
#'     imputation.method = "rf"
#'     )
#'
#' #Get the content of the object created using:
#' names(snowcrab)
#' #To isolate the genlight object (without imputation):
#' genlight.no.imputation <- snowcrab$genlight.no.imputation
#' }

#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011)
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences.
#' G3, 1, 171-182.

#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013)
#' Stacks: an analysis tool set for population genomics.
#' Molecular Ecology, 22, 3124-3140.

#' @references Jombart T (2008) adegenet: a R package for the multivariate
#' analysis of genetic markers. Bioinformatics, 24, 1403-1405.

#' @references Jombart T, Ahmed I (2011) adegenet 1.3-1:
#' new tools for the analysis of genome-wide SNP data.
#' Bioinformatics, 27, 3070-3071.

#' @references Ishwaran H. and Kogalur U.B. (2015). Random Forests for Survival,
#'  Regression and Classification (RF-SRC), R package version 1.6.1.

#' @references Ishwaran H. and Kogalur U.B. (2007). Random survival forests
#' for R. R News 7(2), 25-31.

#' @references Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S. (2008).
#' Random survival forests. Ann. Appl. Statist. 2(3), 841-860.

#' @references Lamy T, Legendre P, Chancerelle Y, Siu G, Claudet J (2015)
#' Understanding the Spatio-Temporal Response of Coral Reef Fish Communities to
#' Natural Disturbances: Insights from Beta-Diversity Decomposition.
#' PLoS ONE, 10, e0138696.

#' @references Danecek P, Auton A, Abecasis G et al. (2011)
#' The variant call format and VCFtools.
#' Bioinformatics, 27, 2156-2158.

#' @references Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR,
#' Bender D, et al.
#' PLINK: a tool set for whole-genome association and population-based linkage
#' analyses.
#' American Journal of Human Genetics. 2007; 81: 559–575. doi:10.1086/519795

#' @references Goudet, J. (1995) FSTAT (Version 1.2): A computer program to
#' calculate F- statistics. Journal of Heredity, 86, 485-486.
#' @references Goudet, J. (2005) hierfstat, a package for r to compute and test hierarchical F-statistics. Molecular Ecology Notes, 5, 184-186.

#' @references Eric Archer, Paula Adams and Brita Schneiders (2016).
#' strataG: Summaries and Population Structure Analyses of
#' Genetic Data. R package version 1.0.5. https://CRAN.R-project.org/package=strataG

#' @references Zheng X, Levine D, Shen J, Gogarten SM, Laurie C, Weir BS.
#' A high-performance computing toolset for relatedness and principal component
#' analysis of SNP data. Bioinformatics. 2012;28: 3326-3328.
#' doi:10.1093/bioinformatics/bts606

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


#' @seealso \code{beta.div} is available on Pierre Legendre web site \url{http://adn.biol.umontreal.ca/~numericalecology/Rcode/}
#' \code{randomForestSRC} is available on CRAN \url{http://cran.r-project.org/web/packages/randomForestSRC/}
#' and github \url{https://github.com/ehrlinger/randomForestSRC}

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com} and
#' Laura Benestan \email{laura.benestan@@icloud.com} (for betadiv)

genomic_converter <- function(
  data,
  output = NULL,
  vcf.metadata = FALSE,
  filename = NULL,
  blacklist.id = NULL,
  blacklist.genotype = NULL,
  whitelist.markers = NULL,
  monomorphic.out = TRUE,
  snp.ld = NULL,
  common.markers = TRUE,
  maf.thresholds = NULL,
  maf.pop.num.threshold = 1,
  maf.approach = "SNP",
  maf.operator = "OR",
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
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE
) {
  if (verbose) {
    cat("#######################################################################\n")
    cat("##################### radiator::genomic_converter #####################\n")
    cat("#######################################################################\n")
    # Check that strataG is installed
    if ("gtypes" %in% output) {
      if (!"strataG" %in% utils::installed.packages()[,"Package"]) {
        stop("Please install strataG for this output option:\n
devtools::install_github('ericarcher/strataG', build_vignettes = TRUE)")
      }
    }

    # Check that snprelate is installed
    # if ("snprelate" %in% output) {
    #   if (!"SNPRelate" %in% utils::installed.packages()[,"Package"]) {
    #     stop("Please install SNPRelate:\n
    #        github::zhengxwen/SNPRelate")
    #   }
    # }

    # Checking for missing and/or default arguments-------------------------------
    if (missing(data)) stop("Input file missing")
    if (!is.null(pop.levels) & is.null(pop.labels)) pop.labels <- pop.levels
    if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")

    message("Function arguments and values:")
    message("Working directory: ", getwd())

    if (is.vector(data)) {
      message("Input file: ", data)
    } else {
      message("Input file: from global environment")
    }

    if (is.null(strata)) {
      message("Strata: no")
    } else {
      message("Strata: ", strata)
    }

    if (is.null(pop.levels)) {
      message("Population levels: no")
    } else {
      message("Population levels: ", stringi::stri_join(pop.levels, collapse = ", "))
    }

    if (is.null(pop.levels)) {
      message("Population labels: no")
    } else {
      message(stringi::stri_join("Population labels: ", stringi::stri_join(pop.labels, collapse = ", ")))
    }

    if (is.null(output)) {
      message("Output format(s): tidy")
    } else {
      message(stringi::stri_join("Output format(s): tidy, ", stringi::stri_join(output, collapse = ", ")))
    }

    if (is.null(filename)) {
      message("Filename prefix: no")
    } else {
      message("Filename prefix: ", filename, "\n")
    }


    message("Filters: ")
    if (is.null(blacklist.id)) {
      message("Blacklist of individuals: no")
    } else {
      message("Blacklist of individuals: ", blacklist.id)
    }

    if (is.null(blacklist.genotype)) {
      message("Blacklist of genotypes: no")
    } else {
      message("Blacklist of genotypes: ", blacklist.genotype)
    }

    if (is.null(whitelist.markers)) {
      message("Whitelist of markers: no")
    } else {
      if (is.vector(whitelist.markers)) {
        message("Whitelist of markers: ", whitelist.markers)
      } else {
        message("Whitelist of markers: ", nrow(whitelist.markers))
      }
    }

    message("monomorphic.out: ", monomorphic.out)
    if (is.null(snp.ld)) {
      message("snp.ld: no")
    } else {
      message("snp.ld: ", snp.ld)
    }
    message("common.markers: ", common.markers)
    if (is.null(max.marker)) {
      message("max.marker: no")
    } else {
      message("max.marker: ", max.marker)
    }

    if (is.null(pop.select)) {
      message("pop.select: no")
    } else {
      message(stringi::stri_join("pop.select: ", stringi::stri_join(pop.select, collapse = ", ")))
    }
    if (is.null(maf.thresholds)) {
      message("maf.thresholds: no")
    } else {
      message(stringi::stri_join("maf.thresholds: ", stringi::stri_join(maf.thresholds, collapse = ", ")))
      message("maf.pop.num.threshold: ", maf.pop.num.threshold)
      message("maf.approach: ", maf.approach)
      message("maf.operator: ", maf.operator)
    }

    message(stringi::stri_join("\n", "Imputations options:"))
    if (is.null(imputation.method)) {
      message("imputation.method: no")
    } else {
      message("imputation.method: ", imputation.method)
      message("hierarchical.levels: ", hierarchical.levels)
    }
    message("\nparallel.core: ", parallel.core, "\n")
    cat("#######################################################################\n")
  }
  timing <- proc.time()
  opt.change <- getOption("width")
  options(width = 70)

  # Filename -------------------------------------------------------------------
  # Get date and time to have unique filenaming
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")

  if (is.null(filename)) {
    filename <- stringi::stri_join("radiator_data_", file.date)

    if (!is.null(imputation.method)) {
      filename.imp <- stringi::stri_join("radiator_data_imputed_", file.date)
    }
  } else {
    filename.problem <- file.exists(filename)
    if (filename.problem) {
      filename <- stringi::stri_join(filename, "_", file.date)
    }
    if (!is.null(imputation.method)) {
      filename.imp <- stringi::stri_join(filename, "_imputed")
    }
  }

  # File type detection --------------------------------------------------------
  data.type <- detect_genomic_format(data = data)

  # Strata argument required for VCF and haplotypes files-----------------------
  if (data.type == "vcf.file" & is.null(strata)) stop("strata argument is required")
  if (data.type == "haplo.file" & is.null(strata)) stop("strata argument is required")

  # Import----------------------------------------------------------------------
  if (verbose) message("\nImporting data\n")
  input <- radiator::tidy_genomic_data(
    data = data,
    vcf.metadata = vcf.metadata,
    blacklist.id = blacklist.id,
    blacklist.genotype = blacklist.genotype,
    whitelist.markers = whitelist.markers,
    monomorphic.out = monomorphic.out,
    max.marker = max.marker,
    snp.ld = snp.ld,
    common.markers = common.markers,
    maf.thresholds = maf.thresholds,
    maf.pop.num.threshold = maf.pop.num.threshold,
    maf.approach = maf.approach,
    maf.operator = maf.operator,
    strata = strata,
    pop.levels = pop.levels,
    pop.labels = pop.labels,
    pop.select = pop.select,
    filename = NULL,
    verbose = FALSE
  )

  if(verbose) message("\nPreparing data for output\n")

  input$GT <- stringi::stri_replace_all_fixed(
    str = input$GT,
    pattern = c("/", ":", "_", "-", "."),
    replacement = c("", "", "", "", ""),
    vectorize_all = FALSE
  )

  # create a strata.df
  # strata.df <- input %>%
  #   distinct(INDIVIDUALS, POP_ID)
  # # strata <- strata.df
  pop.levels <- levels(input$POP_ID)
  pop.labels <- pop.levels


  # prepare output res list
  res <- list()
  res$tidy.data <- input

  # Biallelic detection --------------------------------------------------------
  biallelic <- radiator::detect_biallelic_markers(data = input, verbose = verbose)

  if (!biallelic && "genlight" %in% output || !biallelic && "plink" %in% output) {
    stop("output chosen doesn't work with multi-allelic data")
  }

  # overide genind when marker number > 20K ------------------------------------
  if ("genind" %in% output & biallelic) {
    # detect the number of marker
    marker.number <- dplyr::n_distinct(input$MARKERS)
    if (marker.number > 20000) {

      # When genlight is also selected, remove automatically
      if ("genlight" %in% output) {
        message("Removing the genind output option, the genlight is more suitable with current marker number")
        output <- stringi::stri_replace_all_fixed(
          str = output,
          pattern = "genind",
          replacement = "",
          vectorize_all = FALSE
        )
      } else {
        message("IMPORTANT: you have > 20 000 markers (", marker.number, ")",
                "\nDo you want the more suitable genlight object instead of the current genind? (y/n):")
        overide.genind <- as.character(readLines(n = 1))
        if (overide.genind == "y") {
          output <- stringi::stri_replace_all_fixed(
            str = output,
            pattern = "genind",
            replacement = "genlight",
            vectorize_all = FALSE
          )
        }
      }
    }
  }

  # Imputations-----------------------------------------------------------------
  if (!is.null(imputation.method)) {

    input.imp <- radiator::radiator_imputations_module(
      data = input,
      imputation.method = imputation.method,
      hierarchical.levels = hierarchical.levels,
      num.tree = num.tree,
      pred.mean.matching = pred.mean.matching,
      random.seed = random.seed,
      verbose = verbose,
      parallel.core = parallel.core,
      filename = NULL
    )

    res$tidy.data.imp <- input.imp

  } # End imputations

  # OUTPUT ---------------------------------------------------------------------

  # GENEPOP --------------------------------------------------------------------
  if ("genepop" %in% output) {
    if (verbose) message("Generating genepop file without imputation")
    radiator::write_genepop(
      data = input,
      pop.levels = pop.levels,
      filename = filename
    )

    if (!is.null(imputation.method)) {
      if (verbose) message("Generating genepop file WITH imputations")
      radiator::write_genepop(
        data = input.imp,
        pop.levels = pop.levels,
        filename = filename.imp
      )
    }
  } # end genepop output

  # hierfstat --------------------------------------------------------------------
  if ("hierfstat" %in% output) {
    if (verbose) message("Generating hierfstat file without imputation")
    res$hierfstat.no.imputation <- radiator::write_hierfstat(
      data = input,
      filename = filename
    )

    if (!is.null(imputation.method)) {
      if (verbose) message("Generating hierfstat file WITH imputations")
      res$hierfstat.imputed <- radiator::write_hierfstat(
        data = input.imp,
        filename = filename.imp
      )
    }
  } # end hierfstat output

  # strataG --------------------------------------------------------------------
  if ("gtypes" %in% output) {
    if (verbose) message("Generating strataG gtypes object without imputation")
    res$gtypes.no.imputation <- radiator::write_gtypes(data = input)

    if (!is.null(imputation.method)) {
      if (verbose) message("Generating strataG gtypes object WITH imputations")
      res$gtypes.imputed <- radiator::write_gtypes(data = input.imp)
    }
  } # end strataG output

  # structure --------------------------------------------------------------------
  if ("structure" %in% output) {
    if (verbose) message("Generating structure file without imputation")
    radiator::write_structure(
      data = input,
      pop.levels = pop.levels,
      markers.line = TRUE,
      filename = filename
    )

    if (!is.null(imputation.method)) {
      if (verbose) message("Generating structure file WITH imputations")
      radiator::write_structure(
        data = input.imp,
        pop.levels = pop.levels,
        markers.line = TRUE,
        filename = filename.imp
      )
    }
  } # end structure output

  # betadiv --------------------------------------------------------------------
  if ("betadiv" %in% output) {
    if (!biallelic) stop("betadiv output is currently implemented for biallelic data only")
    if (verbose) message("Generating betadiv object without imputation")
    res$betadiv.no.imputation <- radiator::write_betadiv(data = input)

    if (!is.null(imputation.method)) {
      if (verbose) message("Generating betadiv object WITH imputations")
      res$betadiv.imputed <- radiator::write_betadiv(data = input.imp)
    }
  } # end betadiv output

  # arlequin --------------------------------------------------------------------
  if ("arlequin" %in% output) {
    if (verbose) message("Generating arlequin file without imputation")
    radiator::write_arlequin(
      data = input,
      pop.levels = pop.levels,
      filename = filename
    )

    if (!is.null(imputation.method)) {
      if (verbose) message("Generating arlequin file WITH imputations")
      radiator::write_arlequin(
        data = input.imp,
        pop.levels = pop.levels,
        filename = filename.imp
      )
    }
  } # end arlequin output

  # GENIND ---------------------------------------------------------------------
  if ("genind" %in% output) {
    if (verbose) message("Generating adegenet genind object without imputation")
    res$genind.no.imputation <- radiator::write_genind(data = input)

    if (!is.null(imputation.method)) {
      if (verbose) message("Generating adegenet genind object WITH imputations")
      res$genind.imputed <- radiator::write_genind(data = input.imp)
    }
  } # end genind

  # GENLIGHT -------------------------------------------------------------------
  if ("genlight" %in% output) {
    if (verbose) message("Generating adegenet genlight object without imputation")
    res$genlight.no.imputation <- radiator::write_genlight(data = input,
                                                           biallelic = TRUE)

    if (!is.null(imputation.method)) {
      if (verbose) message("Generating adegenet genlight object WITH imputations")
      res$genlight.imputed <- radiator::write_genlight(data = input.imp,
                                                       biallelic = TRUE)
    }
  } # end genlight output

  # VCF ------------------------------------------------------------------------
  if ("vcf" %in% output) {
    if (!biallelic) stop("vcf output is currently implemented for biallelic data only")
    if (verbose) message("Generating VCF file without imputation")
    radiator::write_vcf(
      data = input,
      filename = filename
    )

    if (!is.null(imputation.method)) {
      if (verbose) message("Generating VCF file WITH imputations")
      radiator::write_vcf(
        data = input.imp,
        filename = filename.imp
      )
    }
  } # end vcf output

  # PLINK ----------------------------------------------------------------------
  if ("plink" %in% output) {
    if (verbose) message("Generating PLINK tped/tfam files without imputation")
    radiator::write_plink(
      data = input,
      filename = filename
    )

    if (!is.null(imputation.method)) {
      if (verbose) message("Generating PLINK tped/tfam files WITH imputations")
      radiator::write_plink(
        data = input.imp,
        filename = filename.imp
      )
    }
  } # end plink output


  # SNPRelate ------------------------------------------------------------------
  # if ("snprelate" %in% output) {
  #   if (verbose) message("Generating SNPRelate object without imputation")
  #   res$snprelate.no.imputation <- radiator::write_snprelate(data = input,
  #                                                            biallelic = TRUE)
  #   if (!is.null(imputation.method)) {
  #     if (verbose) message("Generating SNPRelate object WITH imputations")
  #     res$snprelate.imputed <- radiator::write_snprelate(data = input.imp,
  #                                                        biallelic = TRUE)
  #   }
  # }

  # bayescan -------------------------------------------------------------------
  if ("bayescan" %in% output) {
    if (verbose) message("Generating BayeScan object without imputation")
    res$bayescan.no.imputation <- radiator::write_bayescan(
      data = input,
      pop.select = pop.select,
      snp.ld = snp.ld,
      filename = filename)
    if (!is.null(imputation.method)) {
      if (verbose) message("Generating BayeScan object WITH imputations")
      res$bayescan.imputed <- radiator::write_bayescan(
        data = input.imp,
        pop.select = pop.select,
        snp.ld = snp.ld,
        filename = filename)
    }
  }

  # pcadapt -------------------------------------------------------------------
  if ("pcadapt" %in% output) {
    if (verbose) message("Generating pcadapt file and object without imputation")
    res$pcadapt.no.imputation <- radiator::write_pcadapt(
      data = input,
      pop.select = pop.select,
      snp.ld = snp.ld,
      maf.thresholds = maf.thresholds,
      maf.pop.num.threshold = maf.pop.num.threshold,
      maf.approach = maf.approach,
      maf.operator = maf.operator,
      filename = filename,
      parallel.core = parallel.core
    )

    if (!is.null(imputation.method)) {
      if (verbose) message("Generating pcadapt file and object WITH imputations")
      res$pcadapt.imputed <- radiator::write_pcadapt(
        data = input.imp,
        pop.select = pop.select,
        snp.ld = snp.ld,
        maf.thresholds = maf.thresholds,
        maf.pop.num.threshold = maf.pop.num.threshold,
        maf.approach = maf.approach,
        maf.operator = maf.operator,
        filename = filename,
        parallel.core = parallel.core
      )
    }
  }


  # hzar -------------------------------------------------------------------
  if ("hzar" %in% output) {
    if (verbose) message("Generating HZAR file without imputation")
    res$hzar.no.imputation <- radiator::write_hzar(
      data = input,
      distance = NULL,
      filename = filename,
      parallel.core = parallel.core
    )

    if (!is.null(imputation.method)) {
      if (verbose) message("Generating HZAR file WITH imputations")
      res$hzar.imputed <- radiator::write_hzar(
        data = input.imp,
        distance = NULL,
        filename = filename,
        parallel.core = parallel.core
      )
    }
  }

  # dadi -----------------------------------------------------------------------
  # not yet implemented, use vcf2dadi

  # outout results -------------------------------------------------------------
  n.markers <- dplyr::n_distinct(input$MARKERS)
  if (tibble::has_name(input, "CHROM")) {
    n.chromosome <- dplyr::n_distinct(input$CHROM)
  } else {
    n.chromosome <- "no chromosome info"
  }
  n.individuals <- dplyr::n_distinct(input$INDIVIDUALS)
  n.pop <- dplyr::n_distinct(input$POP_ID)

  if (verbose) {
    cat("############################### RESULTS ###############################\n")
    message("Tidy data in your global environment")
    message("Depending on output selected, check the list in your global environment and your working directory")
    message("Data format of input: ", data.type)
    if (biallelic) {
      message("Biallelic data")
    } else{
      message("Multiallelic data")
    }
    if (common.markers) {
      message("Number of common markers: ", n.markers)
    } else {
      message("Number of markers: ", n.markers)
    }
    message("Number of chromosome/contig/scaffold: ", n.chromosome)
    message("Number of individuals ", n.individuals)
    message("Number of populations ", n.pop)
    timing <- proc.time() - timing
    message("\nComputation time: ", round(timing[[3]]), " sec")
    cat("################ radiator::genomic_converter completed ################\n")
  }
  return(res)
} # end genomic_converter
