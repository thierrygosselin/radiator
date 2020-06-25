# Convert genomic dataset to other useful genomic format with filter and imputation

#' @name genomic_converter

#' @title Conversion tool among several genomic formats

#' @description The arguments in the \code{genomic_converter} function were tailored for the
#' reality of GBS/RADseq data while maintaining a reproducible workflow.
#'
#' \itemize{
#'   \item \strong{Input file:} 14 diploid file formats are supported
#'   (see \code{data} argument below).
#'   \item \strong{Filters:} see \emph{Advance mode} section below for ways to
#'   use blacklist and whitelist related arguments.
#'   For best results with unfiltered datasets, use \code{\link{filter_rad}}
#'   (\code{genomic_converter} is included in that function!).
#'   \item \strong{Imputations:} deprecated module no longer available
#'   in \emph{genomic_converter} (see \emph{Life cycle} section below).
#'   \item \strong{Parallel:} Some parts of the function are designed to be conduncted on multiple CPUs
#'   \item \strong{Output:} 29 output file formats are supported (see \code{output} argument below)
#' }

#' @param output 29 genomic data formats can be exported: tidy (by default),
#' genepop, genind, genlight, vcf (for file format version, see details below),
#' plink, structure, faststructure, arlequin, hierfstat, gtypes (strataG),
#' bayescan, betadiv, pcadapt, hzar, fineradstructure, related, seqarray,
#' snprelate, maverick, genepopedit, rubias, hapmap and dadi.
#' Use a character string,
#' e.g. \code{output = c("genind", "genepop", "structure")}, to have preferred
#' output formats generated. With default, only the tidy format is generated.
#' Default: \code{output = NULL}.

#' @param filename (optional) The filename prefix for the object in the global environment
#' or the working directory. Default: \code{filename = NULL}. A default name will be used,
#' customized with the output file(s) selected.

#' @inheritParams tidy_genomic_data
#' @inheritParams radiator_common_arguments
#' @inheritParams read_strata
#' @inheritParams write_genepop
#' @inheritParams write_genind
#' @inheritParams write_genlight
#' @inheritParams write_structure
#' @inheritParams write_faststructure
#' @inheritParams write_arlequin
#' @inheritParams write_plink
#' @inheritParams write_vcf
#' @inheritParams write_gtypes
#' @inheritParams write_hierfstat
#' @inheritParams write_bayescan
#' @inheritParams write_pcadapt
#' @inheritParams write_hzar
#' @inheritParams write_fineradstructure
#' @inheritParams write_related
#' @inheritParams write_snprelate
#' @inheritParams write_stockr
#' @inheritParams write_genepopedit
#' @inheritParams write_rubias
#' @inheritParams write_hapmap


#' @section Input genomic datasets:
#' \enumerate{
#' \item GDS file or object, must end with \code{.gds} or \code{.rad}:
#' documented in \code{\link{read_vcf}}
#'
#' \item VCF files must end with \code{.vcf}: documented in \code{\link{tidy_vcf}}
#'
#' \item PLINK files must end with \code{.tped} or \code{.bed}: documented in \code{\link{tidy_plink}}
#'
#' \item genind object from
#' \href{https://github.com/thibautjombart/adegenet}{adegenet}:
#' documented in \code{\link{tidy_genind}}.
#'
#' \item genlight object from
#' \href{https://github.com/thibautjombart/adegenet}{adegenet}:
#' documented in \code{\link{tidy_genlight}}.
#'
#' \item gtypes object from
#' \href{https://github.com/EricArcher/strataG}{strataG}:
#' documented in \code{\link{tidy_gtypes}}.
#'
#' \item dart data from \href{http://www.diversityarrays.com}{DArT}:
#' documented in \code{\link{read_dart}}.
#'
#' \item genepop file must end with \code{.gen}, documented in \code{\link{tidy_genepop}}.
#'
#' \item fstat file must end with \code{.dat}, documented in \code{\link{tidy_fstat}}.
#'
#' \item haplotype file created in STACKS (e.g. \code{data = "batch_1.haplotypes.tsv"}).
#' To make the haplotype file population ready, you need the \code{strata} argument.
#'
#' \item Data frames: documented in \code{\link{tidy_wide}}.
#' }

#' @section Advance mode:
#'
#' \emph{dots-dots-dots ...} allows to pass several arguments for fine-tuning the function:
#' \enumerate{
#'
#' \item{path.folder: } use this argument to specify an output folder.
#' Default: \code{path.folder = "radiator_genomic_converter"}.
#'
#' \item \code{vcf.metadata} (optional, logical or string).
#' Default: \code{vcf.metadata = TRUE}. Documented in \code{\link{tidy_vcf}}.
#'
#' \item \code{vcf.stats} (optional, logical).
#' Default: \code{vcf.stats = TRUE}.
#' Documented in \code{\link{tidy_vcf}}.
#'
#' \item \code{whitelist.markers} (optional) Default \code{whitelist.markers = NULL}.
#' Documented in \code{\link{read_whitelist}}.
#'
#' \item \code{filter.common.markers} (optional, logical).
#' Default: \code{filter.common.markers = TRUE}.
#' By defaults, only common markers are kept in the dataset.
#' Documented in \code{\link{filter_common_markers}}.
#'
#' \item \code{filter.monomorphic} (logical, optional)
#' Default: \code{filter.monomorphic = TRUE}.
#' By defaults, only polymorphic markers across strata are kept in the dataset.
#' Documented in \code{\link{filter_monomorphic}}.
#'
#' \item \emph{individuals to blacklist ? } Use the strata file for this.
#' Documented in \code{\link{read_strata}}.
#'
#' \item \code{keep.allele.names} argument used when tidying genind object.
#' Documented in \code{\link{tidy_genind}}.
#' Default: \code{keep.allele.names = FALSE}.
#'
#' \item \code{blacklist.id} (optional) Default: \code{blacklist.id = NULL}.
#' Ideally, managed in the strata file.
#' Documented in \code{\link{read_strata}} and \code{\link{read_blacklist_id}}.
#' }

#' @section Life cycle:
#'
#' Map-independent imputation of missing genotype is avaible in my other R
#' package called \href{https://github.com/thierrygosselin/grur}{grur}.
#'
#' Use \href{https://github.com/thierrygosselin/grur}{grur} to :
#' \enumerate{
#' \item \strong{Visualize your missing data: } before imputing your genotypes,
#' visualize your missing data.
#' Several visual tools are available inside \href{https://github.com/thierrygosselin/grur}{grur} to
#' help you decide the best strategy after.
#' \item \strong{Optimize: }
#' use \href{https://github.com/thierrygosselin/grur}{grur} imputation module
#' and other functions to optimize the imputations of your dataset.
#' You need to test arguments. Failing to conduct tests and adjust imputations arguments
#' will \strong{generate artifacts} and/or \strong{exacerbate bias}.
#' Using defaults is not optional here...
#' \item \strong{genomic_converter: }
#' use the output argument inside \href{https://github.com/thierrygosselin/grur}{grur}
#' imputation module to generate the required formats.
#' }


#' @section VCF file format version:
#'
#' If you need a different VCF file format version than the current one, just change
#' the version inside the newly created VCF, that should do the trick.
#' \href{https://vcftools.github.io/specs.html}{For more
#' information on Variant Call Format specifications}.

#' @return The function returns an object (list). The content of the object
#' can be listed with \code{names(object)} and use \code{$} to isolate specific
#' object (see examples). Some output format will write the output file in the
#' working directory. The tidy genomic data frame is generated automatically.

#' @export
#' @rdname genomic_converter
#' @examples
#' \dontrun{
#' #To verify your file is detected by radiator as the correct format:
#' radiator::detect_genomic_format(data = "populations.snps.vcf")
#'
#' # The simplest form of the function:
#' require(SeqVarTools) # when using vcf as input file
#' snowcrab <- genomic_converter(
#'                    data = "populations.snps.vcf", strata = "snowcrab.strata.tsv",
#'                    output = c("genlight", "genepop"))
#'
#' #Get the content of the object created using:
#' names(snowcrab)
#' #To isolate the genlight object (without imputation):
#' genlight <- snowcrab$genlight
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

#' @references Malinsky M, Trucchi E, Lawson D, Falush D (2018)
#' RADpainter and fineRADstructure: population inference from RADseq data.
#' bioRxiv, 057711.

#' @references Pew J, Muir PH, Wang J, Frasier TR (2015)
#' related: an R package for analysing pairwise relatedness from codominant
#' molecular markers.
#' Molecular Ecology Resources, 15, 557-561.

#' @references Raj A, Stephens M, Pritchard JK (2014)
#' fastSTRUCTURE: Variational Inference of Population Structure in Large SNP
#' Datasets. Genetics, 197, 573-589.

#' @references Verity R, Nichols RA (2016) Estimating the Number of
#' Subpopulations (K) in Structured Populations.
#' Genetics, 203, genetics.115.180992-1839.

#' @references Zheng X, Gogarten S, Lawrence M, Stilp A, Conomos M, Weir BS,
#' Laurie C, Levine D (2017). SeqArray -- A storage-efficient high-performance
#' data format for WGS variant calls.
#' Bioinformatics.

#' @seealso \code{beta.div} is available on Pierre Legendre web site \url{http://adn.biol.umontreal.ca/~numericalecology/Rcode/}
#'
#'
#' \code{\link{detect_genomic_format}}
#'
#'
#' \code{\link{tidy_genomic_data}}
#'
#' \href{https://github.com/thierrygosselin/grur}{grur}
#'
#' \code{\link{read_strata}}

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

genomic_converter <- function(
  data,
  strata = NULL,
  output = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
  ...
) {

  ## Testing
  # data
  # strata = NULL
  # output = "genind"
  # filename = NULL
  # parallel.core = parallel::detectCores() - 1
  # verbose = TRUE
  # # dots dots dots
  # path.folder <- NULL
  # keep.allele.names <- FALSE
  # vcf.metadata = TRUE
  # vcf.stats <- TRUE
  # whitelist.markers = NULL
  # filter.monomorphic = TRUE
  # filter.common.markers <- TRUE
  # internal <- TRUE
  # blacklist.id <- NULL
  # blacklist.genotypes <- NULL
  # parameters <- NULL


  # Check for specific format vs package required-----------------------------
  if ("gtypes" %in% output) radiator_packages_dep(package = "strataG")
  if ("genlight" %in% output) radiator_packages_dep(package = "adegenet")
  if ("seqarray" %in% output) {
    radiator_packages_dep(package = "SeqVarTools", cran = FALSE, bioc = TRUE)
  }
  if ("snprelate" %in% output) {
    radiator_packages_dep(package = "SNPRelate", cran = FALSE, bioc = TRUE)
  }

  # Cleanup-------------------------------------------------------------------
  # obj.keeper <- c(ls(envir = globalenv()), "res", "verbose")
  radiator_function_header(f.name = "genomic_converter", verbose = verbose)
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (verbose) message("Execution date@time: ", file.date)
  old.dir <- getwd()
  opt.change <- getOption("width")
  options(width = 70)
  timing <- radiator_tic()
  #back to the original directory and options
  on.exit(setwd(old.dir), add = TRUE)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(radiator_toc(timing, verbose = verbose), add = TRUE)
  # on.exit(rm(list = setdiff(ls(envir = sys.frame(-1L)), obj.keeper), envir = sys.frame(-1L)))
  on.exit(radiator_function_header(f.name = "genomic_converter", start = FALSE, verbose = verbose), add = TRUE)
  res <- list()

  # Function call and dotslist -------------------------------------------------
  rad.dots <- radiator_dots(
    func.name = as.list(sys.call())[[1]],
    fd = rlang::fn_fmls_names(),
    args.list = as.list(environment()),
    dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE),
    keepers = c("path.folder", "keep.allele.names",
                "whitelist.markers", "filter.common.markers",
                "filter.monomorphic", "vcf.metadata", "vcf.stats",
                "parameters", "blacklist.genotypes", "internal", "blacklist.id"),
    deprecated = c("maf.thresholds", "common.markers",
                   "max.marker","monomorphic.out", "snp.ld", "filter.call.rate",
                   "filter.markers.coverage", "filter.markers.missing",
                   "number.snp.reads",
                   "mixed.genomes.analysis", "duplicate.genomes.analysis",
                   "maf.data",
                   "hierarchical.levels", "imputation.method",
                   "pred.mean.matching", "num.tree",
                   "pop.levels", "pop.labels", "pop.select"
    ),
    verbose = FALSE
  )

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("data is missing")

  # Filename -------------------------------------------------------------------
  if (is.null(filename)) {
    filename <- stringi::stri_join("radiator_data_", file.date)
  } else {
    filename.problem <- file.exists(filename)
    if (filename.problem) {
      filename <- stringi::stri_join(filename, "_", file.date)
    }
  }

  # Folders---------------------------------------------------------------------
  path.folder <- generate_folder(
    f = path.folder,
    rad.folder = "radiator_genomic_converter",
    internal = internal,
    file.date = file.date,
    verbose = verbose)

  # write the dots file
  write_rad(
    data = rad.dots,
    path = path.folder,
    filename = stringi::stri_join("radiator_genomic_converter_args_", file.date, ".tsv"),
    tsv = TRUE,
    internal = internal,
    write.message = "Function call and arguments stored in: ",
    verbose = verbose
  )

  # radiator_parameters---------------------------------------------------------
  filters.parameters <- radiator_parameters(
    generate = TRUE,
    initiate = FALSE,
    update = FALSE,
    parameter.obj = parameters,
    path.folder = path.folder,
    file.date = file.date,
    internal = FALSE,
    verbose = verbose)

  # File type detection --------------------------------------------------------
  data.type <- detect_genomic_format(data = data)

  # Import----------------------------------------------------------------------
  # back.up strata
  strata.bk <- strata


  if (verbose) message("\nImporting data\n")
  input <- tidy_genomic_data(
    data = data,
    strata = strata.bk,
    filename = filename,
    parallel.core = parallel.core,
    whitelist.markers = whitelist.markers,
    blacklist.id = blacklist.id,
    vcf.metadata = vcf.metadata,
    vcf.stats = vcf.stats,
    keep.allele.names = keep.allele.names,
    filter.common.markers = filter.common.markers,
    filter.monomorphic = filter.monomorphic,
    path.folder = path.folder,
    parameters = filters.parameters,
    internal = TRUE,
    verbose = FALSE
  )

  if(verbose) message("\nPreparing data for output\n")

  if (!is.null(strata.bk) || rlang::has_name(input, "POP_ID")) {
    if (is.factor(input$POP_ID)) {
      pop.levels <- levels(input$POP_ID)
    } else {
      pop.levels <- unique(input$POP_ID)
    }
  }

  # Biallelic detection --------------------------------------------------------
  biallelic <- radiator::detect_biallelic_markers(data = input, verbose = verbose)

  if (!biallelic && "genlight" %in% output || !biallelic && "plink" %in% output) {
    rlang::abort("output chosen doesn't work with multi-allelic data")
  }


  # GT requirement -------------------------------------------------------------
  if (TRUE %in% (c("genepop", "hierfstat", "structure", "hzar", "gsi_sim",
                   "genepopedit", "arlequin", "bayescan") %in% output)) {
    input <- calibrate_alleles(
      data = input,
      biallelic = biallelic,
      verbose = FALSE,
      gt = TRUE
    ) %$%
      input
  }
  # overide genind when marker number > 20K ------------------------------------
  if ("genind" %in% output & biallelic) {
    # detect the number of marker
    marker.number <- length(unique(input$MARKERS))
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
        message("\nIMPORTANT genind object related question...")
        message("you have > 20 000 markers (", marker.number, ")",
                "\nDo you want the more suitable genlight object instead ? (y/n):")
        overide.genind <- as.character(readLines(n = 1))
        if (overide.genind == "y") {
          message("Using genlight object...")
          output <- stringi::stri_replace_all_fixed(
            str = output,
            pattern = "genind",
            replacement = "genlight",
            vectorize_all = FALSE
          )
        } else {
          message("Continue working on genind object")
        }
      }
    }
  }

  # Imputations-----------------------------------------------------------------
  # if (!is.null(imputation.method)) {
  #   input.imp <- radiator::radiator_imputations_module(
  #     data = input,
  #     imputation.method = imputation.method,
  #     hierarchical.levels = hierarchical.levels,
  #     num.tree = num.tree,
  #     pred.mean.matching = pred.mean.matching,
  #     random.seed = random.seed,
  #     verbose = verbose,
  #     parallel.core = parallel.core,
  #     filename = NULL
  #   ) %>%
  #     dplyr::arrange(POP_ID, INDIVIDUALS, MARKERS)
  # } # End imputations

  # OUTPUT ---------------------------------------------------------------------

  setwd(path.folder)
  if (!is.null(strata.bk) && is.vector(strata.bk)) {
    strata.bk <- file.path(old.dir, strata.bk)
  }
  # GENEPOP --------------------------------------------------------------------
  if ("genepop" %in% output) {
    if (verbose) message("Generating genepop file")
    radiator::write_genepop(
      data = input,
      pop.levels = pop.levels,
      filename = filename
    )
  } # end genepop output

  # GENEPOPEDIT --------------------------------------------------------------------
  if ("genepopedit" %in% output) {
    if (verbose) message("Generating genepopedit flatten object")
    res$genepopedit <- radiator::write_genepopedit(data = input)
  } # end genepopedit output


  # GSI_SIM --------------------------------------------------------------------
  if ("gsi_sim" %in% output) {
    if (verbose) message("Generating gsi_sim output")
    res$gsi_sim <- radiator::write_gsi_sim(data = input,
                                           pop.levels = pop.levels,
                                           strata = strata.bk,
                                           filename = filename)
  } # end gsi_sim output

  # RUBIAS --------------------------------------------------------------------
  if ("rubias" %in% output) {
    if (verbose) message("Generating rubias output")
    res$rubias <- radiator::write_rubias(data = input,
                                         strata = strata.bk,
                                         filename = filename,
                                         parallel.core = parallel.core)
  } # end rubias output

  # HapMap --------------------------------------------------------------------
  if ("hapmap" %in% output) {
    if (verbose) message("Generating hapmap output")
    radiator::write_hapmap(data = input, filename = filename)
  } # end hapmap output


  # hierfstat --------------------------------------------------------------------
  if ("hierfstat" %in% output) {
    if (verbose) message("Generating hierfstat file")
    res$hierfstat <- radiator::write_hierfstat(
      data = input,
      filename = filename
    )
  } # end hierfstat output

  # strataG --------------------------------------------------------------------
  if ("gtypes" %in% output) {
    if (!requireNamespace("strataG", quietly = TRUE)) {
      rlang::abort("strataG needed for this function to work
                 Install with install.packages('strataG')")
    }
    if (verbose) message("Generating strataG gtypes object")
    res$gtypes <- radiator::write_gtypes(
      data = input,
      write = TRUE,
      filename = filename
    )
  } # end strataG output

  # structure --------------------------------------------------------------------
  if ("structure" %in% output) {
    if (verbose) message("Generating structure file")
    radiator::write_structure(
      data = input,
      pop.levels = pop.levels,
      markers.line = TRUE,
      filename = filename
    )
  } # end structure output

  # faststructure --------------------------------------------------------------------
  if ("faststructure" %in% output) {
    if (verbose) message("Generating faststructure file")
    radiator::write_faststructure(
      data = input,
      pop.levels = pop.levels,
      filename = filename
    )
  } # end faststructure output

  # betadiv --------------------------------------------------------------------
  if ("betadiv" %in% output) {
    if (!biallelic) rlang::abort("betadiv output is currently implemented for biallelic data only")
    if (verbose) message("Generating betadiv object")
    res$betadiv <- radiator::write_betadiv(data = input)

    # if (!is.null(imputation.method)) {
    #   if (verbose) message("Generating betadiv object WITH imputations")
    #   res$betadiv.imputed <- radiator::write_betadiv(data = input.imp)
    # }
  } # end betadiv output

  # arlequin --------------------------------------------------------------------
  if ("arlequin" %in% output) {
    if (verbose) message("Generating arlequin file")
    radiator::write_arlequin(
      data = input,
      pop.levels = pop.levels,
      filename = filename
    )

    # if (!is.null(imputation.method)) {
    #   if (verbose) message("Generating arlequin file WITH imputations")
    #   radiator::write_arlequin(
    #     data = input.imp,
    #     pop.levels = pop.levels,
    #     filename = filename.imp
    #   )
    # }
  } # end arlequin output

  # GENIND ---------------------------------------------------------------------
  if ("genind" %in% output) {
    if (verbose) message("Generating adegenet genind object")
    res$genind <- radiator::write_genind(data = input, write = TRUE)

    # if (!is.null(imputation.method)) {
    #   if (verbose) message("Generating adegenet genind object WITH imputations")
    #   res$genind.imputed <- radiator::write_genind(data = input.imp)
    # }
  } # end genind

  # GENLIGHT -------------------------------------------------------------------
  if ("genlight" %in% output) {
    if (verbose) message("Generating adegenet genlight object")
    res$genlight <- radiator::write_genlight(
      data = input, biallelic = TRUE, write = TRUE)

    # if (!is.null(imputation.method)) {
    #   if (verbose) message("Generating adegenet genlight object WITH imputations")
    #   res$genlight.imputed <- radiator::write_genlight(data = input.imp,
    #                                                    biallelic = TRUE)
    # }
  } # end genlight output

  # VCF ------------------------------------------------------------------------
  if ("vcf" %in% output) {
    if (!biallelic) rlang::abort("vcf output is currently implemented for biallelic data only")
    if (verbose) message("Generating VCF file")
    radiator::write_vcf(
      data = input,
      filename = filename
    )

    # if (!is.null(imputation.method)) {
    #   if (verbose) message("Generating VCF file WITH imputations")
    #   radiator::write_vcf(
    #     data = input.imp,
    #     filename = filename.imp
    #   )
    # }
  } # end vcf output

  # PLINK ----------------------------------------------------------------------
  if ("plink" %in% output) {
    if (verbose) message("Generating PLINK tped/tfam files")
    radiator::write_plink(
      data = input,
      filename = filename
    )

    # if (!is.null(imputation.method)) {
    #   if (verbose) message("Generating PLINK tped/tfam files WITH imputations")
    #   radiator::write_plink(
    #     data = input.imp,
    #     filename = filename.imp
    #   )
    # }
  } # end plink output


  # SNPRelate ------------------------------------------------------------------
  if ("snprelate" %in% output) {
    res$snprelate <- radiator::write_snprelate(
      data = input,
      biallelic = TRUE,
      filename = filename,
      verbose = verbose
    )
    # if (!is.null(imputation.method)) {
    #   if (verbose) message("Generating SNPRelate object WITH imputations")
    #   res$snprelate.imputed <- radiator::write_snprelate(
    #     data = input.imp,
    #     biallelic = TRUE,
    #     filename = filename.imp,
    #     verbose = verbose
    #   )
    # }
  }

  # bayescan -------------------------------------------------------------------
  if ("bayescan" %in% output) {
    if (verbose) message("Generating BayeScan object")
    res$bayescan <- radiator::write_bayescan(
      data = input,
      pop.select = pop.select,
      filename = filename)
  }

  # pcadapt -------------------------------------------------------------------
  if ("pcadapt" %in% output) {
    if (verbose) message("Generating pcadapt file and object")
    res$pcadapt <- radiator::write_pcadapt(
      data = input,
      filename = filename,
      parallel.core = parallel.core
    )
  }


  # hzar -------------------------------------------------------------------
  if ("hzar" %in% output) {
    if (verbose) message("Generating HZAR file")
    res$hzar <- radiator::write_hzar(
      data = input,
      distance = NULL,
      filename = filename,
      parallel.core = parallel.core
    )

  }

  # fineradstructure -----------------------------------------------------------
  if ("fineradstructure" %in% output) {
    if (verbose) message("Generating fineradstructure file")
    message("\n\nNote: This function recently changed")
    message("check radiator news and function documentation:")
    message("?radiator::write_fineradstructure\n\n")
    res$fineradstructure <- radiator::write_fineradstructure(
      data = input, strata = strata.bk, filename = filename)
  }

  # related --------------------------------------------------------------------
  if ("related" %in% output) {
    if (verbose) message("Generating related file")
    res$related <- radiator::write_related(
      data = input, filename = filename, parallel.core = parallel.core)

    # if (!is.null(imputation.method)) {
    #   if (verbose) message("Generating related file WITH imputations")
    #   res$related.imputed <- radiator::write_related(
    #     data = input.imp, filename = filename.imp, parallel.core = parallel.core)
    # }
  }

  # stockr --------------------------------------------------------------------
  if ("stockr" %in% output) {
    if (verbose) message("Generating stockR file")
    res$stockr <- radiator::write_stockr(data = input, verbose = verbose)
  }

  # structure --------------------------------------------------------------------
  if ("maverick" %in% output) {
    if (verbose) message("Generating MavericK files")
    radiator::write_maverick(
      data = input,
      filename = filename
    )

    # if (!is.null(imputation.method)) {
    #   if (verbose) message("Generating MavericK files WITH imputations")
    #   radiator::write_maverick(
    #     data = input.imp,
    #     filename = filename.imp
    #   )
    # }
  } # end MavericK output

  # dadi -----------------------------------------------------------------------
  if ("dadi" %in% output) {
    radiator::write_dadi(
      data = input
    )
    message("\n\nNote: To use an outgroup, use radiator::write_dadi separately\n")
  }



  # Writing tidy on disk -------------------------------------------------------
  # tidy.name <- stringi::stri_join(filename, ".rad")
  # message("\nWriting tidy data set:\n", tidy.name)
  # write_rad(data = input, path = tidy.name)


  # outout results -------------------------------------------------------------
  if (verbose) {
    n.markers <- length(unique(input$MARKERS))
    if (tibble::has_name(input, "CHROM")) {
      n.chromosome <- length(unique(input$CHROM))
    } else {
      n.chromosome <- "no chromosome info"
    }
    n.individuals <- length(unique(input$INDIVIDUALS))
    if (!is.null(strata.bk)) n.pop <- length(unique(input$POP_ID))

    cat("################################### RESULTS ####################################\n")
    message("Data format of input: ", data.type)
    if (biallelic) {
      message("Biallelic data")
    } else{
      message("Multiallelic data")
    }
    message("Number of markers: ", n.markers)
    message("Number of chromosome/contig/scaffold: ", n.chromosome)
    if (!is.null(strata.bk))  message("Number of strata: ", n.pop)
    message("Number of individuals: ", n.individuals)
  }
  res$tidy.data <- input
  return(res)
} # end genomic_converter
