# radiator imputations module

#' @name radiator_imputations_module
#' @title Map-independent imputations of missing genotypes
#'
#' @description Used internally in \href{https://github.com/thierrygosselin/assigner}{assigner} and
#' \href{https://github.com/thierrygosselin/radiator}{radiator} and
#' might be of interest for users.
#' The goal of this module is to provide a simple solution for
#' a complicated problem: missing genotypes in RADseq genomic datasets.
#' This function will performed \strong{map-independent imputations} of missing
#' genotypes.
#'
#' \strong{Key features:}
#'
#' \itemize{
#' \item \strong{Imputation algorithms/methods: } Random forests (on-the-fly-imputation, ),
#' Extreme gradient tree boosting,
#' Multiple Correspondence Analysis (MCA) and
#' the classic Strawman imputation
#' (~ max/mean/mode: the most frequently observed, i.e. non-missing, genotypes is used).
#' \item \strong{Hierarchical level: } Imputations conducted by populations or globally.
#' \item \strong{Haplotype/SNP approach: } Correlation among SNPs is accounted for during
#' rf and tree boosting imputation, i.e. imputation is automatically conducted by haplotype
#' when marker meta-information is avaialble (chromosome, locus and position,
#' usually from VCF files). The alternative, considers all the markers independent
#' and imputation is conducted by SNPs.
#' \item \strong{Genotype likelihood (GL): } The GL info is detected automatically
#' (GL column in FORMAT field of VCF files). Genotypes with higher likelihood
#' will have higher probability during bootstrap samples of trees in Random
#' forests.
#' Notes: (1) option only available with Random forests;
#' (2) the use of genotype likelihoods
#' in the form of normalized, phred-scaled likelihoods (PL, e.g. from GATK)
#' are not recognized, yet... it's still under development.
#' \item \strong{Predictive mean matching: } the rf option uses a fast k-nearest neighbor
#' (KNN) searching algorithms (see argument documentation and details below).
#' \item \strong{Optimized for speed: } the package
#' \href{https://github.com/imbs-hl/ranger}{ranger}
#' (see Wright and Ziegler, 2016) provides a fast C++ version
#' of the original implementation of rf from Breiman (2001).
#' The \href{https://github.com/dmlc/xgboost}{XGBoost} provides the fast C++
#' implementation for the extreme gradient tree boosting algorithm.
#' Imputations of genotypes are conducted in parallel across CPUs.
#' A progress bar is now available to see if you have time for a coffee break!
#' }
#'
#'
#' Before running this function to populate the original dataset with synthetic
#' data I highly recommend you look for patterns of missingness
#' \code{\link[grur]{missing_visualization}}
#' and explore the reasons for their presence.
#' Follow the \href{https://www.dropbox.com/s/4zf032g6yjatj0a/vignette_missing_data_analysis.nb.html?dl=0}{vignette}
#' for more info.


#' @param data A tidy genomic dataset.
#' It can be file in the working directory or
#' an object in the global environment.
#' To get a tidy dataset from various genomic format, see
#' \href{https://github.com/thierrygosselin/radiator}{radiator}
#' \code{\link[radiator]{tidy_genomic_data}}.
#' \emph{See details of this function for more info}.

#' @param imputation.method (character, optional)
#' Methods available for map-independent imputations of missing genotype
#' (see details for more info):
#'
#' \enumerate{
#' \item \code{imputation.method = "max"} Strawman imputation,
#' the most frequently observed genotypes (ties are broken at random).
#'
#' \item \code{imputation.method = "rf"} On-the-fly-imputations using
#' Random Forests algorithm.
#'
#' \item \code{imputation.method = "rf_pred"} Random Forests algorithm is used
#' as a prediction problem.
#'
#' \item \code{imputation.method = "boost"} extreme gradient boosting trees.
#'
#' \item \code{imputation.method = "mca"} Multiple Correspondence Analysis (in devel).
#'
#' \code{imputation.method = NULL} the function will stop.
#' Default: \code{imputation.method = NULL}.
#' }

#' @param hierarchical.levels (character, optional) \code{c("global", "strata")}.
#' Should the imputations be computed by markers globally or by strata.
#' Historically, this was \code{"populations"}.
#'
#' Note that imputing genotype globally in conjunction with
#' \code{imputation.method = "max"} can potentially create huge bias.
#' e.g. by introducing foreign genotypes/haplotypes in some populations
#' (see note for more info).
#' Default: \code{hierarchical.levels = "strata"}.
#' @param num.tree (integer, optional) The number of trees to grow
#' when \code{imputation.method = "rf"} or \code{imputation.method = "rf_pred"}.
#' Default: \code{num.tree = 50}.

#' @param pred.mean.matching (integer, optional) Used in conjunction with
#' random Forests (\code{imputation.method = "rf_pred"}).
#' Number of candidate non-missing
#' value to sample from during the predictive mean matching step.
#' A fast k-nearest neighbor searching algorithms is used with this approach.
#' \code{pred.mean.matching = 3} will use 3 neighbors.
#' Default: \code{pred.mean.matching = 0}, avoids this step.

#' @param random.seed (integer, optional) For reproducibility, set an integer
#' that will be used to initialize the random generator. With default,
#' a random number is generated.
#' Default: \code{random.seed = NULL}.

#' @param verbose (optional, logical) When \code{verbose = TRUE}
#' the function is a little more chatty during execution.
#' Default: \code{verbose = TRUE}.

#' @param parallel.core (optional) The number of core used for parallel
#' execution when \code{imputation.method = "rf"}.
#' Markers are imputed in parallel, populations are processed sequentially.
#' Default: \code{parallel::detectCores() - 1}.

#' @param filename (optional) The function uses \code{\link[fst]{write.fst}},
#' to write the tidy data frame in
#' the working directory. The file extension appended to
#' the \code{filename} provided is \code{.rad}.
#' With default: \code{filename = NULL}, the imputed tidy data frame is
#' in the global environment only (i.e. not written in the working directory...).



#' @param ... (optional) To pass further argument for fine-tuning your
#' imputations. See details below.

#' @return The output in your global environment is the imputed tidy data frame.
#' If \code{filename} is provided, the imputed tidy data frame is also
#' written to the working directory. The original data is returned for markers
#' with \emph{all} or \emph{no} NA.

#' @details
#' \strong{Predictive mean matching:}
#'
#' Random Forests already behave like a nearest neighbor
#' classifier, with adaptive metric. Now we have the option to conduct
#' predictive mean matching on top of the prediction based missing value
#' imputation.PMM tries to raise the variance in the resulting conditional
#' distributions to a realistic level.
#' The closest k predicted values are identified by a fast
#' k-nearest neighbour approach wrapped in the package
#' \href{https://github.com/mayer79/missRanger}{missRanger}
#' Returned value correspond to the mean value.
#'
#'
#' \strong{haplotype/SNP approach:}
#'
#' The \emph{haplotype approach} is automatically used when markers meta-information
#' is detected (chromosome/CHROM, locus/ID and SNP/POS columns, usually from a VCF file).
#' Missing genotypes from SNPs on the same locus or same RADseq read is undertaken
#' simulteneously to account for the correlation of the linked SNPs. When one or
#' more SNPs on the same read/haplotype is missing, the haplotype is deleted and
#' consequently, imputation might results in different genotype for those SNPs
#' that were not missing. This approach is much safer than potentially creating
#' weird chimeras during haplotype imputations.
#' Alternatively, a \emph{snp approach} is used, and the SNP are considered
#' independent. Imputations of genotypes is then conducted for each marker separately.
#'
#'
#' \strong{Imputing globally or by populations ?}
#' \code{hierarchical.levels = "global"} argument will act differently depending
#' on the \code{imputation.method} selected.
#'
#' \strong{Strawman imputations (~ max/mean/mode) considerations: }
#'
#' With \code{imputation.method = "max"} and \code{hierarchical.levels = "global"}
#' \emph{will likely create bias}.
#'
#' \emph{Example 1 (unbalanced sample size):} Consider 2 populations evolving more
#' by drift than selection: pop1 (n = 36) and pop2 (n = 50).
#' You'll likely have a few polymorphic marker, where pop1 and pop2 are
#' monomorphic for different alleles (pop1 is fixed for the minor/ALT allele and
#' pop2 is fixed for the major/REF allele). Missing genotypes in pop1
#' using the most common filling technique in the literature (using mean/mode/max),
#' will result in pop1 having individuals with the REF allele.
#' Not something you want... unless your population membership is not 100% accurate,
#' (e.g. you might have migrants or wrong assignation),
#' which in this case you still don't want to impute with
#' \code{imputation.method = "max"} (see alternative below).
#'
#' \emph{Example 2 (balanced sample size):} pop1 (n = 100) and pop2 (n = 100).
#' For a particular marker, pop1 as 85 individuals genotyped and pop2 100.
#' Again, if the populations are fixed for different alleles
#' (pop1 = ALT and pop2 = REF), you will end up having REF allele in your pop1,
#' not something you want... unless your population membership is not 100% accurate,
#' (e.g. you might have migrants or wrong assignation),
#' which in this case you still don't want to impute with
#' \code{imputation.method = "max"} (see alternative below).
#'
#' \strong{Random Forests imputations: }
#'
#' Random Forests use machine learning and you can take this into account while
#' choosing argument values. Uncertain of the groupings ? Use random forests with
#' \code{hierarchical.levels = "global"}. Random forests will account for the
#' potential linkage and correlation between
#' markers and genotypes to make the best imputation available. This can potentially
#' results in genotypes for a certain combo population/marker with new groupings
#' (e.g. a new allele). This is much more accurate and not the same thing as
#' the \code{imputation.method = "max"} because the imputed genotype was validated
#' after considering all the other genotype values of the individual being imputed.
#' \emph{Test the option and report bug if you find one.}
#'
#' \emph{random forest with on-the-fly-imputation (rf): }the technique is described
#' in Tang and Ishwaran (2017). Non-missing genotypes are used for
#' the split-statistics. Daughter node assignation membership use random
#' non-missing genotypes from the inbag data. Missing genotypes are imputed at
#' terminal nodes using maximal class rule with out-of-bag non-missing genotypes.
#'
#' \emph{random forest as a prediction problem (rf_pred): }markers with
#' missing genotypes are imputed one at a time. The fitted forest is used to
#' predict missing genotypes. Missingness in the response variables are
#' incorporated as attributes for growing the forest.
#'
#' \strong{... :dot dot dot arguments}
#'
#' The argument is available to tailor your imputations using
#' extreme gradient tree boosting and random forest:
#'
#' Available arguments for extreme gradient tree boosting tree method:
#' \emph{eta, gamma, max_depth, min_child_weight, subsample, colsample_bytree,
#' num_parallel_tree, nrounds, save_name, early_stopping_rounds}.
#' Refer to \code{\link[xgboost]{xgboost}} for arguments documentation.
#'
#'
#' Available arguments for Random forests method:
#' \emph{nodesize, nsplit, nimpute}.
#' Refer to \code{\link[randomForestSRC]{impute.rfsrc}} for arguments documentation.
#'
#'
#  Multiple Correspondence Analysis option available (upcomming):
# \emph{ncp}.
# Refer to \code{\link[missMDA]{imputeMCA}} for argument documentation.

#' @note
#'
#' \strong{Reference genome or linkage map available ?}
#'
#' Numerous approaches are available and more appropriate, please search
#' the literature
#' (\href{https://online.papersapp.com/collections/05d6e65a-73c9-49e6-9c75-289a818f76f3/share}{references}).
#'
#'
#' \strong{What's simple imputation message when running the function ?}
#'
#' Before conducting the imputations by populations with random forest or extreme
#' gradient tree boosting, the data is first screened for markers that are
#' monomorphic within populations. Because for those cases, it's clear what the
#' missing genotypes should be, the imputations is very \emph{simple} and missing
#' genotypes are imputed with the only genotype found for the particular population.
#' The small cost in time is worth it, because the random forest or extreme
#' gradient tree boosting model will benefit having more complete and
#' reliable genotypes.
#'
#'
#' \strong{Deprecated arguments:}
#'
#' \itemize{
#' \item \code{imputations.group} is now replaced by \code{hierarchical.levels}
#' \item \code{impute} is no longer available.
#' Imputing using \code{impute = "allele"} option was wrong because it
#' was using F1 genotypes for imputations. Now imputation is only conducted at
#' the genotype level.
#' \item \code{iteration.rf} is no longer used. This argument is now available
#' inside the \code{...} for on-the-fly-imputations (see details). The default
#' is now set to 10 iterations.
#' \item \code{split.number} is automatically set.
#' }

#' @seealso
#' \href{https://github.com/mayer79/missRanger}{missRanger}
#'
#' \href{https://github.com/imbs-hl/ranger}{ranger}
#'
#' \href{https://github.com/stekhoven/missForest}{missForest}
#'
#' \href{https://github.com/kogalur/randomForestSRC}{randomForestSRC}
#'
#' \href{https://github.com/dmlc/xgboost}{XGBoost}

#' @export
#' @rdname radiator_imputations_module
#' @importFrom parallel detectCores
#' @importFrom dplyr distinct group_by ungroup rename arrange tally filter select select_ one_of mutate mutate_all summarise left_join funs bind_rows
#' @importFrom tidyr gather unite drop_na
#' @importFrom purrr map flatten keep discard flatten_chr flatten_dbl flatten_lgl
#' @importFrom purrrlyr invoke_rows
#' @importFrom stringi stri_replace_na
#' @importFrom tibble has_name as_data_frame
#' @importFrom stats predict reformulate as.formula
#' @importFrom rlang .data
#' @importFrom ranger ranger
#' @importFrom xgboost xgb.DMatrix cb.early.stop xgb.train
#' @importFrom randomForestSRC impute.rfsrc
#' @importFrom readr write_lines write_tsv
# @importFrom fst write.fst

#' @examples
#' \dontrun{
#' # The simplest way to run when you have a tidy dataset:
#'
#' wolf.imputed <- radiator::radiator_imputations_module(data = "wolf.tidy.dataset.tsv")
#'
#' # This will impute the missing genotypes by population using random Forests.
#' # The remaining arguments will be the defaults.
#'
#' # When you start with a vcf file you can use magrittr %>% to `pipe` the
#' # result. Below, an example with more arguments offered by the functions:
#'
#' wolf.imp <- radiator::tidy_genomic_data(
#'     data = "batch_1.vcf",
#'     strata = "strata.wolf.10pop.tsv",
#'     vcf.metadata = TRUE,
#'     whitelist.markers = "whitelist.loci.txt",
#'     verbose = TRUE) %>%
#' radiator::radiator_imputations_module(
#'     data = ., imputation.method = "boost", parallel.core = 32)
#' }

#' @references Wright, M. N. & Ziegler, A. (2016).
#' ranger: A Fast Implementation of Random Forests for High Dimensional Data
#' in C++ and R.
#' Journal of Statistical Software, in press. http://arxiv.org/abs/1508.04409.
#' @references Breiman, L. (2001). Random forests. Machine learning, 45(1), 5-32.
#' @references Chen T, Guestrin C. (2016).
#' XGBoost: A scalable tree boosting system. arXivorg. 2016.
#' doi:10.1145/2939672.2939785
#' @references Tang F, Ishwaran H. (2017) Random Forest Missing Data Algorithms.
#' arXiorg: 1–24.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

radiator_imputations_module <- function(
  data,
  imputation.method = NULL,
  hierarchical.levels = "strata",
  # markers.linkage = "multivariate",
  num.tree = 50,
  pred.mean.matching = 0,
  verbose = TRUE,
  parallel.core = parallel::detectCores() - 1,
  random.seed = NULL,
  filename = NULL,
  ...
) {
  timing <- proc.time() #for timing

  if (verbose) {
    cat("\n\n")
    cat("#######################################################################\n")
    cat("####################### grur::grur_imputations ########################\n")
    cat("#######################################################################\n")
  }
  if (is.null(imputation.method)) {
    message("Imputation method: NULL")
    message("Returning the tidy dataset")
    input.imp <- data
  } else {
  message("Imputation method: ", imputation.method)
  message("Hierarchical levels: ", hierarchical.levels)
  # message("Markers linkage: ", markers.linkage, "\n")

  # Capture unevaluated ...
  # Inspired by Eric Anderson and Hadley Wickham codes
  # could also importFrom pryr::named_dots

  dotslist <- list(...)
  unknowned_param <- setdiff(
    names(dotslist),
    c("eta", "gamma", "max_depth", "min_child_weight", "subsample", "colsample_bytree",
      "num_parallel_tree", "nrounds", "save_name", "early_stopping_rounds",
      "nodesize", "nsplit", "nimpute", "ncp"))

  if (length(unknowned_param) > 0) {
    stop("Unknowned \"...\" parameters ",
         stringi::stri_join(unknowned_param, collapse = " "),
         " to radiator imputation module")
  }

  boost.dots <- dotslist[names(dotslist) %in%
                           c("eta", "gamma", "max_depth", "min_child_weight",
                             "subsample", "colsample_bytree",
                             "num_parallel_tree", "nrounds", "save_name",
                             "early_stopping_rounds", "nodesize", "nsplit",
                             "nimpute", "ncp")]

  # learning rate
  if (!is.null(boost.dots[["eta"]])) {
    eta <- boost.dots[["eta"]]
  } else {
    eta = 0.1
  }

  # gamma for minimum loss reduction (larger the number, more conservative the algorithm)
  # prevent overfitting through regularization
  if (!is.null(boost.dots[["gamma"]])) {
    gamma <- boost.dots[["gamma"]]
  } else {
    gamma = 0
  }

  # max_depth
  if (!is.null(boost.dots[["max_depth"]])) {
    max_depth <- boost.dots[["max_depth"]]
  } else {
    max_depth = 6
  }

  # min_child_weight: minimum sum of instance weight (hessian) needed in a child.
  if (!is.null(boost.dots[["min_child_weight"]])) {
    min_child_weight <- boost.dots[["min_child_weight"]]
  } else {
    min_child_weight = 1
  }

  # subsample: subsample ratio of the training instance for growing tree and prevent overfitting
  if (!is.null(boost.dots[["subsample"]])) {
    subsample <- boost.dots[["subsample"]]
  } else {
    subsample = 0.5
  }

  # colsample_bytree: subsample ratio of columns when growing each tree
  if (!is.null(boost.dots[["colsample_bytree"]])) {
    colsample_bytree <- boost.dots[["colsample_bytree"]]
  } else {
    colsample_bytree = 1
  }

  # num_parallel_tree: Experimental parameter.
  # number of trees to grow per round.
  if (!is.null(boost.dots[["num_parallel_tree"]])) {
    num_parallel_tree <- boost.dots[["num_parallel_tree"]]
  } else {
    num_parallel_tree = 1
  }

  # nrounds: the max number of iterations
  if (!is.null(boost.dots[["nrounds"]])) {
    nrounds <- boost.dots[["nrounds"]]
  } else {
    nrounds = 2000
  }
  # save_name: the name or path for periodically saved model file
  if (!is.null(boost.dots[["save_name"]])) {
    save_name <- boost.dots[["save_name"]]
  } else {
    save_name = "imputation.model.temp"
  }

  # early_stopping_rounds: If NULL, the early stopping function is not triggered.
  # If set to an integer k, training with a validation set will stop
  # if the performance doesn't improve for k rounds.
  if (!is.null(boost.dots[["early_stopping_rounds"]])) {
    early_stopping_rounds <- boost.dots[["early_stopping_rounds"]]
  } else {
    early_stopping_rounds = 20
  }


  # randomForestSRC arguments:
  if (!is.null(boost.dots[["nodesize"]])) {
    nodesize <- boost.dots[["nodesize"]]
  } else {
    nodesize = 1
  }
  if (!is.null(boost.dots[["nsplit"]])) {
    nsplit <- boost.dots[["nsplit"]]
  } else {
    nsplit = 10
  }
  if (!is.null(boost.dots[["nimpute"]])) {
    nimpute <- boost.dots[["nimpute"]]
  } else {
    nimpute = 10
  }

  if (!is.null(boost.dots[["ncp"]])) {
    ncp <- boost.dots[["ncp"]]
  } else {
    ncp = 2
  }

  if (verbose) {
    if (imputation.method == "boost") {
      message("Extreme gradient tree boosting options:")
      message("    learning rate: ", eta)
      message("    regularization, minimum loss reduction (gamma): ", gamma)
      message("    maximum depth of tree: ", max_depth)
      message("    minimum sum of instance weight: ", min_child_weight)
      message("    subsample ratio for training when growing trees (prevent overfitting): ", subsample)
      message("    subsample ratio of columns when growing trees: ", colsample_bytree)
      message("    number of trees to grow per round: ", num_parallel_tree)
      message("    max number of iterations: ", nrounds)
      message("    filename for the model for periodical saving: ", save_name)
      message("    early stopping round integer: ", early_stopping_rounds, "\n")
    }

    if (imputation.method == "rf") {
      message("On-the-fly-imputations options:")
      message("    number of trees to grow: ", num.tree)
      message("    minimum terminal node size: ", nodesize)
      message("    non-negative integer value used to specify random splitting: ", nsplit)
      message("    number of iterations: ", nimpute)
    }

    if (imputation.method == "rf_pred") {
      message("Random Forests options:")
      message("    number of trees: ", num.tree)
      message("    minimum terminal node size: ", nodesize)
      message("    non-negative integer value used to specify random splitting: ", nsplit)
      message("    number of iterations: ", nimpute)
      message("    predictive mean matching: ", pred.mean.matching, "\n")
    }

    if (imputation.method == "mca") {
      message("Multiple Correspondence Analysis options:")
      message("    number of dimentions used to predict the missing values: ", ncp)
      message("    MCA algorithm: Regularized")
    }

    message("Number of CPUs: ", parallel.core)
    message("Note: If you have speed issues: follow radiator's vignette on parallel computing\n")
    if (!is.null(filename)) message("Filename: ", filename)
  }
  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file missing")

  # Set seed for sampling reproducibility
  if (is.null(random.seed)) {
    random.seed <- sample(x = 1:1000000, size = 1)
    set.seed(random.seed)
  } else {
    set.seed(random.seed)
  }

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    input <- radiator::tidy_wide(data = data, import.metadata = TRUE)
  } else {
    input <- data
  }
  data <- NULL #unused object

  message("\nNumber of populations: ", dplyr::n_distinct(input$POP_ID))
  message("Number of individuals: ", dplyr::n_distinct(input$INDIVIDUALS))
  message("Number of markers: ", dplyr::n_distinct(input$MARKERS))

  # output the proportion of missing genotypes BEFORE imputations
  na.before <- dplyr::summarise(.data = input, MISSING = round(length(GT[GT == "000000"])/length(GT), 6)) %>%
    purrr::flatten_dbl(.) %>% format(., scientific = FALSE)
  message("\nProportion of missing genotypes before imputations: ", na.before)

  # necessary steps to make sure we work with unique markers and not duplicated LOCUS
  if (!tibble::has_name(input, "MARKERS") && tibble::has_name(input, "LOCUS")) {
    input <- dplyr::rename(.data = input, MARKERS = LOCUS)
  }

  # New simple id for markers, because
  # formula in RF difficult with markers containing separators and numbers
  if (imputation.method %in% c("rf", "boost")) {
    marker.meta <- dplyr::distinct(.data = input, MARKERS) %>%
      dplyr::mutate(NEW_MARKERS = stringi::stri_join("M", seq(1, nrow(.))))
    input <- dplyr::full_join(marker.meta, input, by = "MARKERS")
    if (tibble::has_name(input, "CHROM") && tibble::has_name(input, "LOCUS") && tibble::has_name(input, "POS")) {
      marker.meta <- dplyr::distinct(.data = input, NEW_MARKERS, MARKERS, CHROM, LOCUS, POS)
    }
    input <- dplyr::select(.data = input, -MARKERS) %>%
      dplyr::rename(MARKERS = NEW_MARKERS)
  } else {
    # not sure necessary for max, need checking
    # scan for the columnn CHROM and keep the info to include back after imputations
    if (tibble::has_name(input, "CHROM") && tibble::has_name(input, "LOCUS") && tibble::has_name(input, "POS")) {
      marker.meta <- dplyr::distinct(.data = input, MARKERS, CHROM, LOCUS, POS)
    } else {
      marker.meta <- dplyr::distinct(.data = input, MARKERS)
    }
  }

  # scan for REF allele column
  if (tibble::has_name(input, "REF")) {
    ref.column <- TRUE
  } else {
    ref.column <- FALSE
  }
  biallelic <- radiator::detect_biallelic_markers(data = input)


  # Manage Genotype Likelihood -------------------------------------------------
  if (tibble::has_name(input, "GL")) {
    if (verbose) message("\nGenotype likelihood (GL) column detected")
  }
  have <- colnames(input)

  # For haplotype VCF
  if (!biallelic && ref.column) {
    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "POP_ID", "INDIVIDUALS", "GT_VCF_NUC", "GL")
    selected.columns <- purrr::keep(.x = have, .p = have %in% want)

    input <- dplyr::select(.data = input,
                           dplyr::one_of(selected.columns)) %>%
      dplyr::mutate(GT = replace(GT_VCF_NUC, which(GT_VCF_NUC == "./."), NA)) %>%
      dplyr::select(-GT_VCF_NUC)
  } else {
    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "POP_ID", "INDIVIDUALS", "GT", "GL")
    selected.columns <- purrr::keep(.x = have, .p = have %in% want)

    input <- dplyr::select(.data = input,
                           dplyr::one_of(selected.columns)) %>%
      dplyr::mutate(GT = replace(GT, which(GT == "000000"), NA))
  }
  have <- want <- selected.columns <- NULL
  # keep stratification
  strata.before <- dplyr::distinct(.data = input, INDIVIDUALS, POP_ID)

  # SNP/haplotype approach -----------------------------------------------------
  # detect the presence of SNP/LOCUS info and combine SNPs on the same RADseq locus together

  if (tibble::has_name(input, "CHROM") && tibble::has_name(input, "LOCUS") && imputation.method != "max") {
    input <- tidyr::unite(data = input, col = CHROM_LOCUS, CHROM, LOCUS) %>%
      dplyr::select(-POS) # no longer necessary info in MARKERS

    # Locus with > 1 SNP
    snp.number <- dplyr::distinct(.data = input, MARKERS, CHROM_LOCUS) %>%
      dplyr::count(CHROM_LOCUS)

    if (max(snp.number$n) > 1) {
      if (verbose) message("Encoding SNPs into haplotypes: combining SNPs into groups defined by chromosomes and locus info")

      if (tibble::has_name(input, "GL")) {
        keep.gl <- dplyr::ungroup(input) %>%
          dplyr::filter(!is.na(GL)) %>%
          dplyr::distinct(CHROM_LOCUS, POP_ID, INDIVIDUALS, GL) %>%
          dplyr::group_by(CHROM_LOCUS, POP_ID, INDIVIDUALS) %>%
          dplyr::summarise(GL = mean(GL))
        input <- dplyr::select(.data = input, -GL)
      } else {
        keep.gl <- NULL
      }

      locus.multiple.snp <- dplyr::filter(.data = snp.number, n > 1) %>%
        dplyr::select(CHROM_LOCUS) %>%
        purrr::flatten_chr(.)
      data.multiple.snp <- dplyr::filter(.data = input, CHROM_LOCUS %in% locus.multiple.snp)
      data.one.snp <- dplyr::filter(.data = input, !CHROM_LOCUS %in% locus.multiple.snp)
      if (length(locus.multiple.snp) > 100) {
        # parallel
        input <- list()
        input <- .radiator_parallel_mc(
          X = locus.multiple.snp,
          FUN = radiator::rad_encoding_snp,
          mc.cores = parallel.core,
          data = data.multiple.snp
        ) %>% dplyr::bind_rows(.) %>%
          dplyr::bind_rows(data.one.snp)
      } else {
        input <- purrr::map(.x = locus.multiple.snp,
                            .f = radiator::rad_encoding_snp,
                            data = data.multiple.snp) %>%
          dplyr::bind_rows(.) %>%
          dplyr::bind_rows(data.one.snp)
      }

      # include GL back and use relative measure group_by locus
      if (!is.null(keep.gl)) {
        input <- dplyr::filter(.data = input, !is.na(GT)) %>%
          dplyr::select(CHROM_LOCUS, POP_ID, INDIVIDUALS) %>%
          dplyr::left_join(keep.gl, by = c("CHROM_LOCUS", "POP_ID", "INDIVIDUALS")) %>%
          dplyr::right_join(input, by = c("CHROM_LOCUS", "POP_ID", "INDIVIDUALS")) %>%
          dplyr::select(MARKERS, CHROM_LOCUS, POP_ID, INDIVIDUALS, GT, GL) %>%
          dplyr::group_by(CHROM_LOCUS) %>%
          dplyr::mutate(GL = GL/max(GL, na.rm = TRUE)) %>%
          dplyr::ungroup(.)
      }
      separate.haplo <- TRUE
      data.multiple.snp <- data.one.snp <- snp.number <- keep.gl <- NULL
    } else {
      separate.haplo <- FALSE
    }
    # End of grouping SNPs

    input <- dplyr::select(.data = input, -CHROM_LOCUS)
  } else {
    separate.haplo <- FALSE
  }#End preparing SNP/haplo approach

  # Strawman imputations (max) -------------------------------------------------
  if (imputation.method == "max") {
    if (hierarchical.levels == "strata") {
      if (verbose) message("Using the most observed genotype per marker/strata for imputations")
      if (tibble::has_name(input, "GL")) {
        input.imp <- dplyr::select(input, MARKERS, POP_ID, INDIVIDUALS, GT, GL) %>%
          dplyr::group_by(MARKERS, POP_ID) %>%
          dplyr::mutate(
            GT = stringi::stri_replace_na(GT, replacement = max(GT, na.rm = TRUE)),
            GT = replace(GT, which(GT == "NA"), NA),
            GL = stringi::stri_replace_na(GL, replacement = mean(GL, na.rm = TRUE)),
            GL = replace(GL, which(GL == "NA"), NA),
            GL = as.numeric(GL)) %>%
          dplyr::ungroup(.)
      } else {
        input.imp <- dplyr::select(input, MARKERS, POP_ID, INDIVIDUALS, GT) %>%
          dplyr::group_by(MARKERS, POP_ID) %>%
          dplyr::mutate(GT = stringi::stri_replace_na(GT, replacement = max(GT, na.rm = TRUE)),
                        GT = replace(GT, which(GT == "NA"), NA)) %>%
          dplyr::ungroup(.)
      }
      input <- NULL
      # detect remaining NA
      # e.g. if one strata is missing all GT... when not using common markers
      # if (anyNA(data.one.snp)) {
      #   warning("Missing data is still present in the dataset",
      #           "\n    2 options:",
      #           "\n    run the function again with hierarchical.levels = 'global'",
      #           "\n    use common.markers = TRUE when using hierarchical.levels = 'strata'")
      # }
    }# End imputation max populations

    # global
    if (hierarchical.levels == "global") {
      if (verbose) message("Using the most observed genotype per marker for imputations")
      if (tibble::has_name(input, "GL")) {
        input.imp <- dplyr::select(input, MARKERS, POP_ID, INDIVIDUALS, GT, GL) %>%
          dplyr::group_by(MARKERS) %>%
          dplyr::mutate(
            GT = stringi::stri_replace_na(str = GT, replacement = max(GT, na.rm = TRUE)),
            GL = as.numeric(stringi::stri_replace_na(str = GL, replacement = mean(GL, na.rm = TRUE)))
          ) %>%
          dplyr::ungroup(.)
      } else {
        input.imp <- dplyr::select(input, MARKERS, POP_ID, INDIVIDUALS, GT) %>%
          dplyr::group_by(MARKERS) %>%
          dplyr::mutate(GT = stringi::stri_replace_na(str = GT, replacement = max(GT, na.rm = TRUE))) %>%
          dplyr::ungroup(.)
      }
      input <- NULL
    } # End imputation max global
  }#End imputation max

  # Imputation with Random Forests and tree boosting ---------------------------
  ### Note to myself: Need to add CHROM hierarchy (by markers inside CHROM, one at a time)

  if (imputation.method %in% c("rf", "boost", "mca")) {
    # Vector of markers
    marker.list <- dplyr::distinct(input, MARKERS) %>% purrr::flatten_chr(.)

    # Simple imputation with monomorphic markers -------------------------------

    # Problem encountered:
    # Merged SNPs might end up be missing if one of the SNP on the read is missing
    # It's usually safer to delete the whole genotype and impute it back,
    # instead of creating weird chimeras

    # In some dataset I've seen markers becomming monomorphic after this change
    # These markers have very very very low polymorphism and further test are
    # required to validate this technique. I think that better filtering (MAF, etc.)
    # can remove those problem...

    if (hierarchical.levels == "strata") {
      # First: dont' waist time imputing, some screening first
      # The small cost in time is worth it,
      # because model in RF and xgboost will benefit having more complete and reliable genotypes
      if (verbose) message("Scanning dataset for population(s) with monomorphic marker(s)...")
      # scanning for populations with one genotype group
      scan.pop <- dplyr::group_by(.data = input, MARKERS, POP_ID, GT) %>%
        dplyr::tally(.)

      markers.pop.na <- dplyr::filter(.data = scan.pop, is.na(GT)) %>%
        dplyr::select(MARKERS, POP_ID)

      if (nrow(markers.pop.na) > 1) {
        simple.imputation <- dplyr::filter(.data = scan.pop, !is.na(GT)) %>%
          dplyr::select(-n) %>%
          dplyr::group_by(MARKERS, POP_ID) %>%
          dplyr::tally(.) %>%
          dplyr::filter(n == 1) %>%
          dplyr::select(MARKERS, POP_ID) %>%
          dplyr::semi_join(markers.pop.na, by = c("MARKERS", "POP_ID"))

        simple.imputation.number <- nrow(simple.imputation)

        if (simple.imputation.number > 1) {
          if (verbose) message("    Simple strawman imputations conducted on ", simple.imputation.number, " markers/pops combo")
          # update marker.list
          marker.list <- dplyr::anti_join(markers.pop.na, simple.imputation, by = c("MARKERS", "POP_ID")) %>%
            dplyr::ungroup(.) %>%
            dplyr::distinct(MARKERS) %>%
            purrr::flatten_chr(.)

          data.imp <- dplyr::inner_join(
            input, simple.imputation, by = c("MARKERS", "POP_ID")) %>%
            dplyr::group_by(MARKERS, POP_ID) %>%
            dplyr::mutate(
              GT = stringi::stri_replace_na(
                GT, replacement = max(GT, na.rm = TRUE))) %>%
            dplyr::ungroup(.)

          # if GL is present give the mean value for the imputed genotype
          # note to myself: update doc to say what you're doing here...
          if (tibble::has_name(data.imp, "GL")) {
            data.imp <- data.imp %>%
              dplyr::group_by(MARKERS, POP_ID) %>%
              dplyr::mutate(GL = as.numeric(stringi::stri_replace_na(GL, replacement = mean(GL, na.rm = TRUE)))) %>%
              dplyr::ungroup(.)
          }
          # update dataframe
          input <- dplyr::anti_join(
            input, simple.imputation, by = c("MARKERS", "POP_ID")) %>%
            dplyr::bind_rows(data.imp) %>%
            dplyr::arrange(MARKERS, POP_ID, INDIVIDUALS)

          # data.imp.bk <- dplyr::filter(input, !MARKERS %in% marker.list)

          # if (nrow(data.imp.bk) <= 0) {
          #   data.imp.bk <- NULL
          # }
        }
        # else {
        #   data.imp.bk <- NULL
        # }
      }
      # else {
      #   data.imp.bk <- NULL
      # }
      # removed unused object
      data.imp <- simple.imputation <- simple.imputation.number <- scan.pop <- markers.pop.na <- NULL
    }# End imputation prep by pop

    if (hierarchical.levels == "global") {
      scan.markers.na <- dplyr::filter(input, is.na(GT)) %>%
        dplyr::distinct(MARKERS) %>%
        purrr::flatten_chr(.)

      if (length(scan.markers.na) < length(marker.list)) {

        simple.imputation <- dplyr::group_by(.data = input, MARKERS, GT) %>%
          dplyr::tally(.) %>%
          dplyr::filter(!is.na(GT)) %>%
          dplyr::select(-n) %>%
          dplyr::group_by(MARKERS) %>%
          dplyr::tally(.) %>%
          dplyr::filter(n == 1) %>%
          dplyr::select(MARKERS) %>%
          purrr::flatten_chr(.)

        simple.imputation.number <- length(simple.imputation)
        if (simple.imputation.number >= 1) {
          if (verbose) message("Simple imputations conducted on ", simple.imputation.number, " markers")

          data.imp <- dplyr::filter(input, MARKERS %in% simple.imputation) %>%
            dplyr::group_by(MARKERS) %>%
            dplyr::mutate(
              GT = stringi::stri_replace_na(GT, replacement = max(GT, na.rm = TRUE))) %>%
            dplyr::ungroup(.)

          # if GL is present give the mean value for the imputed genotype
          if (tibble::has_name(data.imp, "GL")) {
            data.imp <- data.imp %>%
              dplyr::group_by(MARKERS) %>%
              dplyr::mutate(
                GL = as.numeric(stringi::stri_replace_na(GL, replacement = mean(GL, na.rm = TRUE)))) %>%
              dplyr::ungroup(.)
          }

          input <- dplyr::filter(input, !MARKERS %in% simple.imputation) %>%
            dplyr::bind_rows(data.imp)

          data.imp <- NULL

          marker.list <- dplyr::filter(input, is.na(GT)) %>%
            dplyr::distinct(MARKERS) %>%
            purrr::flatten_chr(.)

          # data.imp.bk <- dplyr::filter(input, !MARKERS %in% marker.list)
        } else {
          marker.list <- scan.markers.na
          # data.imp.bk <- dplyr::filter(input, !MARKERS %in% marker.list)
        }
      }
      # else {
      # data.imp.bk <- NULL
      # }
      scan.markers.na <- NULL
    }

    # On-the-fly-imputations using Random Forests ------------------------------
    if (imputation.method == "rf") {
      if (verbose) message("On-the-fly-imputations using Random Forests algorithm")
      # Parallel computations options
      options(rf.cores = parallel.core, mc.cores = parallel.core)

      # on-the-fly-imputations with randomForestSRC package
      impute_rf <- function(
        x,
        num.tree = 10,
        nodesize = 1,
        # splitrule = "random",
        nsplit = 10,
        nimpute = 10,
        verbose = FALSE,
        hierarchical.levels = "strata") {

        if (hierarchical.levels == "strata") {
          message("        Imputations for pop: ", unique(x$POP_ID))
          x <- dplyr::select(x, -POP_ID)
        }

        res <- randomForestSRC::impute.rfsrc(
          data = data.frame(x),
          ntree = num.tree,
          nodesize = nodesize,
          # splitrule = "random",
          nsplit = nsplit, #split.number,
          nimpute = nimpute, #iteration.rf,
          do.trace = verbose)
        return(res)
      } # End on-the-fly imputation function

      # Random Forest by pop
      if (hierarchical.levels == "strata") {
        message("    Imputations computed by strata, take a break...")

        input.imp <- dplyr::select(input, MARKERS, POP_ID, INDIVIDUALS, GT) %>%
          dplyr::group_by(POP_ID, INDIVIDUALS) %>%
          tidyr::spread(data = ., key = MARKERS, value = GT) %>%
          dplyr::ungroup(.) %>%
          dplyr::mutate_all(.tbl = ., .funs = factor) %>%
          split(x = ., f = .$POP_ID) %>%
          purrr::map(.x = ., .f = impute_rf,
                     num.tree = num.tree, nodesize = nodesize, nsplit = nsplit,
                     nimpute = nimpute,
                     verbose = FALSE,
                     hierarchical.levels = "strata") %>%
          dplyr::bind_rows(.) %>%
          dplyr::mutate_all(.tbl = ., .funs = as.character) %>%
          tidyr::gather(data = ., key = MARKERS, value = GT, -INDIVIDUALS) %>%
          # Reintroduce the stratification (check if required)
          dplyr::right_join(strata.before, by = "INDIVIDUALS") %>%
          dplyr::arrange(MARKERS, POP_ID, INDIVIDUALS)
      }#End RF by pop

      # Random Forests global
      if (hierarchical.levels == "global") { # Globally/overall
        if (verbose) message("Imputations computed globally, take a break...")

        input.imp <- dplyr::select(input, MARKERS, INDIVIDUALS, GT) %>%
          dplyr::group_by(INDIVIDUALS) %>%
          tidyr::spread(data = ., key = MARKERS, value = GT) %>%
          dplyr::ungroup(.) %>%
          dplyr::mutate_all(.tbl = ., .funs = factor)

        input.imp <- impute_rf(
          x = input.imp,
          num.tree = num.tree, nodesize = nodesize, nsplit = nsplit,
          nimpute = nimpute, verbose = FALSE,
          hierarchical.levels = "global") %>%
          dplyr::mutate_all(.tbl = ., .funs = as.character) %>%
          tidyr::gather(data = ., key = MARKERS, value = GT, -INDIVIDUALS) %>%
          # Reintroduce the stratification (check if required)
          dplyr::right_join(strata.before, by = "INDIVIDUALS") %>%
          dplyr::arrange(MARKERS, POP_ID, INDIVIDUALS)
      } #End RF global


      # separate the haplotypes/snp group
      if (separate.haplo) {
        if (verbose) message("Decoding haplotypes: separating SNPs on the same locus and chromosome, back to original data format")
        input.imp <- radiator::rad_decoding_haplotypes(
          data = input.imp, parallel.core = parallel.core)
      }
    }# End RF

    # Random Forest imputation as a prediction problem -------------------------
    if (imputation.method == "rf_pred") {
      if (verbose) message("Using Random Forests algorithm as a prediction problem, take a break...")

      if (hierarchical.levels == "strata") {
        input.imp <- purrr::map(.x = input,
                                .f = radiator::radiator_imputer,
                                hierarchical.levels = hierarchical.levels,
                                num.tree = num.tree,
                                pred.mean.matching = pred.mean.matching,
                                random.seed = random.seed,
                                parallel.core = parallel.core) %>%
          dplyr::bind_rows(.)
      }

      # Random Forests global
      if (hierarchical.levels == "global") { # Globally/overall
        # if (verbose) message("Imputations computed globally, take a break...")
        input.rf.imp <- list() # to store results
        input.rf.imp <- radiator::radiator_imputer(data = input,
                                           hierarchical.levels = hierarchical.levels,
                                           num.tree = num.tree,
                                           pred.mean.matching = pred.mean.matching,
                                           random.seed = random.seed,
                                           parallel.core = parallel.core)
      } # End imputation RF global


    }# End rf_pred

    # Extreme Gradient Boosting Tree Imputations -------------------------------
    if (imputation.method == "boost") {
      if (verbose) message("Using extreme gradient tree boosting algorithm, take a break...")

      if (hierarchical.levels == "strata") {
        input <- dplyr::ungroup(input) %>%
          dplyr::arrange(MARKERS, POP_ID, INDIVIDUALS) %>%
          dplyr::mutate(
            # GT_N = as.numeric(factor(replace(GT, which(is.na(GT)), 0))),
            POP_ID_N = as.numeric(factor(POP_ID)),
            INDIVIDUALS_N = as.numeric(factor(INDIVIDUALS))
          ) %>%
          dplyr::group_by(MARKERS) %>%
          dplyr::mutate(GT_N = radiator::rad_factorize_gt(GT)) %>%
          dplyr::ungroup(.)

        input.wide <- dplyr::select(
          .data = input,
          MARKERS, POP_ID = POP_ID_N, INDIVIDUALS = INDIVIDUALS_N, GT = GT_N) %>%
          dplyr::arrange(MARKERS, POP_ID, INDIVIDUALS) %>%
          dplyr::group_by(POP_ID, INDIVIDUALS) %>%
          tidyr::spread(data = ., key = MARKERS, value = GT) %>%
          dplyr::ungroup(.)
      }

      if (hierarchical.levels == "global") {
        input <- dplyr::ungroup(input) %>%
          dplyr::arrange(MARKERS, POP_ID, INDIVIDUALS) %>%
          dplyr::select(-POP_ID) %>%
          dplyr::mutate(INDIVIDUALS_N = as.numeric(factor(INDIVIDUALS))) %>%
          dplyr::group_by(MARKERS) %>%
          dplyr::mutate(GT_N = radiator::rad_factorize_gt(GT)) %>%
          dplyr::ungroup(.)

        input.wide <- dplyr::select(
          .data = input,
          MARKERS, INDIVIDUALS = INDIVIDUALS_N, GT = GT_N) %>%
          dplyr::arrange(MARKERS, INDIVIDUALS) %>%
          dplyr::group_by(INDIVIDUALS) %>%
          tidyr::spread(data = ., key = MARKERS, value = GT) %>%
          dplyr::ungroup(.)
      }



      # XGBoost not adapted yet to weight variable by GL...yet
      # if (tibble::has_name(input, "GL")) {
      #   gl.wide <- dplyr::select(
      #     .data = input,
      #     MARKERS, POP_ID = POP_ID_N, INDIVIDUALS = INDIVIDUALS_N, GL) %>%
      #     dplyr::arrange(MARKERS, POP_ID, INDIVIDUALS) %>%
      #     dplyr::group_by(POP_ID, INDIVIDUALS) %>%
      #     tidyr::spread(data = ., key = MARKERS, value = GL) %>%
      #     dplyr::ungroup(.)
      # } else {
      gl.wide <- NULL
      # }

      # with map # serial test
      # input.imp <- purrr::map(.x = marker.list,
      #                         .f = radiator::radiator_boost_imputer,
      #                         input.wide = input.wide,
      #                         gl.wide = gl.wide,
      #                         eta = eta,
      #                         gamma = gamma,
      #                         max_depth = max_depth,
      #                         min_child_weight = min_child_weight,
      #                         subsample = subsample,
      #                         colsample_bytree = colsample_bytree,
      #                         num_parallel_tree = num_parallel_tree,
      #                         parallel.core = parallel.core,
      #                         nrounds = nrounds,
      #                         early_stopping_rounds = early_stopping_rounds,
      #                         save_name = save_name,
      #                         hierarchical.levels = hierarchical.levels
      #  ) %>% dplyr::bind_rows(.))


      # parallel
      input.imp <- list()
      input.imp <- .radiator_parallel_mc(
        X = marker.list,
        FUN = radiator::radiator_boost_imputer,
        mc.cores = parallel.core,
        input.wide = input.wide,
        gl.wide = gl.wide,
        eta = eta,
        gamma = gamma,
        max_depth = max_depth,
        min_child_weight = min_child_weight,
        subsample = subsample,
        colsample_bytree = colsample_bytree,
        num_parallel_tree = num_parallel_tree,
        parallel.core = parallel.core,
        nrounds = nrounds,
        early_stopping_rounds = early_stopping_rounds,
        save_name = save_name,
        hierarchical.levels = hierarchical.levels
      ) %>% dplyr::bind_rows(.)
      input.wide <- NULL# remove unused objects

      # Remove factor/integer type genotype (revert back to original)
      input.imp <- radiator::rad_defactorize_gt(data.to.change = input.imp, data.with.info = input)

      # Reintroduce the stratification (check if required)
      input.imp <- dplyr::left_join(strata.before, input.imp, by = "INDIVIDUALS") %>%
        dplyr::arrange(MARKERS, POP_ID, INDIVIDUALS)
      strata.before <- NULL# remove unused objects

      if (separate.haplo) {
        # separate the haplotypes/snp group
        if (verbose) message("Decoding haplotypes: separating SNPs on the same locus and chromosome, back to original data format")
        input.imp <- radiator::rad_decoding_haplotypes(
          data = input.imp, parallel.core = parallel.core)
      }

      input <- NULL# remove unused objects

      # Impute GL
      if (tibble::has_name(input.imp, "GL") && hierarchical.levels == "strata") {
        message("Imputing GL with mean value per populations")
        input.imp <- dplyr::group_by(.data = input.imp, MARKERS, POP_ID, GT) %>%
          dplyr::mutate(
            GL = stringi::stri_replace_na(GL, replacement = mean(GL, na.rm = TRUE)),
            GL = replace(GL, which(GL %in% c("NA", "NaN")), NA),
            GL = as.numeric(GL)) %>%
          dplyr::group_by(MARKERS, GT) %>%
          dplyr::mutate(
            GL = as.numeric(stringi::stri_replace_na(GL, replacement = mean(GL, na.rm = TRUE)))
          ) %>%
          dplyr::ungroup(.)
      }
      if (tibble::has_name(input.imp, "GL") && hierarchical.levels == "global") {
        message("Imputing GL with mean overall value")
        input.imp <- dplyr::group_by(.data = input.imp, MARKERS, GT) %>%
          dplyr::mutate(
            GL = as.numeric(stringi::stri_replace_na(GL, replacement = mean(GL, na.rm = TRUE)))
          ) %>%
          dplyr::ungroup(.)
      }
    }# End boost

    # Multiple Correspondence Analysis Imputations -----------------------------
    if (imputation.method == "mca") {
      if (verbose) message("Using Multiple Correspondence Analysis algorithm...")

      impute_mca <- function(x, ncp = 2) {
        # res <- missMDA::imputeMCA(data.frame(x), ncp = ncp)$completeObs
        res <- rep(1, nrow(x)) #test
        return(res)
      } # End impute_mca

      input <- dplyr::select(input, MARKERS, POP_ID, INDIVIDUALS, GT) %>%
        dplyr::mutate(GT = replace(GT, which(GT == "000000"), NA)) %>%
        dplyr::group_by(POP_ID, INDIVIDUALS) %>%
        tidyr::spread(data = ., key = MARKERS, value = GT) %>%
        dplyr::ungroup(.) %>%
        dplyr::mutate_all(.tbl = ., .funs = factor)

      # test <- impute_mca(x = input, ncp = 2)


      # Random Forest by pop
      if (hierarchical.levels == "strata") {
        message("    Imputations computed by strata, take a break...")

        input.split <- split(x = input, f = input$POP_ID)
        input.imp <- list()
        input.imp <- .radiator_parallel_mc(
          X = input.split,
          FUN = impute_mca,
          mc.cores = parallel.core,
          ncp = ncp
        ) %>%
          dplyr::bind_rows(.) %>%
          dplyr::mutate_all(.tbl = ., .funs = as.character) %>%
          tidyr::gather(data = ., key = MARKERS, value = GT, -INDIVIDUALS) %>%
          # Reintroduce the stratification (check if required)
          dplyr::right_join(strata.before, by = "INDIVIDUALS") %>%
          dplyr::arrange(MARKERS, POP_ID, INDIVIDUALS)
      }#End RF by pop

      # Random Forests global
      if (hierarchical.levels == "global") { # Globally/overall
        if (verbose) message("Imputations computed globally, take a break...")

        input.imp <- dplyr::select(input, MARKERS, INDIVIDUALS, GT) %>%
          dplyr::group_by(INDIVIDUALS) %>%
          tidyr::spread(data = ., key = MARKERS, value = GT) %>%
          dplyr::ungroup(.) %>%
          dplyr::mutate_all(.tbl = ., .funs = factor)

        input.imp <- impute_rf(
          x = input.imp,
          num.tree = num.tree, nodesize = nodesize, nsplit = nsplit,
          nimpute = nimpute,verbose = FALSE,
          hierarchical.levels = "global") %>%
          dplyr::mutate_all(.tbl = ., .funs = as.character) %>%
          tidyr::gather(data = ., key = MARKERS, value = GT, -INDIVIDUALS) %>%
          # Reintroduce the stratification (check if required)
          dplyr::right_join(strata.before, by = "INDIVIDUALS") %>%
          dplyr::arrange(MARKERS, POP_ID, INDIVIDUALS)
      } #End RF global


      # separate the haplotypes/snp group
      if (separate.haplo) {
        if (verbose) message("Decoding haplotypes: separating SNPs on the same locus and chromosome, back to original data format")
        input.imp <- radiator::rad_decoding_haplotypes(
          data = input.imp, parallel.core = parallel.core)
      }

    }# End mca

  } # End imputation RF,  boost and MCA

  # prep results ---------------------------------------------------------------

  # Replace NA by 000000 in GT column if found
  if (anyNA(input.imp)) {
    warning("Missing data is still present in the dataset",
            "\n    2 options:",
            "\n    run the function again with hierarchical.levels = 'global'",
            "\n    use common.markers = TRUE when using hierarchical.levels = 'strata'")
    if (!biallelic && ref.column) {
      input.imp$GT <- stringi::stri_replace_na(str = input.imp$GT, replacement = "./.")
    } else {
      input.imp$GT <- stringi::stri_replace_na(str = input.imp$GT, replacement = "000000")
    }
  }

  if (!biallelic && ref.column) {
    input.imp <- dplyr::rename(input.imp, GT_VCF_NUC = GT)
  }

  # Compute REF/ALT allele... might have change depending on prop of missing values
  if (verbose) message("Adjusting REF/ALT alleles to account for imputations...")
  input.imp <- radiator::change_alleles(
    data = input.imp,
    biallelic = biallelic,
    parallel.core = parallel.core,
    verbose = verbose)$input

  if (tibble::has_name(input.imp, "POLYMORPHIC.x")) input.imp <- dplyr::select(input.imp, -POLYMORPHIC.x)
  if (tibble::has_name(input.imp, "POLYMORPHIC.y")) input.imp <- dplyr::select(input.imp, -POLYMORPHIC.y)


  # Integrate marker.meta columns and sort
  if (!is.null(marker.meta)) {
    want <- c( "MARKERS", "CHROM", "LOCUS", "POS", "POP_ID",
               "INDIVIDUALS", "REF", "ALT", "GT", "GT_VCF",
               "GT_VCF_NUC", "GT_BIN", "GL")

    if (tibble::has_name(marker.meta, "NEW_MARKERS")) {
      input.imp <- suppressWarnings(dplyr::left_join(
        dplyr::rename(input.imp, NEW_MARKERS = MARKERS),
        marker.meta, by = "NEW_MARKERS") %>%
          dplyr::select(dplyr::one_of(want)))
    } else {
      input.imp <- suppressWarnings(dplyr::left_join(
        input.imp, marker.meta, by = "MARKERS") %>%
          dplyr::select(dplyr::one_of(want)))
    }
    want <- marker.meta <- NULL
  } else {
    input.imp <- dplyr::arrange(.data = input.imp, MARKERS, POP_ID, INDIVIDUALS)
  }

  # Write to working directory
  if (!is.null(filename)) {
    tidy.name <- stringi::stri_join(filename, ".rad")
    if (verbose) message("Writing the imputed tidy data: \n", tidy.name)
    #fst::write.fst(x = input.imp, path = tidy.name, compress = 85)
    # readr::write_tsv(x = input.imp, path = filename, col_names = TRUE)
  }

  # Missing after imputation:
  na.after <- dplyr::summarise(.data = input.imp, MISSING = round(length(GT[GT == "000000"])/length(GT), 6)) %>%
    purrr::flatten_dbl(.) %>% format(., scientific = FALSE)
  message("\nProportion of missing genotypes after imputations: ", na.after)

  # Error notices
  if (imputation.method == "boost" && file.exists("radiator_imputations_error.txt")) {
    message("Error notice: Tree boosting imputations encountered error(s),
            please look in the working directory for a file:
            radiator_imputations_error.txt
            email the problem to the author: thierrygosselin@icloud.com")
  }
}
  if (verbose) {
    # output the proportion of missing genotypes after imputations
    timing <- proc.time() - timing
    message("\nComputation time: ", round(timing[[3]]), " sec")
    cat("################## grur::grur_imputations completed ###################\n")
  }
  return(input.imp)
  } # End imputations

# Internal nested Function -----------------------------------------------------
# radiator_imputer ---------------------------------------------------------------
#' @title radiator_imputer
#' @description imputations using Ranger package and predictive mean matching
#' @rdname radiator_imputer
#' @keywords internal
#' @export

radiator_imputer <- function(
  data,
  num.tree = 100,
  pred.mean.matching = 0,
  random.seed = NULL,
  parallel.core = parallel::detectCores() - 1,
  # markers.linkage = "multivariate",
  hierarchical.levels = "strata",
  marker.list = marker.list,
  verbose = verbose
) {
  # data <- input #test

  data <- data %>%
    dplyr::select(POP_ID, INDIVIDUALS, MARKERS, GT) %>%
    dplyr::group_by(INDIVIDUALS, POP_ID) %>%
    dplyr::mutate(GT = replace(GT, which(is.na(GT)), "missing")) %>%
    tidyr::spread(data = ., key = MARKERS, value = GT) %>%
    dplyr::ungroup(.)

  # data[is.na(data.model)] <- "missing"


  data.gl <- NULL
  # if (tibble::has_name(data.pop, "GL")) {
  #   # not sure useful
  #   data.gl <- data.pop %>%
  #     dplyr::select(POP_ID, INDIVIDUALS, MARKERS, GL_RF) %>%
  #     dplyr::group_by(INDIVIDUALS, POP_ID) %>%
  #     tidyr::spread(data = ., key = MARKERS, value = GL_RF) %>%
  #     dplyr::ungroup(.)
  # } else {
  #   data.gl <- NULL
  # }


  data.na <- NULL # remove after test


  # initiate while loop
  i <- 1
  pred.error <- rep(1, length(marker.list))
  names(pred.error) <- marker.list
  oob.error <- TRUE
  maxiter <- 10000
  data.imp <- tibble::data_frame(character(0))

  # imp <- list()
  while (oob.error && i <= maxiter) {
    data.last <- data
    pred.error.last <- pred.error

    data.rf <- list()
    data.rf <- .radiator_parallel_mc(
      X = marker.list,
      FUN = radiator::rad_impute_genotypes,
      mc.cores = parallel.core,
      data = data,
      data.na = data.na,
      data.gl = data.gl,
      num.tree = num.tree,
      pred.mean.matching = pred.mean.matching,
      random.seed = random.seed,
      parallel.core = parallel.core,
      hierarchical.levels = hierarchical.levels,
      # markers.linkage = markers.linkage,
      pred.error = pred.error
    ) #%>%
    # dplyr::bind_rows(.)

    system.time(test <- purrr::map(
      .x = marker.list, .f = radiator::rad_impute_genotypes,
      data = data,
      data.na = data.na,
      data.gl = data.gl,
      data.imp = data.imp,
      num.tree = num.tree,
      pred.mean.matching = pred.mean.matching,
      random.seed = random.seed,
      parallel.core = parallel.core,
      hierarchical.levels = hierarchical.levels,
      # markers.linkage = markers.linkage,
      pred.error = pred.error))


    test <- dplyr::bind_rows(data.imp)

    # update error
    oob.error <- mean(pred.error) < mean(pred.error.last)
    i <- i + 1 # update iteration
  } # End of loop

  if (i == maxiter && oob.error || i == 2) {
    imputed.dataset <- data
  } else {
    imputed.dataset <- data.last
  }

  imputed.dataset <- dplyr::mutate_all(.tbl = imputed.dataset,
                                       .funs = as.character, exclude = NA)


  # results --------------------------------------------------------------------


  return(data.imp)
} #End radiator_imputer

# rad_impute_genotypes -------------------------------------------------------------
#' @title rad_impute_genotypes
#' @description imputations using Ranger package and predictive mean matching of missRanger
#' @rdname rad_impute_genotypes
#' @keywords internal
#' @export
#' @export

rad_impute_genotypes <- function(
  marker.list,
  data,
  data.na,
  data.gl = NULL,
  data.imp,
  num.tree = 100,
  pred.mean.matching = 0,
  random.seed = NULL,
  parallel.core = parallel::detectCores() - 1,
  hierarchical.levels = "strata",
  # markers.linkage = "multivariate",
  pred.error = pred.error
) {
  # m <- "BINDED_M7114_M7115_M7116_M7117"
  # m <- "BINDED_M1_M2_M3_M4_M5"
  # m <- "BINDED_M101_M102"
  # m <- "M993"
  m <- marker.list
  message("Marker: ", m)# for diagnostic

  # Handling complete and missing data ---------------------------------------
  data.model <- dplyr::filter(.data = data, rlang::.data[[m]] != "missing") %>%
    dplyr::mutate_all(.tbl = ., .funs = factor)
  data.missing <- dplyr::filter(.data = data, rlang::.data[[m]] == "missing") %>%
    dplyr::select(-dplyr::one_of(m))

  # If all missing screening # this should be done with marker.list before all this
  # if (nrow(data.missing) > 0) {

  # GL
  data.gl <- NULL

  data.complete <- NULL # remove after test

  if (!is.null(data.gl)) {
    # mean GL per sample
    case.weights <- suppressWarnings(
      dplyr::select(.data = data.complete, INDIVIDUALS) %>%
        dplyr::left_join(data.gl, by = "INDIVIDUALS") %>%
        dplyr::ungroup(.) %>%
        dplyr::select(-dplyr::one_of(c("POP_ID", "INDIVIDUALS"))) %>%
        purrrlyr::invoke_rows(.f = purrr::lift_vd(mean), .to = "GL", .collate = "cols") %>%
        dplyr::select(GL) %>%
        purrr::flatten_dbl(.)
    )
    # mean GL per markers
    split.select.weights <- suppressWarnings(
      dplyr::select(.data = data.complete, INDIVIDUALS) %>%
        dplyr::left_join(data.gl, by = "INDIVIDUALS") %>%
        dplyr::ungroup(.) %>%
        dplyr::select(-dplyr::one_of(c(m, "POP_ID", "INDIVIDUALS"))) %>%
        dplyr::summarise_all(.tbl = ., .funs = mean) %>%
        purrr::flatten_dbl(.)
    )
  } else {
    case.weights <- NULL
    split.select.weights <- NULL
  }
  message("Data preparation: ok")# for diagnostic

  # Formula ------------------------------------------------------------------
  # if (markers.linkage == "multivariate") {
  # if (hierarchical.levels == "strata") {
  # discard.columns <- c(m, "POP_ID", "INDIVIDUALS")
  # discard.columns <- c(m, "POP_ID")
  # model.columns <- setdiff(colnames(data.complete), discard.columns)
  model.columns <- setdiff(colnames(data.model), m)
  rf.formula <- stats::as.formula(
    stringi::stri_join(m, " ~ ",
                       stringi::stri_join(model.columns, collapse = "+")))
  always.split.variables <- NULL
  always.split.variables <- "POP_ID"
  # } else {
  #   discard.columns <- c(m, "INDIVIDUALS")
  #   model.columns <- setdiff(colnames(data.complete), discard.columns)
  #   model.columns <- setdiff(colnames(data.complete), m)
  #   rf.formula <- stats::as.formula(
  #     stringi::stri_join(m, " ~ ",
  #                        stringi::stri_join(model.columns, collapse = "+")))
  #   # rf.formula <- stats::reformulate(termlabels = "POP_ID", response = m)
  #   always.split.variables <- c("POP_ID")
  # }
  # } else {#univariate (one marker at a time)
  #   if (hierarchical.levels == "strata") {
  #     rf.formula <- stats::reformulate(termlabels = ".", response = m)
  #     always.split.variables <- NULL
  #   } else {
  #     rf.formula <- stats::reformulate(termlabels = "POP_ID", response = m)
  #     always.split.variables <- c("POP_ID")
  #   }
  # }#End composing formula
  message("Formula: ok")# for diagnostic
  # RF -----------------------------------------------------------------------
  system.time(ranger.res <- ranger::ranger(
    formula = rf.formula,
    data = data.model,
    # num.trees = 1000,
    num.trees = num.tree,
    case.weights = case.weights,
    split.select.weights = split.select.weights,
    always.split.variables = always.split.variables,
    num.threads = 1,
    # num.threads = parallel.core,
    seed = random.seed))

  # ranger.res
  predicted <- stats::predict(ranger.res$forest, data.missing)$predictions

  message("RF: ok")# for diagnostic
  # predictive mean matching ---------------------------------------------
  if (pred.mean.matching > 0) {

    # ytrain <- dplyr::select(data.complete, GT) %>%
    #   dplyr::mutate(GT = as.character(GT)) %>%
    #   purrr::flatten_chr(.)

    # wide format
    # ytrain <- dplyr::select(.data = data.complete, dplyr::one_of(m)) %>%
    #   dplyr::ungroup(.) %>%
    #   dplyr::mutate_all(.tbl = ., .funs = as.character) %>%
    #   purrr::flatten_chr(.)

    # long format (doesn't give reliable results with missRanger...)
    # ytrain <- dplyr::select(.data = data.complete, GT_IMP) %>%
    #   dplyr::ungroup(.) %>%
    #   dplyr::mutate_all(.tbl = ., .funs = as.character) %>%
    #   purrr::flatten_chr(.)

    # To pass Travis
    # predicted <- missRanger::pmm(
    #   xtrain = ranger.res$predictions,
    #   xtest = predicted,
    #   ytrain = ytrain,
    #   k = pred.mean.matching)
    predicted <- NULL
  }

  message("pred.mean.matching: ok")# for diagnostic

  # data.missing[,m] <- predicted
  # data <- suppressWarnings(dplyr::bind_rows(data.complete, data.missing) %>% dplyr::arrange(INDIVIDUALS))
  # data.imp <- data[,m]


  # imp <- dplyr::select(.data = data2, dplyr::one_of(c("POP_ID", "INDIVIDUALS", m))) %>%
  #   radiator::rad_decoding_haplotypes(parallel.core = parallel.core)
  # imp[[m]] <- radiator::rad_decoding_haplotypes(data = data.imp, parallel.core = parallel.core)

  data.imp <- dplyr::select(.data = data.missing, INDIVIDUALS) %>%
    dplyr::mutate(
      MARKERS = rep(m, nrow(data.missing)),
      GT = predicted
    )



  pred.error[[m]] <- ranger.res$prediction.error
  if (is.nan(pred.error[[m]])) pred.error[[m]] <- 0

  res <- list(data = data, pred.error = pred.error, data.imp = data.imp)
  # } else {
  #   pred.error[[m]] <- 0
  #   data.imp <- data[,m]
  #   res <- list(data = data, pred.error = pred.error, data.imp = data.imp)
  # }
  # difference with missForest, missRanger and randomForestSRC:
  # other package are updating the dataset with the imputed value and continue
  # looping through the columns, creating a bias or differences between columns
  # imputed first through last as none have the same predictor columns missigness
  # e.g. missRanger: completed <- union(completed, v)
  # I think it's preferable to leave the data frame as is and merge columns in the end

  # if (separate.haplo) {
  #   imputed.dataset <- imputed.dataset %>%
  #     tidyr::separate(data = ., col = GT, into = haplo.meta, sep = "-", remove = TRUE) %>%
  #     tidyr::gather(data = ., key = MARKERS, value = GT, -c(CHROM_LOCUS, POP_ID, INDIVIDUALS))
  # }
  return(res)
  message("results: ok")# for diagnostic
} #End rad_impute_genotypes

# radiator_boost_imputer ---------------------------------------------------------------
#' @title radiator_boost_imputer
#' @description imputations using Ranger package and predictive mean matching
#' @rdname radiator_boost_imputer
#' @keywords internal
#' @export
radiator_boost_imputer <- function(
  marker.list = NULL,
  input.wide = NULL,
  gl.wide = NULL,
  eta = 0.2,
  gamma = 0,
  max_depth = 6,
  min_child_weight = 1,
  subsample = 0.8,
  colsample_bytree = 1,
  num_parallel_tree = 1,
  parallel.core = 1,
  nrounds = 200,
  early_stopping_rounds = 20,
  save_name = "imputation.model.temp",
  hierarchical.levels = "strata"
) {
  # marker.list <- "BINDED_M1003_M1004_M1005_M1006_M1007"
  m <- rlang::sym(marker.list)
  message("Imputation of marker: ", m)
  data.complete <- dplyr::filter(.data = input.wide, !is.na(rlang::UQ(m)))
  data.label <- dplyr::select(.data = data.complete, rlang::UQ(m)) %>% purrr::flatten_dbl(.)
  data.complete <- dplyr::select(.data = data.complete, -rlang::UQ(m)) %>% as.matrix(.)
  data.complete <- xgboost::xgb.DMatrix(data = data.complete, label = data.label, missing = NA)
  data.missing <- dplyr::filter(.data = input.wide, is.na(rlang::UQ(m)))

  if (nrow(data.missing) == 0) stop("code error: email author")

  res <- dplyr::select(.data = data.missing, INDIVIDUALS) %>%
    dplyr::mutate(MARKERS = rep(rlang::quo_text(m), nrow(data.missing)))
  missing.label <- dplyr::select(.data = data.missing, rlang::UQ(m)) %>% purrr::flatten_dbl(.)
  data.missing <- dplyr::select(.data = data.missing, -rlang::UQ(m)) %>% as.matrix(.)
  data.missing <- xgboost::xgb.DMatrix(data = data.missing, label = missing.label, missing = NA)

  params <- list(
    booster = "dart", silent = 0, eta = eta, gamma = gamma, max_depth = max_depth,
    min_child_weight = min_child_weight, subsample = subsample,
    colsample_bytree = colsample_bytree, num_parallel_tree = num_parallel_tree,
    objective = "multi:softmax", num_class = length(unique(data.label)), #max(data.label) + 1,
    nthread = 1)

  watchlist <- list(train = data.complete) # not an argument of xgboost::xgboost

  callbacks <- list(xgboost::cb.early.stop(
    stopping_rounds = early_stopping_rounds,
    metric_name = "train_merror", verbose = FALSE))

  # Catch error while tree boosting
  safe_boost <- purrr::safely(.f = xgboost::xgb.train)

  boost.res <- safe_boost(
    data = data.complete,
    missing = NA,
    weight = NULL,
    params = params,
    nrounds = nrounds,
    verbose = 0,
    watchlist = watchlist,
    callbacks = callbacks,
    save_period = 0,
    save_name = save_name)


  if (is.null(boost.res$error)) {
    boost.res <- boost.res$result
    res$GT <- stats::predict(boost.res, data.missing,
                             ntreelimit = boost.res$best_ntreelimit,
                             missing = NA)
  } else {
    readr::write_lines(x = boost.res$error, path = "radiator_imputations_error.txt", append = TRUE)
    res$GT <- as.numeric(rep(NA, nrow(res)))
  }

  # boost.res <- xgboost::xgb.train(
  #   data = data.complete,
  #   missing = NA,
  #   weight = NULL,
  #   params = params,
  #   nrounds = nrounds,
  #   verbose = 0,
  #   watchlist = watchlist,
  #   callbacks = callbacks,
  #   save_period = 0,
  #   save_name = save_name)

  # test
  # predicted <- stats::predict(boost.res, data.complete, ntreelimit = boost.res$best_ntreelimit, missing = 0)
  # observed <- dplyr::filter(.data = input.wide, rlang::.data[[m]] != 0) %>%
  # dplyr::select(dplyr::one_of(m)) %>% purrr::flatten_dbl(.)
  # predicted
  # observed
  # identical(predicted, observed)
  # predicted <- stats::predict(boost.res, data.missing, ntreelimit = boost.res$best_ntreelimit, missing = 0)
  # predicted
  # if (anyNA(predicted)) stop("Not enough trees were selected...")
  # data.imp[[m]] <- predicted


  # Unused arguments
  m <- data.complete <- data.label <- data.missing <- missing.label <- NULL
  params <- watchlist <- callbacks <- boost.res <- NULL
  return(res)
}#End boost

# rad_encoding_snp --------------------------------------------------------------------
#' @title rad_encoding_snp
#' @description bind snp found on the same locus
#' @rdname rad_encoding_snp
#' @keywords internal
#' @export
rad_encoding_snp <- function(locus.list = NULL, data = NULL) {
  # locus.list <- "1_135"
  res <- dplyr::filter(.data = data, CHROM_LOCUS %in% locus.list)
  binded.markers <- dplyr::distinct(.data = res, MARKERS) %>%
    purrr::flatten_chr(.) %>%
    stringi::stri_join(., collapse = "_")
  binded.markers <- stringi::stri_join("BINDED_", binded.markers)

  res <- res %>%
    dplyr::group_by(CHROM_LOCUS, POP_ID, INDIVIDUALS) %>%
    tidyr::spread(data = ., key = MARKERS, value = GT) %>%
    tidyr::unite(data = ., col = GT, -CHROM_LOCUS, -POP_ID, -INDIVIDUALS, sep = "_", remove = TRUE) %>%
    dplyr::mutate(#Haplotype with combination of SNP and NA = NA (up for an argument?)
      GT = stringi::stri_replace_all_fixed(
        str = GT, pattern = "NA", replacement = NA, vectorize_all = FALSE),
      MARKERS = rep(binded.markers, n())
    ) %>%
    dplyr::ungroup(.) %>%
    dplyr::select(MARKERS, CHROM_LOCUS, POP_ID, INDIVIDUALS, GT)

  return(res)
}#End rad_encoding_snp

# rad_decoding_haplotypes------------------------------------------------------------
#' @title rad_decoding_haplotypes
#' @description separate snp group merged with rad_encoding_snp
#' @rdname rad_decoding_haplotypes
#' @keywords internal
#' @export

rad_decoding_haplotypes <- function(data = NULL, parallel.core = parallel::detectCores() - 1) {
  # data <- data.imp.bk#test
  # data <- input.imp#test

  # find the markers that need sep.
  if (tibble::has_name(data, "MARKERS")) {
    markers.sep <- dplyr::distinct(data, MARKERS) %>%
      dplyr::filter(stringi::stri_detect_regex(str = MARKERS, pattern = "^BINDED")) %>%
      purrr::flatten_chr(.)
  } else {
    col.names.data <- colnames(data)
    markers.sep <- purrr::keep(
      .x = col.names.data,
      .p = stringi::stri_detect_regex(str = col.names.data, pattern = "^BINDED"))
    col.names.data <- NULL
  }
  # nested function required ---------------------------------------------------
  rad_separate_locus <- function(binded.markers = NULL, data = NULL) {
    # binded.markers <- markers.sep[[1]]

    col.replace <- stringi::stri_replace_all_fixed(
      str = binded.markers, pattern = "BINDED_", replacement = "", vectorize_all = FALSE)

    if (tibble::has_name(data, "MARKERS")) {
      data.sep <- dplyr::filter(.data = data, MARKERS %in% binded.markers) %>%
        dplyr::select(-MARKERS) %>%
        tidyr::separate_(
          data = .,
          col = "GT",
          into = stringi::stri_split_fixed(str = col.replace, pattern = "_", simplify = TRUE),
          sep = "_", extra = "drop")

      if (tibble::has_name(data.sep, "GL")) {
        data.sep <- tidyr::gather(data = data.sep, key = MARKERS, value = GT, -POP_ID, -INDIVIDUALS, -GL) %>%
          dplyr::select(POP_ID, INDIVIDUALS, MARKERS, GT, GL) %>%
          dplyr::arrange(POP_ID, INDIVIDUALS, MARKERS, GT, GL)
      } else {
        data.sep <- tidyr::gather(data = data.sep, key = MARKERS, value = GT, -POP_ID, -INDIVIDUALS) %>%
          dplyr::select(POP_ID, INDIVIDUALS, MARKERS, GT) %>%
          dplyr::arrange(POP_ID, INDIVIDUALS, MARKERS, GT)
      }
    } else {
      data.sep <- dplyr::select(.data = data, dplyr::one_of(c("POP_ID", "INDIVIDUALS", binded.markers)))
      colnames(data.sep) <- c("POP_ID", "INDIVIDUALS", col.replace)
      data.sep <- tidyr::separate_(
        data = data.sep,
        col = col.replace,
        into = stringi::stri_split_fixed(str = col.replace, pattern = "_", simplify = TRUE),
        sep = "_", extra = "drop") %>%
        tidyr::gather(data = ., key = MARKERS, value = GT, -POP_ID, -INDIVIDUALS)
    }

    return(data.sep)
  }#End rad_separate_locus

  if (length(markers.sep) > 0) {
    if (length(markers.sep) > 100) {
      data.sep <- list()
      data.sep <- .radiator_parallel_mc(
        X = markers.sep,
        FUN = rad_separate_locus,
        mc.cores = parallel.core,
        data = data
      ) %>% dplyr::bind_rows(.)
    } else {
      data.sep <- purrr::map(.x = markers.sep, .f = rad_separate_locus, data = data) %>%
        dplyr::bind_rows(.)
    }

    # Include markers 1 snp/read
    if (tibble::has_name(data, "MARKERS")) {
      markers.no.sep <- dplyr::distinct(data, MARKERS) %>%
        dplyr::filter(!stringi::stri_detect_regex(str = MARKERS, pattern = "^BINDED")) %>%
        purrr::flatten_chr(.)
    } else {
      col.names.data <- colnames(data)
      markers.no.sep <- purrr::discard(
        .x = col.names.data,
        .p = stringi::stri_detect_regex(str = col.names.data, pattern = "^BINDED"))
      markers.no.sep <- purrr::discard(
        .x = markers.no.sep,
        .p = markers.no.sep %in% c("POP_ID", "INDIVIDUALS"))
      col.names.data <- NULL
    }

    if (length(markers.no.sep) > 0) {
      if (tibble::has_name(data, "MARKERS")) {
        data.no.sep <- suppressWarnings(
          dplyr::filter(.data = data, MARKERS %in% markers.no.sep) %>%
            dplyr::select(dplyr::one_of(c("POP_ID", "INDIVIDUALS", "MARKERS", "GT", "GL")))
        )
      } else {
        data.no.sep <- suppressWarnings(
          dplyr::select(
            .data = data,
            dplyr::one_of(c("POP_ID", "INDIVIDUALS", markers.no.sep))) %>%
            tidyr::gather(data = ., key = MARKERS, value = GT, -POP_ID, -INDIVIDUALS)
        )
      }
    } else {
      data.no.sep <- NULL
    }

    # combined data
    if (!is.null(data.no.sep)) {
      data.sep <- dplyr::bind_rows(data.sep, data.no.sep) %>%
        dplyr::arrange(POP_ID, INDIVIDUALS, MARKERS)
    }
  } else {
    data.sep <- data
  }
  return(data.sep)
}#End rad_decoding_haplotypes

# rad_factorize_gt------------------------------------------------------------
#' @title rad_factorize_gt
#' @description Necessary to factorize by markers in tidy format.
#' XGBoost needs numbering to start at , hence the codes below.
#' @rdname rad_factorize_gt
#' @keywords internal
#' @export

rad_factorize_gt <- function(x) {
  x <- as.numeric(factor(x)) - 1
}#End rad_factorize_gt


# rad_defactorize_gt------------------------------------------------------------
#' @title rad_defactorize_gt
#' @description Function to "defactorize/decode" the imputed data back to original.
#' @rdname rad_defactorize_gt
#' @keywords internal
#' @export

rad_defactorize_gt <- function(data.to.change, data.with.info) {
  clean.id <- dplyr::distinct(.data = data.with.info, INDIVIDUALS, INDIVIDUALS_N)
  clean.gt <- dplyr::distinct(.data = data.with.info, MARKERS, GT, GT_N) %>%
    tidyr::drop_na(.)

  if (tibble::has_name(data.to.change, "POP_ID")) {
    clean.pop <- dplyr::distinct(.data = data.with.info, POP_ID, POP_ID_N)
    res <- suppressWarnings(
      dplyr::arrange(.data = data.to.change, POP_ID, INDIVIDUALS) %>%
        dplyr::rename(POP_ID_N = POP_ID, INDIVIDUALS_N = INDIVIDUALS, GT_N = GT) %>%
        dplyr::inner_join(clean.id, by = "INDIVIDUALS_N") %>%
        dplyr::select(-INDIVIDUALS_N) %>%
        dplyr::inner_join(clean.pop, by = "POP_ID_N") %>%
        dplyr::select(-POP_ID_N) %>%
        dplyr::left_join(clean.gt, by = c("MARKERS", "GT_N")) %>%
        dplyr::select(POP_ID, INDIVIDUALS, MARKERS, GT) %>%
        dplyr::bind_rows(
          tidyr::drop_na(
            data = dplyr::select(
              .data = data.with.info,
              dplyr::one_of(c("POP_ID", "INDIVIDUALS", "MARKERS", "GT", "GL"))
            ))))
  } else {
    res <- suppressWarnings(dplyr::arrange(.data = data.to.change, INDIVIDUALS) %>%
                              dplyr::rename(INDIVIDUALS_N = INDIVIDUALS, GT_N = GT) %>%
                              dplyr::inner_join(clean.id, by = "INDIVIDUALS_N") %>%
                              dplyr::select(-INDIVIDUALS_N) %>%
                              dplyr::left_join(clean.gt, by = c("MARKERS", "GT_N")) %>%
                              dplyr::select(INDIVIDUALS, MARKERS, GT) %>%
                              dplyr::bind_rows(
                                tidyr::drop_na(
                                  data = dplyr::select(
                                    .data = data.with.info,
                                    dplyr::one_of(c("INDIVIDUALS", "MARKERS", "GT", "GL"))
                                  ))))
  }
  return(res)
}#End rad_defactorize_gt
