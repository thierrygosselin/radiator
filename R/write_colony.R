# From STACKS haplotypes file write COLONY input files to the working directory

#' @name write_colony
#' @title Write a \code{COLONY} input file
#' @description Write a \code{COLONY} input file.

#' @inheritParams radiator_common_arguments

#' @param sample.markers (number) \code{COLONY} can take a long time to run,
#' use a random subsample of your markers to speed test \code{COLONY}
#' e.g. \code{sample.markers = 500} to use only 500 randomly chosen markers.
#' Default: \code{sample.markers = NULL}, will use all markers.

#' @param allele.freq (optional, string) Allele frequency can be computed from
#' a select group.
#' e.g. \code{allele.freq = "QUE"} or \code{allele.freq = c("QUE", "ONT")}.
#' Using \code{allele.freq = "overall"} will use all the samples to compute the
#' allele frequency.
#' Default: \code{allele.freq = NULL}, will not compute allele frequency.

#' @param inbreeding (boolean) 0/1 no inbreeding/inbreeding.
#' Default: \code{inbreeding = 0}
#' @param mating.sys.males (boolean) Mating system in males.
#' 0/1 polygyny/monogyny.
#' Default: \code{mating.sys.males = 0}.
#' @param mating.sys.females (boolean) Mating system in females.
#' 0/1 polygyny/monogyny.
#' Default: \code{mating.sys.females = 0}.
#' @param clone (boolean) Should clones and duplicated individuals be inferred.
#' 0/1, yes/no. Default: \code{clone = 0}.
#' @param run.length (integer) Length of run. 1 (short), 2 (medium), 3 (long),
#' 4 (very long). Start with short or medium run and consider longer run if your
#' estimates probability are not stable or really good.
#' Default: \code{run.length = 2}.
#' @param analysis (integer) Analysis method.
#' 0 (Pairwise-Likelihood Score), 1 (Full Likelihood),
#' 2 (combined Pairwise-Likelihood Score and Full Likelihood).
#' Default: \code{analysis = 1}.
#' @param allelic.dropout Locus allelic dropout rate.
#' Default : \code{allelic.dropout = 0}.
#' @param error.rate Locus error rate.
#' Default:\code{error.rate = 0.02}.
#' @param print.all.colony.opt (logical) Should all \code{COLONY} options be printed in the file.
#'
#' \strong{This require manual curation, for the file to work directly with \code{COLONY}}.
#' Default = \code{print.all.colony.opt = FALSE}.

#' @inheritParams read_strata
# @inheritParams radiator_imputations_module
#' @inheritParams tidy_genomic_data

#' @param filename Name of the acronym for filenaming in the working directory.
#' @param random.seed (integer, optional) For reproducibility, set an integer
#' that will be used inside the function that requires randomness. With default,
#' a random number is generated and printed in the appropriate output.
#' Default: \code{random.seed = NULL}.

#' @details \strong{It is highly recommended to read (twice!) the user guide distributed with
#' \code{COLONY} to find out the details for input and output of the software.}
#'
#' Not all options are provided here.
#'
#' But to ease the process, all the required options to properly run \code{COLONY}
#' will be printed in the file written in your working directory.
#' Change the values accordingly and wisely.


#' @return A \code{COLONY} file in your working directory (2 if you selected imputations arguments...)

#' @export
#' @rdname write_colony
#' @references Jones OR, Wang J (2010) COLONY: a program for parentage and
#' sibship inference from multilocus genotype data.
#' Molecular Ecology Resources, 10, 551–555.
#' @references Wang J (2012) Computationally Efficient Sibship and
#' Parentage Assignment from Multilocus Marker Data. Genetics, 191, 183–194.


#' @seealso \code{COLONY} is available on Jinliang Wang web site
#' \url{http://www.zsl.org/science/software/colony}
#'
#'
#' \href{https://thierrygosselin.github.io/radiator/articles/rad_genomics_computer_setup.html#colony}{colony installation instructions}

#' @examples
#' \dontrun{
#' # Simplest way to run the function with a tidy dataset:
#' colony.file <- radiator::write_colony(data = "turtle.data.rad")
#' }


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

write_colony <- function(
  data,
  strata = NULL,
  sample.markers = NULL,
  pop.select = NULL,
  allele.freq = NULL,
  inbreeding = 0,
  mating.sys.males = 0,
  mating.sys.females = 0,
  clone = 0,
  run.length =2,
  analysis = 1,
  allelic.dropout = 0,
  error.rate = 0.02,
  print.all.colony.opt = FALSE,
  random.seed = NULL,
  verbose = FALSE,
  parallel.core = parallel::detectCores() - 1,
  filename = NULL,
  ...
) {
  cat("#######################################################################\n")
  cat("######################## radiator::write_colony #######################\n")
  cat("#######################################################################\n")
  timing <- proc.time()

  # Checking for missing and/or default arguments-------------------------------
  if (missing(data)) rlang::abort("Input file missing")

  # Filename -------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M") # Get date and time

  if (is.null(filename)) {
    filename <- stringi::stri_join("radiator_colony_", file.date)
  }

  # Import----------------------------------------------------------------------
  message("Importing data...")
  # File type detection----------------------------------------------------------
  data.type <- radiator::detect_genomic_format(data)

  # Import data ---------------------------------------------------------------

  if (data.type %in% c("SeqVarGDSClass", "gds.file")) {
    if (data.type == "gds.file") {
      data <- radiator::read_rad(data, verbose = FALSE)
    }
    data <- gds2tidy(gds = data, parallel.core = parallel::detectCores() - 1)
    data.type <- "tbl_df"
  } else {
    if (is.vector(data)) {
      data <- radiator::tidy_wide(data = data, import.metadata = TRUE)
    }
  }

  if (rlang::has_name(data, "STRATA") && !rlang::has_name(data, "POP_ID")) {
    data %<>% dplyr::rename(POP_ID = STRATA)
  }

  if (!is.null(strata)) {
    data <- join_strata(data = data, strata = strata, pop.id = TRUE)
  }

  if (!rlang::has_name(data, "GT")) {
    data <- calibrate_alleles(data = data, verbose = FALSE) %$% input
  }

  if (!is.null(pop.select)) {
    data %<>% dplyr::filter(POP_ID %in% pop.select)
    data %<>% filter_monomorphic(data = .)
  }

  # Subsampling markers --------------------------------------------------------
  if (!is.null(sample.markers)) {
    message("Randomly subsampling ", sample.markers, " markers...")
    # sample.markers <- 500 # test
    markers.list <- dplyr::distinct(data, MARKERS) %>%
      dplyr::sample_n(tbl = ., size = sample.markers, replace = FALSE)

    data <- suppressWarnings(dplyr::semi_join(data, markers.list, by = "MARKERS"))
    markers.list <- NULL
  }

  # Generate colony without imputations ----------------------------------------
  message("Generating COLONY file...\n")
  radiator_colony(data = data, filename = filename)

  # Imputations-----------------------------------------------------------------
  # if (!is.null(imputation.method)) {
  #   message("\nMap-independent imputations...\n\n")
  #   input.imp <- radiator::radiator_imputations_module(
  #     data = dplyr::select(.data = input, MARKERS, POP_ID, INDIVIDUALS, GT),
  #     imputation.method = imputation.method,
  #     hierarchical.levels = hierarchical.levels,
  #     num.tree = num.tree,
  #     pred.mean.matching = pred.mean.matching,
  #     random.seed = random.seed,
  #     verbose = verbose,
  #     parallel.core = parallel.core,
  #     filename = NULL
  #   )
  #   message("\n\nGenerating COLONY file WITH imputations...\n")
  #   radiator_colony(data = input.imp, filename = filename.imp)
  #
  # } # End imputations

  res <- "COLONY file(s) written in the working directory"
  message(res)
  # Imputations: colony with imputed haplotypes using Random Forest ------------
  timing <- proc.time() - timing
  message(stringi::stri_join("Computation time: ", round(timing[[3]]), " sec"))
  cat("############################## completed ##############################\n")
  return(res)
}

# radiator_colony internal function ----------------------------------------------
# @keywords internal
# @param data the tidy data to transform into a colony file
# @param filename name of the file written in the directory
# @param ... other argument passed to the remaining part of the code

radiator_colony <- function(
  data,
  filename,
  allele.freq = NULL,
  inbreeding = 0,
  mating.sys.males = 0,
  mating.sys.females = 0,
  clone = 0,
  run.length =2,
  analysis = 1,
  allelic.dropout = 0,
  error.rate = 0.02,
  print.all.colony.opt = FALSE,
  ...) {
  # data <- input #test

  # Separate the alleles -------------------------------------------------------
  data <- dplyr::select(.data = data, MARKERS, POP_ID, INDIVIDUALS, GT) %>%
    dplyr::mutate(
      A1 = stringi::stri_sub(GT, 1, 3),
      A2 = stringi::stri_sub(GT, 4, 6),
      GT = NULL
    ) %>%
    tidyr::gather(data = ., key = ALLELE_GROUP, value = ALLELES, -c(MARKERS, INDIVIDUALS, POP_ID)) %>%
    dplyr::mutate(ALLELES = as.numeric(ALLELES)) %>%
    dplyr::arrange(MARKERS, POP_ID, INDIVIDUALS)

  # Allele frequency per locus--------------------------------------------------
  if (!is.null(allele.freq)) {
    message("Computing allele frequency")
    if (allele.freq != "overall") {
      input.alleles <- dplyr::filter(.data = data, POP_ID %in% allele.freq) %>%
        dplyr::filter(ALLELES != 0)
    } else {
      input.alleles <- dplyr::filter(.data = data, ALLELES != 0)
    }

    allele.per.locus <- dplyr::distinct(input.alleles, MARKERS, ALLELES) %>%
      dplyr::count(x = ., MARKERS) %>%
      dplyr::arrange(MARKERS) %>%
      dplyr::select(n) %>%
      purrr::flatten_chr(.)

    # If someone finds a beter and faster wayto do all this, I'm all in!
    freq <- input.alleles %>%
      dplyr::group_by(MARKERS, ALLELES) %>%
      dplyr::tally(.) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(FREQ = round(n/sum(n), 2)) %>%
      dplyr::select(MARKERS, ALLELES, FREQ) %>%
      dplyr::arrange(MARKERS) %>%
      dplyr::mutate(GROUP = seq(1, n(), by = 1)) %>%
      dplyr::mutate_all(.tbl = ., .funs = as.character) %>%
      tidyr::gather(data = ., key = ALLELES_FREQ, value = VALUE, -c(MARKERS, GROUP)) %>%
      dplyr::mutate(ALLELES_FREQ = factor(ALLELES_FREQ, levels = c("ALLELES", "FREQ"), ordered = TRUE)) %>%
      dplyr::group_by(MARKERS, ALLELES_FREQ) %>%
      tidyr::spread(data = ., key = GROUP, value = VALUE) %>%
      dplyr::ungroup(.) %>%
      tidyr::unite(data = ., col = INFO, -c(MARKERS, ALLELES_FREQ), sep = " ") %>%
      dplyr::mutate(
        INFO = stringi::stri_replace_all_regex(str = INFO, pattern = "NA", replacement = "", vectorize_all = FALSE),
        INFO = stringi::stri_trim_right(str = INFO, pattern = "\\P{Wspace}")
      ) %>%
      dplyr::select(-MARKERS, -ALLELES_FREQ)
    input.alleles <- NULL
  }

  markers.name <- tibble::as_tibble(t(dplyr::distinct(data, MARKERS) %>% dplyr::arrange(MARKERS)))
  marker.num <- ncol(markers.name)

  data <- tidyr::unite(data = data, col = MARKERS.ALLELE_GROUP, MARKERS, ALLELE_GROUP, sep = ".") %>%
    dplyr::group_by(POP_ID, INDIVIDUALS) %>%
    tidyr::spread(data = ., key = MARKERS.ALLELE_GROUP, value = ALLELES) %>%
    dplyr::arrange(POP_ID, INDIVIDUALS) %>%
    dplyr::ungroup(.) %>%
    dplyr::select(-POP_ID) %>%
    dplyr::mutate_all(.tbl = ., .funs = as.character)

  # Line 1 = Dataset name-------------------------------------------------------
  dataset.opt <- "`My first COLONY run`                ! Dataset name\n"
  readr::write_file(x = dataset.opt, path = filename, append = FALSE)

  # Line 2 = Output filename----------------------------------------------------
  colony.output.filename <- stringi::stri_replace_all_fixed(
    filename,
    pattern = ".dat",
    replacement = "")
  colony.output.filename <- paste(colony.output.filename, "         ! Output file name\n")
  readr::write_file(x = colony.output.filename, path = filename, append = TRUE)

  # Line 3 = Offspring number---------------------------------------------------
  off.num.opt <- paste(nrow(data), "                                  ! Number of offspring in the sample\n", sep = "")
  readr::write_file(x = off.num.opt, path = filename, append = TRUE)

  # Line 4 = Number of loci-----------------------------------------------------
  marker.num.opt <- paste(marker.num, "                                 ! Number of loci\n", sep = "")
  readr::write_file(x = marker.num.opt, path = filename, append = TRUE)

  # Line 5 = Seed random number generator  -------------------------------------
  seed.opt <- "1234                                 ! Seed for random number generator\n"
  readr::write_file(x = seed.opt, path = filename, append = TRUE)

  # Line 6 = Updating allele frequency------------------------------------------
  update.allele.freq.opt <- "0                                    ! 0/1=Not updating/updating allele frequency\n"
  readr::write_file(x = update.allele.freq.opt, path = filename, append = TRUE)

  # Line 7 = Dioecious/Monoecious species---------------------------------------
  dioecious.opt <- "2                                    ! 2/1=Dioecious/Monoecious species\n"
  readr::write_file(x = dioecious.opt, path = filename, append = TRUE)

  # Line 8 = inbreeding---------------------------------------------------------
  inbreeding.opt <- paste(inbreeding, "                                    ! 0/1=No inbreeding/inbreeding\n", sep = "")
  readr::write_file(x = inbreeding.opt, path = filename, append = TRUE)

  # Line 9 = Ploidy------------------------------------------------------------
  ploidy.opt <- "0                                    ! 0/1=Diploid species/HaploDiploid species\n"
  readr::write_file(x = ploidy.opt, path = filename, append = TRUE)

  # Line 10 = Mating system (polygamous: 0, monogamous: 1).---------------------
  mating.opt <- paste(mating.sys.males,"  ", mating.sys.females, "                                 ! 0/1=Polygamy/Monogamy for males & females\n", sep = "")
  readr::write_file(x = mating.opt, path = filename, append = TRUE)

  # Line 11 = Clone inference---------------------------------------------------
  clone.opt <- paste(clone, "                                    ! 0/1=Clone inference =No/Yes\n", sep = "")
  readr::write_file(x = clone.opt, path = filename, append = TRUE)

  # Line 12 = Sibship size scaling----------------------------------------------
  sib.size.scal.opt <- "1                                    ! 0/1=Full sibship size scaling =No/Yes\n"
  readr::write_file(x = sib.size.scal.opt, path = filename, append = TRUE)

  # Line 13 = Sibship prior indicator (Integer), average paternal sibship size (Real, optional), average maternal sibship size (Real, optional)
  sib.prior.opt <- "0 0 0                                ! 0, 1, 2, 3 = No, weak, medium, strong sibship size prior; mean paternal & maternal sibship size\n"
  readr::write_file(x = sib.prior.opt, path = filename, append = TRUE)

  # Line 14 = Population allele frequency indicator-----------------------------
  if (is.null(allele.freq)) {
    allele.freq.ind.opt <- "0                                    ! 0/1=Unknown/Known population allele frequency\n"
  } else {
    allele.freq.ind.opt <- "1                                    ! 0/1=Unknown/Known population allele frequency\n"
  }
  readr::write_file(x = allele.freq.ind.opt, path = filename, append = TRUE)


  # Line 15 and more------------------------------------------------------------
  # Numbers of alleles per locus (Integer, optional).
  # Required when the Population allele frequency indicator is set to 1
  # Alleles and their frequencies per locus (Integer, Real, optional)
  # If someone finds a beter to do all this, I'm all in!

  if (!is.null(allele.freq)) {
    readr::write_file(x = paste0(paste0(allele.per.locus, collapse = " "), "  !Number of alleles per locus\n"), path = filename, append = TRUE)
    utils::write.table(x = freq, file = filename, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  }

  # Number of runs--------------------------------------------------------------
  num.run.opt <- "1                                    ! Number of runs\n"
  readr::write_file(x = num.run.opt, path = filename, append = TRUE)

  # Length of run --------------------------------------------------------------
  # give a value of 1, 2, 3, 4 to indicate short, medium, long, very long run
  run.length.opt <- paste(run.length, "                                    ! Length of run\n", sep = "")
  readr::write_file(x = run.length.opt, path = filename, append = TRUE)

  # Monitor method (Time in second)---------------------------------------------
  monitor.met.opt <- "0                                    ! 0/1=Monitor method by Iterate\n"
  readr::write_file(x = monitor.met.opt, path = filename, append = TRUE)


  # Monitor interval (Time in second)-------------------------------------------
  monitor.int.opt <- "10000                                ! Monitor interval in Iterate\n"
  readr::write_file(x = monitor.int.opt, path = filename, append = TRUE)

  # WindowsGUI/DOS, 0 when running Colony in DOS mode or on other platforms-----
  windows.gui.opt <- "0                                    ! non-Windows version\n"
  readr::write_file(x = windows.gui.opt, path = filename, append = TRUE)

  # Analysis method : ----------------------------------------------------------
  # 0, 1 or 2 for Pairwise-Likelihood score (PLS), full likelihood method (FL),
  # or the FL and PLS combined method (FPLS). More on these methods are explained
  # above in the Windows GUI data input section.
  analysis.opt <- paste(analysis, "                                    ! Analysis 0 (Pairwise-Likelihood Score), 1 (Full Likelihood), 2 (combined Pairwise-Likelihood Score and Full Likelihood)\n", sep = "")
  readr::write_file(x = analysis.opt, path = filename, append = TRUE)

  # Precision-------------------------------------------------------------------
  precision.opt <- "3                                    ! 1/2/3=low/medium/high Precision for Full likelihood\n"
  readr::write_file(x = precision.opt, path = filename, append = TRUE)

  # Marker IDs/Names------------------------------------------------------------
  # snp.id <- seq(from = 1, to = marker.number, by = 1)
  # markers.name.opt <- as.character(c(markers.name, "        ! Marker IDs\n"))
  # readr::write_file(x = markers.name.opt, path = filename, append = TRUE)
  markers.name.opt <- stringi::stri_join(markers.name, collapse = " ")
  markers.name.opt <- stringi::stri_join(markers.name.opt, "! Marker IDs\n")
  readr::write_file(x = markers.name.opt, path = filename, append = TRUE)

  # Marker types ---------------------------------------------------------------
  # marker.type (codominant/dominant)
  marker.type.opt <- stringi::stri_join(rep(0, marker.num), collapse = " ")
  marker.type.opt <- stringi::stri_join(marker.type.opt, "  ! Marker types, 0/1=Codominant/Dominant\n")
  readr::write_file(x = marker.type.opt, path = filename, append = TRUE)

  # Allelic dropout rates-------------------------------------------------------
  dropout <- stringi::stri_join(rep(allelic.dropout, marker.num), collapse = " ")
  dropout <- stringi::stri_join(dropout, "     ! Allelic dropout rate at each locus\n")
  readr::write_file(x = dropout, path = filename, append = TRUE)

  # Error rates-----------------------------------------------------------------
  error <- stringi::stri_join(rep(error.rate, marker.num), collapse = " ")
  error <- stringi::stri_join(error, "     ! False allele rate\n\n")
  readr::write_file(x = error, path = filename, append = TRUE)

  # Offspring IDs and genotype--------------------------------------------------
  utils::write.table(x = data, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)

  # Probabilities that the father and mother of an offspring are included-------
  # in the candidate males and females. The two numbers must be provided even if
  # there are no candidate males or/and females.
  prob.opt <- "\n\n0  0                                 ! Prob. of dad/mum included in the candidates\n"
  readr::write_file(x = prob.opt, path = filename, append = TRUE)

  # Numbers of candidate males and females--------------------------------------
  candidate.opt <- "0 0                                  ! Numbers of candidate males & females\n"
  readr::write_file(x = candidate.opt, path = filename, append = TRUE)

  # PRINT ALL REMAINING COLONY OPTIONS -----------------------------------------
  if (print.all.colony.opt) {
    message("Printing all COLONY options...")
    # Candidate male IDs/names and genotypes------------------------------------
    candidate.male.id.geno.opt <- "!Candidate male ID and genotypes\n"
    readr::write_file(x = candidate.male.id.geno.opt, path = filename, append = TRUE)
    # Candidate female IDs/names and genotypes ----------------------------------
    candidate.female.id.geno.opt <- "!Candidate female ID and genotypes\n"
    readr::write_file(x = candidate.female.id.geno.opt, path = filename, append = TRUE)
  }

  # Number of offspring with known paternity------------------------------------
  known.paternity.opt <- "0                                    ! Number of offspring with known father\n"
  readr::write_file(x = known.paternity.opt, path = filename, append = TRUE)

  # Known offspring-father dyad-------------------------------------------------
  if (print.all.colony.opt) {
    known.father.dyad.opt <- "! Offspring ID and known father ID (Known offspring-father dyad)\n"
    readr::write_file(x = known.father.dyad.opt, path = filename, append = TRUE)
  }

  # Number of offspring with known maternity------------------------------------
  known.maternity.opt <- "0                                    ! Number of offspring with known mother\n"
  readr::write_file(x = known.maternity.opt, path = filename, append = TRUE)

  # Known offspring-mother dyad-------------------------------------------------
  if (print.all.colony.opt) {
    known.mother.dyad.opt <- "! Offspring ID and known mother ID (Known offspring-mother dyad)\n"
    readr::write_file(x = known.mother.dyad.opt, path = filename, append = TRUE)
  }

  # Number of known paternal sibships (Integer)---------------------------------
  known.paternal.sibships.opt <- "0                                    ! Number of known paternal sibships\n"
  readr::write_file(x = known.paternal.sibships.opt, path = filename, append = TRUE)


  # Paternal sibship size and members (Integer, String, optional).--------------
  if (print.all.colony.opt) {
    known.paternal.sibships.size.opt <- "! Paternal sibship size and members\n"
    readr::write_file(x = known.paternal.sibships.size.opt, path = filename, append = TRUE)
  }

  # Number of known maternal sibships (Integer)---------------------------------
  known.maternal.sibships.opt <- "0                                    ! Number of known maternal sibships\n"
  readr::write_file(x = known.maternal.sibships.opt, path = filename, append = TRUE)


  # Maternal sibship size and members (Integer, String, optional).--------------
  if (print.all.colony.opt) {
    known.maternal.sibships.size.opt <- "! Maternal sibship size and members\n"
    readr::write_file(x = known.maternal.sibships.size.opt, path = filename, append = TRUE)
  }

  # Number of offspring with known excluded paternity (Integer). ---------------
  offspring.known.excl.paternity.opt <- "0                                    ! Number of offspring with known excluded fathers\n"
  readr::write_file(x = offspring.known.excl.paternity.opt, path = filename, append = TRUE)


  # Excluded paternity ---------------------------------------------------------
  if (print.all.colony.opt) {
    excl.paternity.opt <- "! Offspring ID, number of excluded fathers, and excluded father IDs\n"
    readr::write_file(x = excl.paternity.opt, path = filename, append = TRUE)
  }

  # Number of offspring with known excluded maternity (Integer). ---------------
  offspring.known.excl.maternity.opt <- "0                                    ! Number of offspring with known excluded mothers\n"
  readr::write_file(x = offspring.known.excl.maternity.opt, path = filename, append = TRUE)


  # Excluded maternity----------------------------------------------------------
  if (print.all.colony.opt) {
    excl.maternity.opt <- "! Offspring ID, number of excluded mothers, and excluded father IDs\n"
    readr::write_file(x = excl.maternity.opt, path = filename, append = TRUE)
  }

  # Number of offspring with known excluded paternal sibships-------------------
  offspring.known.excl.paternal.sibships.opt <- "0                                    ! Number of offspring with known excluded paternal sibships\n"
  readr::write_file(x = offspring.known.excl.paternal.sibships.opt, path = filename, append = TRUE)

  # Excluded paternal siblings--------------------------------------------------
  if (print.all.colony.opt) {
    excluded.paternal.siblings.opt <- "! Excluded paternal siblings\n"
    readr::write_file(x = excluded.paternal.siblings.opt, path = filename, append = TRUE)
  }

  # Number of offspring with known excluded maternal sibships-------------------
  offspring.known.excl.maternal.sibships.opt <- "0                                    ! Number of offspring with known excluded maternal sibships\n"
  readr::write_file(x = offspring.known.excl.maternal.sibships.opt, path = filename, append = TRUE)

  # Excluded maternal siblings--------------------------------------------------
  if (print.all.colony.opt) {
    excluded.maternal.siblings.opt <- "! Excluded maternal siblings\n"
    readr::write_file(x = excluded.maternal.siblings.opt, path = filename, append = TRUE)
  }
}# End radiator_colony





