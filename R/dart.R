# This file is for everything related to DArT
# There is one for VCF and one for GDS


# extract_dart_target_id---------------------------------------------------------
#' @name extract_dart_target_id

#' @title Extract \href{http://www.diversityarrays.com}{DArT} target id

#' @description Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users. The function allows to extract DArT
#' target id from a DArT file. To help prepare the appropriate STRATA file.

#' @inheritParams read_dart

#' @param write With default \code{write = TRUE}, the dart target id column is
#' written in a file in the working directory.

#' @return A tidy dataframe with a \code{TARGET_ID} column. For cleaning, the \code{TARGET_ID}
#' column is treated like the column \code{INDIVIDUALS}. Spaces and \code{,}
#' are removed, \code{_} and \code{:} are changed to a dash \code{-} and
#' UPPER case is used.
#' \href{https://thierrygosselin.github.io/radiator/reference/clean_ind_names.html}{see cleaning doc for logic behind this}.

#' @export
#' @rdname extract_dart_target_id

#' @examples
#' \dontrun{
#' # Built a strata file:
#' strata <- radiator::extract_dart_target_id("mt.dart.file.csv") %>%
#'     dplyr::mutate(
#'         INDIVIDUALS = "new id you want to give",
#'         STRATA = "fill this"
#'     ) %>%
#'     readr::write_tsv(x = ., file = "my.new.dart.strata.tsv")
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com} and Peter Grewe \email{peter.grewe@csiro.au}

extract_dart_target_id <- function(data, write = TRUE) {
  ##TEST
  # write = FALSE
  # write = TRUE

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("Input file missing")

  # Check DArT format file -----------------------------------------------------
  dart.check <- check_dart(data)

  # Import data ----------------------------------------------------------------
  dart.target.id <- readr::read_delim(
    file = data,
    delim = dart.check$tokenizer.dart,
    skip = dart.check$skip.number,
    n_max = 1,
    na = "-",
    col_names = FALSE,
    col_types = readr::cols(.default = readr::col_character())) %>%
    t %>%
    magrittr::set_colnames(x = ., value = "TARGET_ID") %>%
    tibble::as_tibble(.)


  if (dart.check$star.number > 0) {
    dart.target.id %<>%
      dplyr::filter(dplyr::row_number() > dart.check$star.number)
  } else {
    # This is a string of known DArT col header not wanted
    discard <- c(
      "TARGET_ID", "ALLELE_ID", "CLONE_ID", "CLUSTER_TEMP_INDEX", "ALLELE_SEQUENCE",
      "CLUSTER_CONSENSUS_SEQUENCE", "CLUSTER_SIZE", "ALLELE_SEQ_DIST", "SNP",
      "SNP_POSITION", "CALL_RATE", "ONE_RATIO_REF", "ONE_RATIO_SNP", "FREQ_HOM_REF",
      "FREQ_HOM_SNP", "FREQ_HETS", "PIC_REF", "PIC_SNP", "AVG_PIC", "AVG_COUNT_REF",
      "AVG_COUNT_SNP", "RATIO_AVG_COUNT_REF_AVG_COUNT_SNP",
      "FREQ_HETS_MINUS_FREQ_MIN_HOM",
      "ALLELE_COUNTS_CORRELATION", "AGGREGATE_TAGS_TOTAL", "DERIVED_CORR_MINUS_SEED_CORR",
      "REP_REF", "REP_SNP", "REP_AVG", "PIC_REP_REF", "PIC_REP_SNP", "TOTAL_PIC_REP_REF_TEST",
      "TOTAL_PIC_REP_SNP_TEST", "BIN_ID", "BIN_SIZE", "ALLELE_SEQUENCE_REF",
      "ALLELE_SEQUENCE_SNP", "TRIMMED_SEQUENCE_REF", "TRIMMED_SEQUENCE", "ONE_RATIO",
      "PIC", "AVG_READ_DEPTH", "STDEV_READ_DEPTH", "Q_PMR", "REPRODUCIBILITY", "MAF",
      "TOT_COUNTS", "TOTAL_PIC_REP_TEST", "PIC_REP")

    discard.genome <- c("CHROM_|CHROM_POS_|ALN_CNT_|ALN_EVALUE_")

    dart.target.id %<>%
      dplyr::filter(!radiator_snakecase(x = TARGET_ID) %in% discard) %>%
      dplyr::filter(!stringi::stri_detect_regex(
        str = stringi::stri_trans_toupper(TARGET_ID),
        pattern = discard.genome, negate = FALSE)) %>%
      dplyr::filter(!stringi::stri_detect_regex(
        str = radiator_snakecase(TARGET_ID),
        pattern = discard.genome, negate = FALSE))
  }

  dart.target.id %<>%
    dplyr::mutate(
      TARGET_ID = clean_ind_names(x = TARGET_ID),
      TARGET_ID = stringi::stri_trans_toupper(TARGET_ID)
    )
  if (write) readr::write_tsv(x = dart.target.id, file = "dart.target.id.tsv")

  # Check that DArT file as good target id written -----------------------------
  if (nrow(dart.target.id) != length(unique(dart.target.id$TARGET_ID))) {
    rlang::warn(
      "Non unique TARGET_ID or sample names used in the DArT file.
Solution: edit manually")
  }
  return(dart.target.id)
}#End extract_dart_target_id


# read_dart -------------------------------------------------------------------------
# Import, filter and transform a dart output file to different formats

#' @name read_dart

#' @title Read and tidy \href{http://www.diversityarrays.com}{DArT} output files.

#' @description Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users. The function generate a GDS object/file
#' and optionally, a tidy dataset using
#' \href{http://www.diversityarrays.com}{DArT} files.

#' @param data One of the DArT output files. 6 formats used by DArT are recognized
#' by radiator.
#' recognised:
#' \enumerate{
#' \item \code{1row}: Genotypes are in 1 row and coded (0, 1, 2, -).
#' \code{0 for 2 reference alleles REF/REF}, \code{1 for 2 alternate alleles ALT/ALT},
#' \code{2 for heterozygote REF/ALT}, \code{- for missing}.
#' \item \code{2rows}: No genotypes. It's absence/presence, 0/1, of the REF and ALT alleles.
#' Sometimes called binary format.
#' \item \code{counts}: No genotypes, It's counts/read depth for the REF and ALT alleles.
#' Sometimes just called count data.
#' \item \code{silico.dart}: SilicoDArT data. No genotypes, no REF or ALT alleles.
#' It's a file coded as absence/presence, 0/1, for the presence of sequence in
#' the clone id.
#' \item \code{silico.dart.counts}: SilicoDArT data. No genotypes, no REF or ALT alleles.
#' It's a file coded as absence/presence, with counts for the presence of sequence in
#' the clone id.
#' \item \code{dart.vcf}: For DArT VCFs, please use \code{\link{read_vcf}}.
#' }
#'
#' Depending on the number of markers, these format will be recoded similarly to
#' VCF files (dosage of alternate allele, see details).
#'
#' The function can import \code{.csv} or \code{.tsv} files.
#'
#'
#' If you encounter a problem, sent me your data so that I can update
#' the function.

#' @param strata A tab delimited file or object with 3 columns.
#' Columns header is:
#' \code{TARGET_ID}, \code{INDIVIDUALS} and \code{STRATA}.
#' Note: the column \code{STRATA} refers to any grouping of individuals.
#' You need to make sure that
#' the column \code{TARGET_ID} match the id used by DArT. With the \code{counts}
#' format the \code{TARGET_ID} is a series of integer.
#' With \code{1row} and \code{2rows} the \code{TARGET_ID} is actually the sample
#' name submitted to DArT.
#' The column \code{INDIVIDUALS} and \code{STRATA} will be kept in the tidy data.
#' Only individuals in the strata file are kept in the tidy, i.e. that the strata
#' is also used as a whitelist of individuals/strata.
#' Silico DArT data is currently used to detect sex markers, so the \code{STRATA}
#' column should be filed with sex information: \code{M} or \code{F}.
#'
#' See example on how to extract the TARGET_ID of your DArT file.
#'
#' \href{https://www.dropbox.com/s/utq2h6o00v55kep/example.dart.strata.tsv?dl=0)}{example.dart.strata.tsv}.

#' @param tidy.dart (logical, optional) Generate a tidy dataset.
#' Default:\code{tidy.dart = FALSE}.


#' @inheritParams tidy_genomic_data
#' @inheritParams radiator_common_arguments
#' @inheritParams filter_whitelist



#' @param ... (optional) To pass further argument for fine-tuning the function.

#' @return A radiator GDS file and tidy dataframe with several columns depending on DArT file:
#' \code{silico.dart:} A tibble with 5 columns: \code{CLONE_ID, SEQUENCE, VALUE, INDIVIDUALS, STRATA}.
#' This object is also saved in the directory (file ending with .rad).
#'
#' Common to \code{1row, 2rows and counts}: A GDS file is automatically generated.
#' To have a tidy tibble, the argument \code{tidy.dart = TRUE} must be used.
#'
#' \enumerate{
#' \item VARIANT_ID: generated by radiator and correspond the markers in integer.
#' \item MARKERS: generated by radiator and correspond to CHROM + LOCUS + POS separated by 2 underscores.
#' \item CHROM: the chromosome info, for de novo: CHROM_1.
#' \item LOCUS: the locus info.
#' \item POS: the SNP id on the LOCUS.
#' \item COL: the position of the SNP on the short read.
#' \item REF: the reference allele.
#' \item ALT: the alternate allele.
#' \item INDIVIDUALS: the sample name.
#' \item STRATA/POP_ID: populations id of the sample.
#' \item GT_BIN: the genotype based on the number of alternate allele in the genotype
#' (the count/dosage of the alternate allele). \code{0, 1, 2, NA}.
#' \item REP_AVG: the reproducibility average, output specific of DArT.
#' }
#' Other columns potentially in the tidy tibble:
#' \enumerate{
#' \item GT: the genotype in 6 digit format \emph{à la genepop}.
#' \item GT_VCF: the genotype in VCF format \code{0/0, 0/1, 1/1, ./.}.
#' \item GT_VCF_NUC: the genotype in VCF format, but keeping the nucleotide information.
#' \code{A/A, A/T, T/T, ./.}
#' \item AVG_COUNT_REF: the coverage for the reference allele, output specific of DArT.
#' \item AVG_COUNT_SNP: the coverage for the alternate allele, output specific of DArT.
#' \item READ_DEPTH: the number of reads used for the genotype (count data).
#' \item ALLELE_REF_DEPTH: the number of reads of the reference allele (count data).
#' \item ALLELE_ALT_DEPTH: the number of reads of the alternate allele (count data).
#' }
#'
#'
#' Written in the working directory:
#' \itemize{
#' \item The radiator GDS file
#' \item The DArT metadata information
#' \item The tidy DArT data
#' \item The strata file associated with this tidy dataset
#' \item The allele dictionary is a tibble with columns:
#' \code{MARKERS, CHROM, LOCUS, POS, REF, ALT}.
#' }

#' @section Advance mode:
#'
#' \emph{dots-dots-dots ...} allows to pass several arguments for fine-tuning the function:
#' \enumerate{
#' \item \code{whitelist.markers}: detailed in \code{\link[radiator]{filter_whitelist}}.
#' Defautl: \code{whitelist.markers = NULL}.
#' \item \code{missing.memory} (option, path)
#' This argument allows to erase genotypes that have bad statistics.
#' It's the path to a file \code{.rad} file that contains 3 columns:
#' \code{MARKERS, INDIVIDUALS, ERASE}. The file is produced by several radiator
#' functions. For DArT data, \code{\link[radiator]{filter_rad}} generate the file.
#' Defautl: \code{missing.memory = NULL}. Currently not used.
#' \item \code{path.folder}: (optional, path) To write output in a specific folder.
#' Default: \code{path.folder = NULL}. The working directory is used.
#' \item \code{pop.levels}: detailed in \code{\link[radiator]{tidy_genomic_data}}.
#' }


#' @export
#' @rdname read_dart

#' @examples
#' \dontrun{
#' clownfish.dart.tidy <- radiator::read_dart(
#'     data = "clownfish.dart.csv",
#'     strata = "clownfish.strata.tsv"
#'     )
#' }
#' @seealso \code{\link{extract_dart_target_id}}

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

read_dart <- function(
  data,
  strata,
  filename = NULL,
  tidy.dart = FALSE,
  verbose = FALSE,
  parallel.core = parallel::detectCores() - 1,
  ...
) {

  # # for testing
  # filename = NULL
  # verbose = TRUE
  # parallel.core = parallel::detectCores() - 1
  # whitelist.markers = NULL
  # missing.memory <- NULL
  # path.folder = NULL
  # internal <- FALSE
  # tidy.dart = FALSE
  # tidy.check = FALSE
  # gt = NULL
  # gt.bin = NULL
  # gt.vcf = NULL
  # gt.vcf.nuc = NULL
  # pop.levels = NULL

  # Cleanup-------------------------------------------------------------------
  # obj.keeper <- c(ls(envir = globalenv()), "data")
  radiator_function_header(f.name = "read_dart", verbose = verbose)
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (verbose) message("Execution date@time: ", file.date)
  old.dir <- getwd()
  opt.change <- getOption("width")
  options(width = 70)
  timing <- radiator_tic()
  #back to the original directory and options
  on.exit(setwd(old.dir), add = TRUE)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(radiator_toc(timing), add = TRUE)
  # on.exit(rm(list = setdiff(ls(envir = sys.frame(-1L)), obj.keeper), envir = sys.frame(-1L)))
  on.exit(radiator_function_header(f.name = "read_dart", start = FALSE, verbose = verbose), add = TRUE)

  # Function call and dotslist -------------------------------------------------
  rad.dots <- radiator_dots(
    func.name = as.list(sys.call())[[1]],
    fd = rlang::fn_fmls_names(),
    args.list = as.list(environment()),
    dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE),
    keepers = c("whitelist.markers", "missing.memory",
                "path.folder", "internal", "pop.levels",
                "gt", "gt.bin", "gt.vcf", "gt.vcf.nuc",
                "tidy.check"
    ),
    verbose = FALSE
  )
  if (is.null(tidy.check)) tidy.check <- FALSE
  if (tidy.check) tidy.dart <- TRUE
  if (!tidy.dart) gt <- gt.vcf <- gt.vcf.nuc <- FALSE

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("data is missing")

  # Folders---------------------------------------------------------------------
  wf <- path.folder <- generate_folder(
    f = path.folder,
    rad.folder = "read_dart",
    prefix_int = FALSE,
    internal = internal,
    file.date = file.date,
    verbose = verbose)


  # write the dots file
  write_rad(
    data = rad.dots,
    path = path.folder,
    filename = stringi::stri_join("radiator_tidy_dart_args_", file.date, ".tsv"),
    tsv = TRUE,
    internal = internal,
    verbose = verbose
  )

  # Generate filenames ---------------------------------------------------------
  filename.gds <- generate_filename(
    name.shortcut = filename,
    path.folder = path.folder,
    date = TRUE,
    extension = "gds.rad")

  filename <- generate_filename(
    name.shortcut = filename,
    path.folder = path.folder,
    date = TRUE,
    extension = "rad")

  meta.filename <- generate_filename(
    name.shortcut = "radiator_tidy_dart_metadata",
    path.folder = path.folder,
    date = TRUE,
    extension = "rad")

  dic.filename <- generate_filename(
    name.shortcut = "radiator_tidy_dart_allele_dictionary",
    path.folder = path.folder,
    date = TRUE,
    extension = "tsv")

  # Import data ---------------------------------------------------------------
  data <- import_dart(
    data = data,
    strata = strata,
    pop.levels = pop.levels,
    parallel.core = parallel.core,
    verbose = verbose
  )
  strata <- data$strata
  dart.format <- data$dart.format
  data <- data$data

  # Whitelist ------------------------------------------------------------------
  if (!is.null(whitelist.markers)) {
    data %<>% filter_whitelist(data = ., whitelist.markers = whitelist.markers)
  }

  # Silico DArT ----------------------------------------------------------------
  if ("silico.dart" %in% dart.format) {
    want <- c("CLONE_ID", "SEQUENCE", strata$INDIVIDUALS)
    data %<>%
      dplyr::select(tidyselect::any_of(want)) %>%
      radiator::rad_long(
        x = .,
        cols = c("CLONE_ID", "SEQUENCE"),
        names_to = "INDIVIDUALS",
        values_to = "VALUE",
        variable_factor = FALSE
      )

    n.clone <- length(unique(data$CLONE_ID))
    data <- radiator::join_strata(data = data, strata = strata)

    filename <- generate_filename(
      name.shortcut = "radiator.silico.dart",
      path.folder = path.folder,
      date = TRUE,
      extension = "rad")

    write_rad(
      data = data,
      path = filename$filename,
      filename = filename$filename.short,
      write.message = "standard",
      verbose = verbose
    )

    if (verbose) cat("################################### SUMMARY ####################################\n")
    message("\nNumber of clones: ", n.clone)
    summary_strata(strata)
    return(data)
  }


  # STRATEGY tidy dart  --------------------------------------------------------
  # Depending on the number of markers ...
  # All this can be overwritten in ... argument
  # gt.bin is the dosage of ALT allele: 0, 1, 2 NA
  ERASE <- NULL
  n.markers <- length(unique(data$MARKERS))
  gt.bin <- TRUE

  # gt.vcf is genotype coding in the VCF: 0/0, 0/1, 1/1, ./.
  if (is.null(gt.vcf)) {
    if (n.markers < 5000) gt.vcf <- TRUE
    if (n.markers >= 5000 && n.markers < 30000) gt.vcf <- TRUE
    if (n.markers >= 30000) gt.vcf <- FALSE
  }
  # gt.vcf.nuc is genotype coding in the VCF but with nucleotides: A/C, ./.
  if (is.null(gt.vcf.nuc)) {
    if (n.markers < 5000) gt.vcf.nuc <- TRUE
    if (n.markers >= 5000 && n.markers < 30000) gt.vcf.nuc <- FALSE
    if (n.markers >= 30000) gt.vcf.nuc <- FALSE
  }
  # gt is genotype coding a la genepop: 001002, 000000
  if (is.null(gt)) {
    if (n.markers < 5000) gt <- TRUE
    if (n.markers >= 5000 && n.markers < 30000) gt <- FALSE
    if (n.markers >= 30000) gt <- FALSE
  }

  # dart2gds -------------------------------------------------------------------
  # check if some rows are missing...
  if (anyNA(unique(data$REF))) {
    if (length(unique(data$MARKERS[!is.na(data$REF)])) != length(unique(data$MARKERS[is.na(data$REF)]))) {
      rlang::abort("The DArT file is missing rows...")
    }
  }

  # Markers meta
  markers.meta <- dplyr::ungroup(data) %>%
    dplyr::select(-tidyselect::any_of(strata$INDIVIDUALS)) %>%
    dplyr::filter(!is.na(REF) | !is.na(ALT)) %>%
    dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
    dplyr::arrange(VARIANT_ID) %>%
    dplyr::mutate(FILTERS = "whitelist") %>%
    dplyr::select(
      FILTERS, M_SEQ, VARIANT_ID, MARKERS, CHROM, LOCUS, POS, COL, REF, ALT,
      tidyselect::everything(.)
    )

  write_rad(
    data = markers.meta,
    path = meta.filename$filename,
    filename = meta.filename$filename.short,
    write.message = "standard",
    verbose = verbose
  )

  #GDS
  #filename = filename.gds

  data <- suppressWarnings(
    dart2gds(
      data = data,
      strata = strata,
      markers.meta = markers.meta,
      filename = filename.gds,
      dart.format = dart.format,
      gt.vcf = gt.vcf,
      gt.vcf.nuc = gt.vcf.nuc,
      gt = gt,
      parallel.core = parallel.core,
      verbose = verbose
    )
  )

  # Tidy data -----------------------------------------------------------------
  n.chrom <- length(unique(markers.meta$CHROM))
  n.locus <- length(unique(markers.meta$LOCUS))
  n.snp <- length(unique(markers.meta$MARKERS))

  # Tidy check ---------------------------------------------------------------
  if (tidy.dart && tidy.check && n.markers > 20000) {
    cat("\n\n################################## IMPORTANT ###################################\n")
    message("Tidying DArT data with ", n.markers, " SNPs is not optimal:")
    message("    1. a computer with lots of RAM is required")
    message("    2. it's very slow to generate")
    message("    3. it's very slow to run codes after")
    message("    4. for most non model species this number of markers is not realistic...")
    message("\nRecommendation:")
    message("    1. stop here and just use the GDS")
    message("    2. filter your dataset. e.g. with filter_rad")
    message("\nIdeally target a maximum of ~ 10 000 - 20 0000 unlinked SNPs\n")

    if (n.markers > 20000) tidy.dart <- FALSE
    tidy.dart <- radiator_question(
      x = "\nContinue tidying the DArT data (y/n) ?",
      answer.opt = c("Y", "N", "Yes", "No", "YES", "NO", "yes", "no", "y", "n"))
    if (any(c("y", "Y", "Yes", "YES", "yes") %in% tidy.dart)) {
      tidy.dart <- TRUE
      message("Tidying the large DArT file...")
    } else {
      message("\nKeeping only the GDS object/file")
      tidy.dart <- FALSE
    }
  }

  # Print genotypes tidying
  if (verbose) {
    message("\nGenotypes formats generated with ", n.markers, " SNPs: ")
    message("    GT_BIN (the dosage of ALT allele: 0, 1, 2 NA): ", gt.bin)
    message("    GT_VCF (the genotype coding VCFs: 0/0, 0/1, 1/1, ./.): ", gt.vcf)
    message("    GT_VCF_NUC (the genotype coding in VCFs, but with nucleotides: A/C, ./.): ", gt.vcf.nuc)
    message("    GT (the genotype coding 'a la genepop': 001002, 001001, 000000): ", gt)
  }

  if (tidy.dart) {
    notwanted <- c("REF", "ALT")
    tidy.data <- markers.meta %>%
      dplyr::left_join(
        extract_genotypes_metadata(gds = data, whitelist = TRUE) %>%
          dplyr::select(-tidyselect::any_of(notwanted))
        , by = "MARKERS"
      )

    # Check that merging was successful
    if (!identical(tidy.data$VARIANT_ID.x, tidy.data$VARIANT_ID.y)) {
      rlang::abort("Tidying problem: raise an issue on GitHub...")
    } else {
      tidy.data %<>%
        dplyr::select(-VARIANT_ID.y) %>%
        dplyr::rename(VARIANT_ID = VARIANT_ID.x)
    }

    # Calibration of ref/alt alleles ------------------------------------------
    # Now done during coding of genotypes

    # Final strata ---------------------------------------------------------
    strata.filename <- generate_filename(
      name.shortcut = "radiator_tidy_dart_strata",
      path.folder = path.folder,
      date = TRUE,
      extension = "tsv")
    strata <- extract_individuals_metadata(gds = data, whitelist = TRUE)
    readr::write_tsv(x = strata, file = strata.filename$filename)

    if (!is.null(strata)) {
      if (rlang::has_name(tidy.data, "TARGET_ID")) {
        tidy.data %<>%
          dplyr::left_join(strata, by = "TARGET_ID") %>%
          dplyr::select(-TARGET_ID) %>%
          dplyr::rename(POP_ID = STRATA)
      } else {
        tidy.data %<>%
          dplyr::left_join(strata, by = "INDIVIDUALS") %>%
          dplyr::rename(POP_ID = STRATA)
      }
    } else {
      tidy.data %<>% dplyr::mutate(POP_ID = 1L)
    }


    # Erase genotypes
    if (!is.null(missing.memory)) {
      message("Using missing.memory file to erase genotypes")
      missing <- radiator::read_rad(data = missing.memory) %>%
        dplyr::arrange(MARKERS, INDIVIDUALS)

      data %<>% dplyr::arrange(MARKERS, INDIVIDUALS)

      #check identical markers
      same.markers <- identical(unique(missing$MARKERS), unique(data$MARKERS))
      same.individuals <- identical(unique(missing$INDIVIDUALS), unique(data$INDIVIDUALS))
      if (!same.markers || !same.individuals) {
        message("note: data and missing memory don't share all the same markers and/or individuals")
        data %<>% dplyr::left_join(missing, by = c("MARKERS", "INDIVIDUALS"))
        #%>% dplyr::mutate(ERASE = replace(ERASE, which(is.na(ERASE)), FALSE))
        data$ERASE[is.na(data$ERASE)] <- FALSE # faster
        which.missing <- which(data$ERASE)
        data %<>% dplyr::select(-ERASE)
      } else {
        which.missing <- which(missing$ERASE)
      }

      message("Erasing genotypes and genotypes metadata...")
      if (rlang::has_name(data, "GT_BIN")) data$GT_BIN[which.missing] <- NA
      if (rlang::has_name(data, "GT")) data$GT[which.missing] <- "000000"
      if (rlang::has_name(data, "GT_VCF")) data$GT_VCF[which.missing] <- "./."
      if (rlang::has_name(data, "GT_VCF_NUC")) data$GT_VCF_NUC[which.missing] <- "./."
      if (rlang::has_name(data, "READ_DEPTH")) data$READ_DEPTH[which.missing] <- NA
      if (rlang::has_name(data, "ALLELE_REF_DEPTH")) data$ALLELE_REF_DEPTH[which.missing] <- NA
      if (rlang::has_name(data, "ALLELE_ALT_DEPTH")) data$ALLELE_ALT_DEPTH[which.missing] <- NA
    }#End missing.memory


    # write tidy
    write_rad(data = tidy.data, path = filename$filename)

    # Whitelist
    # readr::write_tsv(x = markers.meta, file.path(wf, "whitelist.markers.tidy.dart.tsv"))
    write_rad(data = markers.meta,
              path = wf,
              filename = "whitelist.markers.tidy.dart.tsv",
              tsv = TRUE,
              write.message = "standard",
              verbose = verbose)
    data <- tidy.data
  }#tidy

  # Results --------------------------------------------------------------------
  if (verbose) cat("################################### SUMMARY ####################################\n")
  message("\nNumber of chrom: ", n.chrom)
  message("Number of locus: ", n.locus)
  message("Number of SNPs: ", n.snp)
  summary_strata(strata)
  return(data)
}#End read_dart

# import_dart-------------------------------------------------------------------
#' @title import_dart
#' @description Read DArT file
#' @rdname import_dart
#' @keywords internal
#' @export
import_dart <- function(
  data,
  strata,
  pop.levels = NULL,
  path.folder = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE
) {

  # # TEST
  # pop.levels = NULL
  # path.folder = NULL
  # parallel.core = parallel::detectCores() - 1
  # verbose = TRUE

  message("Reading DArT file...")

  # Check that DArT file as good target id written -----------------------------
  target.id <- extract_dart_target_id(data, write = FALSE)
  n.ind.dart <- nrow(target.id)
  if (verbose) message("    Number of individuals: ", n.ind.dart)

  # Check DArT format file -----------------------------------------------------
  dart.check <- check_dart(data)
  silico.dart <- "silico.dart" %in% dart.check$data.type

  # Strata file ------------------------------------------------------------------
  strata.df <- radiator::read_strata(
    strata = strata,
    pop.levels = pop.levels,
    pop.labels = NULL,
    pop.select = NULL,
    blacklist.id = NULL,
    keep.two = FALSE,
    verbose = verbose
  )
  pop.levels <- strata.df$pop.levels
  pop.labels <- strata.df$pop.labels
  pop.select <- strata.df$pop.select
  blacklist.id <- strata.df$blacklist.id
  strata.df <- strata.df$strata
  n.ind.strata <- nrow(strata.df)

  strata.df %<>%
    dplyr::mutate(
      ID_SEQ = seq_len(length.out = dplyr::n()),
      STRATA_SEQ = as.integer(factor(x = STRATA, levels = levels(STRATA)))
    )

  # Check that TARGET_ID in strata match TARGET_ID in the DArT file ------------
  if (verbose) message("Using individuals in strata file to filter individuals in DArT file")
  if (n.ind.dart != n.ind.strata) {
    strata.id.check <- strata.df %>%
      dplyr::mutate(IN_DART = strata.df$TARGET_ID %in% target.id$TARGET_ID)
    strata.id.pass <- !FALSE %in% (unique(strata.id.check$IN_DART))
    if (!strata.id.pass) {
      problem.filename <- generate_filename(
        name.shortcut = "radiator_tidy_dart_strata_problem",
        path.folder = path.folder,
        date = TRUE,
        extension = "tsv")

      readr::write_tsv(
        x = strata.id.check,
        file = problem.filename$filename)
      rlang::abort("\nSome of the samples in the strata are not found in the DArT file.
                     For more info: ", problem.filename$filename.short)
    }
    if (verbose) message("\nNote: Careful if using DArT statistics generated for all samples...")
    strata.id.check <- NULL
  } else {
    if (!identical(sort(target.id$TARGET_ID), sort(strata.df$TARGET_ID))) {
      rlang::abort("\nThe DArT and strata files don't have the same TARGET_IDs")
    }
  }
  blacklist.id <- c(
    dplyr::filter(strata.df, !TARGET_ID %in% target.id$TARGET_ID) %$% TARGET_ID,
    dplyr::filter(target.id, !TARGET_ID %in% strata.df$TARGET_ID) %$% TARGET_ID
  ) %>% unique
  message("Number of blacklisted samples: ", length(blacklist.id))
  if (length(blacklist.id) == 0) blacklist.id <- NULL
  target.id <- NULL

  # need to check for duplicate names... yes happening all the time
  duplicate.id.strata <- length(strata.df$INDIVIDUALS) - length(unique(strata.df$INDIVIDUALS))

  if (duplicate.id.strata > 0) {
    message("Duplicated individuals names found in the strata.\n   number of duplicate names = ", duplicate.id.strata, "\n")
    rlang::abort("\nFix the strata with unique names and\nverify the DArT file for the same issue, adjust accordingly...")
  }


  # There's no big difference really here if we import everything and filter after...
  data <- suppressWarnings(
    data.table::fread(
      file = data,
      header = TRUE,
      strip.white = TRUE,
      na.strings = c("-", "NA"),
      stringsAsFactors = FALSE,
      skip = "CallRate",
      # drop = blacklist.id,
      # cannot use this because of upper/lower case and all the shity way
      # people name their samples...
      select = NULL,
      showProgress = TRUE,
      nThread = parallel.core,
      verbose = FALSE) %>%
      tibble::as_tibble(.)
  )

  # We want snakecase not camelcase
  # Change the TARGET_ID by INDIVIDUALS...
  data %<>%
    clean_dart_colnames(
      data = .,
      blacklist.id = blacklist.id,
      dart.col.num = dart.check$star.number,
      strata = strata.df
    )

  # keep consensus sequence if found -------------------------------------------
  # or rename TRIMMED_SEQUENCE

  # turns out I think it's better to keep ALLELE_SEQUENCE...
  # Changed 2020-02-06
  if (rlang::has_name(data, "ALLELE_SEQUENCE")) {
    data %<>% dplyr::rename(SEQUENCE = ALLELE_SEQUENCE)
  } else if (rlang::has_name(data, "TRIMMED_SEQUENCE")) {
    data %<>% dplyr::rename(SEQUENCE = TRIMMED_SEQUENCE)
  } else {
    if (rlang::has_name(data, "CLUSTER_CONSENSUS_SEQUENCE")) {
      data %<>% dplyr::rename(SEQUENCE = CLUSTER_CONSENSUS_SEQUENCE)
    }
  }


  # how many columns with sequence in it
  # test <- dplyr::select(.data = data, dplyr::ends_with("_SEQUENCE"))

  # Pre 2020-02-06
  # if (rlang::has_name(data, "CLUSTER_CONSENSUS_SEQUENCE")) {
  #   data %<>% dplyr::rename(SEQUENCE = CLUSTER_CONSENSUS_SEQUENCE)
  # } else if (rlang::has_name(data, "TRIMMED_SEQUENCE")) {
  #   data %<>% dplyr::rename(SEQUENCE = TRIMMED_SEQUENCE)
  # } else {
  #   if (rlang::has_name(data, "ALLELE_SEQUENCE")) {
  #     data %<>% dplyr::rename(SEQUENCE = ALLELE_SEQUENCE)
  #   }
  # }

  if (!silico.dart) {
    colnames(data) %<>%
      stringi::stri_replace_all_fixed(
        str = .,
        pattern = c("ALLELE_ID","SNP_POSITION"),
        replacement = c("LOCUS", "POS"),
        vectorize_all = FALSE)

    if (rlang::has_name(data, "CLONE_ID")) {
      if (rlang::has_name(data, "LOCUS")) {
        data %<>% dplyr::select(-CLONE_ID)
      } else {
        colnames(data) %<>%
          stringi::stri_replace_all_fixed(
            str = .,
            pattern = "CLONE_ID",
            replacement = "LOCUS",
            vectorize_all = FALSE)
      }
    }

    if (!rlang::has_name(data, "LOCUS")) {
      rlang::abort("\nProblem tidying DArT dataset: contact author")
    }

    # necessary steps...observed with DArT file using ref genome ---------------
    data %<>% dplyr::filter(!is.na(LOCUS))
    # Check for duplicate rows (sometimes people combine DArT data...)
    if (rlang::has_name(data, "POS")) {
      data %<>% dplyr::arrange(LOCUS, POS)
      data.dup <- nrow(dplyr::distinct(data, LOCUS, SNP, POS, CALL_RATE, .keep_all = FALSE))
    } else {
      data %<>% dplyr::arrange(LOCUS)
      data.dup <- nrow(dplyr::distinct(data, LOCUS, SNP, CALL_RATE, .keep_all = FALSE))
    }

    # make sure no duplicates
    if (nrow(data) != data.dup) {
      message("Duplicate rows were identified")
      message("    using distinct rows")
      message("    check data if downstream problems")
      data %<>% dplyr::distinct(LOCUS, SNP, POS, CALL_RATE, .keep_all = TRUE)
    }
    data.dup <- NULL

    # Screen for duplicate names -------------------------------------------------
    id <- purrr::discard(
      .x = colnames(data),
      .p = colnames(data) %in% c("LOCUS", "SNP", "POS", "CALL_RATE", "AVG_COUNT_REF",
                                 "AVG_COUNT_SNP", "REP_AVG")
    )
    dup.id <- length(id) - length(unique(id))
    if (dup.id > 0) {
      rlang::abort(stringi::stri_join("Duplicated individual names in the data: ", dup.id))
    }
    # removing unused object
    remove.list <- id <- dup.id <- NULL

    # clean locus, generate MARKERS and VARIANT_ID
    data %<>% clean_dart_locus(.)
  }

  # DArT characteristics--------------------------------------------------------
  if (verbose)  message("\nDArT characteristics:")
  # dart.format:
  # "1row" !binary
  # "2rows" binary
  # "counts" binary
  # "silico.dart"
  dart.format <- detect_dart_format(
    x = data,
    target.id = strata.df$INDIVIDUALS,
    verbose = TRUE)

  return(res = list(data = data, strata = strata.df, dart.format = dart.format))
}#End import_dart

# clean_dart_locus--------------------------------------------------------------
#' @title clean_dart_locus
#' @description Clean LOCUS and generate VARIANT_ID and MARKERS
#' @rdname clean_dart_locus
#' @keywords internal
#' @export
clean_dart_locus <- function(x, fast = TRUE) {

  if (fast) {
    # x <- data
    if (!rlang::has_name(x, "CHROM")) x %<>% dplyr::mutate(CHROM = "CHROM_1")
    want <- c("VARIANT_ID", "M_SEQ", "MARKERS", "CHROM", "LOCUS", "POS", "COL",
              "REF", "ALT", "SEQUENCE",
              "CALL_RATE", "AVG_COUNT_REF", "AVG_COUNT_SNP", "REP_AVG")
    x %<>%
      dplyr::mutate(
        COL = stringi::stri_extract_first_regex(
          str = LOCUS,
          pattern = "[-][0-9]+[\\:]"),
        COL = stringi::stri_replace_all_fixed(
          str = COL,
          pattern = c("-", ":"),
          replacement = c("", ""),
          vectorize_all = FALSE),
        COL = as.integer(COL),
        LOCUS = stringi::stri_extract_first_regex(str = LOCUS, pattern = "^[0-9]+"),
        REF = stringi::stri_extract_first_regex(str = SNP, pattern = "[A-Z]"),
        ALT = stringi::stri_extract_last_regex(str = SNP, pattern = "[A-Z]"),
        MARKERS = stringi::stri_join(CHROM, LOCUS, POS, sep = "__"),
        VARIANT_ID = as.integer(factor(MARKERS)),
        M_SEQ = VARIANT_ID,
        SNP = NULL
      ) %>%
      dplyr::select(tidyselect::any_of(want), tidyselect::everything()) %>%
      dplyr::mutate(
        dplyr::across(
          .cols = c(MARKERS, CHROM, LOCUS, POS),
          .fns = as.character
        )
      ) %>%
      dplyr::arrange(CHROM, LOCUS, POS, REF)

  } else {
    want <- c("VARIANT_ID", "M_SEQ", "MARKERS", "CHROM", "LOCUS", "POS", "COL",
              "REF", "ALT", "SEQUENCE",
              "CALL_RATE", "AVG_COUNT_REF", "AVG_COUNT_SNP", "REP_AVG")
    x %<>%
      tidyr::separate(col = LOCUS,
                      into = c("LOCUS", NA),
                      sep = "\\|",
                      extra = "drop"
      ) %>%
      tidyr::separate(col = SNP,
                      into = c(NA, "KEEPER"),
                      sep = ":",
                      extra = "drop") %>%
      tidyr::separate(col = KEEPER, into = c("REF", "ALT"), sep = ">") %>%
      dplyr::mutate(
        CHROM = rep("CHROM_1", n()),
        MARKERS = stringi::stri_join(CHROM, LOCUS, POS, sep = "__"),
        VARIANT_ID = as.integer(factor(MARKERS)),
        M_SEQ = VARIANT_ID,
        COL = POS
      ) %>%
      dplyr::select(tidyselect::any_of(want), tidyselect::everything()) %>%
      dplyr::mutate(
        dplyr::across(
          .cols = c(MARKERS, CHROM, LOCUS, POS),
          .fns = as.character
        )
      ) %>%
      dplyr::arrange(CHROM, LOCUS, POS, REF)
  }
}#End clean_dart_locus

# detect_dart_format-------------------------------------------------------------
#' @title detect_dart_format
#' @description Detect the dart genotype format: 1row, 2rows or counts
#' @rdname detect_dart_format
#' @keywords internal
#' @export
detect_dart_format <- function(x = NULL, target.id = NULL, verbose = TRUE) {
  # dart.format:
  # silico.dart
  # 1row
  # 2rows
  # counts

  if (rlang::has_name(x, "CLONE_ID") && rlang::has_name(x, "SEQUENCE")) {
    if (verbose) message("DArT SNP format: silico DArT")
    dart.format <- "silico.dart"
  } else {

    binary <- anyNA(x$REF)

    if (!binary) {
      if (verbose) message("DArT SNP format: genotypes in 1 Row")
      dart.format <- "1row"
    } else {
      count.data <- x %>%
        dplyr::select(
          tidyselect::any_of(
            sample(x = target.id, size = min(10, floor(0.1 * length(target.id))))
          )
        ) %>%
        dplyr::mutate(dplyr::across(tidyselect::everything(), .fns = as.numeric)) %>%
        purrr::flatten_dbl(.) %>%
        unique(.)

      count.data <- any(count.data > 1, na.rm = TRUE)

      if (count.data) {
        dart.format <- "counts"
        if (verbose) message("DArT SNP format: alleles coverage in 2 Rows counts")
      } else {
        dart.format <- "2rows"
        if (verbose) message("DArT SNP format: alleles absence/presence in 2 Rows")
      }
    }
  }
  return(dart.format)
}#End detect_dart_format

# switch_allele_count----------------------------------------------------------------
#' @title switch_allele_count
#' @description Switch from dart genotyping 2 presence absence
#' @rdname switch_allele_count
#' @keywords internal
#' @export

switch_allele_count <- function(x, dart.group = FALSE, ref = TRUE) {
  # dart.group == TRUE (=1rows) vs presence/absence coding
  # dart.group == FALSE 2rows: presence/absence

  if (dart.group) {
    # 1 = HOM ALT
    # 2 = HET
    # 0 = HOM REF
    if (ref) {#REF count
      x <- dplyr::case_when(
        x == 1L ~ 0L,
        x == 2L ~ 1L,
        x == 0L ~ 1L
      )
    } else {#ALT count
      x <- dplyr::if_else(x == 2L, 1L, x)
      # the other coding remain the same
    }
  } else {
    # 2rows: presence/absence
    # here we want count of alternate allele instead...
    # x <- as.integer(dplyr::recode(.x = as.character(x), "0" = "1", "1" = "0"))
    # case_when is much faster than recode...
    x <- dplyr::case_when(
      x == 0L ~ 1L,
      x == 1L ~ 0L
    )
  }
  return(x)
}# End switch_allele_count


# dart2gds----------------------------------------------------------------
#' @title dart2gds
#' @description Transform dart to GDS format
#' @rdname dart2gds
#' @keywords internal
#' @export

dart2gds <- function(
  data,
  strata = NULL,
  markers.meta,
  filename,
  dart.format,
  gt.vcf = FALSE,
  gt.vcf.nuc = FALSE,
  gt = FALSE,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE
) {

  # Note to myself:
  # It looks intimidating but you've given a lot of thought about it
  # It's one of the fastest way you've found...


  # data.bk <- data
  if (dart.format == "1row") data.source <- c("dart", "1row")
  if (dart.format == "2rows") data.source <- c("dart", "2rows")
  if (dart.format == "counts") data.source <- c("dart", "counts")
  parallel.core.temp <- max(1, 2 * floor(parallel.core / 2))
  geno.coding <- "alt.dos"
  dp <- gl <- NULL

  want <- c("VARIANT_ID", "MARKERS", "REF", "ALT", strata$INDIVIDUALS)
  data %<>%
    dplyr::select(tidyselect::any_of(want)) %>%
    dplyr::arrange(MARKERS, REF)

  if (gt.vcf.nuc || gt) {
    ref.alt <- TRUE
    want <- c("VARIANT_ID", "MARKERS","REF", "ALT")
  } else {
    ref.alt <- FALSE
    want <- c("VARIANT_ID", "MARKERS")
  }

  message("Generating genotypes and calibrating REF/ALT alleles...")
  #2-rows or counts
  if (TRUE %in% (c("counts", "2rows") %in% data.source)) {
    n.ind <- ncol(data) - 4
    n.snp <- nrow(data) / 2
    # alt <- dplyr::filter(data, !is.na(REF)) %>% dplyr::select(-REF, -ALT, -MARKERS)
    # ref <- dplyr::filter(data, is.na(REF)) %>% dplyr::select(-REF, -ALT, -MARKERS)
    alt <- dplyr::filter(data, !is.na(REF)) %>% dplyr::select(-REF, -ALT)
    ref <- dplyr::filter(data, is.na(REF)) %>% dplyr::select(-REF, -ALT)
  }#2rows genotypes

  #1-row
  if ("1row" %in% data.source) {
    n.ind <- ncol(data) - 4
    n.snp <- nrow(data)

    # generate allele count from the 1 row dart-----------------------------------
    alt <- dplyr::select(data, -REF, -ALT) %>%
      dplyr::mutate(
        dplyr::across(
          .cols = strata$INDIVIDUALS,
          .fns = switch_allele_count,
          dart.group = TRUE,
          ref = FALSE
        )
      )

    ref <- dplyr::select(data, -REF, -ALT) %>%
      dplyr::mutate(
        dplyr::across(
          .cols = strata$INDIVIDUALS,
          .fns = switch_allele_count,
          dart.group = TRUE,
          ref = TRUE
        )
      )
  }#1row genotypes

  # Common approach between all---------- -------------------------------------
  switch <- ref$VARIANT_ID[
    rowSums(x = dplyr::select(ref, -VARIANT_ID, -MARKERS), na.rm = TRUE) <
      rowSums(x = dplyr::select(alt, -VARIANT_ID, -MARKERS), na.rm = TRUE)
  ]
  n.switch <- length(switch)
  if (n.switch > 0) {
    if ("counts" %in% data.source) {
      count.what <- " read depth"
    } else {
      count.what <- ""
    }
    message("Number of markers recalibrated based on counts of allele", count.what,": ", n.switch)
    alt.s <- dplyr::filter(alt, VARIANT_ID %in% switch)
    alt <- dplyr::filter(alt, !VARIANT_ID %in% switch)
    ref.s <- dplyr::filter(ref, VARIANT_ID %in% switch)
    ref <- dplyr::filter(ref, !VARIANT_ID %in% switch)

    alt %<>% dplyr::bind_rows(ref.s) %>%
      dplyr::arrange(VARIANT_ID)
    ref %<>% dplyr::bind_rows(alt.s) %>%
      dplyr::arrange(VARIANT_ID)
    alt.s <- ref.s <- NULL

    # re-calibrate the markers.meta
    markers.meta.s <- dplyr::filter(markers.meta, VARIANT_ID %in% switch)

    switch.ref <- dplyr::select(markers.meta.s, dplyr::contains("REF")) %>%
      dplyr::select(-dplyr::contains("REF_AVG"))
    colnames(switch.ref) <- stringi::stri_replace_all_regex(
      str = colnames(switch.ref), pattern = c("_REF", "^REF"),
      replacement = c("_SNP", "ALT"), vectorize_all = FALSE
    )
    switch.alt <- dplyr::select(markers.meta.s, dplyr::contains("snp")) %>%
      dplyr::select(-dplyr::contains("REF_AVG"))
    colnames(switch.alt) <- stringi::stri_replace_all_fixed(
      str = colnames(switch.alt), pattern = "SNP",
      replacement = "REF", vectorize_all = FALSE
    )
    markers.meta.s %<>%
      dplyr::select(-tidyselect::any_of(unique(c(colnames(switch.ref), colnames(switch.alt))))) %>%
      dplyr::bind_cols(switch.ref) %>%
      dplyr::bind_cols(switch.alt)
    switch.alt <- switch.ref <- NULL

    markers.meta %<>%
      dplyr::filter(!VARIANT_ID %in% switch) %>%
      dplyr::bind_rows(markers.meta.s) %>%
      dplyr::arrange(VARIANT_ID)
    markers.meta.s <- NULL
  }#End if switching alleles

  if (ref.alt) {
    alt %<>% dplyr::left_join(
      dplyr::select(markers.meta, VARIANT_ID, REF, ALT),
      by = "VARIANT_ID")
  }

  genotypes.meta <- suppressMessages(
    radiator::rad_long(
      x = alt,
      cols = want,
      names_to = "INDIVIDUALS",
      values_to = "ALLELE_ALT_DEPTH",
      variable_factor = FALSE
    ) %>%
      dplyr::bind_cols(
        radiator::rad_long(
          x = ref,
          cols = c("VARIANT_ID", "MARKERS"),
          names_to = "INDIVIDUALS",
          values_to = "ALLELE_REF_DEPTH",
          variable_factor = FALSE
        )
      )
  )
  ref <- alt <- NULL

  # Faster to check that the bind_cols worked by checking the variant id
  if (!identical(genotypes.meta$VARIANT_ID...1, genotypes.meta$VARIANT_ID...5)) {
    rlang::abort("Contact author, DArT tiding problem")
  } else {
    genotypes.meta %<>% dplyr::select(-VARIANT_ID...5) %>% dplyr::rename(VARIANT_ID = VARIANT_ID...1)
  }

  # Faster to check that the bind_cols worked by checking markers
  if (!identical(genotypes.meta$MARKERS...2, genotypes.meta$MARKERS...6)) {
    rlang::abort("Contact author, DArT tiding problem")
  } else {
    genotypes.meta %<>% dplyr::select(-MARKERS...6) %>% dplyr::rename(MARKERS = MARKERS...2)
  }

  if (!identical(genotypes.meta$INDIVIDUALS...3, genotypes.meta$INDIVIDUALS...7)) {
    rlang::abort("Contact author, DArT tiding problem")
  } else {
    genotypes.meta %<>% dplyr::select(-INDIVIDUALS...7) %>% dplyr::rename(INDIVIDUALS = INDIVIDUALS...3)
  }

  # M_SEQ and ID_SEQ
  strata.temp <- dplyr::select(.data = strata, tidyselect::any_of(c("INDIVIDUALS", "ID_SEQ")))
  genotypes.meta %<>%
    dplyr::mutate(M_SEQ = VARIANT_ID) %>%
    dplyr::left_join(strata.temp, by = "INDIVIDUALS")


  # modify genotypes meta to generate GT_BIN
  if ("counts" %in% data.source) {
    genotypes.meta %<>%
      dplyr::mutate(
        READ_DEPTH = ALLELE_REF_DEPTH + ALLELE_ALT_DEPTH,
        GT_BIN = dplyr::case_when(
          ALLELE_REF_DEPTH > 0L & ALLELE_ALT_DEPTH == 0L ~ 0L,
          ALLELE_REF_DEPTH > 0L & ALLELE_ALT_DEPTH > 0L ~ 1L,
          ALLELE_REF_DEPTH == 0L & ALLELE_ALT_DEPTH > 0L ~ 2L
        )
      ) %>%
      dplyr::mutate(
        dplyr::across(
          .cols = c(READ_DEPTH, ALLELE_REF_DEPTH, ALLELE_ALT_DEPTH),
          .fns = replace_by_na, what = 0L
        )
      )
  } else {
    genotypes.meta %<>%
      dplyr::mutate(
        dplyr::across(
          .cols = "ALLELE_REF_DEPTH",
          .fns = switch_allele_count, dart.group = FALSE
        )
      ) %>%
      dplyr::mutate(
        GT_BIN = ALLELE_REF_DEPTH + ALLELE_ALT_DEPTH,
        ALLELE_ALT_DEPTH = NULL,
        ALLELE_REF_DEPTH = NULL
      )
  }



  # genotypes array ----------------------------------------------------------
  genotypes <- gt2array(
    gt.bin = genotypes.meta$GT_BIN, n.ind = n.ind, n.snp = n.snp
  )

  if (gt.vcf) {
    genotypes.meta %<>%
      dplyr::mutate(
        GT_VCF = dplyr::case_when(
          GT_BIN == 0 ~ "0/0", GT_BIN == 1 ~ "0/1", GT_BIN == 2 ~ "1/1",
          is.na(GT_BIN) ~ "./.")
      )
  }
  if (gt.vcf.nuc) {
    genotypes.meta %<>%
      dplyr::mutate(
        GT_VCF_NUC = dplyr::case_when(
          GT_BIN == 0 ~ stringi::stri_join(REF, REF, sep = "/"),
          GT_BIN == 1 ~ stringi::stri_join(REF, ALT, sep = "/"),
          GT_BIN == 2 ~ stringi::stri_join(ALT, ALT, sep = "/"),
          is.na(GT_BIN) ~ "./.")
      )
  }
  if (gt) {
    genotypes.meta %<>%
      dplyr::mutate(
        GT = stringi::stri_replace_all_fixed(
          str = GT_VCF_NUC,
          pattern = c("A", "C", "G", "T", "/", ".."),
          replacement = c("001", "002", "003", "004", "", "000000"),
          vectorize_all = FALSE)
      )
  }

  # Generate GDS ---------------------------------------------------------------
  ## test
  # biallelic = TRUE
  # ad = NULL
  # gl = NULL
  # open = TRUE
  data <- radiator_gds(
    strata = strata,
    genotypes = genotypes,
    markers.meta = markers.meta,
    genotypes.meta = genotypes.meta,
    biallelic = TRUE,
    data.source = data.source,
    dp = dp,
    ad = NULL,
    geno.coding = geno.coding,
    filename = filename,
    open = TRUE,
    verbose = verbose
  )

  if (verbose) message("done!")
  return(data)
}# End dart2gds



# tidy_dart_metadata--------------------------------------------------------------
#' @name tidy_dart_metadata

#' @title Import and tidy \href{http://www.diversityarrays.com}{DArT} metadata.

#' @description Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.
#' The function generate a tidy dataset of
#' \href{http://www.diversityarrays.com}{DArT} markers and associated metadata.
#' Usefull to filter before importing the actual dataset.

#' @param data DArT output file. Note that most popular formats used by DArT are
#' recognised (1- and 2- row format, also called binary, and count data.).
#' If you encounter a problem, sent me your data so that I can update
#' the function. The function can import \code{.csv} or \code{.tsv} files.


#' @inheritParams tidy_genomic_data
#' @inheritParams radiator_common_arguments

#' @return A tidy dataframe with these columns:
#' \enumerate{
#' \item MARKERS: generated by radiator and correspond to CHROM + LOCUS + POS
#' separated by 2 underscores.
#' \item CHROM: the chromosome, for de novo: CHROM_1.
#' \item LOCUS: the locus.
#' \item POS: the SNP id on the LOCUS.
#' \item REF: the reference allele.
#' \item ALT: the alternate allele.
#' \item CALL_RATE: call rate output specific of DArT.
#' \item AVG_COUNT_REF: the coverage for the reference allele, output specific of DArT.
#' \item AVG_COUNT_SNP: the coverage for the alternate allele, output specific of DArT.
#' \item REP_AVG: the reproducibility average, output specific of DArT.
#' }

#' @export
#' @rdname tidy_dart_metadata

#' @examples
#' \dontrun{
#' clownfish.dart.tidy <- radiator::tidy_dart_metadata(
#' data = "clownfish.dart.tsv",
#' verbose = TRUE)
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

tidy_dart_metadata <- function(
  data,
  filename = NULL,
  verbose = FALSE,
  parallel.core = parallel::detectCores() - 1
) {
  # Cleanup-------------------------------------------------------------------
  # obj.keeper <- c(ls(envir = globalenv()), "input")
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (verbose) message("Execution date@time: ", file.date)
  old.dir <- getwd()
  opt.change <- getOption("width")
  options(width = 70)
  timing <- radiator_tic()
  # res <- list()
  #back to the original directory and options
  on.exit(setwd(old.dir), add = TRUE)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(radiator_toc(timing), add = TRUE)
  on.exit(radiator_function_header(f.name = "", start = FALSE, verbose = verbose), add = TRUE)
  # on.exit(rm(list = setdiff(ls(envir = sys.frame(-1L)), obj.keeper), envir = sys.frame(-1L)))


  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("Input file missing")
  # Date and Time --------------------------------------------------------------
  if (!is.null(filename)) {
    meta.filename <- stringi::stri_join(filename, "_metadata_", file.date,".rad")
  }

  data.type <- radiator::detect_genomic_format(data)

  if (!data.type %in% c("dart", "fst.file")) {
    rlang::abort("Contact author to show your DArT data, problem duting import")
  }

  if (verbose) message("Importing DArT markers metadata")

  # Import metadata-------------------------------------------------------------
  if (data.type == "dart") {
    dart.delim <- "," # for csv files

    if (stringi::stri_detect_fixed(
      str = stringi::stri_sub(str = data, from = -4, to = -1),
      pattern = ".tsv")) {
      dart.delim <- "\t" # for tsv files
    }

    dart.check <- check_dart(data)
    if (!dart.check$data.type %in% c("dart", "silico.dart")) {
      rlang::abort("Contact author to show your DArT data, problem during import")
    } else {
      skip.number <- dart.check$skip.number
    }

    dart.col.type <- readr::read_delim(
      file = data,
      delim = dart.delim,
      skip = skip.number,
      n_max = 1,
      na = "-",
      col_names = FALSE,
      col_types = readr::cols(.default = readr::col_character()))


    want <- tibble::tibble(
      INFO = c("ALLELEID", "SNP", "SNPPOSITION", "CALLRATE",
               "AVGCOUNTREF", "AVGCOUNTSNP", "REPAVG"),
      COL_TYPE = c("c", "c", "i", "d", "d", "d", "d"))

    dart.col.type %<>%
      radiator::rad_long(
        x = .,
        cols = tidyselect::everything(),
        names_to = "DELETE",
        values_to = "INFO"
      ) %>%
      dplyr::select(-DELETE) %>%
      dplyr::mutate(INFO = stringi::stri_trans_toupper(INFO)) %>%
      dplyr::left_join(want, by = "INFO") %>%
      dplyr::mutate(COL_TYPE = stringi::stri_replace_na(str = COL_TYPE, replacement = "_")) %>%
      dplyr::select(COL_TYPE) %>%
      purrr::flatten_chr(.) %>% stringi::stri_join(collapse = "")

    input <- suppressMessages(
      suppressWarnings(
        readr::read_delim(
          file = data,
          delim = dart.delim,
          skip = skip.number,
          na = "-",
          col_names = TRUE,
          col_types = dart.col.type)
      )
    )
    colnames(input) <- stringi::stri_trans_toupper(colnames(input))
    colnames(input) <- stringi::stri_replace_all_fixed(
      str = colnames(input),
      pattern = c("AVGCOUNTREF", "AVGCOUNTSNP", "REPAVG", "ALLELEID", "SNPPOSITION", "CALLRATE"),
      replacement = c("AVG_COUNT_REF", "AVG_COUNT_SNP", "REP_AVG", "LOCUS", "POS", "CALL_RATE"),
      vectorize_all = FALSE)

    # Check for duplicate rows (sometimes people combine DArT data...)----------
    input.dup <- dplyr::distinct(input, LOCUS, SNP, POS, CALL_RATE, .keep_all = FALSE)

    # make sure no duplicates
    if (nrow(input) != nrow(input.dup)) {
      input.dup <- NULL
      message("Duplicate rows were identified")
      message("    using distinct rows")
      message("    check input data if downstream problems")
      input <- dplyr::distinct(input, LOCUS, SNP, POS, CALL_RATE, .keep_all = TRUE)
    }
    input.dup <- NULL

    # Tidying data ---------------------------------------------------------------
    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "REF", "ALT", "CALL_RATE",
              "AVG_COUNT_REF", "AVG_COUNT_SNP", "REP_AVG")

    input %<>%
        tidyr::separate(col = LOCUS, into = c("LOCUS", "NOT_USEFUL"), sep = "\\|", extra = "drop") %>%
        dplyr::select(-NOT_USEFUL) %>%
        tidyr::separate(col = SNP, into = c("NOT_USEFUL", "KEEPER"), sep = ":", extra = "drop") %>%
        dplyr::select(-NOT_USEFUL) %>%
        tidyr::separate(col = KEEPER, into = c("REF", "ALT"), sep = ">") %>%
        dplyr::mutate(
          CHROM = rep("CHROM_1", n()),
          MARKERS = stringi::stri_join(CHROM, LOCUS, POS, sep = "__")) %>%
        dplyr::select(tidyselect::any_of(want), tidyselect::everything()) %>%
        dplyr::mutate(
          dplyr::across(
            .cols = c(MARKERS, CHROM, LOCUS, POS),
            .fns = as.character
          )
        ) %>%
        dplyr::arrange(CHROM, LOCUS, POS, REF) %>%
        dplyr::filter(!is.na(REF) | !is.na(ALT)) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
        dplyr::arrange(MARKERS)
  }

  if (data.type == "fst.file") {
    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "REF", "ALT", "CALL_RATE",
              "AVG_COUNT_REF", "AVG_COUNT_SNP", "REP_AVG")
    input <- radiator::read_rad(data) %>%
      dplyr::select(tidyselect::any_of(want)) %>%
      dplyr::distinct(MARKERS, .keep_all = TRUE)
  }

  if (!is.null(filename)) {
    radiator::write_rad(data = input, path = meta.filename)
    short.name <- list.files(path = ".", pattern = "metadata")
    if (length(short.name) > 1) {
      short.name <- file.info(short.name) %>%
        tibble::rownames_to_column(df = ., var = "FILE") %>%
        dplyr::filter(mtime == max(mtime))
      short.name <- short.name$FILE
    }
    message("Marker's metadata file written:\n    ", short.name)
  }
  # Results --------------------------------------------------------------------
  if (verbose) {
    n.chrom <- length(unique(input$CHROM))
    n.locus <- length(unique(input$LOCUS))
    n.snp <- length(unique(input$MARKERS))
    message("\nNumber of chrom: ", n.chrom)
    message("Number of locus: ", n.locus)
    message("Number of SNPs: ", n.snp)
  }
  return(input)
}#End dart_markers_metadata


# clean_dart_colnames----------------------------------------------------------------
#' @title clean_dart_colnames
#' @description clean_dart_colnames (only the DArT columns = snakecase...)
#' @rdname clean_dart_colnames
#' @keywords internal
#' @export
clean_dart_colnames <- function(data, blacklist.id = NULL, dart.col.num = NULL, strata = NULL) {

  if (dart.col.num > 0) {
    keeper <- dart.col.num# the most reliable if available
  } else {
    if (is.null(blacklist.id)) {
      keeper <- length(colnames(data)) - length(strata$TARGET_ID)
    } else {
      keeper <- length(colnames(data)) - length(strata$TARGET_ID) - length(blacklist.id)
    }
  }

  # clean the dart header
  # clean the dart target ids
  # here doing all colnames is ok, because keeper above is allready upper caps
  # added 20190528
  colnames(data) <- c(
    radiator::radiator_snakecase(x = colnames(data)[1:keeper]),
    stringi::stri_trans_toupper(clean_ind_names(colnames(data)[-c(1:keeper)]))
  )

  if (!is.null(blacklist.id)) {
    data %<>% dplyr::select(-tidyselect::any_of(blacklist.id))
  }

  colnames(data) <- tibble::tibble(TARGET_ID = colnames(data)) %>%
    dplyr::left_join(strata, by = "TARGET_ID") %>%
    dplyr::mutate(
      INDIVIDUALS = dplyr::if_else(
        is.na(INDIVIDUALS), TARGET_ID, INDIVIDUALS)
    ) %$% INDIVIDUALS

  # Below generate errors when some id are very close... ID-10 and ID-1
  # colnames(data) <- stringi::stri_replace_all_fixed(
  #   str = colnames(data),
  #   pattern = strata$TARGET_ID,
  #   replacement = strata$INDIVIDUALS,
  #   vectorize_all = FALSE
  #   )
  return(data)
}#End clean_dart_colnames




## Merge dart-------------------------------------------------------------------

#' @title Merge DArT files
#' @description This function allows to merge 2 DArT files (filtered of not).

#' @param dart1 Full path of the first DArT file.
#' @param strata1 Full path of the first strata file for dart1.
#' @param dart2 Full path of the second DArT file.
#' @param strata2 Full path of the second strata file for dart2.

# deprecated
# @param keep.rad Unless the default is changed, the function removes the
# temporary tidy dart files generated for each DArT datasets during import.
# Default: \code{keep.rad = FALSE}.

#' @param filename Name of the merged DArT file.
#' By default, the function gives the merged data the filename:\code{merge_dart}
#' with date and time appended. The function will also append the date and time
#' to the filename provided.
#' The data is written in the working directory.
#' Default: \code{filename = NULL}.

#' @inheritParams tidy_genomic_data
#' @inheritParams read_dart
#' @inheritParams radiator_common_arguments
#' @inheritParams filter_common_markers
#' @inheritParams filter_monomorphic

#' @param remove.non.immortalized.dart.markers (logical). By default the function
#' will remove markers starting \code{1000}, those are called \code{non-immortalized}
#' markers by DArT.
#' Default: \code{remove.non.immortalized.dart.markers = TRUE}.

#' @details
#' The function average across markers the columns: CALL_RATE, REP_AVG,
#' AVG_COUNT_REF and AVG_COUNT_SNP, when found in the data.
#' For DArT, theses columns represent:
#' \itemize{
#' \item CALL_RATE: is the proportion of samples for which the genotype was called.
#' \item REP_AVG: is the proportion of technical replicate assay pairs for which
#' the marker score is consistent.
#' \item AVG_COUND_REF and AVG_COUND_SNP: the mean coverage for the reference and
#' alternate alleles, respectively.
#' }
#'
#' The function removes markers with starting with 1000 that are not immortalized by DArT
#'
#'
#' When the argument \strong{common.markers} is kept to TRUE, the function
#' produces an
#' \href{https://github.com/hms-dbmi/UpSetR}{UpSet plot} to visualize the number
#' of markers common or not between populations. The plot is not saved automatically,
#' this as to be done manually by the user.

#' @return The function returns a list in the global environment and 2 data frames
#' in the working directory. The dataframes are the tidy dataset and a strata file
#' of the 2 merged DArT files.

#' @examples
#' \dontrun{
#' # The simplest way to run the function:
#' merged.data <- radiator::merge_dart(
#' dart1 = "bluefin_larvae.tsv", strata1 = "strata1_bft_larvae.tsv",
#' dart2 = "bluefin_adults.csv", strata2 = "strata2_bft_adults.tsv")
#' }
#' @rdname merge_dart
#' @export
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

merge_dart <- function(
  dart1, strata1,
  dart2, strata2,
  filter.monomorphic = TRUE,
  filter.common.markers = TRUE,
  # keep.rad = FALSE,
  filename = NULL,
  remove.non.immortalized.dart.markers = TRUE,
  parallel.core = parallel::detectCores() - 1,
  ...
) {

  # test
  # keep.rad = TRUE
  # filename = NULL
  # filter.monomorphic = TRUE
  # filter.common.markers = TRUE
  # parallel.core = 12L
  # remove.non.immortalized.dart.markers = TRUE

  # Cleanup---------------------------------------------------------------------
  obj.keeper <- c(ls(envir = globalenv()), "res")
  radiator_function_header(f.name = "merge_dart", verbose = TRUE)
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  message("Execution date@time: ", file.date)
  old.dir <- getwd()
  opt.change <- getOption("width")
  options(future.globals.maxSize= Inf)
  options(width = 70)
  timing <- radiator_tic()
  #back to the original directory and options
  on.exit(setwd(old.dir), add = TRUE)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(radiator_toc(timing), add = TRUE)
  on.exit(rm(list = setdiff(ls(envir = sys.frame(-1L)), obj.keeper), envir = sys.frame(-1L)))
  on.exit(radiator_function_header(f.name = "merge_dart", start = FALSE, verbose = TRUE), add = TRUE)

  # manage missing arguments -----------------------------------------------------
  if (missing(dart1)) rlang::abort("dart1 file missing")
  if (missing(dart2)) rlang::abort("dart2 file missing")
  if (missing(strata1)) rlang::abort("strata1 file missing")
  if (missing(strata2)) rlang::abort("strata2 file missing")

  # Filename -------------------------------------------------------------------
  # Get date and time to have unique filenaming
  if (is.null(filename)) {
    filename <- stringi::stri_join("merged_dart_", file.date, ".rad")
    strata.filename <- stringi::stri_join("merged_dart_strata_", file.date, ".tsv")
  } else {
    filename <- stringi::stri_join(filename, "_", file.date, ".rad")
    strata.filename <- stringi::stri_join(filename, "_strata_", file.date, ".tsv")
  }

  # File type detection---------------------------------------------------------
  data.type <- radiator::detect_genomic_format(dart2)
  if (data.type != "dart") rlang::abort("DArT file not recognize by radiator")
  data.type <- radiator::detect_genomic_format(dart1)
  if (data.type != "dart") rlang::abort("DArT file not recognize by radiator")

  # import data ----------------------------------------------------------------
  # generate the folder to store the merged data
  path.folder <- generate_folder(
    f = NULL,
    rad.folder = "merge_dart",
    prefix_int = FALSE,
    internal = FALSE,
    file.date = file.date,
    verbose = TRUE)

  # the DArT data
  message("Importing and tidying dart1...")
  input <- suppressWarnings(
    radiator::read_dart(
      data = dart1,
      strata = strata1,
      filename = "dart1",
      tidy.dart = TRUE,
      parallel.core = parallel.core,
      path.folder = path.folder,
      verbose = FALSE
    )
  )


  message("Importing and tidying dart2...")
  dart2 <- suppressWarnings(
    radiator::read_dart(
      data = dart2,
      strata = strata2,
      filename = "dart2",
      tidy.dart = TRUE,
      parallel.core = parallel.core,
      path.folder = path.folder,
      verbose = FALSE
    )
  )

  # cleaning up non-immortalized markers ---------------------------------------
  markers.before <- length(unique(input$MARKERS)) + length(unique(dart2$MARKERS))
  if (remove.non.immortalized.dart.markers) {
    message("Removing non-immortalized DArT markers...")
    input <- suppressWarnings(dplyr::filter(input, !stringi::stri_detect_regex(str = LOCUS, pattern = "^1000")))
    dart2 <- suppressWarnings(dplyr::filter(dart2, !stringi::stri_detect_regex(str = LOCUS, pattern = "^1000")))
    markers.after <- length(unique(input$MARKERS)) + length(unique(dart2$MARKERS))
    message("    Number of blacklisted markers: ", markers.before - markers.after)
  }
  # merging DArT tidy data -----------------------------------------------------
  # we keep common column
  dart.col <- dplyr::intersect(colnames(input), colnames(dart2))

  input %<>%
    dplyr::select(tidyselect::any_of(dart.col)) %>%
    dplyr::bind_rows(dplyr::select(dart2, tidyselect::any_of(dart.col)))
  dart2 <- NULL



  # Remove weird markers ---------------------------------------------------------
  setwd(path.folder)
  biallelic.data <- NULL # global variables...
  input <- radiator::detect_biallelic_problems(
    data = input,
    verbose = TRUE,
    parallel.core = parallel.core
  ) %$%
    biallelic.data
  setwd(old.dir)

  if (rlang::has_name(input, "POLYMORPHIC")) {
    input %<>% dplyr::select(-POLYMORPHIC)
  }

  # Averaging across markers the call rate and other DArT markers metadata statistics
  if (rlang::has_name(input, "CALL_RATE")) {
    message("Averaging across markers the call rate")
    input %<>%
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(CALL_RATE = mean(CALL_RATE)) %>%
      dplyr::ungroup(.)
  }

  if (rlang::has_name(input, "REP_AVG")) {
    message("Averaging across markers the REP_AVG")
    input %<>%
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(REP_AVG = mean(REP_AVG)) %>%
      dplyr::ungroup(.)
  }

  if (rlang::has_name(input, "AVG_COUNT_REF")) {
    message("Averaging across markers the coverage for the REF allele")
    input %<>%
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(AVG_COUNT_REF = mean(AVG_COUNT_REF)) %>%
      dplyr::ungroup(.)
  }

  if (rlang::has_name(input, "AVG_COUNT_ALT")) {
    message("Averaging across markers the coverage for the ALT allele")
    input %<>%
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(AVG_COUNT_REF = mean(AVG_COUNT_REF)) %>%
      dplyr::ungroup(.)
  }

  # monomorphic markers (the TRUE/FALSE is included inside the function)
  input <- radiator::filter_monomorphic(
    data = input,
    filter.monomorphic = filter.monomorphic,
    path.folder = path.folder,
    verbose = TRUE)

  # common markers (the TRUE/FALSE is included inside the function)
  input <- filter_common_markers(
    data = input,
    filter.common.markers = filter.common.markers,
    path.folder = path.folder,
    verbose = TRUE
  )


  # Write tidy in the working directory
  radiator::write_rad(data = input, path = file.path(path.folder, filename))

  strata <- radiator::generate_strata(data = input, pop.id = FALSE)
  readr::write_tsv(x = strata, file = file.path(path.folder, strata.filename))

  # results --------------------------------------------------------------------
  message("\nMerged DArT file and strata file generated:")
  message("    ", filename)
  message("    ", strata.filename)

  message("\nNumber of chrom / locus / SNP:")
  message(length(unique(input$CHROM)), " / ", length(unique(input$LOCUS)), " / ",
          length(unique(input$MARKERS)), "\n")
  radiator::summary_strata(strata = strata)
  return(res = list(merged.dart = input, strata = strata))
}#End merge_dart




