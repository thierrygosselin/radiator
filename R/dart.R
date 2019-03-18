# extract_dart_target_id---------------------------------------------------------
#' @name extract_dart_target_id

#' @title Extract \href{http://www.diversityarrays.com}{DArT} target id

#' @description Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users. The function allows to extract DArT
#' target id from a DArT file. To help prepare the appropriate STRATA file.

#' @param data DArT output file. Note that most popular formats used by DArT are
#' recognised (1- and 2- rows format, also called binary, and count data.).
#' If you encounter a problem, sent me your data so that I can update
#' the function. The function can import \code{.csv} or \code{.tsv} files.

#' @param write With default \code{write = TRUE}, the dart target id column is
#' written in a file in the working directory.

#' @return A tidy dataframe with a \code{TARGET_ID} column. Spaces are remove and
#' UPPER case is used.

#' @export
#' @rdname extract_dart_target_id
#' @importFrom dplyr group_by select rename filter mutate summarise distinct n_distinct arrange left_join semi_join anti_join inner_join full_join tally bind_rows
#' @importFrom parallel detectCores
#' @importFrom stringi stri_replace_all_fixed stri_join stri_sub stri_replace_na stri_pad_left
#' @importFrom purrr discard
# @importFrom data.table as.data.table dcast.data.table
#' @importFrom readr read_tsv write_tsv read_lines read_table
#' @importFrom tibble as_data_frame data_frame
#' @importFrom tidyr spread gather unite separate

#' @examples
#' \dontrun{
#' clownfish.dart.tidy <- radiator::extract_dart_target_id(
#' data = "clownfish.dart.tsv")
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com} and Peter Grewe \email{peter.grewe@csiro.au}

extract_dart_target_id <- function(data, write = TRUE) {
  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("Input file missing")

  # Check DArT format file -----------------------------------------------------
  dart.check <- check_dart(data)
  if (!dart.check$data.type %in% c("dart", "silico.dart")) {
    rlang::abort("Contact author to show your DArT data, problem during import")
  } else {
    skip.number <- dart.check$skip.number
  }

  # Import data ---------------------------------------------------------------
  if (stringi::stri_detect_fixed(
    str = stringi::stri_sub(str = data, from = -4, to = -1),
    pattern = ".csv")) {
    csv <- TRUE
  } else {
    csv <- FALSE
  }

  if (csv) {
    dart.target.id <- readr::read_csv(
      file = data,
      skip = skip.number, n_max = 1,
      na = "-",
      col_names = FALSE,
      col_types = readr::cols(.default = readr::col_character()))
  } else {
    dart.target.id <- readr::read_tsv(
      file = data,
      skip = skip.number, n_max = 1,
      na = "-",
      col_names = FALSE,
      col_types = readr::cols(.default = readr::col_character()))
  }


  # This is a string of col header not wanted
  discard <- c(
    "TARGET_ID", "ALLELEID", "CLONEID", "CLUSTERTEMPINDEX", "ALLELESEQUENCE",
    "CLUSTERCONSENSUSSEQUENCE", "CLUSTERSIZE", "ALLELESEQDIST", "SNP",
    "SNPPOSITION", "CALLRATE", "ONERATIOREF", "ONERATIOSNP", "FREQHOMREF",
    "FREQHOMSNP", "FREQHETS", "PICREF", "PICSNP", "AVGPIC", "AVGCOUNTREF",
    "AVGCOUNTSNP", "RATIOAVGCOUNTREFAVGCOUNTSNP", "FREQHETSMINUSFREQMINHOM",
    "ALLELECOUNTSCORRELATION", "AGGREGATETAGSTOTAL", "DERIVEDCORRMINUSSEEDCORR",
    "REPREF", "REPSNP", "REPAVG", "PICREPREF", "PICREPSNP", "TOTALPICREPREFTEST",
    "TOTALPICREPSNPTEST", "BINID", "BIN.SIZE", "ALLELESEQUENCEREF",
    "ALLELESEQUENCESNP", "TRIMMEDSEQUENCEREF", "TRIMMEDSEQUENCE", "ONERATIO",
    "PIC", "AVGREADDEPTH", "STDEVREADDEPTH", "QPMR", "REPRODUCIBILITY", "MAF",
    "TOTCOUNTS")

  discard.genome <- c("CHROM_|CHROMPOS_|ALNCNT_|ALNEVALUE_")

  dart.target.id <- tidyr::gather(data = dart.target.id, key = DISCARD,
                                  value = TARGET_ID) %>%
    dplyr::select(-DISCARD) %>%
    dplyr::mutate(TARGET_ID = stringi::stri_trans_toupper(TARGET_ID)) %>%
    dplyr::filter(!TARGET_ID %in% discard) %>%
    dplyr::filter(stringi::stri_detect_regex(
      str = stringi::stri_trans_toupper(TARGET_ID),
      pattern = discard.genome, negate = TRUE)) %>%
    dplyr::mutate(TARGET_ID = stringi::stri_replace_all_fixed(
      TARGET_ID, pattern = c(" ", "_"), replacement = c("", "-"), vectorize_all = FALSE))

  if (write) readr::write_tsv(x = dart.target.id, path = "dart.target.id.tsv")

  # Check that DArT file as good target id written -----------------------------
  if (nrow(dart.target.id) != length(unique(dart.target.id$TARGET_ID))) {
    message("non unique target id are used in the DArT file...
What you want are different target ids at the end of the row that contains AlleleID, AlleleSequence, etc
Edit manually the DArT file before trying the functions: read_dart
If you're still encountering problem, email author for help")
  }

  # check if target id are numericnumeric
  # target.num <- unique(stringi::stri_detect_regex(
  #   str = unique(dart.target.id$TARGET_ID), pattern = "^[[:digit:]]+L"))
  # if (!target.num) {
  #   dart.target.id <- dart.target.id %>%
  #     dplyr::mutate(TARGET_ID = clean_ind_names(TARGET_ID))
  # }
  return(dart.target.id)
}#End extract_dart_target_id



## Merge dart-------------------------------------------------------------------

#' @title Merge DArT files
#' @description This function allows to merge 2 DArT files (filtered of not).

#' @param dart1 Full path of the first DArT file.
#' @param strata1 Full path of the first strata file for dart1.
#' @param dart2 Full path of the second DArT file.
#' @param strata2 Full path of the second strata file for dart2.
#' @param keep.rad Unless the default is changed, the function removes the
#' temporary tidy dart files generated for each DArT datasets during import.
#' Default: \code{keep.rad = FALSE}.
#' @param filename Name of the merged DArT file.
#' By default, the function gives the merged data the filename:\code{merge_dart}
#' with date and time appended. The function will also append the date and time
#' to the filename provided.
#' The data is written in the working directory.
#' Default: \code{filename = NULL}.
#' @inheritParams tidy_genomic_data
#' @inheritParams read_dart
#' @inheritParams radiator_common_arguments


#' @importFrom stringi stri_replace_all_fixed stri_replace_na stri_join stri_count_fixed
#' @importFrom tibble as_data_frame data_frame add_column add_row
#' @importFrom dplyr select rename n_distinct distinct mutate summarise group_by ungroup arrange left_join full_join semi_join anti_join bind_rows bind_cols if_else
#' @importFrom readr write_tsv read_tsv
#' @importFrom tidyr separate gather
#' @importFrom parallel detectCores
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot stat_summary labs theme element_blank element_text geom_jitter scale_colour_manual scale_y_reverse theme_light geom_bar facet_grid stat_smooth
#' @importFrom purrr flatten_chr map_df flatten_dbl

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
#' sum <- radiator::merge_dart(
#' dart1 = "bluefin_larvae.tsv", strata1 = "strata1_bft_larvae.tsv",
#' dart1 = "bluefin_adults.csv", strata2 = "strata2_bft_adults.tsv",
#' filename = "bluefin_combined")
#' }
#' @rdname merge_dart
#' @export
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

merge_dart <- function(
  dart1, strata1,
  dart2, strata2,
  #pop.select = NULL,
  # filter.monomorphic = TRUE,
  # filter.common.markers = TRUE,
  keep.rad = FALSE,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  ...
) {
  cat("#######################################################################\n")
  cat("########################## radiator::merge_dart #######################\n")
  cat("#######################################################################\n")
  timing <- proc.time()
  opt.change <- getOption("width")
  options(width = 70)
  # manage missing arguments -----------------------------------------------------
  if (missing(dart1)) rlang::abort("dart1 file missing")
  if (missing(dart2)) rlang::abort("dart2 file missing")
  if (missing(strata1)) rlang::abort("strata1 file missing")
  if (missing(strata2)) rlang::abort("strata2 file missing")

  # Filename -------------------------------------------------------------------
  # Get date and time to have unique filenaming
  if (is.null(filename)) {
    filename <- stringi::stri_join("merged_dart_", format(Sys.time(), "%Y%m%d@%H%M"), ".rad")
    strata.filename <- stringi::stri_join("merged_dart_strata_", format(Sys.time(), "%Y%m%d@%H%M"), ".tsv")
  } else {
    filename <- stringi::stri_join(filename, "_",format(Sys.time(), "%Y%m%d@%H%M"), ".rad")
    strata.filename <- stringi::stri_join(filename, "_strata_",format(Sys.time(), "%Y%m%d@%H%M"), ".tsv")
  }

  # File type detection---------------------------------------------------------
  data.type <- radiator::detect_genomic_format(dart1)

  # import data ----------------------------------------------------------------
  if (!is.null(pop.select)) {
    pop.select <- clean_pop_names(pop.select)
    message(stringi::stri_join(length(pop.select), "population(s) selected", sep = " "))
  }

  # Columns we want
  want <- c("MARKERS", "CHROM", "LOCUS", "POS", "INDIVIDUALS", "POP_ID", "GT_VCF_NUC",
            "CALL_RATE", "AVG_COUNT_REF", "AVG_COUNT_SNP", "REP_AVG")

  # the DArT data
  if (data.type == "dart") {
    message("Importing and tidying dart1...")
    input <- suppressWarnings(radiator::read_dart(
      data = dart1, strata = strata1,
      filename = "temp.tidy.dart1", verbose = FALSE) %>%
        dplyr::select(dplyr::one_of(want)))
    message("Importing and tidying dart2...")
    dart2 <- suppressWarnings(radiator::read_dart(
      data = dart2, strata = strata2,
      filename = "temp.tidy.dart2", verbose = FALSE) %>%
        dplyr::select(dplyr::one_of(want)))
    if (!keep.rad) {
      message("Removing temporary tidy DArT files...")
      file.remove("temp.tidy.dart1.rad")
      file.remove("temp.tidy.dart2.rad")
    }
  }

  # Import the filtered DArT data
  if (data.type == "fst.file") {
    message("Importing the filtered dart1...")
    input <- suppressWarnings(radiator::read_rad(dart1) %>% dplyr::select(dplyr::one_of(want)))
    message("Importing the filtered dart2...")
    dart2 <- suppressWarnings(radiator::read_rad(dart2) %>% dplyr::select(dplyr::one_of(want)))
  }

  # Keeping selected pop -------------------------------------------------------
  if (!is.null(pop.select)) {
    input <- suppressWarnings(dplyr::filter(input, POP_ID %in% pop.select))
    dart2 <- suppressWarnings(dplyr::filter(dart2, POP_ID %in% pop.select))
  }

  # cleaning up non-immortalized markers ---------------------------------------
  message("Removing markers starting with 1000 (non-immortalized DArT markers)")
  markers.before <- length(unique(input$MARKERS)) + length(unique(dart2$MARKERS))

  input <- suppressWarnings(dplyr::filter(input, !stringi::stri_detect_regex(str = LOCUS, pattern = "^1000")))
  dart2 <- suppressWarnings(dplyr::filter(dart2, !stringi::stri_detect_regex(str = LOCUS, pattern = "^1000")))
  markers.after <- length(unique(input$MARKERS)) + length(unique(dart2$MARKERS))
  message("    Markers removed: ", markers.before - markers.after)

  # merging DArT tidy data -----------------------------------------------------
  # check alt alleles
  # test1 <- dplyr::filter(input, (stringi::stri_count_fixed(str = ALT, pattern = ",") + 1) > 1)
  # test2 <- dplyr::filter(dart2, (stringi::stri_count_fixed(str = ALT, pattern = ",") + 1) > 1)

  # we keep common column
  dart.col <- dplyr::intersect(colnames(input), colnames(dart2))

  input <- suppressWarnings(
    dplyr::select(input, dplyr::one_of(dart.col)) %>%
      dplyr::bind_rows(dplyr::select(dart2, dplyr::one_of(dart.col)))
  )
  dart2 <- NULL

  # Pop select cleanup ---------------------------------------------------------
  if (!is.null(pop.select)) {
    if (is.factor(input$POP_ID)) input$POP_ID <- droplevels(input$POP_ID)
  }

  # Remove weird markers
  input <- radiator::detect_biallelic_problems(data = input, verbose = FALSE, parallel.core = parallel.core)$biallelic.data

  # Averaging across markers the call rate and other DArT markers metadata statistics
  # Note to myself: Might be easier/faster to use mutate_if
  if (rlang::has_name(input, "CALL_RATE")) {
    message("Averaging across markers the call rate")
    input <- input %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(CALL_RATE = mean(CALL_RATE)) %>%
      dplyr::ungroup(.)
  }

  if (rlang::has_name(input, "REP_AVG")) {
    message("Averaging across markers the REP_AVG")
    input <- input %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(REP_AVG = mean(REP_AVG)) %>%
      dplyr::ungroup(.)
  }

  message("Adjusting REF/ALT alleles...")
  input <- radiator::calibrate_alleles(
    data = input, biallelic = NULL,
    parallel.core = parallel.core,
    verbose = TRUE)$input

  if (rlang::has_name(input, "POLYMORPHIC")) {
    input <- dplyr::select(input, -POLYMORPHIC)
  }

  if (rlang::has_name(input, "AVG_COUNT_REF")) {
    message("Averaging across markers the coverage for the REF allele")
    input <- input %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(AVG_COUNT_REF = mean(AVG_COUNT_REF)) %>%
      dplyr::ungroup(.)
  }

  if (rlang::has_name(input, "AVG_COUNT_ALT")) {
    message("Averaging across markers the coverage for the ALT allele")
    input <- input %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(AVG_COUNT_REF = mean(AVG_COUNT_REF)) %>%
      dplyr::ungroup(.)
  }

  if (filter.monomorphic) {
    test <- radiator::filter_monomorphic(data = input, verbose = TRUE)
  }

  if (filter.common.markers) {
    input <- radiator::filter_common_markers(data = input, fig = TRUE, verbose = TRUE) %$% input
    message("\n    **** Manually save the figure ****")
  }

  # Write tidy in the working directory
  radiator::write_rad(data = input, path = filename)

  strata <- dplyr::select(input, INDIVIDUALS, STRATA = POP_ID) %>%
    dplyr::distinct(INDIVIDUALS, STRATA) %>%
    readr::write_tsv(x = ., path = strata.filename)

  # results --------------------------------------------------------------------
  message("\nMerged DArT file and strata file in the working directory:")
  message("    ", filename)
  message("    ", strata.filename)
  timing <- proc.time() - timing
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("############################## completed ##############################\n")
  options(width = opt.change)
  return(res = list(merged.dart = input, strata = strata))
}#End merge_dart


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
#' @importFrom dplyr group_by select rename filter mutate summarise distinct n_distinct arrange left_join semi_join anti_join inner_join full_join tally bind_rows
#' @importFrom parallel detectCores
#' @importFrom stringi stri_replace_all_fixed stri_join stri_sub stri_replace_na stri_pad_left
#' @importFrom purrr discard
# @importFrom data.table as.data.table dcast.data.table
#' @importFrom readr read_tsv write_tsv read_lines read_table
#' @importFrom tibble as_data_frame data_frame
#' @importFrom tidyr spread gather unite separate

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
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (verbose) message("Execution date@time: ", file.date)
  old.dir <- getwd()
  opt.change <- getOption("width")
  options(width = 70)
  timing <- proc.time()# for timing
  # res <- list()
  #back to the original directory and options
  on.exit(setwd(old.dir), add = TRUE)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(timing <- proc.time() - timing, add = TRUE)
  on.exit(if (verbose) message("\nTiming: ", round(timing[[3]]), " sec"), add = TRUE)
  on.exit(if (verbose) cat("############################## completed ##############################\n"), add = TRUE)


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
    if (stringi::stri_detect_fixed(
      str = stringi::stri_sub(str = data, from = -4, to = -1),
      pattern = ".csv")) {
      csv <- TRUE
    } else {
      csv <- FALSE
    }

    dart.check <- check_dart(data)
    if (!dart.check$data.type %in% c("dart", "silico.dart")) {
      rlang::abort("Contact author to show your DArT data, problem during import")
    } else {
      skip.number <- dart.check$skip.number
    }


    if (csv) {
      dart.col.type <- readr::read_csv(
        file = data,
        skip = skip.number, n_max = 1,
        na = "-",
        col_names = FALSE,
        col_types = readr::cols(.default = readr::col_character()))
    } else {
      dart.col.type <- readr::read_tsv(
        file = data,
        skip = skip.number, n_max = 1,
        na = "-",
        col_names = FALSE,
        col_types = readr::cols(.default = readr::col_character()))
    }

    want <- tibble::tibble(
      INFO = c("ALLELEID", "SNP", "SNPPOSITION", "CALLRATE",
               "AVGCOUNTREF", "AVGCOUNTSNP", "REPAVG"),
      COL_TYPE = c("c", "c", "i", "d", "d", "d", "d"))

    dart.col.type <- dart.col.type %>%
      tidyr::gather(data = .,key = DELETE, value = INFO) %>%
      dplyr::select(-DELETE) %>%
      dplyr::mutate(INFO = stringi::stri_trans_toupper(INFO)) %>%
      dplyr::left_join(want, by = "INFO") %>%
      dplyr::mutate(COL_TYPE = stringi::stri_replace_na(str = COL_TYPE, replacement = "_")) %>%
      dplyr::select(COL_TYPE) %>%
      purrr::flatten_chr(.) %>% stringi::stri_join(collapse = "")

    if (csv) {
      input <- suppressMessages(suppressWarnings(
        readr::read_csv(
          file = data,
          skip = skip.number,
          na = "-",
          col_names = TRUE,
          col_types = dart.col.type)
      ))
    } else {
      input <- suppressMessages(suppressWarnings(
        readr::read_tsv(
          file = data,
          skip = skip.number,
          na = "-",
          col_names = TRUE,
          col_types = dart.col.type)
      ))
    }

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

    input <- suppressWarnings(
      input %>%
        tidyr::separate(col = LOCUS, into = c("LOCUS", "NOT_USEFUL"), sep = "\\|", extra = "drop") %>%
        dplyr::select(-NOT_USEFUL) %>%
        tidyr::separate(col = SNP, into = c("NOT_USEFUL", "KEEPER"), sep = ":", extra = "drop") %>%
        dplyr::select(-NOT_USEFUL) %>%
        tidyr::separate(col = KEEPER, into = c("REF", "ALT"), sep = ">") %>%
        dplyr::mutate(
          CHROM = rep("CHROM_1", n()),
          MARKERS = stringi::stri_join(CHROM, LOCUS, POS, sep = "__")) %>%
        dplyr::select(dplyr::one_of(want), dplyr::everything()) %>%
        dplyr::mutate_at(.tbl = ., .vars = c("MARKERS", "CHROM", "LOCUS", "POS"), .funs = as.character) %>%
        dplyr::arrange(CHROM, LOCUS, POS, REF) %>%
        dplyr::filter(!is.na(REF) | !is.na(ALT)) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
        dplyr::arrange(MARKERS))
  }

  if (data.type == "fst.file") {
    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "REF", "ALT", "CALL_RATE",
              "AVG_COUNT_REF", "AVG_COUNT_SNP", "REP_AVG")
    input <- suppressWarnings(
      radiator::read_rad(data) %>%
        dplyr::select(dplyr::one_of(want)) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE))
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

# Tidy dart -------------------------------------------------------------------------
# Import, filter and transform a dart output file to different formats

#' @name read_dart

#' @title Tidy \href{http://www.diversityarrays.com}{DArT} output file.

#' @description Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users. The function generate a tidy dataset using
#' \href{http://www.diversityarrays.com}{DArT} file.

#' @param data DArT output file. Note that most popular formats used by DArT are
#' recognised (1- and 2- row format, also called binary, and count data.).
#' If you encounter a problem, sent me your data so that I can update
#' the function. The function can import \code{.csv} or \code{.tsv} files.

#' @param strata A tab delimited file or object with 3 columns.
#' Columns header is:
#' \code{TARGET_ID}, \code{INDIVIDUALS} and \code{STRATA}.
#' Note: the column \code{STRATA} refers to any grouping of individuals.
#' You need to make sure that
#' the column \code{TARGET_ID} match the id used by DArT.
#' The column \code{INDIVIDUALS} and \code{STRATA} will be kept in the tidy data.
#' Only individuals in the strata file are kept in the tidy, i.e. that the strata
#' is also used as a whitelist of individuals/strata.

#' @inheritParams tidy_genomic_data
#' @inheritParams radiator_common_arguments
#' @inheritParams filter_whitelist

#' @param ... (optional) To pass further argument for fine-tuning the function.


#' @return A radiator GDS file and tidy dataframe with several columns depending on DArT file:
#' \enumerate{
#' \item MARKERS: generated by radiator and correspond to CHROM + LOCUS + POS separated by 2 underscores.
#' \item CHROM: the chromosome, for de novo: CHROM_1.
#' \item LOCUS: the locus.
#' \item POS: the SNP id on the LOCUS.
#' \item REF: the reference allele.
#' \item ALT: the alternate allele.
#' \item INDIVIDUALS: the sample name.
#' \item POP_ID: populations id of the sample.
#' \item GT: the genotype in 6 digit format à la genepop.
#' \item GT_VCF: the genotype in VCF format.
#' \item GT_VCF_NUC: the genotype in VCF format, but keeping the nucleotide information.
#' \item GT_BIN: the genotype in binary format similar to PLINK. The number correspond to the number of alternate allele in the genotype.
#' \item CALL_RATE: call rate output specific of DArT.
#' \item AVG_COUNT_REF: the coverage for the reference allele, output specific of DArT.
#' \item AVG_COUNT_SNP: the coverage for the alternate allele, output specific of DArT.
#' \item REP_AVG: the reproducibility average, output specific of DArT.
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

#' \strong{Advance mode, using \emph{dots-dots-dots ...}}
#' \enumerate{
#' \item \code{whitelist.markers}: detailed in \code{\link[radiator]{filter_whitelist}}.
#' Defautl: \code{whitelist.markers = NULL}.
#' \item \code{missing.memory} (option, path)
#' This argument allows to erase genotypes that have bad statistics.
#' It's the path to a file \code{.rad} file that contains 3 columns:
#' \code{MARKERS, INDIVIDUALS, ERASE}. The file is produced by several radiator
#' functions. For DArT data, \code{\link[radiator]{filter_rad}} generate the file.
#' Defautl: \code{missing.memory = NULL}.
#' \item \code{dart.sequence}: (optional, logical) To keep the sequence of markers in
#' the GDS. Theres no computational cost to keep it in the GDS...
#' Default: \code{dart.sequence = TRUE}.
#' \item \code{path.folder}: (optional, path) To write output in a specific folder.
#' Default: \code{path.folder = NULL}. With defautl, the working directory is used.
#' }


#' @export
#' @rdname read_dart
#' @importFrom dplyr group_by select rename filter mutate summarise distinct n_distinct arrange left_join semi_join anti_join inner_join full_join tally bind_rows
#' @importFrom parallel detectCores
#' @importFrom stringi stri_replace_all_fixed stri_join stri_sub stri_replace_na stri_pad_left
#' @importFrom purrr discard
#' @importFrom data.table as.data.table dcast.data.table melt.data.table
#' @importFrom readr read_tsv write_tsv read_lines read_table
#' @importFrom tibble as_data_frame data_frame
#' @importFrom tidyr spread gather unite separate

#' @examples
#' \dontrun{
#' clownfish.dart.tidy <- radiator::read_dart(
#' data = "clownfish.dart.tsv",
#' strata = "clownfish.strata.tsv",
#' verbose = TRUE)
#'
#' # To get a strip, bare minimal version of the DArT tidy data:
#' clownfish.dart.tidy <- radiator::read_dart(
#' data = "clownfish.dart.tsv",
#' strata = "clownfish.strata.tsv",
#' tidy = "strip")
#' # Only MARKERS, INDIVIDUALS, POP_ID, GT_BIN is generated in the output
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

read_dart <- function(
  data,
  strata,
  filename = NULL,
  verbose = FALSE,
  parallel.core = parallel::detectCores() - 1,
  ...
) {

  # for testing
  # filename = NULL
  # verbose = TRUE
  # parallel.core = parallel::detectCores() - 1
  # whitelist.markers = NULL
  # missing.memory <- NULL
  # dart.sequence <- TRUE
  # path.folder = NULL
  # radiator.pipeline = NULL
  # internal <- FALSE
  # tidy.dart = FALSE
  # tidy.check = FALSE
  # gt = NULL
  # gt.bin = NULL
  # gt.vcf = NULL
  # gt.vcf.nuc = NULL

  if (verbose) {
    cat("################################################################################\n")
    cat("############################## radiator::read_dart #############################\n")
    cat("################################################################################\n")
  }

  # Cleanup-------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (verbose) message("Execution date@time: ", file.date)
  old.dir <- getwd()
  opt.change <- getOption("width")
  options(width = 70)
  timing <- proc.time()# for timing
  #back to the original directory and options
  on.exit(setwd(old.dir), add = TRUE)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(timing <- proc.time() - timing, add = TRUE)
  on.exit(message("\nDArT conversion timing: ", round(timing[[3]]), " sec"), add = TRUE)
  on.exit(if (verbose) cat("############################## completed read_dart #############################\n"), add = TRUE)

  # Function call and dotslist -------------------------------------------------
  rad.dots <- radiator_dots(
    func.name = as.list(sys.call())[[1]],
    fd = rlang::fn_fmls_names(),
    args.list = as.list(environment()),
    dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE),
    keepers = c("whitelist.markers", "missing.memory", "dart.sequence",
                "path.folder", "internal", "pop.levels",
                "gt", "gt.bin", "gt.vcf", "gt.vcf.nuc",
                "tidy.dart", "tidy.check"
    ),
    verbose = verbose
  )
  if (is.null(tidy.dart)) tidy.dart <- FALSE
  if (is.null(tidy.check)) tidy.check <- FALSE
  if (tidy.check) tidy.dart <- TRUE
  if (!tidy.dart) {
    gt <- gt.vcf <- gt.vcf.nuc <- FALSE
  }

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

  # radiator.folder <- generate_folder(
  #   f = path.folder,
  #   rad.folder = "import_dart",
  #   prefix_int = TRUE,
  #   internal = FALSE,
  #   file.date = file.date,
  #   verbose = verbose)

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

  # initialization ------
  dummy1 <- data.frame(
    ERASE = character(0),
    GT_BIN = integer(0),
    GT = character(0),
    GT_VCF = character(0),
    GT_VCF_NUC = character(0),
    READ_DEPTH = integer(0),
    ALLELE_REF_DEPTH = integer(0),
    ALLELE_ALT_DEPTH = integer(0)
  )

  # Import data ---------------------------------------------------------------
  data <- import_dart(
    data = data,
    strata = strata,
    pop.levels = pop.levels,
    verbose = verbose
  )
  strata <- data$strata
  dart.format <- data$dart.format
  data <- data$data

  # Tidying data ---------------------------------------------------------------
  # Whitelist ------------------------------------------------------------------
  if (!is.null(whitelist.markers)) {
    data %<>% filter_whitelist(data = ., whitelist.markers = whitelist.markers)
  }
  # STRATEGY tidy dart  ---------------------------------------------------------
  # Depending on the number of markers ...
  # All this can be overwritten in ... argument
  # gt.bin is the dosage of ALT allele: 0, 1, 2 NA

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

  # 1-row format----------------------------------------------------------------
  if (dart.format == "1row") {
    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "REF", "ALT",
              "CALL_RATE", "AVG_COUNT_REF", "AVG_COUNT_SNP", "REP_AVG",
              "SEQUENCE")
    data <- suppressWarnings(
      data.table::as.data.table(data) %>%
        data.table::melt.data.table(
          data = .,
          id.vars = want,
          variable.name = "TARGET_ID", variable.factor = FALSE,
          value.name = "GT"
        ) %>%
        tibble::as_tibble(.)
    )

    # markers metadata
    suppressWarnings(
      dplyr::ungroup(data) %>%
        dplyr::select(dplyr::one_of(want)) %>%
        dplyr::filter(!is.na(REF) | !is.na(ALT)) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
        dplyr::arrange(MARKERS) %>%
        radiator::write_rad(
          data = ., path = meta.filename$filename))
    if (verbose) message("File written: ", meta.filename$filename.short)

    notwanted <- c("CHROM", "LOCUS", "POS", "CALL_RATE", "AVG_COUNT_REF",
                   "AVG_COUNT_SNP", "REP_AVG", "SEQUENCE")

    suppressWarnings(
      data  %<>%
        dplyr::select(-dplyr::one_of(notwanted)) %>%
        dplyr::arrange(MARKERS) %>%
        dplyr::left_join(
          dplyr::distinct(data, MARKERS) %>%
            dplyr::mutate(SPLIT_VEC = split_vec_row(., 4, parallel.core = parallel.core))
          , by = "MARKERS") %>%
        split(x = ., f = .$SPLIT_VEC)
    )

    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "REF", "ALT", "TARGET_ID",
              "CALL_RATE", "AVG_COUNT_REF", "AVG_COUNT_SNP", "REP_AVG",
              "GT", "GT_VCF", "GT_VCF_NUC", "GT_BIN",
              "READ_DEPTH", "ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH", "SEQUENCE")
    data <- suppressWarnings(
      .radiator_parallel(
        X = data,
        FUN = dart2gt,
        mc.cores = parallel.core,
        dart.format = "1row"
      ) %>%
        dplyr::bind_rows(.) %>%
        dplyr::left_join(radiator::read_rad(meta.filename$filename), by = "MARKERS")
    )
  }#End 1 row format DArT file

  # Binary dart counts and 2-rows format-----------------------------------------
  if (dart.format != "1row") {
    # Markers meta
    want <- c("VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS", "COL",
              "REF", "ALT", "SEQUENCE",
              "CALL_RATE", "AVG_COUNT_REF", "AVG_COUNT_SNP", "REP_AVG")
    markers.meta <- suppressWarnings(
      dplyr::ungroup(data) %>%
        dplyr::select(dplyr::one_of(want)) %>%
        dplyr::filter(!is.na(REF) | !is.na(ALT)) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
        dplyr::arrange(MARKERS) %>%
        dplyr::select(-AVG_COUNT_REF, -AVG_COUNT_SNP) %>%
        dplyr::mutate(FILTERS = "whitelist") %>%
        dplyr::select(FILTERS, VARIANT_ID, MARKERS, CHROM, LOCUS, POS, COL, REF, ALT, dplyr::everything(.))
    )
    write_rad(
      data = markers.meta,
      path = meta.filename$filename,
      filename = meta.filename$filename.short,
      write.message = "standard",
      verbose = verbose
    )
    # if (verbose) message("File written: ", meta.filename$)
    notwanted <- c("FILTERS", "CHROM", "LOCUS", "POS", "COL", "CALL_RATE", "AVG_COUNT_REF",
                   "AVG_COUNT_SNP", "REP_AVG", "SEQUENCE")
    suppressWarnings(
      data %<>%
        dplyr::select(-dplyr::one_of(notwanted)) %>%
        dplyr::arrange(MARKERS, REF))

    #GDS -----------------------------------------------------------------------
    data <- dart2gds(
      data = data,
      strata = strata,
      markers.meta = markers.meta,
      filename = filename.gds$filename,
      dart.format = dart.format,
      gt.vcf = gt.vcf,
      gt.vcf.nuc = gt.vcf.nuc,
      gt = gt,
      parallel.core = parallel.core,
      verbose = verbose
    )
  }#End binary

  # Tidy data -----------------------------------------------------------------
  n.chrom <- length(unique(markers.meta$CHROM))
  n.locus <- length(unique(markers.meta$LOCUS))
  n.snp <- length(unique(markers.meta$MARKERS))


  # message("Closing the GDS connection")
  # SeqArray::seqClose(object = data.gds)

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
          dplyr::select(-dplyr::one_of(notwanted))
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
    #
    # dplyr::distinct(tidy.data, MARKERS, CHROM, LOCUS, POS, REF, ALT) %>%
    #   readr::write_tsv(x = ., path = dic.filename$filename)
    # if (verbose) message("File written: ", dic.filename$filename.short)


    # Final strata ---------------------------------------------------------
    strata.filename <- generate_filename(
      name.shortcut = "radiator_tidy_dart_strata",
      path.folder = path.folder,
      date = TRUE,
      extension = "tsv")
    strata <- extract_individuals_metadata(gds = data, whitelist = TRUE)
    readr::write_tsv(x = strata, path = strata.filename$filename)

    if (!is.null(strata)) {
      tidy.data <- join_strata(
        data = tidy.data,
        strata = dplyr::select(strata, STRATA, INDIVIDUALS),
        pop.id = TRUE,
        verbose = FALSE
      )

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
      same.individuals<- identical(unique(missing$INDIVIDUALS), unique(data$INDIVIDUALS))
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

# INTERNAL FUNCTIONS------------------------------------------------------------
#' @title import_dart
#' @description Read DArT file
#' @rdname import_dart
#' @keywords internal
#' @export
import_dart <- function(data, strata, pop.levels = NULL, path.folder= NULL, verbose = TRUE) {
  message("Reading DArT file...")

  # Check that DArT file as good target id written -----------------------------
  target.id <- extract_dart_target_id(data, write = FALSE)
  n.ind.dart <- nrow(target.id)
  if (verbose) message("    Number of individuals: ", n.ind.dart)
  if (nrow(target.id) != length(unique(target.id$TARGET_ID))) {
    rlang::abort("\nnon unique target id are used in the DArT file...
                 What you want are different target ids at the end of the row that contains AlleleID, AlleleSequence, etc
                 Edit manually before trying again
                 If you're still encountering problem, email author for help")
  }

  # Check DArT format file -----------------------------------------------------
  dart.check <- check_dart(data)
  if (!dart.check$data.type %in% c("dart", "silico.dart")) {
    rlang::abort("\nContact author to show your DArT data, problems during import")
  } else {
    skip.number <- dart.check$skip.number
  }
  dart.check <- NULL

  # Strata file ------------------------------------------------------------------
  strata.df <- radiator::read_strata(
    strata = strata,
    pop.levels = pop.levels, pop.labels = NULL,
    pop.select = NULL, blacklist.id = NULL,
    keep.two = FALSE, verbose = verbose)
  pop.levels <- strata.df$pop.levels
  pop.labels <- strata.df$pop.labels
  pop.select <- strata.df$pop.select
  blacklist.id <- strata.df$blacklist.id
  strata.df <- strata.df$strata
  n.ind.strata <- nrow(strata.df)

  # Check that TARGET_ID in strata match TARGET_ID in the DArT file ------------
  if (verbose) message("Using individuals in strata file to filter individuals in DArT file")
  if (n.ind.dart != n.ind.strata) {
    strata.id.check <- strata.df %>%
      dplyr::mutate(IN_DART = stringi::stri_trans_toupper(strata.df$TARGET_ID)
                    %in% stringi::stri_trans_toupper(target.id$TARGET_ID))
    strata.id.pass <- !FALSE %in% (unique(strata.id.check$IN_DART))
    if (!strata.id.pass) {
      problem.filename <- generate_filename(
        name.shortcut = "radiator_tidy_dart_strata_problem",
        path.folder = path.folder,
        date = TRUE,
        extension = "tsv")

      readr::write_tsv(
        x = strata.id.check,
        path = problem.filename$filename)
      rlang::abort("\nSome of the samples in the strata are not found in the DArT file.
                     For more info: ", problem.filename$filename.short)
    }
    if (verbose) message("\nNote: Careful if using DArT statistics generated for all samples...\n")
    strata.id.check <- NULL
  } else {
    if (!identical(sort(target.id$TARGET_ID), sort(strata.df$TARGET_ID))) {
      rlang::abort("\nThe DArT and strata files don't have the same TARGET_IDs")
    }
  }
  target.id <- NULL

  # need to check for duplicate names... yes happening all the time
  duplicate.id.strata <- length(strata.df$INDIVIDUALS) - length(unique(strata.df$INDIVIDUALS))

  if (duplicate.id.strata > 0) {
    message("Duplicated individuals names found in the strata.\n   number of duplicate names = ", duplicate.id.strata, "\n")
    rlang::abort("\nFix the strata with unique names and\nverify the DArT file for the same issue, adjust accordingly...")
  }

  if (stringi::stri_detect_fixed(
    str = stringi::stri_sub(str = data, from = -4, to = -1),
    pattern = ".csv")) {
    csv <- TRUE
  } else {
    csv <- FALSE
  }
  if (csv) {
    dart.col.type <- readr::read_csv(
      file = data,
      skip = skip.number, n_max = 1,
      na = "-",
      col_names = FALSE,
      col_types = readr::cols(.default = readr::col_character()))
  } else {
    dart.col.type <- readr::read_tsv(
      file = data,
      skip = skip.number, n_max = 1,
      na = "-",
      col_names = FALSE,
      col_types = readr::cols(.default = readr::col_character()))
  }
  dart.sequence <- TRUE
  # if (dart.sequence) {
  info <- c("ALLELEID", "SNP", "SNPPOSITION", "CALLRATE",
            "AVGCOUNTREF", "AVGCOUNTSNP", "REPAVG", "CLONEID", "AVGREADDEPTH",
            "REPRODUCIBILITY", "CLUSTERCONSENSUSSEQUENCE")
  info.type <- c("c", "c", "i", "d", "d", "d", "d", "c", "d", "d", "c")
  want <- tibble::tibble(INFO = info, COL_TYPE = info.type) %>%
    dplyr::bind_rows(
      dplyr::select(strata.df, INFO = TARGET_ID) %>%
        dplyr::mutate(
          COL_TYPE = rep("c", n()),
          INFO = stringi::stri_trans_toupper(INFO)))

  dart.col.type %<>%
    tidyr::gather(data = .,key = DELETE, value = INFO) %>%
    dplyr::select(-DELETE) %>%
    dplyr::mutate(
      INFO = stringi::stri_trans_toupper(INFO),
      INFO = stringi::stri_replace_all_fixed(
        str = INFO, pattern = c(" ", "_"),
        replacement = c("", "-"), vectorize_all = FALSE)
    ) %>%
    dplyr::left_join(want, by = "INFO") %>%
    dplyr::mutate(COL_TYPE = stringi::stri_replace_na(str = COL_TYPE, replacement = "_")) %>%
    dplyr::select(COL_TYPE) %>%
    purrr::flatten_chr(.) %>% stringi::stri_join(collapse = "")
  want <- NULL

  if (csv) {
    data <- suppressMessages(suppressWarnings(
      readr::read_csv(
        file = data,
        skip = skip.number,
        na = c("-", " ", "", "NA"),
        col_names = TRUE,
        col_types = dart.col.type)
    ))
  } else {
    data <- suppressMessages(suppressWarnings(
      readr::read_tsv(
        file = data,
        skip = skip.number,
        na = c("-", " ", "", "NA"),
        col_names = TRUE,
        col_types = dart.col.type)
    ))
  }
  dart.col.type <- NULL
  colnames(data) %<>%
    stringi::stri_trans_toupper(str = .) %>%
    stringi::stri_replace_all_fixed(
      str = .,
      pattern = c(" ", "_"),
      replacement = c("", "-"),
      vectorize_all = FALSE
    ) %>%
    stringi::stri_replace_all_fixed(
      str = .,
      pattern = c("AVGCOUNTREF", "AVGCOUNTSNP", "REPAVG", "ALLELEID",
                  "SNPPOSITION", "CALLRATE", "REPRODUCIBILITY", "AVGREADDEPTH",
                  "CLUSTERCONSENSUSSEQUENCE"),
      replacement = c("AVG_COUNT_REF", "AVG_COUNT_SNP", "REP_AVG", "LOCUS",
                      "POS", "CALL_RATE", "REP_AVG", "AVG_READ_DEPTH",
                      "SEQUENCE"),
      vectorize_all = FALSE)

  if (rlang::has_name(data, "CLONEID")) {
    if (rlang::has_name(data, "LOCUS")) {
      data %<>% dplyr::select(-CLONEID)
    } else {
      colnames(data) %<>%
        stringi::stri_replace_all_fixed(
          str = .,
          pattern = "CLONEID",
          replacement = "LOCUS",
          vectorize_all = FALSE)
    }
  }

  if (!rlang::has_name(data, "LOCUS")) {
    rlang::abort("\nProblem tidying DArT dataset: contact author")
  }
  # necessary steps...observed with DArT file using ref genome -----------------
  data %<>% dplyr::filter(!is.na(LOCUS))
  if (rlang::has_name(data, "POS")) {
    data %<>% dplyr::arrange(LOCUS, POS)
  } else {
    data %<>% dplyr::arrange(LOCUS)
  }

  # Check for duplicate rows (sometimes people combine DArT data...)----------
  data.dup <- nrow(dplyr::distinct(data, LOCUS, SNP, POS, CALL_RATE, .keep_all = FALSE))

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
                               "AVG_COUNT_SNP", "REP_AVG", "SEQUENCE")
  )
  dup.id <- length(id) - length(unique(id))
  if (dup.id > 0) {
    rlang::abort(stringi::stri_join("Duplicated individual names in the data: ", dup.id))
  }
  # removing unused object
  remove.list <- id <- dup.id <- NULL

  # clean locus, generate MARKERS and VARIANT_ID
  data %<>% clean_dart_locus(.)

  # DArT characteristics--------------------------------------------------------
  if (verbose)  message("\nDArT characteristics:")
  # dart.format:
  # "1row" !binary
  # "2rows" binary
  # "counts" binary

  # Determine the type of DArT file: binary (2-row formats) or not 1 (row genotypes)
  dart.format <- detect_dart_format(
    x = data,
    target.id = strata.df$TARGET_ID,
    verbose = TRUE)

  return(res = list(data = data, strata = strata.df, dart.format = dart.format))
}#End import_dart



# clean_dart_locus--------------------------------------------------------------
#' @title clean_dart_locus
#' @description Clean LOCUS and generate VARIANT_ID and MARKERS
#' @rdname clean_dart_locus
#' @keywords internal
#' @export
clean_dart_locus <- function(x) {
  want <- c("VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS", "COL",
            "REF", "ALT", "SEQUENCE",
            "CALL_RATE", "AVG_COUNT_REF", "AVG_COUNT_SNP", "REP_AVG")
  suppressWarnings(
    x %<>%
      tidyr::separate(col = LOCUS,
                      into = c("LOCUS", "NOT_USEFUL"),
                      sep = "\\|",
                      extra = "drop"
      ) %>%
      dplyr::select(-NOT_USEFUL) %>%
      tidyr::separate(col = SNP,
                      into = c("NOT_USEFUL", "KEEPER"),
                      sep = ":",
                      extra = "drop") %>%
      dplyr::select(-NOT_USEFUL) %>%
      tidyr::separate(col = KEEPER, into = c("REF", "ALT"), sep = ">") %>%
      dplyr::mutate(
        CHROM = rep("CHROM_1", n()),
        MARKERS = stringi::stri_join(CHROM, LOCUS, POS, sep = "__"),
        VARIANT_ID = as.integer(factor(MARKERS)),
        COL = POS
      ) %>%
      dplyr::select(dplyr::one_of(want), dplyr::everything()) %>%
      dplyr::mutate_at(
        .tbl = .,
        .vars = c("MARKERS", "CHROM", "LOCUS", "POS"),
        .funs = as.character
      ) %>%
      dplyr::arrange(CHROM, LOCUS, POS, REF)
  )
}#End clean_dart_locus

# detect_dart_format-------------------------------------------------------------
#' @title detect_dart_format
#' @description Detect the dart genotype format: 1row, 2rows or counts
#' @rdname detect_dart_format
#' @keywords internal
#' @export
detect_dart_format <- function(x = NULL, target.id = NULL, verbose = TRUE) {
  binary <- anyNA(x$REF)

  if (!binary) {
    if (verbose) message("DArT SNP format: genotypes in 1 Row")
    dart.format <- "1row"
  } else {
    count.data <- x %>%
      dplyr::select(
        dplyr::one_of(
          sample(x = target.id, size = min(10, floor(0.1 * length(target.id))))
        )
      ) %>%
      dplyr::mutate_all(.tbl = ., .funs = as.numeric) %>%
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
  return(dart.format)


  # prior code:
  # # keep one marker and check if genotypes are count data
  # want <- c("VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS", "REF", "ALT", "CALL_RATE",
  #           "AVG_COUNT_REF", "AVG_COUNT_SNP", "REP_AVG", "SEQUENCE")
  # count.data <- 0 #required to start the while loop
  # while (count.data == 0) {
  #   count.data <- suppressWarnings(dplyr::select(
  #     x, -dplyr::one_of(want)) %>%
  #       dplyr::sample_n(tbl = ., size = 3) %>%
  #       purrr::flatten_chr(.) %>%
  #       as.integer %>%
  #       unique %>%
  #       sum(na.rm = TRUE))
  # }
  # count.data <- count.data > 3
}#End detect_dart_format
# dart2gt----------------------------------------------------------------
#' @title dart2gt
#' @description Transform dart genotypes to radiator genotypes fields
#' @rdname dart2gt
#' @keywords internal
#' @export
dart2gt <- function(x, dart.format, write.tidy = "full") {
  # x <- input
  # 1 row format----------------------------------------------------------------
  if (dart.format == "1row") {
    x <- dplyr::select(x, -SPLIT_VEC)
    if (write.tidy == "full") {
      x <- dplyr::bind_rows(
        dplyr::filter(x, is.na(GT)) %>%
          dplyr::mutate(
            GT = "000000",
            GT_VCF = "./.",
            GT_BIN = NA_integer_,
            GT_VCF_NUC = "./.",
            REF = NULL,
            ALT = NULL
          ),
        dplyr::filter(x, !is.na(GT)) %>%
          dplyr::mutate(
            GT = stringi::stri_replace_all_fixed(
              str = GT, pattern = c("0", "1", "2"),
              replacement = c("RR", "AA", "RA"),
              vectorize_all = FALSE),
            GT_VCF = stringi::stri_replace_all_fixed(
              str = GT, pattern = c("RR", "AA", "RA"),
              replacement = c("0/0", "1/1", "0/1"),
              vectorize_all = FALSE),
            GT_BIN = stringi::stri_replace_all_fixed(
              str = GT, pattern = c("RR", "AA", "RA"),
              replacement = c("0", "2", "1"),
              vectorize_all = FALSE),
            GT_BIN = as.integer(GT_BIN),
            GT_VCF_NUC = dplyr::if_else(GT_BIN == "0", stringi::stri_join(REF, REF, sep = "/"),
                                        dplyr::if_else(GT_BIN == "2", stringi::stri_join(ALT, ALT, sep = "/"),
                                                       stringi::stri_join(REF, ALT, sep = "/")), "./."),
            REF = stringi::stri_replace_all_regex(
              str = REF,
              pattern = c("^A$", "^C$", "^G$", "^T$"),
              replacement = c("001", "002", "003", "004"),
              vectorize_all = FALSE),
            ALT = stringi::stri_replace_all_regex(
              str = ALT,
              pattern = c("^A$", "^C$", "^G$", "^T$"),
              replacement = c("001", "002", "003", "004"),
              vectorize_all = FALSE),
            GT = dplyr::if_else(GT == "RR", stringi::stri_join(REF, REF, sep = ""),
                                dplyr::if_else(GT == "AA", stringi::stri_join(ALT, ALT, sep = ""),
                                               stringi::stri_join(REF, ALT, sep = "")), "000000"),
            REF = NULL,
            ALT = NULL
          ))
    } else {
      x <- dplyr::bind_rows(
        dplyr::filter(x, is.na(GT)) %>%
          dplyr::mutate(GT = NA_integer_),
        dplyr::filter(x, !is.na(GT)) %>%
          dplyr::mutate(
            GT = stringi::stri_replace_all_fixed(
              str = GT, pattern = c("0", "1", "2"),
              replacement = c("RR", "AA", "RA"),
              vectorize_all = FALSE),
            GT = stringi::stri_replace_all_fixed(
              str = GT, pattern = c("RR", "AA", "RA"),
              replacement = c("0", "2", "1"),
              vectorize_all = FALSE),
            GT = as.integer(GT),
            REF = NULL,
            ALT = NULL
          )) %>%
        dplyr::rename(GT_BIN = GT)
    }
  }#End 1row

  # Common between 2 rows and counts -------------------------------------------
  if (dart.format %in% c("2rows", "counts")) {
    # x <- input
    x %<>% dplyr::select(-SPLIT_VEC) %>%
      dplyr::arrange(MARKERS, REF) %>%
      dplyr::mutate(TEMP = rep(1:2, n()/2)) %>%
      dplyr::select(dplyr::one_of(c("TEMP", "MARKERS", "REF", "ALT")), dplyr::everything())

    x.alt <- dplyr::filter(x, TEMP == 1) %>%
      dplyr::arrange(MARKERS) %>%
      dplyr::select(-TEMP) %>%
      data.table::as.data.table(.) %>%
      data.table::melt.data.table(
        data = .,
        id.vars = c("MARKERS", "REF", "ALT"),
        variable.name = "TARGET_ID", variable.factor = FALSE,
        value.name = "ALLELE_ALT_DEPTH"
      ) %>%
      tibble::as_tibble(.) %>%
      dplyr::arrange(MARKERS, TARGET_ID)

    x %<>% dplyr::filter(TEMP == 2) %>%
      dplyr::arrange(MARKERS) %>%
      dplyr::select(-dplyr::one_of(c("TEMP", "REF", "ALT"))) %>%
      data.table::as.data.table(.) %>%
      data.table::melt.data.table(
        data = .,
        id.vars = "MARKERS",
        variable.name = "TARGET_ID", variable.factor = FALSE,
        value.name = "ALLELE_REF_DEPTH"
      ) %>%
      tibble::as_tibble(.) %>%
      dplyr::arrange(MARKERS, TARGET_ID) %>%
      dplyr::select(ALLELE_REF_DEPTH) %>%
      dplyr::bind_cols(x.alt) %>%
      dplyr::mutate_at(
        .tbl = .,
        .vars = c("ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH"),
        .funs = as.integer) %>%
      dplyr::select(
        dplyr::one_of(
          c("MARKERS", "REF", "ALT", "TARGET_ID", "ALLELE_REF_DEPTH",
            "ALLELE_ALT_DEPTH")))

    x.alt <- NULL
  }

  # 2 rows format also called binary--------------------------------------------
  if (dart.format == "2rows") {
    x %<>% dplyr::mutate(MISSING = is.na(ALLELE_REF_DEPTH) & is.na(ALLELE_ALT_DEPTH))

    if (write.tidy == "full") {
      x <- dplyr::bind_rows(
        dplyr::filter(x, MISSING) %>%
          dplyr::mutate(
            MISSING = NULL,
            ALLELE_REF_DEPTH = NULL,
            ALLELE_ALT_DEPTH = NULL,
            REF = NULL,
            ALT = NULL,
            GT_VCF = "./.",
            GT_BIN = NA_integer_,
            GT_VCF_NUC = "./.",
            GT = "000000"
          ),
        dplyr::filter(x, !MISSING) %>%
          dplyr::mutate(
            MISSING = NULL,
            A1 = dplyr::if_else(ALLELE_REF_DEPTH > 0, 0L, 1L),
            A2 = dplyr::if_else(ALLELE_ALT_DEPTH > 0, 1L, 0L),
            GT_VCF = stringi::stri_join(A1, A2, sep = "/"),
            GT_BIN = A1 + A2,
            A1 = dplyr::if_else(ALLELE_REF_DEPTH > 0, REF, ALT),
            A2 = dplyr::if_else(ALLELE_ALT_DEPTH > 0, ALT, REF),
            ALLELE_REF_DEPTH = NULL,
            ALLELE_ALT_DEPTH = NULL,
            REF = NULL,
            ALT = NULL,
            GT_VCF_NUC = stringi::stri_join(A1, A2, sep = "/"),
            A1 = stringi::stri_replace_all_regex(
              str = A1,
              pattern = c("^A$", "^C$", "^G$", "^T$"),
              replacement = c("001", "002", "003", "004"),
              vectorize_all = FALSE),
            A2 = stringi::stri_replace_all_regex(
              str = A2,
              pattern = c("^A$", "^C$", "^G$", "^T$"),
              replacement = c("001", "002", "003", "004"),
              vectorize_all = FALSE),
            GT = stringi::stri_join(A1, A2),
            A1 = NULL,
            A2 = NULL
          )
      )
    } else {
      x <- dplyr::bind_rows(
        dplyr::filter(x, MISSING) %>%
          dplyr::mutate(
            MISSING = NULL,
            REF = NULL,
            ALT = NULL,
            ALLELE_REF_DEPTH = NULL,
            ALLELE_ALT_DEPTH = NULL,
            GT_BIN = NA_integer_
          ),
        dplyr::filter(x, !MISSING) %>%
          dplyr::mutate(
            MISSING = NULL,
            REF = NULL,
            ALT = NULL,
            A1 = dplyr::if_else(ALLELE_REF_DEPTH > 0, 0L, 1L),
            A2 = dplyr::if_else(ALLELE_ALT_DEPTH > 0, 1L, 0L),
            ALLELE_REF_DEPTH = NULL,
            ALLELE_ALT_DEPTH = NULL,
            GT_BIN = A1 + A2,
            A1 = NULL,
            A2 = NULL
          )
      )
    }
  }#End 2rows

  # Count data------------------------------------------------------------------
  if (dart.format == "counts") {
    x<- dplyr::mutate(x, MISSING = ALLELE_REF_DEPTH == 0 & ALLELE_ALT_DEPTH == 0)
    if (write.tidy == "full") {
      x <- dplyr::bind_rows(
        dplyr::filter(x, MISSING) %>%
          dplyr::mutate(
            MISSING = NULL,
            REF = NULL,
            ALT = NULL,
            GT_VCF = "./.",
            GT_VCF_NUC = "./.",
            GT_BIN = NA_integer_,
            GT = "000000",
            ALLELE_REF_DEPTH = NA_integer_,
            ALLELE_ALT_DEPTH = NA_integer_,
            READ_DEPTH = NA_integer_
          ),
        dplyr::filter(x, !MISSING) %>%
          dplyr::mutate(
            MISSING = NULL,
            A1 = dplyr::if_else(ALLELE_REF_DEPTH > 0, 0L, 1L),
            A2 = dplyr::if_else(ALLELE_ALT_DEPTH > 0, 1L, 0L),
            GT_VCF = stringi::stri_join(A1, A2, sep = "/"),
            GT_BIN = A1 + A2,
            A1 = dplyr::if_else(ALLELE_REF_DEPTH > 0, REF, ALT),
            A2 = dplyr::if_else(ALLELE_ALT_DEPTH > 0, ALT, REF),
            REF = NULL,
            ALT = NULL,
            GT_VCF_NUC = stringi::stri_join(A1, A2, sep = "/"),
            A1 = stringi::stri_replace_all_regex(
              str = A1,
              pattern = c("^A$", "^C$", "^G$", "^T$"),
              replacement = c("001", "002", "003", "004"),
              vectorize_all = FALSE),
            A2 = stringi::stri_replace_all_regex(
              str = A2,
              pattern = c("^A$", "^C$", "^G$", "^T$"),
              replacement = c("001", "002", "003", "004"),
              vectorize_all = FALSE),
            GT = stringi::stri_join(A1, A2),
            A1 = NULL,
            A2 = NULL,
            READ_DEPTH = ALLELE_REF_DEPTH + ALLELE_ALT_DEPTH
          )
      )
    } else {
      # now faster:
      x <- dplyr::bind_rows(
        dplyr::filter(x, MISSING) %>%
          dplyr::mutate(
            MISSING = NULL,
            REF = NULL,
            ALT = NULL,
            GT_BIN = NA_integer_,
            ALLELE_REF_DEPTH = NA_integer_,
            ALLELE_ALT_DEPTH = NA_integer_,
            READ_DEPTH = NA_integer_
          ),
        dplyr::filter(x, !MISSING) %>%
          dplyr::mutate(
            MISSING = NULL,
            REF = NULL,
            ALT = NULL,
            A1 = dplyr::if_else(ALLELE_REF_DEPTH > 0, 0L, 1L),
            A2 = dplyr::if_else(ALLELE_ALT_DEPTH > 0, 1L, 0L),
            GT_BIN = A1 + A2,
            A1 = NULL,
            A2 = NULL,
            READ_DEPTH = ALLELE_REF_DEPTH + ALLELE_ALT_DEPTH
          )
      )
    }
  }#End counts data

  return(x)
}#End dart2gt

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
  gt.vcf = NULL,
  gt.vcf.nuc = NULL,
  gt = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE
) {
  # gt.vcf <- gt.vcf.nuc <- gt <- NULL

  # 1 row format----------------------------------------------------------------
  if (dart.format == "1row") {
    data %<>%
      dplyr::select(-dplyr::one_of(c("REF", "ALT", "MARKERS"))) %>%
      data.matrix(.) %>%
      magrittr::inset(. == 1L, 3L) %>%
      magrittr::inset(. == 2L, 1L) %>%
      magrittr::inset(. == 3L, 2L) %>%
      magrittr::set_rownames(x = ., value = variant.id) %>%
      as.matrix(.)
    dp <- genotypes.meta <- NULL
    source <- c("dart", "1row")
  }# End 1row

  # 2rows and counts -----------------------------------------------------------
  if (dart.format != "1row") { # 2rows and counts
    message("Generating genotypes...")
    if (dart.format == "2rows") source <- c("dart", "2rows")
    if (dart.format == "counts") source <- c("dart", "counts")
    parallel.core.temp <- 2 * floor(parallel.core / 2)
    geno.coding <- "alt.dos"
    dp <- NULL

    genotypes.meta <- dplyr::left_join(
      data,
      dplyr::distinct(data, MARKERS) %>%
        dplyr::mutate(SPLIT_VEC = split_vec_row(., 4, parallel.core = parallel.core.temp))
      , by = "MARKERS") %>%
      dplyr::group_split(SPLIT_VEC, keep = FALSE) %>%
      .radiator_parallel_mc(
        X = .,
        FUN = generate_geno,
        mc.cores = parallel.core.temp,
        source = source,
        gt.vcf = gt.vcf,
        gt.vcf.nuc = gt.vcf.nuc,
        gt = gt
      ) %>%
      dplyr::bind_rows(.) %>%
      tibble::add_column(.data = ., FILTERS = "whitelist", .before = 1) %>%
      dplyr::arrange(MARKERS, TARGET_ID)

    data <- genotypes.meta %>%
      dplyr::select(VARIANT_ID, TARGET_ID, GT_BIN) %>%
      data.table::as.data.table(.) %>%
      data.table::dcast.data.table(
        data = .,
        formula = VARIANT_ID ~ TARGET_ID,
        value.var = "GT_BIN"
      ) %>%
      tibble::as_tibble(.) %>%
      dplyr::arrange(VARIANT_ID) %>%
      tibble::column_to_rownames(.data = ., var = "VARIANT_ID")
  }# End counts and 2rows

  # Generate GDS ---------------------------------------------------------------
  data <- radiator_gds(
    genotypes.df = data,
    geno.coding = geno.coding,
    strata = strata,
    biallelic = TRUE,
    markers.meta = markers.meta,
    genotypes.meta = genotypes.meta,
    dp = dp,
    filename = filename,
    source = source,
    open = TRUE,
    verbose = verbose)
  if (verbose) message("done!")
  return(data)
}# End dart2gds

# generate_geno----------------------------------------------------------------
#' @title generate_geno
#' @description Generate the genotypes and prep for GDS
#' @rdname generate_geno
#' @keywords internal
#' @export
generate_geno <- function(
  x,
  source,
  gt = FALSE,
  gt.vcf.nuc = FALSE,
  gt.vcf = FALSE
) {
  # res <- list()
  message("Generating genotypes...")

  #2-rows
  if ("2rows" %in% source) {
    # x <- genotypes.meta[[1]]
    res <- dplyr::filter(x, !is.na(REF)) %>%
      dplyr::arrange(MARKERS) %>%
      data.table::as.data.table(.) %>%
      data.table::melt.data.table(
        data = .,
        id.vars = c("VARIANT_ID", "MARKERS", "REF", "ALT"),
        variable.name = "TARGET_ID",
        value.name = "A2",
        variable.factor = FALSE) %>%
      tibble::as_tibble(.) %>%
      dplyr::bind_cols(
        dplyr::filter(x, is.na(REF)) %>%
          dplyr::select(-REF, -ALT, -VARIANT_ID) %>%
          data.table::as.data.table(.) %>%
          data.table::melt.data.table(
            data = .,
            id.vars = "MARKERS",
            variable.name = "TARGET_ID",
            value.name = "A1",
            variable.factor = FALSE) %>%
          tibble::as_tibble(.)
      ) %>%
      dplyr::mutate_at(.tbl = ., .vars = c("A1", "A2"), .funs = as.integer)
    x <- NULL

    if (!identical(res$MARKERS, res$MARKERS1)) {
      rlang::abort("Contact author, DArT tiding problem")
    } else {
      res %<>% dplyr::select(-MARKERS1)
    }

    if (!identical(res$TARGET_ID, res$TARGET_ID1)) {
      rlang::abort("Contact author, DArT tiding problem")
    } else {
      res %<>% dplyr::select(-TARGET_ID1)
    }

    switch_allele_count <- function(x) {
      # x <- as.integer(dplyr::recode(.x = as.character(x), "0" = "1", "1" = "0"))
      # case_when is much faster than recode...
      x <- dplyr::case_when(
        x == 0 ~ 1,
        x == 1 ~ 0
      )
      # here we want count of alternate allele instead...
    }

    res %<>% dplyr::mutate(GT_BIN = switch_allele_count(A1) + A2, A1 = NULL, A2 = NULL)

    # Counts for allele calibration
    switch <- dplyr::select(res, MARKERS, GT_BIN) %>%
      dplyr::filter(!is.na(GT_BIN)) %>%
      dplyr::count(GT_BIN, MARKERS) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::summarise(
        REF_COUNT = sum((2 * n[GT_BIN == 0]), n[GT_BIN == 1], na.rm = TRUE),
        ALT_COUNT = sum((2 * n[GT_BIN == 2]), n[GT_BIN == 1], na.rm = TRUE)
      ) %>%
      dplyr::ungroup(.) %>%
      dplyr::filter(dplyr::if_else(REF_COUNT < ALT_COUNT, TRUE, FALSE)) %$%
      MARKERS

    n.switch <- length(switch)
    if (n.switch > 0) {
      message("Calibration REF/ALT based on counts of alleles: ", n.switch)
      res <- dplyr::filter(res, !MARKERS %in% switch) %>%
        dplyr::bind_rows(
          dplyr::filter(res, MARKERS %in% switch) %>%
            dplyr::rename(ALT = REF, REF = ALT)
        )
    }
    switch <- NULL
  }#2rows genotypes

  # counts
  if ("counts" %in% source) {
    res <- dplyr::filter(x, !is.na(REF)) %>%
      dplyr::arrange(MARKERS) %>%
      data.table::as.data.table(.) %>%
      data.table::melt.data.table(
        data = .,
        id.vars = c("VARIANT_ID", "MARKERS", "REF", "ALT"),
        variable.name = "TARGET_ID",
        value.name = "ALLELE_ALT_DEPTH",
        variable.factor = FALSE) %>%
      tibble::as_tibble(.) %>%
      dplyr::bind_cols(
        dplyr::filter(x, is.na(REF)) %>%
          dplyr::select(-REF, -ALT, -VARIANT_ID) %>%
          data.table::as.data.table(.) %>%
          data.table::melt.data.table(
            data = .,
            id.vars = "MARKERS",
            variable.name = "TARGET_ID",
            value.name = "ALLELE_REF_DEPTH",
            variable.factor = FALSE) %>%
          tibble::as_tibble(.)
      )
    x <- NULL

    if (!identical(res$MARKERS, res$MARKERS1)) {
      rlang::abort("Contact author, DArT tiding problem")
    } else {
      res %<>% dplyr::select(-MARKERS1)
    }

    if (!identical(res$TARGET_ID, res$TARGET_ID1)) {
      rlang::abort("Contact author, DArT tiding problem")
    } else {
      res %<>% dplyr::select(-TARGET_ID1)
    }

    res %<>%
      dplyr::mutate_at(
        .tbl = .,
        .vars = c("ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH"),
        .funs = as.numeric
      ) %>%
      dplyr::mutate(READ_DEPTH = ALLELE_REF_DEPTH + ALLELE_ALT_DEPTH)

    # Coverage for allele calibration
    want <- c("ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH")
    switch <- dplyr::group_by(res, MARKERS) %>%
      dplyr::summarise_at(.tbl = ., .vars = want, .funs = sum, na.rm = TRUE) %>%
      dplyr::mutate_at(.tbl = ., .vars = want, .funs = round, digits = 0) %>%
      dplyr::mutate_at(.tbl = ., .vars = want, .funs = as.integer) %>%
      dplyr::ungroup(.) %>%
      dplyr::filter(dplyr::if_else(ALLELE_ALT_DEPTH > ALLELE_REF_DEPTH, TRUE, FALSE)) %$%
      MARKERS
    n.switch <- length(switch)
    if (n.switch > 0) {
      message("Calibration REF/ALT based on read depth of alleles: ", n.switch)
      res <- dplyr::filter(res, !MARKERS %in% switch) %>%
        dplyr::bind_rows(
          dplyr::filter(res, MARKERS %in% switch) %>%
            dplyr::rename(
              ALT = REF,
              REF = ALT,
              ALLELE_REF_DEPTH = ALLELE_ALT_DEPTH,
              ALLELE_ALT_DEPTH = ALLELE_REF_DEPTH
            )
        )
    }
    switch <- NULL

    # GT_BIN
    res %<>%
      dplyr::mutate(
        GT_BIN = dplyr::case_when(
          ALLELE_REF_DEPTH > 0 & ALLELE_ALT_DEPTH == 0 ~ 0,
          ALLELE_REF_DEPTH > 0 & ALLELE_ALT_DEPTH > 0 ~ 1,
          ALLELE_REF_DEPTH == 0 & ALLELE_ALT_DEPTH > 0 ~ 2
        )
      ) %>%
      dplyr::mutate_at(
        .tbl = .,
        .vars = c("READ_DEPTH", "ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH"),
        .funs = replace_by_na, what = 0
      )
  }# Counts

  if (gt.vcf) {
    res %<>%
      dplyr::mutate(
        GT_VCF = dplyr::case_when(
          GT_BIN == 0 ~ "0/0", GT_BIN == 1 ~ "0/1", GT_BIN == 2 ~ "1/1",
          is.na(GT_BIN) ~ "./.")
      )
  }

  if (gt.vcf.nuc) {
    res %<>%
      dplyr::mutate(
        GT_VCF_NUC = dplyr::case_when(
          GT_BIN == 0 ~ stringi::stri_join(REF, REF, sep = "/"),
          GT_BIN == 1 ~ stringi::stri_join(REF, ALT, sep = "/"),
          GT_BIN == 2 ~ stringi::stri_join(ALT, ALT, sep = "/"),
          is.na(GT_BIN) ~ "./.")
      )
  }

  if (gt) {
    res %<>%
      dplyr::mutate(
        GT = stringi::stri_replace_all_fixed(
          str = GT_VCF_NUC,
          pattern = c("A", "C", "G", "T", "/", ".."),
          replacement = c("001", "002", "003", "004", "", "000000"),
          vectorize_all = FALSE)
      )
  }

  return(res)
}# End generate_geno


# tidy_silico_dart -------------------------------------------------------------
# Import, filter and transform a silico dart output file to different formats

#' @name tidy_silico_dart

#' @title Tidy silico \href{http://www.diversityarrays.com}{DArT} output file.

#' @description Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users. The function generate a tidy dataset using
#' \href{http://www.diversityarrays.com}{DArT} file.

#' @param data Silico DArT output file.
#' The function can import \code{.csv} or \code{.tsv} files.

#' @param strata A tab delimited file or object with 3 columns.
#' Columns header is:
#' \code{TARGET_ID}, \code{INDIVIDUALS} and \code{STRATA}.
#' Note: the column \code{STRATA} refers to any grouping of individuals.
#' You need to make sure that
#' the column \code{TARGET_ID} match the id used by DArT.
#' The column \code{INDIVIDUALS} and \code{STRATA} will be kept in the tidy data.
#' Only individuals in the strata file are kept in the tidy, i.e. that the strata
#' is also used as a whitelist of individuals/strata.
#' Silico DArT data is currently used to detect sex markers, so the \code{STRATA}
#' column should be filed with sex information: \code{M} or \code{F}.

#' @param verbose (optional, logical) When verbose = TRUE the function is a
#' little more chatty during execution.
#' Default: \code{verbose = FALSE}.

#' @inheritParams tidy_genomic_data


#' @return A tidy dataframe with several columns depending on DArT file:
#' #TODO
#'
#'
#' Written in the working directory:
#' #TODO

#' @export
#' @rdname tidy_silico_dart
#' @importFrom dplyr group_by select rename filter mutate summarise distinct n_distinct arrange left_join semi_join anti_join inner_join full_join tally bind_rows
#' @importFrom parallel detectCores
#' @importFrom stringi stri_replace_all_fixed stri_join stri_sub stri_replace_na stri_pad_left
#' @importFrom purrr discard
#' @importFrom data.table as.data.table dcast.data.table melt.data.table
#' @importFrom readr read_tsv write_tsv read_lines read_table
#' @importFrom tibble as_data_frame data_frame
#' @importFrom tidyr spread gather unite separate

#' @examples
#' \dontrun{
#' clownfish.dart.tidy <- radiator::tidy_silico_dart(
#' data = "clownfish.dart.tsv",
#' strata = "clownfish.strata.tsv",
#' verbose = TRUE)
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

tidy_silico_dart <- function(
  data,
  strata = NULL,
  filename = NULL,
  verbose = TRUE
) {

  # Test
  # data = "Report_DAci18-3647_1_moreOrders_SilicoDArT_2.csv"
  # strata = "sturgeon_sex_strata_no_u.tsv"
  # filename = NULL
  # verbose = TRUE

  opt.change <- getOption("width")
  options(width = 70)
  # for timing
  timing <- proc.time()
  message("Importing silico DArT data...")

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("Input file missing")
  if (missing(strata)) rlang::abort("strata file missing")

  # Check that DArT file as good target id written -----------------------------
  target.id <- radiator::extract_dart_target_id(data, write = FALSE)
  n.ind.dart <- nrow(target.id)
  if (verbose) message("Number of individuals in DArT file: ", n.ind.dart)
  if (nrow(target.id) != length(unique(target.id$TARGET_ID))) {
    rlang::abort("\nnon unique target id are used in the DArT file...
                 What you want are different target ids at the end of the row
                 that contains AlleleID, AlleleSequence.
                 Edit manually before trying again
                 If you're still encountering problem, email author for help")
  }

  # Date & Time and filenames --------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")

  if (is.null(filename)) {
    filename <- stringi::stri_join("radiator_tidy_silico_dart_", file.date, ".tsv")
    strata.filename <- stringi::stri_join("radiator_tidy_silico_dart_strata_", file.date, ".tsv")
  } else {
    strata.filename <- stringi::stri_join(filename, "_strata_", file.date, ".tsv")
    filename <- stringi::stri_join(filename, "_", file.date, ".tsv")
  }

  # Check DArT format file -----------------------------------------------------
  dart.check <- check_dart(data)
  if (!dart.check$data.type %in% c("dart", "silico.dart")) {
    rlang::abort("\nContact author to show your silico DArT data, problems during import")
  } else {
    skip.number <- dart.check$skip.number
  }

  # Strata file ------------------------------------------------------------------
  if (verbose) message("\nMaking DArT data population-wise...")
  strata.df <- radiator::read_strata(
    strata = strata,
    pop.levels = NULL, pop.labels = NULL,
    pop.select = NULL, blacklist.id = NULL,
    keep.two = FALSE, verbose = verbose)
  pop.levels <- strata.df$pop.levels
  pop.labels <- strata.df$pop.labels
  # pop.select <- strata.df$pop.select
  # blacklist.id <- strata.df$blacklist.id
  strata.df <- strata.df$strata
  n.ind.strata <- nrow(strata.df)

  # Check that TARGET_ID in strata match TARGET_ID in the DArT file ------------
  if (n.ind.dart != n.ind.strata) {
    strata.id.check <- strata.df %>%
      dplyr::mutate(IN_DART = stringi::stri_trans_toupper(strata.df$TARGET_ID)
                    %in% stringi::stri_trans_toupper(target.id$TARGET_ID))
    strata.id.pass <- !FALSE %in% (unique(strata.id.check$IN_DART))
    if (!strata.id.pass) {
      problem.filename <- stringi::stri_join("radiator_tidy_dart_strata_problem_", file.date, ".tsv")
      readr::write_tsv(
        x = strata.id.check,
        path = problem.filename)
      rlang::abort("\nSome of the samples in the strata are not found in the DArT file.
                   For more info: ", problem.filename)
    }
    message("\nCaution: you've chosen to tidy a subsample of your DArT file.
            DArT statistics generated for all samples might not apply...\n")
    strata.id.check <- NULL
  } else {
    if (!identical(sort(target.id$TARGET_ID), sort(strata.df$TARGET_ID))) {
      rlang::abort("\nThe DArT and strata files don't have the same TARGET_IDs")
    }
  }
  target.id <- NULL

  # need to check for duplicate names... yes happening all the time
  duplicate.id.strata <- length(strata.df$INDIVIDUALS) - dplyr::n_distinct(strata.df$INDIVIDUALS)

  if (duplicate.id.strata > 0) {
    message("Duplicated individuals names found in the strata.\n   number of duplicate names = ", duplicate.id.strata, "\n")
    rlang::abort("\nFix the strata with unique names and\nverify the DArT file for the same issue, adjust accordingly...")
  }

  # Import data ---------------------------------------------------------------
  if (stringi::stri_detect_fixed(
    str = stringi::stri_sub(str = data, from = -4, to = -1),
    pattern = ".csv")) {
    csv <- TRUE
  } else {
    csv <- FALSE
  }

  if (csv) {
    dart.col.type <- readr::read_csv(
      file = data,
      skip = skip.number, n_max = 1,
      na = "-",
      col_names = FALSE,
      col_types = readr::cols(.default = readr::col_character()))
  } else {
    dart.col.type <- readr::read_tsv(
      file = data,
      skip = skip.number, n_max = 1,
      na = "-",
      col_names = FALSE,
      col_types = readr::cols(.default = readr::col_character()))
  }

  want <- tibble::data_frame(
    INFO = c("CLONEID", "ALLELESEQUENCE", "SEQUENCE", "TRIMMEDSEQUENCE"),
    COL_TYPE = c("c", "c", "c", "c")) %>%
    dplyr::bind_rows(
      dplyr::select(strata.df, INFO = TARGET_ID) %>%
        dplyr::mutate(
          COL_TYPE = rep("c", n()),
          INFO = stringi::stri_trans_toupper(INFO)))

  dart.col.type <- dart.col.type %>%
    tidyr::gather(data = .,key = DELETE, value = INFO) %>%
    dplyr::select(-DELETE) %>%
    dplyr::mutate(
      INFO = stringi::stri_trans_toupper(INFO),
      INFO = stringi::stri_replace_all_fixed(INFO, pattern = " ", replacement = "", vectorize_all = FALSE)
    ) %>%
    dplyr::left_join(want, by = "INFO") %>%
    dplyr::mutate(COL_TYPE = stringi::stri_replace_na(str = COL_TYPE, replacement = "_")) %>%
    dplyr::select(COL_TYPE) %>%
    purrr::flatten_chr(.) %>% stringi::stri_join(collapse = "")
  want <- NULL

  if (csv) {
    input <- suppressMessages(suppressWarnings(
      readr::read_csv(
        file = data,
        skip = skip.number,
        na = c("-", " ", "", "NA"),
        col_names = TRUE,
        col_types = dart.col.type)
    ))
  } else {
    input <- suppressMessages(suppressWarnings(
      readr::read_tsv(
        file = data,
        skip = skip.number,
        na = c("-", " ", "", "NA"),
        col_names = TRUE,
        col_types = dart.col.type)
    ))
  }
  dart.col.type <- NULL
  colnames(input) <- stringi::stri_trans_toupper(colnames(input))
  colnames(input) = stringi::stri_replace_all_fixed(
    str = colnames(input), pattern = " ", replacement = "", vectorize_all = FALSE)

  check.data <- c("ALLELESEQUENCE", "SEQUENCE", "TRIMMEDSEQUENCE")
  check.data <- purrr::keep(.x = colnames(input), .p = colnames(input) %in% check.data)

  if (length(check.data > 1)) {
    keeper <- sample(check.data, size = 1)
    keeper <- c("CLONEID", keeper, strata.df$TARGET_ID)

    input <- dplyr::select(input, dplyr::one_of(keeper))
    colnames(input) <- stringi::stri_replace_all_fixed(
      str = colnames(input),
      pattern = c("ALLELESEQUENCE", "TRIMMEDSEQUENCE"),
      replacement = c("SEQUENCE", "SEQUENCE"),
      vectorize_all = FALSE)
  }

  input <- data.table::as.data.table(input) %>%
    data.table::melt.data.table(
      data = .,
      id.vars = c("CLONEID", "SEQUENCE"),
      variable.name = "TARGET_ID",
      variable.factor = FALSE,
      value.name = "VALUE"
    ) %>%
    tibble::as_data_frame(.)

  # Strata file to include populations ----------------------------------------
  # To make sure target ids match
  input <- dplyr::left_join(input, strata.df, by = "TARGET_ID") %>%
    dplyr::select(-TARGET_ID)
  strata.df <- NULL

  if (tibble::has_name(input, "STRATA")) {
    input <- dplyr::rename(input, POP_ID = STRATA)
  }


  input <- dplyr::mutate(input, VALUE = as.integer(VALUE))

  # write tidy
  input <- suppressWarnings(
    dplyr::ungroup(input) %>%
      dplyr::arrange(CLONEID, POP_ID, INDIVIDUALS)) %>%
    readr::write_tsv(x = ., path = filename)

  # clean...
  input$INDIVIDUALS <- radiator::clean_ind_names(input$INDIVIDUALS)
  input$POP_ID <- radiator::clean_pop_names(input$POP_ID)

  # Generate a new strata file -------------------------------------------------
  strata <- dplyr::distinct(input, INDIVIDUALS, POP_ID) %>%
    dplyr::rename(STRATA = POP_ID) %>%
    readr::write_tsv(x = ., path = strata.filename)
  # Results --------------------------------------------------------------------
  if (verbose) {
    n.pop <- dplyr::n_distinct(strata$STRATA)
    n.ind <- dplyr::n_distinct(strata$INDIVIDUALS)
    n.snp <- dplyr::n_distinct(input$CLONEID)

    message("Number of SNPs: ", n.snp)
    message("Number of populations: ", n.pop)
    message("Number of individuals: ", n.ind)

    strata.stats <- strata %>%
      dplyr::group_by(STRATA) %>%
      dplyr::tally(.) %>%
      dplyr::mutate(POP_IND = stringi::stri_join(STRATA, n, sep = " = "))

    message("\nNumber of ind/pop:\n", stringi::stri_join(strata.stats$POP_IND, collapse ="\n"))
    timing <- proc.time() - timing
    message("\nComputation time: ", round(timing[[3]]), " sec")
  }
  options(width = opt.change)
  return(input)
}#End tidy_silico_dart
