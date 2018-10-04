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
  if (missing(data)) stop("Input file missing")
  if (missing(strata)) stop("strata file missing")

  # Check that DArT file as good target id written -----------------------------
  target.id <- radiator::extract_dart_target_id(data, write = FALSE)
  n.ind.dart <- nrow(target.id)
  if (verbose) message("Number of individuals in DArT file: ", n.ind.dart)
  if (nrow(target.id) != length(unique(target.id$TARGET_ID))) {
    stop("\nnon unique target id are used in the DArT file...
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
    stop("\nContact author to show your silico DArT data, problems during import")
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
      stop("\nSome of the samples in the strata are not found in the DArT file.
For more info: ", problem.filename)
    }
    message("\nCaution: you've chosen to tidy a subsample of your DArT file.
DArT statistics generated for all samples might not apply...\n")
    strata.id.check <- NULL
  } else {
    if (!identical(sort(target.id$TARGET_ID), sort(strata.df$TARGET_ID))) {
      stop("\nThe DArT and strata files don't have the same TARGET_IDs")
    }
  }
  target.id <- NULL

  # need to check for duplicate names... yes happening all the time
  duplicate.id.strata <- length(strata.df$INDIVIDUALS) - dplyr::n_distinct(strata.df$INDIVIDUALS)

  if (duplicate.id.strata > 0) {
    message("Duplicated individuals names found in the strata.\n   number of duplicate names = ", duplicate.id.strata, "\n")
    stop("\nFix the strata with unique names and\nverify the DArT file for the same issue, adjust accordingly...")
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
