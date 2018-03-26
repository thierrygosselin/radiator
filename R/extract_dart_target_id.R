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
  if (missing(data)) stop("Input file missing")

  # Check DArT format file -----------------------------------------------------
  data.type <- readChar(con = data, nchars = 16L, useBytes = TRUE)
  dart.with.header <- TRUE %in% (stringi::stri_detect_fixed(str = data.type, pattern = c("*\t", "*,")))
  if (dart.with.header) {
    temp.file <- suppressWarnings(suppressMessages(readr::read_table(file = data, n_max = 20, col_names = "HEADER")))
    skip.number <- which(stringi::stri_detect_fixed(str = temp.file$HEADER,
                                                    pattern = "AlleleID") |
                           stringi::stri_detect_fixed(str = temp.file$HEADER,
                                                      pattern = "CloneID")) - 1
    data.type <- readr::read_lines(file = data, skip = skip.number, n_max = skip.number + 1)[1] %>%
      stringi::stri_sub(str = ., from = 1, to = 16)
  } else {
    skip.number <- 0
  }
  temp.file <- NULL
  dart.clone.id <- stringi::stri_detect_fixed(str = data.type, pattern = "CloneID")
  dart.allele.id <- stringi::stri_detect_fixed(str = data.type, pattern = "AlleleID")

  if (dart.clone.id || dart.allele.id) {
    data.type <- "dart"
  } else {
    stop("Contact author to show your DArT data, problem duting import")
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
    "PIC", "AVGREADDEPTH", "STDEVREADDEPTH", "QPMR", "REPRODUCIBILITY")

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
      TARGET_ID, pattern = " ", replacement = "", vectorize_all = FALSE))

  if (write) readr::write_tsv(x = dart.target.id, path = "dart.target.id.tsv")

  # Check that DArT file as good target id written -----------------------------
  if (nrow(dart.target.id) != length(unique(dart.target.id$TARGET_ID))) {
    message("non unique target id are used in the DArT file...
What you want are different target ids at the end of the row that contains AlleleID, AlleleSequence, etc
Edit manually the DArT file before trying the functions: tidy_dart and filter_dart
If you're still encountering problem, email author for help")
  }
  return(dart.target.id)
}#End extract_dart_target_id
