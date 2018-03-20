## Merge dart

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
#' @inheritParams tidy_dart

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
#' \item AVG_COUND_REF and AVG_COUND_SNP: the coverage for the reference and
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
  pop.select = NULL,
  monomorphic.out = TRUE,
  common.markers = TRUE,
  keep.rad = FALSE,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1
) {
  cat("#######################################################################\n")
  cat("########################## radiator::merge_dart #######################\n")
  cat("#######################################################################\n")
  timing <- proc.time()

  # manage missing arguments -----------------------------------------------------
  if (missing(dart1)) stop("dart1 file missing")
  if (missing(dart2)) stop("dart2 file missing")
  if (missing(strata1)) stop("strata1 file missing")
  if (missing(strata2)) stop("strata2 file missing")

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

  # the DArT data
  if (data.type == "dart") {
    message("Importing and tidying dart1...")
    input <- radiator::tidy_dart(data = dart1, strata = strata1, filename = "temp.tidy.dart1", verbose = FALSE)
    message("Importing and tidying dart2...")
    dart2 <- radiator::tidy_dart(data = dart2, strata = strata2, filename = "temp.tidy.dart2", verbose = FALSE)
    if (!keep.rad) {
      message("Removing temporary tidy DArT files...")
      file.remove("temp.tidy.dart1.rad")
      file.remove("temp.tidy.dart2.rad")
    }
  }

  # The filtered DArT data
  if (data.type == "fst.file") {
    message("Importing the filtered dart1...")
    input <- radiator::read_rad(dart1)
    message("Importing the filtered dart2...")
    dart2 <- radiator::read_rad(dart2)
  }

  if (!is.null(pop.select)) {
    input <- suppressWarnings(dplyr::filter(input, POP_ID %in% pop.select))
    dart2 <- suppressWarnings(dplyr::filter(dart2, POP_ID %in% pop.select))
  }


  # merging DArT tidy data -----------------------------------------------------
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


  # cleaning up non-immortalized markers
  message("Removing markers starting with 1000 (non-immortalized DArT markers id")
  markers.before <- dplyr::n_distinct(input$MARKERS)

  input <- input %>%
    dplyr::mutate(DISCARD = stringi::stri_detect_regex(str = LOCUS, pattern = "^1000")) %>%
    dplyr::filter(!DISCARD) %>%
    dplyr::select(-DISCARD)

  markers.after <- dplyr::n_distinct(input$MARKERS)
  message("    Markers removed: ", markers.before - markers.after)

  # Averaging across markers the call rate and other DArT markers metadata statistics
  if (tibble::has_name(input, "CALL_RATE")) {
    message("Averaging across markers the call rate")
    input <- input %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(CALL_RATE = mean(CALL_RATE)) %>%
      dplyr::ungroup(.)
  }

  if (tibble::has_name(input, "REP_AVG")) {
    message("Averaging across markers the REP_AVG")
    input <- input %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(REP_AVG = mean(REP_AVG)) %>%
      dplyr::ungroup(.)
  }

  message("Adjusting REF/ALT alleles...")
  input <- radiator::change_alleles(
    data = input,
    parallel.core = parallel.core,
    verbose = TRUE)$input

  if (tibble::has_name(input, "POLYMORPHIC")) {
    input <- dplyr::select(input, -POLYMORPHIC)
  }

  if (tibble::has_name(input, "AVG_COUNT_REF")) {
    message("Averaging across markers the coverage for the REF allele")
    input <- input %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(AVG_COUNT_REF = mean(AVG_COUNT_REF)) %>%
      dplyr::ungroup(.)
  }

  if (tibble::has_name(input, "AVG_COUNT_ALT")) {
    message("Averaging across markers the coverage for the ALT allele")
    input <- input %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(AVG_COUNT_REF = mean(AVG_COUNT_REF)) %>%
      dplyr::ungroup(.)
  }

  if (monomorphic.out) {
    input <- radiator::discard_monomorphic_markers(data = input, verbose = TRUE)$input
  }

  if (common.markers) {
    input <- radiator::keep_common_markers(data = input, plot = TRUE, verbose = TRUE)$input
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
  return(res = list(merged.dart = input, strata = strata))
}#End merge_dart
