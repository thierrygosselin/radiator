## Merge vcf

#' @title Merge VCF files
#' @description This function allows to merge 2 VCF files.

#' @param vcf1 First VCF file.
#' @param strata1 strata file for vcf1.
#' @param vcf2 Second VCF file.
#' @param strata2 strata file for vcf2.
#' @param filename Name of the merged VCF file.
#' With the default, the function gives a filename based on date and time.
#' Default: \code{filename = NULL}.
#' @inheritParams tidy_genomic_data

#' @importFrom stringi stri_replace_all_fixed stri_replace_na stri_join stri_count_fixed
#' @importFrom tibble as_data_frame data_frame add_column add_row
#' @importFrom dplyr select rename n_distinct distinct mutate summarise group_by ungroup arrange left_join full_join semi_join anti_join bind_rows bind_cols if_else
#' @importFrom readr write_tsv read_tsv
#' @importFrom tidyr separate gather
#' @importFrom parallel detectCores
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot stat_summary labs theme element_blank element_text geom_jitter scale_colour_manual scale_y_reverse theme_light geom_bar facet_grid stat_smooth
#' @importFrom purrr flatten_chr map_df flatten_dbl

#' @return The function returns in the global environment a tidy dataset with
#' the merged VCF files and the merged VCF in the working directory.

#' @examples
#' \dontrun{
#' # The simplest way to run the function:
#' sum <- radiator::merge_vcf(
#' vcf1 = "batch_1.vcf", strata1 = "strata1_brook_charr.tsv",
#' vcf1 = "batch_2.vcf", strata2 = "strata2_brook_charr.tsv",
#' pop.select = c("QC", "ON", "NE"),
#' maf.thresholds = c(0.002, 0.001),
#' maf.pop.num.threshold = 1,
#' maf.approach = "SNP",maf.operator = "OR",
#' filename = "my_new_VCF.vcf"
#' }


#' @rdname merge_vcf
#' @export
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

merge_vcf <- function(
  vcf1, strata1,
  vcf2, strata2,
  monomorphic.out = TRUE,
  common.markers = TRUE,
  pop.levels = NULL,
  pop.labels = NULL,
  pop.select = NULL,
  blacklist.id = NULL,
  whitelist.markers = NULL,
  maf.thresholds = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1
) {
  cat("#######################################################################\n")
  cat("########################### radiator::merge_vcf #########################\n")
  cat("#######################################################################\n")
  timing <- proc.time()

  # manage missing arguments -----------------------------------------------------
  if (missing(vcf1)) stop("vcf1 file missing")
  if (missing(vcf2)) stop("vcf2 file missing")
  if (missing(strata1)) stop("strata1 file missing")
  if (missing(strata2)) stop("strata2 file missing")


  if (!is.null(pop.levels) & is.null(pop.labels)) {
    pop.levels <- stringi::stri_replace_all_fixed(
      pop.levels, pattern = " ", replacement = "_", vectorize_all = FALSE)
    pop.labels <- pop.levels
  }
  if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")
  if (!is.null(pop.labels)) {
    if (length(pop.labels) != length(pop.levels)) stop("pop.labels and pop.levels must have the same length (number of groups)")
    pop.labels <- stringi::stri_replace_all_fixed(pop.labels, pattern = " ", replacement = "_", vectorize_all = FALSE)
  }

  # Filename -------------------------------------------------------------------
  # Get date and time to have unique filenaming
  if (is.null(filename)) {
    file.date <- stringi::stri_replace_all_fixed(
      Sys.time(),
      pattern = " EDT",
      replacement = "",
      vectorize_all = FALSE
    )
    file.date <- stringi::stri_replace_all_fixed(
      file.date,
      pattern = c("-", " ", ":"),
      replacement = c("", "@", ""),
      vectorize_all = FALSE
    )
    file.date <- stringi::stri_sub(file.date, from = 1, to = 13)

    filename <- stringi::stri_join("radiator_merged_vcf_", file.date)
  }

  # import data ----------------------------------------------------------------
  message("Importing and tidying the vcf1...")
  input <- suppressMessages(radiator::tidy_genomic_data(
    data = vcf1,
    strata = strata1,
    vcf.metadata = FALSE,
    blacklist.id = blacklist.id,
    whitelist.markers = whitelist.markers,
    monomorphic.out = FALSE,
    common.markers = FALSE,
    filename = NULL,
    verbose = FALSE))

  message("Importing and tidying the vcf2...")
  # Also Using pop.levels and pop.labels info if present
  input <- suppressWarnings(
    dplyr::bind_rows(
      input,
      suppressMessages(
        radiator::tidy_genomic_data(
          data = vcf2,
          strata = strata2,
          vcf.metadata = FALSE,
          blacklist.id = blacklist.id,
          whitelist.markers = whitelist.markers,
          monomorphic.out = FALSE,
          common.markers = FALSE,
          filename = NULL,
          verbose = FALSE))) %>%
      radiator::change_pop_names(
        data = .,
        pop.levels = pop.levels, pop.labels = pop.labels)
  )

  # Pop select
  if (!is.null(pop.select)) {
    message(stringi::stri_join(length(pop.select), "population(s) selected", sep = " "))
    input <- suppressWarnings(input %>% dplyr::filter(POP_ID %in% pop.select))
    if (is.factor(input$POP_ID)) input$POP_ID <- droplevels(input$POP_ID)
  }

  message("Adjusting REF/ALT alleles...")
  input <- radiator::change_alleles(
    data = input,
    parallel.core = parallel.core,
    verbose = TRUE)$input

  if (monomorphic.out) {
    input <- radiator::discard_monomorphic_markers(data = input, verbose = TRUE)$input
  }

  if (common.markers) {
    input <- radiator::keep_common_markers(data = input, verbose = TRUE)$input
  }

  if (!is.null(maf.thresholds)) {
  input <- filter_maf(
  data = input,
  interactive.filter = FALSE,
  maf.thresholds = maf.thresholds,
  parallel.core = parallel.core,
  verbose = FALSE)$tidy.filtered.maf
  }

  # Write VCF in the working directory
  radiator::write_vcf(data = input, pop.info = FALSE, filename = filename)

  # results --------------------------------------------------------------------
  message("Merged VCF in the working directory: ", filename, ".vcf")
  timing <- proc.time() - timing
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("############################## completed ##############################\n")
  return(input)
}#End merge_vcf
