# write_genepopedit-----------------------------------------------------------------
#' @name write_genepopedit

#' @title Write a genepopedit flatten object

#' @description Write a genepopedit object from a tidy data frame or GDS file/object.

#' @inheritParams radiator_common_arguments

#' @return A genepopedit object in the global environment.

#' @export
#' @rdname write_genepopedit
#' @references Stanley RRE, Jeffery NW, Wringe BF, DiBacco C, Bradbury IR (2017)
#' genepopedit: a simple and flexible tool for manipulating multilocus molecular
#' data in R. Molecular Ecology Resources, 17, 12-18.
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_genepopedit <- function(data) {
  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("Input file missing")

  # File type detection----------------------------------------------------------
  data.type <- radiator::detect_genomic_format(data)

  # Import data ---------------------------------------------------------------

  if (data.type %in% c("SeqVarGDSClass", "gds.file")) {
    if (data.type == "gds.file") {
      data <- radiator::read_rad(data, verbose = FALSE)
    }
    data <- gds2tidy(gds = data,
                     pop.id = FALSE,
                     calibrate.alleles = FALSE,
                     parallel.core = parallel::detectCores() - 1) %>%
      dplyr::select(INDIVIDUALS, STRATA, MARKERS, GT_BIN, REF, ALT)
    data.type <- "tbl_df"
  } else {
    if (is.vector(data)) {
      data <- radiator::tidy_wide(data = data, import.metadata = TRUE)
    }
  }

  if (!rlang::has_name(data, "GT")) {
    data <- calibrate_alleles(
      data = data,
      verbose = FALSE,
      gt = TRUE,
      gt.vcf.nuc = FALSE,
      gt.vcf = FALSE,
      gt.bin = FALSE
    ) %$% input
  }

  if (!rlang::has_name(data, "STRATA") && rlang::has_name(data, "POP_ID")) {
    data %<>% dplyr::rename(STRATA = POP_ID)
  }

  # generate genepopedit flatten object
  data %<>%
    dplyr::mutate(SampleNum = INDIVIDUALS) %>%
    dplyr::select(SampleID = INDIVIDUALS, Population = STRATA, SampleNum, MARKERS, GT) %>%
    data.table::as.data.table(.) %>%
    data.table::dcast.data.table(
      data = .,
      formula = SampleID + Population + SampleNum ~ MARKERS,
      value.var = "GT"
    ) %>%
    tibble::as_tibble(.)

  return(data)
}# End write_genepopedit
