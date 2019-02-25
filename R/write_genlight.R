# write a genind file from a tidy data frame

#' @name write_genlight
#' @title Write a genlight object from a tidy data frame
#' @description Write a genlight object from a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @inheritParams radiator_common_arguments

#' @param biallelic (logical, optional) If you already know that the data is
#' biallelic use this argument to speed up the function.
#' Default: \code{biallelic = TRUE}.

#' @param write (logical, optional) To write in the working directory the genlight
#' object. The file is written with \code{radiator_genlight_DATE@TIME.RData} and
#' can be open with load or readRDS.
#' Default: \code{write = FALSE}.

#' @export
#' @rdname write_genlight

#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join
#' @importFrom stringi stri_replace_all_fixed
#' @importFrom methods new
# @importFrom adegenet indNames pop chromosome locNames position
#' @importFrom tibble has_name
#' @importFrom tidyr spread

#' @references Jombart T (2008) adegenet: a R package for the multivariate
#' analysis of genetic markers. Bioinformatics, 24, 1403-1405.
#' @references Jombart T, Ahmed I (2011) adegenet 1.3-1:
#' new tools for the analysis of genome-wide SNP data.
#' Bioinformatics, 27, 3070-3071.


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_genlight <- function(data, biallelic = TRUE, write = FALSE, verbose = FALSE) {

  # Checking for missing and/or default arguments ------------------------------
  if (!requireNamespace("adegenet", quietly = TRUE)) {
    rlang::abort("adegenet needed for this function to work
         Install with install.packages('adegenet')")
  }

  if (missing(data)) rlang::abort("Input file missing")

  # File type detection----------------------------------------------------------
  data.type <- radiator::detect_genomic_format(data)

  # Import data ---------------------------------------------------------------
  if (data.type %in% c("SeqVarGDSClass", "gds.file")) {
    if (data.type == "gds.file") {
      data <- radiator::read_rad(data, verbose = verbose)
    }
    data <- gds2tidy(gds = data, parallel.core = parallel::detectCores() - 1)
    data.type <- "tbl_df"
  } else {
    want <- c("MARKERS", "POP_ID", "INDIVIDUALS", "REF", "ALT", "GT", "GT_BIN")
    data <- suppressWarnings(radiator::tidy_wide(data = data, import.metadata = TRUE) %>%
                               dplyr::select(dplyr::one_of(want)) %>%
                               dplyr::arrange(POP_ID, INDIVIDUALS))
  }

  # Detect if biallelic data
  if (is.null(biallelic)) biallelic <- radiator::detect_biallelic_markers(data = input)
  if (!biallelic) rlang::abort("genlight object requires biallelic genotypes")

  # data = input
  input %<>% dplyr::arrange(MARKERS, CHROM, LOCUS, POS, POP_ID, INDIVIDUALS)
  marker.meta <- dplyr::distinct(.data = input, MARKERS, CHROM, LOCUS, POS)

  if (!tibble::has_name(input, "GT_BIN")) {
    input$GT_BIN <- stringi::stri_replace_all_fixed(
      str = input$GT_VCF,
      pattern = c("0/0", "1/1", "0/1", "1/0", "./."),
      replacement = c("0", "2", "1", "1", NA),
      vectorize_all = FALSE
    )
  }

  input <- dplyr::select(.data = input, MARKERS, POP_ID, INDIVIDUALS, GT_BIN) %>%
    dplyr::mutate(GT_BIN = as.integer(GT_BIN)) %>%
    dplyr::group_by(INDIVIDUALS, POP_ID) %>%
    tidyr::spread(data = ., key = MARKERS, value = GT_BIN) %>%
    dplyr::ungroup(.)

  # Generate genlight
  genlight.object <- methods::new(
    "genlight",
    input[,-(1:2)],
    parallel = FALSE
  )
  adegenet::indNames(genlight.object)   <- input$INDIVIDUALS
  adegenet::pop(genlight.object)        <- input$POP_ID
  adegenet::chromosome(genlight.object) <- marker.meta$CHROM
  adegenet::locNames(genlight.object)   <- marker.meta$LOCUS
  adegenet::position(genlight.object)   <- marker.meta$POS


  # Check
  # genlight.object@n.loc
  # genlight.object@ind.names
  # genlight.object@chromosome
  # genlight.object@position
  # genlight.object@loc.names
  # genlight.object@pop
  # genlight.object@strata
  # adegenet::nLoc(genlight.object)
  # adegenet::popNames(genlight.object)
  # adegenet::indNames(genlight.object)
  # adegenet::nPop(genlight.object)
  # adegenet::NA.posi(genlight.object)

  if (write) {
      filename.temp <- generate_filename(extension = "genlight")
      filename.short <- filename.temp$filename.short
      filename.genlight <- filename.temp$filename
      saveRDS(object = genlight.object, file = filename.genlight)
      if (verbose) message("File written: ", filename.short)
  }
  return(genlight.object)
} # End write_genlight
