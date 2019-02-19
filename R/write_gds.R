# write a GDS object from a tidy data frame

#' @name write_gds
#' @title Write a GDS object from a tidy data frame
#' @description Write a Genomic Data Structure (GDS) file format
#' \href{http://zhengxwen.github.io/gdsfmt}{gdsfmt}
#' and object of class \code{SeqVarGDSClass} from a tidy data frame.
#'
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.

#' @param source (optional, character) The name of the software that
#' generated the data. e.g. \code{source = "Stacks v.2.2"}.
#' Default: \code{source = NULL}.

#' @param filename (optional) The file name of the Genomic Data Structure (GDS) file.
#' radiator will append \code{.gds.rad} to the filename.
#' If filename chosen is already present in the
#' working directory, the default \code{radiator_datetime.gds.rad} is chosen.
#' Default: \code{filename = NULL}.

#' @param open (optional, logical) Open or not the radiator connection.
#' Default: \code{open = TRUE}.

#' @inheritParams tidy_genomic_data

#' @return An object in the global environment of class \code{SeqVarGDSClass} and
#' a file in the working directory.

#' @export
#' @rdname write_gds


#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join
#' @importFrom stringi stri_replace_all_fixed
#' @importFrom tibble has_name
#' @importFrom tidyr spread

#' @references Zheng X, Levine D, Shen J, Gogarten SM, Laurie C, Weir BS.
#' A high-performance computing toolset for relatedness and principal component
#' analysis of SNP data. Bioinformatics. 2012;28: 3326-3328.
#' doi:10.1093/bioinformatics/bts606
#' @references Zheng X, Gogarten S, Lawrence M, Stilp A, Conomos M, Weir BS,
#' Laurie C, Levine D (2017). SeqArray -- A storage-efficient high-performance
#' data format for WGS variant calls.
#' Bioinformatics.

#' @examples
#' \dontrun{
#' require(SeqVarTools)
#' data.gds <- radiator::write_gds(data = "shark.rad")
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_gds <- function(
  data,
  source = NULL,
  filename = NULL,
  open = TRUE,
  verbose = TRUE
  ) {
  timing <- proc.time()

  ## testing
  # filename = NULL
  # verbose = TRUE

  # Check that SeqVarTools is installed (it requires automatically: SeqArray and gdsfmt)
  if (!"SeqVarTools" %in% utils::installed.packages()[,"Package"]) {
    rlang::abort('Please install SeqVarTools for this option:\n
         install.packages("BiocManager")
         BiocManager::install("SeqVarTools")
         ')
  }

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("Input file missing")

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    data <- radiator::tidy_wide(data = data, import.metadata = TRUE)
  }

  biallelic <- radiator::detect_biallelic_markers(data)
  if (!biallelic) rlang::abort("Biallelic data required")

  # GDS prep -------------------------------------------------------------------
  if (verbose) message("Generating radiator GDS format...")
  n.markers <- length(unique(data$MARKERS))
  n.ind <- length(unique(data$INDIVIDUALS))

  strata <-
    suppressWarnings(
      data %>%
        dplyr::select(INDIVIDUALS, dplyr::one_of(c("STRATA", "POP_ID", "TARGET_ID")))
    )
  if (tibble::has_name(strata, "POP_ID")) strata  %<>% dplyr::rename(STRATA = POP_ID)

  want <- c("VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS", "COL", "REF",
  "ALT", "CALL_RATE", "REP_AVG", "AVG_COUNT_REF",
  "AVG_COUNT_SNP", "SEQUENCE")
  markers.meta <-
    suppressWarnings(
      data %>%
    dplyr::select(
      dplyr::one_of(want))
    )
  if (ncol(markers.meta) == 0) markers.meta <- NULL

  want <- c(
    "GT", "GT_VCF", "GT_VCF_NUC",
    "READ_DEPTH", "ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH",
    "GQ",
    "GL_HOM_REF", "GL_HET", "GL_HOM_ALT",
    "DP", "AD", "GL", "PL", "HQ", "GOF", "NR", "NV", "RO", "QR", "AO", "QA")
  genotypes.meta <-
    suppressWarnings(
      data %>%
        dplyr::select(
          dplyr::one_of(want))
    )
  want <- NULL
  if (ncol(genotypes.meta) == 0) genotypes.meta <- NULL

  data.gds <- radiator_gds(
    genotypes.df = data,
    strata = strata,
    biallelic = TRUE,
    markers.meta = markers.meta,
    genotypes.meta = genotypes.meta,
    filename = filename,
    source = source,
    verbose = verbose
  )
  markers.meta <- genotypes.meta <- data <- strata <- NULL

  # RETURN ---------------------------------------------------------------------
  if (verbose) message("Number of individuals: ", n.ind)
  if (verbose) message("Number of markers: ", n.markers)

  timing <- proc.time() - timing
  if (verbose) message("\nConversion time: ", round(timing[[3]]), " sec\n")

  if (open) data.gds <- read_rad(data.gds, verbose = FALSE)
  return(data.gds)
} # End write_gds
