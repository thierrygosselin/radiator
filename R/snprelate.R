# write_snprelate --------------------------------------------------------------
#' @name write_snprelate
#' @title Write a SNPRelate object from a tidy data frame
#' @description Write a \href{https://github.com/zhengxwen/SNPRelate}{SNPRelate}
#' object from a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.
#' \strong{The genotypes are biallelic.}

#' @param biallelic (logical, optional) If you already know that the data is
#' biallelic use this argument to speed up the function.
#' Default: \code{biallelic = TRUE}.

#' @param filename (optional) The file name of the Genomic Data Structure (GDS) file.
#' radiator will append \code{.gds} to the filename.
#' If filename chosen is already present in the
#' working directory, the default \code{radiator_snprelate_datetime.gds} is chosen.
#' Default: \code{filename = NULL}.

#' @inheritParams tidy_genomic_data

#' @export
#' @rdname write_snprelate

#' @references Zheng X, Levine D, Shen J, Gogarten SM, Laurie C, Weir BS.
#' A high-performance computing toolset for relatedness and principal component
#' analysis of SNP data. Bioinformatics. 2012;28: 3326-3328.
#' doi:10.1093/bioinformatics/bts606

#' @seealso \href{https://github.com/zhengxwen/SNPRelate}{SNPRelate}

#' @return An object in the global environment of class
#' \code{"SNPGDSFileClass", "gds.class"} and
#' a file in the working directory.

#' @examples
#' \dontrun{
#' require(SNPRelate)
#' data.gds <- radiator::write_snprelate(data = "shark.rad")
#' }


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_snprelate <- function(data, biallelic = TRUE, filename = NULL, verbose = TRUE) {

  #TEST
  # biallelic = TRUE
  # filename = NULL
  # verbose = TRUE

  # Check that snprelate is installed
  radiator_packages_dep(package = "SNPRelate", cran = FALSE, bioc = TRUE)

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("Input file missing")
  if (verbose) message("Generating SNPRelate object/file...")

  # Filename -------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")

  if (!is.null(filename)) filename <- stringi::stri_join(filename, "_snprelate")
  filename <- generate_filename(name.shortcut = filename, extension = "gds.rad")

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    data <- radiator::tidy_wide(data = data, import.metadata = TRUE)
  }

  biallelic <- radiator::detect_biallelic_markers(data)
  if (!biallelic) rlang::abort("Biallelic data required")

  # Generate GDS format --------------------------------------------------------
  if (!rlang::has_name(data, "POP_ID")) data %<>% dplyr::mutate(POP_ID = "pop")
  data %<>% dplyr::arrange(MARKERS, INDIVIDUALS, POP_ID) %>%
    dplyr::mutate(VARIANT_ID = as.integer(factor(MARKERS))) %>%
    dplyr::ungroup(.)

  # MARKERS META WORK ----------------------------------------------------------
  markers.meta <- separate_markers(
    data = data,
    sep = "__",
    markers.meta.all.only = TRUE,
    biallelic = biallelic,
    verbose = verbose) %>%
    dplyr::arrange(VARIANT_ID)
  n.markers <- nrow(markers.meta)

  notwanted <- c("MARKERS", "CHROM", "LOCUS", "POS", "COL", "REF", "ALT")
  data %<>% dplyr::select(-tidyselect::any_of(notwanted))

  # genotypes metadata
  want <- c("ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH", "READ_DEPTH", "GL_HOM_REF",
            "GL_HET", "GL_HOM_ALT", "GQ", "DP", "AD", "GL", "PL", "GQ", "HQ", "GOF", "NR", "NV")
  if (TRUE %in% rlang::has_name(data, want)) {
    genotypes.meta <- data %>% dplyr::select(tidyselect::any_of(want))
    data %<>% dplyr::select(-tidyselect::any_of(want))
  } else {
    genotypes.meta <- NULL
  }

  # strata
  strata <- dplyr::distinct(data, POP_ID, INDIVIDUALS) %>%
    dplyr::rename(STRATA = POP_ID) %>%
    dplyr::arrange(INDIVIDUALS)

  # prep data
  data %<>%
    dplyr::select(VARIANT_ID, INDIVIDUALS, GT_BIN) %>%
    radiator::rad_wide(
      x = .,
      formula = "VARIANT_ID ~ INDIVIDUALS",
      values_from = "GT_BIN"
    ) %>%
    dplyr::arrange(VARIANT_ID) %>%
    magrittr::set_rownames(x = ., value = .$VARIANT_ID) %>%
    dplyr::select(-VARIANT_ID) %>%
    data.matrix(.) %>%
    # 4 steps for SNPRelate genotype coding (change from ALT to REF dosage)
    magrittr::inset(is.na(.), 3L) %>%
    magrittr::inset(. == 0L, 9L) %>%
    magrittr::inset(. == 2L, 0L) %>%
    magrittr::inset(. == 9L, 2L)


  SNPRelate::snpgdsCreateGeno(
    gds.fn = filename$filename,
    genmat = data,
    sample.id = strata$INDIVIDUALS,
    snp.id = markers.meta$VARIANT_ID,
    snp.rs.id = markers.meta$LOCUS,
    snp.chromosome = markers.meta$CHROM,
    snp.position = markers.meta$POS,
    snp.allele = stringi::stri_join(markers.meta$REF, markers.meta$ALT, sep = "/"),
    snpfirstdim = TRUE,
    compress.annotation = "ZIP_RA",
    compress.geno = ""
  )
  data <- SNPRelate::snpgdsOpen(filename$filename, readonly = FALSE, allow.fork = TRUE)

  if (verbose) message("SNPRelate GDS: ", filename$filename.short)

  # timing <- proc.time() - timing
  # if (verbose) message("\nConversion time: ", round(timing[[3]]), " sec\n")
  return(data)
} # End write_snprelate
