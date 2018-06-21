# write a SeqArray object from a tidy data frame

#' @name write_seqarray
#' @title Write a SeqArray GDS file from a vcf file and generate a connection object.
#' @description Write a SeqArray \href{http://zhengxwen.github.io/SeqArray/}{SeqArray}
#' file (Zheng et al. 2017) from a vcf file and generate a connection object.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param vcf The path to the vcf file.

#' @param filename (optional) The file name of the Genomic Data Structure (GDS) file.
#' radiator will append \code{.gds} to the filename.
#' If filename chosen is already present in the
#' working directory, the default \code{radiator_datetime.gds} is chosen.
#' Default: \code{filename = NULL}.

#' @inheritParams tidy_genomic_data

#' @param write.info (logical) When \code{write.info = TRUE}, the markers metadata
#' (CHROM+LOCUS+POS and missingness) and individuals names and missingness is
#' written in the working directory.
#' Default: \code{write.info = FALSE}.

#' @details
#' A vcf file of 35 GB with ~4 millions SNPs take about ~7 min with 8 CPU.
#' A vcf file of 21 GB with ~2 millions SNPs take about ~5 min with 7 CPU.
#'
#' After the file is generated, it's a matter of sec to open a connection.

#' @export
#' @rdname write_seqarray

#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join
#' @importFrom stringi stri_replace_all_fixed
#' @importFrom tibble has_name
#' @importFrom tidyr spread

#' @references Zheng X, Gogarten S, Lawrence M, Stilp A, Conomos M, Weir BS,
#' Laurie C, Levine D (2017). SeqArray -- A storage-efficient high-performance
#' data format for WGS variant calls.
#' Bioinformatics.
#'
#' @references Danecek P, Auton A, Abecasis G et al. (2011)
#' The variant call format and VCFtools.
#' Bioinformatics, 27, 2156-2158.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_seqarray <- function(
  vcf,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
  write.info = FALSE
) {

  #Test
  # vcf = populations.snps.vcf
  # filename = NULL
  # parallel.core = parallel::detectCores() - 1
  # verbose = TRUE


  res <- list() #store the results
  timing.import <- proc.time()
  # Check that SeqArray is installed
  if (!"SeqArray" %in% utils::installed.packages()[,"Package"]) {
    stop('Please install SeqArray for this option:\n
         devtools::install_github("zhengxwen/SeqArray")')
  }

  # Checking for missing and/or default arguments ------------------------------
  if (missing(vcf)) stop("vcf file missing")
  message("\nReading VCF...")

  # Get file size
  big.vcf <- file.size(vcf)

  if (verbose) {
    if (big.vcf > 500000000) message("Large vcf file may take several minutes...")
    if (big.vcf > 5000000000) message("Actually, you have time for a coffee...")
  }

  # Filename -------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (is.null(filename)) {
    ind.file <- stringi::stri_join("vcf_individuals_info_", file.date, ".tsv")
    markers.file <- stringi::stri_join("vcf_markers_metadata_", file.date, ".tsv")
    filename <- stringi::stri_join("radiator_", file.date, ".gds")
  } else {
    ind.file <- stringi::stri_join(filename, "_vcf_individuals_info_", file.date, ".tsv")
    markers.file <- stringi::stri_join(filename, "_vcf_markers_metadata_", file.date, ".tsv")
    filename.problem <- file.exists(filename)
    if (filename.problem) {
      filename <- stringi::stri_join(filename, "_", file.date, ".gds")
    } else {
      filename <- stringi::stri_join(filename, ".gds")
    }
  }

  # Read vcf -------------------------------------------------------------------
  res$vcf.connection <- SeqArray::seqVCF2GDS(
    vcf.fn = vcf,
    out.fn = filename,
    parallel = parallel.core,
    storage.option = "ZIP_RA",
    verbose = FALSE) %>%
    SeqArray::seqOpen(gds.fn = .)

  if (verbose) {
    message("\nSeqArray GDS file generated: ", filename)
    message("\nTo close the connection use SeqArray::seqClose(OBJECT_NAME$vcf.connection)")
  }

  res$markers.meta <- tibble::tibble(
    CHROM = SeqArray::seqGetData(res$vcf.connection, "chromosome"),
    LOCUS = SeqArray::seqGetData(res$vcf.connection, "variant.id"),
    POS = SeqArray::seqGetData(res$vcf.connection, "position")
  ) %>%
    dplyr::mutate(MARKERS = stringi::stri_join(CHROM, LOCUS, POS, sep = "__")) %>%
    dplyr::bind_cols(
      tibble::tibble(
        MISSING_PROP = SeqArray::seqMissing(
          gdsfile = res$vcf.connection,
          per.variant = TRUE, .progress = TRUE, parallel = parallel.core))) %>%
    dplyr::select(MARKERS, CHROM, LOCUS, POS, MISSING_PROP)
  if (write.info) readr::write_tsv(x = res$markers.meta, path = markers.file)

  res$individuals <- tibble::tibble(
    INDIVIDUALS = SeqArray::seqGetData(res$vcf.connection, "sample.id"),
    MISSING_PROP = SeqArray::seqMissing(
      gdsfile = res$vcf.connection, per.variant = FALSE,
      .progress = TRUE,
      parallel = parallel.core))
  if (write.info) readr::write_tsv(x = res$individuals, path = ind.file)

  # Stats ----------------------------------------------------------------------
  res$n.chromosome <- dplyr::n_distinct(res$markers.meta$CHROM)
  res$n.locus <- dplyr::n_distinct(res$markers.meta$LOCUS)
  res$n.markers <- dplyr::n_distinct(res$markers.meta$MARKERS)
  res$n.individuals <- nrow(res$individuals)
  res$ind.missing <- round(mean(res$individuals$MISSING_PROP, na.rm = TRUE), 2)
  res$markers.missing <- round(mean(res$markers.meta$MISSING_PROP, na.rm = TRUE), 2)

  if (verbose) {
    message("\n\nNumber of chromosome/contig/scaffold: ", res$n.chromosome)
    message("Number of locus: ", res$n.locus)
    message("Number of markers: ", res$n.markers)
    message("Number of individuals: ", res$n.individuals)
    message("\n\nMissing data (averaged): ")
    message("    markers: ", res$markers.missing)
    message("    individuals: ", res$ind.missing)
    timing.import <- proc.time() - timing.import
    message(stringi::stri_join("\nImport time: ", round(timing.import[[3]]), " sec"))
  }
  return(res)
} # End write_SeqArray
