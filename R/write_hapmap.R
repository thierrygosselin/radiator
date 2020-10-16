# write_hapmap---------------------------------------------------------------------
# write a hapmap file from a tidy data frame

#' @name write_hapmap
#' @title Write a HapMap file from a tidy data frame
#' @description Write a \href{https://www.ncbi.nlm.nih.gov/variation/news/NCBI_retiring_HapMap/}{HapMap file} from a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#'
#' The data requires nucleotide \code{A, C, G, T} information for the genotypes.
#'
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.

#' @param filename (optional) The file name prefix for the hapmap file
#' written to the working directory. With default: \code{filename = NULL},
#' the date and time is appended to \code{radiator_hapmap}.

#' @export
#' @rdname write_hapmap

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

write_hapmap <- function(
  data,
  filename = NULL
) {
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    data <- radiator::tidy_wide(data = data, import.metadata = TRUE)
  }

  # REF/ALT Alleles and VCF genotype format ------------------------------------
  if (!tibble::has_name(data, "GT_VCF_NUC")) {
    ref.change <- radiator::calibrate_alleles(data = data)$input
    data <- dplyr::left_join(data, ref.change, by = c("MARKERS", "INDIVIDUALS"))
  }

  # remove duplicate REF/ALT column
  if (tibble::has_name(data, "REF.x")) {
    data <- dplyr::select(.data = data, -c(REF.x, ALT.x)) %>%
      dplyr::rename(REF = REF.y, ALT = ALT.y)
  }

  # Include CHROM, LOCUS, POS --------------------------------------------------
  if (!tibble::has_name(data, "CHROM")) {
    data <- dplyr::mutate(
      .data = data,
      CHROM = rep("1", n()),
      LOCUS = MARKERS,
      POS = MARKERS
    )
  }

  # Order/sort by pop and ind --------------------------------------------------
  if (tibble::has_name(data, "POP_ID")) {
    data <- dplyr::arrange(data, POP_ID, INDIVIDUALS)
  } else {
    data <- dplyr::arrange(data, INDIVIDUALS)
  }

  # Remove the POP_ID column ---------------------------------------------------
  if (tibble::has_name(data, "POP_ID")) {
    data <- dplyr::select(.data = data, -POP_ID)
  }

  # hapmap body  ------------------------------------------------------------------
  data %<>%
    dplyr::select(MARKERS, CHROM, LOCUS, POS, REF, ALT, INDIVIDUALS, GT_VCF_NUC) %>%
    dplyr::mutate(
      alleles = stringi::stri_join(REF, ALT, sep = "/"),
      REF = NULL,
      ALT = NULL,
      LOCUS = NULL,
      GT_VCF_NUC = stringi::stri_replace_all_regex(
        str = GT_VCF_NUC,
        pattern = "\\/",
        replacement = "",
        vectorize_all = FALSE
      ),
      GT_VCF_NUC = stringi::stri_replace_all_regex(
        str = GT_VCF_NUC,
        pattern = "\\.\\.",
        replacement = NA_character_,
        vectorize_all = FALSE
      )
    ) %>%
    data.table::as.data.table(.) %>%
    data.table::dcast.data.table(
      data = .,
      formula = MARKERS + alleles + CHROM + POS ~ INDIVIDUALS,
      value.var = "GT_VCF_NUC"
    ) %>%
    tibble::as_tibble(.) %>%
    dplyr::mutate(
      strand = "+",
      'assembly#' = NA_character_,
      center = "radiator",
      protLSID = NA_character_,
      assayLSID = NA_character_,
      panelLSID = NA_character_,
      QCcode = NA_character_
    ) %>%
    dplyr::select('rs#' = MARKERS, alleles, chrom = CHROM, pos = POS, strand,
                  'assembly#', center, protLSID, assayLSID, panelLSID, QCcode,
                  everything(.))

  # Filename ------------------------------------------------------------------
  if (is.null(filename)) {
    # Get date and time to have unique filenaming
    filename <- stringi::stri_join("radiator_hapmap_", file.date, ".tsv")
  } else {
    filename <- stringi::stri_join(filename, ".hapmap.tsv")
  }

  # Write the prunned vcf to the file ------------------------------------------
  suppressWarnings(readr::write_tsv(x = data, path = filename))
}# end write_vcf
