# write a vcf file from a tidy data frame

#' @name write_vcf
#' @title Write a vcf file from a tidy data frame
#' @description Write a vcf file (file format version 4.3, see details below)
#' from a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.


#' @param pop.info (optional, logical) Should the population information be
#' included in the FORMAT field (along the GT info for each samples ?). To make
#' the VCF population-ready use \code{pop.info = TRUE}. The population information
#' must be included in the \code{POP_ID} column of the tidy dataset.
#' Default: \code{pop.info = FALSE}. Experimental.

#' @param filename (optional) The file name prefix for the vcf file
#' written to the working directory. With default: \code{filename = NULL},
#' the date and time is appended to \code{radiator_vcf_file_}.

#' @details \strong{VCF file format version:}
#'
#' If you need a different file format version than the current one, just change
#' the version inside the newly created VCF, that should do the trick.
#' \href{https://vcftools.github.io/specs.html}{For more
#' information on Variant Call Format specifications}.


#' @export
#' @rdname write_vcf
#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join everything
#' @importFrom stringi stri_join stri_replace_all_fixed stri_extract_all_fixed stri_replace_all_regex stri_sub stri_pad_left stri_count_fixed stri_replace_na
#' @importFrom readr write_tsv write_delim
#' @importFrom utils write.table packageVersion
#' @importFrom tibble has_name data_frame
#' @importFrom tidyr spread

#' @references Danecek P, Auton A, Abecasis G et al. (2011)
#' The variant call format and VCFtools.
#' Bioinformatics, 27, 2156-2158.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

write_vcf <- function(data, pop.info = FALSE, filename = NULL) {

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    data <- radiator::tidy_wide(data = data, import.metadata = TRUE)
  }

  # necessary steps to make sure we work with unique markers and not duplicated LOCUS
  if (tibble::has_name(data, "LOCUS") && !tibble::has_name(data, "MARKERS")) {
    data <- dplyr::rename(.data = data, MARKERS = LOCUS)
  }

  # REF/ALT Alleles and VCF genotype format ------------------------------------
  if (!tibble::has_name(data, "GT_VCF")) {
    ref.change <- radiator::change_alleles(data = data)$input
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
  data <- dplyr::arrange(data, POP_ID, INDIVIDUALS)

  id.string <- unique(data$INDIVIDUALS)# keep to sort vcf columns
  # Remove the POP_ID column ---------------------------------------------------
  if (tibble::has_name(data, "POP_ID") && (!pop.info)) {
    data <- dplyr::select(.data = data, -POP_ID)
  }

  # Info field -----------------------------------------------------------------
  info.field <- suppressWarnings(
    dplyr::select(.data = data, MARKERS, GT_VCF) %>%
      dplyr::filter(GT_VCF != "./.") %>%
      dplyr::count(x = ., MARKERS) %>%
      dplyr::mutate(INFO = stringi::stri_join("NS=", n, sep = "")) %>%
      dplyr::select(-n)
  )

  # VCF body  ------------------------------------------------------------------
  GT_VCF_POP_ID <- NULL
  if (pop.info) {
    output <- suppressWarnings(
      dplyr::left_join(data, info.field, by = "MARKERS") %>%
        dplyr::select(MARKERS, CHROM, LOCUS, POS, REF, ALT, INFO, INDIVIDUALS, GT_VCF, POP_ID) %>%
        dplyr::mutate(GT_VCF_POP_ID = stringi::stri_join(GT_VCF, POP_ID, sep = ":")) %>%
        dplyr::select(-c(GT_VCF, POP_ID)) %>%
        dplyr::group_by(MARKERS, CHROM, LOCUS, POS, INFO, REF, ALT) %>%
        tidyr::spread(data = ., key = INDIVIDUALS, value = GT_VCF_POP_ID) %>%
        dplyr::ungroup(.) %>%
        dplyr::mutate(
          QUAL = rep(".", n()),
          FILTER = rep("PASS", n()),
          FORMAT = rep("GT:POP", n())
        )
    )

  } else {
    output <- suppressWarnings(
      dplyr::left_join(data, info.field, by = "MARKERS") %>%
        dplyr::select(MARKERS, CHROM, LOCUS, POS, REF, ALT, INDIVIDUALS, GT_VCF, INFO) %>%
        dplyr::group_by(MARKERS, CHROM, LOCUS, POS, INFO, REF, ALT) %>%
        tidyr::spread(data = ., key = INDIVIDUALS, value = GT_VCF) %>%
        dplyr::ungroup(.) %>%
        dplyr::mutate(
          QUAL = rep(".", n()),
          FILTER = rep("PASS", n()),
          FORMAT = rep("GT", n())
        )
    )
  }

  # Transform the REF/ALT format back to A/C/G/T if 001, 002, etc is found
  ref.change <- TRUE %in% unique(c("001", "002", "003", "004") %in% unique(output$REF))

  if (ref.change) {
    output <- output %>%
      dplyr::mutate(
        REF = stringi::stri_replace_all_fixed(
          str = REF,
          pattern = c("001", "002", "003", "004"),
          replacement = c("A", "C", "G", "T"),
          vectorize_all = FALSE),
        ALT = stringi::stri_replace_all_fixed(
          str = ALT,
          pattern = c("001", "002", "003", "004"),
          replacement = c("A", "C", "G", "T"),
          vectorize_all = FALSE)
      )
  }

  # Keep the required columns
  output <- dplyr::ungroup(output) %>%
    dplyr::arrange(CHROM, LOCUS, POS) %>%
    dplyr::select(-MARKERS) %>%
    dplyr::select('#CHROM' = CHROM, POS, ID = LOCUS, REF, ALT, QUAL, FILTER, INFO, FORMAT, id.string)
  # dplyr::select('#CHROM' = CHROM, POS, ID = LOCUS, REF, ALT, QUAL, FILTER, INFO, FORMAT, dplyr::everything())


  # Filename ------------------------------------------------------------------
  if (is.null(filename)) {
    # Get date and time to have unique filenaming
    file.date <- format(Sys.time(), "%Y%m%d@%H%M")
    filename <- stringi::stri_join("radiator_vcf_file_", file.date, ".vcf")
  } else {
    filename <- stringi::stri_join(filename, ".vcf")
  }

  # File format ----------------------------------------------------------------
  readr::write_delim(
    x = tibble::data_frame("##fileformat=VCFv4.3"),
    path = filename, delim = " ", append = FALSE, col_names = FALSE)

  # File date ------------------------------------------------------------------
  file.date <- stringi::stri_replace_all_fixed(Sys.Date(), pattern = "-", replacement = "")
  file.date <- stringi::stri_join("##fileDate=", file.date, sep = "")
  readr::write_delim(x = tibble::data_frame(file.date), path = filename, delim = " ", append = TRUE, col_names = FALSE)

  # Source ---------------------------------------------------------------------
  readr::write_delim(x = tibble::data_frame(stringi::stri_join("##source=radiator_v.", as.character(utils::packageVersion("radiator")))), path = filename, delim = " ", append = TRUE, col_names = FALSE)

  # Info field 1 ---------------------------------------------------------------
  info1 <- as.data.frame('##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">')
  utils::write.table(x = info1, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)


  # Format field 1 -------------------------------------------------------------
  format1 <- '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
  format1 <- as.data.frame(format1)
  utils::write.table(x = format1, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)

  # Format field 2 ---------------------------------------------------------------
  if (pop.info) {
    format2 <- as.data.frame('##FORMAT=<ID=POP_ID,Number=1,Type=Character,Description="Population identification of Sample">')
    utils::write.table(x = format2, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  }

  # Write the prunned vcf to the file ------------------------------------------
  suppressWarnings(readr::write_tsv(x = output, path = filename, append = TRUE, col_names = TRUE))
}# end write_vcf
