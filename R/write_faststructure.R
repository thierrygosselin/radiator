# write a faststructure file from a tidy data frame

#' @name write_faststructure
#' @title Write a faststructure file from a tidy data frame
#' @description Write a faststructure file from a tidy data frame.
#' For biallelic dataset only.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @inheritParams radiator_common_arguments

#' @param plink.bed To export in plink BED format. This will write 3 files:
#' \code{.bed}, \code{.bim}, \code{.fam}. For this option to work, the argument
#' data must be a GDS file or object.
#' Default: \code{plink.bed = FALSE}.

#' @inheritParams read_strata

#' @param filename (optional) The file name prefix for the faststructure file
#' written to the working directory. With default: \code{filename = NULL},
#' the date and time is appended to \code{radiator_faststructure_}.

#' @param ... other parameters passed to the function.

#' @return A faststructure file is saved to the working directory.

#' @export
#' @rdname write_faststructure

#' @references Raj A, Stephens M, Pritchard JK (2014)
#' fastSTRUCTURE: Variational Inference of Population Structure in Large SNP
#' Datasets. Genetics, 197, 573-589.


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

write_faststructure <- function(
  data,
  plink.bed = FALSE,
  pop.levels = NULL,
  filename = NULL,
  ...
) {

  # Checking for missing and/or default arguments ******************************
  if (missing(data)) rlang::abort("Input file missing")

  if (plink.bed) {
    # Required packages
    if (!requireNamespace("SNPRelate", quietly = TRUE)) {
      rlang::abort('To install SNPRelate:\n
                   install.packages("BiocManager")
                   BiocManager::install("SNPRelate")
                   ')
    }


    # File type detection----------------------------------------------------------
    data.type <- radiator::detect_genomic_format(data)
    # Import data ---------------------------------------------------------------

    if (data.type %in% c("SeqVarGDSClass", "gds.file")) {
      if (data.type == "gds.file") {
        data <- radiator::read_rad(data, verbose = FALSE)
      }
    } else {
      rlang::abort("plink.bed option requires GDS file or object")
    }

    SeqArray::seqGDS2SNP(gdsfile = data, out.gdsfn = "radiator.temp", verbose = TRUE)
    temp <- SNPRelate::snpgdsOpen(
      filename = "radiator.temp",
      readonly = TRUE,
      allow.duplicate = TRUE,
      allow.fork = TRUE
      )
    SNPRelate::snpgdsGDS2BED(gdsobj = temp, bed.fn = filename, verbose = TRUE)
    # Converting from GDS to PLINK binary PED:
    #   Working space: 208 samples, 8719 SNPs
    # Error in if ((opt$autosome.start == 1) & (opt$autosome.end == 22)) { :
    #     missing value where TRUE/FALSE needed
    file.remove(temp$filename)
  } else {
    # Import data ---------------------------------------------------------------
    if (is.vector(data)) {
      data <- radiator::tidy_wide(data = data, import.metadata = FALSE)
    }

    want <- c("INDIVIDUALS", "POP_ID", "MARKERS", "GT", "GT_BIN")
    suppressWarnings(data %<>% dplyr::select(dplyr::one_of(want)))

    # pop.levels -----------------------------------------------------------------
    if (!is.null(pop.levels)) {
      data <- dplyr::mutate(
        .data = data,
        POP_ID = factor(POP_ID, levels = pop.levels, ordered = TRUE),
        POP_ID = droplevels(POP_ID)
      ) %>%
        dplyr::arrange(POP_ID, INDIVIDUALS, MARKERS)
    } else {
      data <- dplyr::mutate(.data = data, POP_ID = factor(POP_ID)) %>%
        dplyr::arrange(POP_ID, INDIVIDUALS, MARKERS)
    }

    # Create a marker vector  ------------------------------------------------
    # markers <- dplyr::distinct(.data = data, MARKERS) %>%
    #   dplyr::arrange(MARKERS) %>%
    #   purrr::flatten_chr(.)

    # faststructure format ----------------------------------------------------------------
    if (rlang::has_name(data, "GT_BIN")) {
      want <- c("INDIVIDUALS", "POP_ID", "MARKERS", "GT_BIN")
      suppressWarnings(
        data %<>%
          dplyr::select(dplyr::one_of(want)) %>%
          dplyr::mutate(
            A1 = dplyr::case_when(
              GT_BIN == 0 ~ 1L,
              GT_BIN == 1 ~ 1L,
              GT_BIN == 2 ~ 2L,
              is.na(GT_BIN) ~ -9L
            ),
            A2 = dplyr::case_when(
              GT_BIN == 0 ~ 1L,
              GT_BIN == 1 ~ 2L,
              GT_BIN == 2 ~ 2L,
              is.na(GT_BIN) ~ -9L
            ),
            GT_BIN = NULL
          ) %>%
          tidyr::gather(data = ., key = ALLELES, value = GT, -c(POP_ID, INDIVIDUALS, MARKERS))
      )
    } else {
      want <- c("INDIVIDUALS", "POP_ID", "MARKERS", "GT")
      suppressWarnings(
        data %<>%
          dplyr::select(dplyr::one_of(want)) %>%
          tidyr::separate(col = GT, into = c("A1", "A2"), sep = 3, extra = "drop", remove = TRUE) %>%
          tidyr::gather(data = ., key = ALLELES, value = GT, -c(POP_ID, INDIVIDUALS, MARKERS)) %>%
          dplyr::mutate(
            GT = stringi::stri_replace_all_fixed(str = GT, pattern = "000", replacement = "-9", vectorize_all = FALSE),
            GT = as.integer(GT)
          )
      )
    }

    # common to both GT and GT_BIN
    data  %<>%
      dplyr::select(INDIVIDUALS, POP_ID, MARKERS, ALLELES, GT) %>%
      tidyr::spread(data = ., key = MARKERS, value = GT) %>%
      dplyr::mutate(POP_ID = as.integer(POP_ID)) %>%
      dplyr::select(-ALLELES)

    markers.col <- purrr::keep(
      .x = colnames(data),
      .p = !colnames(data) %in% c("POP_ID", "INDIVIDUALS")
    )

    data %<>%
      dplyr::mutate(C3 = 1L, C4 = 1L, C5 = 1L, C6 = 1L) %>%
      dplyr::select(INDIVIDUALS, POP_ID, C3, C4, C5, C6, dplyr::everything(markers.col)) %>%
      dplyr::arrange(POP_ID, INDIVIDUALS)

    # Write the file in faststructure format -----------------------------------------
    # Filename
    if (is.null(filename)) {
      file.date <- format(Sys.time(), "%Y%m%d@%H%M")
      filename <- stringi::stri_join("radiator_", file.date, ".faststructure.str")
    } else {
      filename <- stringi::stri_join(filename, ".faststructure.str")
    }

    # filename.connection <- file(filename, "w") # open the connection to the file
    # writeLines(text = stringi::stri_join(markers, sep = "\t", collapse = "\t"),
    #            con = filename.connection, sep = "\n")
    # close(filename.connection) # close the connection
    readr::write_tsv(x = data, path = filename, append = FALSE, col_names = FALSE)
  }
  return(filename)
} # end write_faststructure
