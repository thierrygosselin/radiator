# write a SNPRelate object from a tidy data frame

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

#' @export
#' @rdname write_snprelate

#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join
#' @importFrom stringi stri_replace_all_fixed
#' @importFrom tibble has_name
#' @importFrom tidyr spread
# @importFrom SNPRelate snpgdsOpen snpgdsClose snpgdsCreateGeno

#' @references Zheng X, Levine D, Shen J, Gogarten SM, Laurie C, Weir BS.
#' A high-performance computing toolset for relatedness and principal component
#' analysis of SNP data. Bioinformatics. 2012;28: 3326-3328.
#' doi:10.1093/bioinformatics/bts606

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_snprelate <- function(data, biallelic = TRUE, filename = NULL) {

  # Check that snprelate is installed
  if (!requireNamespace("SNPRelate", quietly = TRUE)) {
    stop('To install SNPRelate:\n
         source("https://bioconductor.org/biocLite.R")
         biocLite("SNPRelate"')
  }

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file missing")


  # Filename -------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (is.null(filename)) {
    filename <- stringi::stri_join("radiator_snprelate_", file.date, ".gds")
  } else {
    filename.problem <- file.exists(filename)
    if (filename.problem) {
      filename <- stringi::stri_join(filename, "_snprelate_", file.date, ".gds")
    } else {
      filename <- stringi::stri_join(filename, ".gds")
    }
  }


  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    data <- radiator::tidy_wide(data = data, import.metadata = TRUE)
  }

  # check genotype column naming
  colnames(data) <- stringi::stri_replace_all_fixed(
    str = colnames(data),
    pattern = "GENOTYPE",
    replacement = "GT",
    vectorize_all = FALSE
  )

  # necessary steps to make sure we work with unique markers and not duplicated LOCUS
  if (tibble::has_name(data, "LOCUS") && !tibble::has_name(data, "MARKERS")) {
    data <- dplyr::rename(.data = data, MARKERS = LOCUS)
  }

  # SNPRelate prep -------------------------------------------------------------
  if (is.null(biallelic)) biallelic <- radiator::detect_biallelic_markers(data = data)
  if (!biallelic) stop("SNPRelate requires biallelic genotypes")
  # verbose <- FALSE

  message("Generating GDS format...")
  # keep markers in common
  # gds.genotypes <- suppressMessages(radiator::keep_common_markers(data = input)$input)

  strata.df <- dplyr::distinct(data, POP_ID, INDIVIDUALS) %>%
    dplyr::mutate(POP_ID = factor(POP_ID))

  snp.id <- dplyr::distinct(.data = data, MARKERS) %>%
    dplyr::arrange(MARKERS) %>%
    purrr::flatten_chr(.)

  if (tibble::has_name(data, "CHROM")) {
    snp.chromosome <- dplyr::select(data, MARKERS, CHROM) %>%
      dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
      dplyr::arrange(MARKERS) %>%
      dplyr::select(-MARKERS) %>%
      purrr::flatten_chr(.)
  } else {
    snp.chromosome <- NULL
  }


  if (tibble::has_name(data, "REF")) {
    snp.allele <- dplyr::select(data, MARKERS, REF, ALT) %>%
      dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
      dplyr::arrange(MARKERS) %>%
      dplyr::mutate(REF_ALT = stringi::stri_join(REF, ALT, sep = "/")) %>%
      dplyr::select(REF_ALT) %>%
      purrr::flatten_chr(.)
  } else {
    snp.allele <- NULL
  }

  data <- suppressWarnings(
    dplyr::select(.data = data, MARKERS, INDIVIDUALS, GT_BIN) %>%
      dplyr::group_by(MARKERS) %>%
      tidyr::spread(data = ., key = INDIVIDUALS, value = GT_BIN) %>%
      dplyr::arrange(MARKERS) %>%
      tibble::column_to_rownames(df = ., var = "MARKERS") %>%
      data.matrix(.)
  )

  # Generate GDS format --------------------------------------------------------

  SNPRelate::snpgdsCreateGeno(
    gds.fn = filename,
    genmat = data,
    sample.id = strata.df$INDIVIDUALS,
    snp.id = snp.id,
    snp.rs.id = NULL,
    snp.chromosome = snp.chromosome,
    snp.position = NULL,
    snp.allele = snp.allele,
    snpfirstdim = TRUE,
    compress.annotation = "ZIP_RA.max",
    compress.geno = "",
    other.vars = NULL
  )

  gds.file.connection <- SNPRelate::snpgdsOpen(filename)

  message("\nNote: \n
GDS filename: ", filename)

  message("\nTo close the connection use SNPRelate::snpgdsClose(OBJECT_NAME)")
  return(gds.file.connection)
} # End write_snprelate
