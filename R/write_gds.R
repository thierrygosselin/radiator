# write a GDS object from a tidy data frame

#' @name write_gds
#' @title Write a GDS object from a tidy data frame
#' @description Write a Genomic Data Structure (GDS) file
#' \href{http://zhengxwen.github.io/gdsfmt}{gdsfmt}
#' object from a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.

#' @param filename (optional) The file name of the Genomic Data Structure (GDS) file.
#' radiator will append \code{.gds} to the filename.
#' If filename chosen is already present in the
#' working directory, the default \code{radiator_datetime.gds} is chosen.
#' Default: \code{filename = NULL}.

#' @inheritParams tidy_genomic_data

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

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_gds <- function(data, filename = NULL, verbose = TRUE) {
  timing <- proc.time()

  ## testing
  # filename = NULL
  # verbose = TRUE

  # Check that snprelate is installed
  if (!requireNamespace("gdsfmt", quietly = TRUE)) {
    stop('Please install gdsfmt for this option:\n
       install.packages("BiocManager")
       BiocManager::install ("gdsfmt")
       ')
  }

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file missing")

  # Filename -------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (is.null(filename)) {
    filename <- stringi::stri_join("radiator_", file.date, ".gds")
  } else {
    filename.problem <- file.exists(stringi::stri_join(filename, ".gds"))
    if (filename.problem) {
      filename <- stringi::stri_join(filename, "_", file.date, ".gds")
    } else {
      filename <- stringi::stri_join(filename, ".gds")
    }
  }

  radiator.temp.file <- stringi::stri_join("radiator_temp_", file.date)

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    data <- radiator::tidy_wide(data = data, import.metadata = TRUE)
  }

  # GDS prep -------------------------------------------------------------------
  if (verbose) message("Generating GDS format...")

  # Generate GDS format --------------------------------------------------------
  n.markers <- length(unique(data$MARKERS))
  SNPRelate::snpgdsCreateGeno(
    gds.fn = radiator.temp.file,
    genmat = suppressWarnings(
      dplyr::select(.data = data, MARKERS, INDIVIDUALS, GT_BIN) %>%
        dplyr::group_by(MARKERS) %>%
        tidyr::spread(data = ., key = INDIVIDUALS, value = GT_BIN) %>%
        dplyr::arrange(MARKERS) %>%
        tibble::column_to_rownames(df = ., var = "MARKERS") %>%
        data.matrix(.)),
    sample.id = unique(data$INDIVIDUALS),
    snp.id = dplyr::distinct(.data = data, MARKERS) %>%
      dplyr::arrange(MARKERS) %>%
      purrr::flatten_chr(.),
    snp.rs.id = NULL,
    snp.chromosome = if (tibble::has_name(data, "CHROM")) {
      dplyr::select(data, MARKERS, CHROM) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
        dplyr::arrange(MARKERS) %>%
        dplyr::select(-MARKERS) %>%
        purrr::flatten_chr(.)
    } else {
      rep(1L, n.markers)
    },
    snp.position = NULL,
    snp.allele = if (tibble::has_name(data, "REF")) {
      dplyr::select(data, MARKERS, REF, ALT) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
        dplyr::arrange(MARKERS) %>%
        dplyr::mutate(REF_ALT = stringi::stri_join(REF, ALT, sep = "/")) %>%
        dplyr::select(REF_ALT) %>%
        purrr::flatten_chr(.)
    } else {
      rep("NA", n.markers)
    },
    snpfirstdim = TRUE,
    compress.annotation = "ZIP_RA",
    compress.geno = ""
  )

  # check <- SNPRelate::snpgdsOpen(radiator.temp.file)
  # check
  # SNPRelate::snpgdsClose(check)

  data.gds <- SeqArray::seqSNP2GDS(
    gds.fn = radiator.temp.file,
    out.fn = filename,
    storage.option = "ZIP_RA",
    major.ref = TRUE, verbose = FALSE) %>%
    SeqArray::seqOpen(gds.fn = ., readonly = FALSE)

  file.remove(radiator.temp.file)
  n.ind <- length(SeqArray::seqGetData(data.gds, "sample.id"))
  n.markers <- length(SeqArray::seqGetData(data.gds, "variant.id"))

  # generate a radiator folder inside the gds file -----------------------------
  radiator.gds <- gdsfmt::addfolder.gdsn(
    node = data.gds,
    name = "radiator",
    replace = TRUE)
  # Add STRATA to GDS ----------------------------------------------------------
  gdsfmt::add.gdsn(
    node = radiator.gds,
    name = "STRATA",
    val = dplyr::distinct(data, POP_ID, INDIVIDUALS) %>%
      dplyr::select(POP_ID),
    replace = TRUE,
    compress = "ZIP_RA",
    closezip = TRUE)

  # ADD MARKERS META to GDS ----------------------------------------------------
  gdsfmt::add.gdsn(
    node = radiator.gds,
    name = "markers.meta",
    val = tibble::add_column(
      .data = suppressWarnings(
        dplyr::select(
          data,
          dplyr::one_of(c("MARKERS", "CHROM", "LOCUS", "POS", "COL", "REF", "ALT")))) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
        dplyr::mutate(MARKERS = as.character(MARKERS)),
      VARIANT_ID = SeqArray::seqGetData(data.gds, "variant.id"),.before = 1),
    replace = TRUE,
    compress = "ZIP_RA",
    closezip = TRUE)

  suppressWarnings(data %<>% dplyr::select(
    -dplyr::one_of(c("MARKERS", "CHROM", "LOCUS", "POS", "COL", "REF", "ALT"))))

  # Add genotypes metadata  ----------------------------------------------------
  if (TRUE %in% tibble::has_name(data, c("ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH", "READ_DEPTH", "GL_HOM_REF",
                                         "GL_HET", "GL_HOM_ALT", "GQ"))) {
    gdsfmt::add.gdsn(
      node = radiator.gds,
      name = "genotypes.meta",
      val = suppressWarnings(
        dplyr::select(
          data,
          dplyr::one_of(
            c("ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH", "READ_DEPTH", "GL_HOM_REF",
              "GL_HET", "GL_HOM_ALT", "GQ")))),
      replace = TRUE,
      compress = "ZIP_RA",
      closezip = TRUE)
  }
  data <- NULL

  # Add to GDS if we're dealing with de novo data or not -----------------------
  vcf.locus <- SeqArray::seqGetData(data.gds, "annotation/id")
  sample.size <- min(length(unique(vcf.locus)), 100)
  gdsfmt::add.gdsn(
    node = radiator.gds,
    name = "reference.genome",
    val = sample(x = unique(vcf.locus), size = sample.size, replace = FALSE) %>%
      stringi::stri_detect_regex(str = ., pattern = "[^[:alnum:]]+") %>%
      unique,
    replace = FALSE,
    compress = "ZIP_RA",
    closezip = TRUE)
  sample.size <- vcf.locus <- NULL

  # bi- or multi-alllelic VCF ------------------------------------------------
  if (max(unique(SeqArray::seqNumAllele(gdsfile = data.gds))) - 1 > 1) {
    biallelic <- FALSE
    if (verbose) message("data is multi-allelic")
  } else {
    biallelic <- TRUE
    if (verbose) message("data is biallelic")
  }

  # add biallelic info to GDS
  gdsfmt::add.gdsn(
    node = radiator.gds,
    name = "biallelic",
    val = biallelic,
    replace = FALSE,
    compress = "ZIP_RA",
    closezip = TRUE)
  timing <- proc.time() - timing
  if (verbose) message("\nConversion time: ", round(timing[[3]]), " sec\n")
  return(data.gds)
} # End write_gds
