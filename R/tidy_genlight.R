#' @name tidy_genlight
#' @title Tidy a genlight object to a tidy dataframe
#' @description Tidy genlight object to a tidy dataframe.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data A genlight object in the global environment.
#' @inheritParams tidy_genomic_data

#' @export
#' @rdname tidy_genlight

#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join
#' @importFrom stringi stri_replace_all_fixed stri_replace_all_regex stri_join
# @importFrom adegenet indNames pop chromosome locNames position
#' @importFrom tibble rownames_to_column data_frame
#' @importFrom tidyr gather unite

#' @references Jombart T (2008) adegenet: a R package for the multivariate
#' analysis of genetic markers. Bioinformatics, 24, 1403-1405.
#' @references Jombart T, Ahmed I (2011) adegenet 1.3-1:
#' new tools for the analysis of genome-wide SNP data.
#' Bioinformatics, 27, 3070-3071.

#' @details
#' A string of the same dimension is generated when genlight:
#' \enumerate{
#' \item \code{is.null(genlight@pop)}: pop1 will be integrated
#' in the tidy dataset.
#' \item \code{is.null(data@chromosome)}: CHROM1 will be integrated
#' in the tidy dataset.
#' \item \code{is.null(data@loc.names)}: LOCUS1 to dim(genlight)[2]
#' will be integrated in the tidy dataset.
#' \item \code{is.null(data@position)}: an integer string of
#' length = dim(genlight)[2] will be integrated in the tidy dataset.
#' }
#'
#' \strong{Note: that if all CHROM, LOCUS and POS is missing the function will be terminated}


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


tidy_genlight <- function(data, verbose = FALSE) {

  if (!"adegenet" %in% utils::installed.packages()[,"Package"]) {
    stop('Please install SeqArray for this option:\n
         devtools::install_github("zhengxwen/SeqArray")
         or the bioconductor version:
         source("https://bioconductor.org/biocLite.R")
         biocLite("SeqArray")')
  }

  ## test
  # data
  # verbose = TRUE

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("data argument required")
  if (class(data)[1] != "genlight") stop("Input is not a genlight object")

  # Import data ---------------------------------------------------------------

  if (verbose) message("genlight info:")
  # strata ?
  if (!is.null(data@pop)) {
    if (verbose) message("    strata: yes")
    strata.df <- tibble::data_frame(
      INDIVIDUALS = data@ind.names,
      POP_ID = data@pop)
  } else {
    if (verbose) message("    strata: no")
    strata.df <- tibble::data_frame(
      INDIVIDUALS = data@ind.names,
      POP_ID = rep("pop1", dim(data)[1]))
  }

  n.markers <- dim(data)[2]

  # Chromosome ?
  if (is.null(data@chromosome)) {
    if (verbose) message("    Chromosome/contig/scaffold: no")
    data@chromosome <- factor(rep("CHROM1", n.markers))
    chrom.info <- FALSE
  } else {
    if (verbose) message("    Chromosome/contig/scaffold: yes")
    chrom.info <- TRUE
  }

  # Locus ?
  if (is.null(data@loc.names)) {
    if (verbose) message("    Locus: no")
    locus.info <- FALSE
    data@loc.names <- stringi::stri_join("LOCUS", seq(from = 1, to = n.markers, by = 1))
  } else {
    if (verbose) message("    Locus: yes")
    locus.info <- TRUE
  }

  # POS ?
  if (is.null(data@position)) {
    if (verbose) message("    POS: no")
    pos.info <- FALSE
    data@position <- rlang::as_integer(seq(from = 1, to = n.markers, by = 1))
  } else {
    if (verbose) message("    POS: yes")
    pos.info <- TRUE
  }


  if (!chrom.info && !locus.info && !pos.info) {
    stop("Tidying the genlight requires at least one of these 3 markers metadata:
       CHROM (genlight@chromosome), LOCUS (genlight@loc.names) or POS (genlight@position)")
  }


  markers <- tibble::data_frame(
    CHROM = data@chromosome,#adegenet::chromosome(data),
    LOCUS = data@loc.names,#adegenet::locNames(data),
    POS = data@position#adegenet::position(data)
  ) %>%
    dplyr::mutate_all(.tbl = ., .funs = as.character) %>%
    dplyr::mutate_all(.tbl = ., .funs = radiator::clean_markers_names) %>%
    tidyr::unite(data = ., col = MARKERS, CHROM, LOCUS, POS, sep = "__", remove = FALSE) %>%
    dplyr::select(MARKERS, CHROM, LOCUS, POS)

  want <- c("MARKERS", "CHROM", "LOCUS", "POS", "POP_ID", "INDIVIDUALS",
            "GT_VCF", "GT_BIN", "GT")

  if (verbose) message("Generating tidy data...")
  data <- suppressWarnings(
    tibble::as_data_frame(data) %>%
      magrittr::set_colnames(x = ., value = markers$MARKERS) %>%
      tibble::add_column(.data = ., INDIVIDUALS = rownames(.), .before = 1) %>%
      data.table::as.data.table(.) %>%
      data.table::melt.data.table(
        data = .,
        id.vars = "INDIVIDUALS",
        variable.name = "MARKERS",
        value.name = "GT_BIN"
      ) %>%
      tibble::as_data_frame(.) %>%
      dplyr::full_join(markers, by = "MARKERS") %>%
      dplyr::full_join(strata.df, by =  "INDIVIDUALS") %>%
      dplyr::mutate(
        GT_VCF = GT_BIN,
        GT_VCF = stringi::stri_replace_all_regex(
          str = GT_BIN,
          pattern = stringi::stri_join("^", c("0", "2", "1"), "$", sep = ""),
          replacement = c("0/0", "1/1", "0/1"),
          vectorize_all = FALSE
        ),
        GT_VCF = replace(GT_VCF, which(is.na(GT_VCF)), "./."),
        GT = stringi::stri_replace_all_fixed(
          str = GT_VCF,
          pattern = c("0/0", "1/1", "0/1", "./."),
          replacement = c("001001", "002002", "001002", "000000"),
          vectorize_all = FALSE
        ),
        INDIVIDUALS = radiator::clean_ind_names(INDIVIDUALS),
        POP_ID = radiator::clean_pop_names(POP_ID)
      ) %>%
      dplyr::arrange(MARKERS, POP_ID, INDIVIDUALS) %>%
      dplyr::select(dplyr::one_of(want)))

return(data)
} # End tidy_genlight
