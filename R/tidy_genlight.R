#' @name tidy_genlight
#' @title Tidy a genlight object to a tidy dataframe
#' @description Tidy genlight object to a tidy dataframe.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data A genlight object in the global environment.
#' @inheritParams tidy_genomic_data

#' @param tidy (logical) Generate a tidy dataset.
#' Default: \code{tidy = TRUE}.

#' @param gds (optional, logical) To write a radiator gds object.
#' Default: \code{gds = TRUE}.

#' @param write (optional, logical) To write in the working directory the tidy
#' data. The file is written with \code{radiator_genlight_DATE@TIME.rad}.
#' Default: \code{write = FALSE}.

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
#' \item \code{is.null(genlight@pop)}: pop will be integrated
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


tidy_genlight <- function(
  data,
  tidy = TRUE,
  gds = TRUE,
  write = FALSE,
  verbose = FALSE
) {

  if (!"adegenet" %in% utils::installed.packages()[,"Package"]) {
    rlang::abort('Please install adegenet for this option:\n
         install.packages("adegenet")
         ')
  }

  ## test
  # data
  # verbose = TRUE

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("data argument required")
  if (class(data)[1] != "genlight") rlang::abort("Input is not a genlight object")

  # Import data ---------------------------------------------------------------

  if (verbose) message("genlight info:")
  # strata ?
  strata <- tibble::tibble(INDIVIDUALS = data@ind.names)
  if (is.null(data@pop)) {
    if (verbose) message("    strata: no")
    if (verbose) message("    'pop' will be added")
    strata %<>% dplyr::mutate(STRATA = "pop")
    # strata <- tibble::tibble(
    #   INDIVIDUALS = data@ind.names,
    #   POP_ID = rep("pop1", dim(data)[1]))
  } else {
    if (verbose) message("    strata: yes")
    strata$STRATA = data@pop
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
    rlang::abort("Tidying the genlight requires at least one of these 3 markers metadata:
       CHROM (genlight@chromosome), LOCUS (genlight@loc.names) or POS (genlight@position)")
  }


  # markers
  markers <- tibble::tibble(
    CHROM = data@chromosome,#adegenet::chromosome(data),
    LOCUS = data@loc.names,#adegenet::locNames(data),
    POS = data@position#adegenet::position(data)
  ) %>%
    dplyr::mutate_all(.tbl = ., .funs = as.character) %>%
    dplyr::mutate_all(.tbl = ., .funs = radiator::clean_markers_names) %>%
    tidyr::unite(data = ., col = MARKERS, CHROM, LOCUS, POS, sep = "__", remove = FALSE) %>%
    dplyr::select(MARKERS, CHROM, LOCUS, POS)

  if (tidy) {
    if (write) {
      filename.temp <- generate_filename(extension = "rad")
      filename.short <- filename.temp$filename.short
      filename.genlight <- filename.temp$filename
    }

    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "POP_ID", "INDIVIDUALS",
              "GT_VCF", "GT_BIN", "GT")

    if (verbose) message("Generating tidy data...")
    data <- suppressWarnings(
      data.frame(data) %>%
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
        dplyr::full_join(strata, by =  "INDIVIDUALS") %>%
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

    if (write) {
      radiator::write_rad(data = data, path = filename.genlight)
      if (verbose) message("File written: ", filename.short)
    }

  }#End tidy genlight
  if (gds) {

    markers %<>% dplyr::mutate(VARIANT_ID = as.integer(factor(MARKERS)))

    gds.filename <- radiator_gds(
      genotypes.df = tibble::as_tibble(data.frame(data) %>% t) %>%
        tibble::add_column(.data = ., VARIANT_ID = markers$VARIANT_ID, .before = 1) %>%
        dplyr::arrange(VARIANT_ID),
      strata = dplyr::rename(strata, STRATA = POP_ID),
      biallelic = TRUE,
      markers.meta = markers,
      filename = NULL,
      verbose = verbose
    )
    if (verbose) message("Written: GDS filename: ", gds.filename)
  }# End gds genlight
  return(data)
} # End tidy_genlight
