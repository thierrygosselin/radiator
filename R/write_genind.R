# write a genind file from a tidy data frame

#' @name write_genind
#' @title Write a genind object from a tidy data frame

#' @description Write a genind object from a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.

#' @export
#' @rdname write_genind

#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join funs
#' @importFrom tidyr spread gather separate complete
#' @importFrom stringi stri_join stri_replace_all_fixed stri_sub stri_replace_na
#' @importFrom tibble has_name as_data_frame
#' @importFrom adegenet genind

#' @references Jombart T (2008) adegenet: a R package for the multivariate
#' analysis of genetic markers. Bioinformatics, 24, 1403-1405.
#' @references Jombart T, Ahmed I (2011) adegenet 1.3-1:
#' new tools for the analysis of genome-wide SNP data.
#' Bioinformatics, 27, 3070-3071.


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_genind <- function(data) {

  # Checking for missing and/or default arguments ******************************
  if (missing(data)) stop("Input file missing")

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    data <- radiator::tidy_wide(data = data, import.metadata = TRUE)
  }
  want <- c("MARKERS", "POP_ID", "INDIVIDUALS", "GT", "GT_BIN")
  data <- suppressWarnings(dplyr::select(data, dplyr::one_of(want)))

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

  data <- dplyr::arrange(data, POP_ID, INDIVIDUALS)

  if (is.factor(data$POP_ID)) {
    pop.levels <- levels(data$POP_ID)
  } else {
    pop.levels <- unique(data$POP_ID)
  }
  # Make sure that POP_ID and INDIVIDUALS are character
  # data <- dplyr::mutate_at(.tbl = data, .vars = c("POP_ID", "INDIVIDUALS"), .funs = as.character)
  data$INDIVIDUALS <- as.character(data$INDIVIDUALS)
  data$POP_ID <- as.character(data$POP_ID)

  # Isolate the strata
  # we convert pop_id to factor because adegenet does it automatically...
  pop.num <- unique(stringi::stri_detect_regex(str = pop.levels, pattern = "^[0-9]+$"))
  if (length(pop.num) == 1 && pop.num) pop.levels <- as.character(sort(as.numeric(pop.levels)))

  # When VCF data available
  if (tibble::has_name(data, "GT_BIN")) {
    data <- suppressWarnings(
      dplyr::select(.data = data, MARKERS, POP_ID, INDIVIDUALS, GT_BIN) %>%
        dplyr::mutate(A1 = abs(GT_BIN - 2)) %>%
        dplyr::rename(A2 = GT_BIN) %>%
        tidyr::gather(data = ., key = ALLELES, value = n, -c(INDIVIDUALS, POP_ID, MARKERS)) %>%
        dplyr::mutate(MARKERS_ALLELES = stringi::stri_join(MARKERS, ALLELES, sep = ".")) %>%
        dplyr::select(-MARKERS, -ALLELES) %>%
        dplyr::group_by(POP_ID, INDIVIDUALS) %>%
        tidyr::spread(data =., key = MARKERS_ALLELES, value = n) %>%
        dplyr::ungroup(.) %>%
        dplyr::mutate(POP_ID = factor(as.character(POP_ID))) %>%# xvalDapc doesn't accept pop as ordered factor
        dplyr::arrange(POP_ID, INDIVIDUALS))
  } else {
    missing.geno <- dplyr::ungroup(data) %>%
      dplyr::select(MARKERS, INDIVIDUALS, GT) %>%
      dplyr::filter(GT == "000000") %>%
      dplyr::select(MARKERS, INDIVIDUALS) #%>%  dplyr::mutate(MISSING = rep("blacklist", n()))

    data <- suppressWarnings(
      dplyr::ungroup(data) %>%
        dplyr::select(MARKERS, INDIVIDUALS, GT, POP_ID) %>%
        dplyr::filter(GT != "000000") %>%
        dplyr::mutate(
          A1 = stringi::stri_sub(str = GT, from = 1, to = 3),
          A2 = stringi::stri_sub(str = GT, from = 4, to = 6)
        ) %>%
        dplyr::select(-GT) %>%
        tidyr::gather(
          data = .,
          key = ALLELES,
          value = GT,
          -c(MARKERS, INDIVIDUALS, POP_ID)
        ) %>%
        dplyr::arrange(MARKERS, POP_ID, INDIVIDUALS, GT) %>%
        dplyr::count(x = ., POP_ID, INDIVIDUALS, MARKERS, GT) %>%
        dplyr::ungroup(.) %>%
        tidyr::complete(data = ., tidyr::nesting(INDIVIDUALS, POP_ID), tidyr::nesting(MARKERS, GT), fill = list(n = 0)) %>%
        dplyr::mutate(MARKERS_ALLELES = stringi::stri_join(MARKERS, GT, sep = ".")) %>%
        dplyr::anti_join(missing.geno, by = c("MARKERS", "INDIVIDUALS")) %>%
        dplyr::select(-MARKERS, -GT) %>%
        dplyr::mutate(POP_ID = factor(as.character(POP_ID))) %>%# xvalDapc doesn't accept pop as ordered factor
        dplyr::arrange(MARKERS_ALLELES, INDIVIDUALS) %>%
        dplyr::group_by(POP_ID, INDIVIDUALS) %>%
        tidyr::spread(data =., key = MARKERS_ALLELES, value = n) %>%
        dplyr::ungroup(.) %>%
        dplyr::arrange(INDIVIDUALS, POP_ID))
  }

  strata.genind <- dplyr::distinct(.data = data, INDIVIDUALS, POP_ID) %>%
    dplyr::mutate(POP_ID = factor(POP_ID, levels = pop.levels, ordered = FALSE),
                  INDIVIDUALS = factor(INDIVIDUALS))

  # genind arguments common to all data.type
  ind <- data$INDIVIDUALS
  pop <- data$POP_ID
  data <-  dplyr::ungroup(data) %>%
    dplyr::select(-c(INDIVIDUALS, POP_ID))
  suppressWarnings(rownames(data) <- ind)

  # genind constructor
  prevcall <- match.call()
  res <- adegenet::genind(
    tab = data,
    pop = pop,
    prevcall = prevcall,
    ploidy = 2,
    type = "codom",
    strata = strata.genind,
    hierarchy = NULL
  )
  data <- NULL
  return(res)
} # End write_genind
