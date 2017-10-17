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
    input <- radiator::tidy_wide(data = data, import.metadata = TRUE)
  } else {
    input <- data
  }

  # check genotype column naming
  colnames(input) <- stringi::stri_replace_all_fixed(
    str = colnames(input),
    pattern = "GENOTYPE",
    replacement = "GT",
    vectorize_all = FALSE
  )

  # necessary steps to make sure we work with unique markers and not duplicated LOCUS
  if (tibble::has_name(input, "LOCUS") && !tibble::has_name(input, "MARKERS")) {
    input <- dplyr::rename(.data = input, MARKERS = LOCUS)
  }

  strata.genind <- dplyr::distinct(.data = input, INDIVIDUALS, POP_ID)

  # When VCF data available
  if (tibble::has_name(input, "GT_VCF")) {
    input <- dplyr::select(.data = input, MARKERS, POP_ID, INDIVIDUALS, GT_VCF) %>%
      dplyr::mutate(
        A1_A2 = stringi::stri_replace_all_fixed(
          str = GT_VCF,
          pattern = c("0/0", "1/1", "0/1", "1/0", "./."),
          replacement = c("2_0", "0_2", "1_1", "1_1", NA),
          vectorize_all = FALSE
        )
      ) %>%
      dplyr::mutate(POP_ID = factor(as.character(POP_ID))) %>%# xvalDapc doesn't accept pop as ordered factor
      dplyr::mutate(
        A1 = stringi::stri_sub(str = A1_A2, from = 1, to = 1),
        A2 = stringi::stri_sub(str = A1_A2, from = 3, to = 3)
      ) %>%
      dplyr::select(-GT_VCF, -A1_A2) %>%
      tidyr::gather(data = ., key = ALLELES, value = n, -c(INDIVIDUALS, POP_ID, MARKERS)) %>%
      dplyr::mutate(MARKERS_ALLELES = stringi::stri_join(MARKERS, ALLELES, sep = ".")) %>%
      dplyr::select(-MARKERS, -ALLELES) %>%
      dplyr::group_by(POP_ID, INDIVIDUALS) %>%
      tidyr::spread(data =., key = MARKERS_ALLELES, value = n) %>%
      dplyr::ungroup(.)
  } else {

    missing.geno <- dplyr::ungroup(input) %>%
      dplyr::select(MARKERS, INDIVIDUALS, GT) %>%
      dplyr::filter(GT == "000000") %>%
      dplyr::select(MARKERS, INDIVIDUALS) %>%
      dplyr::mutate(MISSING = rep("blacklist", n()))

    input <- dplyr::ungroup(input) %>%
      dplyr::select(MARKERS, INDIVIDUALS, GT) %>%
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
        -c(MARKERS, INDIVIDUALS)
      ) %>%
      dplyr::arrange(MARKERS, INDIVIDUALS, GT) %>%
      dplyr::count(x = ., INDIVIDUALS, MARKERS, GT) %>%
      dplyr::ungroup(.) %>%
      tidyr::complete(data = ., INDIVIDUALS, tidyr::nesting(MARKERS, GT), fill = list(n = 0)) %>%
      dplyr::anti_join(missing.geno, by = c("MARKERS", "INDIVIDUALS")) %>%
      dplyr::mutate(MARKERS_ALLELES = stringi::stri_join(MARKERS, GT, sep = ".")) %>%
      dplyr::select(-MARKERS, -GT) %>%
      dplyr::right_join(strata.genind, by = "INDIVIDUALS") %>%#include strata
      dplyr::mutate(POP_ID = factor(as.character(POP_ID))) %>%# xvalDapc doesn't accept pop as ordered factor
      dplyr::arrange(MARKERS_ALLELES, INDIVIDUALS) %>%
      dplyr::group_by(POP_ID, INDIVIDUALS) %>%
      tidyr::spread(data =., key = MARKERS_ALLELES, value = n) %>%
      dplyr::ungroup(.)
   }

  # genind arguments common to all data.type
  ind <- input$INDIVIDUALS
  pop <- input$POP_ID
  input <-  dplyr::ungroup(input) %>% dplyr::select(-c(INDIVIDUALS, POP_ID))
  suppressWarnings(rownames(input) <- ind)

  # genind constructor
  prevcall <- match.call()
  res <- adegenet::genind(
    tab = input,
    pop = pop,
    prevcall = prevcall,
    ploidy = 2,
    type = "codom",
    strata = strata.genind,
    hierarchy = NULL
  )

  return(res)
} # End write_genind
