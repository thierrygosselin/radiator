#' @name tidy_genlight
#' @title Tidy a genlight object to a tidy dataframe
#' @description Tidy genlight object to a tidy dataframe.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data A genlight object in the global environment.

#' @export
#' @rdname tidy_genlight

#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join
#' @importFrom stringi stri_replace_all_fixed stri_replace_all_regex stri_join
#' @importFrom adegenet indNames pop chromosome locNames position
#' @importFrom tibble rownames_to_column data_frame
#' @importFrom tidyr gather unite

#' @references Jombart T (2008) adegenet: a R package for the multivariate
#' analysis of genetic markers. Bioinformatics, 24, 1403-1405.
#' @references Jombart T, Ahmed I (2011) adegenet 1.3-1:
#' new tools for the analysis of genome-wide SNP data.
#' Bioinformatics, 27, 3070-3071.


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


tidy_genlight <- function(data) {

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("data argument required")
  if (class(data)[1] != "genlight") stop("Input is not a genlight object")


  # Import data ---------------------------------------------------------------
  strata.df <- tibble::data_frame(
    INDIVIDUALS = adegenet::indNames(data),
    POP_ID = adegenet::pop(data))

  markers <- tibble::data_frame(
    CHROM = adegenet::chromosome(data),
    LOCUS = adegenet::locNames(data),
    POS = adegenet::position(data)
  ) %>%
    tidyr::unite(data = ., col = MARKERS, CHROM, LOCUS, POS, sep = "__", remove = FALSE) %>%
    dplyr::select(MARKERS, CHROM, LOCUS, POS)



  data <- as.data.frame(data)
  colnames(data) <- markers$MARKERS
  data <- tibble::rownames_to_column(df = data, var = "INDIVIDUALS") %>%
    tidyr::gather(data = ., key = MARKERS, value = GT_BIN, -INDIVIDUALS) %>%
    dplyr::full_join(strata.df, by =  "INDIVIDUALS") %>%
    dplyr::full_join(markers, by = "MARKERS") %>%
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
      )
    ) %>%
    dplyr::select(MARKERS, CHROM, LOCUS, POS, POP_ID, INDIVIDUALS, GT_VCF, GT_BIN, GT)

  return(data)
} # End tidy_genlight
