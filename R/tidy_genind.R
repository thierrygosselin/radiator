#' @name tidy_genind
#' @title Tidy a genind object to a tidy dataframe
#' @description Tidy genind object to a tidy dataframe.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data A genind object in the global environment.

#' @export
#' @rdname tidy_genind

#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join
#' @importFrom stringi stri_replace_all_fixed stri_replace_all_regex stri_join
# @importFrom adegenet indNames
#' @importFrom tibble rownames_to_column data_frame
#' @importFrom tidyr gather unite

#' @references Jombart T (2008) adegenet: a R package for the multivariate
#' analysis of genetic markers. Bioinformatics, 24, 1403-1405.
#' @references Jombart T, Ahmed I (2011) adegenet 1.3-1:
#' new tools for the analysis of genome-wide SNP data.
#' Bioinformatics, 27, 3070-3071.


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


tidy_genind <- function(data) {

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("data argument required")
  if (class(data)[1] != "genind") stop("Input is not a genind object")
  if (is.null(data@pop)) stop("genind object requires population info in @pop")

  biallelic <- max(unique(data@loc.n.all)) == 2
  A2 <- unique(stringi::stri_detect_fixed(str = sample(colnames(data@tab), 100), pattern = ".A2"))

  if (biallelic) {

    # changed adegenet::indNames to rownames(data@tab) to lower dependencies
    if (A2) {
    data <- tibble::as_data_frame(data@tab) %>%
      tibble::add_column(.data = ., INDIVIDUALS = rownames(data@tab), .before = 1) %>%
      tibble::add_column(.data = ., POP_ID = data@pop) %>%
      dplyr::select(POP_ID, INDIVIDUALS, dplyr::ends_with(match = ".A2")) %>%
      tidyr::gather(data = ., key = MARKERS, value = GT_BIN, -c(POP_ID, INDIVIDUALS)) %>%
      dplyr::mutate(
        MARKERS = stringi::stri_replace_all_fixed(
          str = MARKERS, pattern = ".A2", replacement = "", vectorize_all = FALSE),
        GT_VCF = dplyr::if_else(GT_BIN == 0, "0/0",
                                dplyr::if_else(GT_BIN == 1, "0/1", "1/1"), missing = "./."),
        GT = dplyr::if_else(GT_BIN == 0, "001001", dplyr::if_else(GT_BIN == 1, "001002", "002002") , missing = "000000")
      )
    } else {
      data <- tibble::as_data_frame(data@tab) %>%
        tibble::add_column(.data = ., INDIVIDUALS = rownames(data@tab), .before = 1) %>%
        tibble::add_column(.data = ., POP_ID = data@pop) %>%
        dplyr::select(POP_ID, INDIVIDUALS, dplyr::ends_with(match = ".1")) %>%
        tidyr::gather(data = ., key = MARKERS, value = GT_BIN, -c(POP_ID, INDIVIDUALS)) %>%
        dplyr::mutate(
          MARKERS = stringi::stri_replace_all_fixed(
            str = MARKERS, pattern = ".1", replacement = "", vectorize_all = FALSE),
          GT_VCF = dplyr::if_else(GT_BIN == 0, "0/0",
                                  dplyr::if_else(GT_BIN == 1, "0/1", "1/1"), missing = "./."),
          GT = dplyr::if_else(GT_BIN == 0, "001001", dplyr::if_else(GT_BIN == 1, "001002", "002002") , missing = "000000")
        )
    }



  } else {
    data <- tibble::as_data_frame(data@tab) %>%
      tibble::add_column(.data = ., INDIVIDUALS = rownames(data@tab), .before = 1) %>%
      tibble::add_column(.data = ., POP_ID = data@pop) %>%
      tidyr::gather(data = ., key = MARKERS_ALLELES, value = COUNT, -c(POP_ID, INDIVIDUALS)) %>%
      dplyr::filter(COUNT > 0 | is.na(COUNT)) %>%
      tidyr::separate(data = ., col = MARKERS_ALLELES, into = c("MARKERS", "ALLELES"), sep = "\\.") %>%
      dplyr::mutate(
        ALLELES = as.numeric(factor(ALLELES)),
        ALLELES = stringi::stri_pad_left(str = ALLELES, pad = "0", width = 3)
        )

    # #If the genind was coded with allele 0, this will generate missing data with this code
    # allele.zero <- TRUE %in% stringi::stri_detect_regex(str = unique(data3$ALLELES), pattern = "^0$")
    # if (allele.zero) stop("alleles in this multiallelic genind were coded as 0 and won't work with this script")

    #Isolate missing genotype
    missing <- dplyr::filter(data, is.na(COUNT)) %>%
      dplyr::distinct(POP_ID, INDIVIDUALS, MARKERS) %>%
      dplyr::mutate(GT = rep("000000", n())) %>%
      dplyr::ungroup(.)

    #Isolate all Hom
    hom <- dplyr::filter(data, COUNT == 2) %>%
      dplyr::group_by(POP_ID, INDIVIDUALS, MARKERS) %>%
      dplyr::summarise(GT = stringi::stri_join(ALLELES, ALLELES, sep = "")) %>%
      dplyr::ungroup(.)

    #Isolate all Het and combine
    data <- dplyr::filter(data, COUNT != 2) %>%#this will also remove the missing
      dplyr::group_by(POP_ID, INDIVIDUALS, MARKERS) %>%
      dplyr::summarise(GT = stringi::stri_join(ALLELES, collapse = "")) %>%
      dplyr::ungroup(.) %>%
      dplyr::bind_rows(hom, missing)
    hom <- missing <- NULL
  }#End for multi-allelic
return(data)
} # End tidy_genind
