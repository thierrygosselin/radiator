# tidy_gtypes ------------------------------------------------------------------

#' @name tidy_gtypes
#' @title Tidy a gtypes object to a tidy dataframe
#' @description Transform a [strataG gtypes](https://github.com/EricArcher/strataG) object from
#' to a tidy dataframe.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data A gtypes object (>= v.2.0.2) in the global environment.

#' @export
#' @rdname tidy_gtypes

#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join
#' @importFrom stringi stri_replace_all_fixed stri_replace_all_regex stri_join stri_pad_left
#' @importFrom tibble rownames_to_column data_frame as_data_frame
#' @importFrom tidyr gather unite

#' @references Archer FI, Adams PE, Schneiders BB.
#' strataG: An r package for manipulating, summarizing and analysing population
#' genetic data.
#' Molecular Ecology Resources. 2017; 17: 5-11. doi:10.1111/1755-0998.12559


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

tidy_gtypes <- function(data) {

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file necessary...")
  if (class(data) != "gtypes") stop("Input is not a genlight object")

  # import ---------------------------------------------------------------------
  input <- suppressWarnings(
    tibble::as_data_frame(data@data) %>%
      dplyr::rename(INDIVIDUALS = ids, POP_ID = strata) %>%
      dplyr::mutate(ALLELES = rep(c("A1", "A2"), n() / 2)) %>%
      tidyr::gather(data = ., key = MARKERS, value = GT, -c(POP_ID, INDIVIDUALS, ALLELES)))

  # detect stratg genotype coding ----------------------------------------------
  # For GT = c("A", "C", "G", "T")
  gt.format <- sort(unique(input$GT))

  if (unique(gt.format %in% c("A", "C", "G", "T", NA))) {
    input$GT <- stringi::stri_replace_all_regex(
      str = input$GT,
      pattern = c("A", "C", "G", "T"),
      replacement = c("001", "002", "003", "004"),
      vectorize_all = FALSE
    )
  }

  # For GT = c("1", "2")
  if (unique(gt.format %in% c("1", "2", NA))) {
    input$GT <- stringi::stri_pad_left(str = input$GT, pad = "0", width = 3)
  }

  # For GT coded with only 1 number
  # gtypes.number <- unique(stringi::stri_count_boundaries(str = input$GT))
  # unique(stringi::stri_count_boundaries(str = test))

  # prep tidy ------------------------------------------------------------------
  input <- input %>%
    dplyr::mutate(
      GT = replace(GT, which(is.na(GT)), "000"),
      POP_ID = as.character(POP_ID)) %>%
    dplyr::group_by(POP_ID, INDIVIDUALS, MARKERS) %>%
    tidyr::spread(data = ., key = ALLELES, value = GT) %>%
    dplyr::ungroup(.) %>%
    tidyr::unite(data = ., col = GT, A1, A2, sep = "") %>%
    dplyr::select(POP_ID, INDIVIDUALS, MARKERS, GT)
  return(input)
}#End tidy_gtypes


# write_gtypes -----------------------------------------------------------------

#' @name write_gtypes
#' @title Write a \href{https://github.com/EricArcher/}{strataG} object from a tidy data frame

#' @description Write a\href{https://github.com/EricArcher/}{strataG} object from a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.

#' @return An object of the class \href{https://github.com/EricArcher/}{strataG} is returned.



#' @export
#' @rdname write_gtypes
# @import strataG

#' @importFrom tidyr gather
#' @importFrom methods new
#' @importFrom stringi stri_replace_all_fixed stri_sub
#' @importFrom dplyr select arrange rename mutate
#' @importFrom data.table dcast.data.table as.data.table
#' @importFrom purrr safely
#' @importFrom utils installed.packages


#' @seealso \href{https://github.com/EricArcher/}{strataG}

#' @references Archer FI, Adams PE, Schneiders BB.
#' strataG: An r package for manipulating, summarizing and analysing population
#' genetic data.
#' Molecular Ecology Resources. 2017; 17: 5-11. doi:10.1111/1755-0998.12559

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

write_gtypes <- function(data) {
  # Check that strataG is installed --------------------------------------------
  if (!"strataG" %in% utils::installed.packages()[,"Package"]) {
    rlang::abort("Please install strataG for this output option:\n
                 devtools::install_github('ericarcher/strataG', build_vignettes = TRUE)")
  }
  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("Input file missing")

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    data <- radiator::tidy_wide(data = data, import.metadata = TRUE)
  }


  data <- dplyr::select(.data = data, POP_ID, INDIVIDUALS, MARKERS, GT) %>%
    dplyr::arrange(MARKERS, POP_ID, INDIVIDUALS) %>%
    dplyr::mutate(
      GT = replace(GT, which(GT == "000000"), NA),
      POP_ID = as.character(POP_ID),
      `1` = stringi::stri_sub(str = GT, from = 1, to = 3), # most of the time: faster than tidyr::separate
      `2` = stringi::stri_sub(str = GT, from = 4, to = 6)
    ) %>%
    dplyr::select(-GT) %>%
    tidyr::gather(
      data = .,
      key = ALLELES,
      value = GT,
      -c(MARKERS, INDIVIDUALS, POP_ID)
    ) %>%
    data.table::as.data.table(.) %>%
    data.table::dcast.data.table(
      data = .,
      formula = POP_ID + INDIVIDUALS ~ MARKERS + ALLELES,
      value.var = "GT",
      sep = "."
    ) %>%
    tibble::as_data_frame(.)

  # dcast is slightly faster using microbenchmark then unite and spread
  # tidyr::unite(data = ., MARKERS_ALLELES, MARKERS, ALLELES, sep = ".") %>%
  # tidyr::spread(data = ., key = MARKERS_ALLELES, value = GT)

  safe_gtypes <-  purrr::safely(.f = methods::new)

  res <- suppressWarnings(
    safe_gtypes(
      "gtypes",
      gen.data = data[, -(1:2)],
      ploidy = 2,
      ind.names = data$INDIVIDUALS,
      strata = data$POP_ID,
      schemes = NULL,
      sequences = NULL,
      description = NULL,
      other = NULL
    )
  )

  if (is.null(res$error)) {
    res <- res$result
  } else {
    rlang::abort("strataG package must be installed and loaded: library('strataG')")
  }

  return(res)
  }# End write_gtypes
