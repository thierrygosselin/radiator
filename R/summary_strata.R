#' @title read strata
#' @description Read a strata object or file.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.
#' @param strata (path or object) The strata file or object.
#' @param pop.id (logical) When \code{pop.id = TRUE}, the strata returns
#' the stratification colname \code{POP_ID}.
#' Default: \code{pop.id = FALSE}, returns \code{STRATA}.
#' @inheritParams tidy_genomic_data
#' @param keep.two (optional, logical) The output is limited to 2 columns:
#' \code{INDIVIDUALS, STRATA}.
#' Default: \code{keep.two = TRUE}.
#' @rdname read_strata
#' @export
#' @return
#' \enumerate{
#' \item Number of strata/populations
#' \item Number of individuals
#' \item Number of individuals per populations
#' \item Number of duplicate ids.
#' }
#' @examples
#' \dontrun{
#' radiator::read_strata(strata)
#' }

read_strata <- function(strata, pop.id = FALSE,
                        pop.levels = NULL, pop.labels = NULL,
                        pop.select = NULL, blacklist.id = NULL,
                        keep.two = TRUE, verbose = FALSE) {
  if (verbose) message("Analyzing strata file")
  if (is.vector(strata)) {
    strata <- readr::read_tsv(
      file = strata,
      col_types = readr::cols(.default = readr::col_character()))
  }

  if (tibble::has_name(strata, "POP_ID") && !tibble::has_name(strata, "STRATA")) {
    colnames(strata) <- stringi::stri_replace_all_fixed(
      colnames(strata), "POP_ID", "STRATA",
      vectorize_all = FALSE)
  }

  if (keep.two) strata <- dplyr::select(strata, INDIVIDUALS, STRATA)

  # clean....
  strata$INDIVIDUALS <- radiator::clean_ind_names(strata$INDIVIDUALS)
  strata$STRATA <- radiator::clean_pop_names(strata$STRATA)

  if (verbose) message("Number of individuals: ", dplyr::n_distinct(strata$INDIVIDUALS))
  if (verbose) message("Number of strata: ", dplyr::n_distinct(strata$STRATA))
  #blacklist.id ----------------------------------------------------------------
  if (!is.null(blacklist.id)){
    if (is.vector(blacklist.id)) {
      blacklist.id <- suppressMessages(
        readr::read_tsv(file = blacklist.id,
                        col_names = TRUE,
                        col_types = "c",
                        trim_ws = TRUE))
    }
    blacklist.id$INDIVIDUALS <- clean_ind_names(blacklist.id$INDIVIDUALS)


    # remove potential duplicate id
    dup <- dplyr::distinct(.data = blacklist.id, INDIVIDUALS)
    blacklist.id.dup <- nrow(blacklist.id) - nrow(dup)
    if (blacklist.id.dup >1) {
      message("Duplicate id's in blacklist: ", blacklist.id.dup)
      blacklist.id <- dup
    }
    dup <- blacklist.id.dup <- NULL
    n.ind.blacklist <- length(blacklist.id$INDIVIDUALS)
    if (verbose) message("\nNumber of individuals in blacklist: ", n.ind.blacklist, " ind.")
    n.ind.blacklisted <- length(strata$INDIVIDUALS %in% blacklist.id$INDIVIDUALS)
    strata <- dplyr::filter(strata, !INDIVIDUALS %in% blacklist.id$INDIVIDUALS)
    if (verbose) message("\nBlacklisted individuals: ", n.ind.blacklisted, " ind.")
  }


  # manage levels, labels and pop.select ---------------------------------------
  check <- check_pop_levels(pop.levels = pop.levels,
                            pop.labels = pop.labels,
                            pop.select = pop.select)
  pop.levels <- check$pop.levels
  pop.labels <- check$pop.labels
  pop.select <- check$pop.select

  if (!is.null(pop.select)) {
    n.pop.new <- length(pop.select)
    if (verbose) message("\nPopulations/strata selected: ", stringi::stri_join(pop.select, collapse = ", "), " (", n.pop.new," pops)")
    strata <- suppressWarnings(dplyr::filter(strata, STRATA %in% pop.select))
  }


  if (is.null(pop.levels)) { # no pop.levels
    strata$STRATA <- factor(strata$STRATA)
    pop.levels <- pop.labels <- unique(strata$STRATA)
  } else {# with pop.levels
    n.pop <- unique(strata$STRATA)
    if (length(n.pop) != length(pop.levels)) {
      if (verbose) message("pop.levels and unique STRATA/POP_ID have different length")
      if (verbose) message("    using unique STRATA/POP_ID names to replace pop.levels")
      pop.levels <- pop.labels <- unique(strata$STRATA)
    }
    strata$STRATA <- factor(
      x = strata$STRATA,
      levels = pop.levels,
      labels = pop.labels,
      ordered = FALSE)
  }

  if (!is.null(pop.select) || !is.null(blacklist.id)) {
    if (is.factor(pop.levels)) pop.levels <- droplevels(pop.levels)
    if (is.factor(pop.labels)) pop.labels <- droplevels(pop.labels)
  }

  # If dart file manage TARGET_ID ----------------------------------------------
  if (tibble::has_name(strata, "TARGET_ID")) {
    strata <- strata %>%
      dplyr::mutate(
        TARGET_ID = stringi::stri_trans_toupper(TARGET_ID),
        TARGET_ID = stringi::stri_replace_all_fixed(
          TARGET_ID, pattern = " ", replacement = "", vectorize_all = FALSE),
        TARGET_ID = clean_ind_names(TARGET_ID))
  }

  strata <- dplyr::arrange(strata, STRATA, INDIVIDUALS)

  if (isTRUE(pop.id)) strata <- dplyr::rename(strata, POP_ID = STRATA)

  return(
    res = list(
      strata = strata,
      pop.levels = pop.levels,
      pop.labels = pop.labels,
      pop.select = pop.select,
      blacklist.id = blacklist.id)
  )
}#End read_strata


#' @title summary of strata
#' @description Summarise the information of a strata file or object.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.
#' @param strata (path or object) The strata file or object.
#' @rdname summary_strata
#' @export
#' @return
#' \enumerate{
#' \item Number of strata/populations
#' \item Number of individuals
#' \item Number of individuals per populations
#' \item Number of duplicate ids.
#' }
#' @examples
#' \dontrun{
#' radiator::summary_strata(strata)
#' }
summary_strata <- function(strata) {

  strata <- radiator::read_strata(strata = strata)$strata

  strata.stats <- strata %>%
    dplyr::group_by(STRATA) %>%
    dplyr::tally(.) %>%
    dplyr::mutate(POP_IND = stringi::stri_join(STRATA, n, sep = " = "))

  duplicate.id <- nrow(strata) - length(unique(strata$INDIVIDUALS))

  message("Number of populations: ", dplyr::n_distinct(strata$STRATA))
  message("Number of individuals: ", dplyr::n_distinct(strata$INDIVIDUALS))
  message("\nNumber of ind/pop:\n", stringi::stri_join(strata.stats$POP_IND, collapse ="\n"))
  message("\nNumber of duplicate id: ", duplicate.id)
}#End summary_strata
