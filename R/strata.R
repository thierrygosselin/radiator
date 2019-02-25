# read_strata ----------------------------------------------------------------------------------
#' @name read_strata
#' @title read strata
#' @description Read a strata object or file.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param strata (path or object) The strata file or object.
#' Additional documentation is available in \code{\link{read_strata}}.
#' Use that function to whitelist/blacklist populations/individuals.
#' Option to set \code{pop.levels/pop.labels} is also available.

#' @param pop.id (logical) When \code{pop.id = TRUE}, the strata returns
#' the stratification colname \code{POP_ID}.
#' Default: \code{pop.id = FALSE}, returns \code{STRATA}.

#' @param pop.select (optional, string) Selected list of populations for
#' the analysis. e.g. \code{pop.select = c("QUE", "ONT")} to select \code{QUE}
#' and \code{ONT} population samples (out of 20 pops).
#' Default: \code{pop.select = NULL}

#' @param pop.levels (optional, string) This refers to the levels in a factor. In this
#' case, the id of the pop.
#' Use this argument to have the pop ordered your way instead of the default
#' alphabetical or numerical order. e.g. \code{pop.levels = c("QUE", "ONT", "ALB")}
#' instead of the default \code{pop.levels = c("ALB", "ONT", "QUE")}.
#' White spaces in population names are replaced by underscore.
#' Default: \code{pop.levels = NULL}.


#' @param pop.labels (optional, string) Use this argument to rename/relabel
#' your pop or combine your pop. e.g. To combine \code{"QUE"} and \code{"ONT"}
#' into a new pop called \code{"NEW"}:
#' (1) First, define the levels for your pop with \code{pop.levels} argument:
#' \code{pop.levels = c("QUE", "ONT", "ALB")}.
#' (2) then, use \code{pop.labels} argument:
#' \code{pop.labels = c("NEW", "NEW", "ALB")}.
#' To rename \code{"QUE"} to \code{"TAS"}:
#' \code{pop.labels = c("TAS", "ONT", "ALB")}.
#' Default: \code{pop.labels = NULL}. If you find this too complicated,
#' there is also the \code{strata} argument that can do the same thing,
#' see below.
#' White spaces in population names are replaced by underscore.

#' @inheritParams tidy_genomic_data
#' @inheritParams read_blacklist_id
#' @param keep.two (optional, logical) The output is limited to 2 columns:
#' \code{INDIVIDUALS, STRATA}.
#' Default: \code{keep.two = TRUE}.

#' @details the strata file used in radiator is a tab delimited file with
#' a minimum of 2 columns headers:
#' \code{INDIVIDUALS} and \code{STRATA}.
#' If a \code{strata} file is specified with all file formats that don't
#' require it, the strata argument will have precedence on the population
#' groupings used internally in those file formats. For file formats without
#' population/strata groupings (e.g. vcf, haplotype files) if no strata file is
#' provided, 1 pop/strata grouping will automatically be created.
#' For vcf and haplotypes file, the strata can also be used as a whitelist of id.
#' Samples not in the strata file will be discarded from the data set.
#' The \code{STRATA} column can be any hierarchical grouping.
#' To create a strata file see \code{\link[radiator]{individuals2strata}}.
#' If you have already run
#' \href{http://catchenlab.life.illinois.edu/stacks/}{stacks} on your data,
#' the strata file is similar to a stacks \emph{population map file},
#' make sure you
#' have the required column names (\code{INDIVIDUALS} and \code{STRATA}).
#' The strata column is cleaned of a white spaces that interfere with some
#' packages or codes: space is changed to an underscore \code{_}.

#' @seealso \code{\link{summary_strata}},
#' \code{\link{individuals2strata}}, \code{\link{change_pop_names}},
#' \code{\link{join_strata}}, \code{\link{generate_strata}}

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

read_strata <- function(
  strata,
  pop.id = FALSE,
  pop.levels = NULL,
  pop.labels = NULL,
  pop.select = NULL,
  blacklist.id = NULL,
  keep.two = TRUE,
  verbose = FALSE
) {
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")

  if (is.null(strata)) {
    return(res = NULL)
  } else {
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

    if (keep.two) strata  %<>% dplyr::select(INDIVIDUALS, STRATA)

    # clean....
    strata$INDIVIDUALS <- radiator::clean_ind_names(strata$INDIVIDUALS)
    strata$STRATA <- radiator::clean_pop_names(strata$STRATA)

    if (verbose) message("    Number of strata: ", length(unique(strata$STRATA)))
    if (verbose) message("    Number of individuals: ", length(unique(strata$INDIVIDUALS)))

    #blacklist.id ----------------------------------------------------------------
    blacklist.id <- read_blacklist_id(blacklist.id, verbose)
    if (!is.null(blacklist.id)){
      strata  %<>% dplyr::filter(!INDIVIDUALS %in% blacklist.id$INDIVIDUALS)
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
      strata  %<>% suppressWarnings(dplyr::filter(STRATA %in% pop.select))
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
      strata  %<>%
        dplyr::mutate(
          TARGET_ID = stringi::stri_trans_toupper(TARGET_ID),
          TARGET_ID = stringi::stri_replace_all_fixed(
            TARGET_ID, pattern = " ", replacement = "", vectorize_all = FALSE),
          TARGET_ID = clean_ind_names(TARGET_ID))
    }

    strata  %<>% dplyr::arrange(STRATA, INDIVIDUALS)

    if (isTRUE(pop.id)) strata %<>% dplyr::rename(POP_ID = STRATA)


    if (!is.null(blacklist.id) || !is.null(pop.select)) {
      strata.fn <- stringi::stri_join("strata_radiator_filtered_", file.date, ".tsv")
      readr::write_tsv(x = strata, path = strata.fn)
    }

    res = list(
      strata = strata,
      pop.levels = pop.levels,
      pop.labels = pop.labels,
      pop.select = pop.select,
      blacklist.id = blacklist.id)
  }
  return(res)
}#End read_strata

# Summary strata ---------------------------------------------------------------
#' @title summary of strata
#' @description Summarise the information of a strata file or object.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.
#' @param strata (path or object) The strata file or object.
#' @seealso \code{\link{read_strata}},
#' \code{\link{individuals2strata}}, \code{\link{change_pop_names}},
#' \code{\link{join_strata}}, \code{\link{generate_strata}}
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

  message("Number of populations: ", length(unique(strata$STRATA)))
  message("Number of individuals: ", length(unique(strata$INDIVIDUALS)))
  message("\nNumber of ind/pop:\n", stringi::stri_join(strata.stats$POP_IND, collapse ="\n"))
  message("\nNumber of duplicate id: ", duplicate.id)
}#End summary_strata


# individuals2strata------------------------------------------------------------
# Make strata file from individuals

#' @name individuals2strata
#' @title Create a strata file from a list of individuals
#' @description If your individuals have a consistent naming scheme
#' (e.g. SPECIES-POPULATION-MATURITY-YEAR-ID = CHI-QUE-ADU-2014-020),
#' use this function to rapidly create a strata file.
#' Several functions in \pkg{radiator} and \pkg{assigner} requires
#' a \code{strata} argument, i.e. a data frame with the individuals and
#' associated groupings. If you have already run
#' \href{http://catchenlab.life.illinois.edu/stacks/}{stacks} on your data,
#' the strata file is similar to a stacks `population map file`, make sure you
#' have the required column names  (\code{INDIVIDUALS} and \code{STRATA}).

#' @param data A file or data frame object with individuals in a column. The
#' column name is \code{INDIVIDUALS}.

#' @param strata.start (integer) The start of your strata id. See details for more info.

#' @param strata.end (integer) The end of your strata id. See details for more info.

#' @param filename (optional) The file name for the strata object if you
#' want to save it in the working directory.
#' Default: \code{filename = NULL}, the starta object is in the global
#' environment only (i.e. not written in the working directory).

#' @seealso \code{\link{read_strata}}, \code{\link{summary_strata}},
#' \code{\link{change_pop_names}},
#' \code{\link{join_strata}}, \code{\link{generate_strata}}

#' @return a strata object and file, if requested. The file is tab delimited
#' with 2 columns named:
#' \code{INDIVIDUALS} and \code{STRATA}.
#' The \code{STRATA} column can be any hierarchical grouping.


#' @details
#' \code{strata.start} and \code{strata.end}
#' The info must be found within the name of your individual sample. If not,
#' you'll have to create a strata file by hand, the old fashion way.
#' e.g. if your individuals are identified
#' in this form : SPECIES-POPULATION-MATURITY-YEAR-ID = CHI-QUE-ADU-2014-020,
#' then, to have the population id in the \code{STRATA} column,
#' \code{strata.start = 5} and \code{strata.end = 7}.
#' The \code{STRATA} column can be any hierarchical grouping.

#' @export
#' @rdname individuals2strata
#' @importFrom stringi stri_sub
#' @importFrom dplyr mutate
#' @importFrom readr write_tsv read_tsv
#' @importFrom tibble as_data_frame


#' @examples
#' \dontrun{
#' strata.abalone <- individuals2strata(
#' data = "individuals.abalone.tsv",
#' strata.start = 5,
#' strata.end = 7,
#' filename = "strata.abalone.tsv"
#' )
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

individuals2strata <- function(
  data,
  strata.start,
  strata.end,
  filename = NULL
) {


  # Checking for missing and/or default arguments ******************************
  if (missing(data)) rlang::abort("Input file missing")
  if (missing(strata.start)) rlang::abort("strata.start argument missing")
  if (missing(strata.end)) rlang::abort("strata.end argument missing")
  if (is.vector(data)) data <- readr::read_tsv(file = data)

  data <- tibble::as_data_frame(data) %>%
    dplyr::mutate(
      INDIVIDUALS =  as.character(INDIVIDUALS),
      STRATA = stringi::stri_sub(str = INDIVIDUALS, from = strata.start, to = strata.end)
    )


  # Write to working directory
  if (!is.null(filename)) {
    message("Writing the strata object to the working directory: \n", filename)
    readr::write_tsv(x = data, path = filename, col_names = TRUE)
  }

  return(data)
} # end individuals2strata




# change_pop_names--------------------------------------------------------------

#' @name change_pop_names
#' @title Transform into a factor the POP_ID column, change names and reorder the levels
#' @description Transform into a factor the POP_ID column, change names and
#' reorder the levels. If the data as \code{STRATA} column instead of a
#' \code{POP_ID} column, the function will change the column name.

#' @inheritParams tidy_genomic_data
#' @inheritParams read_strata
#' @seealso \code{\link{read_strata}}, \code{\link{summary_strata}},
#' \code{\link{individuals2strata}},
#' \code{\link{join_strata}}, \code{\link{generate_strata}}

#' @rdname change_pop_names
#' @export
#'
#' @importFrom stringi stri_replace_all_regex stri_join stri_replace_all_fixed
#' @importFrom dplyr n_distinct

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


change_pop_names <- function(data, pop.levels = NULL, pop.labels = NULL) {

  # checks ---------------------------------------------------------------------
  if (missing(data)) rlang::abort("Input file missing")

  # POP_ID in gsi_sim does not like spaces, we need to remove space in everything touching POP_ID...

  if (tibble::has_name(data, "STRATA") & !tibble::has_name(data, "POP_ID")) {
    data %<>% dplyr::rename(POP_ID = STRATA)
  }


  # removing spaces in data$POP_ID, pop.levels and pop.labels
  if (!is.null(pop.levels)) {
    if (is.null(pop.labels)) {
      pop.labels <-pop.levels <- clean_pop_names(pop.levels)
    }
    if (dplyr::n_distinct(data$POP_ID) != length(pop.levels)) {
      rlang::abort("The number of strata/POP_ID in the data is different than the number of pop.levels: check argument and data")
    }
  }

  if (!is.null(pop.labels)) {
    if (is.null(pop.levels)) rlang::abort("pop.levels is required if you use pop.labels")
    if (length(pop.labels) != length(pop.levels)) rlang::abort("pop.levels and pop.labels with different length: check arguments")
    pop.labels <- clean_pop_names(pop.labels)
  }

  # in the data
  data$POP_ID <- clean_pop_names(data$POP_ID)

  # convert POP_ID to factor and change names-----------------------------------

  if (is.null(pop.levels)) { # no pop.levels
    data$POP_ID <- factor(data$POP_ID)
  } else {# with pop.levels
    data$POP_ID <- factor(x = data$POP_ID, levels = pop.levels, ordered = FALSE)
    levels(data$POP_ID) <- pop.labels
  }
  data <- dplyr::arrange(data, POP_ID, INDIVIDUALS)
  return(data)
}# end function change_pop_names





# Check pop.levels and pop.labels --------------------------------------------

#' @name check_pop_levels
#' @title Check the use of pop.levels, pop.labels and pop.select arguments.
#' @description Check that pop.levels and pop.labels and pop.select arguments
#' are used correctly and that the values are cleaned for spaces.

#' @inheritParams tidy_genomic_data
#' @inheritParams read_strata

#' @rdname check_pop_levels
#' @export
#'
#' @importFrom stringi stri_replace_all_regex stri_join stri_replace_all_fixed
#' @importFrom dplyr n_distinct

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


check_pop_levels <- function(
  pop.levels = NULL,
  pop.labels = NULL,
  pop.select = NULL
) {

  # checks ---------------------------------------------------------------------
  # removing spaces in data$POP_ID, pop.levels and pop.labels
  if (!is.null(pop.levels) && is.null(pop.labels)) {
    pop.labels <-pop.levels <- clean_pop_names(pop.levels)
  }

  if (!is.null(pop.labels)) {
    if (is.null(pop.levels)) rlang::abort("pop.levels is required if you use pop.labels")
    if (length(pop.labels) != length(pop.levels)) {
      rlang::abort("pop.levels and pop.labels with different length: check arguments")
    }
    pop.labels <- clean_pop_names(pop.labels)
  }
  if (!is.null(pop.select)) pop.select <- clean_pop_names(pop.select)
  return(res = list(pop.levels = pop.levels, pop.labels = pop.labels, pop.select = pop.select))
}# end function change_pop_names



# join_strata ------------------------------------------------------------------

#' @name join_strata
#' @title Join the strata with the data
#' @description The function first filters individuals in data then include the
#' strata.
#' @param data A tidy dataset object.
#' Documented in \code{\link[radiator]{tidy_genomic_data}}.
#' @param strata A strata object.
#' Documented in \code{\link[radiator]{read_strata}}.
#' @return The data filtered by the strata by individuals.

#' @examples
#' \dontrun{
#' data <- radiator::join_strata(
#'     data = my_tidy_dataset_object,
#'     strata = my_strata_object)
#' }


#' @seealso \code{\link{read_strata}}, \code{\link{summary_strata}},
#' \code{\link{individuals2strata}}, \code{\link{change_pop_names}},
#' \code{\link{generate_strata}}



#' @rdname join_strata
#' @export

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

join_strata <- function(data, strata = NULL, verbose = TRUE) {
  if (is.null(strata)) return(data)
  if (verbose) message("Synchronizing data and strata...")
  if (rlang::has_name(data, "POP_ID")) data %<>% dplyr::select(-POP_ID)
  if (rlang::has_name(data, "STRATA")) data %<>% dplyr::select(-STRATA)
  strata %<>% dplyr::filter(INDIVIDUALS %in% data$INDIVIDUALS)
  if (nrow(strata) == 0) {
    rlang::abort("No more individuals in your data, check data and strata ID names...")
  }

  data %<>% dplyr::filter(INDIVIDUALS %in% strata$INDIVIDUALS)
  if (nrow(data) == 0) {
    rlang::abort("No more individuals in your data, check data and strata ID names...")
  }

  data %<>% dplyr::left_join(strata, by = "INDIVIDUALS")
  if (verbose) {
    if (rlang::has_name(data, "POP_ID")) {
      message("    Number of strat: ", length(unique(data$POP_ID)))
    }
    if (rlang::has_name(data, "STRATA")) {
      message("    Number of strata: ", length(unique(data$STRATA)))
    }
    message("    Number of individuals: ", length(unique(data$INDIVIDUALS)))
  }
  return(data)
}#End join_strata


# generate_strata ------------------------------------------------------------------

#' @name generate_strata
#' @title Generate strata object from the data
#' @description Generate strata object from the data \code{POP_ID} or \code{STRATA}
#' with the \code{INDIVIDUALS}.

#' @inheritParams join_strata
#' @seealso \code{\link{read_strata}}, \code{\link{summary_strata}},
#' \code{\link{individuals2strata}}, \code{\link{change_pop_names}},
#' \code{\link{join_strata}}

#' @rdname generate_strata
#' @export

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

generate_strata <- function(data) {
  data %<>% dplyr::ungroup(.)
  if (rlang::has_name(data, "POP_ID")) data %<>% dplyr::distinct(POP_ID, INDIVIDUALS)
  if (rlang::has_name(data, "STRATA")) data %<>% dplyr::distinct(STRATA, INDIVIDUALS)
  return(data)
}#End join_strata


# strata_haplo -----------------------------------------------------------------

#' @title strata_haplo
#' @description Manage strata
#' @rdname strata_haplo
#' @keywords internal
#' @export
strata_haplo <- function(strata = NULL, data = NULL, blacklist.id = NULL) {

  if (is.null(strata)) {
    message("No strata file provided")
    message("    generating a strata with 1 grouping")
    if (is.null(data)) rlang::abort("data required to generate strata")
    strata.df <- readr::read_tsv(
      file = data,
      n_max = 1,
      na = "-",
      col_names = FALSE,
      col_types = readr::cols(.default = readr::col_character())) %>%
      tidyr::gather(data = .,key = DELETE, value = INDIVIDUALS) %>%
      dplyr::mutate(INDIVIDUALS = clean_ind_names(INDIVIDUALS)) %>%
      dplyr::select(-DELETE) %>%
      dplyr::filter(!INDIVIDUALS %in% c("Catalog ID", "Cnt")) %>%
      dplyr::distinct(INDIVIDUALS) %>%
      dplyr::mutate(STRATA = rep("pop1", n()))
  } else {
    if (is.vector(strata)) {
      suppressMessages(
        strata.df <- readr::read_tsv(
          file = strata, col_names = TRUE,
          # col_types = col.types
          col_types = readr::cols(.default = readr::col_character())
        ))
    } else {
      strata.df <- strata
    }
  }

  colnames(strata.df) <- stringi::stri_replace_all_fixed(
    str = colnames(strata.df),
    pattern = "STRATA",
    replacement = "POP_ID",
    vectorize_all = FALSE
  )
  # Remove potential whitespace in pop_id
  strata.df$POP_ID <- clean_pop_names(strata.df$POP_ID)
  colnames.strata <- colnames(strata.df)

  # clean ids
  strata.df$INDIVIDUALS <- clean_ind_names(strata.df$INDIVIDUALS)

  # filtering the strata if blacklist id available
  if (!is.null(blacklist.id)) {
    strata.df <- dplyr::anti_join(x = strata.df, y = blacklist.id, by = "INDIVIDUALS")
  }

  strata.df <- dplyr::distinct(strata.df, POP_ID, INDIVIDUALS)
  return(strata.df)
}#End strata_haplo
