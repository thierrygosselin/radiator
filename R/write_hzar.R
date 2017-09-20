# write a HZAR file from a tidy data frame

#' @name write_hzar
#' @title Write a HZAR file from a tidy data frame.

#' @description Write a HZAR file (Derryberry et al. 2013), from a tidy data frame.
#' Used internally in
#' \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.

#' @param distances (optional) A file with 2 columns, POP_ID The distance information per populations.
#' With default: \code{distances = NULL}, the column is left empty.


#' @inheritParams tidy_genomic_data

#' @param filename (optional) The file name prefix for the HZAR file
#' written to the working directory. With default: \code{filename = NULL},
#' the date and time is appended to \code{radiator_hzar_}.

#' @details \emph{Integrated filters:} Only markers found in common between
#' populations are used and monomorphic markers are automatically removed
#' before generating HZAR file.


#' @return A HZAR file is written in the working directory.

#' @export
#' @rdname write_hzar
#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join
#' @importFrom stringi stri_join stri_replace_all_fixed stri_sub
#' @importFrom tidyr spread unite complete nesting separate gather
#' @importFrom readr write_delim write_tsv write_file write_lines write_csv
#' @importFrom parallel detectCores

#' @references Derryberry EP, Derryberry GE, Maley JM, Brumfield RT.
#' hzar: hybrid zone analysis using an R software package.
#' Molecular Ecology Resources. 2013;14: 652-663. doi:10.1111/1755-0998.12209

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

write_hzar <- function(
  data,
  distances = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1
) {
  opt.change <- getOption("width")
  options(width = 70)


  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file is missing")

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

  # Keeping common markers -----------------------------------------------------
  input <- radiator::keep_common_markers(data = input, verbose = TRUE)$input

  # Removing monomorphic markers -----------------------------------------------
  input <- radiator::discard_monomorphic_markers(data = input, verbose = TRUE)$input

  # detect biallelic markers ---------------------------------------------------
  # biallelic <- radiator::detect_biallelic_markers(data = input)
  n.ind <- dplyr::n_distinct(input$INDIVIDUALS)
  n.pop <- dplyr::n_distinct(input$POP_ID)
  n.markers <- dplyr::n_distinct(input$MARKERS)

  if (is.null(distances)) {
    distances <- dplyr::distinct(input, POP_ID) %>%
      dplyr::mutate(Distance = rep("", n()))
  } else {
    if (is.vector(distances)) {
      suppressMessages(distances <- readr::read_tsv(distances, col_names = TRUE))
    }
    # check same pop_id
    pop.distance <- dplyr::distinct(distances, POP_ID)
    if (!identical(unique(input$POP_ID), pop.distance$POP_ID)) {
      stop("Populations in `distances` file doesn't match populations in `input` file")
    }
  }
  message("Generating HZAR file...")
  output <- dplyr::filter(input, GT != "000000") %>%
    dplyr::left_join(
      dplyr::distinct(input, MARKERS) %>%
        dplyr::mutate(
          SPLIT_VEC = dplyr::ntile(x = 1:nrow(.), n = parallel.core * 3))
      , by = "MARKERS") %>%
    split(x = ., f = .$SPLIT_VEC) %>%
    .radiator_parallel(
      X = .,
      FUN = generate_hzar,
      mc.cores = parallel.core
    ) %>%
    dplyr::bind_rows(.) %>%
    dplyr::arrange(MARKERS, POP_ID, VALUE) %>%
    dplyr::select(-MARKERS) %>%
    dplyr::group_by(POP_ID) %>%
    tidyr::spread(data = ., key = GROUP, value = VALUE, fill = 0) %>%
    dplyr::ungroup(.) %>%
    dplyr::left_join(distances, by = "POP_ID") %>%
    dplyr::select(POP_ID, Distance, dplyr::everything(.)) %>%
    dplyr::rename(Population = POP_ID)

  input <- NULL

  # writting file to directory  ------------------------------------------------
  # Filename: date and time to have unique filenaming
  file.date <- stringi::stri_replace_all_fixed(
    Sys.time(),
    pattern = " EDT", replacement = "") %>%
    stringi::stri_replace_all_fixed(
      str = .,
      pattern = c("-", " ", ":"), replacement = c("", "@", ""),
      vectorize_all = FALSE) %>%
    stringi::stri_sub(str = ., from = 1, to = 13)

  if (is.null(filename)) {
    filename <- stringi::stri_join("radiator_hzar_", file.date, ".csv")
  } else {
    filename.problem <- file.exists(filename)
    if (filename.problem) {
      filename <- stringi::stri_join(filename, "_hzar_", file.date, ".csv")
    } else {
      filename <- stringi::stri_join(filename, "_hzar", ".csv")
    }
  }
  header.line <- stringi::stri_join("# HZAR v.0.2-5 file generated with radiator v.", packageVersion(pkg = "radiator"), " ", file.date)
  readr::write_lines(x = header.line, path = filename)
  readr::write_csv(x = output, path = filename, append = TRUE, col_names = TRUE)

  message("Writting HZAR file, (", filename,") with:
    Number of populations: ", n.pop, "\n    Number of individuals: ", n.ind,
    "\n    Number of markers: ", n.markers)


  res <- stringi::stri_join("Check your working directory for file: ", filename)
  return(res)
}
# Internal nested Function -----------------------------------------------------

#' @title generate_hzar
#' @description function to generate hzar function per groups of markes (to run in parallel)
#' @rdname generate_hzar
#' @keywords internal
#' @export

generate_hzar <- function(x) {
  freq.info <- x %>%
    dplyr::mutate(
      A1 = stringi::stri_sub(GT, 1, 3),
      A2 = stringi::stri_sub(GT, 4,6)
    ) %>%
    dplyr::select(MARKERS, POP_ID, INDIVIDUALS, A1, A2) %>%
    tidyr::gather(data = ., key = ALLELES_GROUP, value = ALLELES, -c(MARKERS, INDIVIDUALS, POP_ID)) %>%
    dplyr::select(-ALLELES_GROUP) %>%
    dplyr::group_by(MARKERS, POP_ID, ALLELES) %>%
    dplyr::tally(.) %>%
    dplyr::group_by(MARKERS, POP_ID) %>%
    dplyr::mutate(
      NN = sum(n),
      N = NN / 2) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(FREQ = n / NN) %>%
    dplyr::select(-n, -NN)

  sample.n.info <- dplyr::distinct(freq.info, MARKERS, POP_ID, N) %>%
    dplyr::mutate(GROUP = stringi::stri_join(MARKERS, ".N")) %>%
    dplyr::select(MARKERS, POP_ID, GROUP, VALUE = N)

  freq.info <- dplyr::select(freq.info, MARKERS, POP_ID, ALLELES,FREQ) %>%
    dplyr::mutate(GROUP = stringi::stri_join(MARKERS, ALLELES, sep = ".")) %>%
    dplyr::select(MARKERS, POP_ID, GROUP, VALUE = FREQ) %>%
    dplyr::bind_rows(sample.n.info)
  sample.n.info <- NULL
  return(freq.info)
}#End generate_hzar


