# tidy_genepop------------------------------------------------------------------
#' @name tidy_genepop

#' @title Import genepop file and convert to a tidy dataframe

#' @description Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.
#' The function \code{tidy_genepop} reads a file in the
#' \href{https://genepop.curtin.edu.au/help_input.html}{genepop format}
#' (see details and note for convention) and output a data frame in wide or long/tidy format.
#'
#' To manipulate and prune the dataset prior to tidying, use the functions
#' \code{\link[radiator]{tidy_genomic_data}} and
#' \code{\link[radiator]{genomic_converter}}, that uses blacklist and whitelist along
#' several other filtering options.

#' @param data A \href{https://genepop.curtin.edu.au/help_input.html}{genepop}
#' filename with extension \code{.gen}.

#' @param strata (optional) A tab delimited file with 2 columns. Header:
#' \code{INDIVIDUALS} and \code{STRATA}.
#' The \code{STRATA} column can be any hierarchical grouping.
#' To create a strata file see \code{\link[radiator]{individuals2strata}}.
#' Default: \code{strata = NULL}.

#' @param tidy (optional, logical) With \code{tidy = FALSE},
#' the markers are the variables and the genotypes the observations (wide format).
#' With the default: \code{tidy = TRUE}, markers and genotypes are variables
#' with their own columns (long format).

#' @param filename (optional) The file name for the tidy data frame
#' written to the working directory.
#' Default: \code{filename = NULL}, the tidy data is
#' in the global environment only (i.e. not written in the working directory).

#' @return The output in your global environment is a wide or long/tidy data frame.
#' If \code{filename} is provided, the wide or long/tidy data frame is also
#' written to the working directory.

#' @details \href{https://genepop.curtin.edu.au/help_input.html}{genepop format}
#' \enumerate{
#' \item \strong{First line:} This line is used to store information about
#' your data, any characters are allowed.
#' This line is not kept inside \code{tidy_genepop}.
#' \item \strong{Second line:} 2 options i) wide: all the locus are stored on this same
#' line with comma (or a comma+space) separator;
#' ii) long: the name of the first locus
#' (the remaining locus are on separate and subsequent rows with the long format:
#' not recommended with genomic datasets with thousands of markers...).
#'
#' The remaining lines are blocks of population and genotypes,
#' \href{https://genepop.curtin.edu.au/help_input.html}{for genepop format examples}.
#' \item \strong{population identifier:} The population block are separated by
#' the word: \code{POP}, or \code{Pop} or \code{pop}. Flavors of the genepop
#' software uses the first or the last identifier of every sub-population,
#' in all output files, to name populations
#' \href{https://genepop.curtin.edu.au/help_input.html}{(more info)}.
#' This is not very convenient for population naming and prone to errors,
#' this is where the \code{strata} argument inside \code{tidy_genepop}
#' (described above) becomes handy.
#' \item \strong{individual identifier:} After the population identifier
#' the individuals that belong to the same population are found subsequently on
#' separate lines. The genepop format specify you can use any character,
#' including a blank space or tab. Spaces are allowed in the identifier names.
#' You may leave it blank except for a comma if you wish. The comma between
#' the individual identifier and the list of genotypes is required.
#' For good naming habit however the function \code{tidy_genepop} will
#' replace \code{"_", ":"} with \code{"-"}, the comma \code{","} along any white
#' space characters defined as \code{"\t", "\n", "\f", "\r", "\p{Z}"} found in the individual
#' name will be trimmed.
#' \item \strong{genotypes:} For each locus, genotypes are separated by one or
#' more blank spaces or tab. 0101 indicates that this individual is homozygous
#' for the 01 allele at the first locus.
#' An alternative input format exists, where each allele is coded by three digits
#' (instead of two as described above). However, the total number of
#' different alleles, for each locus, should not be higher than 99.
#' Missing data is coded with zeros: \code{0000 or 000000}.
#' }

#' @note
#' \href{https://genepop.curtin.edu.au/help_input.html}{genepop format notes:}
#' \itemize{
#' \item No constraint on blanks separating the various fields.
#' \item tabs or spaces allowed.
#' \item Loci names can appear on separate lines (long), or on one line (wide)
#' if separated by commas.
#' \item Individual identifier may have blanks but must end with a comma.
#' \item Alleles are numbered from 01 to 99 (or 001 to 999).
#' Consecutive numbers to designate alleles are not required.
#' \item Populations are defined by the position of the "Pop" separator.
#' To group various populations, just remove relevant "Pop" separators.
#' \item Missing data should be indicated as 00 (or 000) rather than blanks.
#' There are three possibilities for missing data :
#' no information (0000) or (000000), partial information
#' for first allele (1000) or (010000), partial information
#' for second allele (0010) or (000010).
#' \item The number of locus names should correspond to the number of genotypes
#' in each row.
#' If you remove one or several loci from your input file,
#' you should remove both their names and the corresponding genotypes.
#' \item No empty lines should be found within the file.
#' \item No more than one empty line should be present at the end of file.
#' }
#' \strong{not an ideal genomic format: } The nice thing about RADseq dataset is
#' that you have several important genotypes and markers metadata
#' (chromosome, locus, snp, position, read depth, allele depth, etc.) available,
#' these are all lacking in the genepop format. This format is kept for archival
#' reasons in radiator.

#' @export
#' @rdname tidy_genepop
#' @examples
#' \dontrun{
#' # We will use the genepop dataset provided with adegenet package
#' require("adegenet")
#'
#' # The simplest form of the function:
#' nancycats.tidy <- radiator::tidy_genepop(
#'     data = system.file(
#'         "files/nancycats.gen",
#'         package = "adegenet"
#'         )
#'     )
#'
#' # To output a data frame in wide format, with markers in separate columns:
#' nancycats.wide <- radiator::tidy_genepop(
#'     data = system.file(
#'     "files/nancycats.gen",
#'     package="adegenet"
#' ),
#'     tidy = FALSE
#' )
#' }


#' @references Raymond M. & Rousset F, (1995).
#' GENEPOP (version 1.2): population genetics software for exact tests
#' and ecumenicism.
#' J. Heredity, 86:248-249
#' @references Rousset F.
#' genepop'007: a complete re-implementation of the genepop software
#' for Windows and Linux.
#' Molecular Ecology Resources.
#' 2008, 8: 103-106.
#' doi:10.1111/j.1471-8286.2007.01931.x

#' @seealso \href{https://genepop.curtin.edu.au}{genepop}

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


tidy_genepop <- function(data, strata = NULL, tidy = TRUE, filename = NULL) {

  # Checking for missing and/or default arguments-------------------------------
  if (missing(data)) rlang::abort("genepop file missing")

  # Import data ------------------------------------------------------------------
  if (is.vector(data)) {
    data <- readr::read_delim(
      file = data,
      delim = "?",
      skip = 1,
      trim_ws = TRUE,
      col_names = "data",
      col_types = "c")
  } else {
    V1 <- NULL
    data %<>%
      dplyr::rename(.data = ., data = V1) %>%
      dplyr::slice(-1) %>% # removes genepop header
      tibble::as_tibble()
  }

  # Replace white space with only 1 space
  data$data %<>%
    stringi::stri_replace_all_regex(
      str = .,
      pattern = "\\s+",
      replacement = " ",
      vectorize_all = FALSE
    ) %>%
    # Remove unnecessary spaces
    stringi::stri_trim_right(str = ., pattern = "\\P{Wspace}")

  # Pop indices ----------------------------------------------------------------
  pop.indices <- which(data$data %in% c("Pop", "pop", "POP"))
  npop <- length(pop.indices)

  # Markers --------------------------------------------------------------------
  # With this function, it doesn't matter if the markers are on one row or one per row.
  markers <- dplyr::slice(.data = data, 1:(pop.indices[1] - 1)) %>% purrr::flatten_chr(.x = .)
  markers <- unlist(stringi::stri_split_fixed(str = markers, pattern = ","))
  markers <- stringi::stri_replace_all_fixed(
    str = markers,
    pattern = " ",
    replacement = "",
    vectorize_all = FALSE
  )

  # Remove markers from the dataset --------------------------------------------
  data %<>% dplyr::slice(-(1:(pop.indices[1] - 1)))

  # New pop indices and pop string ---------------------------------------------
  pop.indices <- which(data$data %in% c("Pop", "pop", "POP"))
  pop.indices.lagged <- diff(c(pop.indices, (nrow(data) + 1))) - 1
  pop <- factor(rep.int(1:npop, times = pop.indices.lagged))

  # remove pop indices from dataset---------------------------------------------
  data %<>% dplyr::slice(-pop.indices)

  # Scan for genotypes split on 2 lines ----------------------------------------
  # looks for lines without a comma
  problem <- which(!(1:length(data$data) %in% grep(",", data$data)))
  if (length(problem) > 0) {
    for (i in sort(problem, decreasing = TRUE)) {
      # i <- problem
      data$data[i - 1] <- paste(data$data[i - 1], data$data[i], sep = " ")
    }
    data <- dplyr::slice(.data = data, -problem)
    pop <- pop[-problem]
  }

  # preparing the dataset ------------------------------------------------------
  # separate the individuals based on the comma (",")
  data %<>%
    tidyr::separate(
      data = .,
      col = data,
      into = c("INDIVIDUALS", "GT"),
      sep = ","
    )

  # Individuals
  individuals <- dplyr::select(.data = data, INDIVIDUALS) %>%
    dplyr::mutate(INDIVIDUALS = radiator::clean_ind_names(x = INDIVIDUALS))

  # Check for duplicate individual names
  # some genepop format don't provide individuals
  if (length(unique(individuals$INDIVIDUALS)) < nrow(individuals)) {
    message("WARNING: Individual are not named correctly, please check your data")
    message("radiator will generate unique id")
    message("Conversion file to get back to your original naming scheme:\n", stringi::stri_join(getwd(),"/radiator_genepop_id_conversion.tsv"))

    bad.id <- dplyr::select(individuals, BAD_ID = INDIVIDUALS) %>%
      dplyr::mutate(INDIVIDUALS = stringi::stri_join("radiator-individual-", seq(1, nrow(.)))) %>%
      dplyr::select(INDIVIDUALS, BAD_ID)

    individuals <- dplyr::select(bad.id, INDIVIDUALS)

    readr::write_tsv(x = bad.id, file = "radiator_genepop_id_conversion.tsv")
  }

  # isolate the genotypes
  data %<>%
    dplyr::select(GT) %>%
    dplyr::mutate(
      GT = stringi::stri_replace_all_fixed(
        str = GT,
        pattern = c("\t", ","), # remove potential comma and tab
        replacement = c(" ", ""),
        vectorize_all = FALSE
      ),
      # replace white space character: [\t\n\f\r\p{Z}]
      GT = stringi::stri_replace_all_regex(
        str = GT,
        pattern = "\\s+",
        replacement = " ",
        vectorize_all = FALSE
      ),
      # trim unnecessary whitespaces at start and end of string
      GT = stringi::stri_trim_both(
        str = GT,
        pattern = "\\P{Wspace}"
      )
    )

  # create a data frame --------------------------------------------------------
  # separate the dataset by space
  data <- tibble::as_tibble(
    do.call(rbind, stringi::stri_split_fixed(str = data$GT, pattern = " ")),
    .name_repair = "minimal"
  ) %>%
    magrittr::set_colnames(x = ., markers)
  markers <- NULL

  # Population info ------------------------------------------------------------
  # Strata
  if (!is.null(strata)) {
    #join strata and data
    individuals %<>%
      dplyr::left_join(
        radiator::read_strata(strata = strata, pop.id = TRUE, verbose = FALSE) %$%
          strata,
        by = "INDIVIDUALS"
      )

  } else {
    # add pop based on internal genepop: integer and reorder the columns
    individuals %<>% dplyr::mutate(POP_ID = pop)
  }

  # combine the individuals back to the dataset
  data <- dplyr::bind_cols(individuals, data)
  individuals <- NULL

  # Scan for genotype coding and tidy ------------------------------------------
  gt.coding <- dplyr::select(.data = data, -INDIVIDUALS, -POP_ID) %>%
    purrr::flatten_chr(.) %>%
    unique(.) %>%
    nchar(.) %>%
    unique(.)

  if (length(gt.coding) != 1) {
    rlang::abort("Mixed genotype codings are not supported:
  use 1, 2 or 3 characters/numbers for alleles")
  } else {
    if (gt.coding != 6) {
      if (gt.coding == 4) gt.sep <- 2
      if (gt.coding == 2) gt.sep <- 1


      data <- radiator::rad_long(
        x = data,
        cols = c("POP_ID", "INDIVIDUALS"),
        names_to = "MARKERS",
        values_to = "GT"
      ) %>%
        tidyr::separate(
          data = ., col = GT, into = c("A1", "A2"),
          sep = gt.sep, remove = TRUE, extra = "drop"
        ) %>%
        dplyr::mutate(
          A1 = stringi::stri_pad_left(str = A1, pad = "0", width = 3),
          A2 = stringi::stri_pad_left(str = A2, pad = "0", width = 3),
          GT = stringi::stri_join(A1, A2, sep = ""),
          A1 = NULL, A2 = NULL
        ) %>%
        dplyr::arrange(MARKERS, POP_ID, INDIVIDUALS)

      if (!tidy) {
        data %<>%
          radiator::rad_wide(
            x = .,
            formula = "POP_ID + NDIVIDUALS ~ MARKERS",
            values_from = "GT"
          )
      }
    } else {
      if (tidy) {
        data %<>%
          radiator::rad_long(
            x = .,
            cols = c("POP_ID", "INDIVIDUALS"),
            names_to = "MARKERS",
            values_to = "GT"
          )
      }
    }
  }

  # writing to a file  ---------------------------------------------------------
  if (!is.null(filename)) readr::write_tsv(x = data, file = filename, col_names = TRUE)

  return(data)
} # end tidy_genepop


# write_genepop-----------------------------------------------------------------
#' @name write_genepop

#' @title Write a genepop file

#' @description Write a genepop file from a tidy data frame or GDS file/object.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @inheritParams radiator_common_arguments

#' @param pop.levels (optional, string) A character string with your populations ordered.
#' Default: \code{pop.levels = NULL}. Described in \code{\link{read_strata}}.

#' @param genepop.header The first line of the Genepop file.
#' Default: \code{genepop.header = NULL} will use "radiator genepop with date".

#' @param markers.line (optional, logical) In the genepop and structure
#' file, you can write the markers on a single line separated by
#' commas \code{markers.line = TRUE},
#' or have markers on a separate line, i.e. in one column, for the genepop file
#' (not very useful with thousands of markers) and not printed at all for the
#' structure file.
#' Default: \code{markers.line = TRUE}.

#' @param filename (optional) The file name prefix for the genepop file
#' written to the working directory. With default: \code{filename = NULL},
#' the date and time is appended to \code{radiator_genepop_}.

#' @param ... other parameters passed to the function.

#' @return A genepop file is saved to the working directory.

#' @export
#' @rdname write_genepop
#' @references Raymond M. & Rousset F, (1995).
#' GENEPOP (version 1.2): population genetics software for exact tests
#' and ecumenicism.
#' J. Heredity, 86:248-249
#' @references Rousset F.
#' genepop'007: a complete re-implementation of the genepop software
#' for Windows and Linux.
#' Molecular Ecology Resources.
#' 2008, 8: 103-106.
#' doi:10.1111/j.1471-8286.2007.01931.x

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_genepop <- function(
    data,
    pop.levels = NULL,
    genepop.header = NULL,
    markers.line = TRUE,
    filename = NULL,
    ...
) {


  # # For testing
  # pop.levels = NULL
  # genepop.header = NULL
  # markers.line = TRUE
  # filename = NULL

  options(stringsAsFactors = FALSE)
  cli::cli_progress_step("Reading data")
  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("Input file missing")

  # File type detection----------------------------------------------------------
  data.type <- radiator::detect_genomic_format(data)

  # Import data ---------------------------------------------------------------

  if (data.type %in% c("SeqVarGDSClass", "gds.file")) {
    if (data.type == "gds.file") data %<>% radiator::read_rad(data = .)
    data <- gds2tidy(gds = data, pop.id = FALSE, parallel.core = parallel::detectCores() - 1)
    data.type <- "tbl_df"
  } else {
    if (is.vector(data)) data %<>% radiator::tidy_wide(data = ., )
  }

  if (!rlang::has_name(data, "GT")) {
    cli::cli_progress_step("Recoding genotypes...")
    data <- gt_recoding(x = data, gt = TRUE, gt.bin = FALSE, gt.vcf = FALSE, gt.vcf.nuc = FALSE)
  }

  cli::cli_progress_step("Preparing data")

  data %<>% dplyr::select(STRATA, INDIVIDUALS, MARKERS, GT)

  # pop.levels -----------------------------------------------------------------
  if (!is.null(pop.levels)) {
    data %<>%
      dplyr::mutate(
        STRATA = factor(STRATA, levels = pop.levels, ordered = TRUE),
        STRATA = droplevels(STRATA)
      ) %>%
      dplyr::arrange(STRATA, INDIVIDUALS, MARKERS)
  } else {
    data %<>%
      dplyr::mutate(STRATA = factor(STRATA)) %>%
      dplyr::arrange(STRATA, INDIVIDUALS, MARKERS)
  }

  # Create a marker vector  ------------------------------------------------
  markers <- dplyr::distinct(data, MARKERS) %>% dplyr::arrange(MARKERS) %$% MARKERS

  # Wide format ----------------------------------------------------------------
  data  %<>%
    dplyr::arrange(MARKERS) %>%
    radiator::rad_wide(x = ., formula = "STRATA + INDIVIDUALS ~ MARKERS", values_from = "GT") %>%
    dplyr::arrange(STRATA, INDIVIDUALS) %>%
    dplyr::mutate(INDIVIDUALS = paste(INDIVIDUALS, ",", sep = ""))

  # Write the file in genepop format -------------------------------------------
  # Date and time
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")

  # Filename -------------------------------------------------------------------
  if (is.null(filename)) {
    filename <- stringi::stri_join("radiator_genepop_", file.date, ".gen")
  } else {
    filename.problem <- file.exists(filename)
    if (filename.problem) {
      filename <- stringi::stri_join(filename, "_genepop_", file.date, ".gen")
    } else {
      filename <- stringi::stri_join(filename, "_genepop", ".gen")
    }
  }

  # genepop header  ------------------------------------------------------------
  if (is.null(genepop.header)) {
    genepop.header <- stringi::stri_join("radiator genepop ", file.date)
  }

  # genepop construction
  # could probably use purrr here...
  cli::cli_progress_step("Writing genepop")

  filename.connection <- file(filename, "w") # open the connection to the file
  writeLines(text = genepop.header, con = filename.connection, sep = "\n") # write the genepop header
  if (markers.line) { # write the markers on a single line
    writeLines(text = stringi::stri_join(markers, sep = ",", collapse = ", "), con = filename.connection, sep = "\n")
  } else {# write the markers on a single column (separate lines)
    writeLines(text = stringi::stri_join(markers, sep = "\n"), con = filename.connection, sep = "\n")
  }
  close(filename.connection) # close the connection

  write_gen <- function(x, filename) {
    readr::write_delim(x = as.data.frame("pop"), file = filename, delim = "\n", append = TRUE, col_names = FALSE)
    readr::write_delim(x = x, file = filename, delim = " ", append = TRUE, col_names = FALSE)
  }

  purrr::walk(
    .x = dplyr::group_split(.tbl = data, STRATA, .keep = FALSE),
    .f = write_gen,
    filename = filename
  )
  cli::cli_progress_step(stringi::stri_join("Genepop file: ", filename))

  invisible(filename)
}# End write_genepop
