#' @name tidy_genepop

#' @title genepop to tidy dataframe

#' @description Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.
#' The function \code{tidy_genepop} reads a file in the
#' \href{http://genepop.curtin.edu.au/help_input.html}{genepop format}
#' (see details and note for convention) and output a data frame in wide or long/tidy format.
#'
#' To manipulate and prune the dataset prior to tidying, use the functions
#' \code{\link[radiator]{tidy_genomic_data}} and
#' \code{\link[radiator]{genomic_converter}}, that uses blacklist and whitelist along
#' several other filtering options.

#' @param data A \href{http://genepop.curtin.edu.au/help_input.html}{genepop}
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

#' @details \href{http://genepop.curtin.edu.au/help_input.html}{genepop format}
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
#' \href{http://genepop.curtin.edu.au/help_input.html}{for genepop format examples}.
#' \item \strong{population identifier:} The population block are separated by
#' the word: \code{POP}, or \code{Pop} or \code{pop}. Flavors of the genepop
#' software uses the first or the last identifier of every sub-population,
#' in all output files, to name populations
#' \href{http://genepop.curtin.edu.au/help_input.html}{(more info)}.
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
#' \href{http://genepop.curtin.edu.au/help_input.html}{genepop format notes:}
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

#' @export
#' @rdname tidy_genepop


#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join slice bind_cols bind_rows
#' @importFrom purrr invoke flatten_chr
#' @importFrom stringi stri_split_fixed stri_replace_all_fixed stri_join stri_trim_both stri_trim_right stri_pad_left
#' @importFrom readr write_tsv read_tsv read_delim
#' @importFrom utils count.fields
#' @importFrom tidyr separate gather unite spread
#' @importFrom tibble as_data_frame

#' @examples
#' \dontrun{
#' We will use the genepop dataset provided with adegenet package
#' if (!require("adegenet")) install.packages("adegenet")
#'
#' The simplest form of the function:
#' nancycats.tidy <- radiator::tidy_genepop(
#' data = system.file(
#' "files/nancycats.gen",
#' package = "adegenet"))
#'
#' # To output a data frame in wide format, with markers in separate columns:
#' nancycats.wide <- radiator::tidy_genepop(
#' data = system.file(
#' "files/nancycats.gen",
#' package="adegenet"
#' ), tidy = FALSE)
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

#' @seealso \href{http://genepop.curtin.edu.au}{genepop}

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


tidy_genepop <- function(data, strata = NULL, tidy = TRUE, filename = NULL) {

  # Checking for missing and/or default arguments-------------------------------
  if (missing(data)) stop("genepop file missing")

  # Import data ------------------------------------------------------------------
  if (is.vector(data)) {
    data <- readr::read_delim(file = data, delim = "?", skip = 1, trim_ws = TRUE, col_names = "data", col_types = "c")
  } else {
    data <- data %>%
      dplyr::rename(data = V1) %>%
      dplyr::slice(-1) %>% # removes genepop header
      tibble::as_data_frame()
  }

  # Replace white space with only 1 space
  data$data <- stringi::stri_replace_all_regex(
    str = data$data,
    pattern = "\\s+",
    replacement = " ",
    vectorize_all = FALSE
  )

  # Remove unnecessary spaces
  data$data <- stringi::stri_trim_right(str = data$data, pattern = "\\P{Wspace}")

  # Pop indices ----------------------------------------------------------------
  pop.indices <- which(data$data %in% c("Pop", "pop", "POP"))
  npop <- length(pop.indices)

  # Markers --------------------------------------------------------------------
  # With this function, it doesn't matter if the markers are on one row or one per row.
  markers <- dplyr::slice(.data = data, 1:(pop.indices[1] - 1)) %>% purrr::flatten_chr(.x = .)
  markers <- unlist(stringi::stri_split_fixed(str = markers, pattern = ","))

  # Remove markers from the dataset --------------------------------------------
  data <- dplyr::slice(.data = data, -(1:(pop.indices[1] - 1)))

  # New pop indices and pop string ---------------------------------------------
  pop.indices <- which(data$data %in% c("Pop", "pop", "POP"))
  pop.indices.lagged <- diff(c(pop.indices, (nrow(data) + 1))) - 1
  pop <- factor(rep.int(1:npop, times = pop.indices.lagged))

  # remove pop indices from dataset---------------------------------------------
  data <- dplyr::slice(.data = data, -pop.indices)

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
  data <- tidyr::separate(
    data = data,
    col = data,
    into = c("INDIVIDUALS", "GT"),
    sep = ","
  )

  # Individuals
  # remove "_", ":", "," and white space character.
  # White space is defined as [\t\n\f\r\p{Z}].
  individuals <- data %>%
    dplyr::select(INDIVIDUALS) %>%
    dplyr::mutate(
      INDIVIDUALS = stringi::stri_replace_all_fixed(
        str = INDIVIDUALS,
        pattern = c("_", ",", ":"),
        replacement = c("-", "", "-"),
        vectorize_all = FALSE
      ),
      INDIVIDUALS = stringi::stri_replace_all_regex(
        str = INDIVIDUALS,
        pattern = "\\s+",
        replacement = "",
        vectorize_all = FALSE
      )
    )

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

    readr::write_tsv(x = bad.id, path = "radiator_genepop_id_conversion.tsv")
  }



  # isolate the genotypes
  data <- dplyr::select(.data = data, GT) %>%
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
  data <- tibble::as_data_frame(
    purrr::invoke(#similar to do.call
      rbind, stringi::stri_split_fixed(str = data$GT, pattern = " ")
    )
  )

  # add column names
  colnames(data) <- markers

  # Population info ------------------------------------------------------------
  # Strata
  if (!is.null(strata)) {
    if (is.vector(strata)) {
      number.columns.strata <- max(utils::count.fields(strata, sep = "\t"))
      col.types <- stringi::stri_join(rep("c", number.columns.strata), collapse = "")
      strata.df <- readr::read_tsv(file = strata, col_names = TRUE, col_types = col.types) %>%
        dplyr::rename(POP_ID = STRATA)
    } else {
      # message("strata object: yes")
      colnames(strata) <- stringi::stri_replace_all_fixed(
        str = colnames(strata),
        pattern = "STRATA",
        replacement = "POP_ID",
        vectorize_all = FALSE
      )
      strata.df <- strata
    }

    # Remove potential whitespace in pop_id
    strata.df$POP_ID <- stringi::stri_replace_all_fixed(
      strata.df$POP_ID,
      pattern = " ",
      replacement = "_",
      vectorize_all = FALSE
    )

    # Remove unwanted character in individual names
    strata.df <- strata.df %>%
      dplyr::mutate(
        INDIVIDUALS =  stringi::stri_replace_all_fixed(
          str = INDIVIDUALS,
          pattern = c("_", ",", ":"),
          replacement = c("-", "", "-"),
          vectorize_all = FALSE
        ),
        INDIVIDUALS = stringi::stri_replace_all_regex(
          str = INDIVIDUALS,
          pattern = "\\s+",
          replacement = "",
          vectorize_all = FALSE
        )
      )

    #join strata and data
    individuals <- dplyr::left_join(x = individuals, y = strata.df, by = "INDIVIDUALS")

  } else {
    # add pop based on internal genepop: integer and reorder the columns
    individuals <- individuals %>%
      dplyr::mutate(POP_ID = pop)
  }

  # combine the individuals back to the dataset
  data <- dplyr::bind_cols(individuals, data)

  # Scan for genotype coding and tidy ------------------------------------------
  gt.coding <- dplyr::select(data, -INDIVIDUALS, -POP_ID) %>%
    purrr::flatten_chr(.) %>%
    unique(.) %>%
    nchar(.) %>%
    unique(.)

  if (length(gt.coding) != 1) {
    stop("Mixed genotype codings are not supported:
  use 1, 2 or 3 characters/numbers for alleles")
  }

  if (gt.coding == 6) {
    if (tidy) {
      data <- tidyr::gather(
        data = data, key = MARKERS, value = GT, -c(POP_ID, INDIVIDUALS))
    }
  }

  if (gt.coding == 4) {
    data <- tidyr::gather(
      data = data, key = MARKERS, value = GT, -c(POP_ID, INDIVIDUALS)) %>%
      tidyr::separate(
        data = ., col = GT, into = c("A1", "A2"),
        sep = 2, remove = TRUE, extra = "drop"
      ) %>%
      dplyr::mutate(
        A1 = stringi::stri_pad_left(str = A1, pad = "0", width = 3),
        A2 = stringi::stri_pad_left(str = A2, pad = "0", width = 3)
      ) %>%
      tidyr::unite(data = ., col = GT, A1, A2, sep = "") %>%
      dplyr::arrange(MARKERS, POP_ID, INDIVIDUALS)

    if (!tidy) {
      data <- data %>%
        dplyr::group_by(POP_ID, INDIVIDUALS) %>%
        tidyr::spread(data = ., key = MARKERS, value = GT)
    }
  }

  if (gt.coding == 2) {
    data <- tidyr::gather(data = data, key = MARKERS, value = GT, -c(POP_ID, INDIVIDUALS)) %>%
      tidyr::separate(
        data = ., col = GT, into = c("A1", "A2"),
        sep = 1, remove = TRUE, extra = "drop"
      ) %>%
      dplyr::mutate(
        A1 = stringi::stri_pad_left(str = A1, pad = "0", width = 3),
        A2 = stringi::stri_pad_left(str = A2, pad = "0", width = 3)
      ) %>%
      tidyr::unite(data = ., col = GT, A1, A2, sep = "") %>%
      dplyr::arrange(MARKERS, POP_ID, INDIVIDUALS)

    if (!tidy) {
      data <- data %>%
        dplyr::group_by(POP_ID, INDIVIDUALS) %>%
        tidyr::spread(data = ., key = MARKERS, value = GT)
    }
  }

  # writing to a file  ---------------------------------------------------------
  if (!is.null(filename)) {
    readr::write_tsv(x = data, path = filename, col_names = TRUE)
  }

  return(data)
} # end tidy_genepop
