# write an Arlequin file from a tidy data frame

#' @name write_arlequin
#' @title Write an arlequin file from a tidy data frame
#' @description Write a arlequin file from a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.

#' @param pop.levels (optional, string) A character string with your populations ordered.
#' Default: \code{pop.levels = NULL}.

#' @param filename (optional) The file name prefix for the arlequin file
#' written to the working directory. With default: \code{filename = NULL},
#' the filename generated follow this \code{radiator_arlequin_DATE@TIME.csv}.

#' @param ... other parameters passed to the function.

#' @return An arlequin file is saved to the working directory.

#' @references Excoffier, L.G. Laval, and S. Schneider (2005)
#' Arlequin ver. 3.0: An integrated software package for population genetics
#' data analysis. Evolutionary Bioinformatics Online 1:47-50.

#' @export
#' @rdname write_arlequin
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_arlequin <- function(
  data,
  pop.levels = NULL,
  filename = NULL,
  ...
) {

  # Checking for missing and/or default arguments ******************************
  cli::cli_progress_step("Reading data")
  if (missing(data)) rlang::abort("Input file missing")

  # Import data ---------------------------------------------------------------
  data %<>% radiator::tidy_wide(data = .)

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
  markers <- dplyr::distinct(.data = data, MARKERS) %>%
    dplyr::arrange(MARKERS) %>%
    purrr::flatten_chr(.)
  npop <- length(unique(data$STRATA))

  # arlequin format ----------------------------------------------------------------
  data %<>%
    dplyr::mutate(
      A1 = stringi::stri_sub(str = GT, from = 1, to = 3),
      A2 = stringi::stri_sub(str = GT, from = 4, to = 6),
      GT = NULL
    ) %>%
    radiator::rad_long(x = ., cols = c("STRATA", "INDIVIDUALS", "MARKERS"), names_to = "ALLELES", values_to = "GT") %>%
    dplyr::mutate(
      GT = stringi::stri_replace_all_fixed(str = GT, pattern = "000", replacement = "-9", vectorize_all = FALSE)
    ) %>%
    radiator::rad_wide(x = ., formula = "INDIVIDUALS + STRATA  + ALLELES ~ MARKERS", values_from = "GT") %>%
    dplyr::arrange(STRATA, INDIVIDUALS) %>%
    dplyr::select(-ALLELES)

  # Write the file in arlequin format -----------------------------------------

  # date & time
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")

  # Filename
  if (is.null(filename)) {
    filename.temp <- generate_filename(extension = "arlequin")
    filename.short <- filename.temp$filename.short
    filename <- filename.temp$filename
  } else {
    filename <- stringi::stri_join(filename, "_arlequin.csv")
    filename.problem <- file.exists(filename)
    if (filename.problem) {
      filename.temp <- generate_filename(extension = "arlequin")
      filename.short <- filename.temp$filename.short
      filename <- filename.temp$filename
    }
    filename.short <- filename
  }

  cli::cli_progress_step("Writing arlequin")
  filename.connection <- file(filename, "w") # open the connection to the file
  # Profile section
  write("[Profile]", file = filename.connection)
  write(paste('Title = "radiator_arlequin_export', " ", file.date, '"'), file = filename.connection, append = TRUE)
  write(paste("NbSamples = ", npop), file = filename.connection, append = TRUE) # number of pop
  write(paste("GenotypicData = 1"), file = filename.connection, append = TRUE)
  write(paste("LocusSeparator = WHITESPACE"), file = filename.connection, append = TRUE)
  write(paste("GameticPhase = 0"), file = filename.connection, append = TRUE) # 0 = unknown gametic phase
  write(paste("MissingData = '?'"), file = filename.connection, append = TRUE)
  write(paste("DataType = STANDARD"), file = filename.connection, append = TRUE)
  write(paste("[Data]"), file = filename.connection, append = TRUE)
  write(paste("[[Samples]]"), file = filename.connection, append = TRUE)

  # pop <- as.character(data$STRATA) # Create a population vector
  # data.split <- split(data, pop) # split data by populations
  # for (i in 1:length(data.split)) {
  #   # i <- 1
  #   # i <- 2
  #   pop.data <- data.split[[i]]
  #   pop.name <- unique(as.character(pop.data$STRATA))
  #   n.ind <- dplyr::n_distinct(pop.data$INDIVIDUALS)
  #   write(paste("SampleName = ", pop.name), file = filename.connection, append = TRUE)
  #   write(paste("SampleSize = ", n.ind), file = filename.connection, append = TRUE)
  #   write(paste("SampleData = {"), file = filename.connection, append = TRUE)
  #   pop.data$INDIVIDUALS[seq(from = 2, to = n.ind * 2, by = 2)] <- ""
  #   pop.data <- dplyr::select(.data = pop.data, -STRATA)
  #   ncol.data <- ncol(pop.data)
  #   pop.data <- as.matrix(pop.data)
  #   write(x = t(pop.data), ncolumns = ncol.data, sep = "\t", append = TRUE, file = filename.connection)
  #   write("}", file = filename.connection, append = TRUE)
  # }

  write_arl <- function(x, filename.connection) {
    n.ind <- dplyr::n_distinct(x$INDIVIDUALS)
    write(paste("SampleName = ", unique(as.character(x$STRATA))), file = filename.connection, append = TRUE)
    write(paste("SampleSize = ", n.ind), file = filename.connection, append = TRUE)
    write(paste("SampleData = {"), file = filename.connection, append = TRUE)
    x$INDIVIDUALS[seq(from = 2, to = n.ind * 2, by = 2)] <- ""
    x %<>% dplyr::select(-STRATA)
    ncol.data <- ncol(x)
    x <- as.matrix(x)
    write(x = t(x), ncolumns = ncol.data, sep = "\t", append = TRUE, file = filename.connection)
    write("}", file = filename.connection, append = TRUE)
  }

  data.split <- dplyr::group_split(.tbl = data, STRATA, .keep = TRUE)
  purrr::walk(
    .x = data.split,
    .f = write_arl,
    filename.connection = filename.connection
  )


  write(paste("[[Structure]]"), file = filename.connection, append = TRUE)
  write(paste("StructureName = ", "\"", "One cluster", "\""), file = filename.connection, append = TRUE)
  write(paste("NbGroups = 1"), file = filename.connection, append = TRUE)
  write(paste("Group = {"), file = filename.connection, append = TRUE)

  write_arl_end <- function(x, filename.connection) {
    write(paste("\"", unique(x$STRATA), "\"", sep = ""), file = filename.connection, append = TRUE)
  }
  purrr::walk(
    .x = data.split,
    .f = write_arl_end,
    filename.connection = filename.connection
  )
  # for (i in 1:length(data.split)) {
  #   pop.data <- data.split[[i]]
  #   pop.name <- unique(pop.data$STRATA)
  #   write(paste("\"", pop.name, "\"", sep = ""), file = filename.connection, append = TRUE)
  # }
  write(paste("}"), file = filename.connection, append = TRUE)
  cli::cli_progress_step(stringi::stri_join("Arlequin file: ", filename.short))
  invisible(data)
} # end write_arlequin
