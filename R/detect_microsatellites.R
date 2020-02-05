# microsatellites ---------------------------------------------------------------
#' @name detect_microsatellites
#' @title Detect microsatellites
#' @description Detect Simple Sequence Reapeats (SSR)
#' commonly known as microsatellites...
#' \strong{radiator} is not re-inventing the wheel here, it uses the
#' software \href{https://github.com/XuewenWangUGA/GMATA}{GMATA: Genome-wide Microsatellite Analyzing Toward Application}.
#'
#' @param data (path or object)
#' Object in your global environment or a file in the working directory.
#' The tibble must contain 2 columns named: \code{MARKERS} and \code{SEQUENCE}.
#' When RADseq data from DArT is used, \code{\link{filter_rad}} generates
#' automatically this file under the name \code{whitelist.markers.tsv}.

#' @param gmata.dir (path) For the function to work, the path to the directory
#' with GMATA software needs to be given.
#' If not found or \code{NULL}, the function download GMATA
#' from github in the working directory.
#' Default: \code{gmata.path = NULL}.

#' @return 5 files are returned in the folder: detect_microsatellites:
#' \enumerate{
#' \item ".fa.fms": the fasta file of sequences
#' \item ".fa.fms.sat1": the summary of sequences analysed (not important)
#' \item ".fa.ssr": The microsatellites found per markers (see GMATA doc)
#' \item ".fa.ssr.sat2": Extensive summary (see GMATA doc).
#' \item "blacklist.microsatellites.tsv": The list of markers with microsatellites.
#' }
#' The returned object is the blacklist of markers with microsatellites.


#' @inheritParams radiator_common_arguments

#' @examples
#' \dontrun{
#' # The simplest way to run the function:
#' mic <- radiator::detect_microsatellites(data = "mywhitelist.tsv")
#' }
#' @export
#' @rdname detect_microsatellites
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


detect_microsatellites <- function(data, gmata.dir = NULL, ...) {

  # TESTS
  # gmata.dir = "/Users/thierry/Downloads/GMATA-master"
  # gmata.dir = NULL
  # data = "whitelist.markers.tsv"

  # obj.keeper <- c(ls(envir = globalenv()), "blacklist")
  verbose <- TRUE

  # Cleanup---------------------------------------------------------------------
  radiator_function_header(f.name = "detect_microsatellites", verbose = TRUE)
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (verbose) message("Execution date@time: ", file.date)
  old.dir <- getwd()
  opt.change <- getOption("width")
  options(width = 70)
  timing <- radiator_tic()
  #back to the original directory and options
  on.exit(setwd(old.dir), add = TRUE)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(radiator_toc(timing), add = TRUE)
  on.exit(radiator_function_header(f.name = "detect_microsatellites", start = FALSE, verbose = verbose), add = TRUE)
  # on.exit(rm(list = setdiff(ls(envir = sys.frame(-1L)), obj.keeper), envir = sys.frame(-1L)))

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("data argument value is missing")

  # Function call and dotslist -------------------------------------------------
  rad.dots <- radiator_dots(
    func.name = as.list(sys.call())[[1]],
    fd = rlang::fn_fmls_names(),
    args.list = as.list(environment()),
    dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE),
    keepers = NULL,
    verbose = TRUE
  )

  # Folders---------------------------------------------------------------------
  path.folder <- generate_folder(
    f = NULL,
    rad.folder = "detect_microsatellites",
    internal = FALSE,
    file.date = file.date,
    verbose = TRUE)

  # write the dots file
  write_rad(
    data = rad.dots,
    path = path.folder,
    filename = stringi::stri_join(
      "radiator_detect_microsatellites_args_", file.date, ".tsv"),
    tsv = TRUE,
    internal = FALSE,
    verbose = TRUE
  )


  # GMATA CHECKS ---------------------------------------------------------------
  download.gmata <- FALSE
  if (is.null(gmata.dir)) {
    download.gmata <- TRUE
  } else {
    check.gmata <- file.path(gmata.dir, "gmat.pl")
    if (file.exists(check.gmata)) {
      message("GMATA found")
      download.gmata <- FALSE
      gmat.pl.path <- file.path(gmata.dir, "gmat.pl")
    } else {
      message("GMATA not found, downloading the software")
      download.gmata <- TRUE
    }
  }
  check.gmata <- NULL

  if (download.gmata) {
    gmata.temp.dir <- file.path(path.folder, "radiator_dep")
    dir.create(path = gmata.temp.dir)
    gmata <- file.path(gmata.temp.dir, "gmata.zip")
    file.create(gmata)
    utils::download.file(
      url = "https://github.com/XuewenWangUGA/GMATA/archive/master.zip",
      destfile = gmata
      )
    utils::unzip(zipfile = gmata, exdir = gmata.temp.dir)
    gmata.dir <- file.path(gmata.temp.dir, "GMATA-master")
    gmat.pl.path <- file.path(gmata.dir, "gmat.pl")
  }

  Sys.chmod(gmat.pl.path) # make it executable

  # import file
  if (is.vector(data)) {
    data <- readr::read_tsv(
      file = data,
      col_types = readr::cols(.default = readr::col_character())
    )
  }

  # check for the required columns
  want <- c("MARKERS", "SEQUENCE")
  if (!rlang::has_name(data, "MARKERS") || !rlang::has_name(data, "SEQUENCE")) {
    rlang::abort("data tibble missing required column(s)")
  }

  # Generate seq file
  filename <- file.path(
    path.folder,
    stringi::stri_join("radiator_detect_microsatellites_", file.date, ".fa")
    )

  data %>%
    dplyr::select(dplyr::one_of(want)) %>%
    dplyr::mutate(
      MARKERS = stringi::stri_join(">", MARKERS),
      ID = seq(from = 1, to = n(), by = 1)
    ) %>%
    data.table::as.data.table(.) %>%
    data.table::melt.data.table(
      data = .,
      id.vars = "ID",
      variable.name = "MARKERS_SEQUENCE",
      value.name = "MICROSATELLITES",
      variable.factor = FALSE) %>%
    tibble::as_tibble(.) %>%
    dplyr::arrange(ID) %>%
    dplyr::select(MICROSATELLITES) %>%
    readr::write_tsv(x = ., path = filename, col_names = FALSE)


  # Run GMATA ------------------------------------------------------------------
  setwd(gmata.dir)
  gmata.command <- paste("perl", "gmat.pl", "-i", filename)
  system(command = gmata.command)

  setwd(path.folder)
  # removing duplicate file
  file.remove(filename)

  # importing analysis in R ----------------------------------------------------
  message("Analysing the results...")
  ssr.file <- list.files(path = path.folder, pattern = ".fa.ssr", full.names = TRUE)[1]
  micro <- readr::read_tsv(
    file = ssr.file,
    col_names = c("MARKERS", "SEQUENCE_LENGTH", "START", "END", "REPETITIONS", "MOTIF"),
    skip = 1,
    col_types = "ciiiic"
    ) %>%
    dplyr::mutate(
      MARKERS = stringi::stri_replace_all_fixed(
        str = MARKERS,
        pattern = ">",
        replacement = "",
        vectorize_all = FALSE
        )
    ) %>%
    readr::write_tsv(x = ., path = ssr.file)

  n.micro <- nrow(micro)
  message("Number of microsatellite(s) detected: ", n.micro)
  if (n.micro > 0L) {
    blacklist <- micro %>%
      dplyr::distinct(MARKERS) %>%
      readr::write_tsv(x = ., path = "blacklist.microsatellites.tsv")
    message("Number of unique markers with microsatellite(s): ", nrow(blacklist))
  } else {
    blacklist <- NULL
  }
  mean.rep <- round(mean(x = micro$REPETITIONS, na.rm = TRUE), 0)
  min.rep <- round(min(micro$REPETITIONS, na.rm = TRUE), 0)
  max.rep <- round(max(micro$REPETITIONS, na.rm = TRUE), 0)
  message("Mean number of repetitions [min-max]: ", mean.rep, " [", min.rep, " - ", max.rep, "]")

  return(blacklist)
}#End detect_microsatellites

