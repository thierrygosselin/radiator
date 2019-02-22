# read .rad file

#' @name read_rad
#' @title Read radiator file ending .rad or .gds
#' @description Fast read of .rad or .gds file.
#' The function uses \code{\link[fst]{read_fst}} or
#' CoreArray Genomic Data Structure (\href{https://github.com/zhengxwen/gdsfmt}{GDS})
#' file system.
#'
#'
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.


#' @param data A file in the working directory ending with .rad or .gds,
#' and produced by radiator, assigner or grur.

#' @param columns (optional) For fst file. Column names to read. The default is to read all all columns.
#' Default: \code{columns = NULL}.
#' @param from (optional) For fst file. Read data starting from this row number.
#' Default: \code{from = 1}.
#' @param to (optional) For fst file. Read data up until this row number.
#' The default is to read to the last row of the stored dataset.
#' Default: \code{to = NULL}.
#' @param as.data.table (optional, logical) For fst file. If \code{TRUE},
#' the result will be returned as a \code{data.table} object.
#' Any keys set on dataset \code{x} before writing will be retained.
#' This allows for storage of sorted datasets.
#' Default: \code{as.data.table = TRUE}.
#' @param old.format (optional, logical) For fst file. Use \code{TRUE} to read fst
#' files generated with a fst package version lower than v.0.8.0
#' Default: \code{old.format = FALSE}.
#' @param allow.dup (optional, logical) To allow the opening of a GDS file with
#' read-only mode when it has been opened in the same R session.
#' Default: \code{allow.dup = FALSE}.
#' @param check (optional, logical) Verify that GDS number of samples and markers
#' match.
#' Default: \code{check = TRUE}.
#' @param verbose (optional, logical) \code{verbose = TRUE} to be chatty
#' during execution.
#' Default: \code{verbose = FALSE}.

#' @details For GDS file system, \strong{read_rad} will open the GDS connection file
#' set the filters (variants and samples) based on the info found in the file.

#' @return A radiator tidy data frame or
#' GDS object (with read/write permissions) in the global environment.
#' @export
#' @rdname read_rad
#' @importFrom fst read_fst
#' @importFrom purrr safely
#'
#' @seealso \href{https://github.com/fstpackage/fst}{fst} and \href{https://github.com/zhengxwen/gdsfmt}{GDS}
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @examples
#' \dontrun{
#' require(SeqArray)
#' shark <- radiator::read_rad(data = "data.shark.gds")
#' turtle <- radiator::read_rad(data = "data.turtle.rad")
#' }

read_rad <- function(
  data,
  columns = NULL,
  from = 1, to = NULL,
  as.data.table = FALSE,
  old.format = FALSE,
  allow.dup = FALSE,
  check = TRUE,
  verbose = FALSE) {

  ## TEST
  # columns = NULL
  # from = 1
  # to = NULL
  # as.data.table = FALSE
  # old.format = FALSE

  # detect format---------------------------------------------------------------
  data.type <- detect_genomic_format(data)

  # FST FILE -------------------------------------------------------------------
  if ("fst.file" %in% data.type) {

    # since version 0.8.4 there is a distinction between old and new format...
    # Catch error while reading
    read_rad_fst <- function(
      data) {
      fst::read_fst(
        path = data, columns = NULL, from = 1, to = 2,
        as.data.table = FALSE, old_format = FALSE) %>%
        tibble::as_data_frame(.)
    }


    safe_rad_fst <- purrr::safely(.f = read_rad_fst)
    data.safe <- safe_rad_fst(data)
    if (is.null(data.safe$error)) {
      # return(data.safe$result)
      data <- fst::read_fst(
        path = data, columns = columns, from = from, to = to,
        as.data.table = as.data.table, old_format = old.format) %>%
        tibble::as_data_frame(.)
      return(data)
    } else {
      data.old <- suppressWarnings(fst::read_fst(path = data, old_format = TRUE)) %>%
        radiator::write_rad(data = ., path = data)
      data.old <- NULL
      data <- fst::read_fst(
        path = data, columns = columns, from = from, to = to,
        as.data.table = as.data.table, old_format = old.format) %>%
        tibble::as_data_frame(.)
      message("\nThis .rad file was created with an earlier version of the fst package")
      message("A new version with the same name was written")
    }
  }#End fst.file

  # GDS file -------------------------------------------------------------------
  if ("gds.file" %in% data.type) {
    if (verbose) message("Opening GDS file connection")

    if (allow.dup) {
      data <- SeqArray::seqOpen(gds.fn = data, readonly = TRUE, allow.duplicate = TRUE)
    } else {
      data <- SeqArray::seqOpen(gds.fn = data, readonly = FALSE)
    }


    rad_sample <- purrr::safely(.f = function(x) gdsfmt::read.gdsn(gdsfmt::index.gdsn(node = x, path = "radiator/individuals/INDIVIDUALS")))
    rad_markers <- purrr::safely(.f = function(x) gdsfmt::read.gdsn(gdsfmt::index.gdsn(node = x, path = "radiator/markers.meta/VARIANT_ID")))

    # Note to myself: just use , silent = TRUE and it will return null if doesnt work..
    w.s <-  rad_sample(data)
    if (!is.null(w.s$error)) {
      w.s <- SeqArray::seqGetData(gdsfile = data, var.name = "sample.id")
    } else {
      w.s <- w.s$result
    }

    w.m <- rad_markers(data)
    if (!is.null(w.m$error)) {
      w.m <- SeqArray::seqGetData(gdsfile = data, var.name = "variant.id")
    } else {
      w.m <- w.m$result
    }

    if (verbose) message("Setting filters to:")
    if (verbose) message("    number of samples: ", length(w.s))
    if (verbose) message("    number of markers: ", length(w.m))
    SeqArray::seqSetFilter(object = data,
                           variant.id = w.m,
                           sample.id = as.character(w.s),
                           verbose = FALSE)

    # Checks--------------------------------------------------------------------
    if (check) {
      check <- SeqArray::seqGetFilter(data)
      if (length(check$sample.sel[check$sample.sel]) != length(w.s)) {
        rlang::abort("Number of samples don't match, contact author")
      }
      if (length(check$variant.sel[check$variant.sel]) != length(w.m)) {
        rlang::abort("Number of markers don't match, contact author")
      }
    }
  }#End gds.file
  return(data)
}#End read_rad


#' @name write_rad
#' @title Write tidy genomic data file or close GDS file
#' @description When dataset is a tidy data frame, the function provides a
#' fast way to write a radiator \code{.rad} file.
#' The function uses \code{\link[fst]{write_fst}} with a compression level of 85,
#' that work well with RADseq dataset.
#' When the object is a CoreArray Genomic Data Structure
#' (\href{https://github.com/zhengxwen/gdsfmt}{GDS}) file system, the function
#' set filters (variants and samples) based on the info found in the file and
#' close the connection with the GDS file.
#'
#'
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @param data An object in the global environment: tidy data frame or GDS connection file

#' @param path (optional) For tidy data frame, the path to write the data on disk.
#' Default: \code{path = NULL}.

#' @param filename (optional) Name of the file (when using tsv files).
#' Default: \code{filename = NULL}.

#' @param tsv (optinal, logical) To trigger saving using \code{readr::write_tsv}.
#' Default: \code{tsv = FALSE}.

#' @param internal (optional, logical) This is used inside radiator internal code and it stops
#' from writting the file.
#' Default: \code{internal = FALSE}.

#' @param append (optional, logical) If \code{FALSE}, will overwrite existing file.
#' If \code{TRUE}, will append to existing file.
#' In both cases, if file does not exist a new file is created.
#' Default: \code{append = FALSE}.

#' @param col.names (optional, logical) Write columns names at the top of the file?
#' Must be either TRUE or FALSE.
#' Default: \code{col.names = TRUE}.

#' @param write.message (optional, character) Print a message in the console
#' after writting file.
#' With \code{write.message = NULL}, nothing is printed in the console.
#' Default: \code{write.message = "standard"}. This will print
#' \code{message("File written: ", folder_short(filename))}.

#' @param verbose (optional, logical) \code{verbose = TRUE} to be chatty
#' during execution.
#' Default: \code{verbose = FALSE}.

#' @return A file written in the working directory or nothing if it's a GDS connection file.
#' @export
#' @rdname write_rad
#' @importFrom fst write_fst
#'
#'
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @seealso \href{https://github.com/fstpackage/fst}{fst} and \href{https://github.com/zhengxwen/gdsfmt}{GDS}
#' @examples
#' \dontrun{
#' require(SeqArray)
#' radiator::write_rad(data = tidy.data, path = "data.shark.rad")
#' radiator::write_rad(data = gds.object)
#' }


write_rad <- function(
  data,
  path = NULL,
  filename = NULL,
  tsv = FALSE,
  internal = FALSE,
  append = FALSE,
  col.names = TRUE,
  write.message = "standard",
  verbose = FALSE
) {

  if (!internal) {
    if (tsv) {
      # write the dots file
      if (!is.null(path)) {
        path.filename <- file.path(path, filename)
      } else {
        path.filename <- filename
      }
      readr::write_tsv(x = data, path = path.filename, append = append, col_names = col.names)

      if (!is.null(write.message) && verbose) {
        if (write.message == "standard") {
          message("File written: ", folder_short(filename))
        } else {
          write.message
        }
      }

    } else {
      # detect format---------------------------------------------------------------
      data.type <- class(data)

      if ("SeqVarGDSClass" %in% data.type) {
        rad_sample <- purrr::safely(.f = function(x) gdsfmt::read.gdsn(gdsfmt::index.gdsn(node = x, path = "radiator/individuals/INDIVIDUALS")))
        rad_markers <- purrr::safely(.f = function(x) gdsfmt::read.gdsn(gdsfmt::index.gdsn(node = x, path = "radiator/markers.meta/VARIANT_ID")))

        w.s <-  rad_sample(data)
        if (!is.null(w.s$error)) {
          w.s <- SeqArray::seqGetData(gdsfile = data, var.name = "sample.id")
        } else {
          w.s <- w.s$result
        }

        w.m <- rad_markers(data)
        if (!is.null(w.m$error)) {
          w.m <- SeqArray::seqGetData(gdsfile = data, var.name = "variant.id")
        } else {
          w.m <- w.m$result
        }

        if (verbose) message("Setting filters to:")
        if (verbose) message("    number of samples: ", length(w.s))
        if (verbose) message("    number of markers: ", length(w.m))
        SeqArray::seqSetFilter(object = data,
                               variant.id = w.m,
                               sample.id = as.character(w.s),
                               verbose = FALSE)
        if (verbose) message("Closing connection with GDS file:\n", data$filename)
        SeqArray::seqClose(data)
      } else {
        if (is.null(path)) rlang::abort("The function requires the path of the file")
        tibble::as_data_frame(data) %>%
          fst::write_fst(x = ., path = path, compress = 85)
      }
    }
  }
}#End write_rad
