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
#' @param old_format (optional, logical) For fst file. Use \code{TRUE} to read fst
#' files generated with a fst package version lower than v.0.8.0
#' Default: \code{old_format = FALSE}.

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
  old_format = FALSE) {

  ## TEST
  # columns = NULL
  # from = 1
  # to = NULL
  # as.data.table = FALSE
  # old_format = FALSE

  # detect format---------------------------------------------------------------
  data.type <- detect_genomic_format(data)

  # FST FILE -------------------------------------------------------------------
  if (data.type == "fst.file") {

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
        as.data.table = as.data.table, old_format = old_format) %>%
        tibble::as_data_frame(.)
      return(data)
    } else {
      data.old <- suppressWarnings(fst::read_fst(path = data, old_format = TRUE)) %>%
        radiator::write_rad(data = ., path = data)
      data.old <- NULL
      data <- fst::read_fst(
        path = data, columns = columns, from = from, to = to,
        as.data.table = as.data.table, old_format = old_format) %>%
        tibble::as_data_frame(.)
      message("\nThis .rad file was created with an earlier version of the fst package")
      message("A new version with the same name was written")
    }
  }#End fst.file

  # GDS file -------------------------------------------------------------------
  if (data.type == "gds.file") {
    message("Opening GDS file connection")
    data <- SeqArray::seqOpen(gds.fn = data, readonly = FALSE)

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

    message("Setting filters to:")
    message("    number of samples: ", length(w.s))
    message("    number of markers: ", length(w.m))
    SeqArray::seqSetFilter(object = data,
                           variant.id = w.m,
                           sample.id = w.s,
                           verbose = FALSE)

    # Checks--------------------------------------------------------------------
    check <- SeqArray::seqGetFilter(data)
    if (length(check$sample.sel[check$sample.sel]) != length(w.s)) {
      stop("Number of samples don't match, contact author")
    }
    if (length(check$variant.sel[check$variant.sel]) != length(w.m)) {
      stop("Number of markers don't match, contact author")
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
#'
#' @param path (optional) For tidy data frame, the path to write the data on disk.
#'

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


write_rad <- function(data, path) {

  # detect format---------------------------------------------------------------
  data.type <- class(data)

  if (unique(data.type == "SeqVarGDSClass")) {
    rad_sample <- purrr::safely(.f = function(x) gdsfmt::read.gdsn(gdsfmt::index.gdsn(node = x, path = "radiator/STRATA/INDIVIDUALS")))
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

    message("Setting filters to:")
    message("    number of samples: ", length(w.s))
    message("    number of markers: ", length(w.m))
    SeqArray::seqSetFilter(object = data,
                           variant.id = w.m,
                           sample.id = w.s,
                           verbose = FALSE)
    message("Closing connection with GDS file:\n", data$filename)
    SeqArray::seqClose(data)
  } else {
    if (missing(path)) stop("The function requires the path for the file")
    tibble::as_data_frame(data) %>%
      fst::write_fst(x = ., path = path, compress = 85)
  }



}#End write_rad
