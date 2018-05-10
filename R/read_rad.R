# read .rad file

#' @name read_rad
#' @title Read tidy genomic data file ending .rad
#' @description Fast read .rad file.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.
#' The function uses \code{\link[fst]{read.fst}}.

#' @param data A file in the working directory ending with .rad,
#' a tidy genomic data produced by radiator, assigner or grur.

#' @inheritParams fst::read.fst


#' @return A tidy data frame in the global environment.
#' @export
#' @rdname read_rad
#' @importFrom fst read.fst
#' @importFrom purrr safely
#'
#'
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

read_rad <- function(
  data,
  columns = NULL,
  from = 1, to = NULL,
  as.data.table = FALSE, old_format = FALSE) {
  # since version 0.8.4 there is a distinction between old and new format...
  # Catch error while reading
  read_rad_fst <- function(
    data) {
    fst::read.fst(
      path = data, columns = NULL, from = 1, to = 2,
      as.data.table = FALSE, old_format = FALSE) %>%
      tibble::as_data_frame(.)
  }


  safe_rad_fst <- purrr::safely(.f = read_rad_fst)
  data.safe <- safe_rad_fst(data)
  if (is.null(data.safe$error)) {
    # return(data.safe$result)
    data <- fst::read.fst(
      path = data, columns = columns, from = from, to = to,
      as.data.table = as.data.table, old_format = old_format) %>%
      tibble::as_data_frame(.)
    return(data)
  } else {
    data.old <- suppressWarnings(fst::read.fst(path = data, old_format = TRUE)) %>%
      radiator::write_rad(data = ., path = data)
    data.old <- NULL
    data <- fst::read.fst(
      path = data, columns = columns, from = from, to = to,
      as.data.table = as.data.table, old_format = old_format) %>%
      tibble::as_data_frame(.)
    message("\nThis .rad file was created with an earlier version of the fst package")
    message("A new version with the same name was written")

    # message("    as this format will not be supported in future releases.\n")
    return(data)
  }
}#End read_rad


#' @name write_rad
#' @title Write tidy genomic data file ending .rad
#' @description Fast write .rad file.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.
#' The function uses \code{\link[fst]{write.fst}} with a compression level of 85,
#' that work well with RADseq dataset.

#' @param data A file in the working directory ending with .rad,
#' a tidy genomic data produced by radiator, assigner or grur.
#'
#' @param path Path to write the data on disk

#' @return A tidy data frame in the global environment.
#' @export
#' @rdname write_rad
#' @importFrom fst write.fst
#'
#'
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

write_rad <- function(data, path) {
  tibble::as_data_frame(data) %>%
  fst::write.fst(x = ., path = path, compress = 85)
}#End write_rad
