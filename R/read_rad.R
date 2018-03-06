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

#' @return A tidy data frame in the global environment.
#' @export
#' @rdname read_rad
#' @importFrom fst read.fst
#' @importFrom purrr safely
#'
#'
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

read_rad <- function(data) {
  # since version 0.8.4 there is a distinction between old and new format...
  # Catch error while reading
  safe_fst <- purrr::safely(.f = fst::read.fst)
  data.safe <- safe_fst(data)
  if (is.null(data.safe$error)) {
    return(data.safe$result)
  } else {
    data <- suppressWarnings(fst::read.fst(path = data, old_format = TRUE))
    message("This .rad file was created with an earlier version of the fst package")
    message("Please re-write this tidy data with: radiator::write_rad")
    message("    as this format will not be supported in future releases.")
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
  fst::write.fst(x = data, path = path, compress = 85)
}#End write_rad
