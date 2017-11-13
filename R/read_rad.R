# read .rad file

#' @name read_rad

#' @title Read tidy genomic data file ending .rad

#' @description Fast read .rad file.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @param data A file in the working directory ending with .rad,
#' a tidy genomic data produced by radiator, assigner or grur.

#' @return A tidy data frame in the global environment.
#' @export
#' @rdname read_rad
#' @importFrom fst read.fst
#'
#'
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

read_rad <- function(data) {
  data <- fst::read.fst(path = data)
}#End read_rad
