# Separate a column (markers) into CHROM LOCUS and POS

#' @name separate_markers

#' @title Separate markers column into chrom, locus and pos

#' @description Radiator uses unique marker names by combining
#' \code{CHROM}, \code{LOCUS}, \code{POS} columns, with double underscore
#' separators, into \code{MARKERS = CHROM__LOCUS__POS}. All the info is kept most
#' of the time, but in case user need to get back to the original metadata the
#' function provides an easy wau to do it.

#' @param data An object with a column named \code{MARKERS}.
#' If \code{CHROM}, \code{LOCUS}, \code{POS} are already present,
#' the columns are removed before separating the \code{MARKERS} column.
#' The data can be whitelists and blacklists of markers or datasets.

#' @param sep (optional, character) Separator used to identify the different
#' field in the \code{MARKERS} column.
#' Default: \code{sep = "__"}.

#' @return The same data in the global environment, with 3 new columns:
#' \code{CHROM}, \code{LOCUS}, \code{POS}
#' @rdname separate_markers

#' @examples
#' \dontrun{
#' whitelist <- radiator::separate_markers(data = whitelist.markers)
#' tidy.data <- radiator::separate_markers(data = bluefintuna.data)
#' }#' @export


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

separate_markers <- function(data, sep = "__") {
# data.bk <- data
# sep <- "__"

# check if markers column is present
if (!tibble::has_name(data, "MARKERS")) {
  stop("The data require a column named MARKERS")
}
want <- c("CHROM", "LOCUS", "POS")
data <- suppressWarnings(
  data %>%
    dplyr::select(-dplyr::one_of(want)) %>%
    tidyr::separate(data = ., col = "MARKERS",
                    into = want, sep = sep, remove = FALSE))
return(data)
}#End separate_markers
