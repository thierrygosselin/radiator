# write a betadiv from a tidy data frame

#' @name write_betadiv

#' @title Write a betadiv file from a tidy data frame

#' @description Write a betadiv file from a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.

#' @return A betadiv object is returned.

#' @export
#' @rdname write_betadiv

#' @references Lamy T, Legendre P, Chancerelle Y, Siu G, Claudet J (2015)
#' Understanding the Spatio-Temporal Response of Coral Reef Fish Communities to
#' Natural Disturbances: Insights from Beta-Diversity Decomposition.
#' PLoS ONE, 10, e0138696.

#' @seealso \code{beta.div} is available on Pierre Legendre web site \url{http://adn.biol.umontreal.ca/~numericalecology/Rcode/} \code{randomForestSRC} is available on CRAN \url{http://cran.r-project.org/web/packages/randomForestSRC/} and github \url{https://github.com/ehrlinger/randomForestSRC}

#' @author Laura Benestan \email{laura.benestan@@icloud.com} and
#' Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_betadiv <- function(data) {

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    data <- radiator::tidy_wide(data = data, import.metadata = TRUE)
  }

  # necessary steps to make sure we work with unique markers and not duplicated LOCUS
  if (tibble::has_name(data, "LOCUS") && !tibble::has_name(data, "MARKERS")) {
    data <- dplyr::rename(.data = data, MARKERS = LOCUS)
  }

  # Compute count and Minor Allele Frequency -----------------------------------
  # We split the alleles here to prep for MAF
  # need to compute REF/ALT allele for non VCF file
  if (!tibble::has_name(data, "GT_VCF")) {
    ref.change <- radiator::calibrate_alleles(data = data)$input
    data <- dplyr::left_join(data, ref.change, by = c("MARKERS", "INDIVIDUALS"))
  }

  # MAF
  betadiv <- dplyr::select(.data = data, MARKERS, POP_ID, GT_VCF) %>%
    dplyr::filter(GT_VCF != "./.") %>%
    dplyr::group_by(MARKERS, POP_ID) %>%
    dplyr::summarise(
      N = as.numeric(n()),
      PQ = as.numeric(length(GT_VCF[GT_VCF == "1/0" | GT_VCF == "0/1"])),
      QQ = as.numeric(length(GT_VCF[GT_VCF == "1/1"]))
    ) %>%
    dplyr::mutate(MAF = ((QQ * 2) + PQ) / (2 * N)) %>%
    dplyr::ungroup(.) %>%
    dplyr::select(POP_ID, MARKERS, MAF) %>%
    dplyr::group_by(POP_ID) %>%
    tidyr::spread(data = ., key = MARKERS, value = MAF) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(POP_ID = as.integer(POP_ID))

  return(betadiv)
} # end write_betadiv

