% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/detect_microsatellites.R
\name{detect_microsatellites}
\alias{detect_microsatellites}
\title{Detect microsatellites}
\usage{
detect_microsatellites(data, gmata.dir = NULL, ...)
}
\arguments{
\item{data}{(path or object)
Object in your global environment or a file in the working directory.
The tibble must contain 2 columns named: \code{MARKERS} and \code{SEQUENCE}.
When RADseq data from DArT is used, \code{\link{filter_rad}} generates
automatically this file under the name \code{whitelist.markers.tsv}.}

\item{gmata.dir}{(path) For the function to work, the path to the directory
with GMATA software needs to be given.
If not found or \code{NULL}, the function download GMATA
from github in the working directory.
Default: \code{gmata.path = NULL}.}

\item{...}{(optional) Advance mode that allows to pass further arguments
for fine-tuning the function. Also used for legacy arguments (see details or
special section)}
}
\value{
6 files are returned in the folder: detect_microsatellites:
\enumerate{
\item ".fa.fms": the fasta file of sequences
\item ".fa.fms.sat1": the summary of sequences analysed (not important)
\item ".fa.ssr": The microsatellites found per markers (see GMATA doc)
\item ".fa.ssr.sat2": Extensive summary (see GMATA doc).
\item "blacklist.microsatellites.tsv": The list of markers with microsatellites.
\item "whitelist.microsatellites.tsv": The whitelist of markers with NO microsatellites.
}
In the global environment, the object is a list with the blacklist and
the whitelist.
}
\description{
Detect Simple Sequence Reapeats (SSR)
commonly known as microsatellites...
\strong{radiator} is not re-inventing the wheel here, it uses the
software \href{https://github.com/XuewenWangUGA/GMATA}{GMATA: Genome-wide Microsatellite Analyzing Toward Application}.
}
\note{
Thanks to Peter Grewe for the idea of including this type of filter inside
radiator.
}
\examples{
\dontrun{
# The simplest way to run the function when the raw data was DArT:
mic <- radiator::detect_microsatellites(data = "my_whitelist.tsv")

# With stacks pipeline, the populations module need to be run with --fasta-loci
# You could prepare the file this way (uncomment the function):
#
# prep_stacks_fasta <- function(fasta.file) {
#   fasta <- suppressWarnings(
#     vroom::vroom(
#      file = fasta.file,
#      delim = "\t",
#      col_names = "DATA",
#      col_types = "c",
#      comment = "#"
#    ) \%>\%
#      dplyr::mutate(MARKERS = stringi::stri_sub(str = DATA, from = 3, to = 7)) \%>\%
#      tidyr::separate(data = ., col = DATA, into = c("SEQUENCE", "LOCUS"), sep = "_")
#  )
#
#  fasta <- dplyr::bind_cols(
#   dplyr::filter(fasta, MARKERS == "Locus") \%>\%
#   dplyr::select(LOCUS),
#   dplyr::filter(fasta, MARKERS != "Locus") \%>\%
#   dplyr::select(SEQUENCE)
#  ) \%>\%
#  dplyr::mutate(LOCUS = as.numeric(LOCUS))
#   return(fasta)
#  } #prep_stacks_fasta

}
}
\author{
Thierry Gosselin \email{thierrygosselin@icloud.com}
}
