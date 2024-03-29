% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gds.R
\name{reset_filters}
\alias{reset_filters}
\title{Reset filters (individuals and markers) in radiator GDS object.}
\usage{
reset_filters(
  gds,
  list.filters = TRUE,
  reset.all = FALSE,
  filter.individuals.missing = FALSE,
  filter.individuals.heterozygosity = FALSE,
  filter.individuals.coverage.total = FALSE,
  detect.mixed.genomes = FALSE,
  detect.duplicate.genomes = FALSE,
  filter.reproducibility = FALSE,
  filter.monomorphic = FALSE,
  filter.common.markers = FALSE,
  filter.ma = FALSE,
  filter.mean.coverage = FALSE,
  filter.genotyping = FALSE,
  filter.snp.position.read = FALSE,
  filter.snp.number = FALSE,
  filter.short.ld = FALSE,
  filter.long.ld = FALSE,
  filter.hwe = FALSE,
  filter.whitelist = FALSE
)
}
\arguments{
\item{gds}{(2 options) A Genomic Data Structure (GDS) file or object

\emph{How to get GDS ?}
Look into \code{\link{tidy_genomic_data}},
\code{\link{read_vcf}},
\code{\link{tidy_vcf}} or \code{\link{write_gds}}.}

\item{list.filters}{filters (logical, optional) List current active filters for individuals and markers.
Default: \code{list.filters = TRUE}.}

\item{reset.all}{(logical, optional) Reset all individuals and markers filters.
Default: \code{reset.all = FALSE}.}

\item{filter.individuals.missing}{(logical, optional)
Default: \code{filter.individuals.missing = FALSE}.}

\item{filter.individuals.heterozygosity}{(logical, optional)
Default: \code{filter.individuals.heterozygosity = FALSE}.}

\item{filter.individuals.coverage.total}{(logical, optional)
Default: \code{filter.individuals.coverage.total = FALSE}.}

\item{detect.mixed.genomes}{(logical, optional)
Default: \code{detect.mixed.genomes = FALSE}.}

\item{detect.duplicate.genomes}{(logical, optional)
Default: \code{detect.duplicate.genomes = FALSE}.}

\item{filter.reproducibility}{(logical, optional)
Default: \code{filter.reproducibility = FALSE}.}

\item{filter.monomorphic}{(logical, optional)
Default: \code{filter.monomorphic = FALSE}.}

\item{filter.common.markers}{(logical, optional)
Default: \code{filter.common.markers = FALSE}.}

\item{filter.ma}{(logical, optional)
Default: \code{filter.ma = FALSE}.}

\item{filter.mean.coverage}{(logical, optional)
Default: \code{filter.mean.coverage = FALSE}.}

\item{filter.genotyping}{(logical, optional)
Default: \code{filter.genotyping = FALSE}.}

\item{filter.snp.position.read}{(logical, optional)
Default: \code{filter.snp.position.read = FALSE}.}

\item{filter.snp.number}{(logical, optional)
Default: \code{filter.snp.number = FALSE}.}

\item{filter.short.ld}{(logical, optional)
Default: \code{filter.short.ld = FALSE}.}

\item{filter.long.ld}{(logical, optional)
Default: \code{filter.long.ld = FALSE}.}

\item{filter.hwe}{(logical, optional)
Default: \code{filter.hwe = FALSE}.}

\item{filter.whitelist}{(logical, optional)
Default: \code{filter.whitelist = FALSE}.}
}
\description{
List current active filters in radiator GDS object.
Reset specific filters or all filters at once.
}
\examples{
\dontrun{
# List active filters for individuals and markers
reset_filters(gds)

# You changed your decision concerning the genotyping threshold or
# entered the wrong one:
reset_filters(gds, filter.genotyping = TRUE)
}
}
\seealso{
\code{\link{sync_gds}}, \code{\link{list_filters}}.
}
\author{
Thierry Gosselin \email{thierrygosselin@icloud.com}
}
\keyword{internal}
