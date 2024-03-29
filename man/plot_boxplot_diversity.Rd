% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{plot_boxplot_diversity}
\alias{plot_boxplot_diversity}
\alias{radiator_helper_plot}
\title{Box plot of the diversity (Gene and Haplotypes)}
\usage{
plot_boxplot_diversity(data, aes.x.y, y.title)

radiator_helper_plot(
  data,
  strata = FALSE,
  stats,
  x.axis.title = NULL,
  x.breaks = NULL,
  text.orientation = "horizontal",
  plot.filename = NULL
)
}
\arguments{
\item{data}{containing the stats}

\item{aes.x.y}{The GGPLOT2 aesthetics,
e.g. aes.x.y = aes(x = factor(POP_ID), y = GENE_DIVERSITY, na.rm = T).}

\item{y.title}{Title of the y-axis.}

\item{strata}{Containing strata info for the visual.
Default: \code{strata = FALSE}.}

\item{stats}{e.g. MAC, MAF or MAD.}

\item{x.axis.title}{Title for the x axis}

\item{text.orientation}{Orientation of the x axis text.
Default: \code{text.orientation = "horizontal"horizontal}.}

\item{plot.filename}{To save the figure provide a filename.
Default: \code{plot.filename = NULL}}
}
\description{
GGPLOT2 Box plot of the diversity (Gene and Haplotypes).

GGPLOT2 radiator helper plot used in MAC, MAF and MAD filter.
Lines and points, to help decide filter threshold.
}
\keyword{internal}
