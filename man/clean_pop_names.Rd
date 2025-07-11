% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean.R
\name{clean_pop_names}
\alias{clean_pop_names}
\title{Clean population's names for radiator and other packages}
\usage{
clean_pop_names(x, factor = TRUE)
}
\arguments{
\item{x}{(character string) Population character string.}

\item{factor}{(logical) Default: \code{factor = TRUE}. Will also keep or transform the factor levels.}
}
\description{
Function to clean pop's name
that interfere with some packages
or codes. Space is changed to an underscore \code{_}.
Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
and might be of interest for users.
}
