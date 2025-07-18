% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phip_data-class.R
\name{get_comparisons}
\alias{get_comparisons}
\title{Retrieve the comparisons definition table}
\usage{
get_comparisons(x)
}
\arguments{
\item{x}{A valid \code{phip_data} object.}
}
\value{
A tibble with columns: \code{comparison}, \code{group1}, \code{group2},
and \code{variable}.
}
\description{
Returns the two-way comparison specifications stored inside a
\strong{phip_data} object.
}
