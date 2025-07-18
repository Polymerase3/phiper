% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phip_data-class.R
\name{get_meta}
\alias{get_meta}
\title{Retrieve the metadata list}
\usage{
get_meta(x)
}
\arguments{
\item{x}{A valid \code{phip_data} object.}
}
\value{
A named list.
}
\description{
Accesses the \code{meta} slot, which holds flags such as whether the
table is a full peptide × sample grid, the available outcome columns, etc.
}
