% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/peptide-metadata.R
\name{get_peptide_meta}
\alias{get_peptide_meta}
\title{Retrieve the peptide metadata table into DuckDB, forcing atomic types}
\usage{
get_peptide_meta(force_refresh = FALSE)
}
\arguments{
\item{force_refresh}{Logical. If TRUE, re-downloads and rebuilds the cache.}
}
\value{
A dplyr tbl_dbi pointing to the \code{peptide_meta} table.
}
\description{
\itemize{
\item Downloads the RDS once, sanitizes types (logical, character, numeric),
and writes into a DuckDB cache on disk.
\item Subsequent calls return a lazy tbl_dbi without loading into R memory.
}
}
