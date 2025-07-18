% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phip_data-class.R
\name{new_phip_data}
\alias{new_phip_data}
\title{Construct a \strong{phip_data} object}
\usage{
new_phip_data(
  data_long,
  comparisons,
  backend = c("memory", "duckdb", "arrow"),
  peptide_library = NULL,
  meta = list()
)
}
\arguments{
\item{data_long}{A tidy data frame (or \code{tbl_lazy}) with one row per
\code{peptide_id} x \code{sample_id} combination. \strong{Required.}}

\item{comparisons}{A data frame describing two-way contrasts
(\code{comparison}, \code{group1}, \code{group2}, \code{variable});
defaults to an empty tibble if \code{NULL}.}

\item{backend}{Character string specifying the storage engine to use:
\code{"memory"}, \code{"duckdb"}, or \code{"arrow"}.  If \code{NULL}
the implicit default is \code{"duckdb"}.}

\item{peptide_library}{A data frame with one row per \code{peptide_id}
and its annotations.  If \code{NULL}, the package’s current default
library is used.}

\item{meta}{Optional named list of metadata flags to pre-populate the
\code{meta} slot (rarely needed by users).}
}
\value{
An object of class \code{"phip_data"}.
}
\description{
Creates a fully-validated S3 object that bundles the tidy
PhIP-Seq counts (\code{data_long}), optional comparison definitions, a
peptide-library annotation table, and backend metadata.  The function
performs a minimal sanity check on \emph{comparisons} and normalises the chosen
storage
\emph{backend} before returning the object (validation of the data itself
happens via \code{validate_phip_data()} helper).
}
\examples{
\dontrun{
## minimal constructor call
pd <- new_phip_data(
  data_long = tidy_counts,
  comparisons = NULL,
  backend = "duckdb",
  peptide_library = default_peptides
)
}
## list available backends
c("memory", "duckdb", "arrow")

}
