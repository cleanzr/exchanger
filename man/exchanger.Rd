% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ExchangERModel.R
\name{exchanger}
\alias{exchanger}
\title{Initialize an ExchangERModel}
\usage{
exchanger(
  data,
  attr_params,
  clust_prior,
  file_id_colname = NULL,
  rec_id_colname = NULL
)
}
\arguments{
\item{data}{A data frame or matrix of observed records. Column names are
are used to refer to particular columns in the other arguments.}

\item{attr_params}{A named list of \code{\link{Attribute}} objects. Each entry
in the list specifies the model parameters for an entity attribute,
and should match one of the column names (attributes) in
\code{data}.}

\item{clust_prior}{A \code{\linkS4class{RP}} object which specifies the
prior over the clustering (a.k.a. linkage structure).}

\item{file_id_colname}{Column name in \code{data} that contains file/source
identifiers for the records. If NULL, the records are assumed to be
from a single file/source.}

\item{rec_id_colname}{Column name in \code{data} that contains unique record
identifiers. If NULL, use row names to identify records.}
}
\value{
An \code{\linkS4class{ExchangERModel}} object. This can be used as the starting
point for inference using \code{\link{run_inference}} function.
}
\description{
Initializes an \code{\linkS4class{ExchangERModel}} given observed data and model
hyperparameters.
}
\examples{
## Initialize a model for entity resolution of RLdata500
library(comparator)
distort_prior <- BetaRV(1,5)
clust_prior <- PitmanYorRP(alpha = GammaRV(1.0, 1.0), d = BetaRV(1.0, 1.0))
attr_params <- list(
  fname_c1 = Attribute(Levenshtein(), distort_prior),
  lname_c1 = Attribute(Levenshtein(), distort_prior),
  bd = CategoricalAttribute(distort_prior),
  by = CategoricalAttribute(distort_prior),
  bm = CategoricalAttribute(distort_prior)
)
model <- exchanger(RLdata500, attr_params, clust_prior)

}
