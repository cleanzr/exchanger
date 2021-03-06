% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coreference_probs.R
\name{coreference_probs}
\alias{coreference_probs}
\alias{coreference_probs.matrix}
\alias{coreference_probs.mcmc}
\alias{coreference_probs.ExchangERFit}
\title{Coreference Probabilities}
\usage{
coreference_probs(x, pairs, rec_ids = NULL)

\method{coreference_probs}{matrix}(x, pairs, rec_ids = NULL)

\method{coreference_probs}{mcmc}(x, pairs, rec_ids = NULL)

\method{coreference_probs}{ExchangERFit}(x, pairs, rec_ids = NULL)
}
\arguments{
\item{x}{an \R object containing samples of clusterings. Currently, methods
are defined for \code{\linkS4class{ExchangERFit}} objects, \code{\link[coda:mcmc]{coda::mcmc}} objects
and \code{\link[base:matrix]{base::matrix}} objects (where rows index samples and columns index
cluster assignments for each record).}

\item{pairs}{record pairs for which an estimate of the posterior
coreference probability is required. The record pairs are assumed to be
stored in a matrix or data.frame, where rows correspond to pairs and
columns correspond to the identifiers for each record in the pair.}

\item{rec_ids}{a vector specifying the complete set of record identifiers
in canonical order. If NULL, the record identifiers are extracted from \code{x}
if available. Otherwise integers are used to identify the records.}
}
\value{
Returns a vector of pairwise match probabilities for each record
pair.
}
\description{
Computes coreference probabilities for a given subset of record pairs
based on samples of clusterings.
}
\examples{
## Initialize a model for entity resolution of RLdata500
library(comparator)
distort_prior <- BetaRV(1,5)
clust_prior <- PitmanYorRP(alpha = GammaRV(1.0, 1.0), d = BetaRV(1.0, 1.0))
attr_params <- list(
  fname_c1 = CategoricalAttribute(distort_prior),
  lname_c1 = CategoricalAttribute(distort_prior),
  bd = CategoricalAttribute(distort_prior),
  by = CategoricalAttribute(distort_prior),
  bm = CategoricalAttribute(distort_prior)
)
model <- exchanger(RLdata500, attr_params, clust_prior)

## Run inference
fit <- run_inference(model, n_samples=100, burnin_interval=100)

## Compute posterior probability that the following record pairs refer to 
## the same entity
pairs <- rbind(c(130,147), c(406,483))
coreference_probs(fit, pairs)

}
