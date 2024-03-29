% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ExchangERModel.R
\docType{class}
\name{ExchangERModel-class}
\alias{ExchangERModel-class}
\alias{ExchangERModel}
\title{ExchangER Model Class}
\description{
This class records the state of all model parameters. It also stores
indexing data structures for each attribute (\code{attr_indices}) which
are used to accelerate sampling.
}
\section{Slots}{

\describe{
\item{\code{iteration}}{A non-negative integer. Counts the number of successive
applications of the Markov transition operator. This is not necessarily
equal to the number of samples if thinning is applied.}

\item{\code{attr_params}}{A named list of \code{\link{Attribute}} objects. Each entry
in the list specifies the model parameters for an entity attribute,
and should match one of the row names (attributes) in
\code{rec_attrs}.}

\item{\code{file_ids}}{An integer vector of file identifiers. The i-th entry
specifies the file/source of the record in the i-th column of
\code{rec_attrs}. The file identifiers should be in the set
\eqn{\{1, \ldots, n_files\}}{{1, ..., n_files}}.}

\item{\code{rec_ids}}{A character vector of record identifiers.}

\item{\code{rec_attrs}}{An integer matrix representation of the observed
records. Rows correspond to attributes and columns correspond to records.
The attributes are encoded using integer ids, which index the domain of
each attribute (starting at 1).}

\item{\code{rec_distortions}}{An integer matrix of distortion indicators
of the same dimension as \code{rec_attrs}. A value of 1 (0) means that
the corresponding attribute value in \code{rec_attrs} is distorted
(not distorted).}

\item{\code{ent_attrs}}{An integer matrix representation of the
true entity attributes. Rows correspond to attributes and columns
correspond to entities. The encoding for the attributes (rows) matches the
encoding used in \code{rec_attrs}. The entity identifier of the i-th
column is recorded in the i-th entry of \code{ent_ids}.}

\item{\code{ent_ids}}{An integer vector of unique entity identifiers. These
identifiers are used to encode the \code{links} between records and entities.}

\item{\code{links}}{An integer vector of entity identifiers for each record in
\code{rec_attrs}. The i-th entry gives the linked entity for the
i-th record (column) of \code{rec_attrs}.}

\item{\code{distort_probs}}{A numeric matrix of distortion probabilities.
Rows correspond to files/sources and columns correspond to attributes.
The (i,j)-th entry stores the distortion probability for the i-th file
id and the j-th attribute.}

\item{\code{clust_params}}{A \code{\linkS4class{RP}} object where all parameters
(including random variables) are realized. If all parameters are
fixed, this slot will be identical to \code{clust_prior}.}

\item{\code{attr_indices}}{A named list of \code{\link{AttributeIndex}} objects. Each
entry in the list corresponds to an attribute in \code{attr_params}.}

\item{\code{clust_prior}}{A \code{\linkS4class{RP}} object. Represents the prior
over the clustering of records into entities (a.k.a. linkage structure).}

\item{\code{distort_dist_concs}}{A numeric vector of positive concentration
parameters. Each entry specifies the Dirichlet Process concentration
parameter for one of the row names (attributes) in \code{rec_attrs}.}
}}

\seealso{
The \code{\link{exchanger}} function should be used to initialize an ER model
from data.
}
