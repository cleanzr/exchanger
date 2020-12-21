#' @include ExchangERFit.R
NULL

#' Coreference Probabilities
#'
#' @description
#' Computes coreference probabilities for a given subset of record pairs
#' based on samples of clusterings.
#'
#' @param x an \R object containing samples of clusterings. Currently, methods 
#'   are defined for [`ExchangERFit-class`] objects, [`coda::mcmc`] objects 
#'   and [`base::matrix`] objects (where rows index samples and columns index 
#'   cluster assignments for each record).
#' @param pairs record pairs for which an estimate of the posterior
#'   coreference probability is required. The record pairs are assumed to be
#'   stored in a matrix or data.frame, where rows correspond to pairs and
#'   columns correspond to the identifiers for each record in the pair.
#' @param rec_ids a vector specifying the complete set of record identifiers
#'   in canonical order. If NULL, the record identifiers are extracted from `x` 
#'   if available. Otherwise integers are used to identify the records.
#' @return Returns a vector of pairwise match probabilities for each record 
#'   pair.
#' 
#' @examples 
#' ## Initialize a model for entity resolution of RLdata500
#' library(comparator)
#' distort_prior <- BetaRV(1,5)
#' clust_prior <- PitmanYorRP(alpha = GammaRV(1.0, 1.0), d = BetaRV(1.0, 1.0))
#' attr_params <- list(
#'   fname_c1 = CategoricalAttribute(distort_prior),
#'   lname_c1 = CategoricalAttribute(distort_prior),
#'   bd = CategoricalAttribute(distort_prior),
#'   by = CategoricalAttribute(distort_prior),
#'   bm = CategoricalAttribute(distort_prior)
#' )
#' model <- exchanger(RLdata500, attr_params, clust_prior)
#' 
#' ## Run inference
#' fit <- run_inference(model, n_samples=100, burnin_interval=100)
#' 
#' ## Compute posterior probability that the following record pairs refer to 
#' ## the same entity
#' pairs <- rbind(c(130,147), c(406,483))
#' coreference_probs(fit, pairs)
#' 
#' @rdname coreference_probs
#' @export
coreference_probs <- function(x, pairs, rec_ids=NULL) {
  UseMethod("coreference_probs")
}

#' @rdname coreference_probs
#' @export
coreference_probs.matrix <- function(x, pairs, rec_ids=NULL) {
  if (any(is.na(x))) stop("`x` must not contain NAs")

  if (!is.null(rec_ids)) {
    # Check consistent dimensions
    if (length(rec_ids) != ncol(x))
      stop("length of `rec_ids` must match number of columns in `x`")
  } else {
    # Try to get rec_ids from colnames
    if (!is.null(colnames(x))) {
      rec_ids <- colnames(x)
      mode(pairs) <- "character"
    } else {
      rec_ids <- seq_len(ncol(x))
    }
  }

  original_dim <- dim(pairs)
  pairs <- match(pairs, rec_ids)
  dim(pairs) <- original_dim

  num_pairs <- nrow(pairs)
  out <- numeric(length=num_pairs)
  for (i in seq_len(num_pairs)) {
    out[i] = mean(x[,pairs[i,1]] == x[,pairs[i,2]])
  }

  return(out)
}

#' @rdname coreference_probs
#' @export
coreference_probs.mcmc <- function(x, pairs, rec_ids=NULL) {
  coreference_probs(unclass(x), pairs, rec_ids)
}

#' @rdname coreference_probs
#' @export
coreference_probs.ExchangERFit <- function(x, pairs, rec_ids=NULL) {
  samples <- x@history$links
  if (is.null(samples))
    stop("ExchangERFit does not contain history of the linkage structure along the chain")
  rec_ids <- x@state@rec_ids
  coreference_probs(samples, pairs, rec_ids)
}
