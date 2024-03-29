#' @include utils.R
NULL

setClass("AttributeIndex", 
         slots = c(
           domain = "character", 
           probs = "numeric", 
           close_values = "list", 
           max_exp_factor = "numeric", 
           n_missing = "numeric"
          ), 
          validity = function(object) {
            errors <- character()
            domSize <- length(object@domain)
            if (!is.vector(object@probs)) {
              errors <- c(errors, "probs must be a vector")
            }
            if (!approx_equal(sum(object@probs), 1.)) {
              errors <- c(errors, "probs must be normalized")
            }
            if (!is.vector(object@domain)) {
              errors <- c(errors, "domain must be a vector")
            }
            if (length(object@probs) != domSize) {
              errors <- c(errors, "probs and domain must be the same length")
            }
            if (length(object@close_values) != 0 & length(object@close_values) != domSize) {
              errors <- c(errors, "close_values and domain must be the same length")
            }
            if (length(object@max_exp_factor) != 0 & length(object@close_values) != domSize) {
              errors <- c(errors, "max_exp_factor and domain must be the same length")
            }
            if (!all(object@max_exp_factor >= 0.0 & object@max_exp_factor <= 1.0)) {
              errors <- c(errors, "entries in max_exp_factor must be in the unit interval [0,1]")
            }
            if (!is_scalar(object@n_missing) | object@n_missing < 0) {
              errors <- c(errors, "n_missing must be a non-negative scalar")
            }
            if (length(errors)==0) TRUE else errors
          })

#' Attribute Index
#' 
#' @description
#' Builds an `AttributeIndex` object for a collection of observed attribute 
#' values. 
#' 
#' @details
#' The index includes quantities that are required for running inference 
#' such as the empirical distribution over the domain and an index for 
#' conducting range queries. This is an internal function and should not be 
#' required unless constructing an [`ExchangERModel-class`] object manually.
#' 
#' @param values A vector of attribute values (e.g. a column from the source 
#'   data) represented as a factor.
#' @param dist_fn A vectorized distance function or NULL (for a constant 
#'   distance function).
#' @param exclude_entity_value A boolean specifying whether the entity 
#'   value \eqn{y} should be excluded from the support of the distortion 
#'   distribution for record value \eqn{x}.
#' @return An `AttributeIndex` object.
#' 
#' @note This function supports parallel evaluation via the `future` package. 
#'   To configure the evaluation strategy, see `future::plan`.
#' 
#' @keywords internal
#' @export
AttributeIndex <- function(values, dist_fn, exclude_entity_value) {
  if (!is.null(dim(values))) stop("values must be a 1D array")
  if (!is.factor(values)) stop("values must be a factor")
  # TODO: check for duplicate levels?
  
  val_ids <- unclass(values)
  domain <- attr(val_ids, "levels")
  attributes(val_ids) <- NULL
  probs <- as.numeric(table(values)) + 1 # additive smoothing in case a value is not observed
  probs <- probs/sum(probs)
  n_missing <- sum(is.na(val_ids))
  
  if (is.null(dist_fn)) {
    new("AttributeIndex", domain=domain, probs=probs, close_values=list(), 
        max_exp_factor=numeric(), n_missing=n_missing)
  } else {
    res <- compute_close_values(domain, probs, dist_fn, !exclude_entity_value)
    new("AttributeIndex", domain=domain, probs=probs, close_values=res$close_values, 
        max_exp_factor=res$max_exp_factor, n_missing=n_missing)
  }
}

#' @param x an \R object
#' @param parallel Whether to parallelize the computation
#' @return `is.AttributeIndex` returns TRUE if the argument is an 
#' `AttributeIndex` object and FALSE otherwise.
#' 
#' @noRd
#' @keywords internal
is.AttributeIndex <- function(x) inherits(x, "AttributeIndex")


#' Compute index of close values for each value in the domain
#' 
#' @param domain A character vector containing the unique values in the domain.
#' @param probs A numeric vector containing the probability mass associated 
#'   with each value in the domain.
#' @param dist_fn A vectorized distance function.
#' @param include_self Whether to include the value itself in the list of 
#'   close values.
#' @return A list containing two entries:
#'   - `close_values`: a list of close values for each value in `domain`
#'   - `max_exp_factor`: a vector containing max_x exp(-dist_fn(y, x)) for each 
#'     value in the domain y
#' 
#' @importFrom future.apply future_lapply
#' 
#' @keywords internal
#' @noRd
compute_close_values <- function(domain, probs, dist_fn, include_self) {
  close_values <- future_lapply(seq_along(domain), function(i) {
    dists <- dist_fn(domain[i], domain)
    if (!include_self) { 
      dists[i] <- Inf # exclude i from the distortion distribution
    }
    min_dist <- min(dists)
    finite_dist_val_ids <- which(is.finite(dists))
    dists <- dists[finite_dist_val_ids]
    max_dist <- if (length(dists)) max(dists) else NA
    #distort_probs <- probs[finite_dist_val_ids] * exp(-dists)
    distort_probs <- exp(-dists)
    distort_probs <- distort_probs / sum(distort_probs)
    return(list(valIds=finite_dist_val_ids, probs=distort_probs, min_dist=min_dist, max_dist=max_dist))
  }, future.seed = TRUE)
  min_dists <- sapply(close_values, function(x) x$min_dist) 
  max_dists <- sapply(close_values, function(x) x$max_dist)
  max_dist <- if (all(is.na(max_dists))) 0 else max(max_dists, na.rm=TRUE)
  max_exp_factor <- ifelse(is.finite(min_dists), exp(-0.5*min_dists/max_dist), 0)
  return(list(close_values=close_values, max_exp_factor=max_exp_factor))
}