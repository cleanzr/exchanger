#' Transformed Distance Function
#' 
#' @description
#' Generates a transformed distance function where scaling and thresholding 
#' are applied.
#' 
#' @details 
#' The transformed distance function returns `Inf` if 
#' \eqn{s \times dist_fn(y, x) > \tau}{s * dist_fn(y, x) > τ} where 
#' \eqn{s} is the `scaling_factor` and \eqn{\tau}{τ} is the `threshold`. 
#' Otherwise, the transformed distance function returns 
#' \eqn{s \times dist_fn(y, x)}{s * dist_fn(y, x)}.
#' 
#' @param dist_fn A vectorized distance function that takes a pair of 
#'   vectors y, x and returns their elementwise distances.
#' @param threshold a non-negative numeric. This threshold is applied to 
#'   scaled distances, so that any distances above the threshold are set to 
#'   `Inf`.
#' @param scaling_factor a non-negative numeric. Raw distances are multiplied by 
#'   this scaling factor. Defaults to 1.0.
#' @return a transformed distance function
#' 
#' @seealso 
#' This is useful for preparing distance functions for [`Attribute`].
#' 
#' @export
transform_dist_fn <- function(dist_fn, threshold, scaling_factor = 1.0) {
  if (!is_numeric_scalar(scaling_factor)) 
    stop("scaling_factor must be a scalar numeric")
  if (!is_numeric_scalar(threshold)) 
    stop("threshold must be a scalar numeric")
  if (scaling_factor < 0) 
    stop("scaling_factor must be non-negative")
  if (threshold < 0) 
    stop("threshold must be non-negative")
  tdist_fn <- function(y, x) {
    dist_trunc <- scaling_factor * dist_fn(y, x)
    return(ifelse(dist_trunc <= threshold, dist_trunc, Inf))
  }
  return(tdist_fn)
}