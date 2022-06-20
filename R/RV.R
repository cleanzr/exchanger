# Define a generic random variable class

#' Random Variable
#' 
#' @description 
#' Virtual class for a random variable
#' 
#' @export
setClass("RV")

#' @param x An \R object.
#' @return `is.RV` returns TRUE if the argument is a [RV-class] object and 
#' FALSE otherwise.
#' 
#' @keywords internal
#' @noRd
is.RV <- function(x) inherits(x, "RV")

#' Mean of a Random Variable
#' 
#' @description 
#' Computes the mean (expected value) of a random variable
#' 
#' @param x a [`RV-class`] object
#' @param ... other parameters passed on to other methods.
#' 
#' @return 
#' Returns the mean of the random variable
#' 
#' @examples
#' x <- BetaRV(1, 1)
#' mean(x)
#' 
#' @export
setGeneric("mean", function(x, ...) standardGeneric("mean"), 
           signature = c("x"))

#' Standard Deviation of a Random Variable
#' 
#' @description 
#' Computes the standard deviation of a random variable
#' 
#' @param x a [`RV-class`] object
#' @param ... other parameters passed on to other methods.
#' 
#' @return 
#' Returns the mean of the random variable
#' 
#' @examples
#' x <- BetaRV(1, 1)
#' sd(x)
#' 
#' @export
setGeneric("sd", function(x, ...) {
  sigma <- var(x)
  if (is_numeric_scalar(sigma))
    return(sqrt(sigma))
  else
    stop("sd not defined")
}, signature = c("x"))

#' Variance-Covariance of a Random Variable
#' 
#' @description 
#' Computes the variance-covariance matrix of a random variable
#' 
#' @param x a [`RV-class`] object
#' @param ... other parameters passed on to other methods.
#' 
#' @return 
#' Returns the variance if the random variable is a scalar, or a 
#' variance-covariance matrix if the random variable is a vector.
#' 
#' @examples
#' x <- BetaRV(1, 1)
#' var(x)
#' 
#' @export
setGeneric("var", function(x, ...) standardGeneric("var"), 
           signature = c("x"))

#' Density Function Associated with a Random Variable
#' 
#' @description 
#' Computes the density of a random variable at a point
#' 
#' @param x a vector of quantiles
#' @param rv a [`RV-class`] object
#' @param log a logical; if TRUE, probabilities are given as log(p).
#' @param ... other parameters passed on to other methods.
#' 
#' @return 
#' Returns the density of the random variable at `x`.
#' 
#' @examples
#' x <- BetaRV(1, 1)
#' drv(x)
#' 
#' @export
setGeneric("drv", function(x, rv, log, ...) standardGeneric("drv"), 
           signature = c("x", "rv"))

#' Generate a Random Variable
#'
#' @param rv a [`RV-class`] object
#' @param ... other parameters passed on to other methods.
#' 
#' @return 
#' Returns the generated value of the random variable.
#' 
#' @examples
#' x <- BetaRV(1, 1)
#' rrv(x)
#' 
#' @export
setGeneric("rrv", function(rv, ...) standardGeneric("rrv"), 
           signature = c("rv"))