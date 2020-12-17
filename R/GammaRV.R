#' @include RV.R
NULL

.check_GammaRV <- function(object) {
  errors = character()
  if (!is_scalar(object@shape)) {
    errors <- c(errors, "shape must be a scalar")
  }
  if (!is_scalar(object@rate)) {
    errors <- c(errors, "rate must be a scalar")
  }
  if (object@shape <= 0) {
    errors <- c(errors, "shape must be positive")
  }
  if (object@rate <= 0) {
    errors <- c(errors, "rate must be positive")
  }
  if (!is.finite(object@shape)) {
    errors <- c(errors, "shape must be finite")
  }
  if (!is.finite(object@rate)) {
    errors <- c(errors, "rate must be finite")
  }
  if (length(errors)==0) TRUE else errors
}

setClass("GammaRV", slots = c(shape = "numeric", rate = "numeric"), 
         validity=.check_GammaRV, contains = "RV")

#' Gamma-distributed Random Variable
#' 
#' @description
#' Represents a random variable with a Gamma distribution.
#' 
#' @details `GammaRV` constructs a `GammaRV` object.
#' The parameterization used here is in terms of a shape parameter 
#' \eqn{shape = a > 0} and rate parameter \eqn{rate = b > 0}.
#' The density is
#' \deqn{p(X = x) = \frac{b^a}{\Gamma(a)} x^{a-1} \exp(- b x).}{p(X = x) = (b^a)/(Γ(a)) x^{a-1} exp(- b x).}
#' for \eqn{x \geq 0}{x >= 0}.
#' 
#' @param shape Positive real-valued shape hyperparameter
#' @param rate Positive real-valued rate hyperparameter
#' @return The `GammaRV` constructor returns a `GammaRV` object, 
#' which is a subclass of [`RV-class`].
#' 
#' @examples
#' x <- GammaRV(1, 1)
#' mean(x)
#' 
#' @seealso Other random variables defined in this package 
#' include [`BetaRV`], [`DirichletRV`] and [`ShiftedNegBinomRV`].
#' 
#' @export
GammaRV <- function(shape, rate) {
  new("GammaRV", shape=shape, rate=rate)
}

#' @param x An \R object
#' @return \code{is.GammaRV} returns TRUE if the argument is a `BetaRV` 
#' object and FALSE otherwise.
#' 
#' @noRd
#' @keywords internal
is.GammaRV <- function(x) inherits(x, "GammaRV")

#' @describeIn mean Specialization for [`GammaRV`]
#' @export
setMethod("mean", signature(x = "GammaRV"), 
  function(x) x@shape/x@rate)

#' @describeIn var Specialization for [`GammaRV`]
#' @export
setMethod("var", signature(x = "GammaRV"), 
  function(x) x@shape/(x@rate)^2)

#' @describeIn drv Specialization for [`GammaRV`]
#' @importFrom stats dgamma
#' @export
setMethod("drv", signature(x = "numeric", rv = "GammaRV"), 
  function(x, rv, log=FALSE, ...) dgamma(x, shape=rv@shape, rate=rv@rate, log=log))