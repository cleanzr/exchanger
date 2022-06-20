#' @include RV.R
NULL

.check_BetaRV <- function(object) {
  errors = character()
  if (!is_scalar(object@shape1)) {
    errors <- c(errors, "shape1 must be a scalar")
  }
  if (!is_scalar(object@shape2)) {
    errors <- c(errors, "shape2 must be a scalar")
  }
  if (object@shape1 <= 0) {
    errors <- c(errors, "shape1 must be positive")
  }
  if (object@shape2 <= 0) {
    errors <- c(errors, "shape2 must be positive")
  }
  if (!is.finite(object@shape1)) {
    errors <- c(errors, "shape1 must be finite")
  }
  if (!is.finite(object@shape2)) {
    errors <- c(errors, "shape2 must be finite")
  }
  if (length(errors)==0) TRUE else errors
}

setClass("BetaRV", slots = c(shape1 = "numeric", shape2 = "numeric"), 
         validity=.check_BetaRV, contains = "RV")

#' Beta Random Variable
#'
#' @description
#' Represents a random variable with a Beta distribution.
#'
#' @details 
#' The density of the Beta distribution with parameters \eqn{shape1 = a} and 
#' \eqn{shape2 = b} is given by
#' \deqn{p(X = x) = \frac{\Gamma(a+b)}{\Gamma(a) \Gamma(b)} x^{a-1} (1-x)^{b-1}.}{p(X = x) = Gamma(a+b)/(Gamma(a) Gamma(b)) x^{a-1} (1-x)^{b-1}.}
#' for \eqn{0 < x < 1}, \eqn{a > 0} and \eqn{b > 0}.
#'
#' @param shape1,shape2 Positive shape parameters of the Beta distribution
#' @return The `BetaRV` constructor returns a `BetaRV` object, 
#' which is a subclass of [`RV-class`].
#' 
#' @examples
#' # A uniform distribution on the unit interval 
#' x <- BetaRV(1, 1)
#' mean(x)
#' 
#' @seealso Other random variables defined in this package 
#' include [`GammaRV`], [`DirichletRV`] and [`ShiftedNegBinomRV`].
#' 
#' @export
BetaRV <- function(shape1, shape2) {
  new("BetaRV", shape1=shape1, shape2=shape2)
}

setMethod("format", "BetaRV", function(x, ...) {
  paste0("BetaRV(shape1=", x@shape1, ", shape2=", x@shape2, ")")
})

setMethod("show", "BetaRV", function(object) {
  writeLines(format(object))
  invisible(object)
})

#' @param x An \R object
#' @return `is.BetaRV` returns TRUE if the argument is a `BetaRV` 
#' object and FALSE otherwise
#' 
#' @noRd
#' @keywords internal
is.BetaRV <- function(x) inherits(x, "BetaRV")

#' @describeIn mean Specialization for [`BetaRV`]
#' @export
setMethod("mean", signature(x = "BetaRV"), 
  function(x) x@shape1/(x@shape1 + x@shape2))

#' @describeIn var Specialization for [`BetaRV`]
#' @export
setMethod("var", signature(x = "BetaRV"), 
  function(x) {
    total <- x@shape1 + x@shape2
    x@shape1 * x@shape2/(total^2 * (total + 1))
  })

#' @describeIn drv Specialization for [`BetaRV`]
#' @importFrom stats dbeta
#' @export
setMethod("drv", signature(x = "numeric", rv = "BetaRV"), 
  function(x, rv, log=FALSE, ...) dbeta(x, rv@shape1, rv@shape2, log=log))

#' @describeIn rrv Specialization for [`BetaRV`]
#' @importFrom stats rbeta
#' @export
setMethod("rrv", signature(rv = "BetaRV"), 
          function(rv, ...) rbeta(1, shape1=rv@shape1, shape2=rv@shape2))