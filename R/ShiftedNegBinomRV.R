#' @include RV.R
NULL

.check_ShiftedNegBinomRV <- function(object) {
  errors = character()
  if (!is_scalar(object@size)) {
    errors <- c(errors, "size must be a scalar")
  }
  if (!is_scalar(object@prob)) {
    errors <- c(errors, "prob must be a scalar")
  }
  if (object@size <= 0) {
    errors <- c(errors, "size must be positive")
  }
  if (object@prob <= 0 || object@prob > 1) {
    errors <- c(errors, "prob must be on the interval (0, 1]")
  }
  if (!is.finite(object@size)) {
    errors <- c(errors, "size must be finite")
  }
  if (length(errors)==0) TRUE else errors
}

setClass("ShiftedNegBinomRV", slots = c(size = "numeric", prob = "numeric"), 
         validity=.check_ShiftedNegBinomRV, contains = "RV")

#' Shifted Negative Binomial Random Variable
#' 
#' @description
#' Represents a random variable with a negative binomial distribution, 
#' shifted by 1 so that the support is the positive integers.
#' 
#' @details
#' The parameterization adopted here matches [`stats::NegBinomial`]. 
#' The distribution is over the number of failures which occur in a sequence 
#' of Bernoulli trials, before a target number of successes occur, given that 
#' the number of failures is at least one.
#' If \eqn{x} denotes the number of failures, \eqn{size = r} denotes the 
#' number of success by and \eqn{prob = p} denotes the probability of success 
#' in each trial, then the density of \eqn{x} is
#' \deqn{p(X = x) = \frac{\Gamma(x + r - 1)}{(x - 1)! \Gamma(r)} p^n (1 - p)^{x - 1}}{p(X = x) = (Γ(x + r - 1))/((x - 1)! Γ(r)) p^n (1 - p)^{x - 1}}
#' for \eqn{x = 1, 2, 3, ...}, \eqn{n > 0} and \eqn{0 < p \leq 1}{0 < p <= 1}.
#' 
#' @param size Target number of successes before experiment is stopped.
#' @param prob Probability of success in each trial. Must be on the unit 
#' interval \eqn{(0, 1]}
#' @return The `ShiftedNegBinomRV` constructor returns a `ShiftedNegBinomRV` 
#' object, which is a subclass of [`RV-class`].
#' 
#' @examples
#' x <- ShiftedNegBinomRV(100, 0.1)
#' mean(x)
#' 
#' @seealso Other random variables defined in this package 
#' include [`BetaRV`], [`GammaRV`] and [`DirichletRV`].
#' 
#' @export
ShiftedNegBinomRV <- function(size, prob) {
  new("ShiftedNegBinomRV", size=size, prob=prob)
}

setMethod("format", "ShiftedNegBinomRV", function(x, ...) {
  paste0("ShiftedNegBinomRV(size=", x@size, ", prob=", x@prob, ")")
})

setMethod("show", "ShiftedNegBinomRV", function(object) {
  writeLines(format(object))
  invisible(object)
})

#' @param x an \R object
#' @return \code{is.ShiftedNegBinomRV} returns TRUE if the argument is a 
#' `ShiftedNegBinomRV` object and FALSE otherwise.
#' 
#' @noRd
#' @keywords internal
is.ShiftedNegBinomRV <- function(x) inherits(x, "ShiftedNegBinomRV")

#' @describeIn mean Specialization for [`ShiftedNegBinomRV`]
#' @export
setMethod("mean", signature(x = "ShiftedNegBinomRV"), 
  function(x) x@size * (1 - x@prob) / x@prob + 1)

#' @describeIn var Specialization for [`ShiftedNegBinomRV`]
#' @export
setMethod("var", signature(x = "ShiftedNegBinomRV"), 
  function(x) x@size * (1 - x@prob) / (x@prob)^2)

#' @describeIn drv Specialization for [`ShiftedNegBinomRV`]
#' @importFrom stats dnbinom
#' @export
setMethod("drv", signature(x = "numeric", rv = "ShiftedNegBinomRV"), 
  function(x, rv, log=FALSE, ...) dnbinom(x - 1, rv@size, rv@prob, log=log))