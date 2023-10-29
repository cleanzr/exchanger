#' @include RP.R
NULL

.check_ESCD <- function(object) {
  errors = character()
  if (!is_scalar(object@p) && !is.BetaRV(object@p)) {
    errors <- c(errors, "p must be a scalar or a BetaRV")
  }
  if (!is_scalar(object@n) & !is.GammaRV(object@n)) {
    errors <- c(errors, "n must be a scalar or a GammaRV")
  }
  if (!is_scalar(object@alpha)) {
    errors <- c(errors, "alpha must be a scalar")
  }
  if (is_numeric_scalar(object@p)) {
    if (object@p < 0 || object@p > 1)
      errors <- c(errors, "p must be on the interval [0, 1]")
  }
  if (is_numeric_scalar(object@n)) {
    if (object@n <= 0) {
      errors <- c(errors, "n must be strictly positive")
    }  
  }
  if (is_numeric_scalar(object@alpha)) {
    if (object@alpha <= 0) {
      errors <- c(errors, "alpha must be strictly positive")
    }
  }
  if (length(errors)==0) TRUE else errors
}

setClass("ESCDRP", 
         slots = c(p="ANY", n="ANY", alpha="numeric"), 
         validity=.check_ESCD, contains = "RP")

#' ESC-D Random Partition
#' 
#' @description
#' An ESC model where the cluster sizes are drawn from a distribution on the 
#' positive integers, which is itself drawn from a Dirichlet process with a 
#' shifted negative binomial base distribution with concentration parameter 
#' \eqn{`alpha`}. 
#' The parameterization of the shifted negative binomial distribution matches 
#' [`stats::NegBinomial`]. 
#' For \eqn{`size` = n} and \eqn{`prob` = p}, the probability mass associated 
#' with cluster size \eqn{x} is:
#' \deqn{p(X = x) = \frac{\Gamma(x + n - 1)}{(x - 1)! \Gamma(n)} p^n (1 - p)^{x - 1}}{p(X = x) = (Γ(x + n - 1))/((x - 1)! Γ(n)) p^n (1 - p)^{x - 1}}
#' for \eqn{x = 1, 2, 3, ...}, \eqn{n > 0} and \eqn{0 < p \leq 1}{0 < p <= 1}.
#' Hyperpriors on the `p` and `n` parameters are supported.
#' 
#' @param p Probability of success. Must be a numeric scalar on the unit 
#'   interval \eqn{(0, 1]} or a [`BetaRV`] object.
#' @param n Target for number of successful trials, or dispersion parameter 
#'   (the shape parameter of the gamma mixing distribution). Must be a positive 
#'   numeric scalar or a [`GammaRV`] object.
#' @param alpha Concentration parameter for the Dirichlet process. Must be a 
#'   positive numeric scalar.
#' @return The constructor `ESCDRP` returns a `ESCDRP` object which is a 
#'   subclass of [`RP-class`].
#' 
#' @seealso 
#' Other random partitions defined in this package include [`PitmanYorRP`] and 
#' [`EwensRP`].
#' 
#' @references
#' Brenda Betancourt, Giacomo Zanella & Rebecca C. Steorts (2020): Random
#' Partition Models for Microclustering Tasks, _Journal of the American 
#' Statistical Association_, DOI: 10.1080/01621459.2020.1841647.
#' 
#' @export
ESCDRP <- function(p, n, alpha) {
  new("ESCDRP", p=p, n=n, alpha=alpha)
}

setMethod("format", "ESCDRP", function(x, ...) {
  p <- format(x@p, ...)
  n <- format(x@n, ...)
  alpha <- format(x@alpha, ...)
  paste0("ESCDRP(p=", p, ", n=", n, ", alpha=", alpha, ")")
})

setMethod("show", "ESCDRP", function(object) {
  writeLines(format(object))
  invisible(object)
})

#' @param x An \R object.
#' @return `is.ESCDRP` returns TRUE if the argument is a `ESCDRP` object and 
#'   FALSE otherwise.
#' 
#' @noRd
#' @keywords internal
is.ESCDRP <- function(x) inherits(x, "ESCDRP")