#' @include RP.R
NULL

.check_NBD <- function(object) {
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
  if (is_numeric_scalar(object@a)) {
    if (object@a <= 0) 
      errors <- c(errors, "a must be strictly positive")
  }
  if (is_numeric_scalar(object@q)) {
    if (object@q < 0 || object@q > 1)
      errors <- c(errors, "q must be on the interval [0, 1]")
  }
  if (length(errors)==0) TRUE else errors
}

setClass("NBDRP", 
         slots = c(p="ANY", n="ANY", alpha="numeric", a="numeric", q="numeric"), 
         validity=.check_NBD, contains = "RP")

#' NBD Random Partition
#' 
#' @description
#' TODO: update description.
#' A Kolchin partition model where the number of clusters is drawn from a 
#' shifted negative binomial distribution, and the cluster sizes are drawn from 
#' a distribution which is itself drawn from a Dirichlet process prior with 
#' concentration `alpha` and a shifted negative binomial base distribution. 
#' Following the paramterization of the negative binomial distribution in 
#' [`stats::NegBinomial`]: 
#' * the shifted negative binomial distribution on the number of clusters has 
#'   parameters \eqn{`size` = n} and \eqn{`prob` = p}. 
#' * the shifted negative binomial distribution for the cluster sizes has
#'   parameters \eqn{`size` = a} and \eqn{`prob` = q}. 
#' Hyperpriors on the `p` and `n` parameters are supported.
#' 
#' @param p Probability of success for the shifted negative binomial base 
#'   distribution on the cluster sizes. Must be a numeric scalar on the unit 
#'   interval \eqn{(0, 1]} or a [`BetaRV`] object.
#' @param n Target for number of successful trials for the shifted negative 
#'   binomial base distribution on the cluster sizes. Must be a positive 
#'   numeric scalar or a [`GammaRV`] object.
#' @param alpha Concentration parameter for the Dirichlet process prior on 
#'   the cluster sizes. Must be a positive numeric scalar.
#' @param a Target for number of successful trials for the shifted negative 
#'   binomial prior on the number of clusters. Must be a positive numeric 
#'   scalar.
#' @param q Probability of success for the shifted negative binomial prior on 
#'   the number of clusters. Must be a numeric scalar on the unit interval.
#' @return The constructor `NBDRP` returns a `NBDRP` object which is a 
#'   subclass of [`RP-class`].
#' 
#' @seealso 
#' Other random partitions defined in this package include [`PitmanYorRP`] and 
#' [`EwensRP`].
#' 
#' @references
#' Brenda Betancourt, Giacomo Zanella, Hannah Wallach, Jeffrey W. Miller, Abbas 
#' Zaidi & Rebecca C. Steorts (2016): Flexible Models for Microclustering with 
#' Application to Entity Resolution, _Advances in Neural Information Processing Systems_, 29 
#' URL: https://proceedings.neurips.cc/paper_files/paper/2016/file/670e8a43b246801ca1eaca97b3e19189-Paper.pdf
#' 
#' @export
NBDRP <- function(p, n, alpha, a, q) {
  new("NBDRP", p=p, n=n, alpha=alpha, a=a, q=q)
}

setMethod("format", "NBDRP", function(x, ...) {
  p <- format(x@p, ...)
  n <- format(x@n, ...)
  alpha <- format(x@alpha, ...)
  a <- format(x@a, ...)
  q <- format(x@q, ...)
  paste0("NBDRP(p=", p, ", n=", n, ", alpha=", alpha, ", a=", a, ", q=", q, ")")
})

setMethod("show", "NBDRP", function(object) {
  writeLines(format(object))
  invisible(object)
})

#' @param x An \R object.
#' @return `is.NBDRP` returns TRUE if the argument is a `NBDRP` object and 
#'   FALSE otherwise.
#' 
#' @noRd
#' @keywords internal
is.NBDRP <- function(x) inherits(x, "NBDRP")