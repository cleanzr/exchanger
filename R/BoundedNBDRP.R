#' @include RP.R
NULL

.check_BoundedNBD <- function(object) {
  errors = character()
  if (!is_scalar(object@shape1)) {
    errors <- c(errors, "shape1 must be a scalar")
  }
  if (!is_scalar(object@shape2)) {
    errors <- c(errors, "shape2 must be a scalar")
  }
  if (!is_scalar(object@n)) {
    errors <- c(errors, "n must be a scalar")
  }
  if (!is_scalar(object@alpha)) {
    errors <- c(errors, "alpha must be a scalar")
  }
  if (is_numeric_scalar(object@shape1)) {
    if (object@shape1 <= 0)
      errors <- c(errors, "shape1 must be positive")
  }
  if (is_numeric_scalar(object@shape2)) {
    if (object@shape2 <= 0)
      errors <- c(errors, "shape2 must be positive")
  }
  if (is_numeric_scalar(object@n)) {
    if (object@n < 0) {
      errors <- c(errors, "n must be non-negative")
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

setClass("BoundedNBDRP", 
         slots = c(shape1="numeric", shape2="numeric", n="numeric", alpha="numeric", a="numeric", q="numeric"), 
         validity=.check_BoundedNBD, contains = "RP")

#' Bounded NBD Random Partition
#' 
#' @description
#' TODO: update description.
#' A Kolchin partition model where the number of clusters is drawn from a 
#' shifted negative binomial distribution, and the cluster sizes are drawn from 
#' a distribution which is itself drawn from a Dirichlet process prior with 
#' concentration `alpha` and a shifted beta-binomial base distribution. 
#' Following the paramterization of the negative binomial distribution in 
#' [`stats::NegBinomial`], the shifted negative binomial distribution for the 
#' cluster sizes has parameters \eqn{`size` = a} and \eqn{`prob` = q}. 
#' The shifted beta-binomial base distribution has parameters 
#' \eqn{`shape1` = s_1} and \eqn{`shape2` = s_2} with probability mass function 
#' \deqn{p(K = k) = {n \choose k - 1} \frac{\mathrm{B}(k - 1 + s_1, n - k + 1 + s_2)}{\mathrm{B}(s_1, s_2)}}
#' over the integers \eqn{k = 1, 2, \ldots, n + 1}.
#' 
#' @param shape1,shape2 Shape parameters for the shifted beta-binomial base 
#'  distribution on the cluster sizes. Must be positive numeric scalars.
#' @param n Number of trials for the binomial base distribution on the cluster 
#'   sizes. Must be a non-negative integer scalar.
#' @param alpha Concentration parameter for the Dirichlet process prior on 
#'   the cluster sizes. Must be a positive numeric scalar.
#' @param a Target for number of successful trials for the shifted negative 
#'   binomial prior on the number of clusters. Must be a positive numeric 
#'   scalar.
#' @param q Probability of success for the shifted negative binomial prior on 
#'   the number of clusters. Must be a numeric scalar on the unit interval.
#' @return The constructor `BoundedNBDRP` returns a `BoundedNBDRP` object which 
#'   is a subclass of [`RP-class`].
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
BoundedNBDRP <- function(shape1, shape2, n, alpha, a, q) {
  new("BoundedNBDRP", shape1=shape1, shape2=shape2, n=n, alpha=alpha, a=a, q=q)
}

setMethod("format", "BoundedNBDRP", function(x, ...) {
  shape1 <- format(x@shape1, ...)
  shape2 <- format(x@shape2, ...)
  n <- format(x@n, ...)
  alpha <- format(x@alpha, ...)
  a <- format(x@a, ...)
  q <- format(x@q, ...)
  paste0("BoundedNBDRP(shape1=", shape1, ", shape2=", shape2, ", n=", n, ", alpha=", alpha, ", a=", a, ", q=", q, ")")
})

setMethod("show", "BoundedNBDRP", function(object) {
  writeLines(format(object))
  invisible(object)
})

#' @param x An \R object.
#' @return `is.BoundedNBDRP` returns TRUE if the argument is a `BoundedNBDRP` 
#'   object and FALSE otherwise.
#' 
#' @noRd
#' @keywords internal
is.BoundedNBDRP <- function(x) inherits(x, "BoundedNBDRP")