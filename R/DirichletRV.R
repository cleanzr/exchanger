#' @include RV.R
NULL

.check_DirichletRV <- function(object) {
  errors = character()
  if (any(object@alpha <= 0)) {
    errors <- c(errors, "all entries in alpha must be positive")
  }
  if (!any(is.finite(object@alpha))) {
    errors <- c(errors, "all entries in alpha must be finite")
  }
  if (length(errors)==0) TRUE else errors
}

setClass("DirichletRV", slots = c(alpha = "numeric"), 
         validity=.check_DirichletRV, contains = "RV")

#' Dirichlet-distributed Random Vector
#' 
#' @description
#' Represents a random vector with a Dirichlet distribution.
#' 
#' @details
#' With concentration parameter vector \eqn{alpha = a} of length \eqn{K}, 
#' the density is
#' \deqn{p(X = x) = \frac{\prod_{i=1}^{K} \Gamma(\alpha_i)}{\Gamma(\sum_{i=1}^{K} \alpha_i)} \prod_{i=1}^{K} x_i^{\alpha_i-1} \exp(- b x).}{p(X = x) = (prod_i Γ(ɑ_i))/(Γ(sum_i ɑ_i)) prod_i x_i^{ɑ_i-1} exp(- b x).}
#' for vector \eqn{x} in the \eqn{K - 1} simplex.
#' 
#' @param alpha Positive concentration parameter vector. If a scalar is 
#' passed, this value is repeated for each entry in the vector.
#' 
#' @return The `DirichletRV` constructor returns a `DirichletRV` object, 
#' which is a subclass of [`RV-class`].
#' 
#' @examples
#' x <- DirichletRV(c(1,1,1))
#' mean(x)
#' 
#' @seealso Other random variables defined in this package 
#' include [`BetaRV`], [`GammaRV`] and [`ShiftedNegBinomRV`].
#' 
#' @export
DirichletRV <- function(alpha) {
  new("DirichletRV", alpha=alpha)
}

setMethod("format", "DirichletRV", function(x, ...) {
  alpha <- toString(x@alpha, width=20)
  paste0("DirichletRV(alpha=", alpha, ")")
})

setMethod("show", "DirichletRV", function(object) {
  writeLines(format(object))
  invisible(object)
})

#' @param x An \R object
#' @return `is.DirichletRV` returns TRUE if the argument is a `BetaRV` 
#' object and FALSE otherwise.
#' 
#' @noRd
#' @keywords internal
is.DirichletRV <- function(x) inherits(x, "DirichletRV")

#' @describeIn mean Specialization for [`DirichletRV`]
#' @export
setMethod("mean", signature(x = "DirichletRV"), 
  function(x) {
    if (is_scalar(x@alpha)) {
      stop("cannot compute mean with scalar alpha")
    }
    x@alpha/sum(x@alpha)
  })

#' @describeIn var Specialization for [`DirichletRV`]
#' @export
setMethod("var", signature(x = "DirichletRV"), 
  function(x) {
    if (is_scalar(x@alpha)) {
      stop("cannot compute variance with scalar alpha")
    }
    s <- sum(x@alpha)
    alpha_norm <- x@alpha / s
    sigma <- -outer(alpha_norm, alpha_norm)
    diag(sigma) <- alpha_norm * (1.0 - alpha_norm)
    sigma <- sigma / (s + 1.0)
    return(sigma)
  })