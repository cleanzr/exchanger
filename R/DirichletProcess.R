#' @include GammaRV.R
NULL

.check_DirichletProcess <- function(object) {
  errors = character()
  if (!is_scalar(object@alpha) & !is.GammaRV(object@alpha)) {
    errors <- c(errors, "alpha must be a constant scalar or gamma random variable")
  }
  if (is_numeric_scalar(object@alpha)) {
    if (object@alpha <= 0) {
      errors <- c(errors, "alpha must be > 0")
    }
  }
  if (!is.null(object@H)) {
    errors <- c(errors, "H must be NULL")
  }
  if (length(errors)==0) TRUE else errors
}

setClass("DirichletProcess", 
         slots = c(alpha="ANY", H="ANY"), validity=.check_DirichletProcess)

#' Dirichlet Process
#' 
#' @description
#' Represents a Dirichlet process.
#' 
#' @param alpha Concentration parameter. Must be a strictly positive real 
#' number or a [`GammaRV`] hyperprior.
#' @param H Base distribution. Currently only NULL values are supported.
#' @return A `DirichletProcess` object.
#' @rdname DirichletProcess
#' 
#' @export
DirichletProcess <- function(alpha, H=NULL) {
  new("DirichletProcess", alpha=alpha, H=H)
}

#' @param x An \R object
#' @return \code{is.DirichletProcess} returns TRUE if the argument is a 
#' `DirichletProcess` object and FALSE otherwise.
#' 
#' @noRd
#' @keywords internal
is.DirichletProcess <- function(x) inherits(x, "DirichletProcess")