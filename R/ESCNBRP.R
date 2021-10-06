#' @include RP.R
NULL

.check_ESCNB <- function(object) {
  errors = character()
  if (!is_scalar(object@p) && !is.BetaRV(object@p)) {
    errors <- c(errors, "p must be a scalar or a BetaRV")
  }
  if (!is_scalar(object@r) & !is.GammaRV(object@r)) {
    errors <- c(errors, "r must be a scalar or a GammaRV")
  }
  if (is_numeric_scalar(object@p)) {
    if (object@p < 0 || object@p > 1)
      errors <- c(errors, "p must be on the interval [0, 1]")
  }
  if (is_numeric_scalar(object@r)) {
    if (object@r <= 0) {
      errors <- c(errors, "r must be strictly positive")
    }  
  }
  if (length(errors)==0) TRUE else errors
}

setClass("ESCNBRP", 
         slots = c(p="ANY", r="ANY"), 
         validity=.check_ESCNB, contains = "RP")

#' ESC-NB Random Partition
#' 
#' @description
#' Represents an ESC-NB (Exchangeable Sequences of Clusters-Negative Binomial) 
#' random partition.
#' 
#' @details 
#' TODO
#' 
#' @param p Probability of success. Must be on the unit interval \eqn{(0, 1]}
#'   or a [`BetaRV`] object.
#' @param r Target number of successes before negative binomial experiment 
#'   is stopped. Must be strictly a positive real or a [`GammaRV`] object.
#' @return The constructor `ESCNBRP` returns a `ESCNBRP` object which is a 
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
ESCNBRP <- function(p, r) {
  new("ESCNBRP", p=p, r=r)
}

setMethod("format", "ESCNBRP", function(x, ...) {
  p <- format(x@p, ...)
  r <- format(x@r, ...)
  paste0("ESCNBRP(p=", p, ", r=", r, ")")
})

setMethod("show", "ESCNBRP", function(object) {
  writeLines(format(object))
  invisible(object)
})

#' @param x An \R object.
#' @return `is.ESCNBRP` returns TRUE if the argument is a `ESCNBRP` object and 
#'   FALSE otherwise.
#' 
#' @noRd
#' @keywords internal
is.ESCNBRP <- function(x) inherits(x, "ESCNBRP")