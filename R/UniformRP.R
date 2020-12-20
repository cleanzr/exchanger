#' @include RP.R
NULL

setClass("UniformRP", contains = "RP")

#' Uniform Partition
#' 
#' @description
#' Specifies a uniform prior for the clustering---i.e. every possible 
#' partitioning of the records into clusters is equally likely.
#' 
#' @return A `UniformRP` object which is a sub-class of 
#' [`RP-class`].
#' 
#' @seealso \code{\link{PitmanYorRP}}, 
#' \code{\link{EwensRP}}, \code{\link{GeneralizedCouponRP}}
#' 
#' @export
UniformRP <- function() {
  new("UniformRP")
}

setMethod("format", "UniformRP", function(x, ...) {
  "UniformRP()"
})

setMethod("show", "UniformRP", function(object) {
  writeLines(format(object))
  invisible(object)
})

#' @param x An \R object.
#' @return `is.UniformRP` returns TRUE if the argument is a 
#'   `UniformRP` object and FALSE otherwise.
#' 
#' @noRd
#' @keywords internal
is.UniformRP <- function(x) inherits(x, "UniformRP")