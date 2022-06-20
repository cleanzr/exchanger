#' @include RV.R
NULL

#' Random Partition
#' 
#' @description
#' A virtual class for a random partition, which is a subclass of 
#' [`RV-class`].
#' 
#' @details
#' Three types of exchangeable random partitions are currently supported.
#' 
#' * _Pitman-Yor partition._ This yields a random partition where the  
#'   number of blocks scales sub-linearly (as a power law) in the 
#'   number of items. [`PitmanYorRP`] specifies a random partition of this 
#'   type.
#' * _Ewens partition._ This yields a random partition where the number 
#'   of blocks scales logarithmically in the number of items. [`EwensRP`] 
#'   specifies a random partition of this type.
#' * _Generalized coupon partition._ This yields a random partition where 
#'   the number of blocks approaches a constant asymptotically. 
#'   [`GeneralizedCouponRP`] specifies a random partition of this type.
setClass("RP", contains = "RV")


#' @param x An \R object.
#' @return `is.RP` returns TRUE if the argument is a `RP` object and FALSE 
#' otherwise.
#' 
#' @keywords internal
#' @noRd
is.RP <- function(x) inherits(x, "RP")

#' Generate a Random Partition
#'
#' @param rp a [`RP-class`] object
#' @param n number of integers in partition
#' @param ... other parameters passed on to other methods.
#' 
#' @return 
#' Returns a partition of the integers \eqn{\{1, 2, ..., n\}}{{1, 2, ..., n}} 
#' as a list of blocks.
#' 
#' @examples
#' rp <- GeneralizedCouponRP(1, 1)
#' rrp(rp, 100)
#' 
#' @export
setGeneric("rrp", function(rp, n, ...) standardGeneric("rrp"), 
           signature = c("rp"))