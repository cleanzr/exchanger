is_scalar <- function(x) is.atomic(x) && length(x) == 1L

is_numeric_scalar <- function(x) is_scalar(x) && is.numeric(x)

hmean <- function(x) 1.0 / mean(1.0 / x)

approx_equal <- function(x, y, tol = .Machine$double.eps^0.5) abs(x - y) < tol

is_wholenumber <- function(x, tol = .Machine$double.eps^0.5)  approx_equal(x, round(x), tol = tol)
