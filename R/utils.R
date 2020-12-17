is_scalar <- function(x) is.atomic(x) && length(x) == 1L

is_numeric_scalar <- function(x) is_scalar(x) && is.numeric(x)

hmean <- function(x) 1.0 / mean(1.0 / x)