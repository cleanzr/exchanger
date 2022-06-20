#' @include utils.R
NULL

#' Sample Ewens-Pitman partition
#' 
#' @description This function uses a Chinese Restaurant Process construction to 
#'   generate a sample.
#' 
#' @param n Number of elements to partition.
#' @param alpha,sigma Ewens-Pitman parameters as described in Pitman (2006).
#' 
#' @return A partition of the integers {1, ..., n} represented as a list of 
#'     blocks. 
sample_ewens_pitman <- function(n, alpha, sigma) {
  if (sigma >= 0 && sigma <= 1) {
    if (alpha <= -sigma) {
      stop("`alpha` must be strictly greater than `sigma` when `sigma` is in the interval [0, 1]")
    }
  } else if (sigma < 0) {
    if (alpha < 0) {
      stop ("`alpha` must positive")
    }
    if (!is_wholenumber(alpha / sigma)) {
      stop("`alpha` must be an integer multiple of negative `sigma`")
    }
  }
  
  clusters <- list()
  for (elem_id in seq_len(n)) {
    # Weights associated with existing entities
    weights <- as.numeric(sapply(clusters, function(members) abs(length(members) - sigma)))
    # Add weight associated with new entity
    weights <- c(weights, abs(alpha + length(clusters) * sigma))
    ent_id <- sample.int(length(weights), size = 1, replace = TRUE, prob = weights)
    if (ent_id > length(clusters)) {
      clusters[[ent_id]] <- elem_id
    } else {
      members <- clusters[[ent_id]]
      clusters[[ent_id]] <- c(members, elem_id)
    }
  }
  
  return(clusters)
}