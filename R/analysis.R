#' @include utils.R
NULL

#' Most Probable Clusters
#'
#' @description
#' Given samples of clusterings, this function returns a list containing the
#' _most probable cluster_ for each record.
#'
#' @note
#' The collection of most probable clusters is not a valid clustering, as 
#' records may appear in multiple clusters. To obtain a clustering that 
#' one can apply the [`smp_clusters`] function to the output of 
#' this function. 
#' 
#' In Steorts et al. (2016), the most probable clusters are referred to as 
#' _most probable maximal matching sets_.
#'
#' @param x an \R object containing samples of clusterings. Currently, methods 
#'   are defined for [`MCMCResult-class`] objects, [`coda::mcmc`] objects 
#'   and [`base::matrix`] objects (where rows index samples and columns index 
#'   cluster assignments for each record).
#' @param ... further arguments passed to or from other methods.
#' 
#' @return 
#' Returns a `MPClusters` object, which is a list containing the following
#' entries:
#' \item{clusters}{a list of vectors. The i-th entry specifies the most 
#'   frequently occurring cluster involving the i-th record among the samples 
#'   of clusterings. The cluster is represented as a vector of record 
#'   identifiers.}
#' \item{probabilities}{a vector of probabilities. The i-th entry specifies
#'   the probability (frequency of occurrence) of the most frequently 
#'   occurring cluster involving the i-th record.}
#'
#' @references Steorts, R. C., Hall, R. & Fienberg, S. E. A Bayesian Approach
#' to Graphical Record Linkage and Deduplication. _JASA_ **111**,
#' 1660–1672 (2016).
#' 
#' @rdname mp_clusters
#' @importFrom coda mcmc
#' @export
mp_clusters <- function(x, ...) {
  UseMethod("mp_clusters")
}

#' @rdname mp_clusters
#' @param rec_ids an optional vector of record identifiers. If NULL, the
#'   record identifiers are extracted from `x` if available.  Otherwise 
#'   integers are used to identify the records.
#'   
#' @export
mp_clusters.matrix <- function(x, rec_ids = NULL, ...)
{
  if (any(is.na(x))) stop("x must not contain NAs")
  
  if (!is.null(rec_ids)) {
    # Check consistent dimensions
    if (length(rec_ids) != ncol(x))
      stop("length of rec_ids must match number of columns in `x`")
  } else {
    # Try to get rec_ids from colnames
    if (!is.null(colnames(x))) {
      rec_ids <- colnames(x)
    }
  }
  
  # Represent cluster identifiers as integers (starting at 0), for
  # compatibility with the C++ function
  if (!is.integer(x))
  {
    original_dim <- dim(x)
    x <- unclass(factor(x))
    dim(x) <- original_dim
  }
  
  mpc <- .mp_clusters(x)
  
  # Replace integer rec ids with provided ones
  if (!is.null(rec_ids)) {
    mpc$clusters <- lapply(mpc$clusters, function(cid) rec_ids[cid])
    names(mpc$clusters) <- as.character(rec_ids)
    names(mpc$probabilities) <- as.character(rec_ids)
  }
  
  attr(mpc, "rec_ids") <- rec_ids
  class(mpc) <- c("MPClusters", class(mpc))
  return(mpc)
}

#' @rdname mp_clusters
#' @export
mp_clusters.MCMCResult <- function(x, ...) {
  samples <- x@history$links
  if (is.null(samples))
    stop("MCMCResult does not contain history of the linkage structure along the chain")
  rec_ids <- x@state@rec_ids
  mp_clusters(samples, rec_ids, ...)
}

#' @rdname mp_clusters
#' @export
mp_clusters.mcmc <- function(x, rec_ids=NULL, ...) {
  mp_clusters(unclass(x), rec_ids, ...)
}

#' Shared Most Probable Clusters
#'
#' @description
#' Computes a point estimate of the most likely clustering that obeys
#' transitivity constraints based on samples of clusterings. The method was
#' proposed by Steorts et al. (2016), where it is referred to as the
#' method of _shared most probable maximal matching sets_.
#'
#' @param x an \R object containing samples of clusterings or a summary 
#'   of the most probable clusters. Currently, methods are defined for 
#'   `MPClusters` objects (as output by [`mp_clusters`]), 
#'   [`MCMCResult-class`] objects, [`coda::mcmc`] objects 
#'   and [`base::matrix`] objects (where rows index samples and columns index 
#'   cluster assignments for each record).
#' @param ... further arguments passed to or from other methods.
#' 
#' @return Returns the shared most probable clustering, represented as a list
#'   of vectors. Each vector contains the identifiers of records belonging to
#'   a cluster.
#'
#' @references Steorts, R. C., Hall, R. & Fienberg, S. E. A Bayesian Approach
#' to Graphical Record Linkage and Deduplication. _JASA_ **111**,
#' 1660–1672 (2016).
#'
#' @seealso
#' The [`clevr::clusters_to_pairs`] and [`clevr::clusters_to_membership`] 
#' functions are useful for converting the output of this function to different
#' clustering representations.
#'
#' @rdname smp_clusters
#' @export
smp_clusters <- function(x, ...)
  UseMethod("smp_clusters")


#' @rdname smp_clusters
#' @importFrom utils relist
#' @export
smp_clusters.MPClusters <- function(x, ...) {
  # Special case
  if (length(x$clusters) == 0) {
    return(list())
  }
  
  # Express clusters in terms of integer record identifiers (starting at 1)
  # for compatibility with C++
  unlisted_clusters <- unlist(x$clusters)
  rec_ids <- attr(x, "rec_ids")
  unlisted_clusters <- unclass(factor(unlisted_clusters, rec_ids))
  clusters <- relist(unlisted_clusters, x$clusters)
  
  smpc <- .smp_clusters(clusters)
  
  # Express in terms of original rec_ids
  unlisted_smpc <- unlist(smpc)
  unlisted_smpc <- rec_ids[unlisted_smpc]
  smpc <- relist(unlisted_smpc, smpc)
  
  return(smpc)
}

#' @rdname smp_clusters
#' @param rec_ids an optional vector of record identifiers. This argument is 
#'   ignored when an `MPClusters` object is passed for `x`. If NULL, the
#'   record identifiers are extracted from `x` if available.  Otherwise 
#'   integers are used to identify the records.
#'
#' @export
smp_clusters.matrix <- function(x, rec_ids=NULL, ...) {
  mpc <- mp_clusters(x, rec_ids = rec_ids, ...)
  smp_clusters.MPClusters(mpc, ...)
}

#' @rdname smp_clusters
#' @export
smp_clusters.MCMCResult <- function(x, ...)
{
  mpc <- mp_clusters(x, ...)
  smp_clusters(mpc, ...)
}

#' @rdname smp_clusters
#' @export
smp_clusters.mcmc <- function(x, rec_ids=NULL, ...) {
  smp_clusters(unclass(x), rec_ids, ...)
}

# #' Posterior match probabilities
# #' 
# #' @description 
# #' Computes posterior match probabilities for a given subset of record pairs 
# #' based on MCMC samples.
# #' 
# #' @param result a [`MCMCResult-class`] object
# #' @param rec_ids.x,rec_ids.y vectors of identifiers which specify the record 
# #'   pairs of interest. The vectors must be the same length, and they 
# #'   conform to the record identifier format used in `result`.
# #' @return a vector of pairwise match probabilities of the same length as 
# #'   `rec_ids.x`.
# pairwiseMatchProbability <- function(result, rec_ids.x, rec_ids.y) {
#   if (!inherits(result, "MCMCResult")) 
#     stop("result must be a 'MCMCResult' object")
#   if (length(rec_ids.x) != length(rec_ids.y))
#     stop("rec_ids.x and rec_ids.y must be vectors of the same length")
#   if (any(rec_ids.x == rec_ids.y))
#     stop("corresponding entries in rec_ids.x and rec_ids.y must not be equal")
#   
#   rec_ids.x <- as.character(rec_ids.x)
#   rec_ids.y <- as.character(rec_ids.y)
#   
#   linksHist <- result@history$links
#   numPairs <- length(rec_ids.x)
#   out <- numeric(length=numPairs)
#   for (i in seq_len(numPairs)) {
#     out[i] = mean(linksHist[,rec_ids.x[i]] == linksHist[,rec_ids.y[i]])
#   }
#   
#   return(out)
# }
