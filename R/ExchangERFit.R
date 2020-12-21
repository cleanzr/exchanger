#' @include ExchangERModel.R utils.R
NULL

#' Fitted Model
#' 
#' @description 
#' Stores data associated with a fitted [`ExchangERModel-class`] object, 
#' including posterior samples, diagnostics, etc.
#' 
#' @slot history a list containing model parameters of interest and summary 
#'   statistics along the Markov chain. See below for details.
#' @slot state an [`ExchangERModel-class`] object which represents the state 
#'   of the model at the end of the Markov chain.
#' 
#' @details 
#' The `history` may include the following objects:
#' \describe{
#'   \item{clust_params}{a [`coda::mcmc`] matrix object recording samples of 
#'     the clustering prior hyperparameters (if they are random). Rows index 
#'     samples along the Markov chain and columns index parameters.}
#'   \item{links}{a [`coda::mcmc`] matrix object recording samples of the 
#'     linkage/clustering structure. Rows index samples along the Markov chain 
#'     and columns index records. Each sample of the linkage/clustering 
#'     structure is encoded as an integer membership vector: records with the 
#'     same integer value are assigned to the same entity.}
#'   \item{distort_probs}{a [`coda::mcmc`] matrix object recording 
#'     the distortion probabilities for each attribute/file. Rows index samples 
#'     along the Markov chain and columns index attributes/files.}
#'   \item{n_linked_ents}{a [`coda::mcmc`] vector object recording the 
#'     total number of entities (clusters) that are linked to at least one 
#'     record.}
#' }
#' They can be accessed using the [`extract`] method.
setClass("ExchangERFit", 
         slots = c(history = "list",
                   state = "ExchangERModel"))

#' @importFrom mcmcse multiESS ess
setMethod("show", "ExchangERFit", function(object) {
  # Extract info about burn-in, thinning, number of samples
  mcpar <- attr(object@history, "mcpar")
  start_iter <- mcpar[1]
  end_iter <- mcpar[2]
  thin <- mcpar[3]
  n_samples <- end_iter - start_iter + thin
  
  # Compute effective sample size for particular variables, if present
  valid_ess_varnames <- c("n_linked_ents", "distort_probs", "clust_params")
  ess <- list()
  for (varname in intersect(valid_ess_varnames, names(object@history))) {
    var <- object@history[[varname]]
    if (ncol(var) > 1)
      ess[[varname]] <- mcmcse::multiESS(var)
    else
      ess[[varname]] <- mcmcse::ess(var)
  }
  
  cat("Fitted ExchangERModel\n",
      n_samples, if (n_samples > 1) " samples" else " sample", 
      " after thinning (interval=", thin, 
      ") with burn-in (interval=", start_iter, ")\n",
      "Recorded parameters/summary statistics: \n", 
      "  ",  toString(names(object@history)), "\n", sep="")
  if (length(ess) > 0) {
    cat("Effective sample size:\n")
    for (varname in names(ess)) {
      cat("  ", varname, ": ", ess[[varname]], "\n", sep="")
    }
  }
})

#' Function to concatenate two ExchangERFit objects
#' 
#' @param resultA,resultB [`ExchangERFit-class`] objects. `resultB` must occur 
#'   after `resultA` along the Markov chain. Additionally, the results must use 
#'   the same `thin_interval`.
#' @return an [`ExchangERFit-class`] object
#' 
#' @importFrom coda mcmc
#' @noRd
combine_results <- function(resultA, resultB) {
  # Assume that both inputs are instances of type 'ExchangERFit' and that 
  # they refer to the same data set & model
  
  mcparA <- attr(resultA@history, 'mcpar')
  mcparB <- attr(resultB@history, 'mcpar')
  
  # Check validity of inputs
  if (diff(mcparA[1:2]) == 0) stop("resultA has no history")
  if (diff(mcparB[1:2]) == 0) stop("resultB has no history")
  if (mcparA[2] >= mcparB[1] ) {
    stop("resultB must occur chronologically after resultA")
  }
  if (mcparA[3] != mcparB[3]) 
    stop("results must have the same thin_interval")
  
  if ( !setequal(names(resultA@history), names(resultB@history)) ) {
    stop("results must have the same history variables")
  }
  
  # Combine history matrices/vectors
  hist_combined <- vector(mode="list", length = length(resultA@history))
  attr(hist_combined, 'mcpar') <- c(mcparA[1], resultB@state@iteration, mcparA[3])
  names(hist_combined) <- names(resultA@history)
  for (varname in names(resultA@history)) {
    varA <- resultA@history[[varname]]
    varB <- resultB@history[[varname]]
    if (is.matrix(varA)) {
      history <- rbind(varA, varB)
    } else {
      history <- c(varA, varB)
    }
    hist_combined[[varname]] <- mcmc(history, start=mcparA[1], thin = mcparA[3])
  }
  
  return(new("ExchangERFit", history=hist_combined, state=resultB@state))
}