#' @include ExchangERModel.R utils.R
NULL

#' Inference Result
#' 
#' @description 
#' Represents the result of running inference on an [`ExchangERModel-class`] 
#' object.
#' 
#' @slot history a list containing parameters of interest and summary 
#'   statistics along the Markov chain. See below for details.
#' @slot state an [`ExchangERModel-class`] object which represents the state 
#'   of the model at the end of the Markov chain.
#' 
#' @details 
#' The `history` may include the following objects:
#' \item{clust_params}{a [`coda::mcmc`] matrix object recording samples of 
#' the clustering prior hyperparameters (if they are random). Rows index 
#' samples along the Markov chain and columns index parameters.}
#' \item{links}{a [`coda::mcmc`] matrix object recording samples of the 
#' linkage/clustering structure. Rows index samples along the Markov chain 
#' and columns index records. Each sample of the linkage/clustering structure 
#' is encoded as an integer membership vector: records with the same integer 
#' value are assigned to the same entity.}
#' \item{distort_probs}{a [`coda::mcmc`] matrix object recording 
#' the distortion probabilities for each attribute/file. Rows index samples 
#' along the Markov chain and columns index attributes/files.}
#' \item{n_linked_ents}{a [`coda::mcmc`] vector object recording the 
#' total number of entities (clusters) that are linked to at least one 
#' record.}
#' 
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


#' Extract Samples
#' 
#' @description 
#' Extract samples from a fitted model
#' 
#' @param x an [`ExchangERFit-class`] object
#' @param params an optional character vector specifying the names of
#'   parameters/summary statistics to extract. If not specified, all 
#'   parameters and summary statistics are extracted.
#' @param include a logical scalar indicating whether the parameters in 
#'   `params` should be included (TRUE) or excluded (FALSE).
#'  
#' @return
#' If params is of length 1, the samples for the requested parameter are 
#' returned as a [`coda::mcmc`] object.
#' 
#' If params is NULL or of length greater than 1, the parameters are returned 
#' in a named list. The samples for each parameter are represented as a
#' [`coda::mcmc`] object.
#' 
#' @examples
#' #' ## Initialize a model for ER of RLdata500
#' library(comparator)
#' distortionPrior <- BetaRV(1,5)
#' clust_prior <- PitmanYorRP(alpha = GammaRV(1.0, 1.0), d = BetaRV(1.0, 1.0))
#' attr_params <- list(
#'   "fname_c1" = Attribute(Levenshtein(), distortionPrior),
#'   "lname_c1" = Attribute(Levenshtein(), distortionPrior),
#'   "bd" = CategoricalAttribute(distortionPrior),
#'   "by" = CategoricalAttribute(distortionPrior),
#'   "bm" = CategoricalAttribute(distortionPrior)
#' )
#' model <- exchanger(RLdata500, attr_params, clust_prior)
#' 
#' ## Run inference
#' fit <- run_inference(model, n_samples=2000, burnin_interval=1000)
#' 
#' ## Extract samples of the linkage structure
#' links <- extract(fit, "links")
#' 
#' @export
setGeneric("extract", function(x, params = NULL, include = TRUE) standardGeneric("extract"))


#' @export
setMethod("extract", signature = c(x = "ExchangERFit"), 
          function (x, params = NULL, include = TRUE) {
            if (is.null(params)) return(x@history)
            avail_params <- names(x@history)
            sel_params <- if (include) {
              intersect(params, avail_params)
            } else {
              setdiff(avail_params, params)
            }
            if (length(params) == 1 & length(sel_params) == 1) {
              x@history[[sel_params]]
            } else {
              x@history[sel_params]
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

#' Run inference for a model
#' 
#' @description
#' Runs approximate inference for a model using Markov chain Monte 
#' Carlo (MCMC).
#' 
#' @param x an \R object representing a model or a model with sampling 
#'   history. Currently, there are methods defined for [`ExchangERModel-class`] 
#'   objects and [`ExchangERFit-class`] objects. Passing a model begins 
#'   sampling from scratch, in which case a burn-in period is recommended. 
#'   Passing a model with sampling history will resume inference, adding 
#'   samples to the existing history.
#' @param n_samples a positive integer. The number of samples to generate 
#'   after applying burn-in and thinning.
#' @param thin_interval a positive integer. The period for saving samples 
#'   (intermediate samples are discared). The default value is 1, which 
#'   means no thinning is applied.
#' @param burnin_interval a non-negative integer. The number of initial samples 
#'   to discard as burn-in. The default value is 0, which means no burn-in is 
#'   applied.
#' @param ... further arguments passed to or from other methods. 
#' @return 
#' Returns an [`ExchangERFit-class`] object with the following slots:
#' \item{history}{a list containing the sampling history for 
#'   parameters/summary statistics along the Markov chain.}
#' \item{state}{a model object of the same class as the input, which 
#'   represents the state of all parameters in the last step of the Markov 
#'   chain.}
#' 
#' @rdname run_inference
#' @export
setGeneric("run_inference", 
           function(x, n_samples, thin_interval=1, burnin_interval=0, ...) {
             standardGeneric("run_inference")
           })

#' @describeIn run_inference Specialization for [`ExchangERModel-class`]
#' @importFrom coda mcmc
#' @export
setMethod("run_inference", signature(x = "ExchangERModel"), 
  function(x, n_samples, thin_interval=1, burnin_interval=0, ...) {
    # Check validity of arguments
    if (!is_scalar(n_samples)) stop("n_samples must be a scalar")
    if (!is_scalar(thin_interval)) stop("thin_interval must be a scalar")
    if (!is_scalar(burnin_interval)) stop("burnin_interval must be a scalar")
    if (n_samples < 0) stop("n_samples must be positive")
    if (thin_interval < 1) stop("thin_interval must be positive")
    if (burnin_interval < 0) stop("burnin_interval must be non-negative")
    if (as.integer(n_samples) != n_samples) 
      stop("n_samples must be an integer")
    if (as.integer(thin_interval) != thin_interval) 
      stop("thin_interval must be an integer")
    if (as.integer(burnin_interval) != burnin_interval) 
      stop("burnin_interval must be an integer")
    
    t1 <- Sys.time()
    if (x@iteration == 0) {
      # Not resuming a Markov chain
      start_iter <- ifelse(burnin_interval, burnin_interval, thin_interval)
    } else {
      # Resuming a Markov chain
      start_iter <- x@iteration + thin_interval
    }
    
    # Call C++ function
    result <- .sample(x, n_samples, thin_interval, burnin_interval)
    
    # These attributes are needed for the resuming logic
    attr(result@history, 'mcpar') <- c(start_iter, result@state@iteration, thin_interval)
    
    # Convert history to coda::mcmc objects
    for (i in seq_along(result@history)) {
      result@history[[i]] <- mcmc(result@history[[i]], start=start_iter, 
                                  thin=thin_interval)
    }
    t2 <- Sys.time()
    message("Completed sampling in ", format(t2-t1))
  
    return(result)
  })

#' @describeIn run_inference Specialization for [`ExchangERFit-class`]
#' @export
setMethod("run_inference", signature(x = "ExchangERFit"), 
  function(x, n_samples, thin_interval=1, burnin_interval=0, ...) {
    # Resuming a previous chain
    if (burnin_interval != 0) 
      warning("ignoring burnin_interval=", burnin_interval, " when resuming")
    old_thin_interval <- attr(x@history, 'mcpar')[3]
    if (thin_interval != old_thin_interval) 
      warning("using thin_interval=", old_thin_interval, " from previous run")
    result <- run_inference(x@state, n_samples, thin_interval, burnin_interval)
    result <- combine_results(x, result)
    return(result)
  })
