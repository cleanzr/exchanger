#' @include EXERModel.R utils.R
NULL

#' Class for MCMC output
#' 
#' @slot history a list containing parameters of interest and summary 
#'   statistics along the Markov chain.
#' @slot state an [`EXERModel-class`] object which represents the state 
#'   of the model at the end of the Markov chain.
setClass("MCMCResult", 
         slots = c(history = "list",
                   state = "EXERModel"))

#' Function to concatenate two MCMCResult objects
#' 
#' @param resultA,resultB [`MCMCResult-class`] objects. `resultB` must occur 
#'   after `resultA` along the Markov chain. Additionally, the results must use 
#'   the same `thin_interval`.
#' @return an [`MCMCResult-class`] object
#' 
#' @importFrom coda mcmc
#' @noRd
combine_results <- function(resultA, resultB) {
  # Assume that both inputs are instances of type 'MCMCResult' and that 
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
  
  return(new("MCMCResult", history=hist_combined, state=resultB@state))
}

#' Run inference for a model
#' 
#' @description
#' Runs approximate inference for a model using Markov chain Monte 
#' Carlo (MCMC).
#' 
#' @param x an \R object representing a model or a model with sampling 
#'   history. Currently, there are methods defined for [`EXERModel-class`] 
#'   objects and [`MCMCResult-class`] objects. Passing a model begins 
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
#' Returns an [`MCMCResult-class`] object with the following slots:
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

#' @describeIn run_inference Specialization for [`EXERModel-class`]
#' @importFrom coda mcmc
#' @export
setMethod("run_inference", signature(x = "EXERModel"), 
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

#' @describeIn run_inference Specialization for [`MCMCResult-class`]
#' @export
setMethod("run_inference", signature(x = "MCMCResult"), 
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
