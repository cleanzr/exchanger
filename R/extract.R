#' @include ExchangERFit.R
NULL

setGeneric("extract", function(x, params = NULL, include = TRUE) standardGeneric("extract"))

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
#' ## Initialize a model for entity resolution of RLdata500
#' library(comparator)
#' distort_prior <- BetaRV(1,5)
#' clust_prior <- PitmanYorRP(alpha = GammaRV(1.0, 1.0), d = BetaRV(1.0, 1.0))
#' attr_params <- list(
#'   fname_c1 = CategoricalAttribute(distort_prior),
#'   lname_c1 = CategoricalAttribute(distort_prior),
#'   bd = CategoricalAttribute(distort_prior),
#'   by = CategoricalAttribute(distort_prior),
#'   bm = CategoricalAttribute(distort_prior)
#' )
#' model <- exchanger(RLdata500, attr_params, clust_prior)
#' 
#' ## Run inference
#' fit <- run_inference(model, n_samples=100, burnin_interval=100)
#' 
#' ## Extract samples of the linkage structure
#' links <- extract(fit, "links")
#' 
#' @aliases extract
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