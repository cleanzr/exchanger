#' @include BetaRV.R utils.R DirichletProcess.R
NULL

setClass("Attribute", slots=c(dist_fn="function",
                              distort_prob_prior="BetaRV",
                              distort_dist_prior="DirichletProcess",
                              exclude_entity_value="logical",
                              entity_dist_prior="ANY"), 
         validity = function(object) {
           errors <- character()
           if (length(object@exclude_entity_value) != 1) {
             errors <- c(errors, "`exclude_entity_value` must be a logical vector of length 1")
           }
           if (!(is.DirichletRV(object@entity_dist_prior) || is.null(object@entity_dist_prior) || is.numeric(object@entity_dist_prior))) {
             errors <- c(errors, "`entity_dist_prior` must be a DirichletRV, numeric vector or NULL")
           }
           if (is.numeric(object@entity_dist_prior)) {
             if (any(object@entity_dist_prior <= 0) || sum(object@entity_dist_prior) != 1.0) {
               errors <- c(errors, "numeric vector passed to `entity_dist_prior` must be a normalized pmf")
             }
           }
           if (length(errors)==0) TRUE else errors
         })

#' Attribute specification
#' 
#' @description
#' Constructs an `Attribute` object, which stores model parameters associated 
#' with an attribute (e.g. 'name of company', 'zip code', 'date of birth' etc.).
#' 
#' @details 
#' The model parameters specified by this function are: 
#' * The prior on the _entity attribute distribution_, denoted 
#'   `entity_dist_prior` in the parameter list. The entity attribute distribution 
#'   is an important component of the entity model. When instantiating a new 
#'   entity, this distribution selects a value for the corresponding 
#'   attribute---e.g. "July 10 1973" for a "date of birth" attribute.
#' * The prior on the _distortion probability_ for the attribute, denoted 
#'   `distort_prob_prior` in the parameter list. This controls the likelihood 
#'   that values of the attribute are distorted when they are recorded in the 
#'   data.
#' * The prior on the _distortion distribution_ for the attribute, denoted 
#'   `distort_dist_prior` in the parameter list. The distortion distribution 
#'   is responsible for generating a distorted attribute value \eqn{x} 
#'   (recorded in the data) given the latent entity attribute value \eqn{y}. 
#'   It is parameterized by a distance function (see next bullet point).
#' * The _attribute distance function_, denoted `dist_fn` in the parameter 
#'   list. It takes arguments \eqn{y} and \eqn{x} and returns a distance, 
#'   which quantifies the likelihood that \eqn{x} is a distorted alternative 
#'   to \eqn{y}. The function need not satisfy the properties of a distance 
#'   metric.
#' 
#' @param dist_fn A vectorized distance function that takes two vectors of 
#'   attribute values \eqn{y} and \eqn{x} as arguments, and returns a vector 
#'   containing the elementwise distances between the values. The function 
#'   is not assumed to be symmetric in it's arguments. The first argument 
#'   \eqn{y} corresponds to a (latent) entity attribute value and the second
#'   argument \eqn{x} corresponds to an (observed) record attribute value. 
#'   This function parameterizes the distortion distribution (see below). 
#'   The efficiency/scalability of inference can be improved by returning a 
#'   distance of `Inf` whenever \eqn{x} is unlikely to be a distortion of 
#'   \eqn{y}. This can be achieved by wrapping a distance function in 
#'   [`transform_dist_fn`].
#' @param distort_prob_prior A [`BetaRV`] object. Specifies the prior 
#'   on the distortion probability for the attribute. 
#' @param distort_dist_prior A [`DirichletProcess`] object. Specifies 
#'   the prior on the distortion distribution for a record value \eqn{x} 
#'   with entity value \eqn{y}. A Dirichlet Process with a small 
#'   concentration parameter yields a more concentrated distribution, which 
#'   may be useful for large attribute domains. If the concentration parameter 
#'   is set to `Inf`, the Dirichlet Process collapses to the base distribution 
#'   (i.e. the prior is a point mass on the base distribution). 
#'   The base distribution is set to \eqn{p(x|y) \propto \exp(-\mathrm{dist_fn}(y,x))}{p(x|y) ∝ exp(-dist_fn(y,x))}. 
#'   Defaults to `DirichletProcess(Inf)`.
#' @param exclude_entity_value A boolean specifying whether the entity 
#'   value \eqn{y} should be excluded from the support of the distortion 
#'   distribution for record value \eqn{x}. If TRUE, the distance returned 
#'   by \eqn{\mathrm{dist_fn}(y,y)}{dist_fn(y,y)} is replaced by `Inf` for all 
#'   values \eqn{y}. Defaults to TRUE.
#' @param entity_dist_prior A [`DirichletRV`] object, numeric vector or 
#'   NULL. This parameter specifies the prior on the entity attribute 
#'   distribution. 
#'   
#'   If a [`DirichletRV`] object is passed, a Dirichlet prior is placed on the 
#'   entity attribute distribution.
#'   
#'   An arbitrary fixed distribution can be specified by passing a numeric 
#'   vector. It must have the same length as the domain of the attribute and 
#'   be normalized to 1.
#'   
#'   If NULL, the entity attribute distribution is fixed to the empirical 
#'   distribution of the attribute values, as computed from the observed data.
#'   
#'   
#' @return The `Attribute` constructor returns an `Attribute` object.
#' 
#' @seealso 
#' A named list of `Attribute` objects is required when initializing the model 
#' (see [`exchanger`]).
#' 
#' @export
Attribute <- function(dist_fn, distort_prob_prior, 
                      distort_dist_prior=DirichletProcess(Inf), 
                      exclude_entity_value=TRUE, entity_dist_prior=NULL) {
  new("Attribute", dist_fn=dist_fn, distort_prob_prior=distort_prob_prior,
      distort_dist_prior=distort_dist_prior, 
      exclude_entity_value=exclude_entity_value,
      entity_dist_prior=entity_dist_prior)
}

setClass("CategoricalAttribute", contains = "Attribute")

#' @rdname Attribute
#' @return The `CategoricalAttribute` constructor returns a 
#' `CategoricalAttribute` object. 
#' It is intended for modeling categorical attributes, and uses a 
#' a constant distance function.
#' 
#' @importFrom comparator Constant
#' @export
CategoricalAttribute <- function(distort_prob_prior, 
                                 distort_dist_prior=DirichletProcess(Inf), 
                                 exclude_entity_value=TRUE,
                                 entity_dist_prior=NULL) {
  new("CategoricalAttribute", dist_fn=Constant(), 
      distort_prob_prior=distort_prob_prior, 
      distort_dist_prior=distort_dist_prior, 
      exclude_entity_value=exclude_entity_value,
      entity_dist_prior=entity_dist_prior)
}