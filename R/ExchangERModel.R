#' @include RP.R AttributeIndex.R Attribute.R
NULL

.check_ExchangERModel <- function(object) {
  errors <- character()
  # Type checks
  if (!all(sapply(object@attr_params, function(x) inherits(x, "Attribute")))) {
    errors <- c(errors, "`attr_params` must be a list of Attribute instances")
  }
  if (!all(sapply(object@attr_indices, function(x) inherits(x, "AttributeIndex")))) {
    errors <- c(errors, "`attr_indices` must be a list of AttributeIndex instances")
  }
  if (!is.integer(object@rec_attrs)) {
    errors <- c(errors, "`rec_attrs` must be an integer matrix")
  }
  if (!is.integer(object@ent_attrs)) {
    errors <- c(errors, "`ent_attrs` must be an integer matrix")
  }
  if (!is.integer(object@rec_distortions)) {
    errors <- c(errors, "`rec_distortions` must be an integer matrix")
  }
  if (!is.double(object@distort_probs)) {
    errors <- c(errors, "`distort_probs` must be a numeric matrix")
  }
  if (!is.atomic(object@rec_ids)) {
    errors <- c(errors, "`rec_ids` must be atomic")
  }
  # Consistency checks
  n_attrs <- length(object@attr_params)
  if (n_attrs == 0) {
    errors < c(errors, "`attr_params` cannot be empty")
  }
  if (n_attrs != nrow(object@ent_attrs)) {
    errors <- c(errors, paste("`ent_attrs` must have",n_attrs,"rows"))
  }
  if (n_attrs != nrow(object@rec_attrs)) {
    errors <- c(errors, paste("`rec_attrs` must have",n_attrs,"rows"))
  }
  if (n_attrs != nrow(object@rec_distortions)) {
    errors <- c(errors, paste("`rec_distortions` must have",n_attrs,"rows"))
  }
  if (n_attrs != ncol(object@distort_probs)) {
    errors <- c(errors, paste("`distort_probs` must have",n_attrs,"columns"))
  }
  n_entities <- length(object@ent_ids)
  if (n_entities == 0) {
    errors <- c(errors, "`ent_ids` must not be empty")
  }
  if (n_entities != ncol(object@ent_attrs)) {
    errors <- c(errors, paste("`ent_attrs` must have",n_entities,"rows"))
  }
  if (length(unique(object@ent_ids)) != n_entities) {
    errors <- c(errors, paste("`ent_ids` must be unique"))
  }
  if (any(is.na(object@ent_ids))) {
    errors <- c(errors, "`ent_ids` must not contain missing values")
  }
  if (length(setdiff(object@links, object@ent_ids)) != 0) {
    errors <- c(errors, "`links` contains invalid entity identifiers")
  }
  n_records <- ncol(object@rec_attrs)
  if (n_records == 0) {
    errors <- c(errors, "`rec_attrs` must contain at least one row")
  }
  if (n_records != length(object@links)) {
    errors <- c(errors, paste("`links` must be of length",n_records))
  }
  if (n_records != ncol(object@rec_distortions)) {
    errors <- c(errors, paste("`rec_distortions` must have",n_records,"rows"))
  }
  if (length(setdiff(object@rec_distortions, c(0, 1))) != 0) {
    errors <- c(errors, "`rec_distortions` must only contain 0 or 1 values")
  }
  if (n_records != length(object@file_ids)) {
    errors <- c(errors, paste("`file_ids` must be of length",n_records))
  }
  file_ids_unique <- unique(object@file_ids)
  n_files <- length(file_ids_unique)
  file_ids_max <- max(file_ids_unique)
  file_ids_min <- min(file_ids_unique)
  if (file_ids_min != 1 | file_ids_max != n_files) {
    errors <- c(errors, "`file_ids` must contain consecutive integer indices starting at 1")
  }
  if (nrow(object@distort_probs) != n_files) {
    errors <- c(errors, paste("`distort_probs` must have",n_files,"rows"))
  }
  if (any(object@distort_probs < 0) | any(object@distort_probs > 1)) {
    errors <- c(errors, "`distort_probs` contains invalid probabilities")
  }
  # TODO: check consistency between val_ids and attr_indices
  if (any(is.na(object@ent_attrs))) {
    errors <- c(errors, "`ent_attrs` cannot have missing values")
  }
  if (class(object@clust_params)[1] != class(object@clust_prior)[1]) {
    errors <- c(errors, "`clust_params` and `clust_prior` are incompatible")
  }
  if (length(errors)==0) TRUE else errors
}

#' ExchangER Model Class
#' 
#' @description
#' This class records the state of all model parameters. It also stores 
#' indexing data structures for each attribute (`attr_indices`) which 
#' are used to accelerate sampling. 
#' 
#' @slot iteration A non-negative integer. Counts the number of successive 
#'   applications of the Markov transition operator. This is not necessarily 
#'   equal to the number of samples if thinning is applied.
#' @slot attr_params A named list of [`Attribute`] objects. Each entry 
#'   in the list specifies the model parameters for an entity attribute, 
#'   and should match one of the row names (attributes) in 
#'   `rec_attrs`.
#' @slot file_ids An integer vector of file identifiers. The i-th entry 
#'   specifies the file/source of the record in the i-th column of 
#'   `rec_attrs`. The file identifiers should be in the set 
#'   \eqn{\{1, \ldots, n_files\}}{{1, ..., n_files}}.
#' @slot rec_ids A character vector of record identifiers. 
#' @slot rec_attrs An integer matrix representation of the observed 
#'   records. Rows correspond to attributes and columns correspond to records. 
#'   The attributes are encoded using integer ids, which index the domain of 
#'   each attribute (starting at 1). 
#' @slot rec_distortions An integer matrix of distortion indicators
#'   of the same dimension as `rec_attrs`. A value of 1 (0) means that 
#'   the corresponding attribute value in `rec_attrs` is distorted 
#'   (not distorted).
#' @slot ent_attrs An integer matrix representation of the 
#'   true entity attributes. Rows correspond to attributes and columns 
#'   correspond to entities. The encoding for the attributes (rows) matches the 
#'   encoding used in `rec_attrs`. The entity identifier of the i-th 
#'   column is recorded in the i-th entry of `ent_ids`. 
#' @slot ent_ids An integer vector of unique entity identifiers. These 
#'   identifiers are used to encode the `links` between records and entities. 
#' @slot links An integer vector of entity identifiers for each record in 
#'   `rec_attrs`. The i-th entry gives the linked entity for the 
#'   i-th record (column) of `rec_attrs`.
#' @slot distort_probs A numeric matrix of distortion probabilities. 
#'   Rows correspond to files/sources and columns correspond to attributes. 
#'   The (i,j)-th entry stores the distortion probability for the i-th file 
#'   id and the j-th attribute.
#' @slot clust_params A [`RP-class`] object where all parameters 
#'   (including random variables) are realized. If all parameters are 
#'   fixed, this slot will be identical to `clust_prior`.
#' @slot attr_indices A named list of [`AttributeIndex`] objects. Each 
#'   entry in the list corresponds to an attribute in `attr_params`.
#' @slot clust_prior A [`RP-class`] object. Represents the prior 
#'   over the clustering of records into entities (a.k.a. linkage structure).
#' 
#' @seealso 
#' The [`exchanger`] function should be used to initialize an ER model 
#' from data.
#' 
#' @export
ExchangERModel <- setClass("ExchangERModel", 
                      slots = c(
                        iteration = "integer",
                        attr_params = "list",
                        file_ids = "integer",
                        rec_ids = "vector",
                        rec_attrs = "matrix",
                        rec_distortions = "matrix",
                        ent_attrs = "matrix",
                        ent_ids = "integer",
                        links = "integer",
                        distort_probs = "matrix",
                        clust_params = "RP",
                        attr_indices = "list",
                        clust_prior = "RP"
                      ),
                      validity = .check_ExchangERModel)

#' @importFrom utils capture.output
setMethod("show", "ExchangERModel", function(object) {
  n_records <- length(object@rec_ids)
  n_files <- nlevels(object@file_ids)
  cat("ExchangERModel\n",
      "Defined on ", n_records, " records from ", n_files, 
      if (n_files > 1) " files" else " file", "\n",
      "Clustering prior:  ", format(object@clust_prior), "\n",
      "Attributes used for matching:\n", sep="")
  attr_params <- object@attr_params
  for (a in names(attr_params)) {
    attr_str <- capture.output(print(attr_params[[a]]))
    cat("  * ", a, ": ", attr_str[1], "\n", sep="")
    writeLines(paste0("  ", attr_str[-1]))
  }
})

#' Initialize an ExchangERModel
#' 
#' @description
#' Initializes an [`ExchangERModel-class`] given observed data and model 
#' hyperparameters.
#' 
#' @param data A data frame or matrix of observed records. Column names are 
#'   are used to refer to particular columns in the other arguments.
#' @param attr_params A named list of [`Attribute`] objects. Each entry 
#'   in the list specifies the model parameters for an entity attribute, 
#'   and should match one of the column names (attributes) in 
#'   `data`.
#' @param clust_prior A [`RP-class`] object which specifies the 
#'   prior over the clustering (a.k.a. linkage structure).
#' @param file_id_colname Column name in `data` that contains file/source 
#'   identifiers for the records. If NULL, the records are assumed to be 
#'   from a single file/source.
#' @param rec_id_colname Column name in `data` that contains unique record 
#'   identifiers. If NULL, use row names to identify records.
#' @return An [`ExchangERModel-class`] object. This can be used as the starting 
#'   point for inference using [`run_inference`] function.
#' 
#' @examples 
#' ## Initialize a model for ER of RLdata500
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
#' @export
exchanger <- function(data, attr_params, clust_prior, 
                      file_id_colname=NULL, rec_id_colname=NULL) {
  # Input validation
  data_colnames <- colnames(data)
  if (is.null(data_colnames))
    stop("data must have colnames")
  n_records <- nrow(data)
  if (n_records < 2)
    stop("data must contain at least two records for ER to be feasible")
  
  if (!is.list(attr_params)) stop("attr_params must be a list")
  attr_names <- names(attr_params)
  if (is.null(attr_names)) stop("attr_params is missing names")
  # Sort attributes so that categorical attributes appear first. The C++ 
  # code for generating candidate matches works faster this way.
  order_attrs <- order(sapply(attr_params, function(x) class(x)), decreasing = TRUE)
  attr_params <- attr_params[order_attrs]
  attr_names <- attr_names[order_attrs]
  n_attrs <- length(attr_params)
  
  missing_attrs <- setdiff(attr_names, data_colnames)
  if (length(missing_attrs) != 0) {
    stop("`data` is missing the following attributes: ", paste(missing_attrs))
  }
  
  if (!is.null(file_id_colname)) {
    if (!is_scalar(file_id_colname) || !is.character(file_id_colname)) 
      stop("`file_id_colname` must be a string")
    if (!is.element(file_id_colname, data_colnames))
      stop("`file_id_colname` does not match any columns in data")
    file_ids <- data[[file_id_colname]]
    if (any(is.na(file_ids))) 
      stop("file identifiers cannot contain missing values")
  } else {
    file_ids <- integer(length = n_records)
  }
  
  if (!is.null(rec_id_colname)) {
    if (!is_scalar(rec_id_colname) || !is.character(rec_id_colname)) 
      stop("rec_id_colname must be a string")
    if (!is.element(rec_id_colname, data_colnames))
      stop("rec_id_colname does not match any columns in data")
    rec_ids <- data[[rec_id_colname]]
    if (length(unique(rec_ids)) != n_records)
      stop("record identifiers must be unique")
  } else {
    if (!is.null(rownames(data))) {
      rec_ids <- rownames(data)
    } else {
      rec_ids <- seq.int(length.out = n_records)  
    }
  }
  
  if (!is.RP(clust_prior)) 
    stop("clust_prior must be a RP")
  
  if (is.GeneralizedCouponRP(clust_prior) && is.numeric(clust_prior@m) && 
      clust_prior@m < n_records)
    stop("Coupon collecting parameter m < n_records is currently not supported.")
  
  if (is.GeneralizedCouponRP(clust_prior) && is.numeric(clust_prior@kappa) && 
      is.infinite(clust_prior@kappa) && !is.numeric(clust_prior@m))
    stop("Random parameter m is currently not supported when kappa = Inf")
  
  # Encode attributes in data as factors. The levels for a factor represent 
  # the domain of that attribute, and are used later on when building the 
  # attribute index.
  data <- data_to_factor(data, attr_params)
  # Extract only the relevant attributes, transform to integer codes and take 
  # the transpose
  # TODO: avoid transposing. Would need to update variables by column, 
  # rather than by row in the Gibbs sampler
  rec_attrs <- t(data.matrix(data[attr_names]))  
  
  # Represent file_ids as integer codes (original values preserved in levels)
  file_ids <- unclass(factor(file_ids))
  files <- levels(file_ids)
  
  # Initialize all random variables
  rec_distortions <- matrix(0L, nrow=n_attrs, ncol=n_records)
  links <- seq_len(n_records)
  distort_probs <- distort_prob_from_prior(attr_params, files)
  clust_params <- clust_params_from_prior(clust_prior, n_records)
  attr_indices <- buildAttributeIndices(data, attr_params)
  
  # Initialize so that each record is linked to a separate entity
  ent_ids <- seq_len(n_records)
  ent_attrs <- rec_attrs
  # Fill missing values with sample drawn from empirical distribution
  for (i in seq_along(attr_indices)) {
    index = attr_indices[[i]]
    if (index@n_missing > 0) {
      ent_attrs[i,is.na(ent_attrs[i,])] <- sample.int(length(index@probs), size=index@n_missing, prob=index@probs, replace=TRUE)
    }
  }
  
  # Fill in names
  rec_ids_char <- as.character(rec_ids)
  names(file_ids) <- rec_ids_char
  colnames(rec_attrs) <- rec_ids_char
  rownames(rec_distortions) <- rownames(rec_attrs)
  colnames(rec_distortions) <- rec_ids_char
  names(links) <- rec_ids_char
  
  ExchangERModel(iteration=0L, attr_params=attr_params, file_ids=file_ids, 
            rec_ids=rec_ids, rec_attrs=rec_attrs, 
            rec_distortions=rec_distortions, 
            ent_attrs=ent_attrs, ent_ids=ent_ids, 
            links=links, distort_probs=distort_probs, 
            clust_params=clust_params, attr_indices=attr_indices, 
            clust_prior=clust_prior)
}


#' Build attribute indices
#' 
#' @param data A data frame or matrix of observed records. Column names are 
#'   are used to refer to particular columns in the other arguments.
#' @param attr_params A named list of [`Attribute`] objects. Each entry 
#'   in the list specifies the model parameters for an entity attribute, 
#'   and should match one of the column names (attributes) in 
#'   `data`.
#' @return A list of attribute indices
#' 
#' @keywords internal
#' @noRd
buildAttributeIndices <- function(data, attr_params) {
  # Preallocate list
  attr_indices <- vector(mode = "list", length=length(attr_params))
  names(attr_indices) <- names(attr_params)
  
  for (i in seq_along(attr_params)) {
    a <- attr_params[[i]]
    attr_name <- names(attr_params)[i]
    message("Building index for \'", attr_name,"\' attribute.")
    if (inherits(a, "CategoricalAttribute")) {
      attr_indices[[i]] <- AttributeIndex(data[[attr_name]], NULL, a@exclude_entity_value)
    } else {
      attr_indices[[i]] <- AttributeIndex(data[[attr_name]], a@dist_fn, a@exclude_entity_value)
    }
  }
  
  return(attr_indices)
}

#' Initialize distortion probabilities based on the prior
#' 
#' @param attr_params A named list of [`Attribute`] objects. Each entry 
#'   in the list specifies the model parameters for an entity attribute, 
#'   and should match one of the column names (attributes) in 
#'   `data`.
#' @param files A character vector containing the unique file identifiers as 
#'   represented in the source `data`.
#' @return A matrix containing an instantiated distortion probability for each 
#'   attribute (indexed by row) and file (indexed by column)
#' 
#' @keywords internal
#' @noRd
distort_prob_from_prior <- function(attr_params, files) {
  n_files <- length(files)
  probs <- sapply(attr_params, function(a) mean(a@distort_prob_prior))
  distort_probs <- t(replicate(probs, n = n_files))
  rownames(distort_probs) <- files
  colnames(distort_probs) <- names(attr_params)
  return(distort_probs)
}


#' Initialize cluster parameters based on the prior
#' 
#' @param clust_prior A [`RP-class`] object which specifies the 
#'   prior over the clustering (a.k.a. linkage structure).
#' @param n_records Number of records in `data`.
#' @return a [`RP-class`] object where any random variables are instantiated
#' 
#' @keywords internal
#' @noRd
clust_params_from_prior <- function(clust_prior, n_records) {
  clust_params <- clust_prior
  
  # Special treatment for the generalized coupon partition
  if (is.GeneralizedCouponRP(clust_params)) {
    if (is.RV(clust_params@m)) clust_params@m <- mean(clust_params@m)
    
    # Prevent infinite values (can happen for some hyperparameter settings)
    if (!is.finite(clust_params@m)) clust_params@m <- n_records
    
    # Prevent number of records exceeding the number of coupons, as it may be 
    # impossible (depending on the data)
    clust_params@m <- max(clust_params@m, n_records)
  }
  
  # Initialize any random variables with mean
  for (i in slotNames(clust_params)) { 
    if (is.RV(slot(clust_params, i))) 
      slot(clust_params,i) <- mean(slot(clust_params,i)) 
    if (!is.finite(slot(clust_params,i)))
      warning("cluster parameter '", i, "' initialized to infinite value")
  }
  
  return(clust_params)
}


#' Convert relevant attributes in data to a factor representation
#' 
#' @param data A data frame or matrix of observed records. Column names are 
#'   are used to refer to particular columns in the other arguments.
#' @param attr_params A named list of [`Attribute`] objects. Each entry 
#'   in the list specifies the model parameters for an entity attribute, 
#'   and should match one of the column names (attributes) in 
#'   `data`.
#' @return A copy of `data`, where the attributes in `attr_params` are 
#'   encoded as factors, and the levels represent the unique values in the 
#'   domain.
#'   
#' @keywords internal
#' @noRd
data_to_factor <- function(data, attr_params) {
  for (i in seq_along(attr_params)) {
    attr_name <- names(attr_params)[i]
    entity_dist_prior <- attr_params[[i]]@entity_dist_prior
    emp_domain <- unique(as.character(data[[attr_name]]))
    
    # Check whether domain of the attribute is specified in entity_dist_prior
    domainGiven <- FALSE
    if (is.numeric(entity_dist_prior)) {
      # Distribution on entity values is a fixed parameter
      domain <- names(entity_dist_prior)
      domainGiven <- TRUE
    } else if (is.DirichletRV(entity_dist_prior) && length(entity_dist_prior@alpha) > 1) {
      # Distribution on entity values is random with a Dirichlet prior and 
      # concentration parameters (alpha) are not symmetric
      domain <- names(entity_dist_prior@alpha)
      domainGiven <- TRUE
    }
    
    if (domainGiven) {
      # Check validity of given domain
      if (is.null(domain))
        stop(paste0("entity_dist_prior for attribute '", attr_name, "' is missing names for domain values"))
      if (length(setdiff(emp_domain, domain)) != 0)
        stop(paste0("attribute '", attr_name, "' in data contains values not in entity_dist_prior"))
    } else {
      # Use empirical domain
      domain <- emp_domain
    }
    if (length(emp_domain) <= 1) 
      stop("attribute '", attr_name,"' has a domain of size 1 and cannot be used for entity resolution")
    
    data[[attr_name]] <- factor(data[[attr_name]], domain)
  }
  
  return(data)
}