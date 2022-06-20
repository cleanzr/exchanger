#pragma once

#include "types.h"
#include "rv.h"
#include "entities.h"
#include "records.h"
#include "links.h"
#include "cache.h"
#include "attribute_index.h"
#include <boost/optional.hpp>

class Records;
class Entities;
class Links;
class Cache;

/**
 * Represents the distortion distribution concentration parameters for each 
 * attribute
 */
class DistortDistConcs {
private:
  /**
   * Initialize the Gamma priors from the R representation
   * @param attr_params  a list of `Attribute` S4 instances.
   * @return a vector of GammaRV structs---one for each attribute.
   */
  static std::vector<boost::optional<GammaRV>> init_prior(const Rcpp::List &attr_params);

  /**
   * A vector of concentration parameters for each attribute.
   */
  arma::vec values_;

  /**
   * A vector of Beta priors for each attribute
   */
  std::vector<boost::optional<GammaRV>> priors_;

public:
  /**
   * Getter for a concentration parameter
   * @param aid attribute identifier
   * @return the value of the parameter
   */
  double get(attr_id aid) const;

  /**
   * Perform a Gibbs update for the distortion probabilities
   * @param recs links container
   * @param recs records container
   * @param ents entities container
   * @param cache cache container
   */
  void update(const Links &links, const Records &recs, const Entities &ents, const Cache &cache);

  /**
   * Get an R representation of the concentration parameters
   * @return an R numeric vector
   */
  Rcpp::NumericVector to_R() const;

  /**
   * Constructor from the R representation
   * @param probs a vector of distortion distribution concentration 
   *   parameters---one for each attribute.
   * @param attr_params a list of `Attribute` S4 instances.
   */
  DistortDistConcs(const arma::vec &values, const Rcpp::List &attr_params) 
    : values_(values),
      priors_{init_prior(attr_params)}
  {}
};