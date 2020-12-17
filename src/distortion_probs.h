#pragma once

#include "types.h"
#include <RcppArmadillo.h>
#include "records.h"
#include "rv.h"

class Records;

/**
 * Represents the distortion probabilities for each file/attribute
 */
class DistortProbs {
private:
  /**
   * Initialize the Beta priors from the R representation
   * @param attr_params  a list of `Attribute` S4 instances.
   * @return a vector of BetaRV structs---one for each attribute.
   */
  static std::vector<BetaRV> init_prior(const Rcpp::List &attr_params);

  /**
   * A matrix of distortion probabilities. Rows correspond to files and 
   * columns correspond to attributes.
   */
  arma::mat probs_;

  /**
   * A vector of Beta priors for each attribute
   */
  std::vector<BetaRV> distort_prob_priors_;

public:
  /**
   * Getter for a distortion probability
   * @param fid file identifier
   * @param aid attribute identifier
   * @return the value of the probability
   */
  double get(file_id fid, attr_id aid) const;

  /**
   * Perform a Gibbs update for the distortion probabilities
   * @param recs records container
   * @param corr_nondistort_counts a matrix containing "corrected" distortion 
   * counts for each file/attribute
   */
  void update(const Records& recs, const arma::imat &corr_nondistort_counts);

  /**
   * Get an R representation of the distortion probabilities
   * @return an R numeric matrix
   */
  Rcpp::NumericMatrix to_R() const;

  /**
   * Constructor from the R representation
   * @param probs a matrix of distortion probabilities. Rows correspond to 
   * files and columns correspond to attributes.
   * @param attr_params a list of `Attribute` S4 instances.
   */
  DistortProbs(const arma::mat &probs, const Rcpp::List &attr_params) 
    : probs_(probs),
      distort_prob_priors_{init_prior(attr_params)}
  {}
};


