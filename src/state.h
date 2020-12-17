#pragma once

#include <RcppArmadillo.h>
#include "attribute_index.h"
#include "types.h"
#include <iostream>
#include "records.h"
#include "entities.h"
#include "links.h"
#include "distortion_probs.h"
#include "clust_params.h"

/**
 * Container class for the state of the model parameters
 */
struct State {
private:
  int iteration_;
public:
  /**
   * Getter for the iteration counter
   */
  int get_iteration() const { return iteration_; }

  /**
   * Setter for the iteration counter
   */
  void set_iteration(const int i) { iteration_ = i; }
  
  Cache cache_;
  Records recs_;
  Entities ents_;
  DistortProbs distort_probs_;
  Links links_;  
  std::shared_ptr<ClustParams> clust_params_;

  /**
   * Constructor moves all components of the model into this class
   */
  State(int iteration, Entities entities, DistortProbs distort_probs, Links links, Records records, Cache cache, 
    std::shared_ptr<ClustParams> clust_params)
    : iteration_(iteration),
      cache_(std::move(cache)),
      recs_(std::move(records)),
      ents_(std::move(entities)),
      distort_probs_(std::move(distort_probs)),
      links_(std::move(links)),
      clust_params_(std::move(clust_params))
  {
    cache_.update(ents_);
  }

  /**
   * Peform a full Gibbs update (application of the Markov transition operator). 
   * Note that the update order is important, as we are using partially-collapsed Gibbs sampling.
   */
  void update() {
    //Rcpp::Rcout << "Updating links" << std::endl;
    clust_params_->update(links_);
    links_.update(ents_, recs_, clust_params_, cache_.attr_indices_);
    //Rcpp::Rcout << "Updating entities" << std::endl;
    ents_.update_attributes(links_, recs_, distort_probs_, cache_);
    ents_.update_distributions();
    cache_.update(ents_);
    //Rcpp::Rcout << "Updating distortion indicators" << std::endl;
    arma::imat corr_nondistort_counts = recs_.update_distortion(ents_, links_, distort_probs_, cache_);
    //Rcpp::Rcout << "Updating distortion probs" << std::endl;
    distort_probs_.update(recs_, corr_nondistort_counts);
    iteration_++;
  }
};
