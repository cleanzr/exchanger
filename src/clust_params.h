#pragma once

#include "types.h"
#include "links.h"
#include "rv.h"
#include <boost/optional.hpp>

class Links;

class ClustParams {
public:
  /**
   * Get the prior weight associated with an item spawning a new cluster. This is used in the CRP representation 
   * of the random partition. 
   * @param n_clusters the current number of clusters
   * @return a positive weight (not normalized)
   */
  virtual double prior_weight_new(int n_clusters) const = 0;

  /**
   * Get the prior weight associated with an item joining an existing cluster of a particular size. This is used in 
   * the CRP representation of the random partition. 
   * @param clusterSize the size of the existing cluster
   * @return a positive weight (not normalized)
   */
  virtual double prior_weight_existing(int clusterSize) const = 0;
  
  /**
   * Convert to an R representation
   * @return a `RP` S4 object
   */
  virtual Rcpp::S4 to_R() const = 0;

  /**
   * Convert to an R vector representation
   * @param includeFixed whether to include clustering parameters that are fixed
   * @return a named numeric vector, where the names correspond to clustering parameters
   */
  virtual Rcpp::NumericVector to_R_vec(bool includeFixed=false) const = 0;
  
  /**
   * Update the clustering parameters (if they're random)
   * @param links reference to links container
   */
  virtual void update(Links &links) = 0;

  /**
   * Get the number of clustering parameters that are random (not fixed)
   */
  virtual int num_random() const = 0;
};

class UnifClustParams : public ClustParams {
public:
  virtual double prior_weight_new(int n_clusters) const;

  /**
   * Get the prior weight associated with an item joining an existing cluster of a particular size. This is used in 
   * the CRP representation of the random partition. 
   * @param clusterSize the size of the existing cluster
   * @return a positive weight (not normalized)
   */
  virtual double prior_weight_existing(int clustSize) const;

  /**
   * Constructor doesn't require R objects, as there are no parameters
   */
  UnifClustParams() {}

  /**
   * Convert to an R representation
   * @return a `RP` S4 object
   */
  Rcpp::S4 to_R() const;

  /**
   * Convert to an R vector representation
   * @param includeFixed whether to include clustering parameters that are fixed
   * @return a named numeric vector, where the names correspond to clustering parameters
   */
  virtual Rcpp::NumericVector to_R_vec(bool includeFixed=false) const;

  /**
   * Update the clustering parameters (if they're random)
   * @param links reference to links container
   */
  virtual void update(Links &links);

  /**
   * Get the number of clustering parameters that are random (not fixed)
   */
  virtual int num_random() const;
};


class CouponClustParams : public UnifClustParams {
private:
   int m_;
public:
   /**
    * Get the prior weight associated with an item spawning a new cluster. This is used in the CRP representation 
    * of the random partition. 
    * @param n_clusters the current number of clusters
    * @return a positive weight (not normalized)
    */
   double prior_weight_new(int n_clusters) const;
   
   /**
    * Get the prior weight associated with an item joining an existing cluster of a particular size. This is used in 
    * the CRP representation of the random partition. 
    * @param clusterSize the size of the existing cluster
    * @return a positive weight (not normalized)
    */
   double prior_weight_existing(int clusterSize) const;
   
   /**
    * Constructor from R objects
    * @param clust_prior a `RP` S4 object representing the priors over the cluster parameters
    * @return a CouponClustParams instance
    */
   CouponClustParams(const Rcpp::S4 &clust_prior);
   
   /**
    * Convert to an R representation
    * @return a `RP` S4 object
    */
   Rcpp::S4 to_R() const;
   
   /**
    * Convert to an R vector representation
    * @param includeFixed whether to include clustering parameters that are fixed
    * @return a named numeric vector, where the names correspond to clustering parameters
    */
   virtual Rcpp::NumericVector to_R_vec(bool includeFixed=false) const;
   
   /**
    * Update the clustering parameters (if they're random)
    * @param links reference to links container
    */
   virtual void update(Links &links);
   
   /**
    * Get the number of clustering parameters that are random (not fixed)
    */
   int num_random() const;
};


class GenCouponClustParams : public UnifClustParams {
private:
  boost::optional<ShiftedNegBinomRV> m_prior_;
  boost::optional<GammaRV> kappa_prior_;
  int m_;
  double kappa_;
public:
  /**
   * Get the prior weight associated with an item spawning a new cluster. This is used in the CRP representation 
   * of the random partition. 
   * @param n_clusters the current number of clusters
   * @return a positive weight (not normalized)
   */
  double prior_weight_new(int n_clusters) const;

  /**
   * Get the prior weight associated with an item joining an existing cluster of a particular size. This is used in 
   * the CRP representation of the random partition. 
   * @param clusterSize the size of the existing cluster
   * @return a positive weight (not normalized)
   */
  double prior_weight_existing(int clusterSize) const;

  /**
   * Constructor from R objects
   * @param clust_params a `RP` S4 object representing the instantiated cluster parameters
   * @param clust_prior a `RP` S4 object representing the priors over the cluster parameters
   * @return a GenCouponClustParams instance
   */
  GenCouponClustParams(const Rcpp::S4 &clust_params, const Rcpp::S4 &clust_prior);

  /**
   * Convert to an R representation
   * @return a `RP` S4 object
   */
  Rcpp::S4 to_R() const;

  /**
   * Convert to an R vector representation
   * @param includeFixed whether to include clustering parameters that are fixed
   * @return a named numeric vector, where the names correspond to clustering parameters
   */
  virtual Rcpp::NumericVector to_R_vec(bool includeFixed=false) const;

  /**
   * Update the clustering parameters (if they're random)
   * @param links reference to links container
   */
  virtual void update(Links &links);

  /**
   * Get the number of clustering parameters that are random (not fixed)
   */
  int num_random() const;
};

class PitmanYorClustParams : public UnifClustParams {
private:
  boost::optional<GammaRV> alpha_prior_;
  boost::optional<BetaRV> d_prior_;
  double alpha_;
  double d_;
public:
  /**
   * Get the prior weight associated with an item spawning a new cluster. This is used in the CRP representation 
   * of the random partition. 
   * @param n_clusters the current number of clusters
   * @return a positive weight (not normalized)
   */
  double prior_weight_new(int n_clusters) const;

  /**
   * Get the prior weight associated with an item joining an existing cluster of a particular size. This is used in 
   * the CRP representation of the random partition. 
   * @param clusterSize the size of the existing cluster
   * @return a positive weight (not normalized)
   */
  double prior_weight_existing(int clusterSize) const;

  /**
   * Constructor from R objects
   * @param clust_params a `RP` S4 object representing the instantiated cluster parameters
   * @param clust_prior a `RP` S4 object representing the priors over the cluster parameters
   * @return a PitmanYorClustParams instance
   */
  PitmanYorClustParams(const Rcpp::S4 &clust_params, const Rcpp::S4 &clust_prior);

  /**
   * Convert to an R representation
   * @return a `RP` S4 object
   */
  Rcpp::S4 to_R() const;

  /**
   * Convert to an R vector representation
   * @param includeFixed whether to include clustering parameters that are fixed
   * @return a named numeric vector, where the names correspond to clustering parameters
   */
  Rcpp::NumericVector to_R_vec(bool includeFixed=false) const;

  /**
   * Update the clustering parameters (if they're random)
   * @param links reference to links container
   */
  void update(Links &links);

  /**
   * Get the number of clustering parameters that are random (not fixed)
   */
  int num_random() const;
};

/**
 * Build a ClustParams instance from R `RP` objects
 * @param clust_params a `RP` S4 object representing the instantiated cluster parameters
 * @param clust_prior a `RP` S4 object representing the priors over the cluster parameters
 * @return a shared pointer to a ClustParams instance
 */
std::shared_ptr<ClustParams> read_clust_params(const Rcpp::S4 &clust_params, const Rcpp::S4 &clust_prior);