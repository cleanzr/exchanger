#pragma once

#include "types.h"
#include "discrete_dist.h"
#include <memory>

/**
 * Abstract class for an attribute index. Serves as a container for model 
 * parameters related to an attribute, and contains various utility functions 
 * that are needed to accelerate sampling.
 */
class AbstractAttributeIndex {
public:
  /**
   * Size of the attribute domain
   */
  virtual int domain_size() const = 0;

  /**
   * Get the empirical frequency of a value in the domain
   * @param v value index
   * @return a probability
   */
  virtual double get_empirical_prob(val_id v) const = 0;

  /**
   * Perform a range query for a record value in the domain
   * @param w value index
   * @return a map containing values that fall below the distance threshold 
   * (as keys) and their base distortion probabilities as values.
   */
  virtual const std::unordered_map<val_id, double>& get_close_values(val_id w) const = 0;

  /**
   * Get the base distortion probability for "v distorts to w". 
   * @param v value index for the entity attribute
   * @param w value index for the record attribute
   * @return a probability
   */
  virtual const double get_distortion_prob(val_id v, val_id w) const = 0;

  /**
   * Get the value of max_w exp(-dist(v,w))
   * @param v value index
   * @return a probability
   */
  virtual double get_max_exp_factor(val_id v) const = 0;
  
  /**
   * Whether entity value v is excluded from the support of the distortion 
   * distribution for record value v.
   */
  virtual bool exclude_entity_value() const = 0;
  
  /**
   * Get expected distortion probability for a given record value w.r.t. the 
   * given distribution over the entity values.
   * @param v value index for the record attribute
   * @param a distribution over entity values
   * @return a probability
   */
  virtual double get_exp_distortion_prob(val_id v, DiscreteDist<val_id> *distribution) const = 0;

  /**
   * Get maximum distortion probability for a given record value (over all 
   * possible entity values).
   * @param v value index for the record attribute
   * @return
   */
  virtual double get_max_distortion_prob(val_id v) const = 0;

  /**
   * Update the index based on the latest distribution for the entity 
   * attributes.
   * @param distribution a pointer to a discrete distribution
   */
  virtual void update(DiscreteDist<val_id> *distribution) = 0;

  /**
   * Get a shared pointer to the empirical distribution instance
   */
  virtual std::shared_ptr<DiscreteDist<val_id> > get_empirical_dist() const = 0;

  /**
   * Dirichlet concentration parameter for the distortion distribution
   * @return a ?positive concentration parameter
   */
  virtual double dirichlet_concentration() const = 0;

  virtual ~AbstractAttributeIndex() {}
};


/**
 * Class for a categorical attribute index (with a constant distance function). 
 */
class ConstantAttributeIndex : public AbstractAttributeIndex {
protected:
  /**
   * Represents the empirical distribution
   */
  std::shared_ptr<DiscreteDist<val_id> > dist_;

  /** Dirichlet concentration parameter for the distortion distribution */
  double dirichlet_concentration_;

  /** Whether to exclude entity values in the distortion distributions */
  bool exclude_entity_value_;
  
  /** Expected normalization constant */
  double expected_norm_constant_;
  
public:
  int domain_size() const;
  double get_empirical_prob(val_id v) const;
  std::shared_ptr<DiscreteDist<val_id> > get_empirical_dist() const;
  virtual const std::unordered_map<val_id, double>& get_close_values(val_id w) const;
  virtual const double get_distortion_prob(val_id v, val_id w) const;
  virtual double get_max_exp_factor(val_id v) const;
  virtual double get_max_distortion_prob(val_id v) const;
  bool exclude_entity_value() const;
  virtual void update(DiscreteDist<val_id>* distribution);
  virtual double get_exp_distortion_prob(val_id v, DiscreteDist<val_id> *distribution) const;
  double dirichlet_concentration() const;
  ConstantAttributeIndex(const arma::vec &probs, double dirichlet_concentration, bool exclude_entity_value);
  ConstantAttributeIndex(const ConstantAttributeIndex&) = delete;
  ConstantAttributeIndex& operator=(const ConstantAttributeIndex&) = delete;
};


/**
 * Class for a general attribute index (with a non-constant distance function). 
 */
class NonConstantAttributeIndex : public ConstantAttributeIndex {
protected:
  /**
   * Index used for answering range queries
   */
  const std::vector<std::unordered_map<val_id, double>> close_values_;

  /**
   * Vector containing max_w exp(-dist(v,w)) for each value v
   */
  std::vector<double> max_exp_factor_;
public:
  const std::unordered_map<val_id, double>& get_close_values(val_id w) const;
  const double get_distortion_prob(val_id v, val_id w) const;
  double get_max_exp_factor(val_id v) const;
  void update(DiscreteDist<val_id>* distribution);
  double get_max_distortion_prob(val_id v) const;
  double get_exp_distortion_prob(val_id v, DiscreteDist<val_id> *distribution) const;
  NonConstantAttributeIndex(const arma::vec &probs, const std::vector<std::unordered_map<val_id, double>> &close_values, 
    const Rcpp::NumericVector &max_exp_factor, double dirichlet_concentration, bool exclude_entity_value);
  NonConstantAttributeIndex(const NonConstantAttributeIndex&) = delete;
  NonConstantAttributeIndex& operator=(const NonConstantAttributeIndex&) = delete;
};

/**
 * Builds an attribute index from an `AttributeIndex` R object.
 * @param R_index a reference to the S4 class representation of the index
 * @param dirichlet_concentration 
 * @return a shared pointer to an AbstractAttributeIndex
 */ 
std::shared_ptr<AbstractAttributeIndex> readAttributeIndex(const Rcpp::S4 &R_index, double dirichlet_concentration, bool exclude_entity_value);
