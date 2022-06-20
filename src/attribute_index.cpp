#include "attribute_index.h"
#include <memory>
#include <iostream>
#include <RcppArmadillo.h>

ConstantAttributeIndex::ConstantAttributeIndex(const arma::vec &probs, bool exclude_entity_value) 
  : dist_(std::shared_ptr<DiscreteDist<val_id> >(new IndexNonUniformDiscreteDist(probs.cbegin(), probs.cend(), true))), 
    exclude_entity_value_(exclude_entity_value),
    expected_norm_constant_(1.0)
{}

bool ConstantAttributeIndex::exclude_entity_value() const { return exclude_entity_value_; }

int ConstantAttributeIndex::domain_size() const { return dist_->num_items(); }

double ConstantAttributeIndex::get_empirical_prob(val_id v) const { return dist_->get_probability(v); }

std::shared_ptr<DiscreteDist<val_id>> ConstantAttributeIndex::get_empirical_dist() const { return dist_; }

const std::unordered_map<val_id, double>& ConstantAttributeIndex::get_close_values(val_id w) const 
{
  // TODO: this should return all values in the domain, excluding w. Right now, we don't use this?
  return emptyval<std::unordered_map<val_id, double>>::value;
}

void ConstantAttributeIndex::update(DiscreteDist<val_id> *distribution) 
{
  if (exclude_entity_value_) {
    expected_norm_constant_ = 0.0;
    for (val_id v = 0; v < domain_size(); v++) {
      expected_norm_constant_ += distribution->get_probability(v) / (1.0 - get_empirical_prob(v));
    }  
  }
}

const double ConstantAttributeIndex::get_distortion_prob(val_id v, val_id w) const 
{
  if (exclude_entity_value_) {
    if (v == w) {
      return 0.0;
    } else {
      return get_empirical_prob(w) / (1.0 - get_empirical_prob(v));
    }
  } else {
    return get_empirical_prob(w);
  }
}

double ConstantAttributeIndex::get_max_exp_factor(val_id v) const { return 1.0; }

double ConstantAttributeIndex::get_max_distortion_prob(val_id v) const { 
  return 1.0 / (1.0 - dist_->max_probability()); 
}

double ConstantAttributeIndex::get_exp_distortion_prob(val_id v, DiscreteDist<val_id> *distribution) const { 
  double e_prob = get_empirical_prob(v);
  if (exclude_entity_value_) {
    return expected_norm_constant_ * e_prob - distribution->get_probability(v) * e_prob / (1.0 - e_prob);  
  } else {
    return e_prob;
  }
}

NonConstantAttributeIndex::NonConstantAttributeIndex(const arma::vec &probs, 
  const std::vector<std::unordered_map<val_id, double>> &close_values, const Rcpp::NumericVector &max_exp_factor, 
  bool exclude_entity_value) 
  : ConstantAttributeIndex(probs, exclude_entity_value),
    close_values_(close_values)
{
  max_exp_factor_ = std::vector<double>(max_exp_factor.begin(), max_exp_factor.end());
}

const std::unordered_map<val_id, double>& NonConstantAttributeIndex::get_close_values(val_id w) const {
  return close_values_[w];
}

const double NonConstantAttributeIndex::get_distortion_prob(val_id v, val_id w) const {
  if (v == w) {
    if (exclude_entity_value_) {
      return 0.0;
    } else {
      return 1.0;
    }
  } else {
    const auto &cv = close_values_[w];
    const auto got = cv.find(v);
    if ( got != cv.end() ) {
      return got->second;
    }
    return 0.0;
  }
}

double NonConstantAttributeIndex::get_exp_distortion_prob(val_id v, DiscreteDist<val_id> *distribution) const {
  const auto close_values = get_close_values(v);
  double exp = 0.0;
  for (auto const cv : close_values) {
    exp += distribution->get_probability(cv.first) * cv.second;
  }
  return exp;
}

double NonConstantAttributeIndex::get_max_distortion_prob(val_id w) const { 
  const auto &cv = close_values_[w];
  double max = 0.0;
  for (const auto &x : cv) {
    if (x.second > max) max = x.second;
  }
  return max;
}

void NonConstantAttributeIndex::update(DiscreteDist<val_id>* distribution) {}

double NonConstantAttributeIndex::get_max_exp_factor(val_id v) const {
  return max_exp_factor_[v];
}

std::shared_ptr<AbstractAttributeIndex> readAttributeIndex(const Rcpp::S4 &R_index, bool exclude_entity_value) {
  const arma::vec &probs = R_index.slot("probs");
  const Rcpp::List &close_values = R_index.slot("close_values");
  const Rcpp::NumericVector &max_exp_factor = R_index.slot("max_exp_factor");

  AbstractAttributeIndex* index;
  if (close_values.size() > 0) {
    // Invert close_values data structure
    // In R the set {w : expDist(v, w) > 0} of close values to v is located 
    // in entry v of close_values.
    // But for inference we need fast access to {v : expDist(v, w) > 0} for a 
    // given value w (these sets are different in general when the similarity 
    // function is asymmetric).
    std::vector<std::unordered_map<val_id, double>> close_values_(close_values.size());
    for (val_id v = 0; v < close_values.size(); v++) {
      const Rcpp::List& cv = close_values[v];
      const Rcpp::IntegerVector& close_val_ids = cv["valIds"];
      const Rcpp::NumericVector& close_probs = cv["probs"];
      for (int i = 0; i < close_val_ids.size(); i++) {
        val_id w = close_val_ids[i] - 1;
        close_values_[w][v] = close_probs[i];
      }
    }
    index = new NonConstantAttributeIndex(probs, close_values_, max_exp_factor, exclude_entity_value);
  } else {
    index = new ConstantAttributeIndex(probs, exclude_entity_value);
  }

  std::shared_ptr<AbstractAttributeIndex> out(index);

  return out;
}