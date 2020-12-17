#pragma once

#include <RcppArmadillo.h>
#include "alias_sampler.h"
#include "types.h"
#include <tuple>
#include <limits>

/**
 * Abstract class for a discrete distribution over items of arbitrary type
 */
template <class T> 
class DiscreteDist {
public:
  /**
   * Get the number of items in the sample space
   */
  virtual int num_items() const = 0;

  /**
   * Get the total unnormalized weight (effectively the normalization constant)
   */
  virtual double total_weight() const = 0;
  
  /**
   * Get the probability mass associated with an item
   * @param item 
   * @return a probability
   */
  virtual double get_probability(T item) const = 0;

  /**
   * Get the smallest probability mass associated with an item
   */
  virtual double min_probability() const = 0;

  /**
   * Get the largest probability mass associated with an item
   */
  virtual double max_probability() const = 0;

  /**
   * Draw a random variate from the sample space
   * @return a random item
   */
  virtual T draw() const = 0;

  virtual ~DiscreteDist() {}
};


/**
 * Class for a non-uniform discrete distribution over items of arbitrary type.
 * Random variate generation is backed by an alias sampler.
 */
template <class T>
class NonUniformDiscreteDist : public DiscreteDist<T> {
private:
  std::unordered_map<T, double> valuesProbs_;
  std::vector<T> values_;
  AliasSampler alias_;
  double min_prob_;
  double max_prob_;

  /**
   * Initializer for the members
   * @param itemsWeights map from items in the sample space to weights
   * @param normalized boolean that indicates whether the weights are already normalized
   * @return a tuple of members
   */
  static std::tuple<std::unordered_map<T, double>, std::vector<T>, AliasSampler, double, double> 
    init_members(const std::unordered_map<T, double> &itemsWeights, bool normalized);

  /**
   * Constructor based on a tuple of members
   */
  NonUniformDiscreteDist(std::tuple<std::unordered_map<T, double>, std::vector<T>, AliasSampler, double, double> members);
public:
  /**
   * Get the number of items in the sample space
   */
  int num_items() const;

  /**
   * Get the total unnormalized weight (effectively the normalization constant)
   */
  double total_weight() const;

  /**
   * Get the probability mass associated with an item
   * @param item 
   * @return a probability
   */
  double get_probability(T item) const;

  /**
   * Get the smallest probability mass associated with an item
   */
  double min_probability() const;

  /**
   * Get the largest probability mass associated with an item
   */
  double max_probability() const;

  /**
   * Draw a random variate from the sample space
   * @return a random item
   */
  T draw() const;

  /**
   * Constructor
   * @param itemsWeights map from items in the sample space to weights
   * @param normalized boolean that indicates whether the weights are already normalized
   * @return an NonUniformDiscreteDist instance.
   */
  NonUniformDiscreteDist(const std::unordered_map<T, double> &itemsWeights, bool normalized);
};

template <class T>
std::tuple<std::unordered_map<T, double>, std::vector<T>, AliasSampler, double, double> 
  NonUniformDiscreteDist<T>::init_members(const std::unordered_map<T, double> &itemsWeights, bool normalized) 
{
  double total_weight = 0.0;
  std::vector<T> values(itemsWeights.size());
  std::vector<double> probs(itemsWeights.size());  
  std::unordered_map<T, double> itemsProbs;
  itemsProbs.reserve(itemsWeights.size());
  double min_weight = std::numeric_limits<double>::infinity();
  double max_weight = 0.0;

  auto itValues = values.begin();
  auto itProbs = probs.begin();
  for (auto const iw : itemsWeights) 
  {
    *itValues = iw.first;
    *itProbs = iw.second;
    itemsProbs[iw.first] = iw.second;
    total_weight += iw.second;
    if (iw.second > max_weight) { max_weight = iw.second; }
    if (iw.second < min_weight) { min_weight = iw.second; }
    ++itProbs;
    ++itValues;
  }

  if (!normalized) { 
    for (auto ip : itemsProbs) { ip.second *= 1.0/total_weight; } 
    min_weight = min_weight / total_weight;
    max_weight = max_weight / total_weight;
  }

  AliasSampler alias(probs.cbegin(), probs.cend());

  return std::make_tuple(itemsProbs, values, std::move(alias), min_weight, max_weight);
}

template <class T> 
NonUniformDiscreteDist<T>::NonUniformDiscreteDist(std::tuple<std::unordered_map<T, double>, std::vector<T>, AliasSampler, double, double> members) 
  : valuesProbs_{std::get<0>(members)},
    values_{std::get<1>(members)},
    alias_{std::move(std::get<2>(members))},
    min_prob_{std::get<3>(members)},
    max_prob_{std::get<4>(members)}
{}

template <class T> 
NonUniformDiscreteDist<T>::NonUniformDiscreteDist(const std::unordered_map<T, double> &itemsWeights, bool normalized) 
  : NonUniformDiscreteDist{init_members(itemsWeights, normalized)}
{}

template <class T> 
double NonUniformDiscreteDist<T>::total_weight() const { return alias_.total_weight(); }

template <class T> 
int NonUniformDiscreteDist<T>::num_items() const { return values_.size(); }

template <class T> 
double NonUniformDiscreteDist<T>::get_probability(T item) const { return valuesProbs_.at(item); }

template <class T>
double NonUniformDiscreteDist<T>::min_probability() const { return min_prob_; }

template <class T>
double NonUniformDiscreteDist<T>::max_probability() const { return max_prob_; }

template <class T> 
T NonUniformDiscreteDist<T>::draw() const { 
  int idx = alias_.draw();
  return values_[idx];
}

/**
 * Class for a non-uniform discrete distribution over integer-index items. 
 * Random variate generation is backed by an alias sampler.
 */
class IndexNonUniformDiscreteDist : public DiscreteDist<int> {
private:
  std::vector<double> probs_;
  AliasSampler alias_;
  double min_prob_;
  double max_prob_;
  
  /**
   * Initializer for the members
   * @param begin,end input iterator over the weights of the items
   * @param normalized boolean that indicates whether the weights are already normalized
   * @return a tuple of members
   */
  template <class InputIterator>
  static std::tuple<std::vector<double>, AliasSampler, double, double> 
    init_members(InputIterator begin, InputIterator end, bool normalized);

  /**
   * Constructor based on a tuple of members
   */
  IndexNonUniformDiscreteDist(std::tuple<std::vector<double>, AliasSampler, double, double> members);
public:
  /**
   * Get the number of items in the sample space
   */
  int num_items() const;

  /**
   * Get the total unnormalized weight (effectively the normalization constant)
   */
  double total_weight() const;

  /**
   * Get the probability mass associated with an item
   * @param item integer item identifier
   * @return a probability
   */
  double get_probability(int id) const;

  /**
   * Get the smallest probability mass associated with an item
   */
  double min_probability() const;

  /**
   * Get the largest probability mass associated with an item
   */
  double max_probability() const;

  /**
   * Draw a random variate from the sample space
   * @return a random item identifier
   */
  int draw() const;

  /**
   * Constructor
   * @param begin,end input iterator over the weights of the items
   * @param normalized boolean that indicates whether the weights are already normalized
   * @return an IndexNonUniformDiscreteDist instance.
   */
  template <class InputIterator>
  IndexNonUniformDiscreteDist(InputIterator begin, InputIterator end, bool normalized = false);
};

template <class InputIterator>
std::tuple<std::vector<double>, AliasSampler, double, double> 
  IndexNonUniformDiscreteDist::init_members(InputIterator begin, InputIterator end, bool normalized) 
{
  std::vector<double> probs(std::distance(begin, end));
  double min_weight = std::numeric_limits<double>::infinity();
  double max_weight = 0.0;
  double total_weight = 0.0;

  for (auto &p : probs) 
  {
    p = *begin;
    if (p < min_weight) { min_weight = p;}
    if (p > max_weight) { max_weight = p;}
    total_weight += p;
    ++begin;
  }

  if (!normalized) { 
    for (auto &p : probs) { p *= 1.0/total_weight; } 
    min_weight = min_weight / total_weight;
    max_weight = max_weight / total_weight;
  }
  
  AliasSampler alias(probs.cbegin(), probs.cend());

  return std::make_tuple(probs, std::move(alias), min_weight, max_weight);
}

template <class InputIterator>
IndexNonUniformDiscreteDist::IndexNonUniformDiscreteDist(InputIterator begin, InputIterator end, bool normalized) 
  : IndexNonUniformDiscreteDist{init_members(begin, end, normalized)}
{}
