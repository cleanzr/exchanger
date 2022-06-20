#pragma once

#include "attribute_index.h"
#include "types.h"
#include <RcppArmadillo.h>
#include "records.h"
#include "links.h"
#include "distortion_probs.h"
#include "distort_dist_concs.h"
#include "cache.h"
#include <memory>

class Records;
class Cache;
class DistortProbs;
class DistortDistConcs;
class Links;

/**
 * Represents the prior over an entity attribute generative distribution. 
 * May be a Dirichlet distribution or fixed.
 */
struct AttributePrior {
  public:
    std::vector<double> param_vec_;
    bool isDirichlet_;

    template <class InputIterator>
    AttributePrior(InputIterator begin, InputIterator end, bool isDirichlet) 
      : isDirichlet_(isDirichlet)
    {
      param_vec_.resize(std::distance(begin, end));
      for (auto &p : param_vec_) 
      {
        p = *begin;
        ++begin;
      }
    }
};

/**
 * Represents the entity attribute values and their generative distributions
 */
class Entities {
private:
  ent_id next_id_;
  std::list<ent_id> unused_ids_;
  int total_;
  int n_attributes_;
  std::unordered_map<ent_id, attr_vec> forward_index_;
  std::vector<std::vector<std::unordered_set<ent_id> > > inverse_index_;
  std::vector<std::shared_ptr<IndexNonUniformDiscreteDist> > distributions_;
  std::vector<std::shared_ptr<AttributePrior> > attr_priors_;
  ent_id next_id();
public:
  /**
   * Add an entity to the population
   * @param avec a vector containing the entity's attribute values
   * @return the identifier of the entity
   */
  ent_id add(const attr_vec& avec);

  /**
   * Remove an entity from the population
   * @param eid identifier of the entity
   */
  void remove(ent_id eid);

  /**
   * Getter for the complete set of entity identifiers
   * @return a vector containing all entity identifiers
   */
  std::vector<ent_id> get_entity_ids() const;

  /**
   * Query the set of entities that have a given value for an attribute
   * @param aid attribute identifier
   * @param vid value identifier
   * @return a set containing identifiers for entities that satisfy the query
   */
  const std::unordered_set<ent_id>& get_entity_ids(attr_id aid, val_id vid) const;

  /**
   * Getter for an entity's attribute values
   * @param eid entity's identifier
   * @return a vector of attribute values
   */
  const attr_vec& get_attributes(ent_id eid) const;

  /**
   * Getter for an attribute value of an entity
   * @param eid entity's identifier
   * @param aid attribute identifier
   * @return a value identifier
   */
  val_id get_attribute_value(ent_id eid, attr_id aid) const;

  /**
   * Return the current size of the entity population
   */
  int n_entities() const;

  /**
   * Convert to an R representation
   */
  std::pair<Rcpp::IntegerVector, Rcpp::IntegerMatrix> to_R() const;

  /** 
   * Peform a Gibbs update for the distributions over the entity attributes
   */
  void update_distributions();

  /** 
   * Peform a Gibbs update for the entity attributes
   * @param links reference to the linkage structure
   * @param recs reference to the records
   * @param distort_probs reference to the distortion probabilities
   * @param distort_dist_concs reference to the distortion dist concentration parameters
   * @param cache reference to the cache
   */
  void update_attributes(const Links &links, const Records &recs, const DistortProbs &distort_probs, 
                         const DistortDistConcs &distort_dist_concs, const Cache &cache);

  /** 
   * Getter for an entity attribute distribution 
   * @param aid attribute identifier
   * @return a discrete distribution for the specified attribute
   */
  IndexNonUniformDiscreteDist* get_distribution(attr_id aid) const;

private:
  /**
   * Setter for an entity attribute value
   */
  void set_attribute_value(ent_id eid, attr_id aid, val_id vid);

  /**
   * Add to the forward index
   */
  void add_forward(ent_id eid, attr_vec avec);

  /**
   * Add to the inverted index
   */
  void add_inverted(attr_id aid, val_id vid, ent_id eid);

  /**
   * Remove from the inverted index
   */
  void remove_inverted(attr_id aid, val_id vid, ent_id eid);
public:
  //Entities(const arma::ivec &eids, const arma::imat &attr_vals, const Cache &cache);

  /**
   * Constructor from R objects
   * @param eids vector of entity identifiers
   * @param attr_vals matrix of entity attribute values
   * @param attr_params list of `Attribute` S4 objects
   * @param attr_indices list of `AttributeIndex` S4 objects
   */
  Entities(const Rcpp::IntegerVector &eids, const arma::imat &attr_vals, 
           const Rcpp::List &attr_params, const Rcpp::List &attr_indices);
};