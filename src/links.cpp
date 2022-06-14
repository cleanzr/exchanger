#include "links.h"
#include "clust_params.h"
#include <RcppArmadillo.h>

void Links::set_link(rec_id rid, ent_id eid) {
  ent_id &old_eid = forward_index_[rid];
  if (eid != old_eid) {
    remove_inverted(rid, old_eid);
    add_inverted(rid, eid);
    old_eid = eid;
  }
}

void Links::add_inverted(rec_id rid, ent_id eid) {
  // Try to insert new key-value pair if eid doesn't already exist in the index
  std::pair<std::unordered_map<ent_id, std::unordered_set<rec_id>>::iterator, bool> result = 
    inverse_index_.insert(std::make_pair(eid, std::unordered_set<rec_id> {rid}));
  
  // If we didn't insert anything, key already exists. Need to insert in existing set.
  if (!result.second) { (result.first)->second.insert(rid); }
}

void Links::remove_inverted(rec_id rid, ent_id eid) {
  // Lookup eid in index
  std::unordered_map<ent_id, std::unordered_set<rec_id>>::iterator got = 
    inverse_index_.find(eid);
  
  if ( got != inverse_index_.end() ) {
    // Entry exists in map for eid
    // Now remove rid from set
    got->second.erase(rid);
    
    // If the set becomes empty, remove entry for eid in map
    if (got->second.empty()) { inverse_index_.erase(got); }
  } else {
    Rcpp::warning("Links: Trying to remove entry from inverted index that doesn't exist.");
  }
}

ent_id Links::get_entity_id(rec_id rid) const { return forward_index_[rid]; }

const std::unordered_set<rec_id>& Links::get_record_ids(ent_id eid) const {
  std::unordered_map<ent_id, std::unordered_set<rec_id>>::const_iterator got = inverse_index_.find(eid);

  if ( got != inverse_index_.end() ) { 
    return got->second;
  }
  
  return emptyval<std::unordered_set<rec_id> >::value;
}

Rcpp::IntegerVector Links::to_R() const {
  return Rcpp::wrap(forward_index_);
}

int Links::n_linked_ents() const {
  return inverse_index_.size();
}

int Links::n_records(ent_id eid) const {
  std::unordered_map<ent_id, std::unordered_set<rec_id>>::const_iterator got = 
    inverse_index_.find(eid);
  
  if ( got != inverse_index_.end() ) { return got->second.size(); }

  return 0;
}

int Links::n_records() const {
  return forward_index_.size();
}

Links::Links(const Rcpp::IntegerVector &forward_index) 
{
  forward_index_.reserve(forward_index.size());
  
  rec_id rid = 0;
  for (auto const eid : forward_index) {
    forward_index_.push_back(eid);
    add_inverted(rid, eid);
    rid++;
  }
}

template <class T> 
bool cmpSets(const std::unordered_set<T>* a, const std::unordered_set<T>* b) { return (a->size() < b->size()); }

std::vector<ent_id> getPossibleLinks(rec_id rid, const Records &recs, 
  const Entities &ents, const std::vector<std::shared_ptr<AbstractAttributeIndex> > &attr_indices) 
{ 
  // Initialize an array which will store sets of entity ids that are 
  // candidate matches for the given record. The intersection of all sets 
  // will be the final set of candidates.
  std::vector<const std::unordered_set<ent_id> * > sets;
  int smallestSetSize = ents.n_entities();
  
#ifdef FILTER_CANDIDATES_BY_DIST
  std::vector<std::unordered_set<ent_id> > new_sets;
  // Keep track of the smallest set. If it reaches zero, we can break early.
  const AbstractAttributeIndex* index;
  new_sets.reserve(recs.n_attributes());
#endif
  sets.reserve(recs.n_attributes());
  val_id r_vid;
  
  for (attr_id aid = 0; aid < recs.n_attributes(); aid++) 
  {
    if (smallestSetSize == 0) {
      // Multi-set intersection will yield no candidates. Return now.
      return std::vector<ent_id>();
    }
    
    r_vid = recs.get_attribute(rid, aid);
    if (r_vid != NA_INTEGER) 
    {
      if (!recs.get_attribute_distortion(rid, aid)) 
      {
        // Record attribute is not distorted. Only entities that match exactly 
        // on this attribute are possible.
        const std::unordered_set<ent_id> &matchEntIds = ents.get_entity_ids(aid, r_vid);
        if (matchEntIds.size() < smallestSetSize) {
          smallestSetSize = matchEntIds.size();
        }
        sets.push_back(&matchEntIds);
      } 
      else
      {
#ifdef FILTER_CANDIDATES_BY_DIST
        // Record attribute is distorted. Assuming a thresholded distance funtion,
        // only entities that are a "close" match on this attribute are possible.
        index = attr_indices[aid].get();

        // It's only worth proceeding here for a non-constant attribute, specifically
        // one for which the sets of "close" values are small. If this is not the
        // case, then we'll end up wasting time preparing a large set of entity ids,
        // which won't help with pruning.
        if (typeid(*index) != typeid(ConstantAttributeIndex)) {
          const std::unordered_map<val_id, double> &close_values = index->get_close_values(r_vid);
          
          bool abandonedLoop = false;
          std::unordered_set<ent_id> matchEntIds;
          for (const auto &vw : close_values) {
            const std::unordered_set<ent_id> &eIds = ents.get_entity_ids(aid, vw.first);
            if (eIds.size() > smallestSetSize) {
            // if (eIds.size() > smallestSetSize || matchEntIds.size() > smallestSetSize) {
              // Not efficient to continue building the set.
              abandonedLoop = true;
              break;
            }
            matchEntIds.insert(eIds.begin(), eIds.end());
          }
          
          if (abandonedLoop) {
            // Didn't finish building matchEntIds. Skip to next attribute
            continue; 
          }
          
          if (matchEntIds.size() < smallestSetSize) {
            smallestSetSize = matchEntIds.size();
          }
          
          new_sets.push_back(matchEntIds);
          sets.push_back(&new_sets.back());
        }
#endif
      }
    }
  }
  
  // Now compute the multiple set intersection, but first handle special cases
  if (sets.empty()) 
  {
    // Haven't been able to prune the possible links at all (possibly because 
    // all attributes were unobserved, or it wasn't feasible to prune for some 
    // distorted attributes).
    // Have to return the set of all entities.
    std::vector<ent_id> allEntityIds = ents.get_entity_ids();
    return allEntityIds;
  }

  // Sort sets in increasing order of size to improve the efficiency of the multiple set
  // intersection algorithm
  std::sort(sets.begin(), sets.end(), cmpSets<ent_id>);

  std::vector<const std::unordered_set<ent_id> * >::iterator itSets = sets.begin();
  std::vector<const std::unordered_set<ent_id> * >::iterator itSetsEnd = sets.end();
  std::vector<ent_id> result( (*itSets)->begin(), (*itSets)->end() );
  ++itSets;

  if (itSets == itSetsEnd) { return result; }

  std::vector<ent_id> buffer;
  for (; itSets != itSetsEnd; ++itSets) 
  {
    buffer.clear();
    for (auto const &eid : result) 
    {
      // Look up in second set
      std::unordered_set<ent_id>::const_iterator got = (*itSets)->find(eid);
      if (got != (*itSets)->end()) { buffer.push_back(eid); } // entity is in set
    }
    std::swap(result, buffer);
  }
  
  return result;
}



double likelihoodWeightExisting(rec_id rid, ent_id eid, const std::unordered_set<rec_id> &linked_rids, 
  const Records &recs, const Entities &ents, const std::vector<std::shared_ptr<AbstractAttributeIndex> > &attr_indices) 
{
  double weight = 1.0;
  val_id e_vid, r_vid, linked_rvid;
  attr_id aid;
  const AbstractAttributeIndex* index;

  // Factor for each attribute
  for (aid = 0; aid < recs.n_attributes(); aid++) 
  {
    if (weight == 0.0) break;
    
    index = attr_indices[aid].get();
    e_vid = ents.get_attribute_value(eid, aid);
    
    // Counts number of observed, distorted record attribute values linked to 
    // this entity
    int ctrDistorted = 0; 
    // Same as above, but separates counts per unique record attribute value
    std::unordered_map<val_id, int> ctrDistortedValues;

    // Contribution associated with records linked to the entity, with 
    // attribute values that are observed and distorted
    double b = index->dirichlet_concentration();
    if (b != R_PosInf) 
    {
      for (auto const &linked_rid : linked_rids) 
      {
        // Need to filter out rid, since we didn't remove it from the links 
        // data structure
        if (linked_rid != rid) { 
          linked_rvid = recs.get_attribute(linked_rid, aid);
          if (recs.get_attribute_distortion(linked_rid, aid) && linked_rvid != NA_INTEGER) {
            // Linked record attribute value is distorted and observed
            ctrDistorted++;
            if (ctrDistorted == 1) {
              // Guard against floating point rounding errors
              weight *= 1.0 / b;
            } else {
              weight *= 1.0 / (1.0 + (b - 1.0)/ctrDistorted);
            }
            auto got = ctrDistortedValues.find(linked_rvid);
            if (got == ctrDistortedValues.end()) {
              ctrDistortedValues[linked_rvid] = 1;
              weight *= b * index->get_distortion_prob(e_vid, linked_rvid);
            } else {
              got->second++;
              weight *= 1.0 + (b * index->get_distortion_prob(e_vid, linked_rvid) - 1.0) / got->second;
            }
          }
        }
      }
    }
    
    // Contribution associated with the record, assuming it were reassigned to 
    // this entity
    r_vid = recs.get_attribute(rid, aid);
    if (recs.get_attribute_distortion(rid, aid) && r_vid != NA_INTEGER) {
      // Record attribute value is distorted and observed
      if (b != R_PosInf) {
        if (ctrDistorted == 0) {
          // Guard against floating point rounding errors
          weight *= 1.0 / b;
        } else {
          weight *= 1.0 / (1.0 - 1.0 / (ctrDistorted + 1) + b);
        }
        auto got = ctrDistortedValues.find(r_vid);
        if (got == ctrDistortedValues.end()) {
          weight *= b * index->get_distortion_prob(e_vid, r_vid);
        } else {
          weight *= 1.0 - 1.0 / (got->second + 1) + b * index->get_distortion_prob(e_vid, r_vid);
        }
      } else {
        weight *= index->get_distortion_prob(e_vid, r_vid);
      }
    }
  }
  #ifdef CHECK_WEIGHTS
  if (weight < 0.0 || !std::isfinite(weight)) { 
    std::string message = "invalid weight " +  std::to_string(weight) + " on entity " + std::to_string(eid) + " for record " + std::to_string(rid); 
    throw std::runtime_error(message);
  }
  #endif
  return weight;
}

val_id rejectionSampleConstant(val_id r_vid, const IndexNonUniformDiscreteDist* distribution, 
  const AbstractAttributeIndex* index) 
{
  val_id e_vid;
  double ratio;
  bool accepted = false;
  while (!accepted) {
    e_vid = distribution->draw();
    if (e_vid == r_vid) {
      // ratio = 0.0
      continue;
    } else {
      ratio = index->get_distortion_prob(e_vid, r_vid) / index->get_max_distortion_prob(r_vid);
    }
    if (R::unif_rand() < ratio) { accepted = true; }
  }
  return e_vid;
}

double likelihoodWeightNew(rec_id rid, const Records &recs, 
  const std::vector<std::shared_ptr<AbstractAttributeIndex> > &attr_indices, const Entities &ents) 
{
  double weight = 1.0;
  val_id r_vid;
  const AbstractAttributeIndex *index;
  IndexNonUniformDiscreteDist *distribution;
  for (attr_id aid = 0; aid < recs.n_attributes(); aid++)
  {
    r_vid = recs.get_attribute(rid, aid);
    if (r_vid != NA_INTEGER)
    {
      index = attr_indices[aid].get();
      distribution = ents.get_distribution(aid);
      if (recs.get_attribute_distortion(rid, aid)) 
      {
        // Record attribute value is distorted and observed
        weight *= index->get_exp_distortion_prob(r_vid, distribution);
      } 
      else 
      {
        // Record attribute value is observed, but not distorted
        weight *= distribution->get_probability(r_vid);
      }
    }
  }
  #ifdef CHECK_WEIGHTS
  if (weight <= 0.0 || !std::isfinite(weight)) { 
    std::string message = "invalid weight " +  std::to_string(weight) + " on new entity for record " + std::to_string(rid); 
    throw std::runtime_error(message);
  }
  #endif
  return weight;
}


void Links::update_link(rec_id rid, Entities &ents, const Records &recs, std::shared_ptr<ClustParams> clust_params, 
  const std::vector<std::shared_ptr<AbstractAttributeIndex> > &attr_indices) 
{
  ent_id old_linked_eid = get_entity_id(rid);
  if (n_records(old_linked_eid) <= 1) { ents.remove(old_linked_eid); } // remove corresponding entity values

  std::vector<ent_id> pLinks = getPossibleLinks(rid, recs, ents, attr_indices);

  // First weights correspond to items in pLinks (which are existing entities).
  // Last weight corresponds to a new entity.
  std::vector<double> weights(pLinks.size() + 1);
  size_t i = 0;

  // Compute weights vector
  // Weights associated with existing entities
  for (auto const &eid : pLinks) 
  {
    weights[i] = clust_params->prior_weight_existing(n_records(eid) - (eid == old_linked_eid));
    const auto &linked_rids = get_record_ids(eid);
    weights[i] *= likelihoodWeightExisting(rid, eid, linked_rids, recs, ents, attr_indices);
    i++;
  }
  // Weight associated with "new" entity.
  weights[i] = clust_params->prior_weight_new(ents.n_entities());
  weights[i] *= likelihoodWeightNew(rid, recs, attr_indices, ents);
  
  ent_id newLink;
  IndexNonUniformDiscreteDist dist(weights.cbegin(), weights.cend(), false);
  int randIdx = dist.draw();

  if (randIdx == pLinks.size())
  {
    // Selected "new" entity. Below we draw attributes for this new entity according to the posterior distribution.
    const AbstractAttributeIndex* index;
    IndexNonUniformDiscreteDist* distribution;
    val_id r_vid;
    std::vector<val_id> attributeVec(recs.n_attributes());
    for (attr_id aid = 0; aid < recs.n_attributes(); aid++) {
      r_vid = recs.get_attribute(rid, aid);
      distribution = ents.get_distribution(aid);
      if ( r_vid == NA_INTEGER ) {
        // Record value is unobserved, so draw entity value according to prior
        index = attr_indices[aid].get();
        attributeVec[aid] = distribution->draw();
      } else {
        // Record value is observed, so draw entity value according to posterior
        if ( !recs.get_attribute_distortion(rid, aid) ) {
          // Non-distorted, so entity value is fixed to record value
          attributeVec[aid] = r_vid;
        } else {
          // Distorted
          index = attr_indices[aid].get();
          if (typeid(*index) == typeid(ConstantAttributeIndex)) {
            if (index->exclude_entity_value()) {
              attributeVec[aid] = rejectionSampleConstant(r_vid, distribution, index);
            } else {
              std::shared_ptr<DiscreteDist<val_id> > e_dist = index->get_empirical_dist();
              attributeVec[aid] = e_dist->draw();
            }
          } else {
            std::unordered_map<val_id, double> pmf;
            const auto close_values = index->get_close_values(r_vid);
            for (auto const cv : close_values) 
            {
              pmf[cv.first] = distribution->get_probability(cv.first) * cv.second;
            }
            NonUniformDiscreteDist<val_id> dist(pmf, false);
            attributeVec[aid] = dist.draw();
          }
        }
      }
    }
    // Update index and get entity id for new link
    newLink = ents.add(attributeVec);
  }
  else 
  {
    // Selected existing entity
    newLink = pLinks[randIdx];
  }
  set_link(rid, newLink);
}

void Links::update(Entities &ents, const Records &recs, std::shared_ptr<ClustParams> clust_params,
  const std::vector<std::shared_ptr<AbstractAttributeIndex> > &attr_indices)
{
  for (rec_id rid=0; rid < recs.n_records(); rid++) 
  {
    update_link(rid, ents, recs, clust_params, attr_indices);
  }
}