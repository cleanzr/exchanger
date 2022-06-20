#include "entities.h"
#include "dirichlet.h"
#include "distort_dist_concs.h"
#include <memory>

void Entities::add_inverted(attr_id aid, val_id vid, ent_id eid) 
{
  std::unordered_set<ent_id> &eids = inverse_index_[aid][vid];
  eids.insert(eid);
}

void Entities::remove_inverted(attr_id aid, val_id vid, ent_id eid) 
{
  std::unordered_set<ent_id> &eids = inverse_index_[aid][vid];
  eids.erase(eid);
}

ent_id Entities::add(const attr_vec& avec) 
{
  ent_id eid = next_id();
  
  // Add to forward_index
  add_forward(eid, avec);

  // Add to inverse_index
  attr_id aid = 0;
  for (auto const vid : avec) {
    add_inverted(aid, vid, eid);
    aid++;
  }

  total_++;
  return eid;
}

void Entities::add_forward(ent_id eid, attr_vec avec) {
  forward_index_[eid] = avec;
}

ent_id Entities::next_id() {
  ent_id returnId;
  if (unused_ids_.empty()) {
    returnId = next_id_;
    next_id_++;
  } else {
    returnId = unused_ids_.front();
    unused_ids_.pop_front();
  }
  return returnId;
}

void Entities::remove(ent_id eid) {
  // Use forward_index to remove entries in inverse_index efficiently
  std::unordered_map<ent_id, attr_vec>::iterator got = forward_index_.find(eid);
  
  if ( got != forward_index_.end() ) {
    attr_vec& attributes = got->second;
    attr_id aid = 0;
    for (auto const vid : attributes) {
      remove_inverted(aid, vid, eid);
      aid++;
    }

    // Finally remove from forward_index
    forward_index_.erase(got);
    total_--;
    unused_ids_.push_back(eid);
  }
}

std::vector<ent_id> Entities::get_entity_ids() const {
  // Construct set of entity ids based on forward_index
  std::vector<ent_id> allEntityIds;
  allEntityIds.reserve(total_);
  for (auto const &pair : forward_index_) { allEntityIds.push_back(pair.first); }
  return allEntityIds;
}

const std::unordered_set<ent_id>& Entities::get_entity_ids(attr_id aid, val_id vid) const {
  const std::unordered_set<ent_id> &eids = inverse_index_[aid][vid];
  return eids;
}

const attr_vec& Entities::get_attributes(ent_id eid) const {
  std::unordered_map<ent_id, attr_vec>::const_iterator got = forward_index_.find(eid);

  if ( got != forward_index_.end() ) { 
    return got->second; 
  }
  
  return emptyval<attr_vec>::value;
}

val_id Entities::get_attribute_value(ent_id eid, attr_id aid) const {
  const attr_vec& avec = get_attributes(eid);
  return avec[aid];
}

void Entities::set_attribute_value(ent_id eid, attr_id aid, val_id vid) {
  std::unordered_map<ent_id, attr_vec>::iterator got = forward_index_.find(eid);

  if ( got != forward_index_.end() ) {
    val_id &old_vid = got->second[aid];
    if (vid != old_vid) {
      // Need to update inverse_index and forward_index
      remove_inverted(aid, old_vid, eid);
      add_inverted(aid, vid, eid);
      old_vid = vid;
    }
  } else {
    Rcpp::warning("Entities: Trying to set entity attribute that doesn't exist");
  }
}

int Entities::n_entities() const { return total_; }

std::pair<Rcpp::IntegerVector, Rcpp::IntegerMatrix> Entities::to_R() const {
  Rcpp::IntegerVector eids(total_);
  arma::imat attr_vals(n_attributes_, total_);

  int i = 0;
  for (auto const &entry : forward_index_) {
    eids[i] = entry.first;
    attr_vals.col(i) = entry.second + 1;
    i += 1;
  }

  Rcpp::IntegerMatrix Rmat = Rcpp::wrap(attr_vals);

  return std::make_pair(eids, Rmat);
}

void Entities::update_distributions() {
  attr_id aid = 0;
  for (auto const &v : inverse_index_) {
    AttributePrior* attributePrior = attr_priors_[aid].get();
    if (attributePrior->isDirichlet_) 
    {
      // Need to update distribution for this attribute
      std::vector<double> concentration(v.size()); // hard-coded prior
    
      auto priorParamPtr = attributePrior->param_vec_.begin();
      auto concentrationPtr = concentration.begin();
      for (auto const &eids : v) {
        *concentrationPtr = eids.size() + *priorParamPtr;
        concentrationPtr++;
        priorParamPtr++;
      }

      std::vector<double> pmf = rdirichlet(concentration.cbegin(), concentration.cend());
      std::shared_ptr<IndexNonUniformDiscreteDist> &distPtr = distributions_[aid];
      distPtr.reset(new IndexNonUniformDiscreteDist(pmf.cbegin(), pmf.cend(), true));
    }
    aid++;
  }
}

IndexNonUniformDiscreteDist* Entities::get_distribution(attr_id aid) const {
  const std::shared_ptr<IndexNonUniformDiscreteDist> &ptr = distributions_[aid];
  return ptr.get();
}

Entities::Entities(const Rcpp::IntegerVector &eids, const arma::imat &attr_vals, 
           const Rcpp::List &attr_params, const Rcpp::List &attr_indices) 
  : next_id_(0),
    total_(0),
    n_attributes_(attr_vals.n_rows)
{
  inverse_index_.resize(n_attributes_); 
  distributions_.resize(n_attributes_);
  attr_priors_.resize(n_attributes_);

  attr_id aid = 0;
  for ( ; aid < n_attributes_; aid++ ) {
    const Rcpp::S4 &index_R = attr_indices[aid];
    const Rcpp::CharacterVector &domain_R = index_R.slot("domain");

    // Read domain size
    inverse_index_[aid].resize(domain_R.length());

    // Read attribute prior
    const Rcpp::S4 &spec_R = attr_params[aid];
    const SEXP &entity_dist_prior_SEXP = spec_R.slot("entity_dist_prior");
    Rcpp::NumericVector param_vec;
    bool isDirichlet = false;
    if (TYPEOF(entity_dist_prior_SEXP) == S4SXP) 
    {
      const Rcpp::S4 &entity_dist_prior = Rcpp::as<Rcpp::S4>(entity_dist_prior_SEXP);
      if (!entity_dist_prior.is("DirichletRV")) Rcpp::stop("Unrecognized entity_dist_prior");
      isDirichlet = true;
      param_vec = entity_dist_prior.slot("alpha");
      if (param_vec.length() != domain_R.length() && param_vec.length() == 1) {
        param_vec = Rcpp::rep(param_vec[0], domain_R.length());
      } else {
        Rcpp::stop("Invalid entity_dist_prior");
      }
    } 
    else if (TYPEOF(entity_dist_prior_SEXP) == REALSXP)
    {
      param_vec = Rcpp::as<Rcpp::NumericVector>(entity_dist_prior_SEXP);
    } 
    else if (TYPEOF(entity_dist_prior_SEXP) == NILSXP) 
    {
      param_vec = index_R.slot("probs");
    }
    attr_priors_[aid].reset(new AttributePrior(param_vec.begin(), param_vec.end(), isDirichlet));

    // Instantiate distribution for non-Dirichlet (Dirichlet is handled by update method)
    if (!isDirichlet) {
      distributions_[aid].reset(new IndexNonUniformDiscreteDist(param_vec.begin(), param_vec.end(), true));
    }
  }

  // Populate forward and inverse indexes given initial state
  for (unsigned int i=0; i < attr_vals.n_cols; i++) {
    const arma::ivec &avec = attr_vals.col(i);
    ent_id eid = eids[i];

    add_forward(eid, avec);

    aid = 0;
    for (auto const vid : avec) {
      add_inverted(aid, vid, eid);
      aid++;
    }

    total_++;
    if (eid >= next_id_) { next_id_ = eid + 1; }
  }

  // Update distributions for attributes with Dirichlet prior (this must come 
  // after indexes are populated)
  update_distributions();
}

double logWeightEntityValue(
  val_id e_vid, 
  attr_id aid, 
  const Records &recs, 
  const DistortProbs &distort_probs, 
  const DistortDistConcs &distort_dist_concs,
  const std::unordered_set<rec_id> &linked_rids, 
  const AbstractAttributeIndex *index, 
  IndexNonUniformDiscreteDist *distribution
) 
{
  double log_weight, eff_dist_prob;
  val_id r_vid;
  file_id fid;
  int ctr_disagree = 0;
  std::unordered_map<val_id, int> ctr_disagree_value;
  
  double max_exp_factor = index->get_max_exp_factor(e_vid);
  log_weight = std::log(distribution->get_probability(e_vid));
  
  // Loop over records linked to this entity
  for (auto const &linked_rid: linked_rids) 
  {
    fid = recs.get_file_id(linked_rid);
    r_vid = recs.get_attribute(linked_rid, aid);
    if (r_vid != NA_INTEGER) 
    {
      // Record attribute value is observed
      eff_dist_prob = distort_probs.get(fid, aid) * max_exp_factor;
      if (r_vid == e_vid) 
      {
        if (index->exclude_entity_value())
          log_weight += std::log(1.0 - eff_dist_prob);
        else
          log_weight += std::log(1.0 - eff_dist_prob + eff_dist_prob * index->get_distortion_prob(e_vid, r_vid));
      } 
      else 
      {
        log_weight += std::log(eff_dist_prob);
        double b = distort_dist_concs.get(aid);
        if (R_FINITE(b)) { // This also guarantees index->exclude_entity_value() is TRUE, checked elsewhere
          ctr_disagree++;
          log_weight += -std::log(ctr_disagree - 1 + b);
          auto got = ctr_disagree_value.find(r_vid);
          if (got != ctr_disagree_value.end()) {
            // Seen this r_vid before
            got->second++;
            log_weight += std::log(got->second - 1 + b * index->get_distortion_prob(e_vid, r_vid));
          } else {
            // Seeing this r_vid for the first time
            ctr_disagree_value[r_vid] = 1;
            log_weight += std::log(b * index->get_distortion_prob(e_vid, r_vid));
          }
        } else {
          log_weight += std::log(index->get_distortion_prob(e_vid, r_vid));
        }
      }
    }
  }
  return log_weight;
}

val_id drawEntityValueSequential(
  attr_id aid, 
  const Records &recs, 
  const DistortProbs &distort_probs, 
  const DistortDistConcs &distort_dist_concs,
  const std::unordered_set<rec_id> &linked_rids, 
  const AbstractAttributeIndex *index, 
  IndexNonUniformDiscreteDist *distribution
) 
{
  val_id e_vid;
  double log_weight;
  double max_log_weight = R_NegInf;
#ifdef CHECK_WEIGHTS
  double total_weight = 0.0;
#endif
  
  std::vector<double> weights(index->domain_size());
  for (e_vid = 0; e_vid < index->domain_size(); e_vid++) 
  {
    log_weight = logWeightEntityValue(e_vid, aid, recs, distort_probs, distort_dist_concs, linked_rids, index, distribution);
#ifdef CHECK_WEIGHTS
    if (log_weight > 0 && !std::isfinite(log_weight)) { 
      std::string message = "invalid log-weight " + std::to_string(log_weight) + " detected when updating entity attribute " + std::to_string(aid);
      throw std::runtime_error(message); 
    }
#endif
    if (log_weight > max_log_weight) {
      max_log_weight = log_weight;
    }
    weights[e_vid] = log_weight;
  }
  for (auto &weight : weights) {
    weight = std::exp(weight - max_log_weight);
#ifdef CHECK_WEIGHTS
    if (weight < 0.0 || !std::isfinite(weight)) {
      std::string message = "invalid weight " + std::to_string(weight) + " detected when updating entity attribute " + std::to_string(aid);
      throw std::runtime_error(message);
    }
    total_weight += weight;
#endif
  }
#ifdef CHECK_WEIGHTS
  if (total_weight <= 0.0 || !std::isfinite(total_weight)) {
    std::string message = "no valid entity values for attribute " + std::to_string(aid);
    throw std::runtime_error(message);
  }
#endif
  IndexNonUniformDiscreteDist dist(weights.cbegin(), weights.cend(), false);
  return dist.draw();
}


val_id drawEntityValueIndex(
  attr_id aid, 
  const Records &recs, 
  const DistortProbs &distort_probs, 
  const DistortDistConcs &distort_dist_concs,
  const std::unordered_set<rec_id> &linked_rids, 
  const AbstractAttributeIndex *index, 
  IndexNonUniformDiscreteDist *distribution
) 
{
  // Exploit the following observation to avoid computing the pmf on the entire 
  // domain:
  //   pmf(v) > 0 iff one of the linked records has value w == v or 
  //   dist(v,w) <= threshold

  std::unordered_map<val_id, double> vid_weights;
  
  // Initialize vid_weights with necessary support (use 0.0 since working with log-probs)
  for (auto const &linked_rid : linked_rids) 
  {
    val_id r_vid = recs.get_attribute(linked_rid, aid);
    if (r_vid != NA_INTEGER) 
    {
      // Record attribute value is observed
      vid_weights[r_vid] = 0.0;
      
      const auto close_vids = index->get_close_values(r_vid);
      for (auto const &close_vid : close_vids) 
      {
        vid_weights[close_vid.first] = 0.0;
      }
    }
  }
  
  // If vid_weights is empty, then all linked record values are missing. Draw 
  // from entity distribution.
  if (vid_weights.empty()) {
    return distribution->draw();
  }
  
  double log_weight;
  double max_log_weight = R_NegInf;
#ifdef CHECK_WEIGHTS
  double total_weight = 0.0;
#endif
  for (auto &vw : vid_weights) 
  {
    log_weight = logWeightEntityValue(vw.first, aid, recs, distort_probs, distort_dist_concs, linked_rids, index, distribution);
#ifdef CHECK_WEIGHTS
    if (log_weight > 0 && !std::isfinite(log_weight)) { 
      std::string message = "invalid log-weight " + std::to_string(log_weight) + " detected when updating entity attribute " + std::to_string(aid);
      throw std::runtime_error(message); 
    }
#endif
    if (log_weight > max_log_weight) {
      max_log_weight = log_weight;
    }
    vw.second = log_weight;
  }
  for (auto &vw : vid_weights) {
    vw.second = std::exp(vw.second - max_log_weight);
#ifdef CHECK_WEIGHTS
    if (vw.second < 0.0 || !std::isfinite(vw.second)) {
      std::string message = "invalid weight " + std::to_string(vw.second) + " detected when updating entity attribute " + std::to_string(aid);
      throw std::runtime_error(message);
    }
    total_weight += vw.second;
#endif
  }
#ifdef CHECK_WEIGHTS
  if (total_weight <= 0.0 || !std::isfinite(total_weight)) { 
    std::string message = "no valid entity values for attribute " + std::to_string(aid);
    throw std::runtime_error(message);
  }
#endif

  if (vid_weights.empty()) return distribution->draw();

  NonUniformDiscreteDist<val_id> dist(vid_weights, false);
  return dist.draw();
}


void Entities::update_attributes(
  const Links &links, 
  const Records &recs, 
  const DistortProbs &distort_probs, 
  const DistortDistConcs &distort_dist_concs,
  const Cache &cache
) 
{
  IndexNonUniformDiscreteDist *distribution;
  for (auto const &pair : forward_index_) 
  { 
    const ent_id &eid = pair.first;
    const std::unordered_set<rec_id> &linked_rids = links.get_record_ids(eid);

    val_id e_vid;
    const AbstractAttributeIndex *index;
    for (attr_id aid=0; aid < recs.n_attributes(); aid++) 
    {
      index = cache.attr_indices_[aid].get();
      distribution = get_distribution(aid);
      if (typeid(*index) == typeid(ConstantAttributeIndex)) 
      {
        // TODO: replace sequential with method based on rejection sampling? 
        // Or some method that scales sublinearly in the size of the domain.
        e_vid = drawEntityValueSequential(aid, recs, distort_probs, distort_dist_concs, linked_rids, index, distribution);
      } 
      else 
      {
        //e_vid = drawEntityValueSequential(aid, recs, distort_probs, distort_dist_concs, linked_rids, index, distribution);
        e_vid = drawEntityValueIndex(aid, recs, distort_probs, distort_dist_concs, linked_rids, index, distribution);
      }
      
      set_attribute_value(eid, aid, e_vid);
    }
  }
}