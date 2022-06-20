#include "distort_dist_concs.h"

std::vector<boost::optional<GammaRV>> DistortDistConcs::init_prior(const Rcpp::List &attr_params) {
  
  std::vector<boost::optional<GammaRV>> priors;
  priors.reserve(attr_params.size());
  for (const Rcpp::S4 &spec : attr_params) 
  {
    Rcpp::S4 distort_dist_prior = spec.slot("distort_dist_prior");
    SEXP alpha_SEXP = distort_dist_prior.slot("alpha");
    if (TYPEOF(alpha_SEXP) == S4SXP) {
      Rcpp::S4 alpha_S4 = Rcpp::as<Rcpp::S4>(alpha_SEXP);
      GammaRV prior_(alpha_S4);
      priors.push_back(prior_);
    } else {
      priors.push_back(boost::none);
    }
  }
  return priors;
}

double DistortDistConcs::get(attr_id aid) const {
  return values_(aid);
}

Rcpp::NumericVector DistortDistConcs::to_R() const {
  return Rcpp::wrap(values_);
}

void DistortDistConcs::update(
  const Links &links, 
  const Records &recs, 
  const Entities &ents,
  const Cache &cache
)
{
  const AbstractAttributeIndex* index;
  for (attr_id aid = 0; aid < recs.n_attributes(); aid++) 
  {
    const boost::optional<GammaRV> &prior_ = priors_[aid];

    if (prior_.has_value()) {
      index = cache.attr_indices_[aid].get();
      // Update auxiliary variables first
      ent_id eid;
      val_id r_vid, e_vid;
      int ctr_disagree;
      double sum_zeta_ei = 0.0;
      double sum_log_t_e = 0.0;
      double prob_ei;
      std::unordered_map<val_id, int> ctr_disagree_value;
      
      //map vid to counts
      for (auto const &clust : links.inverse_index_) 
      {
        ctr_disagree = 0;
        
        eid = clust.first;
        e_vid = ents.get_attribute_value(eid, aid);
        for (auto const &rid : clust.second)
        {
          r_vid = recs.get_attribute(rid, aid);
          if (r_vid != NA_INTEGER) 
          {
            if (r_vid != e_vid) 
            {
              auto got = ctr_disagree_value.find(r_vid);
              if (got != ctr_disagree_value.end()) {
              // Seen this r_vid before
              double z = values_[aid] * index->get_distortion_prob(e_vid, r_vid);
              int &i = got->second;
              i++;
              prob_ei = z / (z + i - 1);
              } else {
              // Seeing this r_vid for the first time
              ctr_disagree_value[r_vid] = 1;
              prob_ei = 1.0;
              }
              ctr_disagree += 1;
              sum_zeta_ei += (R::unif_rand() < prob_ei);
            }
          }
        }
      
        sum_log_t_e += std::log(R::rbeta(values_[aid], ctr_disagree));
      }
      
      double shape = prior_.value().get_shape() + sum_zeta_ei;
      double rate = prior_.value().get_rate() - sum_log_t_e;
      values_[aid] = R::rgamma(shape, 1.0 / rate);
    }
  }
}