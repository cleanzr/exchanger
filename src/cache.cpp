#include "cache.h"
#include <RcppArmadillo.h>

std::vector<std::shared_ptr<AbstractAttributeIndex> > Cache::init_members(const Rcpp::List &attr_indices, 
  const Rcpp::List &attr_params) 
{
  std::vector<std::shared_ptr<AbstractAttributeIndex> > attr_indices_;
  attr_indices_.reserve(attr_indices.size());
  int attrId = 0;
  for (const Rcpp::S4 &R_index : attr_indices) 
  {
    const Rcpp::S4 &spec = attr_params[attrId];
    const Rcpp::S4 &distort_dist_prior = spec.slot("distort_dist_prior");
    double dirichlet_concentration = distort_dist_prior.slot("alpha");
    bool exclude_entity_value = spec.slot("exclude_entity_value");
    attr_indices_.push_back(readAttributeIndex(R_index, dirichlet_concentration, exclude_entity_value));
    attrId++;
  }
  return attr_indices_;
}

void Cache::update(const Entities &ents) 
{
  attr_id aid = 0;
  for (auto const &index : attr_indices_) {
    IndexNonUniformDiscreteDist *distribution = ents.get_distribution(aid);
    index->update(distribution);
    aid++;
  }
}