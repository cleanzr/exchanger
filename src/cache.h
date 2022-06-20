#pragma once

#include "types.h"
#include "attribute_index.h"
#include "entities.h"
#include <RcppArmadillo.h>
#include <memory>

class Entities;

class Cache {
private:
  static std::vector<std::shared_ptr<AbstractAttributeIndex> > init_members(const Rcpp::List &attr_indices, const Rcpp::List &attr_params);
public:
  const std::vector<std::shared_ptr<AbstractAttributeIndex> > attr_indices_;
public:
  void update(const Entities &ents);
  Cache(std::vector<std::shared_ptr<AbstractAttributeIndex> > attr_indices)
    : attr_indices_(std::move(attr_indices)) 
  {}
  Cache(const Rcpp::List &attr_indices, const Rcpp::List &attr_params) 
    : Cache{init_members(attr_indices, attr_params)} 
  {}
};