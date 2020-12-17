#include "discrete_dist.h"
#include <RcppArmadillo.h>

IndexNonUniformDiscreteDist::IndexNonUniformDiscreteDist(std::tuple<std::vector<double>, AliasSampler, double, double> members) 
  : probs_{std::get<0>(members)},
    alias_{std::move(std::get<1>(members))},
    min_prob_{std::get<2>(members)},
    max_prob_{std::get<3>(members)}
{}

double IndexNonUniformDiscreteDist::total_weight() const { return alias_.total_weight(); }

int IndexNonUniformDiscreteDist::num_items() const { return probs_.size(); }

double IndexNonUniformDiscreteDist::get_probability(int item) const { return probs_[item]; }

int IndexNonUniformDiscreteDist::draw() const { return alias_.draw(); }

double IndexNonUniformDiscreteDist::min_probability() const { return min_prob_; }

double IndexNonUniformDiscreteDist::max_probability() const { return max_prob_; }