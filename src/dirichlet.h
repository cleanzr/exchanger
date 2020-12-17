#pragma once

#include <RcppArmadillo.h>

/**
 * Generate Dirichlet random variates
 * 
 * @param begin,end input iterators over the concentration parameters vector. 
 * All concentration parameters must be positive.
 * @return a pseudorandom probability vector
 */
template <class InputIterator>
std::vector<double> rdirichlet(InputIterator begin, InputIterator end) {
  size_t n_elem = std::distance(begin, end);
  std::vector<double> out(n_elem);
  double total = 0.0;

  // Draw independent gamma variables and compute sum
  for (auto &p : out) 
  {
    if (*begin <= 0.0) { throw std::range_error("rdirichlet: concentration parameter must be positive"); }
    p = R::rgamma(*begin, 1.0);
    total += p;
    ++begin;
  }
  
  // Normalizing results in draw from Dirichlet(alpha)
  for (auto &p : out) {
    p = p / total;
    if (p < 0.0) { Rcpp::stop("rdirichlet: generated negative probability mass"); }
  }

  return out;
}

