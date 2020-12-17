// Walker Alias method for discrete non-uniform sampling
//
// Copyright (C) 2019  Neil Marchant
// Copyright (C) 2018  Australian Bureau of Statistics
// 
// Author: Neil Marchant
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#pragma once

#include <RcppArmadillo.h>
#include "types.h"
#include <iostream>


/**
 * Vose's alias method for generating discrete random variates. 
 */
class AliasSampler
{
private:
  /** Number of items in the discrete sample space */
  size_t n_items_; 
  
  /** Sum of the unnormalized weights */
  double total_weight_;
  
  /** The probability table. Note: this is not the pmf. */
  std::vector<double> probs_;
  
  /** The alias table */
  std::vector<size_t> alias_table_;
  
  /** 
   * Builds the tables assuming probs_ is initialized with the unnormalized 
   * weights 
   */
  void build_tables();
public:
  /**
   * Draw a discrete random variate
   * @return integer index corresponding to an item in the sample space
   */
  int draw() const;
  
  /**
   * @return sum of the unnormalized weights
   */
  double total_weight() const;
  
  /**
   * Construct an alias sampler from an iterator over unnormalized 
   * weights for items in the sample space. 
   * Items are identified by their position in the iterator.
   * @param begin,end input iterators to the initial and final positions 
   * in a range
   * @return an AliasSampler
   */
  template <class InputIterator>
  AliasSampler(InputIterator begin, InputIterator end) 
  {
    total_weight_ = 0;
    n_items_ = std::distance(begin, end);
    if (n_items_ == 0) {
        throw std::runtime_error("AliasSampler: weights must have non-zero length");
    }

    alias_table_.resize(n_items_);
    probs_.resize(n_items_);
    
    // Fill probs_ with weights while checking validity
    for (auto &p : probs_) 
    {
      p = *begin;
      if (p < 0 || !std::isfinite(p)) { 
        std::string message = "AliasSampler: invalid weight " + std::to_string(p) + " detected";
        throw std::runtime_error(message);
      }
      total_weight_ += p;
      ++begin;
    }
    if (total_weight_ <= 0 || !std::isfinite(total_weight_)) { 
      std::string message = "AliasSampler: total weight " + std::to_string(total_weight_) + " is invalid";
      throw std::runtime_error(message);
    }

    build_tables();
  }
};
