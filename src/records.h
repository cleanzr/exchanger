#pragma once

#include "types.h"
#include <RcppArmadillo.h>
#include "entities.h"
#include "links.h"
#include "distortion_probs.h"
#include "cache.h"

class DistortProbs;
class Links;
class Cache;
class Entities;

/**
 * Represents the observed records
 */
class Records {
private:
  const arma::imat attr_vals_;
  const arma::ivec fids_;
  arma::imat attr_distortions_;
  const arma::ivec fsizes_;
  arma::imat distortion_counts_;

  /**
   * Setter for the distortion indicator of a record attribute
   * @param rid a record identifier
   * @param aid an attribute identifier
   * @param distorted the distortion indicator
   */
  void set_attribute_distortion(rec_id rid, attr_id aid, bool distorted);
public:
  /**
   * Get the total number of records
   */
  int n_records() const;
  
  /**
   * Get the number of records from a file/source
   * @param fid a file identifier
   * @return the number of records
   */
  int n_records(file_id fid) const;
  
  /**
   * Get the number of attributes
   */
  int n_attributes() const;
  
  /**
   * Get the number of files/sources
   */
  int n_files() const;
  
  /**
   * Getter for an attribute value of a record
   * @param rid a record identifier
   * @param aid an attribute identifier
   * @return the attribute value 
   */
  val_id get_attribute(rec_id rid, attr_id aid) const;
  
  /**
   * Getter for the file/source of a record
   * @param rid a record identifier
   * @return the file/source of the record
   */
  file_id get_file_id(rec_id rid) const;
  
  /**
   * Getter for the distortion indicator of a record attribute
   * @param rid a record identifier
   * @param aid an attribute identifier
   * @return a boolean distortion indicator
   */
  bool get_attribute_distortion(rec_id rid, attr_id aid) const;
  
  /**
   * Get the number of distorted values for a particular file and attribute
   * @param fid a file identifier
   * @param aid an attribute identifier
   * @return the number of aggregate distortions
   */
  int n_distorted(file_id fid, attr_id aid) const;
  
  /**
   * Perform a Gibbs update for the distortion indicators
   * @param ents a reference to the entities container
   * @param links a reference to the linkage structure container
   * @param distort_probs a reference to the distortion probabilities container
   * @param cache a reference to the cache
   */
  arma::imat update_distortion(const Entities& ents, const Links& links, const DistortProbs& distort_probs, const Cache& cache);
  
  /**
   * Convert the distortion indicators to an R representation
   * @return a distortion indicator matrix. Rows correspond to attributes and 
   * columns correspond to records.
   */
  Rcpp::IntegerMatrix R_rec_distortions() const;

  /**
   * Convert the aggregate number of distortions to an R representation 
   * @return a matrix containing the total number of distortions for each 
   * file/attribute. Rows correspond to files and columns correspond to 
   * attributes.
   */
  Rcpp::IntegerVector R_n_distorted_per_attr() const;

  /**
   * TODO
   */
  Rcpp::IntegerVector R_n_distorted_per_rec() const;
  
  /**
   * Constructor from R objects
   * @param attr_vals a matrix of record attribute values. Rows 
   * correspond to attributes and columns correspond to records.
   * @param fid a vector of file identifiers for each record.
   * @param attr_distortions a matrix of distortion indicators. 
   * Rows correspond to attribute and columns correspond to records.
   * @param fsizes a vector counting the number of records in each file.
   * @return a Records instance.
   */
  Records(const arma::imat &attr_vals, const arma::ivec &fids, 
          const arma::imat &attr_distortions, const arma::ivec &fsizes) 
    : attr_vals_(attr_vals), 
      fids_(fids), 
      attr_distortions_(attr_distortions), 
      fsizes_(fsizes)
  {
    distortion_counts_ = arma::zeros<arma::imat>(n_files(), n_attributes());
    for (rec_id rid=0; rid < n_records(); rid++) {
      file_id fid = fids_[rid];
      for (attr_id aid=0; aid < n_attributes(); aid++) {
        if (get_attribute_distortion(rid, aid)) { distortion_counts_(fid, aid) += 1; }
      }
    }
  }
};