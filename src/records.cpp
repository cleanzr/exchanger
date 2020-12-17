#include "records.h"
#include <RcppArmadillo.h>

int Records::n_records() const { return attr_vals_.n_cols; }

int Records::n_records(file_id fid) const { return fsizes_[fid]; }

int Records::n_attributes() const { return attr_vals_.n_rows; }

int Records::n_files() const { return fsizes_.size(); }

val_id Records::get_attribute(rec_id rid, attr_id aid) const 
{
  return attr_vals_(aid, rid);
}

int Records::get_file_id(rec_id rid) const { return fids_[rid]; }

bool Records::get_attribute_distortion(rec_id rid, attr_id aid) const 
{
  return attr_distortions_(aid, rid);
}

int Records::n_distorted(file_id fid, attr_id aid) const 
{
  return distortion_counts_(fid, aid);
}

void Records::set_attribute_distortion(rec_id rid, attr_id aid, bool distorted) 
{
  int &oldDist = attr_distortions_(aid, rid);
  file_id fid = get_file_id(rid);
  distortion_counts_(fid, aid) += distorted - oldDist;
  oldDist = distorted;
}

Rcpp::IntegerMatrix Records::R_rec_distortions() const 
{
  return Rcpp::wrap(attr_distortions_);
}

Rcpp::IntegerVector Records::R_n_distorted_per_attr() const 
{
  // Sum over files
  arma::irowvec n_distorted_per_attr = arma::sum(distortion_counts_, 0);
  return Rcpp::wrap(n_distorted_per_attr);
}

Rcpp::IntegerVector Records::R_n_distorted_per_rec() const 
{
  // Sum over files
  arma::ivec rec_distort_counts = arma::zeros<arma::ivec>(n_attributes() + 1);
  
  int n_distorted_attr; // number of distorted attributes for a given record
  
  // Loop over records (columns)
  for (size_t c = 0; c < attr_distortions_.n_cols; c++)
  {
    auto it_end = attr_distortions_.end_col(c);
    n_distorted_attr = 0;
    // Loop over attributes for record, counting number of distortions
    for (auto it = attr_distortions_.begin_col(c); it != it_end; ++it)
    {
      n_distorted_attr += (*it);
    }
    rec_distort_counts(n_distorted_attr) += 1;
  }
  return Rcpp::wrap(rec_distort_counts);
}


arma::imat Records::update_distortion(const Entities& ents, const Links& links, 
  const DistortProbs& distort_probs, const Cache& cache) 
{ 
  file_id fid;
  ent_id linked_eid;
  val_id r_vid, e_vid;
  bool newDistortion;
  bool distortedAlternative;
  double distort_prob;
  double pr1, pr0;
  const AbstractAttributeIndex* index;
  arma::imat corr_nondistort_counts = arma::zeros<arma::imat>(this->n_files(), this->n_attributes());

  for (rec_id rid=0; rid < this->n_records(); rid++)
  {
    fid = this->get_file_id(rid);
    for (attr_id aid=0; aid < this->n_attributes(); aid++) 
    {
      linked_eid = links.get_entity_id(rid);
      r_vid = this->get_attribute(rid, aid);
      e_vid = ents.get_attribute_value(linked_eid, aid);
      distort_prob = distort_probs.get(fid, aid);
      index = cache.attr_indices_[aid].get();
      
      // First update the distortion indicator for the attribute. 
      if (r_vid == NA_INTEGER) 
      {
        // If the record attribute value is unobserved, draw indicator from 
        // the prior.
        newDistortion = (R::unif_rand() < distort_prob);
      } 
      else 
      {
        // Indicator is deterministic conditional on the observed entity and 
        // record attribute values.
        if (r_vid == e_vid) {
          if (index->exclude_entity_value()) {
            newDistortion = false;
          } else {
            // Record value can take on entity value, while being distorted
            double pr1 = distort_prob * index->get_distortion_prob(e_vid, r_vid);
            double pr0 = 1.0 - distort_prob;
            newDistortion = (R::unif_rand() < (pr1 / (pr1 + pr0)));
          }
        } else {
          newDistortion = true;
        }
      }
      this->set_attribute_distortion(rid, aid, newDistortion);

      // Next update an auxiliary variable which is used to update the 
      // distortion probabilities
      pr1 = index->get_max_exp_factor(e_vid);
      if (pr1 == 1.0) {
        distortedAlternative = true;
      } else {
        pr0 = 1.0 - pr1;
        if (newDistortion) {
          pr1 *= distort_prob;
        } else {
          pr1 *= 1.0 - distort_prob;
        }
        distortedAlternative = (R::unif_rand() < (pr1 / (pr1 + pr0)));
      }
      corr_nondistort_counts(fid, aid) += distortedAlternative & !newDistortion;
    }
  }
  return corr_nondistort_counts;
}

