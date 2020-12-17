#include "distortion_probs.h"
#include <RcppArmadillo.h>

std::vector<BetaRV> DistortProbs::init_prior(const Rcpp::List &attr_params) {
  std::vector<BetaRV> distort_prob_priors;
  distort_prob_priors.reserve(attr_params.size());
  for (const Rcpp::S4 &spec : attr_params) 
  {
    const Rcpp::S4 &prior = spec.slot("distort_prob_prior");
    distort_prob_priors.push_back(BetaRV(prior));
  }
  return distort_prob_priors;
}

double DistortProbs::get(file_id fid, attr_id aid) const {
  return probs_(fid, aid);
}

Rcpp::NumericMatrix DistortProbs::to_R() const {
  return Rcpp::wrap(probs_);
}

void DistortProbs::update(const Records& recs, const arma::imat &corr_nondistort_counts) 
{
  double n_distorted;
  double eff_n_distorted, eff_n_nondistorted, prob;
  for (file_id fid = 0; fid < recs.n_files(); fid++) 
  {
    for (attr_id aid = 0; aid < recs.n_attributes(); aid++) 
    {
      const BetaRV &prior = distort_prob_priors_[aid];
      n_distorted = recs.n_distorted(fid, aid);
      eff_n_distorted = n_distorted + prior.get_shape1();
      eff_n_nondistorted = corr_nondistort_counts(fid, aid) + prior.get_shape2();
      prob = R::rbeta(eff_n_distorted, eff_n_nondistorted);
      probs_(fid, aid) = prob;
    }
  }
}