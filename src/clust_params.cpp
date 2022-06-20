#include "clust_params.h"
#include <RcppArmadillo.h>

double CouponClustParams::prior_weight_existing(int cluster_size) const
{ 
  return 1.0; 
}

double CouponClustParams::prior_weight_new(int n_clusters) const
{ 
  return std::abs(m_  - n_clusters); 
}

double GenCouponClustParams::prior_weight_existing(int cluster_size) const
{ 
  return cluster_size + kappa_; 
}

double GenCouponClustParams::prior_weight_new(int n_clusters) const
{ 
  return kappa_ * std::abs(m_  - n_clusters); 
}

double PitmanYorClustParams::prior_weight_existing(int cluster_size) const
{ 
  return std::abs(cluster_size - d_);
}

double PitmanYorClustParams::prior_weight_new(int n_clusters) const { 
  return alpha_ + d_ * n_clusters; 
}

void CouponClustParams::update(Links &links) {
  // Nothing to do
}

void GenCouponClustParams::update(Links &links) {
    
  if (!kappa_prior_.has_value() && !m_prior_.has_value()) {
    // Nothing to update
    return;
  }

  // Generate w auxiliary variable
  double w = R::rbeta(m_ * kappa_ + 1, links.n_records() - 1);

  if (kappa_prior_.has_value()) {
    // Generate v auxiliary variables
    int clust_size;
    double prob_j;
    int num_v_success = 0;
    for (auto const &clust : links.inverse_index_) {
      clust_size = clust.second.size();
      for (int j = 1; j <= clust_size - 1; j++) {
        prob_j = (kappa_)/(kappa_ + j);
        num_v_success += (R::unif_rand() < prob_j);
      }
    }

    // Update kappa
    double shape = kappa_prior_.value().get_shape() + links.n_linked_ents() - 1 + num_v_success;
    double rate = kappa_prior_.value().get_rate() - m_ * std::log(w);
    kappa_ = R::rgamma(shape, 1.0 / rate);
  }
  
  if (m_prior_.has_value()) {
    // Update m
    double sz = m_prior_.value().get_size() + links.n_linked_ents() - 1;
    double pb = 1 - (1 - m_prior_.value().get_prob()) * std::pow(w, kappa_);
    m_ = R::rnbinom(sz, pb) + links.n_linked_ents();
  }
}

void PitmanYorClustParams::update(Links &links)
{
  if (!d_prior_.has_value() || !alpha_prior_.has_value()) {
    // Nothing to update
    return;
  }
  
  // Generate u auxiliary variables
  bool u_i;
  double prob_i;
  int num_u_success = 0;
  int num_u_fail = 0;
  for (int i = 1; i <= links.n_linked_ents() - 1; i++) {
    prob_i = alpha_/(alpha_ + d_ * i);
    u_i = (R::unif_rand() < prob_i);
    num_u_success += u_i;
    num_u_fail += 1 - u_i;
  }

  if (alpha_prior_.has_value()) {
    // Generate w auxiliary variable
    double w = R::rbeta(alpha_ + 1, links.n_records() - 1);

    // Update alpha
    double shape = alpha_prior_.value().get_shape() + num_u_success;
    double rate = alpha_prior_.value().get_rate() - std::log(w);
    alpha_ = R::rgamma(shape, 1.0 / rate);
  }

  if (d_prior_.has_value()) {
    // Generate v auxiliary variables
    int clust_size;
    double prob_j;
    int num_v_fail = 0;
    for (auto const &clust : links.inverse_index_) {
      clust_size = clust.second.size();
      for (int j = 1; j <= clust_size - 1; j++) {
        prob_j = (j - 1)/(j - d_);
        num_v_fail += (R::unif_rand() >= prob_j);
      }
    }

    // Update d
    double shape1 = d_prior_.value().get_shape1() + num_u_fail;
    double shape2 = d_prior_.value().get_shape2() + num_v_fail;
    d_ = R::rbeta(shape1, shape2);
  }
}

int CouponClustParams::num_random() const {
  return 0;
}

int GenCouponClustParams::num_random() const {
  return m_prior_.has_value() + kappa_prior_.has_value();
}

int PitmanYorClustParams::num_random() const {
  return alpha_prior_.has_value() + d_prior_.has_value();
}

Rcpp::S4 CouponClustParams::to_R() const
{
  Rcpp::S4 out("GeneralizedCouponRP");
  out.slot("m") = m_;
  out.slot("kappa") = R_PosInf;
  return out;
}

Rcpp::S4 GenCouponClustParams::to_R() const
{
  Rcpp::S4 out("GeneralizedCouponRP");
  out.slot("m") = m_;
  out.slot("kappa") = kappa_;
  return out;
}

Rcpp::S4 PitmanYorClustParams::to_R() const {
  Rcpp::S4 out("PitmanYorRP");
  out.slot("alpha") = alpha_;
  out.slot("d") = d_;
  return out;
}

Rcpp::NumericVector CouponClustParams::to_R_vec(bool include_fixed) const 
{
  std::vector<std::string> names;
  std::vector<double> values;
  if (include_fixed) {
    names.push_back("kappa");
    values.push_back(R_PosInf);
    names.push_back("m");
    values.push_back(m_);
  }
  
  Rcpp::NumericVector out(values.begin(), values.end());
  Rcpp::CharacterVector out_names(names.begin(), names.end());
  out.attr("names") = out_names;
  return out; 
}

Rcpp::NumericVector GenCouponClustParams::to_R_vec(bool include_fixed) const 
{
  std::vector<std::string> names;
  std::vector<double> values;
  if (kappa_prior_.has_value() || include_fixed) {
    names.push_back("kappa");
    values.push_back(kappa_);
  }
  if (m_prior_.has_value() || include_fixed) {
    names.push_back("m");
    values.push_back(m_);
  }

  Rcpp::NumericVector out(values.begin(), values.end());
  Rcpp::CharacterVector out_names(names.begin(), names.end());
  out.attr("names") = out_names;
  return out; 
}

Rcpp::NumericVector PitmanYorClustParams::to_R_vec(bool include_fixed) const 
{
  std::vector<std::string> names;
  std::vector<double> values;
  if (alpha_prior_.has_value() || include_fixed) {
    names.push_back("alpha");
    values.push_back(alpha_);
  }
  if (d_prior_.has_value() || include_fixed) {
    names.push_back("d");
    values.push_back(d_);
  }

  Rcpp::NumericVector out(values.begin(), values.end());
  Rcpp::CharacterVector out_names(names.begin(), names.end());
  out.attr("names") = out_names;
  return out; 
}

CouponClustParams::CouponClustParams(const Rcpp::S4 &clust_prior) 
{
  // Read parameter values
  m_ = clust_prior.slot("m");
}

PitmanYorClustParams::PitmanYorClustParams(const Rcpp::S4 &clust_params, const Rcpp::S4 &clust_prior) 
{
  // Read prior
  SEXP alpha_SEXP = clust_prior.slot("alpha");
  if (TYPEOF(alpha_SEXP) == S4SXP) {
    Rcpp::S4 alpha_S4 = Rcpp::as<Rcpp::S4>(alpha_SEXP);
    alpha_prior_ = GammaRV(alpha_S4);
  }

  SEXP d_SEXP = clust_prior.slot("d");
  if (TYPEOF(d_SEXP) == S4SXP) {
    Rcpp::S4 d_S4 = Rcpp::as<Rcpp::S4>(d_SEXP);
    d_prior_ = BetaRV(d_S4);
  }

  // Read parameter values
  alpha_ = clust_params.slot("alpha");
  d_ = clust_params.slot("d");
}

GenCouponClustParams::GenCouponClustParams(const Rcpp::S4 &clust_params, const Rcpp::S4 &clust_prior) 
{
  // Read prior
  SEXP m_SEXP = clust_prior.slot("m");
  if (TYPEOF(m_SEXP) == S4SXP) {
    Rcpp::S4 m_S4 = Rcpp::as<Rcpp::S4>(m_SEXP);
    m_prior_ = ShiftedNegBinomRV(m_S4);
  }

  SEXP kappa_SEXP = clust_prior.slot("kappa");
  if (TYPEOF(kappa_SEXP) == S4SXP) {
    Rcpp::S4 kappa_S4 = Rcpp::as<Rcpp::S4>(kappa_SEXP);
    kappa_prior_ = GammaRV(kappa_S4);
  }

  // Read parameter values
  m_ = clust_params.slot("m");
  kappa_ = clust_params.slot("kappa");
}

std::shared_ptr<ClustParams> read_clust_params(const Rcpp::S4 &clust_params, const Rcpp::S4 &clust_prior) 
{
  std::shared_ptr<ClustParams> clust_params_;
  if (clust_prior.is("PitmanYorRP")) {
    clust_params_.reset(new PitmanYorClustParams(clust_params, clust_prior));
  } else if (clust_prior.is("GeneralizedCouponRP")) {
    // Sampling implementation differs depending on whether kappa is Inf
    SEXP kappa_SEXP = clust_prior.slot("kappa");
    if (TYPEOF(kappa_SEXP) == REALSXP) {
      Rcpp::NumericVector kappa_vec = Rcpp::as<Rcpp::NumericVector>(kappa_SEXP);
      if (kappa_vec[0] == R_PosInf) {
        clust_params_.reset(new CouponClustParams(clust_prior));  
      } else {
        clust_params_.reset(new GenCouponClustParams(clust_params, clust_prior));
      }
    } else {
      clust_params_.reset(new GenCouponClustParams(clust_params, clust_prior));  
    }
  } else {
    Rcpp::stop("Unrecognized clust_prior");
  }
  return clust_params_;
}
