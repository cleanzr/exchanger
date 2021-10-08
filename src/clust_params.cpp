#include "clust_params.h"
#include "dirichlet.h"
#include <RcppArmadillo.h>
#include <map>

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

double ESCNBClustParams::prior_weight_existing(int cluster_size) const
{ 
  return std::abs(cluster_size + size_);
}

double ESCNBClustParams::prior_weight_new(int n_clusters) const { 
  return (n_clusters + 1) * size_ * std::pow(1 - prob_, size_); 
}

double ESCDClustParams::prior_weight_existing(int cluster_size) const
{ 
  return (cluster_size + 1) * mu_[cluster_size] / mu_[cluster_size - 1];
}

double ESCDClustParams::prior_weight_new(int n_clusters) const { 
  return (n_clusters + 1) * mu_[0]; 
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

void ESCNBClustParams::update(Links &links)
{
  auto n_clusters = links.n_linked_ents();

  // Update prob
  if (p_prior_.has_value()) {
    double shape1 = n_clusters * size_ + p_prior_.value().get_shape1();
    double shape2 = links.n_records() + p_prior_.value().get_shape2();
    prob_ = R::rbeta(shape1, shape2);
  }

  if (r_prior_.has_value()) {
    // Generate v auxiliary variables
    int clust_size;
    double prob_j;
    int num_v_success = 0;
    for (auto const &clust : links.inverse_index_) {
      clust_size = clust.second.size();
      for (int j = 1; j <= clust_size - 2; j++) {
        prob_j = (size_)/(size_ + j);
        num_v_success += (R::unif_rand() < prob_j);
      }
    }

    // Update size
    double shape = r_prior_.value().get_shape() + n_clusters + num_v_success;
    double rate = r_prior_.value().get_rate() - n_clusters * std::log(prob_);
    size_ = R::rgamma(shape, 1.0 / rate);
  }
}

void ESCDClustParams::update(Links &links)
{ 
  // Frequency of each cluster size
  std::map<int, int> clust_size_freq;
  int clust_size;
  int max_clust_size = 0;
  for (auto const &clust : links.inverse_index_) {
    clust_size = clust.second.size();
    auto got = clust_size_freq.find(clust_size);
    if (got != clust_size_freq.end()) {
      got->second += 1;
    } else {
      clust_size_freq.insert(std::make_pair(clust_size, 1));
    }
    
    if (clust_size > max_clust_size)
      max_clust_size = clust_size;
  }
  
  if (p_prior_.has_value() || r_prior_.has_value()) {
    double prod_si;  // \prod_{j}^{s-1} size^{u_{sij}} * (j - 1)^{1 - u_{sij}} / j
    int u_si;        // \sum_{j}^{s-1} u_{sij}
    double z_si;     // alpha * prob^size * (1 - prob)^{s-1} * prod_si
    double prob_si;  // z_{si} / (z_{si} + i - 1)
    double prob_j;   // size / (size + j - 1)
    bool v_si;       // auxiliary variable drawn from Bernoulli(prob_si)

    // Store delta updates to prob_ and size_ hyperparameters
    double shape_prime = 0.0;
    double rate_prime = 0.0;
    double shape1_prime = 0.0;
    double shape2_prime = 0.0;

    for (auto const &x : clust_size_freq) {
      const int &sz = x.first;
      const int &count = x.second;
      for (int i = 1; i <= count; i++) {
        prod_si = 1.0;
        u_si = 0;
        if (r_prior_.has_value()) {
          // Generate u auxiliary variables
          bool u_sij;
          for (int j = 1; j <= sz - 1; j++) {
            prob_j = (size_)/(size_ + j - 1);
            u_sij = (R::unif_rand() < prob_j);
            u_si += u_sij;
            prod_si *= size_/j * u_sij + (1 - 1/j) * (1 - u_sij);
          }
        }
        // Generate v auxiliary variable
        z_si = alpha_ * std::pow(prob_, size_) * std::pow(1 - prob_, sz - 1) * prod_si;
        prob_si = z_si / (z_si + i - 1);
        v_si = (R::unif_rand() < prob_si);

        // Update prob_ and size_ hyperparameters
        shape_prime += v_si * u_si;
        rate_prime += -v_si;
        shape1_prime += v_si;
        shape2_prime += (sz - 1) * v_si;
      }
    }
    
    rate_prime *= std::log(prob_);
    shape1_prime *= size_;
    
    // Update prob
    if (p_prior_.has_value()) {
      double shape1 = p_prior_.value().get_shape1() + shape1_prime;
      double shape2 = p_prior_.value().get_shape2() + shape2_prime;
      prob_ = R::rbeta(shape1, shape2);
    }
    
    // Update size
    if (r_prior_.has_value()) {
      double shape = r_prior_.value().get_shape() + shape_prime;
      double rate = r_prior_.value().get_rate() + rate_prime;
      size_ = R::rgamma(shape, 1.0 / rate);
    }
  }
  
  // Update mu
  std::vector<double> c_vec(max_clust_size + 1);  // concentration vector
  clust_size = 1;
  double c_prior_sum = 0.0;
  for (auto &x : c_vec) {
    // Prior contribution
    if (clust_size > max_clust_size) {
      // Element corresponding to remaining cluster sizes is computed differently
      x = alpha_ - c_prior_sum;
      break; // No contribution from observations
    } else {
      // alpha times negative binomial pmf
      x = alpha_ * R::dnbinom(clust_size, size_, prob_, 0);
      c_prior_sum += x;
    }

    // Contribution from observations
    auto got = clust_size_freq.find(clust_size);
    if (got != clust_size_freq.end()) {
      x += got->second;
    }

    clust_size++;
  }
  mu_ = rdirichlet(c_vec.begin(), c_vec.end());
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

int ESCNBClustParams::num_random() const {
  return p_prior_.has_value() + r_prior_.has_value();
}

int ESCDClustParams::num_random() const {
  return p_prior_.has_value() + r_prior_.has_value();
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

Rcpp::S4 ESCNBClustParams::to_R() const {
  Rcpp::S4 out("ESCNBRP");
  out.slot("size") = size_;
  out.slot("prob") = prob_;
  return out;
}

Rcpp::S4 ESCDClustParams::to_R() const {
  Rcpp::S4 out("ESCDRP");
  out.slot("size") = size_;
  out.slot("prob") = prob_;
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

Rcpp::NumericVector ESCNBClustParams::to_R_vec(bool include_fixed) const 
{
  std::vector<std::string> names;
  std::vector<double> values;
  if (p_prior_.has_value() || include_fixed) {
    names.push_back("prob");
    values.push_back(prob_);
  }
  if (r_prior_.has_value() || include_fixed) {
    names.push_back("size");
    values.push_back(size_);
  }

  Rcpp::NumericVector out(values.begin(), values.end());
  Rcpp::CharacterVector out_names(names.begin(), names.end());
  out.attr("names") = out_names;
  return out; 
}

Rcpp::NumericVector ESCDClustParams::to_R_vec(bool include_fixed) const 
{
  std::vector<std::string> names;
  std::vector<double> values;
  if (p_prior_.has_value() || include_fixed) {
    names.push_back("prob");
    values.push_back(prob_);
  }
  if (r_prior_.has_value() || include_fixed) {
    names.push_back("size");
    values.push_back(size_);
  }
  if (include_fixed) {
    names.push_back("alpha");
    values.push_back(alpha_);
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

ESCNBClustParams::ESCNBClustParams(const Rcpp::S4 &clust_params, const Rcpp::S4 &clust_prior) 
{
  // Read prior
  SEXP p_SEXP = clust_prior.slot("prob");
  if (TYPEOF(p_SEXP) == S4SXP) {
    Rcpp::S4 p_S4 = Rcpp::as<Rcpp::S4>(p_SEXP);
    p_prior_ = BetaRV(p_S4);
  }

  SEXP r_SEXP = clust_prior.slot("size");
  if (TYPEOF(r_SEXP) == S4SXP) {
    Rcpp::S4 r_S4 = Rcpp::as<Rcpp::S4>(r_SEXP);
    r_prior_ = GammaRV(r_S4);
  }

  // Read parameter values
  prob_ = clust_params.slot("prob");
  size_ = clust_params.slot("size");
}

ESCDClustParams::ESCDClustParams(const Rcpp::S4 &clust_params, const Rcpp::S4 &clust_prior) 
{
  // Read prior
  SEXP p_SEXP = clust_prior.slot("prob");
  if (TYPEOF(p_SEXP) == S4SXP) {
    Rcpp::S4 p_S4 = Rcpp::as<Rcpp::S4>(p_SEXP);
    p_prior_ = BetaRV(p_S4);
  }
  
  SEXP r_SEXP = clust_prior.slot("size");
  if (TYPEOF(r_SEXP) == S4SXP) {
    Rcpp::S4 r_S4 = Rcpp::as<Rcpp::S4>(r_SEXP);
    r_prior_ = GammaRV(r_S4);
  }
  
  // Read parameter values
  prob_ = clust_params.slot("prob");
  size_ = clust_params.slot("size");
  alpha_ = clust_params.slot("alpha");
}

std::shared_ptr<ClustParams> read_clust_params(const Rcpp::S4 &clust_params, const Rcpp::S4 &clust_prior) 
{
  std::shared_ptr<ClustParams> clust_params_;
  if (clust_prior.is("PitmanYorRP")) {
    clust_params_.reset(new PitmanYorClustParams(clust_params, clust_prior));
  } else if (clust_prior.is("ESCNBRP")) {
    clust_params_.reset(new ESCNBClustParams(clust_params, clust_prior));
  } else if (clust_prior.is("ESCDRP")) {
    clust_params_.reset(new ESCDClustParams(clust_params, clust_prior));
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