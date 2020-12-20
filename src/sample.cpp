#include "state.h"
#include <progress.hpp>
#include <progress_bar.hpp>
#include <RcppArmadillo.h>
#include <algorithm>

/**
 * Convert a C++ representation of the model state to an `ExchangERModel` S4 object 
 * @param init_state reference to an `ExchangERModel` S4 object prior to 
 * sampling
 * @param state reference to a State instance.
 * @return an `ExchangERModel` S4 object
 */
Rcpp::S4 buildFinalS4State(const Rcpp::S4 &init_state, State &state) 
{
  Rcpp::S4 final_state("ExchangERModel");
  
  // Some slots are unchanged from init_state
  final_state.slot("attr_params") = init_state.slot("attr_params");
  final_state.slot("rec_attrs") = init_state.slot("rec_attrs");
  final_state.slot("file_ids") = init_state.slot("file_ids");
  final_state.slot("rec_ids") = init_state.slot("rec_ids");
  final_state.slot("attr_indices") = init_state.slot("attr_indices");
  final_state.slot("clust_prior") = init_state.slot("clust_prior");
  
  const Rcpp::IntegerMatrix &rec_attrs = final_state.slot("rec_attrs");
  const Rcpp::CharacterVector &rec_ids = final_state.slot("rec_ids");
  const Rcpp::IntegerVector &file_ids = init_state.slot("file_ids");
  
  // Other slots must be updated based on C++ State
  final_state.slot("iteration") = state.get_iteration();
  
  Rcpp::IntegerMatrix rec_distortions = state.recs_.R_rec_distortions();
  rownames(rec_distortions) = rownames(rec_attrs);
  colnames(rec_distortions) = rec_ids;
  final_state.slot("rec_distortions") = rec_distortions;
  
  std::pair<Rcpp::IntegerVector, Rcpp::IntegerMatrix> final_ents = state.ents_.to_R();
  Rcpp::IntegerMatrix &ent_attrs = final_ents.second;
  rownames(ent_attrs) = rownames(rec_attrs);
  final_state.slot("ent_attrs") = ent_attrs;
  
  final_state.slot("ent_ids") = final_ents.first;
  
  Rcpp::IntegerVector links = state.links_.to_R();
  links.attr("names") = rec_ids;
  final_state.slot("links") = links;
  
  Rcpp::NumericMatrix distort_probs = state.distort_probs_.to_R();
  colnames(distort_probs) = rownames(rec_attrs);
  const Rcpp::CharacterVector &fileNames = file_ids.attr("levels");
  rownames(distort_probs) = fileNames;
  final_state.slot("distort_probs") = distort_probs;
  
  Rcpp::S4 clust_params = state.clust_params_.get()->to_R();
  final_state.slot("clust_params") = clust_params;
  
  return final_state;
}

/**
 * Convert an `ExchangERModel` S4 object to a C++ representation of the model state
 * @param init_state reference to an `ExchangERModel` S4 object which 
 * represents the initial state of the model parameters
 * @return a State instance
 */
State readS4State(Rcpp::S4 init_state)
{
  int iteration = init_state.slot("iteration");
  const Rcpp::List &attr_params = init_state.slot("attr_params");
  const Rcpp::IntegerMatrix &rec_attrs_R = init_state.slot("rec_attrs");
  const arma::imat &rec_distortions = init_state.slot("rec_distortions");
  arma::imat ent_attrs = init_state.slot("ent_attrs");
  const Rcpp::IntegerVector &ent_ids = init_state.slot("ent_ids");
  const Rcpp::IntegerVector &file_ids_R = init_state.slot("file_ids");
  const Rcpp::IntegerVector &links = init_state.slot("links");
  const arma::mat &distort_probs = init_state.slot("distort_probs");
  const Rcpp::List &attr_indices = init_state.slot("attr_indices");
  const Rcpp::S4 &clust_prior = init_state.slot("clust_prior");
  const Rcpp::S4 &clust_params = init_state.slot("clust_params");
  
  arma::ivec file_sizes = Rcpp::as<arma::ivec>(Rcpp::table(file_ids_R));
  arma::ivec file_ids = Rcpp::as<arma::ivec>(file_ids_R);
  
  // Shift val_ids and file_ids to start at zero
  file_ids -= 1;
  // Do this using Rcpp to ensure NAs are preserved
  arma::imat rec_attrs = Rcpp::as<arma::imat>(rec_attrs_R - 1);
  ent_attrs -= 1;
  
  std::shared_ptr<ClustParams> clust_params_ = read_clust_params(clust_params, clust_prior);

  return State( iteration, 
    Entities(ent_ids, ent_attrs, attr_params, attr_indices),
    DistortProbs(distort_probs, attr_params),
    Links(links),
    Records(rec_attrs, file_ids, rec_distortions, file_sizes),
    Cache(attr_indices, attr_params),
    std::move(clust_params_) );
}

/**
 * Run Markov chain Monte Carlo on the ER model
 * 
 * @param init_state reference to an `ExchangERModel` S4 object which 
 * represents the initial state of the model parameters
 * @param n_samples number of posterior samples to generate after applying
 * thinning and burnin.
 * @param thin_interval period at which to save samples
 * @param burnin_interval number of iterations to discard as burnin
 * @return an `ExchangERResult` S4 object
 */
//[[Rcpp::export(.sample)]]
Rcpp::S4 sample(Rcpp::S4 init_state, int n_samples, int thin_interval=1, 
int burnin_interval=0) 
{
  #ifdef CHECK_WEIGHTS
  Rcpp::Rcout << "Running with CHECK_WEIGHTS" << std::endl;
  #endif
  State state = readS4State(init_state);
  
  // Convenience variables
  ClustParams *clust_params = state.clust_params_.get();
  
  int n_iter = burnin_interval + thin_interval * n_samples;
  Progress p(n_iter, true);
  
  // Preallocate R arrays to store history of variables along the chain
  Rcpp::CharacterVector rec_ids = init_state.slot("rec_ids");
  Rcpp::IntegerMatrix hist_links(n_samples, state.recs_.n_records());
  colnames(hist_links) = rec_ids;
  
  const Rcpp::IntegerMatrix &rec_attrs = init_state.slot("rec_attrs");
  Rcpp::CharacterVector attrs_names = rownames(rec_attrs);
  const Rcpp::IntegerVector& file_ids = init_state.slot("file_ids");
  Rcpp::CharacterVector file_ids_names = file_ids.attr("levels");
  Rcpp::NumericMatrix hist_distort_probs(n_samples, state.recs_.n_attributes() * state.recs_.n_files());
  Rcpp::CharacterVector hist_distort_probs_nms(state.recs_.n_attributes() * state.recs_.n_files());
  for (int i = 0; i < state.recs_.n_files(); ++i ) {
    for (int j = 0; j < state.recs_.n_attributes(); ++j) {
      std::string name = Rcpp::as<std::string>(attrs_names[j]) + "[" + Rcpp::as<std::string>(file_ids_names[i]) + "]"; 
      hist_distort_probs_nms[i + j * state.recs_.n_files()] = name;
    }
  }
  colnames(hist_distort_probs) = hist_distort_probs_nms;
  
  Rcpp::IntegerMatrix hist_n_linked_ents(n_samples, 1);
  colnames(hist_n_linked_ents) = Rcpp::CharacterVector::create("n_linked_ents");
  
  Rcpp::NumericMatrix hist_clust_params;
  if (clust_params->num_random() > 0) {
    Rcpp::CharacterVector clust_param_names = clust_params->to_R_vec().names();
    hist_clust_params = Rcpp::NumericMatrix(n_samples, clust_param_names.size());
    colnames(hist_clust_params) = clust_param_names;
  }
  
  int sample_ctr = 0;
  int initialIter = state.get_iteration();
  int completedIter;  // number of completed iterations (may be different from state's iteration counter if resuming)
  while (sample_ctr < n_samples) 
  {
    // Update state
    state.update();
    
    completedIter = state.get_iteration() - initialIter;

    if (completedIter >= burnin_interval) 
    {
      // Finished burn-in, so start saving samples
      if ((completedIter - burnin_interval) % thin_interval == 0)
      {
        // Update history using this sample
        Rcpp::IntegerVector links_R = state.links_.to_R();
        Rcpp::IntegerVector n_distorted = state.recs_.R_n_distorted_per_attr();
        hist_links.row(sample_ctr) = links_R;
        hist_n_linked_ents(sample_ctr, 0) = state.links_.n_linked_ents();
        Rcpp::NumericVector distort_prob = state.distort_probs_.to_R();
        hist_distort_probs.row(sample_ctr) = distort_prob;
        if (hist_clust_params.ncol() > 0) {
          hist_clust_params.row(sample_ctr) = clust_params->to_R_vec();
        }
        sample_ctr++;
      }
    }
    
    if (Progress::check_abort()) { break; }
    p.increment();
  }
  
  // Build final S4 state
  Rcpp::S4 final_state = buildFinalS4State(init_state, state);

  // Return output as a list
  Rcpp::S4 result("ExchangERResult");
  Rcpp::List history;
  history["links"] = hist_links;
  history["distort_probs"] = hist_distort_probs;
  history["n_linked_ents"] = hist_n_linked_ents;
  if (hist_clust_params.ncol() > 0) history["clust_params"] = hist_clust_params;
  result.slot("history") = history;
  result.slot("state") = final_state;
  return result;
}