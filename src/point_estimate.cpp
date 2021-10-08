#include <Rcpp.h>
using namespace Rcpp;
#include <unordered_map>

namespace std {
template <>
struct hash<set<int>> {
  size_t operator()(const set<int>& s) const
  {
    size_t hash=0;
    for (auto const item : s) {
      hash += item;
    }
    return hash;
  }
};
}

template<class K, class V>
using map = std::unordered_map<K, V>;

template<class T>
using set = std::set<T>;

// Build an inverted index from entityIds to recordIds
template <class InputIterator>
map<int, set<int> > membership_to_clusters(InputIterator begin, InputIterator end)
{
  map<int, set<int> > out;
  
  int rec_id = 1; // start rec_id at 1 for compatibility with R indexing
  for (; begin != end; ++begin)
  {
    const int &clust_id = *begin;
    if (clust_id != NA_INTEGER) out[clust_id].insert(rec_id);
    rec_id++;
  }
  
  return out;
}

Rcpp::List membership_to_clusters(const Rcpp::IntegerVector &membership)
{
  map<int, set<int> > clusters = membership_to_clusters(membership.cbegin(), membership.cend());
  
  // Convert to List of IntegerVectors: one for each cluster
  Rcpp::List out(clusters.size());
  size_t i = 0;
  for (auto const &cluster : clusters) {
    Rcpp::IntegerVector c = Rcpp::wrap(cluster.second);
    out[i] = c;
    i++;
  }
  
  return out;
}

//[[Rcpp::export(.mp_clusters)]]
Rcpp::List mp_clusters(const Rcpp::IntegerMatrix &samples) {
  // Convenience variables
  int n_records = samples.ncol();
  int n_samples = samples.nrow();
  
  Rcpp::List out(n_records);
  map<set<int>, int> cluster_freqs;
  
  // Each row represents a sample clustering as a membership vector
  // Iterate through samples, counting frequency of each cluster
  for (int r = 0; r < n_samples; r++) {
    Rcpp::IntegerMatrix::ConstRow memb = samples.row(r);
    auto clusters = membership_to_clusters(memb.cbegin(), memb.cend());
    
    for (auto const &cluster : clusters) { cluster_freqs[cluster.second]++; }
  }
  
  // Find most probable clusters (pointer to cluster_freqs) for each record
  std::vector<map<set<int>, int>::iterator > most_probable(n_records, cluster_freqs.end());
  
  map<set<int>, int>::iterator it_cluster_freqs = cluster_freqs.begin();
  for (; it_cluster_freqs != cluster_freqs.end(); ++it_cluster_freqs)
  {
    for (auto const r : it_cluster_freqs->first) {
      map<set<int>, int>::iterator &it_most_probable = most_probable[r - 1];
      if (it_most_probable == cluster_freqs.end() ||
          (it_cluster_freqs->second > it_most_probable->second) ) {
        it_most_probable = it_cluster_freqs;
      }
    }
  }
  
  // Convert to Rcpp::List
  Rcpp::List clusters(n_records);
  Rcpp::NumericVector probabilities(n_records);
  for (int r0 = 0; r0 < n_records; r0++)
  {
    Rcpp::IntegerVector cluster = Rcpp::wrap(most_probable[r0]->first);
    clusters[r0] = cluster;
    probabilities[r0] = 1.0 / n_samples * most_probable[r0]->second ;
  }
  
  return Rcpp::List::create(Rcpp::Named("clusters") = clusters,
                            Rcpp::Named("probabilities") = probabilities);
}

// [[Rcpp::export(.smp_clusters)]]
Rcpp::List smp_clusters(const Rcpp::List& mp_clusters) {
  int n_records = mp_clusters.size();
  Rcpp::List out;
  
  // rec_ids are integers {1, ..., n_records}
  
  // TODO: make key an IntegerVector to avoid copying below
  // Define map to store the clusters (keys) and the corresponding
  // set of record ids (values) for which the cluster is most probable.
  map<set<int>, set<int> > smp_clusters;
  
  // Fill smp_clusters map
  for (int r0  = 0; r0 < n_records; r0++)
  {
    const Rcpp::IntegerVector &cluster = mp_clusters[r0];
    
    set<int> key(cluster.cbegin(), cluster.cend());
    smp_clusters[key].insert(r0 + 1); // rec_ids start at 1
  }
  
  out = Rcpp::List(smp_clusters.size());
  size_t i = 0;
  for (auto const &smpc : smp_clusters) {
    Rcpp::IntegerVector cluster = Rcpp::wrap(smpc.second);
    out[i] = cluster;
    i++;
  }
  
  // Convert to List of IntegerVectors: one for each cluster
  return out;
}