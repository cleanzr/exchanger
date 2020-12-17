#include "alias_sampler.h"

void AliasSampler::build_tables() 
{  
  // At this point, probs_ vector contains valid weights, but is not normalized.
  
  size_t ii, jj, kk; // indices for loops
  
  // index tables, fill with zeros
  std::vector<size_t> HL_dat(n_items_, 0.0);
  std::vector<size_t>::iterator H, L, H0, L0;
  
  //HL0 = HL_dat.begin();
  H0 = H = HL_dat.begin();
  L0 = L = HL_dat.end();
  // fill HL_dat from beginning (small prob) and end (large prob) with indices
  for ( ii = 0; ii < n_items_; ii++ ) 
  {
    probs_[ii] *= n_items_ / total_weight_;
    if( probs_[ii] < 1.) {
      *(H++) = ii;
    } else {
      *(--L) = ii;
    }
  }
  
  // some of both large and small
  if ( (H > H0) && (L < L0) ) 
  {
    for ( kk = 0; kk < n_items_; kk++ ) 
    {
      ii = HL_dat[kk];
      jj = *L;
      alias_table_[ii] = jj;
      probs_[jj] += (probs_[ii] - 1.);
      if (probs_[jj] < 1.) L++;
      if(L == L0) break; // now all probs_ >= 1
    }
  }
  
  for ( ii = 0; ii < n_items_; ii++ )  probs_[ii] += ii;
}

int AliasSampler::draw() const 
{
  double rU = R::unif_rand() * n_items_;
  int kk = (int) rU;
  return ( rU < probs_[kk] ) ? kk : alias_table_[kk];
}

double AliasSampler::total_weight() const 
{
  return total_weight_;
}