// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <algorithm>
using namespace Rcpp;

// compute_delta2(): the pooled estimate of the within-subject intraclass
// correlation rho. The per-peptide score statistics (cpp_delta2_score) live
// with the permutation kernel in delta2_contrast.cpp, so both use one
// implementation.

#include "delta2_common.h"
using namespace delta2;

// Pooled one-way random-effects ANOVA estimate of the intraclass
// correlation of the sample bits within subjects, summed over peptides.
// Also returns each subject's total positives (breadth) over pep_cols and
// each peptide's total positives over the samples.
//
// [[Rcpp::export]]
Rcpp::List cpp_delta2_rho(const Rcpp::RawVector& bitset_raw,
                          int n_words,
                          const Rcpp::IntegerVector& pep_cols,
                          const Rcpp::IntegerVector& sample_rows,
                          const Rcpp::IntegerVector& sample_subject) {
  if (n_words <= 0) stop("n_words must be positive.");
  IntegerVector time_code(sample_rows.size(), 1);
  Contrast cx;
  cx.init(n_words, sample_rows, sample_subject, time_code, nullptr, 1);
  const int k = cx.n_sub;
  const double Ns = cx.n_samples;

  int max_sub = 0;
  for (int s : cx.sub_orig) max_sub = std::max(max_sub, s);
  NumericVector hits(max_sub, 0.0);

  double SSB = 0.0, SSW = 0.0;
  Vec St(1);
  NumericVector pep_tot(pep_cols.size());
  for (int j = 0; j < pep_cols.size(); ++j) {
    const uint8_t* col = column_ptr(bitset_raw, n_words, pep_cols[j]);
    St[0] = 0.0;
    const double C = cx.count(col, n_words, St);
    pep_tot[j] = C;
    double ssb = 0.0;
    for (int i : cx.touched) {
      const double c = cx.cnt[i];
      const double ni = cx.n_i[i];
      ssb += c * c / ni;
      SSW += c * (ni - c) / ni;
      hits[cx.sub_orig[i] - 1] += c;
    }
    SSB += ssb - C * C / Ns;
    cx.reset();
  }

  double sum_n2 = 0.0;
  for (int i = 0; i < k; ++i) sum_n2 += cx.n_i[i] * cx.n_i[i];

  double rho_hat = NA_REAL, MSB = NA_REAL, MSW = NA_REAL, n0 = NA_REAL;
  if (k >= 2 && Ns > k) {
    MSB = SSB / (k - 1.0);
    MSW = SSW / (Ns - k);
    n0 = (Ns - sum_n2 / Ns) / (k - 1.0);
    const double den = MSB + (n0 - 1.0) * MSW;
    rho_hat = (den > 0.0) ? (MSB - MSW) / den : 0.0;
    rho_hat = std::min(1.0, std::max(0.0, rho_hat));
  }
  return List::create(
    _["rho_hat"] = rho_hat,
    _["MSB"] = MSB,
    _["MSW"] = MSW,
    _["n0"] = n0,
    _["subject_hits"] = hits,
    _["pep_tot"] = pep_tot,
    _["n_samples"] = Ns
  );
}
