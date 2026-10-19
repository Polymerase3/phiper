// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <vector>
#include <random>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <algorithm>
using namespace Rcpp;

// compute_delta2() resampling: one shared object per resample (a subject
// labelling, or the F[t][u] time masks), never anything per peptide.
//
//   group   between subjects; serves the group term and the interaction.
//           K = 2 with propensities: rejection sampler on independent
//           Bernoulli(p_i) draws, conditioned on the observed group size.
//           Otherwise: Metropolis-Hastings pairwise swaps (Berrett et al.
//           2020) in the Besag-Clifford parallel form, so the observed
//           labelling is exchangeable with the draws however short the chain.
//           Constant propensities: a uniform random split.
//   time    within subjects; a coin per two-sample subject (biased by
//           covariates when given), otherwise a uniform permutation of the
//           subject's time labels.

namespace {

typedef std::mt19937_64 Rng;

inline double u01(Rng& rng) {
  return std::uniform_real_distribution<double>(0.0, 1.0)(rng);
}

inline int uint_below(Rng& rng, int n) {
  return std::uniform_int_distribution<int>(0, n - 1)(rng);
}

// One Metropolis-Hastings swap step on labels (0-based) with propensity
// matrix q (N x K, column-major as from R). Returns 1 if a swap was made.
inline int swap_step(std::vector<int>& lab, const NumericMatrix& q, Rng& rng) {
  const int N = lab.size();
  const int i = uint_below(rng, N);
  int j = uint_below(rng, N - 1);
  if (j >= i) ++j;
  const int xi = lab[i], xj = lab[j];
  if (xi == xj) return 0;
  // the uniform pair proposal is symmetric, so this is the full MH ratio
  const double num = q(i, xj) * q(j, xi);
  const double den = q(i, xi) * q(j, xj);
  const bool accept = (den > 0.0) ? (num >= den || u01(rng) * den < num)
                                  : (num > 0.0);
  if (!accept) return 0;
  lab[i] = xj;
  lab[j] = xi;
  return 1;
}

} // namespace

// Resample subject labels. labels are 1..K (observed); prob is N x K
// (row i = subject i's propensity over levels), or 0 x 0 for constant
// propensities. Returns list(labels = N x B integer matrix, method,
// n_swaps_accepted).
//
// [[Rcpp::export]]
Rcpp::List cpp_delta2_resample_group(const Rcpp::IntegerVector& labels,
                                     const Rcpp::NumericMatrix& prob,
                                     int K, int B, int seed,
                                     int n_swap_steps, int max_tries) {
  const int N = labels.size();
  if (N < 2) stop("Need at least two subjects.");
  if (B < 1) stop("B must be positive.");
  std::vector<int> obs(N), sizes(K, 0);
  for (int i = 0; i < N; ++i) {
    const int l = labels[i];
    if (l == NA_INTEGER || l < 1 || l > K) stop("labels must be in 1..K.");
    obs[i] = l - 1;
    ++sizes[l - 1];
  }
  const bool constant = (prob.nrow() == 0);
  if (!constant && (prob.nrow() != N || prob.ncol() != K)) {
    stop("prob must be N x K.");
  }

  Rng rng((uint64_t)seed);
  IntegerMatrix out(N, B);
  std::string method;
  double n_acc = 0.0;

  if (constant) {
    // exact uniform split with the observed sizes: shuffle the labels
    method = "uniform";
    std::vector<int> lab = obs;
    for (int b = 0; b < B; ++b) {
      for (int i = N - 1; i > 0; --i) std::swap(lab[i], lab[uint_below(rng, i + 1)]);
      for (int i = 0; i < N; ++i) out(i, b) = lab[i] + 1;
    }
  } else {
    bool use_chain = (K != 2);
    if (!use_chain) {
      // rejection: independent Bernoulli(p_i), keep draws with |B| = n_B
      method = "rejection";
      std::vector<int> lab(N);
      for (int b = 0; b < B && !use_chain; ++b) {
        int tries = 0;
        int cnt;
        do {
          if (++tries > max_tries) {
            use_chain = true;
            break;
          }
          cnt = 0;
          for (int i = 0; i < N; ++i) {
            lab[i] = (u01(rng) < prob(i, 1)) ? 1 : 0;
            cnt += lab[i];
          }
        } while (cnt != sizes[1]);
        if (!use_chain) {
          for (int i = 0; i < N; ++i) out(i, b) = lab[i] + 1;
        }
      }
    }
    if (use_chain) {
      // Besag-Clifford parallel scheme: S steps from the observed labelling
      // to a hub (the chain is reversible, so forward = backward), then S
      // independent steps from the hub for each draw
      method = "swap";
      const int S = std::max(n_swap_steps, 1);
      std::vector<int> hub = obs;
      for (int s = 0; s < S; ++s) swap_step(hub, prob, rng);
      std::vector<int> lab(N);
      for (int b = 0; b < B; ++b) {
        lab = hub;
        for (int s = 0; s < S; ++s) n_acc += swap_step(lab, prob, rng);
        for (int i = 0; i < N; ++i) out(i, b) = lab[i] + 1;
      }
    }
  }

  return List::create(
    _["labels"] = out,
    _["method"] = method,
    _["n_swaps_accepted"] = n_acc
  );
}

// Resample time labels within subjects. Per sample: sample_rows (1-based
// bitset rows), sample_subject, time_code (1..T). keep_prob is per subject
// (indexed by sample_subject): for a subject with exactly two samples at two
// different timepoints, the probability of keeping the observed orientation
// (NA = fair coin); other subjects get a uniform permutation of their time
// labels. Returns list(time_code = n_samples x B, masks = raw of B blocks of
// T*T columns with n_words uint64 words, F[t][u] at column t*T + u: samples
// observed at t that the draw relabels u).
//
// [[Rcpp::export]]
Rcpp::List cpp_delta2_resample_time(const Rcpp::IntegerVector& sample_rows,
                                    const Rcpp::IntegerVector& sample_subject,
                                    const Rcpp::IntegerVector& time_code,
                                    const Rcpp::NumericVector& keep_prob,
                                    int n_time_levels, int n_words,
                                    int B, int seed, bool masks) {
  const int ns = sample_rows.size();
  const int T = n_time_levels;
  if (sample_subject.size() != ns || time_code.size() != ns) {
    stop("sample_rows, sample_subject and time_code must have equal length.");
  }
  if (T < 1 || n_words < 1 || B < 1) stop("Invalid dimensions.");

  int max_sub = 0;
  for (int s = 0; s < ns; ++s) {
    if (sample_subject[s] == NA_INTEGER || sample_subject[s] < 1) stop("Invalid sample_subject.");
    if (time_code[s] == NA_INTEGER || time_code[s] < 1 || time_code[s] > T) stop("Invalid time_code.");
    if (sample_rows[s] == NA_INTEGER || sample_rows[s] < 1 || sample_rows[s] > n_words * 64) {
      stop("sample_rows out of range.");
    }
    max_sub = std::max(max_sub, (int)sample_subject[s]);
  }
  if (keep_prob.size() > 0 && keep_prob.size() < max_sub) {
    stop("keep_prob must have one value per subject.");
  }

  // subject blocks: the sample positions of each subject
  std::vector<std::vector<int>> blocks(max_sub);
  for (int s = 0; s < ns; ++s) blocks[sample_subject[s] - 1].push_back(s);

  Rng rng((uint64_t)seed);
  IntegerMatrix out(ns, B);
  const size_t block_words = (size_t)T * T * n_words;
  RawVector mask_raw(masks ? (R_xlen_t)(block_words * B * 8) : 0);
  std::vector<uint64_t> F(masks ? block_words : 0);
  std::vector<int> lab;

  for (int b = 0; b < B; ++b) {
    for (int k = 0; k < max_sub; ++k) {
      const std::vector<int>& blk = blocks[k];
      const int n = blk.size();
      if (n == 0) continue;
      lab.resize(n);
      for (int a = 0; a < n; ++a) lab[a] = time_code[blk[a]];
      if (n == 2 && lab[0] != lab[1]) {
        double pk = 0.5;
        if (keep_prob.size() > 0 && !NumericVector::is_na(keep_prob[k])) pk = keep_prob[k];
        if (u01(rng) >= pk) std::swap(lab[0], lab[1]);
      } else {
        for (int a = n - 1; a > 0; --a) std::swap(lab[a], lab[uint_below(rng, a + 1)]);
      }
      for (int a = 0; a < n; ++a) out(blk[a], b) = lab[a];
    }
    if (masks) {
      std::fill(F.begin(), F.end(), 0ULL);
      for (int s = 0; s < ns; ++s) {
        const int t = time_code[s] - 1, u = out(s, b) - 1;
        const int z = sample_rows[s] - 1;
        F[(size_t)(t * T + u) * n_words + (z >> 6)] |= (1ULL << (z & 63));
      }
      std::memcpy(RAW(mask_raw) + block_words * 8 * b, F.data(), block_words * 8);
    }
  }

  return List::create(_["time_code"] = out, _["masks"] = mask_raw);
}

// Per-peptide pass over one resample's time masks: positives relabelled to
// each timepoint, count[j, u] = sum_t popcount(col_j & F[t][u]). The masks
// are shared by all peptides.
//
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_delta2_time_counts(const Rcpp::RawVector& bitset_raw,
                                           int n_words,
                                           const Rcpp::IntegerVector& pep_cols,
                                           const Rcpp::RawVector& masks,
                                           int n_time_levels, int b) {
  const int T = n_time_levels;
  const int m = pep_cols.size();
  const size_t col_bytes = (size_t)n_words * 8;
  const size_t block_bytes = (size_t)T * T * col_bytes;
  if (b < 1 || block_bytes * b > (size_t)masks.size()) stop("b out of range of masks.");
  const uint8_t* F = reinterpret_cast<const uint8_t*>(RAW(masks)) + block_bytes * (b - 1);
  const uint8_t* bits = reinterpret_cast<const uint8_t*>(RAW(bitset_raw));

  NumericMatrix out(m, T);
  for (int j = 0; j < m; ++j) {
    const size_t off = (size_t)(pep_cols[j] - 1) * col_bytes;
    if (pep_cols[j] < 1 || off + col_bytes > (size_t)bitset_raw.size()) {
      stop("pep_cols out of range of the bitset.");
    }
    for (int t = 0; t < T; ++t) {
      for (int u = 0; u < T; ++u) {
        const uint8_t* f = F + (size_t)(t * T + u) * col_bytes;
        double c = 0.0;
        for (int w = 0; w < n_words; ++w) {
          uint64_t x, y;
          std::memcpy(&x, bits + off + (size_t)w * 8, 8);
          std::memcpy(&y, f + (size_t)w * 8, 8);
          c += __builtin_popcountll(x & y);
        }
        out(j, u) += c;
      }
    }
  }
  return out;
}
