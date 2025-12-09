// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <cstring>

using namespace Rcpp;

// ---- utilities ---------------------------------------------------------------

static inline uint32_t pc64(uint64_t x) {
#if defined(__GNUC__) || defined(__clang__)
  return (uint32_t)__builtin_popcountll(x);
#else
  // Fallback popcount
  x = x - ((x >> 1) & 0x5555555555555555ULL);
  x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
  return (uint32_t)((((x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL) * 0x0101010101010101ULL) >> 56);
#endif
}

// Build a dense mask (n_words uint64) with bits set for the given 1-based rows
static inline std::vector<uint64_t> build_mask(const IntegerVector& rows, int n_words) {
  std::vector<uint64_t> mask(n_words, 0ULL);
  for (int i = 0; i < rows.size(); ++i) {
    const int r = rows[i];
    if (r <= 0) continue;
    const int z = r - 1;
    mask[(z >> 6)] |= (1ULL << (z & 63));
  }
  return mask;
}

// Count intersections for one peptide column with a given mask (unpaired path)
static inline int count_hits_col(const uint8_t* bitset_raw, int n_words, int col0, const std::vector<uint64_t>& mask) {
  const uint64_t* colp = reinterpret_cast<const uint64_t*>(bitset_raw) + (size_t)col0 * (size_t)n_words;
  uint32_t s = 0;
  for (int w = 0; w < n_words; ++w) s += pc64(colp[w] & mask[w]);
  return (int)s;
}

// Build per-peptide bit masks from 1-based index lists (paired path)
static inline std::vector< std::vector<uint64_t> >
  lists_to_masks(const List& hits_list, const int n_words_p) {
    const int m = hits_list.size();
    std::vector< std::vector<uint64_t> > out; out.reserve(m);
    for (int i = 0; i < m; ++i) {
      IntegerVector idx = hits_list[i];
      std::vector<uint64_t> v(n_words_p, 0ULL);
      for (int k = 0; k < idx.size(); ++k) {
        const int r = idx[k];
        if (r <= 0) continue;
        const int z = r - 1;
        v[(z >> 6)] |= (1ULL << (z & 63));
      }
      out.push_back(std::move(v));
    }
    return out;
  }

// Compute z, weights, T_obs exactly as in R helper
struct CombineOut {
  double T_obs;
  std::vector<double> w_norm;
  std::vector<double> delta;
};

static CombineOut combine_T_internal(const std::vector<double>& p1,
                                     const std::vector<double>& p2,
                                     const std::vector<double>& n1,
                                     const std::vector<double>& n2,
                                     double winsor_z,
                                     const std::string& weight_mode,
                                     const std::string& stat_mode,
                                     const std::string& prev_strat) {
  const int m = (int)p1.size();
  std::vector<double> z(m), w(m), delta(m);

  // z & delta
  for (int i = 0; i < m; ++i) {
    const double d = p2[i] - p1[i];
    delta[i] = d;
    if (stat_mode == "asin") {
      const double den = std::sqrt( 1.0/(4.0*std::max(n1[i],1.0)) + 1.0/(4.0*std::max(n2[i],1.0)) );
      z[i] = (std::asin(std::sqrt(p2[i])) - std::asin(std::sqrt(p1[i]))) / den;
    } else { // diff
      const double se = std::sqrt(p1[i]*(1.0-p1[i])/std::max(n1[i],1.0) + p2[i]*(1.0-p2[i])/std::max(n2[i],1.0));
      z[i] = d / std::max(se, 1e-12);
    }
    if (z[i] >  winsor_z) z[i] =  winsor_z;
    if (z[i] < -winsor_z) z[i] = -winsor_z;
  }

  // weights
  if (weight_mode == "se_invvar" && stat_mode == "asin") {
    for (int i = 0; i < m; ++i) {
      w[i] = 1.0 / std::sqrt( 1.0/(4.0*std::max(n1[i],1.0)) + 1.0/(4.0*std::max(n2[i],1.0)) );
    }
  } else if (weight_mode == "se_invvar") {
    for (int i = 0; i < m; ++i) {
      const double se = std::sqrt(p1[i]*(1.0-p1[i])/std::max(n1[i],1.0) + p2[i]*(1.0-p2[i])/std::max(n2[i],1.0));
      w[i] = 1.0 / std::max(se, 1e-6);
    }
  } else if (weight_mode == "n_eff_sqrt") {
    for (int i = 0; i < m; ++i) {
      w[i] = std::sqrt(std::max(n1[i]*p1[i] + n2[i]*p2[i], 1e-12));
    }
  } else {
    std::fill(w.begin(), w.end(), 1.0);
  }

  // Stouffer combine
  double T_obs = 0.0;
  if (prev_strat == "decile") {
    // pooled prevalence
    std::vector<double> pp(m);
    for (int i = 0; i < m; ++i) pp[i] = (n1[i]*p1[i] + n2[i]*p2[i]) / std::max(n1[i] + n2[i], 1.0);

    // quantile breaks (approximate type 7)
    std::vector<double> sorted = pp;
    std::sort(sorted.begin(), sorted.end());
    std::vector<double> brks; brks.reserve(11);
    for (int k = 0; k <= 10; ++k) {
      double q = k/10.0;
      double pos = q * (m - 1);
      int lo = (int)std::floor(pos);
      int hi = (int)std::ceil(pos);
      double v = (lo == hi) ? sorted[lo] : (sorted[lo] + (pos - lo) * (sorted[hi] - sorted[lo]));
      if (brks.empty() || v > brks.back()) brks.push_back(v); // unique
    }
    if (brks.size() < 2) brks = {sorted.front(), sorted.back()};

    // bin & per-bin Stouffer
    int nb = (int)brks.size() - 1;
    std::vector<double> bin_z;
    bin_z.reserve(nb);
    for (int b = 0; b < nb; ++b) {
      const double L = brks[b], R = brks[b+1];
      double num = 0.0, den2 = 0.0; bool any=false;
      for (int i = 0; i < m; ++i) {
        const double v = pp[i];
        const bool in = ( (b==0 ? (v >= L) : (v > L)) && ( (b==nb-1 ? (v <= R) : (v < R)) ) );
        if (in) { num += w[i]*z[i]; den2 += w[i]*w[i]; any=true; }
      }
      bin_z.push_back(any ? (num / std::sqrt(std::max(den2, 1e-300))) : 0.0);
    }
    // mean of bin-level z's
    double s=0.0; for (double v: bin_z) s += v;
    T_obs = (nb>0) ? s / nb : 0.0;
  } else {
    double num = 0.0, den2 = 0.0;
    for (int i = 0; i < m; ++i) { num += w[i]*z[i]; den2 += w[i]*w[i]; }
    T_obs = num / std::sqrt(std::max(den2, 1e-300));
  }

  // normalized weights
  double sw = 0.0; for (double wi: w) sw += wi;
  std::vector<double> w_norm(m);
  const double denom = std::max(sw, 1e-12);
  for (int i = 0; i < m; ++i) w_norm[i] = w[i] / denom;

  return {T_obs, std::move(w_norm), std::move(delta)};
}

// ---- main exported helper ----------------------------------------------------

/**
 * Compute contrast using:
 *  - UNPAIRED: subject×peptide bitset + group row sets (preserve group sizes).
 *  - PAIRED: per-group per-peptide hit lists aligned to subjects 1..P; flip labels 50/50.
 *
 * @param bitset_raw RawVector: column-major, each column has n_words uint64 (UNPAIRED).
 * @param n_words int: number of 64-bit words per column (UNPAIRED).
 * @param pep_cols IntegerVector (1-based): peptide columns to use (UNPAIRED).
 * @param g1_rows, g2_rows IntegerVector (1-based): subject row indices per group (UNPAIRED).
 * @param hits_g1_paired, hits_g2_paired List: per-peptide IntegerVector of 1..P subject indices (PAIRED).
 * @param P int: number of paired subjects (PAIRED).
 * @param B, seed, smooth_eps_num, smooth_eps_den, min_max_prev, winsor_z numeric.
 * @param weight_mode, stat_mode, prev_strat, design strings.
 *
 * @return List with fields: n_peptides_used, m_eff, T_obs, T_null_sd, b, p_perm,
 *         mean_delta, frac_delta_pos, mean_delta_w, frac_delta_pos_w.
 */
// [[Rcpp::export]]
Rcpp::List cpp_shift_contrast(const Rcpp::RawVector& bitset_raw,
                              const int n_words,
                              const Rcpp::IntegerVector& pep_cols,
                              const Rcpp::IntegerVector& g1_rows,
                              const Rcpp::IntegerVector& g2_rows,
                              const Rcpp::List& hits_g1_paired,
                              const Rcpp::List& hits_g2_paired,
                              const int P,
                              const int B,
                              const int seed,
                              const double smooth_eps_num,
                              const double smooth_eps_den,
                              const double min_max_prev,
                              const std::string& weight_mode,
                              const std::string& stat_mode,
                              const std::string& prev_strat,
                              const double winsor_z,
                              const std::string& design) {

  // --------------------------- PAIRED PATH -----------------------------------
  if (design == "paired") {
    if (P <= 0) stop("cpp_shift_contrast(paired): P must be > 0.");
    if (hits_g1_paired.size() != hits_g2_paired.size())
      stop("cpp_shift_contrast(paired): hits_g1_paired and hits_g2_paired must have the same length.");

    const int m = hits_g1_paired.size();
    const int n_words_p = (P + 63) / 64;

    // Build per-peptide masks for g1 and g2 (aligned to subject order 1..P)
    std::vector< std::vector<uint64_t> > M1 = lists_to_masks(hits_g1_paired, n_words_p);
    std::vector< std::vector<uint64_t> > M2 = lists_to_masks(hits_g2_paired, n_words_p);

    // Observed counts x1=|S1|, x2=|S2|; n1=n2=P (constant)
    std::vector<double> p1; p1.reserve(m);
    std::vector<double> p2; p2.reserve(m);
    std::vector<double> n1(m, (double)P), n2(m, (double)P);
    std::vector<int> keep_idx; keep_idx.reserve(m);

    for (int i = 0; i < m; ++i) {
      // popcount masks
      uint32_t s1 = 0, s2 = 0;
      for (int w = 0; w < n_words_p; ++w) { s1 += pc64(M1[i][w]); s2 += pc64(M2[i][w]); }
      const double p1i = ( (double)s1 + smooth_eps_num ) / ( (double)P + smooth_eps_den * smooth_eps_num );
      const double p2i = ( (double)s2 + smooth_eps_num ) / ( (double)P + smooth_eps_den * smooth_eps_num );
      if (std::max(p1i, p2i) >= min_max_prev) {
        keep_idx.push_back(i);
        p1.push_back(p1i); p2.push_back(p2i);
      }
    }

    const int mu = (int)keep_idx.size();
    if (mu == 0) {
      return Rcpp::List::create(
        _["n_peptides_used"]  = 0L,
        _["m_eff"]            = NA_REAL,
        _["T_obs"]            = NA_REAL,
        _["T_null_sd"]        = NA_REAL,
        _["b"]                = 0L,
        _["p_perm"]           = NA_REAL,
        _["mean_delta"]       = NA_REAL,
        _["frac_delta_pos"]   = NA_REAL,
        _["mean_delta_w"]     = NA_REAL,
        _["frac_delta_pos_w"] = NA_REAL
      );
    }

    // align n1/n2 to kept length
    n1.assign(mu, (double)P);
    n2.assign(mu, (double)P);

    // Observed T
    CombineOut obs = combine_T_internal(p1, p2, n1, n2, winsor_z, weight_mode, stat_mode, prev_strat);

    // Moments
    double mean_delta = 0.0, frac_pos = 0.0, mean_delta_w = 0.0, frac_pos_w = 0.0;
    for (int i = 0; i < mu; ++i) {
      mean_delta += obs.delta[i];
      if (obs.delta[i] > 0) frac_pos += 1.0;
      mean_delta_w += obs.w_norm[i] * obs.delta[i];
      frac_pos_w   += obs.w_norm[i] * (obs.delta[i] > 0 ? 1.0 : 0.0);
    }
    mean_delta   /= (double)mu;
    frac_pos     /= (double)mu;

    // m_eff = 1 / sum(w_norm^2)
    double sum_w2 = 0.0;
    for (double wn : obs.w_norm) sum_w2 += wn * wn;
    const double m_eff = 1.0 / std::max(sum_w2, 1e-12);

    // Permutations: flip labels per subject with prob 0.5
    std::mt19937_64 rng((uint64_t)seed);
    std::bernoulli_distribution coin(0.5);

    int b_hits = 0;
    std::vector<uint64_t> Fmask(n_words_p, 0ULL); // flip mask
    double sum_T = 0.0, sum_T2 = 0.0;

    for (int b = 0; b < B; ++b) {
      // Build flip mask F over P subjects
      std::fill(Fmask.begin(), Fmask.end(), 0ULL);
      for (int j = 0; j < P; ++j) {
        if (coin(rng)) {
          const int w = (j >> 6), r = (j & 63);
          Fmask[w] |= (1ULL << r);
        }
      }

      // For each peptide, counts after flip:
      // x1' = pop(S1 & ~F) + pop(S2 & F)
      // x2' = pop(S2 & ~F) + pop(S1 & F)
      std::vector<double> p1b; p1b.reserve(mu);
      std::vector<double> p2b; p2b.reserve(mu);

      for (int k = 0; k < mu; ++k) {
        const int i = keep_idx[k];
        uint32_t a = 0, b2 = 0, c = 0, d = 0; // temp counts
        for (int w = 0; w < n_words_p; ++w) {
          const uint64_t S1 = M1[i][w], S2 = M2[i][w], F = Fmask[w];
          a += pc64(S1 & ~F);
          b2 += pc64(S2 &  F);
          c += pc64(S2 & ~F);
          d += pc64(S1 &  F);
        }
        const double x1p = (double)a + (double)b2;
        const double x2p = (double)c + (double)d;
        const double pA  = (x1p + smooth_eps_num) / ((double)P + smooth_eps_den * smooth_eps_num);
        const double pB  = (x2p + smooth_eps_num) / ((double)P + smooth_eps_den * smooth_eps_num);
        if (std::max(pA, pB) >= min_max_prev) {
          p1b.push_back(pA); p2b.push_back(pB);
        }
      }

      double Tb = 0.0;
      if (!p1b.empty()) {
        std::vector<double> n1b(p1b.size(), (double)P), n2b(p2b.size(), (double)P);
        Tb = combine_T_internal(p1b, p2b, n1b, n2b, winsor_z, weight_mode, stat_mode, prev_strat).T_obs;
      }

      sum_T  += Tb;
      sum_T2 += Tb * Tb;

      if (std::fabs(Tb) >= std::fabs(obs.T_obs)) ++b_hits;
    }

    const double p_perm = (1.0 + (double)b_hits) / (1.0 + (double)B);

    double T_null_sd = NA_REAL;
    if (B > 0) {
      const double Bd = (double)B;
      const double mean_T = sum_T / Bd;
      double var_T = (sum_T2 / Bd) - mean_T * mean_T;
      if (var_T > 0.0) T_null_sd = std::sqrt(var_T);
    }

    return Rcpp::List::create(
      _["n_peptides_used"]  = mu,
      _["m_eff"]            = m_eff,
      _["T_obs"]            = obs.T_obs,
      _["T_null_sd"]        = T_null_sd,
      _["b"]                = b_hits,
      _["p_perm"]           = p_perm,
      _["mean_delta"]       = mean_delta,
      _["frac_delta_pos"]   = frac_pos,
      _["mean_delta_w"]     = mean_delta_w,
      _["frac_delta_pos_w"] = frac_pos_w
    );
  }

  // --------------------------- UNPAIRED PATH ----------------------------------
  {
    // Build group masks
    std::vector<uint64_t> mask_g1 = build_mask(g1_rows, n_words);
    std::vector<uint64_t> mask_g2 = build_mask(g2_rows, n_words);

    // Observed counts per peptide
    const int m = pep_cols.size();
    std::vector<double> x1(m), x2(m);
    const uint8_t* raw = reinterpret_cast<const uint8_t*>(RAW(bitset_raw));

    for (int i = 0; i < m; ++i) {
      const int col0 = pep_cols[i] - 1;
      x1[i] = (double)count_hits_col(raw, n_words, col0, mask_g1);
      x2[i] = (double)count_hits_col(raw, n_words, col0, mask_g2);
    }

    const double n1_const = (double)g1_rows.size();
    const double n2_const = (double)g2_rows.size();
    if (n1_const <= 0 || n2_const <= 0) {
      return Rcpp::List::create(
        _["n_peptides_used"] = 0L,
        _["m_eff"] = NA_REAL,
        _["T_obs"] = NA_REAL,
        _["T_null_sd"] = NA_REAL,
        _["b"] = NA_INTEGER,
        _["p_perm"] = NA_REAL,
        _["mean_delta"] = NA_REAL,
        _["frac_delta_pos"] = NA_REAL,
        _["mean_delta_w"] = NA_REAL,
        _["frac_delta_pos_w"] = NA_REAL
      );
    }

    // Smoothed prevalences and filter by min_max_prev
    std::vector<double> p1, p2, n1, n2;
    p1.reserve(m); p2.reserve(m); n1.reserve(m); n2.reserve(m);
    std::vector<int> keep_idx; keep_idx.reserve(m);

    for (int i = 0; i < m; ++i) {
      const double p1i = (x1[i] + smooth_eps_num) / (n1_const + smooth_eps_den * smooth_eps_num);
      const double p2i = (x2[i] + smooth_eps_num) / (n2_const + smooth_eps_den * smooth_eps_num);
      if (std::max(p1i, p2i) >= min_max_prev) {
        keep_idx.push_back(i);
        p1.push_back(p1i); p2.push_back(p2i);
        n1.push_back(n1_const); n2.push_back(n2_const);
      }
    }

    const int mu = (int)keep_idx.size();
    if (mu == 0) {
      return Rcpp::List::create(
        _["n_peptides_used"]  = 0L,
        _["m_eff"]            = NA_REAL,
        _["T_obs"]            = NA_REAL,
        _["T_null_sd"]        = NA_REAL,
        _["b"]                = 0L,
        _["p_perm"]           = NA_REAL,
        _["mean_delta"]       = NA_REAL,
        _["frac_delta_pos"]   = NA_REAL,
        _["mean_delta_w"]     = NA_REAL,
        _["frac_delta_pos_w"] = NA_REAL
      );
    }

    // Observed T
    CombineOut obs = combine_T_internal(p1, p2, n1, n2, winsor_z, weight_mode, stat_mode, prev_strat);

    // Moments
    double mean_delta = 0.0, frac_pos = 0.0, mean_delta_w = 0.0, frac_pos_w = 0.0;
    for (int i = 0; i < mu; ++i) {
      mean_delta += obs.delta[i];
      if (obs.delta[i] > 0) frac_pos += 1.0;
      mean_delta_w += obs.w_norm[i] * obs.delta[i];
      frac_pos_w   += obs.w_norm[i] * (obs.delta[i] > 0 ? 1.0 : 0.0);
    }
    mean_delta   /= (double)mu;
    frac_pos     /= (double)mu;

    // m_eff = 1 / sum(w_norm^2)
    double sum_w2 = 0.0;
    for (double wn : obs.w_norm) sum_w2 += wn * wn;
    const double m_eff = 1.0 / std::max(sum_w2, 1e-12);

    // Permutations (two-sided add-one)
    std::mt19937_64 rng((uint64_t)seed);
    std::vector<int> idx_all(g1_rows.size() + g2_rows.size());
    for (size_t i = 0; i < idx_all.size(); ++i) idx_all[i] = (int)i;
    const int nA = (int)g1_rows.size();
    const int nB = (int)g2_rows.size();
    const int N  = nA + nB;

    int b_hits = 0;
    std::vector<int> choose_idx(nA);
    std::vector<uint64_t> mask_A(n_words), mask_B(n_words);
    double sum_T = 0.0, sum_T2 = 0.0;

    for (int b = 0; b < B; ++b) {
      // random split: sample nA indices for group A
      for (int i = 0; i < N; ++i) idx_all[i] = i;
      for (int i = 0; i < nA; ++i) {
        std::uniform_int_distribution<int> U(i, N-1);
        int j = U(rng);
        std::swap(idx_all[i], idx_all[j]);
        choose_idx[i] = idx_all[i];
      }
      // build masks from the chosen subject IDs
      std::fill(mask_A.begin(), mask_A.end(), 0ULL);
      std::fill(mask_B.begin(), mask_B.end(), 0ULL);
      std::vector<char> isA(N, 0);
      for (int i = 0; i < nA; ++i) {
        isA[ choose_idx[i] ] = 1;
        int subj_row = (choose_idx[i] < nA) ? g1_rows[ choose_idx[i] ] : g2_rows[ choose_idx[i]-nA ];
        if (subj_row > 0) {
          int z = subj_row - 1; mask_A[(z >> 6)] |= (1ULL << (z & 63));
        }
      }
      for (int t = 0; t < N; ++t) if (!isA[t]) {
        int subj_row = (t < nA) ? g1_rows[t] : g2_rows[t - nA];
        if (subj_row > 0) {
          int z = subj_row - 1; mask_B[(z >> 6)] |= (1ULL << (z & 63));
        }
      }

      // counts → smoothed p → filter → T_b
      std::vector<double> p1b; p1b.reserve(mu);
      std::vector<double> p2b; p2b.reserve(mu);
      std::vector<double> n1b(mu, (double)nA), n2b(mu, (double)nB);

      for (int k = 0; k < mu; ++k) {
        const int i_keep = keep_idx[k]; // index in original pep list
        const int col0 = pep_cols[i_keep] - 1;
        const double xA = (double)count_hits_col(reinterpret_cast<const uint8_t*>(RAW(bitset_raw)), n_words, col0, mask_A);
        const double xB = (double)count_hits_col(reinterpret_cast<const uint8_t*>(RAW(bitset_raw)), n_words, col0, mask_B);
        const double pA = (xA + smooth_eps_num) / (n1b[0] + smooth_eps_den * smooth_eps_num);
        const double pB = (xB + smooth_eps_num) / (n2b[0] + smooth_eps_den * smooth_eps_num);
        if (std::max(pA, pB) >= min_max_prev) {
          p1b.push_back(pA); p2b.push_back(pB);
        }
      }

      double Tb = 0.0;
      if (!p1b.empty()) {
        n1b.assign(p1b.size(), (double)nA);
        n2b.assign(p2b.size(), (double)nB);
        Tb = combine_T_internal(p1b, p2b, n1b, n2b, winsor_z, weight_mode, stat_mode, prev_strat).T_obs;
      }

      sum_T  += Tb;
      sum_T2 += Tb * Tb;

      if (std::fabs(Tb) >= std::fabs(obs.T_obs)) ++b_hits;
    }

    const double p_perm = (1.0 + (double)b_hits) / (1.0 + (double)B);

    double T_null_sd = NA_REAL;
    if (B > 0) {
      const double Bd = (double)B;
      const double mean_T = sum_T / Bd;
      double var_T = (sum_T2 / Bd) - mean_T * mean_T;
      if (var_T > 0.0) T_null_sd = std::sqrt(var_T);
    }

    return Rcpp::List::create(
      _["n_peptides_used"]  = mu,
      _["m_eff"]            = m_eff,
      _["T_obs"]            = obs.T_obs,
      _["T_null_sd"]        = T_null_sd,
      _["b"]                = b_hits,
      _["p_perm"]           = p_perm,
      _["mean_delta"]       = mean_delta,
      _["frac_delta_pos"]   = frac_pos,
      _["mean_delta_w"]     = mean_delta_w,
      _["frac_delta_pos_w"] = frac_pos_w
    );
  }
}
