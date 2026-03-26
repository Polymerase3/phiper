// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>


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

// Welford's online algorithm for running mean and variance
struct WelfordAccum {
  int64_t n = 0;
  double mean = 0.0;
  double M2   = 0.0;
  void digest(double x) {
    ++n;
    const double d  = x - mean;
    mean += d / static_cast<double>(n);
    M2   += d * (x - mean);
  }
  double variance() const {
    return (n > 1) ? M2 / static_cast<double>(n - 1) : NA_REAL;
  }
};

// Compute z, weights, T_obs exactly as in R helper
struct CombineOut {
  double T_obs;
  std::vector<double> w_norm;
  std::vector<double> delta;
};

static inline Rcpp::List make_perm_result(
    int mu, double m_eff, double T_obs,
    const WelfordAccum& acc, int b_hits, int B,
    double max_delta, double frac_pos, double frac_pos_w)
{
  const double p_perm = (1.0 + (double)b_hits) / (1.0 + (double)B);
  double T_null_mean = NA_REAL;
  double T_null_sd   = NA_REAL;
  if (acc.n > 0) {
    T_null_mean = acc.mean;
    double var_T = acc.variance();
    if (var_T < 0.0 && var_T > -1e-15) var_T = 0.0;
    if (var_T >= 0.0) T_null_sd = std::sqrt(var_T);
  }
  return Rcpp::List::create(
    Rcpp::_["n_peptides_used"]  = mu,
    Rcpp::_["m_eff"]            = m_eff,
    Rcpp::_["T_obs"]            = T_obs,
    Rcpp::_["T_null_mean"]      = T_null_mean,
    Rcpp::_["T_null_sd"]        = T_null_sd,
    Rcpp::_["b"]                = b_hits,
    Rcpp::_["p_perm"]           = p_perm,
    Rcpp::_["max_delta"]        = max_delta,
    Rcpp::_["frac_delta_pos"]   = frac_pos,
    Rcpp::_["frac_delta_pos_w"] = frac_pos_w
  );
}

static inline double ll_binom(double x, double n, double p) {
  if (p <= 0.0) {
    return (x > 0.0) ? -std::numeric_limits<double>::infinity() : 0.0;
  }
  if (p >= 1.0) {
    return ((n - x) > 0.0) ? -std::numeric_limits<double>::infinity() : 0.0;
  }
  double out = 0.0;
  if (x > 0.0) out += x * std::log(p);
  const double nx = n - x;
  if (nx > 0.0) out += nx * std::log1p(-p);
  return out;
}

static inline bool score_boundary(double x1, double n1, double x2, double n2) {
  const double total = x1 + x2;
  const double N = n1 + n2;
  return (total <= 0.0) || (N > 0.0 && total >= N);
}

static inline double score_denom_from_counts(double x1, double n1,
                                             double x2, double n2) {
  const double N = n1 + n2;
  if (N <= 0.0 || n1 <= 0.0 || n2 <= 0.0) return 0.0;
  const double p0 = (x1 + x2) / N;
  const double v = p0 * (1.0 - p0) * (1.0 / n1 + 1.0 / n2);
  return (v > 0.0) ? std::sqrt(v) : 0.0;
}

static CombineOut combine_T_internal(const std::vector<double>& p1,
                                     const std::vector<double>& p2,
                                     const std::vector<double>& x1,
                                     const std::vector<double>& x2,
                                     const std::vector<double>& n1,
                                     const std::vector<double>& n2,
                                     const std::vector<double>& b,
                                     const std::vector<double>& c,
                                     double winsor_z,
                                     const std::string& weight_mode,
                                     const std::string& stat_mode,
                                     const std::string& aggregate_stat,
                                     const std::vector<double>& strat_bins) {
  const int m = (int)p1.size();
  std::vector<double> z(m), w(m), delta(m);

  // z & delta
  for (int i = 0; i < m; ++i) {
    const double d = p2[i] - p1[i];
    delta[i] = d;
    if (stat_mode == "mcnemar") {
      if ((int)b.size() != m || (int)c.size() != m) stop("mcnemar: b/c length mismatch.");
      const double m_dc = b[i] + c[i];
      z[i] = (m_dc > 0.0) ? ((b[i] - c[i]) / std::sqrt(m_dc)) : 0.0;
    } else if (stat_mode == "srlr_paired") {
      if ((int)b.size() != m || (int)c.size() != m) stop("srlr_paired: b/c length mismatch.");
      const double m_dc = b[i] + c[i];
      if (m_dc <= 0.0) {
        z[i] = 0.0;
      } else {
        const double t1 = (b[i] > 0.0) ? (b[i] * std::log((2.0 * b[i]) / m_dc)) : 0.0;
        const double t2 = (c[i] > 0.0) ? (c[i] * std::log((2.0 * c[i]) / m_dc)) : 0.0;
        const double LR = 2.0 * (t1 + t2);
        const double sign = (b[i] > c[i]) ? 1.0 : (b[i] < c[i]) ? -1.0 : 0.0;
        z[i] = (LR > 0.0) ? (sign * std::sqrt(LR)) : 0.0;
      }
    } else if (stat_mode == "asin") {
      const double den = std::sqrt( 1.0/(4.0*std::max(n1[i],1.0)) + 1.0/(4.0*std::max(n2[i],1.0)) );
      z[i] = (std::asin(std::sqrt(p2[i])) - std::asin(std::sqrt(p1[i]))) / den;
    } else if (stat_mode == "score" || stat_mode == "srlr") {
      const double n1i = n1[i];
      const double n2i = n2[i];
      const double x1i = x1[i];
      const double x2i = x2[i];

      if (stat_mode == "score" && score_boundary(x1i, n1i, x2i, n2i)) {
        z[i] = 0.0;
      } else {
        const double p1m = (n1i > 0.0) ? (x1i / n1i) : 0.0;
        const double p2m = (n2i > 0.0) ? (x2i / n2i) : 0.0;
        const long double lhs = (long double)x2i * (long double)n1i;
        const long double rhs = (long double)x1i * (long double)n2i;
        const double sign = (lhs > rhs) ? 1.0 : (lhs < rhs) ? -1.0 : 0.0;

        if (stat_mode == "score") {
          const double den = score_denom_from_counts(x1i, n1i, x2i, n2i);
          z[i] = (den > 0.0) ? ((p2m - p1m) / den) : 0.0;
        } else { // "srlr"
          const double p0 = (n1i + n2i > 0.0) ? ((x1i + x2i) / (n1i + n2i)) : 0.0;
          const double ll_alt  = ll_binom(x1i, n1i, p1m) + ll_binom(x2i, n2i, p2m);
          const double ll_null = ll_binom(x1i, n1i, p0)  + ll_binom(x2i, n2i, p0);
          double LR = 2.0 * (ll_alt - ll_null);
          if (LR < 0.0 && LR > -1e-12) LR = 0.0;
          z[i] = (LR > 0.0) ? (sign * std::sqrt(LR)) : 0.0;
        }
      }
    } else { // diff
      const double se = std::sqrt(p1[i]*(1.0-p1[i])/std::max(n1[i],1.0) + p2[i]*(1.0-p2[i])/std::max(n2[i],1.0));
      z[i] = d / std::max(se, 1e-12);
    }
    if (z[i] >  winsor_z) z[i] =  winsor_z;
    if (z[i] < -winsor_z) z[i] = -winsor_z;
  }

  // weights
  if (weight_mode == "se_invvar" && stat_mode == "srlr") {
    std::fill(w.begin(), w.end(), 1.0);
  } else if (weight_mode == "se_invvar" && stat_mode == "score") {
    for (int i = 0; i < m; ++i) {
      if (score_boundary(x1[i], n1[i], x2[i], n2[i])) {
        w[i] = 0.0;
      } else {
        const double den = score_denom_from_counts(x1[i], n1[i], x2[i], n2[i]);
        w[i] = 1.0 / std::max(den, 1e-6);
      }
    }
  } else if (weight_mode == "se_invvar" && stat_mode == "asin") {
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

  auto maxmean_stat = [&](const std::vector<int>& idx) {
    double sp = 0.0, sn = 0.0;
    int cp = 0, cn = 0;
    for (int i : idx) {
      if (z[i] > 0.0) { sp += z[i]; ++cp; }
      else if (z[i] < 0.0) { sn += -z[i]; ++cn; }
    }
    const double Splus = (cp > 0) ? (sp / (double)cp) : 0.0;
    const double Sminus = (cn > 0) ? (sn / (double)cn) : 0.0;
    return (Splus >= Sminus) ? Splus : -Sminus;
  };

  // Aggregate combine
  double T_obs = 0.0;
  bool use_strat = true;
  if (strat_bins.empty()) use_strat = false;
  if (strat_bins.size() == 1 && strat_bins[0] == 0.0) use_strat = false;

  if (use_strat) {
    // pooled prevalence
    std::vector<double> pp(m);
    for (int i = 0; i < m; ++i) pp[i] = (n1[i]*p1[i] + n2[i]*p2[i]) / std::max(n1[i] + n2[i], 1.0);

    // fixed cutpoints (0 < v < 1), with 0/1 boundaries added
    std::vector<double> cuts;
    cuts.reserve(strat_bins.size());
    for (double v : strat_bins) {
      if (std::isfinite(v) && v > 0.0 && v < 1.0) cuts.push_back(v);
    }
    if (cuts.empty()) {
      use_strat = false;
    } else {
      std::sort(cuts.begin(), cuts.end());
      cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
    }

    if (use_strat) {
      std::vector<double> brks;
      brks.reserve(cuts.size() + 2);
      brks.push_back(0.0);
      brks.insert(brks.end(), cuts.begin(), cuts.end());
      brks.push_back(1.0);

      // bin & per-bin Stouffer — only non-empty bins contribute to the mean
      int nb = (int)brks.size() - 1;
      double s = 0.0; int n_nonempty = 0;
      for (int b = 0; b < nb; ++b) {
        const double L = brks[b], R = brks[b+1];
        double num = 0.0, den2 = 0.0; bool any=false;
        std::vector<int> idx;
        for (int i = 0; i < m; ++i) {
          const double v = pp[i];
          const bool in = ((b == 0 ? (v >= L) : (v > L)) && (v <= R));
          if (in) {
            if (aggregate_stat == "maxmean") {
              idx.push_back(i);
            } else {
              num += w[i]*z[i];
              den2 += w[i]*w[i];
              any=true;
            }
          }
        }
        if (aggregate_stat == "maxmean") {
          if (!idx.empty()) { s += maxmean_stat(idx); ++n_nonempty; }
        } else {
          if (any) { s += num / std::sqrt(std::max(den2, 1e-300)); ++n_nonempty; }
        }
      }
      // mean over non-empty bins only
      T_obs = (n_nonempty > 0) ? s / n_nonempty : 0.0;
    }
  }
  if (!use_strat) {
    if (aggregate_stat == "maxmean") {
      std::vector<int> idx(m);
      for (int i = 0; i < m; ++i) idx[i] = i;
      T_obs = maxmean_stat(idx);
    } else {
      double num = 0.0, den2 = 0.0;
      for (int i = 0; i < m; ++i) { num += w[i]*z[i]; den2 += w[i]*w[i]; }
      T_obs = num / std::sqrt(std::max(den2, 1e-300));
    }
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
 * @param B, seed, winsor_z numeric.
 * @param weight_mode, stat_mode, aggregate_stat, design strings.
 * @param strat_bins Numeric vector of pooled prevalence cutpoints or 0 for no stratification.
 *
 * @return List with fields: n_peptides_used, m_eff, T_obs, T_null_mean, T_null_sd,
 *         b, p_perm, max_delta, frac_delta_pos, frac_delta_pos_w.
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
                              const std::string& weight_mode,
                              const std::string& stat_mode,
                              const std::string& aggregate_stat,
                              const Rcpp::NumericVector& strat_bins,
                              const double winsor_z,
                              const std::string& design) {
  std::vector<double> strat_bins_vec = Rcpp::as<std::vector<double> >(strat_bins);
  const bool paired_only_stat = (stat_mode == "mcnemar" || stat_mode == "srlr_paired");

  if (paired_only_stat && design != "paired") {
    stop("stat_mode requires paired design: %s", stat_mode);
  }

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
    std::vector<double> x1; x1.reserve(m);
    std::vector<double> x2; x2.reserve(m);
    std::vector<double> p1; p1.reserve(m);
    std::vector<double> p2; p2.reserve(m);
    std::vector<double> b_vec; if (paired_only_stat) b_vec.reserve(m);
    std::vector<double> c_vec; if (paired_only_stat) c_vec.reserve(m);
    std::vector<double> n1(m, (double)P), n2(m, (double)P);
    std::vector<int> keep_idx; keep_idx.reserve(m);

    for (int i = 0; i < m; ++i) {
      // popcount masks
      uint32_t s1 = 0, s2 = 0, b = 0, c = 0;
      for (int w = 0; w < n_words_p; ++w) {
        const uint64_t S1 = M1[i][w];
        const uint64_t S2 = M2[i][w];
        s1 += pc64(S1);
        s2 += pc64(S2);
        if (paired_only_stat) {
          b += pc64(S1 & ~S2);
          c += pc64(~S1 & S2);
        }
      }
      const double x1i = (double)s1;
      const double x2i = (double)s2;
      const double p1i = x1i / (double)P;
      const double p2i = x2i / (double)P;
      keep_idx.push_back(i);
      x1.push_back(x1i); x2.push_back(x2i);
      p1.push_back(p1i); p2.push_back(p2i);
      if (paired_only_stat) {
        b_vec.push_back((double)b);
        c_vec.push_back((double)c);
      }
    }

    const int mu = (int)keep_idx.size();
    if (mu == 0) {
      return Rcpp::List::create(
        _["n_peptides_used"]  = 0L,
        _["m_eff"]            = NA_REAL,
        _["T_obs"]            = NA_REAL,
        _["T_null_mean"]      = NA_REAL,
        _["T_null_sd"]        = NA_REAL,
        _["b"]                = 0L,
        _["p_perm"]           = NA_REAL,
        _["max_delta"]       = NA_REAL,
        _["frac_delta_pos"]   = NA_REAL,
        _["frac_delta_pos_w"] = NA_REAL
      );
    }

    // align n1/n2 to kept length
    n1.assign(mu, (double)P);
    n2.assign(mu, (double)P);

    // Observed T
    CombineOut obs = combine_T_internal(p1, p2, x1, x2, n1, n2, b_vec, c_vec, winsor_z, weight_mode, stat_mode, aggregate_stat, strat_bins_vec);

    // Moments
    double max_delta = 0.0;  // max absolute delta
    double frac_pos = 0.0, frac_pos_w = 0.0;
    for (int i = 0; i < mu; ++i) {
      const double d = obs.delta[i];
      const double ad = std::fabs(d);

      if (ad > max_delta) max_delta = ad;
      if (d > 0) frac_pos += 1.0;
      frac_pos_w += obs.w_norm[i] * (d > 0 ? 1.0 : 0.0);
    }
    frac_pos /= (double)mu;

    // m_eff = 1 / sum(w_norm^2)
    double sum_w2 = 0.0;
    for (double wn : obs.w_norm) sum_w2 += wn * wn;
    const double m_eff = 1.0 / std::max(sum_w2, 1e-12);

    // Permutations: flip labels per subject with prob 0.5
    std::mt19937_64 rng((uint64_t)seed);
    std::bernoulli_distribution coin(0.5);

    int b_hits = 0;
    std::vector<uint64_t> Fmask(n_words_p, 0ULL); // flip mask
    WelfordAccum acc;

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
      std::vector<double> x1b; x1b.reserve(mu);
      std::vector<double> x2b; x2b.reserve(mu);
      std::vector<double> p1b; p1b.reserve(mu);
      std::vector<double> p2b; p2b.reserve(mu);
      std::vector<double> bb; if (paired_only_stat) bb.reserve(mu);
      std::vector<double> cc; if (paired_only_stat) cc.reserve(mu);

      for (int k = 0; k < mu; ++k) {
        const int i = keep_idx[k];
        uint32_t a = 0, b2 = 0, c = 0, d = 0; // temp counts
        uint32_t b_disc = 0, c_disc = 0;
        for (int w = 0; w < n_words_p; ++w) {
          const uint64_t S1 = M1[i][w], S2 = M2[i][w], F = Fmask[w];
          a += pc64(S1 & ~F);
          b2 += pc64(S2 &  F);
          c += pc64(S2 & ~F);
          d += pc64(S1 &  F);
          if (paired_only_stat) {
            const uint64_t S1p = (S1 & ~F) | (S2 & F);
            const uint64_t S2p = (S2 & ~F) | (S1 & F);
            b_disc += pc64(S1p & ~S2p);
            c_disc += pc64(~S1p & S2p);
          }
        }
        const double x1p = (double)a + (double)b2;
        const double x2p = (double)c + (double)d;
        const double pA  = x1p / (double)P;
        const double pB  = x2p / (double)P;
        x1b.push_back(x1p); x2b.push_back(x2p);
        p1b.push_back(pA); p2b.push_back(pB);
        if (paired_only_stat) {
          bb.push_back((double)b_disc);
          cc.push_back((double)c_disc);
        }
      }

      double Tb = 0.0;
      if (!p1b.empty()) {
        std::vector<double> n1b(p1b.size(), (double)P), n2b(p2b.size(), (double)P);
        Tb = combine_T_internal(p1b, p2b, x1b, x2b, n1b, n2b, bb, cc, winsor_z, weight_mode, stat_mode, aggregate_stat, strat_bins_vec).T_obs;
      }

      acc.digest(Tb);

      if (std::fabs(Tb) >= std::fabs(obs.T_obs)) ++b_hits;
    }

    return make_perm_result(mu, m_eff, obs.T_obs, acc, b_hits, B,
                            max_delta, frac_pos, frac_pos_w);
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
        _["T_null_mean"] = NA_REAL,
        _["T_null_sd"] = NA_REAL,
        _["b"] = NA_INTEGER,
        _["p_perm"] = NA_REAL,
        _["max_delta"] = NA_REAL,
        _["frac_delta_pos"] = NA_REAL,
        _["frac_delta_pos_w"] = NA_REAL
      );
    }

    // Prevalences
    std::vector<double> p1, p2, n1, n2;
    p1.reserve(m); p2.reserve(m); n1.reserve(m); n2.reserve(m);
    std::vector<int> keep_idx; keep_idx.reserve(m);

    for (int i = 0; i < m; ++i) {
      const double p1i = x1[i] / n1_const;
      const double p2i = x2[i] / n2_const;
      keep_idx.push_back(i);
      p1.push_back(p1i); p2.push_back(p2i);
      n1.push_back(n1_const); n2.push_back(n2_const);
    }

    const int mu = (int)keep_idx.size();
    if (mu == 0) {
      return Rcpp::List::create(
        _["n_peptides_used"]  = 0L,
        _["m_eff"]            = NA_REAL,
        _["T_obs"]            = NA_REAL,
        _["T_null_mean"]      = NA_REAL,
        _["T_null_sd"]        = NA_REAL,
        _["b"]                = 0L,
        _["p_perm"]           = NA_REAL,
        _["max_delta"]        = NA_REAL,
        _["frac_delta_pos"]   = NA_REAL,
        _["frac_delta_pos_w"] = NA_REAL
      );
    }

    // Observed T
    CombineOut obs = combine_T_internal(p1, p2, x1, x2, n1, n2, std::vector<double>(), std::vector<double>(), winsor_z, weight_mode, stat_mode, aggregate_stat, strat_bins_vec);

    // Moments
    double max_delta = 0.0;  // max absolute delta
    double frac_pos = 0.0, frac_pos_w = 0.0;
    for (int i = 0; i < mu; ++i) {
      const double d = obs.delta[i];
      const double ad = std::fabs(d);

      if (ad > max_delta) max_delta = ad;
      if (d > 0) frac_pos += 1.0;
      frac_pos_w += obs.w_norm[i] * (d > 0 ? 1.0 : 0.0);
    }
    frac_pos /= (double)mu;

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
    WelfordAccum acc;

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

      // counts --> p --> T_b
      std::vector<double> x1b; x1b.reserve(mu);
      std::vector<double> x2b; x2b.reserve(mu);
      std::vector<double> p1b; p1b.reserve(mu);
      std::vector<double> p2b; p2b.reserve(mu);
      std::vector<double> n1b(mu, (double)nA), n2b(mu, (double)nB);

      for (int k = 0; k < mu; ++k) {
        const int i_keep = keep_idx[k]; // index in original pep list
        const int col0 = pep_cols[i_keep] - 1;
        const double xA = (double)count_hits_col(reinterpret_cast<const uint8_t*>(RAW(bitset_raw)), n_words, col0, mask_A);
        const double xB = (double)count_hits_col(reinterpret_cast<const uint8_t*>(RAW(bitset_raw)), n_words, col0, mask_B);
        const double pA = xA / n1b[0];
        const double pB = xB / n2b[0];
        x1b.push_back(xA); x2b.push_back(xB);
        p1b.push_back(pA); p2b.push_back(pB);
      }

      double Tb = 0.0;
      if (!p1b.empty()) {
        n1b.assign(p1b.size(), (double)nA);
        n2b.assign(p2b.size(), (double)nB);
        Tb = combine_T_internal(p1b, p2b, x1b, x2b, n1b, n2b, std::vector<double>(), std::vector<double>(), winsor_z, weight_mode, stat_mode, aggregate_stat, strat_bins_vec).T_obs;
      }

      acc.digest(Tb);

      if (std::fabs(Tb) >= std::fabs(obs.T_obs)) ++b_hits;
    }

    return make_perm_result(mu, m_eff, obs.T_obs, acc, b_hits, B,
                            max_delta, frac_pos, frac_pos_w);
  }
}
