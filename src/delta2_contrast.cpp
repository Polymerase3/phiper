// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <algorithm>
#include <string>
using namespace Rcpp;

#include "delta2_common.h"
using namespace delta2;

// compute_delta2() statistics of one (rank, feature, group_col) stratum.
//
// Stratum computes the per-peptide statistic of each term under a subject
// labelling (group term, interaction, pairwise tests) or a within-subject
// time relabelling (time term):
//   group        weighted least-squares score for group: subject values
//                u_i = w_i (ybar_i - Z_i alpha), alpha the w-weighted fit of
//                the subject means on the covariates Z (an intercept, plus
//                the propensity covariates when given); standardised by the
//                exact permutation variance of sum_i (g_i - gbar) u_i
//   time         within-subject deviations, exact within-subject
//                permutation variance
//   interaction  within-subject deviations residualised against the pooled
//                time main effect, compared across groups and standardised
//                by a Welch-type variance: one residual variance per group,
//                design-shaped, with the estimated time effect's pull on the
//                score, recomputed for every labelling
// A df = 1 term gives a signed z; an omnibus term the vector Cov^{+1/2} U.
//
// cpp_delta2_score() returns these per peptide for the observed labels.
// cpp_delta2_contrast() pools them over peptides (strata of pooled
// prevalence, Stouffer / maxmean / AF per component, T = sum of squared
// components for df > 1) and repeats that for every draw, draws outside and
// peptides inside: each draw builds its shared object once (a subject
// labelling with its per-level sums, or the per-timepoint sample masks) and
// nothing per peptide.

namespace {

enum class AggStat { stouffer, maxmean, af };

AggStat parse_agg(const std::string& s) {
  if (s == "maxmean") return AggStat::maxmean;
  if (s == "af") return AggStat::af;
  return AggStat::stouffer;
}

// strata of pooled prevalence: same breaks and bin rule as compute_delta()
std::vector<int> prevalence_bins(const NumericVector& strat_bins,
                                 const Vec& p0, int& n_bins) {
  Vec cuts;
  if (!(strat_bins.size() == 1 && strat_bins[0] == 0.0)) {
    for (double v : strat_bins) {
      if (std::isfinite(v) && v > 0.0 && v < 1.0) cuts.push_back(v);
    }
  }
  std::sort(cuts.begin(), cuts.end());
  cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
  const int m = p0.size();
  std::vector<int> bin(m, 0);
  if (cuts.empty()) {
    n_bins = 1;
    return bin;
  }
  Vec brks;
  brks.push_back(0.0);
  brks.insert(brks.end(), cuts.begin(), cuts.end());
  brks.push_back(1.0);
  n_bins = (int)brks.size() - 1;
  for (int j = 0; j < m; ++j) {
    bin[j] = -1;
    for (int b = 0; b < n_bins; ++b) {
      const double L = brks[b], R = brks[b + 1], v = p0[j];
      if ((b == 0 ? (v >= L) : (v > L)) && (v <= R)) {
        bin[j] = b;
        break;
      }
    }
  }
  return bin;
}

// Pools one component's peptide z, winsorised at `winsor`, into a scalar:
// the statistic per non-empty prevalence bin, averaged over those bins.
struct Pooler {
  AggStat stat;
  double winsor = 4.0;
  int n_bins = 1;
  std::vector<int> bin;
  Vec w;
  // scratch
  Vec num, den2, sp, sn, ppos, pneg;
  std::vector<int> cnt, cp, cn;
  std::vector<std::vector<double>> zs;

  static double af_one(Vec& pv) {
    std::sort(pv.begin(), pv.end());
    double s = 0.0, hk = 0.0, best = 0.0;
    for (int j = 0; j < (int)pv.size(); ++j) {
      s += -std::log(pv[j]);
      hk += 1.0 / (j + 1);
      best = std::max(best, s / hk);
    }
    return best;
  }

  double operator()(const double* z) {
    const int m = bin.size();
    num.assign(n_bins, 0.0);
    den2.assign(n_bins, 0.0);
    cnt.assign(n_bins, 0);
    if (stat == AggStat::maxmean) {
      sp.assign(n_bins, 0.0);
      sn.assign(n_bins, 0.0);
      cp.assign(n_bins, 0);
      cn.assign(n_bins, 0);
    }
    if (stat == AggStat::af) {
      zs.resize(n_bins);
      for (auto& v : zs) v.clear();
    }
    for (int j = 0; j < m; ++j) {
      const int b = bin[j];
      if (b < 0) continue;
      ++cnt[b];
      const double zj = std::max(-winsor, std::min(winsor, z[j]));
      if (stat == AggStat::stouffer) {
        num[b] += w[j] * zj;
        den2[b] += w[j] * w[j];
      } else if (stat == AggStat::maxmean) {
        if (zj > 0.0) { sp[b] += zj; ++cp[b]; }
        else if (zj < 0.0) { sn[b] += -zj; ++cn[b]; }
      } else {
        zs[b].push_back(zj);
      }
    }
    double s = 0.0;
    int n_nonempty = 0;
    for (int b = 0; b < n_bins; ++b) {
      if (cnt[b] == 0) continue;
      ++n_nonempty;
      if (stat == AggStat::stouffer) {
        s += num[b] / std::sqrt(std::max(den2[b], 1e-300));
      } else if (stat == AggStat::maxmean) {
        const double Sp = cp[b] > 0 ? sp[b] / cp[b] : 0.0;
        const double Sn = cn[b] > 0 ? sn[b] / cn[b] : 0.0;
        s += (Sp >= Sn) ? Sp : -Sn;
      } else {
        constexpr double P_MIN = 1e-300;
        const int k = zs[b].size();
        ppos.resize(k);
        pneg.resize(k);
        for (int a = 0; a < k; ++a) {
          ppos[a] = std::max(0.5 * std::erfc(zs[b][a] * M_SQRT1_2), P_MIN);
          pneg[a] = std::max(0.5 * std::erfc(-zs[b][a] * M_SQRT1_2), P_MIN);
        }
        const double Tp = af_one(ppos), Tn = af_one(pneg);
        s += (Tp >= Tn) ? Tp : -Tn;
      }
    }
    return n_nonempty > 0 ? s / n_nonempty : 0.0;
  }
};

// |T_b| against |T_obs|: draws that tie in exact arithmetic can differ by
// rounding once subject values are not integers (rho weights), so ties are
// taken within 1e-9 on the z scale (relative above 1)
inline int compare_abs(double tb, double tobs) {
  const double ax = std::fabs(tb), ao = std::fabs(tobs);
  if (std::fabs(ax - ao) <= 1e-9 * std::max(1.0, std::max(ax, ao))) return 0;
  return (ax > ao) ? 1 : -1;
}

struct NullAcc {
  double T_obs = 0.0;
  int64_t n = 0;
  double mean = 0.0, M2 = 0.0;
  int b_hits = 0, b_equal = 0;
  void digest(double x) {
    ++n;
    const double d = x - mean;
    mean += d / (double)n;
    M2 += d * (x - mean);
    const int c = compare_abs(x, T_obs);
    if (c > 0) ++b_hits;
    else if (c == 0) ++b_equal;
  }
};

double perm_p(int b_hits, int b_equal, int B, bool mid_p) {
  if (mid_p) {
    double p = ((double)b_hits + 0.5 * (double)b_equal) / (double)B;
    return (p <= 0.0) ? 0.5 / (double)B : p;
  }
  return (1.0 + (double)(b_hits + b_equal)) / (1.0 + (double)B);
}

// effect-size summary over peptides; frac columns only for signed effects
List delta_summary(const Vec& d, const Vec& d_adj, const Vec& w_norm,
                   bool is_signed) {
  const int m = d.size();
  double mx = 0.0, mx_adj = 0.0, fp = 0.0, fpw = 0.0;
  for (int j = 0; j < m; ++j) {
    mx = std::max(mx, std::fabs(d[j]));
    mx_adj = std::max(mx_adj, std::fabs(d_adj[j]));
    if (d[j] > 0.0) {
      fp += 1.0;
      fpw += w_norm[j];
    }
  }
  return List::create(
    _["max_delta"] = (m > 0) ? mx : NA_REAL,
    _["max_delta_adj"] = (m > 0) ? mx_adj : NA_REAL,
    _["frac_delta_pos"] = (is_signed && m > 0) ? fp / m : NA_REAL,
    _["frac_delta_pos_w"] = (is_signed && m > 0) ? fpw : NA_REAL
  );
}

// signed difference of the top- and bottom-scored levels for a 1-column
// design, else the largest range across levels of any component
double level_contrast(const Vec& M, int K, int L, const Vec& design,
                      int Ld) {
  if (Ld == 1 && L == 1) {
    int top = 0, bot = 0;
    for (int k = 1; k < K; ++k) {
      if (design[k] > design[top]) top = k;
      if (design[k] < design[bot]) bot = k;
    }
    return M[top] - M[bot];
  }
  double best = 0.0;
  for (int l = 0; l < L; ++l) {
    double lo = M[l], hi = M[l];
    for (int k = 1; k < K; ++k) {
      lo = std::min(lo, M[(size_t)k * L + l]);
      hi = std::max(hi, M[(size_t)k * L + l]);
    }
    best = std::max(best, hi - lo);
  }
  return best;
}


// orthonormal basis (n x (n - 1), row-major) of the vectors orthogonal to 1
Vec centred_basis(int n) {
  Vec C((size_t)n * (n - 1), 0.0);
  for (int j = 0; j < n - 1; ++j) {
    const double norm = std::sqrt((double)(j + 1) * (j + 2));
    for (int i = 0; i <= j; ++i) C[(size_t)i * (n - 1) + j] = 1.0 / norm;
    C[(size_t)(j + 1) * (n - 1) + j] = -(double)(j + 1) / norm;
  }
  return C;
}

// y = M x for a row-major r x c matrix
inline void matvec(const double* M, int r, int c, const double* x, double* y) {
  for (int a = 0; a < r; ++a) {
    double s = 0.0;
    for (int b = 0; b < c; ++b) s += M[(size_t)a * c + b] * x[b];
    y[a] = s;
  }
}

inline double quad(const double* M, int n, const double* x) {
  double s = 0.0;
  for (int a = 0; a < n; ++a)
    for (int b = 0; b < n; ++b) s += x[a] * M[(size_t)a * n + b] * x[b];
  return s;
}

// Shared objects of one subject labelling: per-level sums of the weights
// times the covariates (group term, pairs) and of the interaction design
struct LabShared {
  Vec WZ;                 // R x K x p:    sum_{i in k} w_i Z_i
  Vec W2ZZ;               // R x K x p x p: sum_{i in k} w_i^2 Z_i Z_i'
  Vec Mab;                // R x P x p x p: pinv of the pair's sum w Z Z'
  Vec nI, cred, SAr;      // interaction: size, centred reduced group design,
                          // sum of reduced design spreads
  Vec Q;                  // (Lg2 Lr) x Lr: (sum_k c_k (x) SA_k) SA^+, the
                          // pull of the estimated time effect on S
  Vec SArp, dfk;          // per level: SA_k^+ and the residual degrees of
                          // freedom sum_i rank(A_i) - rank(SA_k)
  bool inter_ok = false;
};

struct Stratum {
  bool do_group = false, do_time = false, do_inter = false, residualize = true;
  int K = 1, T = 1, N = 0, m = 0, Lg = 1, Lt = 1, LL = 1, Lg2 = 1, Lr = 1,
      LR = 1, P = 0, p = 1, R = 1, N_I = 0;
  Vec Gd, Xd, Cg, Ct;
  std::vector<int> pair_a, pair_b;
  Contrast cx;
  std::vector<int> obs_lab;
  Vec n_k, Z;
  std::vector<double> rho_vals;
  std::vector<int> rho_of;
  Vec wr, Mz, W2ZZ_tot;
  Vec gbar_b, Wg_b;
  Vec A, Ar, Arp, SumA_pinv, SumAr_pinv, SumA_full;
  std::vector<int> Ar_rank;
  std::vector<char> spread;

  // per-peptide CSR over positive subjects, and fixed quantities
  std::vector<int> off, csr_sub;
  Vec csr_v, csr_w, csr_D, csr_Dr, csr_dAd;
  std::vector<char> csr_spread;
  Vec tot, p0, alpha, SS, Wt, cst, beta, beta_r;
  std::vector<char> g_ok;
  std::vector<Vec> St_obs;

  // effect sizes
  Vec d_group, d_group_adj, d_time, d_inter, d_inter_adj, d_pair, d_pair_adj;

  // scratch
  Vec A1, A2, A1Z, A3Z, U, SDr, DAD, RS, Sred, Vred, rred, tmp, al, Mk, MS;

  void setup(const IntegerVector& sample_rows, const IntegerVector& sample_subject,
             const IntegerVector& group_code, const IntegerVector& time_code,
             const NumericVector& group_scores, const NumericVector& time_scores,
             const CharacterVector& terms, int n_group_levels, int n_time_levels,
             int n_words, int n_pep, const NumericVector& rho,
             const IntegerMatrix& pairs, const NumericMatrix* covariates,
             bool resid) {
    for (R_xlen_t k = 0; k < terms.size(); ++k) {
      const std::string tm = Rcpp::as<std::string>(terms[k]);
      if (tm == "group") do_group = true;
      else if (tm == "time") do_time = true;
      else if (tm == "interaction") do_inter = true;
      else stop("Unknown term: %s", tm);
    }
    residualize = resid;
    K = std::max(n_group_levels, 1);
    m = n_pep;
    if (m == 0) stop("No peptides.");
    if (rho.size() != 1 && rho.size() != m) stop("rho must have length 1 or one per peptide.");
    for (R_xlen_t i = 0; i < rho.size(); ++i) {
      if (!(rho[i] >= 0.0 && rho[i] <= 1.0)) stop("rho must be in [0, 1].");
    }
    for (R_xlen_t i = 0; i < group_code.size(); ++i) {
      const int g = group_code[i];
      if (g != NA_INTEGER && (g < 1 || g > K)) stop("group_code out of range.");
    }

    cx.init(n_words, sample_rows, sample_subject, time_code, &group_code, n_time_levels);
    T = cx.T;
    N = cx.n_sub;
    Gd = level_design(K, group_scores, Lg);
    Xd = level_design(T, time_scores, Lt);
    LL = Lt * Lt;
    // reduced coordinates for the interaction: an orthonormal basis of the
    // centred level space (identity for a one-column design)
    if (Lg == 1) { Cg = Vec{1.0}; Lg2 = 1; } else { Cg = centred_basis(K); Lg2 = K - 1; }
    if (Lt == 1) { Ct = Vec{1.0}; Lr = 1; } else { Ct = centred_basis(T); Lr = T - 1; }
    LR = Lr * Lr;

    obs_lab.assign(N, 0);
    n_k.assign(K, 0.0);
    for (int i = 0; i < N; ++i) {
      obs_lab[i] = group_code[cx.sub_orig[i] - 1] - 1;
      n_k[obs_lab[i]] += 1.0;
    }
    P = pairs.nrow();
    for (int q = 0; q < P; ++q) {
      pair_a.push_back(pairs(q, 0) - 1);
      pair_b.push_back(pairs(q, 1) - 1);
      if (pair_a.back() < 0 || pair_a.back() >= K || pair_b.back() < 0 || pair_b.back() >= K) {
        stop("pairs out of range.");
      }
    }

    // covariates of the group term: an intercept, plus the given columns
    p = 1 + (covariates ? covariates->ncol() : 0);
    if (covariates && covariates->nrow() != group_code.size()) {
      stop("covariates must have one row per subject.");
    }
    Z.assign((size_t)N * p, 1.0);
    for (int i = 0; i < N; ++i)
      for (int c = 1; c < p; ++c) {
        const double z = (*covariates)(cx.sub_orig[i] - 1, c - 1);
        if (!std::isfinite(z)) stop("covariates must be finite.");
        Z[(size_t)i * p + c] = z;
      }

    // between-subject weights w_i = n_i / (1 + (n_i - 1) rho), one set per
    // distinct rho (rho_by_stratum gives one per prevalence stratum), and the
    // weighted least-squares projection on Z for each
    rho_of.assign(m, 0);
    for (int j = 0; j < m; ++j) {
      const double rj = (rho.size() == 1) ? rho[0] : rho[j];
      int r = 0;
      while (r < (int)rho_vals.size() && rho_vals[r] != rj) ++r;
      if (r == (int)rho_vals.size()) rho_vals.push_back(rj);
      rho_of[j] = r;
    }
    R = rho_vals.size();
    wr.assign((size_t)R * N, 0.0);
    Mz.assign((size_t)R * p * p, 0.0);
    W2ZZ_tot.assign((size_t)R * p * p, 0.0);
    for (int r = 0; r < R; ++r) {
      Vec WZZ((size_t)p * p, 0.0);
      for (int i = 0; i < N; ++i) {
        const double w = cx.n_i[i] / (1.0 + (cx.n_i[i] - 1.0) * rho_vals[r]);
        wr[(size_t)r * N + i] = w;
        const double* zi = &Z[(size_t)i * p];
        for (int a = 0; a < p; ++a)
          for (int b = 0; b < p; ++b) {
            WZZ[(size_t)a * p + b] += w * zi[a] * zi[b];
            W2ZZ_tot[(size_t)r * p * p + a * p + b] += w * w * zi[a] * zi[b];
          }
      }
      const Vec M = psd_pinv_pow(WZZ, p, -1.0);
      std::copy(M.begin(), M.end(), Mz.begin() + (size_t)r * p * p);
    }

    // between-subject group design: fixed, the label multiset is fixed
    gbar_b.assign(Lg, 0.0);
    Wg_b.assign((size_t)Lg * Lg, 0.0);
    if (N >= 2) {
      for (int i = 0; i < N; ++i)
        for (int l = 0; l < Lg; ++l) gbar_b[l] += Gd[(size_t)obs_lab[i] * Lg + l] / N;
      Vec Gb((size_t)Lg * Lg, 0.0);
      for (int i = 0; i < N; ++i)
        for (int a = 0; a < Lg; ++a)
          for (int b = 0; b < Lg; ++b)
            Gb[a * Lg + b] += (Gd[(size_t)obs_lab[i] * Lg + a] - gbar_b[a]) *
              (Gd[(size_t)obs_lab[i] * Lg + b] - gbar_b[b]) / (N - 1.0);
      Wg_b = psd_pinv_pow(Gb, Lg, -0.5);
    }

    // within-subject design spread A_i = X_i' C X_i, and its reduced form
    A.assign((size_t)N * LL, 0.0);
    Ar.assign((size_t)N * LR, 0.0);
    Arp.assign((size_t)N * LR, 0.0);
    Ar_rank.assign(N, 0);
    spread.assign(N, 0);
    SumA_full.assign(LL, 0.0);
    Vec SumAr(LR, 0.0);
    for (int i = 0; i < N; ++i) {
      Vec xs(Lt, 0.0);
      for (int t = 0; t < T; ++t)
        for (int l = 0; l < Lt; ++l) xs[l] += Xd[(size_t)t * Lt + l] * cx.n_it[(size_t)i * T + t];
      double tr = 0.0;
      double* Ai = &A[(size_t)i * LL];
      for (int a = 0; a < Lt; ++a) {
        for (int b = 0; b < Lt; ++b) {
          double s = 0.0;
          for (int t = 0; t < T; ++t)
            s += Xd[(size_t)t * Lt + a] * Xd[(size_t)t * Lt + b] * cx.n_it[(size_t)i * T + t];
          Ai[a * Lt + b] = s - xs[a] * xs[b] / cx.n_i[i];
        }
        tr += Ai[a * Lt + a];
      }
      spread[i] = (tr > 1e-12) ? 1 : 0;
      // Ar = Ct' A Ct
      double* Ari = &Ar[(size_t)i * LR];
      for (int a = 0; a < Lr; ++a)
        for (int b = 0; b < Lr; ++b) {
          double s = 0.0;
          for (int u = 0; u < Lt; ++u)
            for (int v = 0; v < Lt; ++v)
              s += Ct[(size_t)u * Lr + a] * Ai[u * Lt + v] * Ct[(size_t)v * Lr + b];
          Ari[a * Lr + b] = s;
        }
      // Ar_i^+ and its rank: subject i's information per time direction
      {
        const Vec Ari_v(Ari, Ari + LR);
        const Vec pinv = psd_pinv_pow(Ari_v, Lr, -1.0);
        std::copy(pinv.begin(), pinv.end(), Arp.begin() + (size_t)i * LR);
        Vec vals, vecs;
        if (Lr == 1) {
          Ar_rank[i] = Ari[0] > 1e-12 ? 1 : 0;
        } else {
          sym_eigen(Ari_v, Lr, vals, vecs);
          double vmax = 0.0;
          for (double v : vals) vmax = std::max(vmax, v);
          int rk = 0;
          for (double v : vals) rk += (vmax > 0.0 && v > 1e-10 * vmax);
          Ar_rank[i] = rk;
        }
      }
      if (!spread[i]) continue;
      ++N_I;
      for (int e = 0; e < LL; ++e) SumA_full[e] += Ai[e];
      for (int e = 0; e < LR; ++e) SumAr[e] += Ari[e];
    }
    SumA_pinv = psd_pinv_pow(SumA_full, Lt, -1.0);
    SumAr_pinv = psd_pinv_pow(SumAr, Lr, -1.0);

    A1.assign(K, 0.0);
    A2.assign(K, 0.0);
    A1Z.assign((size_t)K * p, 0.0);
    A3Z.assign((size_t)K * p, 0.0);
    U.assign(std::max(Lg, Lt), 0.0);
    SDr.assign((size_t)K * Lr, 0.0);
    DAD.assign(K, 0.0);
    RS.assign(Lr, 0.0);
    Sred.assign((size_t)Lg2 * Lr, 0.0);
    Vred.assign((size_t)Lg2 * Lr * Lg2 * Lr, 0.0);
    rred.assign((size_t)Lg2 * Lr, 0.0);
    tmp.assign(std::max(p, Lr), 0.0);
    al.assign(p, 0.0);
    Mk.assign((size_t)Lg2 * Lr * Lr, 0.0);
    MS.assign((size_t)Lg2 * Lr * Lr, 0.0);
  }

  // per-peptide fixed quantities from the observed data, and effect sizes
  // (adj: inverse-propensity weight per subject, or empty)
  void observed(const RawVector& bitset_raw, int n_words, const IntegerVector& pep_cols,
                const NumericVector& adj_weight, const NumericVector& time_scores,
                const IntegerVector& group_code) {
    off.assign(m + 1, 0);
    tot.assign(m, 0.0);
    p0.assign(m, 0.0);
    alpha.assign((size_t)m * p, 0.0);
    SS.assign(m, 0.0);
    g_ok.assign(m, 0);
    Wt.assign((size_t)m * LL, 0.0);
    cst.assign((size_t)m * Lt, 0.0);
    beta.assign((size_t)m * Lt, 0.0);
    beta_r.assign((size_t)m * Lr, 0.0);
    St_obs.assign(m, Vec(T, 0.0));
    d_group.assign(m, 0.0);
    d_group_adj.assign(m, 0.0);
    d_time.assign(m, 0.0);
    d_inter.assign(m, 0.0);
    d_inter_adj.assign(m, 0.0);
    d_pair.assign((size_t)m * P, 0.0);
    d_pair_adj.assign((size_t)m * P, 0.0);

    const bool has_adj = adj_weight.size() > 0;
    if (has_adj && adj_weight.size() != group_code.size()) stop("adj_weight must have one value per subject.");
    Vec a_i(N, 1.0);
    if (has_adj) for (int i = 0; i < N; ++i) a_i[i] = adj_weight[cx.sub_orig[i] - 1];
    Vec asum_k(K, 0.0);
    for (int i = 0; i < N; ++i) asum_k[obs_lab[i]] += a_i[i];
    // per-group pinv of the (weighted) spread, for the interaction effect
    std::vector<Vec> SumA_k_pinv(K), SumA_k_adj_pinv(K);
    if (do_inter) {
      std::vector<Vec> sa(K, Vec(LL, 0.0)), saa(K, Vec(LL, 0.0));
      for (int i = 0; i < N; ++i) {
        if (!spread[i]) continue;
        for (int e = 0; e < LL; ++e) {
          sa[obs_lab[i]][e] += A[(size_t)i * LL + e];
          saa[obs_lab[i]][e] += a_i[i] * A[(size_t)i * LL + e];
        }
      }
      for (int k = 0; k < K; ++k) {
        SumA_k_pinv[k] = psd_pinv_pow(sa[k], Lt, -1.0);
        SumA_k_adj_pinv[k] = psd_pinv_pow(saa[k], Lt, -1.0);
      }
    }
    double tscale = 1.0;
    if (time_scores.size() > 0) {
      tscale = *std::max_element(time_scores.begin(), time_scores.end()) -
        *std::min_element(time_scores.begin(), time_scores.end());
    }

    Vec D(Lt), Dr(Lr), sumD(Lt), sumDr(Lr), Vt(LL), bz(p), b3z(p), mk(K), mk_adj(K),
        Dk((size_t)K * Lt), Dk_adj((size_t)K * Lt), bk((size_t)K * Lt), bk_adj((size_t)K * Lt);
    for (int j = 0; j < m; ++j) {
      const uint8_t* col = column_ptr(bitset_raw, n_words, pep_cols[j]);
      const int r = rho_of[j];
      const double* wj = &wr[(size_t)r * N];
      tot[j] = cx.count(col, n_words, St_obs[j]);
      p0[j] = cx.n_samples > 0 ? tot[j] / cx.n_samples : 0.0;

      double sv2 = 0.0;
      std::fill(bz.begin(), bz.end(), 0.0);
      std::fill(b3z.begin(), b3z.end(), 0.0);
      std::fill(sumD.begin(), sumD.end(), 0.0);
      std::fill(sumDr.begin(), sumDr.end(), 0.0);
      std::fill(Vt.begin(), Vt.end(), 0.0);
      std::fill(mk.begin(), mk.end(), 0.0);
      std::fill(mk_adj.begin(), mk_adj.end(), 0.0);
      std::fill(Dk.begin(), Dk.end(), 0.0);
      std::fill(Dk_adj.begin(), Dk_adj.end(), 0.0);
      for (int i : cx.touched) {
        const double* ci = &cx.cnt[(size_t)i * T];
        double c = 0.0;
        for (int t = 0; t < T; ++t) c += ci[t];
        const double ni = cx.n_i[i];
        const double ybar = c / ni;
        const double v = wj[i] * ybar;
        sv2 += v * v;
        const double* zi = &Z[(size_t)i * p];
        for (int a = 0; a < p; ++a) {
          bz[a] += v * zi[a];
          b3z[a] += wj[i] * v * zi[a];
        }
        mk[obs_lab[i]] += ybar;
        mk_adj[obs_lab[i]] += a_i[i] * ybar;
        csr_sub.push_back(i);
        csr_v.push_back(v);
        csr_w.push_back(wj[i]);
        csr_spread.push_back(spread[i]);
        for (int l = 0; l < Lt; ++l) {
          double s = 0.0, q = 0.0;
          for (int t = 0; t < T; ++t) {
            s += Xd[(size_t)t * Lt + l] * (ci[t] - cx.n_it[(size_t)i * T + t] * ybar);
            q += Xd[(size_t)t * Lt + l] * cx.n_it[(size_t)i * T + t] * ybar;
          }
          D[l] = spread[i] ? s : 0.0;
          // time const: sum_t X[t,l] n_it ybar_i, over every positive subject
          cst[(size_t)j * Lt + l] += q;
        }
        matvec_t(Ct.data(), Lt, Lr, D.data(), Dr.data());
        csr_D.insert(csr_D.end(), D.begin(), D.end());
        csr_Dr.insert(csr_Dr.end(), Dr.begin(), Dr.end());
        csr_dAd.push_back(spread[i] ? quad(&Arp[(size_t)i * LR], Lr, Dr.data()) : 0.0);
        if (!spread[i]) continue;
        const double f = c * (ni - c) / ni / (ni - 1.0);
        for (int e = 0; e < LL; ++e) Vt[e] += f * A[(size_t)i * LL + e];
        for (int l = 0; l < Lt; ++l) {
          sumD[l] += D[l];
          Dk[(size_t)obs_lab[i] * Lt + l] += D[l];
          Dk_adj[(size_t)obs_lab[i] * Lt + l] += a_i[i] * D[l];
        }
        for (int a = 0; a < Lr; ++a) sumDr[a] += Dr[a];
      }
      off[j + 1] = csr_sub.size();

      // group term: alpha = (sum w Z Z')^+ sum w ybar Z, and
      // SS = sum_i u_i^2 = sum v^2 - 2 alpha' sum w v Z + alpha' sum w^2 Z Z' alpha
      double* aj = &alpha[(size_t)j * p];
      matvec(&Mz[(size_t)r * p * p], p, p, bz.data(), aj);
      double ab3 = 0.0;
      for (int a = 0; a < p; ++a) ab3 += aj[a] * b3z[a];
      const double aWa = quad(&W2ZZ_tot[(size_t)r * p * p], p, aj);
      SS[j] = (N >= 2) ? sv2 - 2.0 * ab3 + aWa : 0.0;
      g_ok[j] = SS[j] > 1e-12 * (sv2 + aWa);

      // time whitening and main effect (within estimator)
      if (do_time) {
        const Vec W = psd_pinv_pow(Vt, Lt, -0.5);
        std::copy(W.begin(), W.end(), Wt.begin() + (size_t)j * LL);
      }
      double* bj = &beta[(size_t)j * Lt];
      matvec(SumA_pinv.data(), Lt, Lt, sumD.data(), bj);
      if (residualize) matvec(SumAr_pinv.data(), Lr, Lr, sumDr.data(), &beta_r[(size_t)j * Lr]);
      if (do_time) {
        if (Lt == 1) {
          d_time[j] = bj[0] * tscale;
        } else {
          double lo = bj[0], hi = bj[0];
          for (int a = 1; a < Lt; ++a) { lo = std::min(lo, bj[a]); hi = std::max(hi, bj[a]); }
          d_time[j] = hi - lo;
        }
      }
      if (do_inter) {
        // within estimate per group, contrasted across groups
        for (int k = 0; k < K; ++k)
          for (int a = 0; a < Lt; ++a) {
            double s = 0.0, sa = 0.0;
            for (int b = 0; b < Lt; ++b) {
              s += SumA_k_pinv[k][a * Lt + b] * Dk[(size_t)k * Lt + b];
              sa += SumA_k_adj_pinv[k][a * Lt + b] * Dk_adj[(size_t)k * Lt + b];
            }
            bk[(size_t)k * Lt + a] = s * tscale;
            bk_adj[(size_t)k * Lt + a] = sa * tscale;
          }
        d_inter[j] = level_contrast(bk, K, Lt, Gd, Lg);
        d_inter_adj[j] = level_contrast(bk_adj, K, Lt, Gd, Lg);
      }
      // group effect: subject-mean prevalence per level
      for (int k = 0; k < K; ++k) {
        mk[k] = n_k[k] > 0 ? mk[k] / n_k[k] : 0.0;
        mk_adj[k] = asum_k[k] > 0 ? mk_adj[k] / asum_k[k] : 0.0;
      }
      d_group[j] = level_contrast(mk, K, 1, Gd, Lg);
      d_group_adj[j] = level_contrast(mk_adj, K, 1, Gd, Lg);
      for (int q = 0; q < P; ++q) {
        d_pair[(size_t)q * m + j] = mk[pair_b[q]] - mk[pair_a[q]];
        d_pair_adj[(size_t)q * m + j] = mk_adj[pair_b[q]] - mk_adj[pair_a[q]];
      }
      cx.reset();
    }
  }

  // y = M' x for a row-major r x c matrix M (y has c entries)
  static void matvec_t(const double* M, int r, int c, const double* x, double* y) {
    for (int b = 0; b < c; ++b) {
      double s = 0.0;
      for (int a = 0; a < r; ++a) s += M[(size_t)a * c + b] * x[a];
      y[b] = s;
    }
  }

  void shared(const int* lab, LabShared& sh) const {
    if (do_group || P > 0) {
      sh.WZ.assign((size_t)R * K * p, 0.0);
      if (P > 0) sh.W2ZZ.assign((size_t)R * K * p * p, 0.0);
      Vec WZZ(P > 0 ? (size_t)R * K * p * p : 0, 0.0);
      for (int r = 0; r < R; ++r)
        for (int i = 0; i < N; ++i) {
          const double w = wr[(size_t)r * N + i];
          const double* zi = &Z[(size_t)i * p];
          const size_t k = lab[i];
          double* wz = &sh.WZ[((size_t)r * K + k) * p];
          for (int a = 0; a < p; ++a) wz[a] += w * zi[a];
          if (P == 0) continue;
          double* wzz = &WZZ[((size_t)r * K + k) * p * p];
          double* w2zz = &sh.W2ZZ[((size_t)r * K + k) * p * p];
          for (int a = 0; a < p; ++a)
            for (int b = 0; b < p; ++b) {
              wzz[a * p + b] += w * zi[a] * zi[b];
              w2zz[a * p + b] += w * w * zi[a] * zi[b];
            }
        }
      if (P > 0) {
        sh.Mab.assign((size_t)R * P * p * p, 0.0);
        Vec S((size_t)p * p);
        for (int r = 0; r < R; ++r)
          for (int q = 0; q < P; ++q) {
            const double* wa = &WZZ[((size_t)r * K + pair_a[q]) * p * p];
            const double* wb = &WZZ[((size_t)r * K + pair_b[q]) * p * p];
            for (int e = 0; e < p * p; ++e) S[e] = wa[e] + wb[e];
            const Vec M = psd_pinv_pow(S, p, -1.0);
            std::copy(M.begin(), M.end(), sh.Mab.begin() + ((size_t)r * P + q) * p * p);
          }
      }
    }
    if (do_inter) {
      sh.inter_ok = false;
      sh.nI.assign(K, 0.0);
      sh.cred.assign((size_t)K * Lg2, 0.0);
      sh.SAr.assign((size_t)K * LR, 0.0);
      Vec gbar(Lg, 0.0);
      for (int i = 0; i < N; ++i) {
        if (!spread[i]) continue;
        sh.nI[lab[i]] += 1.0;
        for (int l = 0; l < Lg; ++l) gbar[l] += Gd[(size_t)lab[i] * Lg + l];
        const double* Ari = &Ar[(size_t)i * LR];
        double* sa = &sh.SAr[(size_t)lab[i] * LR];
        for (int e = 0; e < LR; ++e) sa[e] += Ari[e];
      }
      int levels = 0;
      for (int k = 0; k < K; ++k) levels += sh.nI[k] > 0;
      if (N_I < 2 || levels < 2) return;
      // per level: SA_k^+ and sum_i rank(A_i) - rank(SA_k)
      sh.SArp.assign((size_t)K * LR, 0.0);
      sh.dfk.assign(K, 0.0);
      for (int i = 0; i < N; ++i) if (spread[i]) sh.dfk[lab[i]] += Ar_rank[i];
      for (int k = 0; k < K; ++k) {
        if (sh.nI[k] <= 0) continue;
        const Vec SAk(sh.SAr.begin() + (size_t)k * LR, sh.SAr.begin() + (size_t)(k + 1) * LR);
        const Vec pinv = psd_pinv_pow(SAk, Lr, -1.0);
        std::copy(pinv.begin(), pinv.end(), sh.SArp.begin() + (size_t)k * LR);
        int rk = 0;
        if (Lr == 1) {
          rk = SAk[0] > 1e-12 ? 1 : 0;
        } else {
          Vec vals, vecs;
          sym_eigen(SAk, Lr, vals, vecs);
          double vmax = 0.0;
          for (double v : vals) vmax = std::max(vmax, v);
          for (double v : vals) rk += (vmax > 0.0 && v > 1e-10 * vmax);
        }
        sh.dfk[k] -= rk;
      }
      for (int l = 0; l < Lg; ++l) gbar[l] /= N_I;
      Vec c(Lg);
      for (int k = 0; k < K; ++k) {
        for (int l = 0; l < Lg; ++l) c[l] = Gd[(size_t)k * Lg + l] - gbar[l];
        matvec_t(Cg.data(), Lg, Lg2, c.data(), &sh.cred[(size_t)k * Lg2]);
      }
      // S = sum_i (c_i (x) I - Q) D_i with Q = GA SA^+, GA = sum_k c_k (x) SA_k:
      // zero when every group has the same visit pattern
      const int dI = Lg2 * Lr;
      sh.Q.assign((size_t)dI * Lr, 0.0);
      if (residualize) {
        Vec GA((size_t)dI * Lr, 0.0);
        for (int k = 0; k < K; ++k)
          for (int l = 0; l < Lg2; ++l) {
            const double ck = sh.cred[(size_t)k * Lg2 + l];
            if (ck == 0.0) continue;
            for (int a = 0; a < Lr; ++a)
              for (int b = 0; b < Lr; ++b)
                GA[(size_t)(l * Lr + a) * Lr + b] += ck * sh.SAr[(size_t)k * LR + a * Lr + b];
          }
        for (int e = 0; e < dI; ++e)
          for (int b = 0; b < Lr; ++b) {
            double s = 0.0;
            for (int a = 0; a < Lr; ++a) s += GA[(size_t)e * Lr + a] * SumAr_pinv[a * Lr + b];
            sh.Q[(size_t)e * Lr + b] = s;
          }
      }
      sh.inter_ok = true;
    }
  }

  // group term (Lg x m), interaction (Lg Lt x m) and pairs (P x m) under a
  // labelling, component-major
  void between(const int* lab, const LabShared& sh, Vec& r_group, Vec& r_inter,
               Vec& z_pair) {
    const int dI = Lg2 * Lr;
    for (int j = 0; j < m; ++j) {
      const int r = rho_of[j];
      const double* aj = &alpha[(size_t)j * p];
      if (do_group || P > 0) {
        std::fill(A1.begin(), A1.end(), 0.0);
        std::fill(A2.begin(), A2.end(), 0.0);
        if (P > 0) {
          std::fill(A1Z.begin(), A1Z.end(), 0.0);
          std::fill(A3Z.begin(), A3Z.end(), 0.0);
        }
        for (int e = off[j]; e < off[j + 1]; ++e) {
          const int i = csr_sub[e], k = lab[i];
          const double v = csr_v[e];
          A1[k] += v;
          A2[k] += v * v;
          if (P == 0) continue;
          const double* zi = &Z[(size_t)i * p];
          for (int a = 0; a < p; ++a) {
            A1Z[(size_t)k * p + a] += v * zi[a];
            A3Z[(size_t)k * p + a] += csr_w[e] * v * zi[a];
          }
        }
      }
      if (do_group) {
        // U_l = sum_k (G[k,l] - gbar_l) sum_{i in k} u_i,
        // sum_{i in k} u_i = A1_k - (sum_{i in k} w_i Z_i)' alpha
        for (int l = 0; l < Lg; ++l) {
          double s = 0.0;
          for (int k = 0; k < K; ++k) {
            const double* wz = &sh.WZ[((size_t)r * K + k) * p];
            double fit = 0.0;
            for (int a = 0; a < p; ++a) fit += wz[a] * aj[a];
            s += (Gd[(size_t)k * Lg + l] - gbar_b[l]) * (A1[k] - fit);
          }
          U[l] = s;
        }
        for (int a = 0; a < Lg; ++a) {
          double z = 0.0;
          if (g_ok[j]) {
            for (int b = 0; b < Lg; ++b) z += Wg_b[a * Lg + b] * U[b];
            z /= std::sqrt(SS[j]);
          }
          r_group[(size_t)a * m + j] = z;
        }
      }
      // pairs: the group statistic refitted on the pair's own subjects
      for (int q = 0; q < P; ++q) {
        const int a = pair_a[q], b = pair_b[q];
        const double na = n_k[a], nb = n_k[b], Nab = na + nb;
        double z = 0.0;
        if (na > 0 && nb > 0 && Nab >= 2) {
          for (int c = 0; c < p; ++c) tmp[c] = A1Z[(size_t)a * p + c] + A1Z[(size_t)b * p + c];
          matvec(&sh.Mab[((size_t)r * P + q) * p * p], p, p, tmp.data(), al.data());
          const double* wzb = &sh.WZ[((size_t)r * K + b) * p];
          double fitb = 0.0, a3 = 0.0;
          for (int c = 0; c < p; ++c) {
            fitb += wzb[c] * al[c];
            a3 += al[c] * (A3Z[(size_t)a * p + c] + A3Z[(size_t)b * p + c]);
          }
          const double* w2a = &sh.W2ZZ[((size_t)r * K + a) * p * p];
          const double* w2b = &sh.W2ZZ[((size_t)r * K + b) * p * p];
          const double aWa = quad(w2a, p, al.data()) + quad(w2b, p, al.data());
          const double s2 = A2[a] + A2[b];
          const double ss = s2 - 2.0 * a3 + aWa;
          if (ss > 1e-12 * (s2 + aWa)) {
            z = (A1[b] - fitb) / std::sqrt(na * nb / (Nab * (Nab - 1.0)) * ss);
          }
        }
        z_pair[(size_t)q * m + j] = z;
      }
      if (do_inter) interaction(j, lab, sh, r_inter, dI);
    }
  }

  // Robust interaction of peptide j: residuals R_i = D_i - A_i beta in
  // reduced time coordinates. S = sum_k c_k (x) RS_k = sum_i M_k(i) D_i with
  // M_k = c_k (x) I - Q, so V = sum_k M_k Cov(sum_{i in k} D_i) M_k'. Within a
  // group, Cov(D_i) = s2_k A_i (the subject's own design spread), with s2_k
  // the group's residual variance per unit of within-subject information
  // around its own time effect: a Welch-type variance (one scale per group)
  // whose shape comes from the design, so it stays stable when subjects see
  // different visits. r = (Cg (x) Ct) V^{+1/2} S.
  void interaction(int j, const int* lab, const LabShared& sh, Vec& r_inter, int dI) {
    const int Lfull = Lg * Lt;
    if (!sh.inter_ok) {
      for (int e = 0; e < Lfull; ++e) r_inter[(size_t)e * m + j] = 0.0;
      return;
    }
    const double* br = &beta_r[(size_t)j * Lr];
    std::fill(SDr.begin(), SDr.end(), 0.0);
    std::fill(DAD.begin(), DAD.end(), 0.0);
    for (int e = off[j]; e < off[j + 1]; ++e) {
      if (!csr_spread[e]) continue;
      const int k = lab[csr_sub[e]];
      const double* De = &csr_Dr[(size_t)e * Lr];
      for (int a = 0; a < Lr; ++a) SDr[(size_t)k * Lr + a] += De[a];
      DAD[k] += csr_dAd[e];
    }
    std::fill(Sred.begin(), Sred.end(), 0.0);
    std::fill(Vred.begin(), Vred.end(), 0.0);
    double scale = 0.0;
    for (int k = 0; k < K; ++k) {
      if (sh.nI[k] <= 0) continue;
      const double* sa = &sh.SAr[(size_t)k * LR];
      const double* sd = &SDr[(size_t)k * Lr];
      // RS_k = sum_{i in k} R_i = SD_k - SA_k beta
      for (int a = 0; a < Lr; ++a) {
        double s = 0.0;
        for (int b = 0; b < Lr; ++b) s += sa[a * Lr + b] * br[b];
        RS[a] = sd[a] - s;
      }
      const double* ck = &sh.cred[(size_t)k * Lg2];
      for (int l = 0; l < Lg2; ++l)
        for (int a = 0; a < Lr; ++a) Sred[(size_t)l * Lr + a] += ck[l] * RS[a];
      // s2_k = (sum_i D_i' A_i^+ D_i - SD_k' SA_k^+ SD_k) / df_k
      const double df = sh.dfk[k];
      if (df <= 0.0) continue;
      const double s2 = std::max(DAD[k] - quad(&sh.SArp[(size_t)k * LR], Lr, sd), 0.0) / df;
      if (s2 <= 0.0) continue;
      // V += s2 M_k SA_k M_k'
      double m2 = 0.0, tr = 0.0;
      for (int a = 0; a < Lr; ++a) tr += sa[a * Lr + a];
      for (int l = 0; l < Lg2; ++l)
        for (int a = 0; a < Lr; ++a)
          for (int b = 0; b < Lr; ++b) {
            const size_t e = (size_t)(l * Lr + a);
            const double mv = (a == b ? ck[l] : 0.0) - sh.Q[e * Lr + b];
            Mk[e * Lr + b] = mv;
            m2 += mv * mv;
          }
      scale += s2 * m2 * tr;
      for (int e = 0; e < dI; ++e)
        for (int b = 0; b < Lr; ++b) {
          double s = 0.0;
          for (int a = 0; a < Lr; ++a) s += Mk[(size_t)e * Lr + a] * sa[a * Lr + b];
          MS[(size_t)e * Lr + b] = s;
        }
      for (int e = 0; e < dI; ++e)
        for (int f = 0; f < dI; ++f) {
          double s = 0.0;
          for (int b = 0; b < Lr; ++b) s += MS[(size_t)e * Lr + b] * Mk[(size_t)f * Lr + b];
          Vred[(size_t)e * dI + f] += s2 * s;
        }
    }
    const Vec W = psd_pinv_pow(Vred, dI, -0.5, scale);
    matvec(W.data(), dI, dI, Sred.data(), rred.data());
    // back to the full level coordinates
    for (int l = 0; l < Lg; ++l)
      for (int a = 0; a < Lt; ++a) {
        double s = 0.0;
        for (int l2 = 0; l2 < Lg2; ++l2)
          for (int a2 = 0; a2 < Lr; ++a2)
            s += Cg[(size_t)l * Lg2 + l2] * Ct[(size_t)a * Lr + a2] * rred[(size_t)l2 * Lr + a2];
        r_inter[(size_t)(l * Lt + a) * m + j] = s;
      }
  }

  // time term of peptide j from per-timepoint counts (observed or relabelled)
  void time_r(int j, const Vec& counts, Vec& r_time) {
    for (int l = 0; l < Lt; ++l) {
      double s = -cst[(size_t)j * Lt + l];
      for (int u = 0; u < T; ++u) s += Xd[(size_t)u * Lt + l] * counts[u];
      U[l] = s;
    }
    const double* W = &Wt[(size_t)j * LL];
    for (int a = 0; a < Lt; ++a) {
      double r = 0.0;
      for (int b = 0; b < Lt; ++b) r += W[a * Lt + b] * U[b];
      r_time[(size_t)a * m + j] = r;
    }
  }
};

} // namespace

// Per-peptide score z of each term for the observed labels.
//
// sample_rows, sample_subject and time_code are per sample (same length);
// group_code is per subject (indexed by sample_subject). Samples with an NA
// code, or whose subject has an NA group_code, are left out. Empty
// group_scores / time_scores mean omnibus over the levels (a single 0/1
// contrast when there are two). rho sets the between-subject weights
// w_i = n_i / (1 + (n_i - 1) rho). covariates (one row per subject, as
// group_code) adjust the group term; NULL means an intercept only.
//
// [[Rcpp::export]]
Rcpp::List cpp_delta2_score(const Rcpp::RawVector& bitset_raw,
                            int n_words,
                            const Rcpp::IntegerVector& pep_cols,
                            const Rcpp::IntegerVector& sample_rows,
                            const Rcpp::IntegerVector& sample_subject,
                            const Rcpp::IntegerVector& group_code,
                            const Rcpp::IntegerVector& time_code,
                            const Rcpp::NumericVector& group_scores,
                            const Rcpp::NumericVector& time_scores,
                            const Rcpp::CharacterVector& terms,
                            int n_group_levels, int n_time_levels,
                            double rho,
                            Rcpp::Nullable<Rcpp::NumericMatrix> covariates = R_NilValue) {
  if (n_words <= 0) stop("n_words must be positive.");
  NumericMatrix cov_m;
  if (covariates.isNotNull()) cov_m = NumericMatrix(covariates.get());
  Stratum st;
  st.setup(sample_rows, sample_subject, group_code, time_code, group_scores,
           time_scores, terms, n_group_levels, n_time_levels, n_words,
           pep_cols.size(), NumericVector::create(rho), IntegerMatrix(0, 2),
           covariates.isNotNull() ? &cov_m : nullptr, true);
  st.observed(bitset_raw, n_words, pep_cols, NumericVector(0), time_scores, group_code);
  const int m = st.m;

  Vec r_group((size_t)st.Lg * m), r_time((size_t)st.Lt * m),
      r_inter((size_t)st.Lg * st.Lt * m), z_pair;
  if (st.do_group || st.do_inter) {
    LabShared sh;
    st.shared(st.obs_lab.data(), sh);
    st.between(st.obs_lab.data(), sh, r_group, r_inter, z_pair);
  }
  if (st.do_time) for (int j = 0; j < m; ++j) st.time_r(j, st.St_obs[j], r_time);

  auto as_matrix = [&](const Vec& r, int L) {
    NumericMatrix out(m, L);
    for (int l = 0; l < L; ++l)
      for (int j = 0; j < m; ++j) out(j, l) = r[(size_t)l * m + j];
    return out;
  };
  NumericMatrix St_out(m, st.T);
  for (int j = 0; j < m; ++j)
    for (int t = 0; t < st.T; ++t) St_out(j, t) = st.St_obs[j][t];

  List out = List::create(
    _["tot"] = NumericVector(st.tot.begin(), st.tot.end()),
    _["S_time"] = St_out,
    _["n_subjects"] = st.N,
    _["n_samples"] = st.cx.n_samples,
    _["rho"] = rho
  );
  if (st.do_group) out["group"] = as_matrix(r_group, st.Lg);
  if (st.do_time) out["time"] = as_matrix(r_time, st.Lt);
  if (st.do_inter) {
    out["interaction"] = as_matrix(r_inter, st.Lg * st.Lt);
    out["n_subjects_interaction"] = st.N_I;
  }
  return out;
}

// Pooled statistics, permutation null and effect sizes of one stratum.
//
// Subjects are indexed by sample_subject into group_code; group_draws holds
// one resampled labelling per column (rows as group_code) and time_draws one
// within-subject relabelling of the time codes per column (rows as
// sample_rows). Draws must keep the observed label multiset. pairs lists the
// level pairs (1-based, group1 < group2) for the pairwise tests; adj_weight
// is each subject's inverse-propensity weight for the adjusted effect sizes
// (empty = unweighted). rho is one value or one per peptide. covariates
// (one row per subject, or NULL) adjust the group term and pairs.
// residualize = false drops the interaction's time main-effect correction
// (validation only: it shows what the correction is for). With
// standardize = "two_pass" each peptide's value is standardised by its mean
// and SD over the observed labels and the draws (a first pass over the same
// draws) before pooling.
//
// [[Rcpp::export]]
Rcpp::List cpp_delta2_contrast(const Rcpp::RawVector& bitset_raw,
                               int n_words,
                               const Rcpp::IntegerVector& pep_cols,
                               const Rcpp::IntegerVector& sample_rows,
                               const Rcpp::IntegerVector& sample_subject,
                               const Rcpp::IntegerVector& group_code,
                               const Rcpp::IntegerVector& time_code,
                               const Rcpp::NumericVector& group_scores,
                               const Rcpp::NumericVector& time_scores,
                               const Rcpp::CharacterVector& terms,
                               int n_group_levels, int n_time_levels,
                               const Rcpp::NumericVector& rho,
                               const Rcpp::IntegerMatrix& group_draws,
                               const Rcpp::IntegerMatrix& time_draws,
                               const Rcpp::IntegerMatrix& pairs,
                               const Rcpp::NumericVector& adj_weight,
                               const std::string& weight_mode,
                               const std::string& aggregate_stat,
                               const std::string& perm_method,
                               const Rcpp::NumericVector& strat_bins,
                               double winsor_z,
                               double min_m_eff,
                               Rcpp::Nullable<Rcpp::NumericMatrix> covariates = R_NilValue,
                               bool residualize = true,
                               const std::string& standardize = "none") {
  if (standardize != "none" && standardize != "two_pass") {
    stop("standardize must be \"none\" or \"two_pass\".");
  }
  const bool two_pass = (standardize == "two_pass");
  NumericMatrix cov_m;
  if (covariates.isNotNull()) cov_m = NumericMatrix(covariates.get());
  Stratum st;
  st.setup(sample_rows, sample_subject, group_code, time_code, group_scores,
           time_scores, terms, n_group_levels, n_time_levels, n_words,
           pep_cols.size(), rho, pairs, covariates.isNotNull() ? &cov_m : nullptr,
           residualize);
  const int m = st.m, P = st.P, N = st.N, T = st.T, K = st.K;
  const bool need_group_draws = st.do_group || st.do_inter || P > 0;

  int B = 0;
  if (need_group_draws) {
    if (group_draws.nrow() != group_code.size()) stop("group_draws must have one row per subject.");
    B = group_draws.ncol();
  }
  if (st.do_time) {
    if (time_draws.nrow() != sample_rows.size()) stop("time_draws must have one row per sample.");
    if (B > 0 && time_draws.ncol() != B) stop("group_draws and time_draws must have the same B.");
    B = time_draws.ncol();
  }
  if (B < 1) stop("Need at least one draw.");

  st.observed(bitset_raw, n_words, pep_cols, adj_weight, time_scores, group_code);

  // time-label masks: kept samples (bitset row, position)
  std::vector<int> kept_pos, kept_row;
  for (int s = 0; s < sample_rows.size(); ++s) {
    const int r = sample_rows[s];
    if (r == NA_INTEGER || r < 1 || r > n_words * 64) continue;
    if (st.cx.slot_sub[r - 1] < 0) continue;
    kept_pos.push_back(s);
    kept_row.push_back(r - 1);
  }

  // --- peptide weights and m_eff ---------------------------------------------
  Pooler pool;
  pool.stat = parse_agg(aggregate_stat);
  pool.winsor = winsor_z;
  pool.bin = prevalence_bins(strat_bins, st.p0, pool.n_bins);
  pool.w.assign(m, 1.0);
  for (int j = 0; j < m; ++j) {
    if (weight_mode == "n_eff_sqrt") {
      pool.w[j] = std::sqrt(std::max(st.tot[j], 1e-12));
    } else if (weight_mode == "se_invvar") {
      const double q = st.p0[j] * (1.0 - st.p0[j]);
      pool.w[j] = (q > 0.0) ? 1.0 / std::max(std::sqrt(q), 1e-6) : 0.0;
    }
  }
  double sw = 0.0;
  for (double x : pool.w) sw += x;
  Vec w_norm(m);
  double sw2 = 0.0;
  for (int j = 0; j < m; ++j) {
    w_norm[j] = pool.w[j] / std::max(sw, 1e-12);
    sw2 += w_norm[j] * w_norm[j];
  }
  const double m_eff = 1.0 / std::max(sw2, 1e-12);
  if (min_m_eff > 0.0 && m_eff < min_m_eff) {
    return List::create(_["skipped"] = true, _["m_eff"] = m_eff);
  }

  // --- statistic of one draw --------------------------------------------------
  const int Lg = st.Lg, Lt = st.Lt;
  Vec r_group((size_t)Lg * m), r_time((size_t)Lt * m), r_inter((size_t)Lg * Lt * m),
      z_pair((size_t)P * m), cnt_u(T);

  // pooled T of a term: pooled component(s) E, T = E (df 1) or sum E^2
  auto pooled = [&](const Vec& r, int L, Vec* E) {
    double t2 = 0.0, e0 = 0.0;
    if (E) E->assign(L, 0.0);
    for (int l = 0; l < L; ++l) {
      const double e = pool(&r[(size_t)l * m]);
      if (E) (*E)[l] = e;
      if (l == 0) e0 = e;
      t2 += e * e;
    }
    return (L == 1) ? e0 : t2;
  };

  // raw per-peptide values of one draw (b = -1: the observed labels) into
  // r_group / r_inter / z_pair and r_time
  LabShared sh;
  std::vector<int> lab(N);
  std::vector<uint64_t> G((size_t)T * n_words);
  const uint8_t* bits = reinterpret_cast<const uint8_t*>(RAW(bitset_raw));
  auto fill = [&](int b) {
    if (need_group_draws) {
      if (b < 0) {
        lab = st.obs_lab;
      } else {
        for (int i = 0; i < N; ++i) {
          const int l = group_draws(st.cx.sub_orig[i] - 1, b);
          if (l == NA_INTEGER || l < 1 || l > K) stop("group_draws out of range.");
          lab[i] = l - 1;
        }
      }
      st.shared(lab.data(), sh);
      st.between(lab.data(), sh, r_group, r_inter, z_pair);
    }
    if (!st.do_time) return;
    if (b < 0) {
      for (int j = 0; j < m; ++j) st.time_r(j, st.St_obs[j], r_time);
      return;
    }
    // shared object: the mask of samples relabelled to each timepoint
    std::fill(G.begin(), G.end(), 0ULL);
    for (size_t q = 0; q < kept_pos.size(); ++q) {
      const int u = time_draws(kept_pos[q], b);
      if (u == NA_INTEGER || u < 1 || u > T) stop("time_draws out of range.");
      const int z = kept_row[q];
      G[(size_t)(u - 1) * n_words + (z >> 6)] |= (1ULL << (z & 63));
    }
    for (int j = 0; j < m; ++j) {
      const uint8_t* col = bits + (size_t)(pep_cols[j] - 1) * n_words * 8;
      for (int u = 0; u < T; ++u) {
        double c = 0.0;
        for (int w = 0; w < n_words; ++w) {
          c += __builtin_popcountll(load_word(col + (size_t)w * 8) & G[(size_t)u * n_words + w]);
        }
        cnt_u[u] = c;
      }
      st.time_r(j, cnt_u, r_time);
    }
  };

  // two_pass: each entry's mean and SD over the observed labels and the B
  // draws (first pass), then every draw is replayed standardised; an entry
  // constant over the draws carries no information and is set to 0
  std::vector<Vec*> bufs;
  if (st.do_group) bufs.push_back(&r_group);
  if (st.do_inter) bufs.push_back(&r_inter);
  if (P > 0) bufs.push_back(&z_pair);
  if (st.do_time) bufs.push_back(&r_time);
  std::vector<Vec> mom_mean(bufs.size()), mom_sd(bufs.size());
  if (two_pass) {
    for (size_t q = 0; q < bufs.size(); ++q) {
      mom_mean[q].assign(bufs[q]->size(), 0.0);
      mom_sd[q].assign(bufs[q]->size(), 0.0);
    }
    for (int b = -1; b < B; ++b) {
      fill(b);
      const double n = (double)(b + 2);
      for (size_t q = 0; q < bufs.size(); ++q) {
        const Vec& r = *bufs[q];
        Vec& mu = mom_mean[q];
        Vec& M2 = mom_sd[q];
        for (size_t e = 0; e < r.size(); ++e) {
          const double d = r[e] - mu[e];
          mu[e] += d / n;
          M2[e] += d * (r[e] - mu[e]);
        }
      }
    }
    for (size_t q = 0; q < bufs.size(); ++q) {
      for (double& v : mom_sd[q]) v = std::sqrt(std::max(v, 0.0) / (double)B);
    }
  }
  auto standardise = [&]() {
    if (!two_pass) return;
    for (size_t q = 0; q < bufs.size(); ++q) {
      Vec& r = *bufs[q];
      const Vec& mu = mom_mean[q];
      const Vec& sd = mom_sd[q];
      for (size_t e = 0; e < r.size(); ++e) {
        r[e] = (sd[e] > 1e-12 * std::max(1.0, std::fabs(mu[e]))) ? (r[e] - mu[e]) / sd[e] : 0.0;
      }
    }
  };

  NullAcc acc_group, acc_time, acc_inter;
  std::vector<NullAcc> acc_pair(P);
  Vec E_group, E_time, E_inter;

  // observed
  fill(-1);
  standardise();
  if (st.do_group) acc_group.T_obs = pooled(r_group, Lg, &E_group);
  if (st.do_inter) acc_inter.T_obs = pooled(r_inter, Lg * Lt, &E_inter);
  for (int q = 0; q < P; ++q) acc_pair[q].T_obs = pool(&z_pair[(size_t)q * m]);
  if (st.do_time) acc_time.T_obs = pooled(r_time, Lt, &E_time);

  // draws (with two_pass, a replay of the draws of the first pass)
  Vec max_pair(B, 0.0);
  for (int b = 0; b < B; ++b) {
    fill(b);
    standardise();
    if (need_group_draws) {
      if (st.do_group) acc_group.digest(pooled(r_group, Lg, nullptr));
      if (st.do_inter) acc_inter.digest(pooled(r_inter, Lg * Lt, nullptr));
      double mx = 0.0;
      for (int q = 0; q < P; ++q) {
        const double tp = pool(&z_pair[(size_t)q * m]);
        acc_pair[q].digest(tp);
        mx = std::max(mx, std::fabs(tp));
      }
      max_pair[b] = mx;
    }
    if (st.do_time) acc_time.digest(pooled(r_time, Lt, nullptr));
  }

  const bool mid_p = (perm_method == "mid_p");
  auto term_out = [&](const NullAcc& acc, const Vec& E, const Vec& d, const Vec& d_adj,
                      bool is_signed) {
    double sd = NA_REAL;
    if (acc.n > 1) {
      double v = acc.M2 / (double)(acc.n - 1);
      if (v < 0.0 && v > -1e-15) v = 0.0;
      if (v >= 0.0) sd = std::sqrt(v);
    }
    List out = List::create(
      _["T_obs"] = acc.T_obs,
      _["T_null_mean"] = acc.n > 0 ? acc.mean : NA_REAL,
      _["T_null_sd"] = sd,
      _["b"] = acc.b_hits,
      _["p_perm"] = perm_p(acc.b_hits, acc.b_equal, B, mid_p),
      _["E"] = NumericVector(E.begin(), E.end())
    );
    const List ds = delta_summary(d, d_adj, w_norm, is_signed);
    for (const char* nm : {"max_delta", "max_delta_adj", "frac_delta_pos", "frac_delta_pos_w"}) {
      out[nm] = ds[nm];
    }
    return out;
  };

  List term_list;
  if (st.do_group) term_list["group"] = term_out(acc_group, E_group, st.d_group, st.d_group_adj, Lg == 1);
  if (st.do_time) term_list["time"] = term_out(acc_time, E_time, st.d_time, st.d_time, Lt == 1);
  if (st.do_inter) {
    term_list["interaction"] = term_out(acc_inter, E_inter, st.d_inter, st.d_inter_adj,
                                        Lg == 1 && Lt == 1);
  }

  List pair_out;
  if (P > 0) {
    NumericVector T_obs(P), p_perm(P), p_maxT(P), mx(P), mx_adj(P), fp(P), fpw(P);
    IntegerVector bb(P);
    for (int q = 0; q < P; ++q) {
      const NullAcc& a = acc_pair[q];
      T_obs[q] = a.T_obs;
      bb[q] = a.b_hits;
      p_perm[q] = perm_p(a.b_hits, a.b_equal, B, mid_p);
      int gt = 0, eq = 0;
      for (int b = 0; b < B; ++b) {
        const int c = compare_abs(max_pair[b], a.T_obs);
        if (c > 0) ++gt;
        else if (c == 0) ++eq;
      }
      p_maxT[q] = perm_p(gt, eq, B, mid_p);
      const Vec d(st.d_pair.begin() + (size_t)q * m, st.d_pair.begin() + (size_t)(q + 1) * m);
      const Vec da(st.d_pair_adj.begin() + (size_t)q * m, st.d_pair_adj.begin() + (size_t)(q + 1) * m);
      const List ds = delta_summary(d, da, w_norm, true);
      mx[q] = ds["max_delta"];
      mx_adj[q] = ds["max_delta_adj"];
      fp[q] = ds["frac_delta_pos"];
      fpw[q] = ds["frac_delta_pos_w"];
    }
    pair_out = List::create(
      _["T_obs"] = T_obs, _["p_perm"] = p_perm, _["p_maxT"] = p_maxT, _["b"] = bb,
      _["max_delta"] = mx, _["max_delta_adj"] = mx_adj,
      _["frac_delta_pos"] = fp, _["frac_delta_pos_w"] = fpw
    );
  }

  return List::create(
    _["skipped"] = false,
    _["n_peptides_used"] = m,
    _["m_eff"] = m_eff,
    _["n_subjects"] = N,
    _["n_samples"] = st.cx.n_samples,
    _["n_subjects_interaction"] = st.N_I,
    _["B"] = B,
    _["terms"] = term_list,
    _["pairs"] = pair_out
  );
}

// Per-peptide prevalence of one contrast for the plotters: per group level the
// mean over subjects of their share of positive samples (and, with
// adj_weight, its inverse-propensity weighted version, Hajek-normalised
// within level), and per timepoint the share of positive samples.
//
// [[Rcpp::export]]
Rcpp::List cpp_delta2_prevalence(const Rcpp::RawVector& bitset_raw,
                                 int n_words,
                                 const Rcpp::IntegerVector& pep_cols,
                                 const Rcpp::IntegerVector& sample_rows,
                                 const Rcpp::IntegerVector& sample_subject,
                                 const Rcpp::IntegerVector& group_code,
                                 const Rcpp::IntegerVector& time_code,
                                 int n_group_levels, int n_time_levels,
                                 const Rcpp::NumericVector& adj_weight) {
  const int K = std::max(n_group_levels, 1);
  for (R_xlen_t i = 0; i < group_code.size(); ++i) {
    const int g = group_code[i];
    if (g != NA_INTEGER && (g < 1 || g > K)) stop("group_code out of range.");
  }
  const bool has_adj = adj_weight.size() > 0;
  if (has_adj && adj_weight.size() != group_code.size()) stop("adj_weight must have one value per subject.");

  Contrast cx;
  cx.init(n_words, sample_rows, sample_subject, time_code, &group_code, n_time_levels);
  const int T = cx.T;
  const int N = cx.n_sub;
  const int m = pep_cols.size();

  std::vector<int> lab(N);
  Vec a_i(N, 1.0), n_k(K, 0.0), asum_k(K, 0.0), n_t(T, 0.0);
  for (int i = 0; i < N; ++i) {
    lab[i] = group_code[cx.sub_orig[i] - 1] - 1;
    if (has_adj) a_i[i] = adj_weight[cx.sub_orig[i] - 1];
    n_k[lab[i]] += 1.0;
    asum_k[lab[i]] += a_i[i];
    for (int t = 0; t < T; ++t) n_t[t] += cx.n_it[(size_t)i * T + t];
  }

  NumericMatrix g_prev(m, K), g_prev_adj(m, K), t_prev(m, T);
  NumericVector tot(m);
  Vec St(T);
  for (int j = 0; j < m; ++j) {
    const uint8_t* col = column_ptr(bitset_raw, n_words, pep_cols[j]);
    std::fill(St.begin(), St.end(), 0.0);
    tot[j] = cx.count(col, n_words, St);
    for (int i : cx.touched) {
      double c = 0.0;
      for (int t = 0; t < T; ++t) c += cx.cnt[(size_t)i * T + t];
      const double ybar = c / cx.n_i[i];
      g_prev(j, lab[i]) += ybar;
      g_prev_adj(j, lab[i]) += a_i[i] * ybar;
    }
    for (int k = 0; k < K; ++k) {
      g_prev(j, k) = n_k[k] > 0 ? g_prev(j, k) / n_k[k] : NA_REAL;
      g_prev_adj(j, k) = asum_k[k] > 0 ? g_prev_adj(j, k) / asum_k[k] : NA_REAL;
    }
    for (int t = 0; t < T; ++t) t_prev(j, t) = n_t[t] > 0 ? St[t] / n_t[t] : NA_REAL;
    cx.reset();
  }

  return List::create(
    _["group_prev"] = g_prev,
    _["group_prev_adj"] = g_prev_adj,
    _["time_prev"] = t_prev,
    _["n_group"] = NumericVector(n_k.begin(), n_k.end()),
    _["n_time"] = NumericVector(n_t.begin(), n_t.end()),
    _["tot"] = tot
  );
}
