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

// compute_delta2() permutation kernel for one (rank, feature, group_col)
// stratum. The per-peptide statistics are those of cpp_delta2_score(); this
// pools them over peptides (strata of pooled prevalence, Stouffer / maxmean /
// AF per component, T = sum of squared components for df > 1) and repeats
// that for every draw. The loop runs over draws outside and peptides inside:
// each draw builds its shared object once (a subject labelling, or the
// per-timepoint sample masks) and nothing per peptide.

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

// Pools one component's peptide z (already winsorised) into a scalar: the
// statistic per non-empty prevalence bin, averaged over those bins.
struct Pooler {
  AggStat stat;
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
      if (stat == AggStat::stouffer) {
        num[b] += w[j] * z[j];
        den2[b] += w[j] * w[j];
      } else if (stat == AggStat::maxmean) {
        if (z[j] > 0.0) { sp[b] += z[j]; ++cp[b]; }
        else if (z[j] < 0.0) { sn[b] += -z[j]; ++cn[b]; }
      } else {
        zs[b].push_back(z[j]);
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

} // namespace

// Pooled statistics, permutation null and effect sizes of one stratum.
//
// Subjects are indexed by sample_subject into group_code; group_draws holds
// one resampled labelling per column (rows as group_code) and time_draws one
// within-subject relabelling of the time codes per column (rows as
// sample_rows). Draws must keep the observed label multiset. pairs lists the
// level pairs (1-based, group1 < group2) for the pairwise tests; adj_weight
// is each subject's inverse-propensity weight for the adjusted effect sizes
// (empty = unweighted). rho is one value or one per peptide.
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
                               double min_m_eff) {
  bool do_group = false, do_time = false, do_inter = false;
  for (R_xlen_t k = 0; k < terms.size(); ++k) {
    const std::string tm = Rcpp::as<std::string>(terms[k]);
    if (tm == "group") do_group = true;
    else if (tm == "time") do_time = true;
    else if (tm == "interaction") do_inter = true;
    else stop("Unknown term: %s", tm);
  }
  const int P = pairs.nrow();
  const bool need_group_draws = do_group || do_inter || P > 0;
  const int K = std::max(n_group_levels, 1);
  const int m = pep_cols.size();
  if (m == 0) stop("No peptides.");
  if (rho.size() != 1 && rho.size() != m) stop("rho must have length 1 or one per peptide.");
  for (R_xlen_t i = 0; i < rho.size(); ++i) {
    if (!(rho[i] >= 0.0 && rho[i] <= 1.0)) stop("rho must be in [0, 1].");
  }
  for (R_xlen_t i = 0; i < group_code.size(); ++i) {
    const int g = group_code[i];
    if (g != NA_INTEGER && (g < 1 || g > K)) stop("group_code out of range.");
  }

  int B = 0;
  if (need_group_draws) {
    if (group_draws.nrow() != group_code.size()) stop("group_draws must have one row per subject.");
    B = group_draws.ncol();
  }
  if (do_time) {
    if (time_draws.nrow() != sample_rows.size()) stop("time_draws must have one row per sample.");
    if (B > 0 && time_draws.ncol() != B) stop("group_draws and time_draws must have the same B.");
    B = time_draws.ncol();
  }
  if (B < 1) stop("Need at least one draw.");

  Contrast cx;
  cx.init(n_words, sample_rows, sample_subject, time_code, &group_code, n_time_levels);
  const int T = cx.T;
  const int N = cx.n_sub;

  int Lg = 1, Lt = 1;
  const Vec Gd = level_design(K, group_scores, Lg);
  const Vec Xd = level_design(T, time_scores, Lt);
  const int LL = Lt * Lt;

  // observed labels (0-based) per compact subject; draws read through sub_orig
  std::vector<int> obs_lab(N);
  std::vector<double> n_k(K, 0.0);
  for (int i = 0; i < N; ++i) {
    obs_lab[i] = group_code[cx.sub_orig[i] - 1] - 1;
    n_k[obs_lab[i]] += 1.0;
  }
  const bool has_adj = adj_weight.size() > 0;
  if (has_adj && adj_weight.size() != group_code.size()) stop("adj_weight must have one value per subject.");
  Vec a_i(N, 1.0);
  if (has_adj) for (int i = 0; i < N; ++i) a_i[i] = adj_weight[cx.sub_orig[i] - 1];

  // time-label masks: kept samples (bitset row, original time, position)
  std::vector<int> kept_pos, kept_row;
  for (int s = 0; s < sample_rows.size(); ++s) {
    const int r = sample_rows[s];
    if (r == NA_INTEGER || r < 1 || r > n_words * 64) continue;
    if (cx.slot_sub[r - 1] < 0) continue;
    kept_pos.push_back(s);
    kept_row.push_back(r - 1);
  }

  // --- per-subject design spread A_i -----------------------------------------
  Vec A((size_t)N * LL, 0.0);
  std::vector<char> spread(N, 0);
  for (int i = 0; i < N; ++i) {
    Vec xs(Lt, 0.0);
    for (int t = 0; t < T; ++t)
      for (int l = 0; l < Lt; ++l) xs[l] += Xd[(size_t)t * Lt + l] * cx.n_it[(size_t)i * T + t];
    double tr = 0.0;
    for (int a = 0; a < Lt; ++a) {
      for (int b = 0; b < Lt; ++b) {
        double s = 0.0;
        for (int t = 0; t < T; ++t)
          s += Xd[(size_t)t * Lt + a] * Xd[(size_t)t * Lt + b] * cx.n_it[(size_t)i * T + t];
        A[(size_t)i * LL + a * Lt + b] = s - xs[a] * xs[b] / cx.n_i[i];
      }
      tr += A[(size_t)i * LL + a * Lt + a];
    }
    spread[i] = (tr > 1e-12) ? 1 : 0;
  }
  int N_I = 0;
  Vec SumA(LL, 0.0), Q4((size_t)LL * LL, 0.0);
  for (int i = 0; i < N; ++i) {
    if (!spread[i]) continue;
    ++N_I;
    const double* Ai = &A[(size_t)i * LL];
    for (int e = 0; e < LL; ++e) SumA[e] += Ai[e];
    for (int e = 0; e < LL; ++e)
      for (int f = 0; f < LL; ++f) Q4[(size_t)e * LL + f] += Ai[e] * Ai[f];
  }
  const Vec SumA_pinv = psd_pinv_pow(SumA, Lt, -1.0);

  // between-subject whitening: fixed, the label multiset is fixed
  Vec gbar_b(Lg, 0.0), Wg_b((size_t)Lg * Lg, 0.0);
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

  // interaction shared object of a labelling: centred group design over the
  // spread subjects, its whitening, and sum_i (g_i - gbar) (x) A_i
  struct InterShared { Vec gbar, Wg, GA; };
  auto inter_shared = [&](const int* lab, InterShared& sh) {
    sh.gbar.assign(Lg, 0.0);
    sh.Wg.assign((size_t)Lg * Lg, 0.0);
    sh.GA.assign((size_t)Lg * LL, 0.0);
    if (N_I < 2) return;
    for (int i = 0; i < N; ++i) {
      if (!spread[i]) continue;
      for (int l = 0; l < Lg; ++l) sh.gbar[l] += Gd[(size_t)lab[i] * Lg + l];
    }
    for (int l = 0; l < Lg; ++l) sh.gbar[l] /= N_I;
    Vec GI((size_t)Lg * Lg, 0.0);
    for (int i = 0; i < N; ++i) {
      if (!spread[i]) continue;
      const double* Ai = &A[(size_t)i * LL];
      for (int a = 0; a < Lg; ++a) {
        const double ga = Gd[(size_t)lab[i] * Lg + a] - sh.gbar[a];
        for (int b = 0; b < Lg; ++b)
          GI[a * Lg + b] += ga * (Gd[(size_t)lab[i] * Lg + b] - sh.gbar[b]) / (N_I - 1.0);
        for (int e = 0; e < LL; ++e) sh.GA[(size_t)a * LL + e] += ga * Ai[e];
      }
    }
    sh.Wg = psd_pinv_pow(GI, Lg, -0.5);
  };

  // --- observed pass: per-peptide fixed quantities -----------------------------
  // CSR over each peptide's positive subjects
  std::vector<int> off(m + 1, 0), csr_sub;
  Vec csr_v, csr_D;
  std::vector<char> csr_spread;
  Vec tot(m), p0(m), sumv(m), SS(m), Wt((size_t)m * LL, 0.0), cst((size_t)m * Lt, 0.0),
      beta((size_t)m * Lt, 0.0), WH((size_t)m * LL, 0.0);
  std::vector<char> g_ok(m, 0);
  std::vector<Vec> St_obs(m, Vec(T, 0.0));

  // effect sizes per peptide
  const int n_pairs = P;
  Vec d_group(m, 0.0), d_group_adj(m, 0.0), d_time(m, 0.0), d_inter(m, 0.0), d_inter_adj(m, 0.0);
  Vec d_pair((size_t)m * n_pairs, 0.0), d_pair_adj((size_t)m * n_pairs, 0.0);
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

  Vec D(Lt), sumD(Lt), HDD(LL), H(LL), Ab(Lt), Vt(LL), Rbar(Lt), mk(K), mk_adj(K),
      Dk((size_t)K * Lt), Dk_adj((size_t)K * Lt), bk((size_t)K * Lt), bk_adj((size_t)K * Lt);
  for (int j = 0; j < m; ++j) {
    const uint8_t* col = column_ptr(bitset_raw, n_words, pep_cols[j]);
    const double rj = (rho.size() == 1) ? rho[0] : rho[j];
    tot[j] = cx.count(col, n_words, St_obs[j]);
    p0[j] = cx.n_samples > 0 ? tot[j] / cx.n_samples : 0.0;

    double sv = 0.0, sv2 = 0.0;
    std::fill(sumD.begin(), sumD.end(), 0.0);
    std::fill(HDD.begin(), HDD.end(), 0.0);
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
      const double v = ni / (1.0 + (ni - 1.0) * rj) * ybar;
      sv += v;
      sv2 += v * v;
      mk[obs_lab[i]] += ybar;
      mk_adj[obs_lab[i]] += a_i[i] * ybar;
      csr_sub.push_back(i);
      csr_v.push_back(v);
      csr_spread.push_back(spread[i]);
      for (int l = 0; l < Lt; ++l) {
        double s = 0.0;
        for (int t = 0; t < T; ++t)
          s += Xd[(size_t)t * Lt + l] * (ci[t] - cx.n_it[(size_t)i * T + t] * ybar);
        D[l] = spread[i] ? s : 0.0;
        csr_D.push_back(D[l]);
        // time const: sum_t X[t,l] n_it ybar_i, over every positive subject
        double q = 0.0;
        for (int t = 0; t < T; ++t) q += Xd[(size_t)t * Lt + l] * cx.n_it[(size_t)i * T + t] * ybar;
        cst[(size_t)j * Lt + l] += q;
      }
      if (!spread[i]) continue;
      const double Sy = c * (ni - c) / ni;
      const double f = Sy / (ni - 1.0);
      for (int e = 0; e < LL; ++e) Vt[e] += f * A[(size_t)i * LL + e];
      for (int l = 0; l < Lt; ++l) {
        sumD[l] += D[l];
        Dk[(size_t)obs_lab[i] * Lt + l] += D[l];
        Dk_adj[(size_t)obs_lab[i] * Lt + l] += a_i[i] * D[l];
      }
      for (int a = 0; a < Lt; ++a)
        for (int b = 0; b < Lt; ++b) HDD[a * Lt + b] += D[a] * D[b];
    }
    off[j + 1] = csr_sub.size();
    sumv[j] = sv;
    SS[j] = (N >= 2) ? sv2 - sv * sv / N : 0.0;
    g_ok[j] = SS[j] > 1e-12 * sv2;

    // time whitening
    if (do_time) {
      const Vec W = psd_pinv_pow(Vt, Lt, -0.5);
      std::copy(W.begin(), W.end(), Wt.begin() + (size_t)j * LL);
    }
    // time main effect (within estimator) and interaction residual spread
    for (int a = 0; a < Lt; ++a) {
      double s = 0.0;
      for (int b = 0; b < Lt; ++b) s += SumA_pinv[a * Lt + b] * sumD[b];
      beta[(size_t)j * Lt + a] = s;
    }
    if (do_time) {
      const double* bj = &beta[(size_t)j * Lt];
      if (Lt == 1) {
        d_time[j] = bj[0] * tscale;
      } else {
        double lo = bj[0], hi = bj[0];
        for (int a = 1; a < Lt; ++a) { lo = std::min(lo, bj[a]); hi = std::max(hi, bj[a]); }
        d_time[j] = hi - lo;
      }
    }
    if (do_inter && N_I >= 2) {
      const double* bj = &beta[(size_t)j * Lt];
      for (int e = 0; e < LL; ++e) H[e] = HDD[e];
      double h_scale = 0.0;
      for (int l = 0; l < Lt; ++l) h_scale += HDD[l * Lt + l];
      for (int e = off[j]; e < off[j + 1]; ++e) {
        if (!csr_spread[e]) continue;
        const double* De = &csr_D[(size_t)e * Lt];
        const double* Ai = &A[(size_t)csr_sub[e] * LL];
        for (int a = 0; a < Lt; ++a) {
          double s = 0.0;
          for (int b = 0; b < Lt; ++b) s += Ai[a * Lt + b] * bj[b];
          Ab[a] = s;
        }
        for (int a = 0; a < Lt; ++a)
          for (int b = 0; b < Lt; ++b) H[a * Lt + b] -= De[a] * Ab[b] + Ab[a] * De[b];
      }
      for (int l = 0; l < Lt; ++l)
        for (int q = 0; q < Lt; ++q) {
          double s = 0.0;
          for (int p = 0; p < Lt; ++p)
            for (int mm = 0; mm < Lt; ++mm)
              s += bj[p] * bj[mm] * Q4[(size_t)(l * Lt + p) * LL + (mm * Lt + q)];
          H[l * Lt + q] += s;
          if (l == q) h_scale += s;
        }
      for (int l = 0; l < Lt; ++l) {
        double s = sumD[l];
        for (int p = 0; p < Lt; ++p) s -= SumA[l * Lt + p] * bj[p];
        Rbar[l] = s / N_I;
      }
      for (int a = 0; a < Lt; ++a)
        for (int b = 0; b < Lt; ++b) H[a * Lt + b] -= N_I * Rbar[a] * Rbar[b];
      const Vec W = psd_pinv_pow(H, Lt, -0.5, h_scale);
      std::copy(W.begin(), W.end(), WH.begin() + (size_t)j * LL);

      // interaction effect: within estimate per group, contrasted across groups
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
    for (int p = 0; p < n_pairs; ++p) {
      const int a = pairs(p, 0) - 1, b = pairs(p, 1) - 1;
      d_pair[(size_t)p * m + j] = mk[b] - mk[a];
      d_pair_adj[(size_t)p * m + j] = mk_adj[b] - mk_adj[a];
    }
    cx.reset();
  }

  // --- peptide weights and m_eff ---------------------------------------------
  Pooler pool;
  pool.stat = parse_agg(aggregate_stat);
  pool.bin = prevalence_bins(strat_bins, p0, pool.n_bins);
  pool.w.assign(m, 1.0);
  for (int j = 0; j < m; ++j) {
    if (weight_mode == "n_eff_sqrt") {
      pool.w[j] = std::sqrt(std::max(tot[j], 1e-12));
    } else if (weight_mode == "se_invvar") {
      const double q = p0[j] * (1.0 - p0[j]);
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
  // r buffers: component-major, m entries per component
  Vec r_group((size_t)Lg * m), r_time((size_t)Lt * m), r_inter((size_t)Lg * Lt * m),
      z_pair((size_t)n_pairs * m), S1(K), S2(K), U(std::max(Lg, Lt)), SD((size_t)Lg * Lt),
      Sint((size_t)Lg * Lt), cnt_u(T);
  auto clamp = [&](double z) { return std::max(-winsor_z, std::min(winsor_z, z)); };

  // between-subject terms and pairs under a labelling
  auto between_r = [&](const int* lab, const InterShared* sh) {
    for (int j = 0; j < m; ++j) {
      std::fill(S1.begin(), S1.end(), 0.0);
      std::fill(S2.begin(), S2.end(), 0.0);
      for (int e = off[j]; e < off[j + 1]; ++e) {
        const int k = lab[csr_sub[e]];
        S1[k] += csr_v[e];
        S2[k] += csr_v[e] * csr_v[e];
      }
      if (do_group) {
        for (int l = 0; l < Lg; ++l) {
          double s = -gbar_b[l] * sumv[j];
          for (int k = 0; k < K; ++k) s += Gd[(size_t)k * Lg + l] * S1[k];
          U[l] = s;
        }
        for (int a = 0; a < Lg; ++a) {
          double r = 0.0;
          if (g_ok[j]) {
            for (int b = 0; b < Lg; ++b) r += Wg_b[a * Lg + b] * U[b];
            r /= std::sqrt(SS[j]);
          }
          r_group[(size_t)a * m + j] = clamp(r);
        }
      }
      for (int p = 0; p < n_pairs; ++p) {
        const int a = pairs(p, 0) - 1, b = pairs(p, 1) - 1;
        const double na = n_k[a], nb = n_k[b], Nab = na + nb;
        double z = 0.0;
        if (na > 0 && nb > 0 && Nab >= 2) {
          const double s1 = S1[a] + S1[b], s2 = S2[a] + S2[b];
          const double ss = s2 - s1 * s1 / Nab;
          if (ss > 1e-12 * s2) {
            z = (S1[b] - nb * s1 / Nab) / std::sqrt(na * nb / (Nab * (Nab - 1.0)) * ss);
          }
        }
        z_pair[(size_t)p * m + j] = clamp(z);
      }
      if (do_inter) {
        if (N_I < 2) {
          for (int e = 0; e < Lg * Lt; ++e) r_inter[(size_t)e * m + j] = 0.0;
          continue;
        }
        std::fill(SD.begin(), SD.end(), 0.0);
        for (int e = off[j]; e < off[j + 1]; ++e) {
          if (!csr_spread[e]) continue;
          const int k = lab[csr_sub[e]];
          const double* De = &csr_D[(size_t)e * Lt];
          for (int l = 0; l < Lg; ++l) {
            const double gl = Gd[(size_t)k * Lg + l] - sh->gbar[l];
            for (int a = 0; a < Lt; ++a) SD[(size_t)l * Lt + a] += gl * De[a];
          }
        }
        const double* bj = &beta[(size_t)j * Lt];
        for (int l = 0; l < Lg; ++l)
          for (int a = 0; a < Lt; ++a) {
            double s = 0.0;
            for (int p = 0; p < Lt; ++p) s += sh->GA[(size_t)l * LL + a * Lt + p] * bj[p];
            Sint[(size_t)l * Lt + a] = SD[(size_t)l * Lt + a] - s;
          }
        const double* Wh = &WH[(size_t)j * LL];
        for (int l = 0; l < Lg; ++l)
          for (int a = 0; a < Lt; ++a) {
            double r = 0.0;
            for (int b = 0; b < Lg; ++b)
              for (int q = 0; q < Lt; ++q)
                r += sh->Wg[l * Lg + b] * Wh[a * Lt + q] * Sint[(size_t)b * Lt + q];
            r_inter[(size_t)(l * Lt + a) * m + j] = clamp(r);
          }
      }
    }
  };

  // time term from relabelled per-timepoint counts
  auto time_r_from = [&](int j, const Vec& counts) {
    for (int l = 0; l < Lt; ++l) {
      double s = -cst[(size_t)j * Lt + l];
      for (int u = 0; u < T; ++u) s += Xd[(size_t)u * Lt + l] * counts[u];
      U[l] = s;
    }
    const double* W = &Wt[(size_t)j * LL];
    for (int a = 0; a < Lt; ++a) {
      double r = 0.0;
      for (int b = 0; b < Lt; ++b) r += W[a * Lt + b] * U[b];
      r_time[(size_t)a * m + j] = clamp(r);
    }
  };

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

  NullAcc acc_group, acc_time, acc_inter;
  std::vector<NullAcc> acc_pair(n_pairs);
  Vec E_group, E_time, E_inter;

  // observed
  {
    InterShared sh;
    if (do_inter) inter_shared(obs_lab.data(), sh);
    if (need_group_draws) between_r(obs_lab.data(), &sh);
    if (do_group) acc_group.T_obs = pooled(r_group, Lg, &E_group);
    if (do_inter) acc_inter.T_obs = pooled(r_inter, Lg * Lt, &E_inter);
    for (int p = 0; p < n_pairs; ++p) {
      acc_pair[p].T_obs = pool(&z_pair[(size_t)p * m]);
    }
    if (do_time) {
      for (int j = 0; j < m; ++j) time_r_from(j, St_obs[j]);
      acc_time.T_obs = pooled(r_time, Lt, &E_time);
    }
  }

  // draws
  Vec max_pair(B, 0.0);
  std::vector<int> lab(N);
  std::vector<uint64_t> G((size_t)T * n_words);
  const uint8_t* bits = reinterpret_cast<const uint8_t*>(RAW(bitset_raw));
  for (int b = 0; b < B; ++b) {
    if (need_group_draws) {
      for (int i = 0; i < N; ++i) {
        const int l = group_draws(cx.sub_orig[i] - 1, b);
        if (l == NA_INTEGER || l < 1 || l > K) stop("group_draws out of range.");
        lab[i] = l - 1;
      }
      InterShared sh;
      if (do_inter) inter_shared(lab.data(), sh);
      between_r(lab.data(), &sh);
      if (do_group) acc_group.digest(pooled(r_group, Lg, nullptr));
      if (do_inter) acc_inter.digest(pooled(r_inter, Lg * Lt, nullptr));
      double mx = 0.0;
      for (int p = 0; p < n_pairs; ++p) {
        const double tp = pool(&z_pair[(size_t)p * m]);
        acc_pair[p].digest(tp);
        mx = std::max(mx, std::fabs(tp));
      }
      max_pair[b] = mx;
    }
    if (do_time) {
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
        time_r_from(j, cnt_u);
      }
      acc_time.digest(pooled(r_time, Lt, nullptr));
    }
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
  if (do_group) term_list["group"] = term_out(acc_group, E_group, d_group, d_group_adj, Lg == 1);
  if (do_time) term_list["time"] = term_out(acc_time, E_time, d_time, d_time, Lt == 1);
  if (do_inter) {
    term_list["interaction"] = term_out(acc_inter, E_inter, d_inter, d_inter_adj, Lg == 1 && Lt == 1);
  }

  List pair_out;
  if (n_pairs > 0) {
    NumericVector T_obs(n_pairs), p_perm(n_pairs), p_maxT(n_pairs), mx(n_pairs),
        mx_adj(n_pairs), fp(n_pairs), fpw(n_pairs);
    IntegerVector bb(n_pairs);
    for (int p = 0; p < n_pairs; ++p) {
      const NullAcc& a = acc_pair[p];
      T_obs[p] = a.T_obs;
      bb[p] = a.b_hits;
      p_perm[p] = perm_p(a.b_hits, a.b_equal, B, mid_p);
      int gt = 0, eq = 0;
      for (int b = 0; b < B; ++b) {
        const int c = compare_abs(max_pair[b], a.T_obs);
        if (c > 0) ++gt;
        else if (c == 0) ++eq;
      }
      p_maxT[p] = perm_p(gt, eq, B, mid_p);
      const Vec d(d_pair.begin() + (size_t)p * m, d_pair.begin() + (size_t)(p + 1) * m);
      const Vec da(d_pair_adj.begin() + (size_t)p * m, d_pair_adj.begin() + (size_t)(p + 1) * m);
      const List ds = delta_summary(d, da, w_norm, true);
      mx[p] = ds["max_delta"];
      mx_adj[p] = ds["max_delta_adj"];
      fp[p] = ds["frac_delta_pos"];
      fpw[p] = ds["frac_delta_pos_w"];
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
    _["n_samples"] = cx.n_samples,
    _["n_subjects_interaction"] = N_I,
    _["B"] = B,
    _["terms"] = term_list,
    _["pairs"] = pair_out
  );
}
