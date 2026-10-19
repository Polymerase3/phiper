// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <algorithm>
using namespace Rcpp;

// compute_delta2() score engine: per-peptide score z for the group, time and
// interaction terms of a sample-level bitset (see R/delta2_design.R).
//
// Every term is a linear statistic of the sample bits standardised by its
// exact permutation moments:
//   group        between-subject: subject values u_i = w_i (ybar_i - mu_w),
//                mu_w the w-weighted mean of the subject means, group
//                labels resampled across subjects with sizes fixed.
//   time         within-subject: the subject's deviations from its own mean,
//                timepoints reshuffled inside each subject.
//   interaction  the within-subject deviations, residualised against the
//                pooled time main effect, compared across group labels.
// A df = 1 term (two levels or a score trend) gives a signed z; an omnibus
// term gives the vector r = Cov^{+1/2} U, with sum(r^2) the quadratic form.

#include "delta2_common.h"
using namespace delta2;

// Score z per peptide per term.
//
// sample_rows, sample_subject and time_code are per sample (same length);
// group_code is per subject (indexed by sample_subject). Samples with an NA
// code, or whose subject has an NA group_code, are left out. Empty
// group_scores / time_scores mean omnibus over the levels (a single 0/1
// contrast when there are two). rho sets the between-subject weights
// w_i = n_i / (1 + (n_i - 1) rho).
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
                            double rho) {
  if (n_words <= 0) stop("n_words must be positive.");
  if (!(rho >= 0.0 && rho <= 1.0)) stop("rho must be in [0, 1].");
  bool do_group = false, do_time = false, do_inter = false;
  for (R_xlen_t k = 0; k < terms.size(); ++k) {
    const std::string tm = Rcpp::as<std::string>(terms[k]);
    if (tm == "group") do_group = true;
    else if (tm == "time") do_time = true;
    else if (tm == "interaction") do_inter = true;
    else stop("Unknown term: %s", tm);
  }
  const int K = std::max(n_group_levels, 1);
  for (R_xlen_t i = 0; i < group_code.size(); ++i) {
    const int g = group_code[i];
    if (g != NA_INTEGER && (g < 1 || g > K)) stop("group_code out of range.");
  }

  Contrast cx;
  cx.init(n_words, sample_rows, sample_subject, time_code, &group_code, n_time_levels);
  const int T = cx.T;
  const int N = cx.n_sub;
  const int m = pep_cols.size();

  int Lg = 1, Lt = 1;
  const Vec Gd = level_design(K, group_scores, Lg);
  const Vec Xd = level_design(T, time_scores, Lt);

  // --- per-subject fixed quantities ----------------------------------------
  Vec g((size_t)N * Lg), w(N);
  for (int i = 0; i < N; ++i) {
    const int k = group_code[cx.sub_orig[i] - 1] - 1;
    for (int l = 0; l < Lg; ++l) g[(size_t)i * Lg + l] = Gd[(size_t)k * Lg + l];
    w[i] = cx.n_i[i] / (1.0 + (cx.n_i[i] - 1.0) * rho);
  }

  // A_i = X_i' C X_i, the subject's within-subject design spread (Lt x Lt)
  const int LL = Lt * Lt;
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

  // between-subject design: centred group scores over all subjects
  Vec gbar_b(Lg, 0.0), Wg_b((size_t)Lg * Lg, 0.0);
  if (do_group && N >= 2) {
    for (int i = 0; i < N; ++i)
      for (int l = 0; l < Lg; ++l) gbar_b[l] += g[(size_t)i * Lg + l] / N;
    Vec Gb((size_t)Lg * Lg, 0.0);
    for (int i = 0; i < N; ++i)
      for (int a = 0; a < Lg; ++a)
        for (int b = 0; b < Lg; ++b)
          Gb[a * Lg + b] += (g[(size_t)i * Lg + a] - gbar_b[a]) *
            (g[(size_t)i * Lg + b] - gbar_b[b]) / (N - 1.0);
    Wg_b = psd_pinv_pow(Gb, Lg, -0.5);
  }
  // centring at the weighted mean mu_w = sum w ybar / sum w: these sums over
  // all subjects turn the sparse sums over positive subjects into
  // sum_i (g_i - gbar) u_i and sum_i u_i^2
  double Wsum = 0.0, W2sum = 0.0;
  Vec GW(Lg, 0.0);
  for (int i = 0; i < N; ++i) {
    Wsum += w[i];
    W2sum += w[i] * w[i];
    for (int l = 0; l < Lg; ++l) GW[l] += (g[(size_t)i * Lg + l] - gbar_b[l]) * w[i];
  }

  // interaction: subjects with within-subject spread, residualised against
  // the pooled time main effect beta = (sum A_i)^+ sum D_i
  int N_I = 0;
  Vec gbar_I(Lg, 0.0), Wg_I((size_t)Lg * Lg, 0.0), SumA(LL, 0.0), SumA_pinv(LL, 0.0);
  Vec GA((size_t)Lg * LL, 0.0);         // sum_i (g_i - gbar) (x) A_i
  Vec Q4((size_t)LL * LL, 0.0);         // [l,p,m,q] = sum_i A_i[l,p] A_i[m,q]
  if (do_inter) {
    for (int i = 0; i < N; ++i) {
      if (!spread[i]) continue;
      ++N_I;
      for (int l = 0; l < Lg; ++l) gbar_I[l] += g[(size_t)i * Lg + l];
    }
    if (N_I >= 2) {
      for (int l = 0; l < Lg; ++l) gbar_I[l] /= N_I;
      Vec GI((size_t)Lg * Lg, 0.0);
      for (int i = 0; i < N; ++i) {
        if (!spread[i]) continue;
        const double* Ai = &A[(size_t)i * LL];
        for (int a = 0; a < Lg; ++a) {
          const double ga = g[(size_t)i * Lg + a] - gbar_I[a];
          for (int b = 0; b < Lg; ++b)
            GI[a * Lg + b] += ga * (g[(size_t)i * Lg + b] - gbar_I[b]) / (N_I - 1.0);
          for (int e = 0; e < LL; ++e) GA[(size_t)a * LL + e] += ga * Ai[e];
        }
        for (int e = 0; e < LL; ++e) SumA[e] += Ai[e];
        for (int e = 0; e < LL; ++e)
          for (int f = 0; f < LL; ++f) Q4[(size_t)e * LL + f] += Ai[e] * Ai[f];
      }
      Wg_I = psd_pinv_pow(GI, Lg, -0.5);
      SumA_pinv = psd_pinv_pow(SumA, Lt, -1.0);
    }
  }

  // --- outputs ---------------------------------------------------------------
  NumericMatrix z_group(do_group ? m : 0, Lg);
  NumericMatrix z_time(do_time ? m : 0, Lt);
  NumericMatrix z_inter(do_inter ? m : 0, Lg * Lt);
  NumericVector tot(m);
  NumericMatrix St_out(m, T);

  Vec St(T), Ub(Lg), Ut(Lt), Vt(LL), D(Lt), sumD(Lt), SD((size_t)Lg * Lt),
      HDD(LL), beta(Lt), Ab(Lt), Sint((size_t)Lg * Lt), H(LL), Rbar(Lt);

  for (int j = 0; j < m; ++j) {
    const uint8_t* col = column_ptr(bitset_raw, n_words, pep_cols[j]);
    std::fill(St.begin(), St.end(), 0.0);
    tot[j] = cx.count(col, n_words, St);
    for (int t = 0; t < T; ++t) St_out(j, t) = St[t];

    double sumv = 0.0, sumv2 = 0.0, sumwv = 0.0;
    std::fill(Ub.begin(), Ub.end(), 0.0);
    std::fill(Ut.begin(), Ut.end(), 0.0);
    std::fill(Vt.begin(), Vt.end(), 0.0);
    std::fill(sumD.begin(), sumD.end(), 0.0);
    std::fill(SD.begin(), SD.end(), 0.0);
    std::fill(HDD.begin(), HDD.end(), 0.0);

    // pass 1: everything that only needs the touched subjects
    for (int i : cx.touched) {
      const double* ci = &cx.cnt[(size_t)i * T];
      double c = 0.0;
      for (int t = 0; t < T; ++t) c += ci[t];
      const double ni = cx.n_i[i];
      const double ybar = c / ni;

      if (do_group) {
        const double v = w[i] * ybar;
        sumv += v;
        sumv2 += v * v;
        sumwv += w[i] * v;
        for (int l = 0; l < Lg; ++l) Ub[l] += (g[(size_t)i * Lg + l] - gbar_b[l]) * v;
      }
      if (!(do_time || do_inter) || !spread[i]) continue;

      for (int l = 0; l < Lt; ++l) {
        double s = 0.0;
        for (int t = 0; t < T; ++t)
          s += Xd[(size_t)t * Lt + l] * (ci[t] - cx.n_it[(size_t)i * T + t] * ybar);
        D[l] = s;
      }
      if (do_time) {
        // exact within-subject variance: S_x,i * S_y,i / (n_i - 1)
        const double Sy = c * (ni - c) / ni;
        const double f = Sy / (ni - 1.0);
        for (int l = 0; l < Lt; ++l) Ut[l] += D[l];
        for (int e = 0; e < LL; ++e) Vt[e] += f * A[(size_t)i * LL + e];
      }
      if (do_inter) {
        for (int l = 0; l < Lt; ++l) sumD[l] += D[l];
        for (int a = 0; a < Lg; ++a) {
          const double ga = g[(size_t)i * Lg + a] - gbar_I[a];
          for (int l = 0; l < Lt; ++l) SD[(size_t)a * Lt + l] += ga * D[l];
        }
        for (int a = 0; a < Lt; ++a)
          for (int b = 0; b < Lt; ++b) HDD[a * Lt + b] += D[a] * D[b];
      }
    }

    // group: U = sum_i (g_i - gbar) u_i with u_i = w_i (ybar_i - mu_w);
    // Cov(U) = SS_u * Gb, so r = Gb^{-1/2} U / sqrt(SS_u). Centring at mu_w
    // keeps U free of sum_{i in group} w_i, which moves with the labels and
    // would otherwise shift every peptide together when the w_i differ.
    // Gb carries the exact permutation 1/(N-1); compute_delta()'s "score"
    // uses the binomial 1/N, so with one sample per subject this z is v1's
    // times sqrt((N-1)/N), the same factor for every peptide in a contrast.
    if (do_group) {
      const double mu = (Wsum > 0.0) ? sumv / Wsum : 0.0;
      const double SS = (N >= 2) ? sumv2 - 2.0 * mu * sumwv + mu * mu * W2sum : 0.0;
      const bool ok = SS > 1e-12 * (sumv2 + mu * mu * W2sum);
      for (int b = 0; b < Lg; ++b) Ub[b] -= mu * GW[b];
      for (int a = 0; a < Lg; ++a) {
        double r = 0.0;
        if (ok) {
          for (int b = 0; b < Lg; ++b) r += Wg_b[a * Lg + b] * Ub[b];
          r /= std::sqrt(SS);
        }
        z_group(j, a) = r;
      }
    }

    // time: r = Vt^{+1/2} Ut
    if (do_time) {
      const Vec W = psd_pinv_pow(Vt, Lt, -0.5);
      for (int a = 0; a < Lt; ++a) {
        double r = 0.0;
        for (int b = 0; b < Lt; ++b) r += W[a * Lt + b] * Ut[b];
        z_time(j, a) = r;
      }
    }

    // interaction: residuals R_i = D_i - A_i beta, Cov = G_I (x) H
    if (do_inter) {
      if (N_I >= 2) {
        for (int a = 0; a < Lt; ++a) {
          double s = 0.0;
          for (int b = 0; b < Lt; ++b) s += SumA_pinv[a * Lt + b] * sumD[b];
          beta[a] = s;
        }
        // sum_i (g_i - gbar) (x) R_i
        for (int a = 0; a < Lg; ++a)
          for (int l = 0; l < Lt; ++l) {
            double s = 0.0;
            for (int p = 0; p < Lt; ++p) s += GA[(size_t)a * LL + l * Lt + p] * beta[p];
            Sint[(size_t)a * Lt + l] = SD[(size_t)a * Lt + l] - s;
          }
        // H = sum R R' = sum DD' - sum (D (A b)' + (A b) D') + sum A b b' A
        for (int e = 0; e < LL; ++e) H[e] = HDD[e];
        double h_scale = 0.0;
        for (int l = 0; l < Lt; ++l) h_scale += HDD[l * Lt + l];
        for (int i : cx.touched) {
          if (!spread[i]) continue;
          const double* ci = &cx.cnt[(size_t)i * T];
          double c = 0.0;
          for (int t = 0; t < T; ++t) c += ci[t];
          const double ybar = c / cx.n_i[i];
          for (int l = 0; l < Lt; ++l) {
            double s = 0.0;
            for (int t = 0; t < T; ++t)
              s += Xd[(size_t)t * Lt + l] * (ci[t] - cx.n_it[(size_t)i * T + t] * ybar);
            D[l] = s;
          }
          const double* Ai = &A[(size_t)i * LL];
          for (int a = 0; a < Lt; ++a) {
            double s = 0.0;
            for (int b = 0; b < Lt; ++b) s += Ai[a * Lt + b] * beta[b];
            Ab[a] = s;
          }
          for (int a = 0; a < Lt; ++a)
            for (int b = 0; b < Lt; ++b) H[a * Lt + b] -= D[a] * Ab[b] + Ab[a] * D[b];
        }
        for (int l = 0; l < Lt; ++l)
          for (int q = 0; q < Lt; ++q) {
            double s = 0.0;
            for (int p = 0; p < Lt; ++p)
              for (int mm = 0; mm < Lt; ++mm)
                s += beta[p] * beta[mm] * Q4[(size_t)(l * Lt + p) * LL + (mm * Lt + q)];
            H[l * Lt + q] += s;
            if (l == q) h_scale += s;
          }
        // centre: Rbar = (sum D - SumA beta) / N_I, zero up to rounding
        for (int l = 0; l < Lt; ++l) {
          double s = sumD[l];
          for (int p = 0; p < Lt; ++p) s -= SumA[l * Lt + p] * beta[p];
          Rbar[l] = s / N_I;
        }
        for (int a = 0; a < Lt; ++a)
          for (int b = 0; b < Lt; ++b) H[a * Lt + b] -= N_I * Rbar[a] * Rbar[b];
        const Vec WH = psd_pinv_pow(H, Lt, -0.5, h_scale);
        for (int a = 0; a < Lg; ++a)
          for (int l = 0; l < Lt; ++l) {
            double r = 0.0;
            for (int b = 0; b < Lg; ++b)
              for (int q = 0; q < Lt; ++q)
                r += Wg_I[a * Lg + b] * WH[l * Lt + q] * Sint[(size_t)b * Lt + q];
            z_inter(j, a * Lt + l) = r;
          }
      } else {
        for (int e = 0; e < Lg * Lt; ++e) z_inter(j, e) = 0.0;
      }
    }

    cx.reset();
  }

  List out = List::create(
    _["tot"] = tot,
    _["S_time"] = St_out,
    _["n_subjects"] = N,
    _["n_samples"] = cx.n_samples,
    _["rho"] = rho
  );
  if (do_group) out["group"] = z_group;
  if (do_time) out["time"] = z_time;
  if (do_inter) {
    out["interaction"] = z_inter;
    out["n_subjects_interaction"] = N_I;
  }
  return out;
}

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
