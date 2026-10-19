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
//   group        between-subject: subject values v_i = w_i * ybar_i, group
//                labels resampled across subjects with sizes fixed.
//   time         within-subject: the subject's deviations from its own mean,
//                timepoints reshuffled inside each subject.
//   interaction  the within-subject deviations, residualised against the
//                pooled time main effect, compared across group labels.
// A df = 1 term (two levels or a score trend) gives a signed z; an omnibus
// term gives the vector r = Cov^{+1/2} U, with sum(r^2) the quadratic form.

namespace {

typedef std::vector<double> Vec;

inline uint64_t load_word(const uint8_t* p) {
  uint64_t w;
  std::memcpy(&w, p, sizeof(w));
  return w;
}

// Eigen decomposition of a small symmetric n x n matrix (row-major), cyclic
// Jacobi. vecs holds the eigenvectors as columns.
void sym_eigen(Vec a, int n, Vec& vals, Vec& vecs) {
  vecs.assign((size_t)n * n, 0.0);
  for (int i = 0; i < n; ++i) vecs[i * n + i] = 1.0;
  for (int sweep = 0; sweep < 100; ++sweep) {
    double off = 0.0, diag = 0.0;
    for (int i = 0; i < n; ++i) {
      diag += a[i * n + i] * a[i * n + i];
      for (int j = i + 1; j < n; ++j) off += a[i * n + j] * a[i * n + j];
    }
    if (off <= 1e-30 * (diag + 1e-300)) break;
    for (int p = 0; p < n; ++p) {
      for (int q = p + 1; q < n; ++q) {
        const double apq = a[p * n + q];
        if (apq == 0.0) continue;
        const double theta = (a[q * n + q] - a[p * n + p]) / (2.0 * apq);
        const double t = (theta >= 0.0 ? 1.0 : -1.0) /
          (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
        const double c = 1.0 / std::sqrt(t * t + 1.0);
        const double s = t * c;
        for (int k = 0; k < n; ++k) {
          const double akp = a[k * n + p], akq = a[k * n + q];
          a[k * n + p] = c * akp - s * akq;
          a[k * n + q] = s * akp + c * akq;
        }
        for (int k = 0; k < n; ++k) {
          const double apk = a[p * n + k], aqk = a[q * n + k];
          a[p * n + k] = c * apk - s * aqk;
          a[q * n + k] = s * apk + c * aqk;
        }
        for (int k = 0; k < n; ++k) {
          const double vkp = vecs[k * n + p], vkq = vecs[k * n + q];
          vecs[k * n + p] = c * vkp - s * vkq;
          vecs[k * n + q] = s * vkp + c * vkq;
        }
      }
    }
  }
  vals.resize(n);
  for (int i = 0; i < n; ++i) vals[i] = a[i * n + i];
}

// Symmetric pseudo-inverse of a PSD matrix raised to `power` (-1 or -1/2).
// Eigenvalues below 1e-10 * max(largest eigenvalue, scale) are treated as
// zero; `scale` is the size of the terms m was summed from, so a matrix that
// is rounding noise after cancellation counts as zero.
Vec psd_pinv_pow(const Vec& m, int n, double power, double scale = 0.0) {
  Vec out((size_t)n * n, 0.0);
  if (n == 1) {
    if (m[0] > 1e-10 * scale && m[0] > 0.0) out[0] = std::pow(m[0], power);
    return out;
  }
  Vec vals, vecs;
  sym_eigen(m, n, vals, vecs);
  double vmax = 0.0;
  for (double v : vals) vmax = std::max(vmax, v);
  if (vmax <= 1e-10 * scale || vmax <= 0.0) return out;
  const double tol = 1e-10 * std::max(vmax, scale);
  for (int e = 0; e < n; ++e) {
    if (vals[e] <= tol) continue;
    const double f = std::pow(vals[e], power);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        out[i * n + j] += f * vecs[i * n + e] * vecs[j * n + e];
  }
  return out;
}

// Level design: two levels or scores give one column, otherwise one
// indicator column per level. Row-major, n_levels x L.
Vec level_design(int n_levels, const NumericVector& scores, int& L) {
  if (scores.size() > 0) {
    if (scores.size() != n_levels) stop("scores must have one value per level.");
    L = 1;
    return Vec(scores.begin(), scores.end());
  }
  if (n_levels == 2) {
    L = 1;
    return Vec{0.0, 1.0};
  }
  L = n_levels;
  Vec d((size_t)n_levels * n_levels, 0.0);
  for (int k = 0; k < n_levels; ++k) d[k * n_levels + k] = 1.0;
  return d;
}

// Shared walk over a contrast: samples kept, their subject and timepoint,
// per-subject sample counts, and per-peptide positive counts per
// (subject, timepoint) collected sparsely through a touched-subject list.
struct Contrast {
  int T = 1;
  int n_sub = 0;                 // subjects with >= 1 sample in the contrast
  std::vector<int> slot_sub;     // bitset row (0-based) -> compact subject, -1
  std::vector<int> slot_time;    // bitset row (0-based) -> timepoint 0..T-1
  std::vector<int> sub_orig;     // compact subject -> 1-based input index
  std::vector<double> n_i;       // samples per subject
  std::vector<double> n_it;      // samples per (subject, timepoint)
  double n_samples = 0.0;

  // per-peptide scratch
  std::vector<double> cnt;       // positives per (subject, timepoint)
  std::vector<int> touched;
  std::vector<char> is_touched;

  void init(int n_words, const IntegerVector& sample_rows,
            const IntegerVector& sample_subject, const IntegerVector& time_code,
            const IntegerVector* group_code, int n_time_levels) {
    T = std::max(n_time_levels, 1);
    const int ns = sample_rows.size();
    if (sample_subject.size() != ns) stop("sample_subject must match sample_rows.");
    if (time_code.size() != ns) stop("time_code must match sample_rows.");
    const int n_rows = n_words * 64;
    slot_sub.assign(n_rows, -1);
    slot_time.assign(n_rows, 0);

    int max_sub = 0;
    for (int s = 0; s < ns; ++s) {
      if (sample_subject[s] != NA_INTEGER) max_sub = std::max(max_sub, (int)sample_subject[s]);
    }
    std::vector<int> compact(max_sub + 1, -1);
    for (int s = 0; s < ns; ++s) {
      const int r = sample_rows[s], sub = sample_subject[s], t = time_code[s];
      if (r == NA_INTEGER || sub == NA_INTEGER || t == NA_INTEGER) continue;
      if (group_code != nullptr) {
        if (sub > group_code->size()) stop("sample_subject out of range of group_code.");
        if ((*group_code)[sub - 1] == NA_INTEGER) continue;
      }
      if (r < 1 || r > n_rows) stop("sample_rows out of range.");
      if (t < 1 || t > T) stop("time_code out of range.");
      if (slot_sub[r - 1] >= 0) stop("Duplicated sample_rows.");
      if (compact[sub] < 0) {
        compact[sub] = n_sub++;
        sub_orig.push_back(sub);
      }
      slot_sub[r - 1] = compact[sub];
      slot_time[r - 1] = t - 1;
    }
    n_i.assign(n_sub, 0.0);
    n_it.assign((size_t)n_sub * T, 0.0);
    for (int r = 0; r < n_rows; ++r) {
      const int i = slot_sub[r];
      if (i < 0) continue;
      n_i[i] += 1.0;
      n_it[(size_t)i * T + slot_time[r]] += 1.0;
      n_samples += 1.0;
    }
    cnt.assign((size_t)n_sub * T, 0.0);
    is_touched.assign(n_sub, 0);
  }

  // Count positives of one bitset column; returns the number of positives.
  double count(const uint8_t* col, int n_words, std::vector<double>& St) {
    double tot = 0.0;
    for (int w = 0; w < n_words; ++w) {
      uint64_t bits = load_word(col + (size_t)w * 8);
      while (bits) {
        const int b = __builtin_ctzll(bits);
        bits &= bits - 1;
        const int r = w * 64 + b;
        const int i = slot_sub[r];
        if (i < 0) continue;
        const int t = slot_time[r];
        cnt[(size_t)i * T + t] += 1.0;
        St[t] += 1.0;
        tot += 1.0;
        if (!is_touched[i]) {
          is_touched[i] = 1;
          touched.push_back(i);
        }
      }
    }
    return tot;
  }

  void reset() {
    for (int i : touched) {
      is_touched[i] = 0;
      std::fill(cnt.begin() + (size_t)i * T, cnt.begin() + (size_t)(i + 1) * T, 0.0);
    }
    touched.clear();
  }
};

const uint8_t* column_ptr(const RawVector& bitset_raw, int n_words, int pep_col) {
  const size_t col_bytes = (size_t)n_words * 8;
  const size_t off = (size_t)(pep_col - 1) * col_bytes;
  if (pep_col < 1 || off + col_bytes > (size_t)bitset_raw.size()) {
    stop("pep_cols out of range of the bitset.");
  }
  return reinterpret_cast<const uint8_t*>(RAW(bitset_raw)) + off;
}

} // namespace

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

    double sumv = 0.0, sumv2 = 0.0;
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

    // group: Cov(U) = SS_v * Gb, so r = Gb^{-1/2} U / sqrt(SS_v)
    // Gb carries the exact permutation 1/(N-1); compute_delta()'s "score"
    // uses the binomial 1/N, so with one sample per subject this z is v1's
    // times sqrt((N-1)/N), the same factor for every peptide in a contrast.
    if (do_group) {
      const double SS = (N >= 2) ? sumv2 - sumv * sumv / N : 0.0;
      const bool ok = SS > 1e-12 * sumv2;
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
// Also returns each subject's total positives (breadth) over pep_cols.
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
  for (int j = 0; j < pep_cols.size(); ++j) {
    const uint8_t* col = column_ptr(bitset_raw, n_words, pep_cols[j]);
    St[0] = 0.0;
    const double C = cx.count(col, n_words, St);
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
    _["subject_hits"] = hits
  );
}
