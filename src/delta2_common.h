// Shared helpers of the compute_delta2() C++ code (delta2_score.cpp,
// delta2_contrast.cpp).
#ifndef PHIPER_DELTA2_COMMON_H
#define PHIPER_DELTA2_COMMON_H

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <algorithm>

namespace delta2 {

using namespace Rcpp;

typedef std::vector<double> Vec;

inline uint64_t load_word(const uint8_t* p) {
  uint64_t w;
  std::memcpy(&w, p, sizeof(w));
  return w;
}

// Eigen decomposition of a small symmetric n x n matrix (row-major), cyclic
// Jacobi. vecs holds the eigenvectors as columns.
inline void sym_eigen(Vec a, int n, Vec& vals, Vec& vecs) {
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
inline Vec psd_pinv_pow(const Vec& m, int n, double power, double scale = 0.0) {
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
inline Vec level_design(int n_levels, const NumericVector& scores, int& L) {
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

inline const uint8_t* column_ptr(const RawVector& bitset_raw, int n_words, int pep_col) {
  const size_t col_bytes = (size_t)n_words * 8;
  const size_t off = (size_t)(pep_col - 1) * col_bytes;
  if (pep_col < 1 || off + col_bytes > (size_t)bitset_raw.size()) {
    stop("pep_cols out of range of the bitset.");
  }
  return reinterpret_cast<const uint8_t*>(RAW(bitset_raw)) + off;
}

} // namespace delta2

#endif
