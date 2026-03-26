# Global Shift in Peptide-level Prevalence via Subject-level Permutation

Test for a **global** (antigen-wide) shift in peptide-level prevalence
between two groups within each `(rank, feature, group_col)` stratum by
aggregating per-peptide effects into a single Stouffer-type statistic
and assessing significance via **label permutation**.

## Usage

``` r
compute_delta(
  x,
  rank_cols,
  group_cols,
  exist_col = "exist",
  paired_by = NULL,
  interaction = FALSE,
  combine_cols = NULL,
  interaction_sep = "::",
  B_permutations = 2000L,
  weight_mode = c("equal", "se_invvar", "n_eff_sqrt"),
  stat_mode = c("srlr", "diff", "asin", "score", "mcnemar", "srlr_paired"),
  aggregate_stat = c("stouffer", "maxmean"),
  strat_bins = c(0.002, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50),
  winsor_z = 4,
  rank_feature_keep = NULL,
  peptide_library = NULL,
  log = FALSE,
  log_file = "compute_delta.log"
)
```

## Arguments

- x:

  A `phip_data` object (recommended to attach a `peptide_library` when
  `rank_cols` includes non-`"peptide_id"` and not providing
  `peptide_library` argument) or a data frame with the necessary
  columns.

- rank_cols:

  Character vector of rank columns (e.g. `"peptide_id"`, `"species"`,
  ...). For non-peptide ranks, annotations are joined from a peptide
  library (see Details).

- group_cols:

  Character vector of grouping columns that define universes; pairwise
  group contrasts are constructed within each
  `(rank, feature, group_col)` stratum.

- exist_col:

  Name of the 0/1 presence column. Default `"exist"`.

- paired_by:

  Optional single column name identifying paired samples (e.g. subject
  ID). When provided, paired tests are used where applicable.

- interaction:

  Logical; reserved for future use (currently ignored). If `TRUE`, an
  interaction of the first two `group_cols` may be used as an additional
  `group_col` in future versions.

- combine_cols:

  Optional length-2 character vector; reserved for future use. If
  non-`NULL`, only this interaction would be used as `group_col` in
  future versions.

- interaction_sep:

  Separator for interaction values. Default `"::"`.

- B_permutations:

  Number of permutations `B` used for the permutation p-value. Default
  `2000L`. Must be at least 100.

- weight_mode:

  One of `c("equal", "se_invvar", "n_eff_sqrt")`, controlling how
  peptide-level z-scores are weighted in the Stouffer combination (see
  Details).

- stat_mode:

  One of `c("diff", "asin", "score", "srlr", "mcnemar", "srlr_paired")`,
  determining the peptide-level test statistic (see Details). The
  paired-only modes require paired designs.

- aggregate_stat:

  One of `c("stouffer", "maxmean")`, controlling how peptide-level
  z-scores are aggregated into a single test statistic (see Details).

- strat_bins:

  Numeric vector of pooled prevalence cutpoints in (0, 1). If `0`, no
  stratification is used (single global Stouffer statistic). Otherwise,
  peptides are binned by pooled prevalence and the mean of bin-level
  z-scores is used as the global test statistic.

- winsor_z:

  Winsorization threshold applied to peptide-level z-scores. Values
  beyond `±winsor_z` are truncated. Default `4.0`.

- rank_feature_keep:

  Optional **named list** mapping `rank` to a vector of `feature` values
  to keep. Only rank–feature strata in this list are tested; others are
  dropped after the peptide-level pivot.

- peptide_library:

  Optional data frame providing peptide annotations for non-peptide
  ranks. Must at least contain `peptide_id` and all requested
  `rank_cols` besides `"peptide_id"`. If `NULL`, the function falls back
  to `x$peptide_library` or
  [`phiper::get_peptide_meta()`](https://polymerase3.github.io/phiper/reference/get_peptide_meta.md)
  if available.

- log:

  Logical; if `TRUE`, write progress messages (per contrast and overall)
  using the package's logging helpers.

- log_file:

  Path to a log file used by the logging helpers if `log` is `TRUE`.
  Default `"compute_delta.log"`.

## Value

A tibble with one row per tested stratum:

- `rank`, `feature`: rank and feature identifiers.

- `group_col`, `group1`, `group2`: grouping column and the ordered pair
  of groups compared.

- `design`: `"paired"` or `"unpaired"`.

- `n_subjects_paired`: number of paired subjects used in a paired design
  (or `NA` for unpaired).

- `n_peptides_used`: number of peptides contributing to the test.

- `m_eff`: effective number of peptides after weighting (as returned by
  the C++ helper).

- `T_obs`: observed Stouffer-type test statistic.

- `p_perm`: two-sided permutation p-value.

- `b`: number of permutations actually used (may be `< B_permutations`
  if early stopping is implemented in the C++ helper).

- `max_delta`, `frac_delta_pos`, `frac_delta_pos_w`: maximum absolute
  peptide-level prevalence difference and unweighted/weighted fractions
  of positive peptide-level deltas.

## Details

For each stratum `(rank, feature, group_col)` that defines an ordered
pair of groups `(g1, g2)`, the procedure is:

1.  **Peptide-level prevalences.** For each peptide, prevalences in `g1`
    and `g2` are \$\$ \hat p_k = \frac{x_k}{n_k},\quad k \in \\g1, g2\\,
    \$\$ where `x_k` is the number of samples with presence
    (`exist > 0`) and `n_k` is the number of distinct samples (or paired
    subjects).

2.  **Per-peptide effect and z-score (`stat_mode`).**

    - `"diff"`: effect \\\Delta = \hat p\_{g2} - \hat p\_{g1}\\ with
      \$\$ SE = \sqrt{\hat p\_{g1}(1 - \hat p\_{g1}) / n\_{g1} + \hat
      p\_{g2}(1 - \hat p\_{g2}) / n\_{g2}} , \$\$ and z-score \\z =
      \Delta / SE\\ (with a small floor on `SE`).

    - `"asin"`: variance-stabilized difference of angles, \$\$ z =
      \frac{\arcsin \sqrt{\hat p\_{g2}} - \arcsin \sqrt{\hat p\_{g1}}}
      {\sqrt{1 / (4 n\_{g1}) + 1 / (4 n\_{g2})}} . \$\$

    - `"score"`: pooled score z for equality of proportions (2×2 score
      test), using raw proportions \\\tilde p_k = x_k / n_k\\ and pooled
      \\\tilde p = (x_1 + x_2)/(n_1 + n_2)\\: \$\$ z = \frac{\tilde
      p\_{g2} - \tilde p\_{g1}} {\sqrt{\tilde p(1-\tilde
      p)\left(1/n\_{g1} + 1/n\_{g2}\right)}} . \$\$

    - `"srlr"`: signed root likelihood ratio for \\H_0:
      p\_{g1}=p\_{g2}\\: \$\$ LR = 2\left\[\ell(\hat
      p\_{g1};x_1,n_1)+\ell(\hat p\_{g2};x_2,n_2) -\ell(\hat
      p;x_1,n_1)-\ell(\hat p;x_2,n_2)\right\], \$\$ \$\$ z =
      \mathrm{sign}(\hat p\_{g2}-\hat p\_{g1})\sqrt{LR}, \$\$ where
      \\\ell(p;x,n)=x\log p+(n-x)\log(1-p)\\ and \\\hat
      p\_{gk}=x_k/n_k\\, \\\hat p=(x_1+x_2)/(n_1+n_2)\\.

    - `"mcnemar"` (paired-only): signed McNemar z based on discordant
      pairs \\b = \\(g1=1,g2=0)\\, \\c = \\(g1=0,g2=1)\\: \$\$ z =
      \frac{b - c}{\sqrt{b + c}} . \$\$ If \\b+c=0\\, z is defined as 0.

    - `"srlr_paired"` (paired-only): signed root deviance for the paired
      binomial problem with \\m=b+c\\: \$\$ z =
      \mathrm{sign}(b-c)\sqrt{2\left\[b\log\frac{2b}{m} +
      c\log\frac{2c}{m}\right\]} . \$\$ Terms with \\0\log 0\\ are
      defined as 0; if \\m=0\\, z is 0.

    Each \\z_i\\ is then **winsorized** to \\\pm\\ `winsor_z`.

3.  **Weights (`weight_mode`).**

    - `"equal"`: all peptides get weight 1.

    - `"se_invvar"`: weights proportional to the inverse standard error
      (or the analogous denominator for `"asin"`).

    - `"n_eff_sqrt"`: \$\$ w_i = \sqrt{ n\_{g1} \hat p\_{g1} + n\_{g2}
      \hat p\_{g2} }, \$\$ i.e. the square root of the expected number
      of positives across both groups.

4.  **Combine into a single test statistic (`aggregate_stat`).** Let
    \\w_i\\ be the weights and \\z_i\\ the peptide-level z-scores.

    - If `aggregate_stat = "stouffer"` and `strat_bins = 0`: \$\$
      T\_{\mathrm{obs}} = \frac{\sum_i w_i z_i}{\sqrt{\sum_i w_i^2}} .
      \$\$

    - If `aggregate_stat = "stouffer"` and `strat_bins` is a numeric
      vector of pooled prevalence cutpoints (between 0 and 1), peptides
      are binned by pooled prevalence \$\$ \frac{n\_{g1}\hat p\_{g1} +
      n\_{g2}\hat p\_{g2}}{n\_{g1} + n\_{g2}}, \$\$ a Stouffer statistic
      \\T_b\\ is computed within each bin. The **mean** of the bin-level
      z-scores is reported as \\T\_{\mathrm{obs}}\\. For example, if
      `strat_bins` includes `0.10` and `0.20`, then peptides with pooled
      prevalence in (0.10, 0.20\] share a bin.

    - If `aggregate_stat = "maxmean"`, the signed z-scores are split
      into positive and negative parts, averaged within each set, and
      the dominant direction is selected. If `strat_bins` is provided,
      the maxmean statistic is computed per bin and the **mean** across
      bins is reported.

5.  **Permutation scheme and p-value.** A **two-sided** permutation
    p-value for the global shift is computed:

    - **Paired design**: if both groups have measurements for the same
      `subject_id`, each subject’s labels (`g1` \\\leftrightarrow\\
      `g2`) are independently flipped with probability 1/2 and steps
      (1–4) are recomputed to obtain \\T_b\\.

    - **Unpaired design**: sample labels are shuffled while preserving
      group sizes (random split into `n1` / `n2`), and steps (1–4) are
      recomputed to obtain \\T_b\\.

    With \\B\\ permutations, the p-value is \$\$ p = \frac{1 +
    \sum\_{b=1}^{B} \mathbf{1}\\\|T_b\| \ge \|T\_{\mathrm{obs}}\|\\}
    {1 + B}, \$\$ the standard add-one estimator that remains non-zero
    under the global null.

6.  **Multiplicity.** No multiple-testing correction is performed;
    `p_perm` is returned as computed.

**Input requirements.**

- `exist_col` is treated as 0/1 presence.

- There must be **at most one positive** per (`subject_id`,
  `peptide_id`, `group_col`, `group_value`); paired designs can have up
  to two positives across the two group levels. Violations trigger an
  error. Example (group levels A/B): for a single subject and peptide,
  you may have A=1 and B=0 (or A=0 and B=1, or A=1 and B=1), but you
  cannot have two rows both with A=1 (or two rows both with B=1).

- Non-peptide ranks specified in `rank_cols` must be resolvable from a
  peptide library (see `peptide_library` below).

**Peptide library resolution.**

When `rank_cols` includes non-`"peptide_id"` ranks, a peptide library is
required. It is resolved in the following order:

1.  The explicit `peptide_library` argument.

2.  `x$peptide_library` if `x` is a `phip_data` with an attached
    library.

3.  [`phiper::get_peptide_meta()`](https://polymerase3.github.io/phiper/reference/get_peptide_meta.md)
    if available (i.e. `phiper` is installed and exports this function).

If none of these are available, the function aborts with an informative
error.

## Parallelization

The permutation contrasts are evaluated either sequentially or in
parallel using the **current** `future` backend:

- The function does **not** change the global
  [`future::plan()`](https://future.futureverse.org/reference/plan.html).

- If `future` is installed, the number of workers is inferred from
  [`future::nbrOfWorkers()`](https://future.futureverse.org/reference/nbrOfWorkers.html).
  To run in parallel, set up a `future` plan (e.g.
  `future::plan(multisession, workers = 8)`) **outside** the function.

- Internal C++ code is constrained to a single thread per R worker
  (`RcppParallel::setThreadOptions(numThreads = 1)` and BLAS/OpenMP
  limits), to avoid oversubscription when using multiple workers.

Reproducibility of permutations is controlled via the global RNG: call
[`set.seed()`](https://rdrr.io/r/base/Random.html) before
`compute_delta()` to obtain reproducible results; each contrast draws
its own seed from the global RNG state.

## Examples

``` r
# \donttest{
# Load example PhIP-Seq data shipped with the package
pd <- phip_load_example_data()

# Small unpaired subset with a mock peptide library
pd_filt <- pd |>
  dplyr::filter(
    peptide_id %in% c("16627", "5243", "24799", "16196", "18003"),
    timepoint == "T1"
  ) |>
  dplyr::collect()

mock_peplib <- data.frame(
  peptide_id = c("16627", "5243", "24799", "16196", "18003"),
  species    = rep("mock_species", 5),
  stringsAsFactors = FALSE
)

res <- compute_delta(
  x                  = pd_filt,
  exist_col          = "exist",
  rank_cols          = "species",
  group_cols         = "group",
  peptide_library    = mock_peplib,
  B_permutations     = 500L,  # smaller for speed
  weight_mode        = "n_eff_sqrt",
  stat_mode          = "asin",
  strat_bins         = 0,
  winsor_z           = Inf,
  rank_feature_keep  = list(species = NULL),
  log                = FALSE
)
res
#> # A tibble: 1 × 19
#>   rank  feature group_col group1 group2 design n_subjects_paired n_peptides_used
#>   <chr> <chr>   <chr>     <chr>  <chr>  <chr>              <int>           <int>
#> 1 spec… mock_s… group     A      B      unpai…                NA               5
#> # ℹ 11 more variables: m_eff <dbl>, T_obs <dbl>, T_null_mean <dbl>,
#> #   T_null_sd <dbl>, T_obs_stand <dbl>, Z_from_p <dbl>, p_perm <dbl>, b <int>,
#> #   max_delta <dbl>, frac_delta_pos <dbl>, frac_delta_pos_w <dbl>
# }
```
