# Changelog

## phiper 0.3.1

### Changes to `compute_alpha`

- Added `pielou_evenness` and `berger_parker_dominance` to the output
  (NA for samples with richness ≤ 1 and richness = 0 respectively).
- New `metrics` parameter: request any subset of the five indices;
  defaults to all five.
- New `mode` parameter (`"binary"`, `"threshold"`, `"abundance"`)
  replaces `fc_threshold`. `"abundance"` mode uses raw values from
  `abundance_col` with optional `abundance_agg` (`"mean"`, `"sum"`,
  `"max"`) at higher ranks.
- `shannon_log` renamed to `shannon_base`; old name still works with a
  deprecation warning.
- Performance: all-samples roster now collected once before the rank
  loop instead of re-queried per rank.
- Hardening: all-invalid ranks now aborts with an informative error
  instead of silently returning empty output; `n_samples` attribute
  added to the result.
- Validation: `mode = "threshold"` now requires `threshold` to be finite
  and `abundance_col` (when supplied) to be a character scalar.

## phiper 0.3.0

### Major changes

- Extracted all data-import, class construction, and low-level utility
  code into the new **phiperio** package. `phiper` now declares
  `phiperio` as a hard dependency and re-exports its user-facing
  functions (`load_example_data`, `get_example_path`) so existing
  workflows require no changes.
- Removed all functions that moved to phiperio: `phip_convert`,
  `phip_convert_legacy`, `new_phip_data`, `expand_phip_data`,
  `phip_data_join`, `validate_phip_data`, `disconnect`,
  `get_comparisons`, `phip_example_path`, `phip_load_example_data`, and
  the full logging / validation helper suite (`add_quotes`, `word_list`,
  `.chk_*`).
- Internal logging and validation now use the unified phiperio helpers
  (`.ph_abort`, `.ph_warn`, `.ph_log_info`, `.ph_with_timing`,
  `.ph_check_cond`, `.ph_add_quotes`, `.ph_word_list`, `.ph_check_path`,
  `.ph_check_extension`, `.ph_check_null_default`).
- `get_peptide_meta()` renamed to `get_peptide_library()` throughout, in
  line with the phiperio API.
- `compute_alpha`: restored efficient same-connection peptide library
  handling via `.ph_peplib_on_main()` (DuckDB ATTACH fast path with
  [`copy_to()`](https://dplyr.tidyverse.org/reference/copy_to.html)
  fallback).

### Internal

- R source files renamed to follow the new `<domain>_compute` /
  `<domain>_plots` convention: `binary_alpha` → `alpha_compute`,
  `binary_alpha_plots` → `alpha_plots`, `binary_beta` → `beta_compute`,
  `binary_beta_plots` → `beta_plots`, `prevalence-DELTA_test` →
  `delta_compute`, `prevalence-DELTA_plots` → `delta_plots`,
  `prevalence-POP-test` → `pop_compute`, `prevalence-POP_plots` →
  `pop_plots`, `plot-utils` → `plot_utils`. Test files renamed
  accordingly.
- Naming conventions for source files and functions documented in
  `CONTRIBUTING.md`: exported functions use plain `snake_case`; all
  internal helpers use the `.ph_` prefix.

## phiper 0.2.7

- compute_delta: added `maxmean` (Efron-type) as a test statistic
  option.
- compute_delta: added prevalence bins for the test statistic
  computation.
- compute_delta: `srlr` is now the default test statistic.
- Implemented McNemar test and paired signed root likelihood ratio
  statistic for paired designs.
- Welford’s online algorithm refactored into its own class; added a
  wrapper for post-permutation output generation.
- Docs and R CMD CHECK fixes.

## phiper 0.2.6

- compute_delta: added stat_mode options “score” (pooled score z) and
  “srlr” (signed root likelihood ratio) using raw counts.
- compute_delta: removed smoothing, prevalence filtering, BH adjustment,
  and fold_change/cross_prev summaries; outputs now include T_null_mean
  and T_null_sd and standardized T_obs uses the null mean and sample SD.
- shift_computing: null variance uses Welford’s algorithm with sample
  variance in permutation loops.
- prevalence-DELTA plots/tests updated to drop BH/categorical
  dependencies and use numeric-only filtering where applicable.

## phiper 0.2.5

- Removed 10 unused package dependencies (data.table, fs, htmltools,
  purrr, forcats, filelock, arrow from suggests)
- Moved duckdb and dbplyr to imports; moved knitr to suggests for proper
  dependency management
- Removed 360 lines of unused CAP/dispersion plotting functions and
  associated documentation

## phiper 0.2.4

- Updated the peptide library with the new annotations from Sasha
- in the compute_distance –\> fallback to collecting before pivoting

## phiper 0.2.3

- Added `small_mixture` to `phip_load_example_data`
- Changed all examples and tests to use this `small_mixture` instead of
  redining it
- Added cache to `phip_load_example_data` to make tests/examples faster

## phiper 0.2.2

- Separated the feature associations to PCoA vectors from `compute_pcoa`
  to a separate function.

## phiper 0.2.1

### Minor changes

- Removed superassignments, changed all assignments to “\<-” style
- Stated all dependencies, removed unstated dependencies

### Bug fixes

- Fixed typos in the examples from binary_beta.R
- fixed R CMD check (0 errors, 0 warnings, 0 notes)

### Documentation

- Documented empty params, got rid of all documentation-related
  wearnings,
- removed all non-ASCII chars
- Added CONTRIBUTING.md to .Rbuildignore (it caused a note in R CMD
  CHECK),

## phiper 0.2.0

### Minor changes

- Removed generic and S3 methods for the expand_phip_data
- Renamed the internal helpers involved in the data import to match the
  naming convention .ph\_
- documented the internal helpers involved in the data import stage

### Major changes

- Removed the backend argument entirely from the phip_convert,
  phip_convert_legacy, new_phip_data, .resolve_paths and other, less
  imporatant helpers. DuckDB is now the only supported backend
- Changed the structure of the repo. Now all functions related to the
  standard/legacy import workflows live in the standard-conver.R or
  legacy-convert.R
- removed the resolve-paths.R file. Moved the functions to utils.R

## phiper 0.1.1

### Minor changes

- Removed old, dead and undocumented code: R/binary-analysis_peptides.R,
  R/vinary-analysis_stability.R, R/binary-analysis_stability_cattime.R,
  R/fold_change-analysis.R

- i left other files, even if they were undocumented/untested, as they
  were essential for phiper to work (e.g. the data import paths)

## phiper 0.1.1

### Minor changes

- Removed old, dead and undocumented code: R/binary-analysis_peptides.R,
  R/vinary-analysis_stability.R, R/binary-analysis_stability_cattime.R,
  R/fold_change-analysis.R

- i left other files, even if they were undocumented/untested, as they
  were essential for phiper to work (e.g. the data import paths)

### Documentation

- Regenerated the docs, removed the old functions from NAMESPACE

## phiper 0.1.0

### Minor changes

- None.

### Major changes

- None.

### Features

- None.

### Bug fixes

- None.

### Documentation

- None.

### Internal

- None.
