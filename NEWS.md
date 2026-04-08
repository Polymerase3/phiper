# phiper 0.3.3

## Tests

- New `test-beta_plots.R`: tests for `plot_pcoa()`, `plot_cap()`, `plot_scree()`,
  `plot_dispersion()`, and `plot_tsne()` (2-D and 3-D), covering basic output
  type, grouping/time aesthetics, centroid and ellipse options, axis selection,
  variance-explained labels, and input-validation errors.
- New `test-plot_utils.R`: tests for colour helpers (`phip_palette`,
  `scale_colour_phip()`, `scale_fill_phip()`, `theme_phip()`,
  `phip_use_montserrat()`), internal colour utilities (`.hex2rgb()`,
  `.rgb2hex()`, `.mix_cols()`, `.tint()`, `.blend_hex()`, `.make_shades()`,
  `.build_shaded_map()`), and ordination helpers (`.pick_axes()`,
  `.axis_labels_with_pct()`, `.shaded_colors()`, `.make_point_fills()`,
  `.first_subview_name()`).
- New `test-shift_computing.R`: tests for `compute_delta()` covering all
  `stat_mode` options (`diff`, `score`, `srlr`, `mcnemar`, `srlr_paired`), all
  `weight_mode` options, stratified bins (`strat_bins`), winsorisation, and
  paired designs.
- New `test-zzz.R`: tests for `.onLoad()` idempotency, `load_example_data()`,
  and `get_example_path()`.
- New `test-utils.R`: tests for internal utilities including `.ph_check_cond()`,
  `.ph_check_extension()`, `.ph_check_null_default()`, `.ph_check_path()`,
  `%nin%`, `%||%`, `.ph_opt()`, `.ph_now()`, `.ph_base_prefix()`,
  `.ph_wrap()`, `.ph_compose_lines()`, `.ph_log_info()`, `.ph_log_ok()`,
  `.ph_warn()`, `.ph_abort()`, `.ph_with_timing()`, `.ph_check_pd()`, and
  `.ph_resolve_paths()`.

# phiper 0.3.2

## New functions

- `compute_alpha_significance()`: runs global (Kruskal-Wallis or one-way ANOVA)
  and pairwise (Wilcoxon or Tukey HSD) hypothesis tests on every `(rank, metric)`
  combination returned by `compute_alpha()`. Returns a `"phip_alpha_significance"`
  list with `$global` (statistic, p-value, test) and `$pairwise` (p_raw, p_adj,
  Cohen's d, significance stars) tibbles. Supports `p_adjust_method` (default
  `"BH"`) and subsetting via `metric` / `rank` arguments. Group column is inferred
  automatically when not supplied.
- `plot_alpha_significance()`: visualises pairwise results either as a filtered
  tibble (`type = "table"`) or a symmetric Cohen's d heatmap with significance
  stars (`type = "heatmap"`). Accepts `metric`, `rank`, and `p_threshold`
  arguments to focus on a single comparison.

## Changes to `plot_alpha()` and `plot_alpha_interactive()`

- `metric` now also accepts `"pielou_evenness"` and `"berger_parker_dominance"`,
  matching all five indices added to `compute_alpha()` in 0.3.1.
- New significance-bracket parameters: `significance`, `show_significance`,
  `sig_p_threshold`, `sig_step_increase`, `sig_tip_length`. Pass a
  `"phip_alpha_significance"` object and set `show_significance = TRUE` to overlay
  pairwise brackets via `ggsignif` (optional dependency). Non-significant pairs
  are automatically omitted.
- Added `...` (reserved; ignored) for forward-compatibility.

## Font bundling

- Montserrat Regular, Bold, and Italic TTF files are now shipped in `inst/fonts/`,
  so `phip_use_montserrat()` works offline without a Google Fonts connection.
- `phip_use_montserrat()` tries the local bundle first; falls back to
  `sysfonts::font_add_google()` only when the local files are absent.
- Package registers the font at load time via `.onLoad()` (showtext rendering
  remains opt-in; call `phip_use_montserrat()` explicitly to enable it).

## Tests

- New `test-alpha_significance.R`: 56 tests covering `compute_alpha_significance()`
  (global/pairwise tests, p-adjustment, Cohen's d, metric/rank subsetting, group
  inference, single-group edge case) and `plot_alpha_significance()` (table and
  heatmap modes, error paths). Line coverage: 89.91% → 99.54%.
- New `test-alpha_plots.R`: 63 tests covering all five metrics (static and
  interactive), faceting, group/rank filtering, custom colours, y-range,
  x-ordering, factor group columns, significance brackets, and
  `plot_enrichment_counts()`. Line coverage: 14.88% → 94.57%.

## Documentation

- New vignette `alpha-diversity` demonstrating the full alpha diversity pipeline:
  loading data, `compute_alpha()` (binary, threshold, and abundance modes),
  `plot_alpha()`, `compute_alpha_significance()`,
  `plot_alpha_significance()`, and significance brackets on box plots.

## Dependencies

- Added `ggsignif` and `rmarkdown` to `Suggests`.
- Added `VignetteBuilder: knitr`.

# phiper 0.3.1

## Changes to `compute_alpha`

- Added `pielou_evenness` and `berger_parker_dominance` to the output (NA for
  samples with richness ≤ 1 and richness = 0 respectively).
- New `metrics` parameter: request any subset of the five indices; defaults to
  all five.
- New `mode` parameter (`"binary"`, `"threshold"`, `"abundance"`) replaces
  `fc_threshold`. `"abundance"` mode uses raw values from `abundance_col` with
  optional `abundance_agg` (`"mean"`, `"sum"`, `"max"`) at higher ranks.
- `shannon_log` renamed to `shannon_base`; old name still works with a
  deprecation warning.
- Performance: all-samples roster now collected once before the rank loop
  instead of re-queried per rank.
- Hardening: all-invalid ranks now aborts with an informative error instead of
  silently returning empty output; `n_samples` attribute added to the result.
- Validation: `mode = "threshold"` now requires `threshold` to be finite and
  `abundance_col` (when supplied) to be a character scalar.

# phiper 0.3.0

## Major changes

- Extracted all data-import, class construction, and low-level utility code
  into the new **phiperio** package. `phiper` now declares `phiperio` as a
  hard dependency and re-exports its user-facing functions
  (`load_example_data`, `get_example_path`) so existing workflows require no
  changes.
- Removed all functions that moved to phiperio: `phip_convert`,
  `phip_convert_legacy`, `new_phip_data`, `expand_phip_data`,
  `phip_data_join`, `validate_phip_data`, `disconnect`, `get_comparisons`,
  `phip_example_path`, `phip_load_example_data`, and the full logging /
  validation helper suite (`add_quotes`, `word_list`, `.chk_*`).
- Internal logging and validation now use the unified phiperio helpers
  (`.ph_abort`, `.ph_warn`, `.ph_log_info`, `.ph_with_timing`,
  `.ph_check_cond`, `.ph_add_quotes`, `.ph_word_list`, `.ph_check_path`,
  `.ph_check_extension`, `.ph_check_null_default`).
- `get_peptide_meta()` renamed to `get_peptide_library()` throughout, in line
  with the phiperio API.
- `compute_alpha`: restored efficient same-connection peptide library
  handling via `.ph_peplib_on_main()` (DuckDB ATTACH fast path with
  `copy_to()` fallback).

## Internal

- R source files renamed to follow the new `<domain>_compute` /
  `<domain>_plots` convention: `binary_alpha` → `alpha_compute`,
  `binary_alpha_plots` → `alpha_plots`, `binary_beta` → `beta_compute`,
  `binary_beta_plots` → `beta_plots`, `prevalence-DELTA_test` →
  `delta_compute`, `prevalence-DELTA_plots` → `delta_plots`,
  `prevalence-POP-test` → `pop_compute`, `prevalence-POP_plots` →
  `pop_plots`, `plot-utils` → `plot_utils`. Test files renamed accordingly.
- Naming conventions for source files and functions documented in
  `CONTRIBUTING.md`: exported functions use plain `snake_case`; all internal
  helpers use the `.ph_` prefix.

# phiper 0.2.7

- compute_delta: added `maxmean` (Efron-type) as a test statistic option.
- compute_delta: added prevalence bins for the test statistic computation.
- compute_delta: `srlr` is now the default test statistic.
- Implemented McNemar test and paired signed root likelihood ratio statistic for
  paired designs.
- Welford's online algorithm refactored into its own class; added a wrapper for
  post-permutation output generation.
- Docs and R CMD CHECK fixes.

# phiper 0.2.6

- compute_delta: added stat_mode options "score" (pooled score z) and "srlr"
  (signed root likelihood ratio) using raw counts.
- compute_delta: removed smoothing, prevalence filtering, BH adjustment, and
  fold_change/cross_prev summaries; outputs now include T_null_mean and
  T_null_sd and standardized T_obs uses the null mean and sample SD.
- shift_computing: null variance uses Welford's algorithm with sample variance
  in permutation loops.
- prevalence-DELTA plots/tests updated to drop BH/categorical dependencies and
  use numeric-only filtering where applicable.

# phiper 0.2.5

- Removed 10 unused package dependencies (data.table, fs, htmltools, purrr, 
forcats, filelock, arrow from suggests)
- Moved duckdb and dbplyr to imports; moved knitr to suggests for proper 
dependency management
- Removed 360 lines of unused CAP/dispersion plotting functions and associated 
documentation

# phiper 0.2.4

- Updated the peptide library with the new annotations from Sasha
- in the compute_distance --> fallback to collecting before pivoting

# phiper 0.2.3

- Added `small_mixture` to `phip_load_example_data`
- Changed all examples and tests to use this `small_mixture` instead of redining it
- Added cache to `phip_load_example_data` to make tests/examples faster

# phiper 0.2.2

- Separated the feature associations to PCoA vectors from `compute_pcoa` to a separate function.

# phiper 0.2.1

## Minor changes

- Removed superassignments, changed all assignments to "<-" style
- Stated all dependencies, removed unstated dependencies

## Bug fixes

- Fixed typos in the examples from binary_beta.R
- fixed R CMD check (0 errors, 0 warnings, 0 notes)

## Documentation

- Documented empty params, got rid of all documentation-related wearnings, 
- removed all non-ASCII chars
- Added CONTRIBUTING.md to .Rbuildignore (it caused a note in R CMD CHECK),

# phiper 0.2.0

## Minor changes

- Removed generic and S3 methods for the expand_phip_data
- Renamed the internal helpers involved in the data import to match the naming
convention .ph_<function_name>
- documented the internal helpers involved in the data import stage

## Major changes

- Removed the backend argument entirely from the phip_convert, 
phip_convert_legacy, new_phip_data, .resolve_paths and other, less imporatant 
helpers. DuckDB is now the only supported backend
- Changed the structure of the repo. Now all functions related to the 
standard/legacy import workflows live in the standard-conver.R or
legacy-convert.R
- removed the resolve-paths.R file. Moved the functions to utils.R


# phiper 0.1.1

## Minor changes

- Removed old, dead and undocumented code: R/binary-analysis_peptides.R, 
R/vinary-analysis_stability.R, R/binary-analysis_stability_cattime.R,
R/fold_change-analysis.R

- i left other files, even if they were undocumented/untested, as they were 
essential for phiper to work (e.g. the data import paths)

# phiper 0.1.1

## Minor changes

- Removed old, dead and undocumented code: R/binary-analysis_peptides.R, 
R/vinary-analysis_stability.R, R/binary-analysis_stability_cattime.R,
R/fold_change-analysis.R

- i left other files, even if they were undocumented/untested, as they were 
essential for phiper to work (e.g. the data import paths)

## Documentation

- Regenerated the docs, removed the old functions from NAMESPACE

# phiper 0.1.0

## Minor changes

- None.

## Major changes

- None.

## Features

- None.

## Bug fixes

- None.

## Documentation

- None.

## Internal

- None.
