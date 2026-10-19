# phiper 0.5.0 (2026-10-19)

## New features

- `compute_delta2()` tests, per `(rank, feature, group_col)` stratum, a
  between-subject group term, a within-subject time term (`time_col`) and
  their interaction, on sample-level data: a subject's samples are never
  collapsed, so the number of samples per subject no longer acts as a hidden
  covariate. Each term is a score statistic pooled over peptides like
  `compute_delta()`, with a permutation null drawn once per `group_col` and
  shared by all strata. Groups with more than two levels get an omnibus test
  and pairwise tests with max-T adjustment (`$pairwise`); ordered levels can
  be tested as a trend (`group_scores`, `time_scores`).
- `propensity_formula` makes the group term and pairwise tests valid when
  group membership depends on covariates: labels are drawn from their
  conditional law given the covariates (rejection sampling for two groups,
  a pairwise-swap chain for more), and the statistic is adjusted for the
  same covariates. `time_propensity_formula` does the same for
  within-subject time labels. Effect sizes are reported raw and
  inverse-propensity weighted (`max_delta_adj_*`).
- Subject means in the group term are weighted by the within-subject
  correlation `rho`, estimated from all peptides by default.
- The result is a `delta2_result` with a print method: `$terms` (one row per
  stratum, term columns suffixed `_group`, `_time`, `_interaction`),
  `$pairwise` and `$meta` (design, propensity, dropout and sampler
  diagnostics, `rho_hat`).
- `inst/validation/` holds the simulation studies used as the release gate
  for `compute_delta2()` (calibration, power, unbalanced designs, `rho`) and
  their latest results.
- `standardize = "two_pass"` in `compute_delta2()` standardises each
  peptide's statistic by its mean and SD over the draws before pooling: a
  first pass records them over the observed labels and the `B_permutations`
  draws, a second replays the same draws. Every peptide then has mean 0
  and SD 1 under the null that was actually sampled (propensity draws, swap
  chain, biased time coins), at twice the permutation cost.
- The result gains `$peptides`, the per-level prevalence of every peptide
  (subject-mean prevalence per group level, raw and inverse-propensity
  weighted; share of positive samples per timepoint). `return_peptides =
  FALSE` leaves it out.
- `forestplot()`, `deltaplot()`, `ecdf_plot()` and their interactive
  versions take a `delta2_result` directly, with `term` (and `group_col`,
  `adjusted`, `group_pair_values`) choosing what to plot.
- The documentation and a new vignette, `vignette("delta2-analysis")`,
  explain how `weight_mode` chooses the effect scale (`equal`: arcsine,
  `n_eff_sqrt`: log odds ratio, `se_invvar`: prevalence difference), how
  prevalence strata temper it, and how to read the propensity diagnostics.

## Changes

- `compute_delta()`'s peptide-library resolution, `rank_feature_keep` filter
  and parallel runner moved, unchanged, into helpers shared with
  `compute_delta2()`. Results are unchanged.
- `nnet` is suggested, for the propensity model with three or more groups.

# phiper 0.4.22 (2026-10-19)

## Performance

- The `compute_delta()` permutation kernel parses its mode arguments once
  instead of comparing strings for every peptide, builds the `strat_bins`
  breaks once, and reuses its buffers across permutations instead of
  allocating them for each one. On x86 CPUs that support it, hits are
  counted with the POPCNT instruction, chosen at run time so generic builds
  still run on any CPU. Results are unchanged.
- Unpaired permutations are drawn once per group contrast and shared by all
  of its rank/feature strata, instead of being drawn again for every
  stratum. With 200 peptides in 20 species (220 strata) at the default 2000
  permutations, a cohort of 1,000 subjects now takes 0.29 s instead of
  3.1 s, and 10,000 subjects 2.2 s instead of 30 s. The vignette data at
  2000 permutations takes 1.2 s instead of 3.3 s.
- Each unpaired permutation now counts a peptide's hits in the first group
  only and takes the second group's count as the peptide's total minus
  that, which halves the counting work; the stored permutation splits take
  half the memory. One stratum of 200 peptides at 2000 permutations takes
  42 ms instead of 64 ms with 10,000 subjects, and 23 ms instead of 25 ms
  with 1,000, where computing the statistic dominates. Results are
  unchanged.
- `compute_delta()` makes fewer passes over the input table before the
  permutations. The duplicate-positives check compares the number of
  distinct keys with the number of positive rows instead of counting every
  key, an unused per-subject summary is gone, and the group levels and the
  subject rows of each group level come from one distinct subject x group
  table instead of two pivots over every row. With an in-memory table of
  10,000 subjects and 200 peptides (2 million rows) at the default 2000
  permutations, a run takes 0.46 s instead of 2.3 s. Results are unchanged.
- In unpaired strata with a single peptide (every `peptide_id` stratum),
  the statistic of a permutation depends only on the first group's hit
  count, so it is computed once per distinct count and reused. The
  vignette data at 2000 permutations, where 2,225 of the 2,604 strata hold
  one peptide, takes 0.85 s instead of 1.2 s. Results are unchanged.

## Changes

- Under the same seed, unpaired `p_perm`, `T_null_mean`, `T_null_sd` and
  `b` differ from 0.4.20, because all strata of a group contrast now share
  one seed and the same permutations. Each stratum's p-value is still a
  valid permutation p-value; `T_obs` and the other observed statistics are
  unchanged, and so are runs with only paired contrasts. The regression
  snapshot is updated for the unpaired cases.
- New test in `test-delta_compute.R`: a species holding a single peptide
  gets the same null as that peptide's `peptide_id` stratum.
- The `paired_by` documentation of `compute_delta()` now states that a
  contrast is paired as soon as one `paired_by` unit is observed in both
  groups, and that the paired test then uses only the units present in
  both groups.

## Bug fixes

- `compute_delta()` results are now reproducible with `set.seed()` for
  `phip_data` input. DuckDB returned subjects and library peptides in a
  different order on each run, and since the permutations draw subjects by
  position, the same seed gave a different null distribution (`p_perm`,
  `T_null_mean`, `T_null_sd`, `b`). Subject rows and peptide columns are
  now sorted before the permutation test, so results also no longer depend
  on the row order of `x` or `peptide_library`. Under the same seed,
  results can differ from 0.4.20 for input that was not already in
  subject and peptide order.
- New test in `test-delta_compute.R`: shuffling the rows of `x` and
  `peptide_library` leaves the output identical.
- `compute_delta()` now errors when the same `subject_id` is in both
  groups of an unpaired contrast, e.g. `group_cols = "timepoint"` without
  `paired_by`. A subject's hits are pooled into one row, so the subject
  was counted in both groups and a permutation could place it on both
  sides. The error lists up to 10 such subjects and points to `paired_by`.
- New tests in `test-delta_compute.R`: a within-subject timepoint contrast
  errors without `paired_by` and runs as paired with
  `paired_by = "subject_id"`. The zero-variance toy test, which reused its
  subjects in both groups, now has one subject per sample.

# phiper 0.4.21 (2026-10-19)

## Performance

- `compute_delta()` spends far less time in R per stratum. It now builds
  the peptide columns of every rank/feature and the subject rows of every
  group level once, before the loop, instead of filtering a table for each
  stratum. Each stratum's result row is a plain list, and all rows are
  bound into one tibble at the end. On the typical-workflow vignette data
  (2,604 strata) a run at `B_permutations = 100` drops from 8.5 s to 0.75 s,
  and at the default 2000 from 10.7 s to 3.3 s. Results are unchanged.

# phiper 0.4.20 (2026-10-19)

## Bug fixes

- `compute_delta()` no longer warns about an unexpected many-to-many
  relationship on paired designs. The join that checks which subjects appear
  in both groups is many-to-many by design and now says so.

## Tests

- New `test-delta_regression.R` snapshots the full `compute_delta()` output
  for every combination of `stat_mode`, `weight_mode`, `aggregate_stat`,
  `strat_bins` (0 or the default) and `perm_method`, on a fixed synthetic
  dataset of 300 subjects and 50 peptides (unpaired and paired). Later
  changes to `compute_delta()` must leave this snapshot unchanged. The
  dataset and the script that generates it are in `tests/testthat/fixtures/`.
- The test of the bitset/peptide dimension check in `test-delta_compute.R`
  now restores the mocked `build_bitset_unpaired()` when it ends. The mock
  used to stay in place for the rest of the test run.

# phiper 0.4.19 (2026-10-18)

## Bug fixes

- `scatter_interactive()`, `volcano_interactive()` and
  `plot_alpha_interactive()` showed the literal `%{text}` instead of the hover
  text for a point that was alone in its trace, e.g. the only BH-significant
  peptide of a contrast or the only peptide matching a `color_by` value.
  plotly writes the text of a one-point trace as a single string, which the
  `%{text}` hover template cannot resolve. The hover now uses
  `hoverinfo = "text"`, as the Delta plots already do.

# phiper 0.4.18 (2026-10-18)

## New features

- `profile_plot()` plots the peptides of selected proteins at their position
  along the protein and at their prevalence per group, one panel per protein,
  from `compute_pop()` output. It draws any number of groups. A shaded band
  marks the longest run of overlapping peptides with BH q <= `alpha`, per
  group pair and direction, and `delta_tbl` adds the `compute_delta()`
  statistics (`T_obs`, `p_perm`) to the panel titles. A section in the Delta
  vignette shows it on the example data.

## Dependencies

- phiper now requires phiperio >= 0.6.2, whose peptide libraries carry the
  `protein_seq_id` column that `profile_plot()` uses.

# phiper 0.4.17 (2026-10-17)

## Bug fixes

- `compute_delta()` no longer fails when `rank_cols` mixes column types, such
  as a character rank like `genus` and a logical flag like `is_auto` from the
  combined library ("Can't combine `genus` <character> and `is_auto`
  <logical>"). The rank columns are now converted to character before they are
  stacked, so a logical flag yields the features `"TRUE"` and `"FALSE"`.

# phiper 0.4.16 (2026-10-15)

## Peptide libraries

- `library-metadata/` now hosts the updated combined, human proteome and ICAM
  peptide libraries served by phiperio 0.6.2. All three carry the
  protein-level columns `protein_seq_id` and
  `taxon_with_protein_unique_annotation`, and the combined library is now the
  one of the phipflow pipeline, with its revised taxonomy. The previous files
  stay for phiperio 0.6.1 and earlier.

# phiper 0.4.15 (2026-10-08)

## Continuous integration

- The phipflow compatibility check now calls the workflow in its new home,
  `Vogl-Lab-Research/phipflow`.

# phiper 0.4.14 (2026-10-06)

## Reproducibility

- Plot functions no longer reset the global random seed (#65).
  `scatter_static()`, `scatter_interactive()` and `plot_alpha_interactive()`
  used to call `set.seed(1)` internally, which overwrote any seed the user
  had set. Jitter now follows the user's `set.seed()` and leaves the global
  random state unchanged, the same as `deltaplot()` and
  `deltaplot_interactive()`. Without a seed, jitter differs between sessions.
- The background subsample in `scatter_interactive()` still uses
  `background_seed`, but no longer changes the global random state.

# phiper 0.4.13 (2026-10-06)

## Dependencies

- phiper now imports 13 packages instead of 20, and installing it pulls in 41
  packages instead of 77 (#25). `cli`, `tidyselect`, `scales`, `showtext` and
  `sysfonts` are no longer needed. `plotly` and `Rtsne` moved to Suggests:
  install `plotly` to use the `*_interactive()` plots and the 3D t-SNE view,
  and `Rtsne` to use `compute_tsne()`. phiper now requires ggplot2 >= 4.0.0.

## Plotting

- Loading phiper no longer changes global settings: it does not set the
  default ggplot2 theme, the discrete colour options, or showtext font
  rendering. Plots still use the phiper colours, which `theme_phip()` now
  sets.
- The bundled Montserrat font was removed and plots default to the `"sans"`
  font family.
- Static plots now draw Δ, Σ and subscripts (e.g. in `deltaplot()` and
  `volcano_static()`) with plotmath, so they can be saved to PDF with the
  standard `pdf()` device.

# phiper 0.4.12 (2026-10-06)

## Documentation

- The "Typical workflow" vignette has clearer code comments, notes that
  `ggsignif` is needed for significance brackets in alpha diversity plots,
  and shows how to save a `phip_data` object with `export_parquet()`.
  Written by Nikolas Basler.

# phiper 0.4.11 (2026-10-05)

## Documentation

- The human proteome library article now documents the library's metadata
  fields: general peptide and source-protein information, taxonomic
  information, and the `is_*` annotation flags (proteome, mitochondrial, HLA,
  HLA eplet, control, neoantigen, cryptic peptide, transposable element ORF
  and therapeutic antibody sources). Written by Nicolai Hörstke.
- The function reference on the website is now grouped by module (alpha
  diversity, beta diversity, POP analysis, delta analysis, plot styling and
  example data).

## Continuous integration

- `pkgcheck` passes again under pkgcheck 0.2.0.44. The workflow now waits for
  the latest `R-CMD-check` and `test-coverage` push runs to finish before
  running pkgcheck, which previously sampled them mid-run and reported a
  failing CI. Push-triggered pkgcheck runs on different branches no longer
  cancel each other. Alon Alexander is now listed as copyright holder only, as
  pkgcheck requires an ORCID for every contributor.

# phiper 0.4.10 (2026-10-05)

## Continuous integration

- New `phipflow-compat` workflow runs the end-to-end smoke test of the
  [phipflow](https://github.com/csReynaB/phipflow) Nextflow pipeline with the
  phiper version under test, on every pull request and push to `main` that
  touches the package. Changes that break phipflow now fail before they are
  merged.

# phiper 0.4.9 (2026-10-05)

## Bug fixes

- `compute_alpha()` no longer corrupts phiperio's peptide-library cache. It
  attached the cache's DuckDB file to the data connection in read-write mode,
  so when `get_peptide_library()` later wrote another library to the cache in
  the same R session, that table became unreadable ("INTERNAL Error: invalid
  fsst_symbol_table_offset"). This broke the human proteome vignette in
  `R CMD check` on Linux and macOS. The cache is now attached read-only. An
  already corrupted cache is rebuilt with
  `phiperio::get_peptide_library(..., force_refresh = TRUE)`, or by deleting
  `phip_cache.duckdb` from `tools::R_user_dir("phiperio", "cache")`.

# phiper 0.4.8 (2026-10-05)

## Documentation

- The "Peptide library" navbar item is now a dropdown with one article per
  library. The existing article now covers the Agilent, Twist and Corona2
  library, and two articles are new: the ICAM library, documenting its
  general, taxonomic, annotation-flag and fused-protein fields (written by
  Gabriel Innocenti), and the human proteome library (curated by Nicolai
  Hörstke), whose field documentation is coming soon.

# phiper 0.4.7 (2026-10-05)

## Bug fixes

- When no peptide library is supplied, `scatter_interactive()`, the
  `color_by` highlighting of the scatter and volcano plots, and
  `compute_delta()` now fetch the phiperio libraries the peptide IDs belong
  to, detected from their prefixes with `phiperio::detect_peptide_libraries()`.
  They previously always fell back to the combined Agilent, Twist and Corona2
  library, so human proteome and ICAM peptides got no hover metadata and could
  not be highlighted or grouped by taxonomy. Peptides that match no known
  library get no library: plots leave every point uncoloured, and
  `compute_delta()` asks for `peptide_library` instead of reporting missing
  columns. Requires phiperio 0.6.1.

# phiper 0.4.6 (2026-10-05)

## Peptide libraries

- Added two peptide libraries to `library-metadata/`, next to the combined
  Agilent, Twist and Corona2 library: the human proteome library
  (`human_proteome_library_16.09.26.rds`, 300,000 `humanProteome_*` peptides)
  and the ICAM library (`icam_library_01.10.26.rds`, 420,000 `icam_*`
  peptides). phiperio's `get_peptide_library()` serves them as
  `"human_proteome"` and `"icam"`. Both omit the full protein sequence
  (`full_aa_seq`), which kept the files above GitHub's 100 MB limit; the
  curation scripts are in phiperio's `data-raw/`.

# phiper 0.4.5 (2026-09-09)

## Bug fixes

- `compute_delta()` no longer ignores `paired_by` in the strict hits guard
  (#56). The guard rejected duplicate positives per `subject_id` even when a
  different pairing column was supplied, so a subject contributing two samples
  to the same group aborted the call although the pairing unit was unique
  within that group. The guard now keys on `paired_by` when given, and its
  error message names the column it actually checked.

# phiper 0.4.4 (2026-09-08)

## New features

- `compute_delta()` gains a `min_m_eff` argument (#54). Strata whose effective
  number of peptides (`m_eff`) falls below the threshold are skipped **before**
  any permutation is drawn and are dropped from the returned tibble, since the
  permutation test is only reliable for `m_eff > 5`. Defaults to `0`, which
  preserves the previous behaviour of testing every stratum.

# phiper 0.4.3 (2026-07-14)

## New vignette

- Added "Peptide Library Metadata" vignette documenting the Agilent, Twist, and Corona2 peptide library metadata columns (documentation by Gabriel Innocenti, based on prior annotation work by Sasha Zhernakova and Carlos Reyna-Blanco).

## Documentation

- pkgdown navbar and Articles page now include a dedicated "Peptide library" entry.

# phiper 0.4.2 (2026-07-06)

## New vignette

- Added "Tutorial: Typical phiper workflow" vignette by Nikolas Basler, covering a full cross-sectional analysis (alpha/beta diversity, POP, DELTA) with a bundled dummy dataset (`inst/extdata/typical-workflow/`).

## Data

- Bundled tutorial dummy dataset (40 enrichment CSVs + metadata) under `inst/extdata/typical-workflow/`.
- Updated peptide library to `combined_library_06.07.26.rds` with `protein_id` merged in from `peptide_to_protein_map.csv`.

## Contributors

- Added Nikolas Basler as contributor.

## Documentation

- pkgdown articles section now groups vignettes into "Get started" (tutorial) and "Module vignettes".
- README links directly to the tutorial vignette.

# phiper 0.4.1 (2026-04-15)

## Bug fixes

- Fixed `n_peptides` bug in `compute_pop()`.

## Vignettes

- New vignettes: beta diversity, delta plots, and POP plots.
- Added calibration vignette to `.Rbuildignore`.

## Plotting

- Exported beta diversity plot functions.

## Documentation

- Updated README.

## CI

- Updated pkgcheck workflow to mirror phiperio behaviour.

# phiper 0.4.0

## New functions

- `compute_pop()` replaces the old `ph_prevalence_compare()`. Computes
  per-feature prevalence comparisons (Fisher / McNemar tests) across group
  pairs. Unused arguments were removed, internal filtering now delegates to the
  shared `.ph_filter_pairs()` helper, and imports were trimmed.
- `scatter_static()` — static ggplot2 prevalence scatter (percent1 vs percent2)
  with BH-corrected significance colouring and optional `color_by` highlighting.
- `volcano_static()` — static ggplot2 volcano (log2 ratio vs −log10 p) with
  configurable fold-change and p-value cutoffs and raw / BH p-value modes.
- `volcano_interactive()` — plotly equivalent of `volcano_static()`.

## Changes to `scatter_interactive()`

- Fixed a hover-text bug that caused incorrect peptide labels to appear in
  some multi-rank datasets.
- Updated peptide-library joining to follow the phiperio conventions used
  elsewhere in the package.
- Removed arguments that were no longer used after the peptide-library refactor.

## Plotting

- All plots now use `theme_phip()` as their base theme and the phiper discrete
  colour / fill scales, ensuring a consistent visual style across the package.
- `phip_use_montserrat()` has been removed. Font registration is handled
  automatically in `.onLoad()`.

## Examples and documentation

- Added self-contained `@examples` blocks to `scatter_static()`,
  `volcano_static()`, and `volcano_interactive()`.
- Fixed examples for `deltaplot()`, `deltaplot_interactive()`, `ecdf_plot()`,
  and `ecdf_plot_interactive()`: the old examples called `ph_prevalence_compare()`
  which no longer exists; they now use a minimal inline `data.frame` and require
  no external data.
- Updated prose descriptions in `scatter_static()` and `scatter_interactive()`
  that still referenced `ph_prevalence_compare()`.

## Tests

- New `test-pop_plots.R` covering `scatter_static()`, `scatter_interactive()`,
  `volcano_static()`, and `volcano_interactive()` (return types, pair / rank
  filtering, BH colouring, `color_by` interface, background overlay, and the
  internal `.build_color_group()` helper).
- Removed the two `phip_use_montserrat()` tests from `test-plot_utils.R` for
  the function that was deleted.
- Updated vdiffr reference snapshots for `deltaplot`, `forestplot`, and
  `ecdf_plot` to reflect the new `theme_phip()` styling.

## R CMD CHECK

- Added `delta_ratio`, `n01`, and `n10` to `utils::globalVariables()` in
  `zzz.R` to silence the "no visible binding" notes from `compute_pop()`.

## Vignettes

- Removed the pre-built `.html` vignette that was accidentally committed to
  the repository.

# phiper 0.3.4

## Changes to `compute_delta()`

- New `perm_method = "mid_p"` option: computes the mid-p corrected permutation
  p-value, which halves the contribution of ties
  (\eqn{p_{\mathrm{mid}} = (\#\{|T_b| > |T_{\mathrm{obs}}|\} + 0.5 \cdot \#\{|T_b| = |T_{\mathrm{obs}}|\}) / B}),
  reducing the conservative bias of the standard test for discrete statistics.
  `perm_method` now defaults to `"mid_p"` (previously `"standard"`).
- New `aggregate_stat = "af"` option: implements an adaptive Fisher aggregation
  statistic. Per-peptide z-scores are converted to one-sided p-values, sorted,
  and the optimal truncation point is selected by maximising a
  harmonic-mean-weighted cumulative \eqn{-\log p} sum. Positive and negative
  tails are aggregated separately and the dominant direction is returned.

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
