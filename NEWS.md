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
