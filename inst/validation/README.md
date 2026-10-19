# compute_delta2() validation

The release gate for `compute_delta2()`: simulation studies of calibration,
power and robustness, plus the automated brute-force and equivalence tests.
Nothing here runs on CRAN or during `R CMD check`.

## Running it

```sh
# from this directory; phiper installed, or point PHIPER_SOURCE at a checkout
PHIPER_SOURCE=/path/to/phiper VAL_CORES=8 Rscript run_validation.R
Rscript run_validation.R --quick                 # minutes, for a smoke test
Rscript run_validation.R --sections 1,3 --out my_results
```

Options: `--sections` (any of `1,2,3,4,5,tests`), `--quick` (few replicates),
`--out` (results directory, default `results/`), `--mixture` (path to
`phip_mixture.parquet` or a CSV export of it, for 6.5c; found automatically
in the package when `arrow` or `nanoparquet` is installed). Environment:
`PHIPER_SOURCE` loads phiper with `pkgload` from a source checkout (needed
for the 6.6 / 6.7 tests, which live in `tests/testthat`), `VAL_CORES` sets
the number of forked workers.

A full run (1000 datasets x B = 500 for the calibration studies) takes about
an hour on 4 cores. The run writes `report.md` (gates, tables, figures),
`gates_*.csv`, the figures as PNG and `results.rds` with every replicate.

## Data

`data/phip_mixture_prevalence.csv` holds the per-peptide prevalence of the
phiper example data (`inst/extdata/phip_mixture.parquet`, 1950 peptides). The
simulator draws baseline prevalences from it, or from a sparse log-uniform
profile on [0.01, 0.3] where a section asks for realistic sparsity.

The simulator (`val_sim()` in `R/helpers.R`) gives each subject a covariate
`age`, assigns groups at fixed sizes or from a multinomial logit in `age`,
and gives each subject-peptide a persistent state that a sample copies with
probability `sqrt(rho)`, so repeat samples of a subject correlate at `rho`.
Visits can be complete, depend on group, or stop after a high-breadth visit
(informative dropout).

## Sections and gates

| Section | What | Gate |
|---|---|---|
| 6.1 `R/s1_calibration.R` | {2, 3 groups} x {1, 2, 4 timepoints} under a null where group labels and positivity both depend on `age`; with and without `propensity_formula = ~ age` | with the propensity model, every term's rejection rate at 0.05 within the Monte Carlo band; without it, time and interaction within the band (the group term's over-rejection is reported, not gated) |
| 6.2 `R/s2_pooled_power.R` | score z vs srlr, standardised srlr, Wald t: same data, same permutations, 20/20 and 10/30 groups, sparse prevalence | report only (the null rate of each statistic is gated, since all are permutation tests) |
| 6.3 `R/s3_unbalanced.R` | (a) all 4 visits in cases vs 2 random visits in controls, with a time main effect (and, reported, controls seen at T1-T2 only); (b) informative dropout | (a) residualised interaction within the band, unresidualised over-rejects, the T1-T2-only design warned about; (b) no-dropout time term within the band, dropout flagged by the diagnostics in >= 80% of runs, diagnostics silent without dropout |
| 6.4 `R/s3_unbalanced.R` | one group with higher turnover, very unequal sizes | interaction within the band |
| 6.5 `R/s5_rho.R` | rho recovery, power of `rho` options, rho per prevalence stratum on the example data | \|bias\| <= 0.05 for 0 < rho < 1 |
| 6.6 / 6.7 `tests/testthat/` | swap sampler vs brute-force enumeration; the 3.2 equivalence table at 1e-10; constant propensities = uniform shuffles | the test files pass |

The Monte Carlo band is `0.05 +/- 1.96 * sqrt(0.05 * 0.95 / n)`, about
0.036-0.064 for 1000 replicates. With 50 gates, one or two band misses in a
full run are expected by chance alone; a systematic miss (one design or term
off in every run, or far outside the band) is a failure.
