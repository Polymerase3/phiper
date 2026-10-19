# Notes on the current results

`report.md` is generated; this file is the reading of it. Full run, B = 500,
1000 datasets per calibration cell, 4 cores, about 62 minutes. **45 of 50
gates pass, 1 is accepted as a documented limitation, 4 fail by chance (see
below).** This is the third full run; the fixes it confirms are listed at
the end. It predates the switch of the resampler to platform-independent
draws (after merging the same fix for `compute_delta()`), so a re-run with
the same seeds gives different replicates, equal in distribution.

Decisions (2026-10-19): the 6.3a miss is accepted and documented in
`?compute_delta2`; `rho_by_stratum` stays `FALSE` by default.

## The misses

**6.3a, controls missing 2 of 4 visits at random: residualised interaction
0.075 (band 0.0365-0.0635). Accepted.** Mildly liberal. Two further independent runs
of the same scenario gave 0.070 (600 datasets) and 0.051 (800 datasets),
about 0.066 pooled over 2400 datasets. Before the fixes below this scenario
was conservative (0.024). The likely remaining cause is that the
interaction's variance uses one residual variance per group, while
positivity here changes over visits and the groups see different visits;
a per-timepoint variance would address it, but a run without a time trend
(0.049) and one with it (0.051) did not separate the two cleanly. Accepted
as a known limitation rather than pursued further.

**6.2 null at 20/20 (0.074-0.080, four gates).** Chance: the four
statistics share each dataset and miss together, and a re-run with fresh
seeds gave 0.048 (0.008 at alpha = 0.01, 0.102 at 0.10). The same seeds
give the same miss in every run. These permutation tests are exact by
construction; the seeds are left as they are rather than re-drawn until the
gate passes. With 50 gates at 95% bands, two or three chance misses per run
are expected; a family-wise band (e.g. 99.9%, about +/- 3.3 SE) would make a
full run pass or fail as a whole, but that policy is not adopted yet.

## Reported, not gated

**6.3a, controls seen at T1-T2 only.** The interaction at T3-T4 is not
estimable (no control data there), so the observed statistic has fewer
degrees of freedom than its relabelled reference and the test is valid but
conservative (0 of 200). `compute_delta2()` warns about the gap in every
run (gated, passes).

**Three-group interaction: 0.058-0.062 in all four designs** (all within
the band; about 0.06 on average, so slightly above nominal).

## Passing

- 6.1: with the propensity model every term is calibrated in all six
  designs (group 0.039-0.059); without it the group term over-rejects
  (0.62-0.92) while time and interaction stay calibrated. This is the
  figure for the paper (`calibration.png`).
- 6.3a: without residualisation the interaction rejects in 99-100% of runs.
- 6.3b: informative dropout pushes the time term to 0.42; the diagnostics
  flag 100% of those runs and 0% of the clean ones.
- 6.4: unequal turnover with 10 vs 50 subjects, 0.055 / 0.057.
- 6.5: rho recovered with |bias| < 0.002; clamped at 0 in half the runs at
  rho = 0.
- 6.6 / 6.7: the test files pass (356 expectations).

## 6.2 and 6.5 conclusions

- Score stays the default: srlr's pooled null is offset by about 2.2 at
  10/30 and its power collapses in one direction (0.002 vs 0.63); the Wald t
  is lopsided the other way; score and standardised srlr are symmetric.
- rho = "pooled" is best at true rho = 0.2 (0.64 vs 0.56 "one", 0.62 "zero")
  and matches "one" at rho = 0.8 (0.32 vs 0.33). Per-stratum rho on the
  example data ranges 0.27-0.76 around a pooled 0.49, rising with
  prevalence; whether `rho_by_stratum = TRUE` gains power is not yet tested,
  so the default stays `FALSE`.

## Fixes found by this validation

1. The group term summed `w_i * ybar_i` (as specified in 3.2/3.5); with
   unequal weights that sum moves with every relabelling and shifts all
   peptides together (pooled null SD 1.53 instead of 0.45 at rho = 0). Now
   centred at the weighted mean.
2. With `propensity_formula` the group test never rejected (0.000-0.001):
   the propensity model, fitted on the labels it resamples, pins the
   statistic's covariate-driven part at the centre of its null. The group
   term and pairs are now adjusted for the propensity covariates (weighted
   least-squares score): 0.039-0.059.
3. The interaction used one pooled permutation variance, wrong when groups
   differ in within-subject variability: 0.109 / 0.025 with unequal
   turnover. It now uses one residual variance per group, shaped by each
   subject's visits and accounting for the estimated time effect: 0.055 /
   0.057.
