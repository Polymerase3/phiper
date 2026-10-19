# Notes on the 2026-10-19 run

`report.md` is generated; this file is the reading of it. Full run, B = 500,
1000 datasets per calibration cell, 4 cores, about 52 minutes. 35 of 49 gates
pass. The misses, by cause:

**Fixed during validation.** The group term's subject values were `w_i * ybar_i`
(as specified in 3.2/3.5). With unequal weights (unequal sample counts,
rho < 1), `sum over group of w_i` moves with every relabelling and shifts all
peptides together: the pooled null SD was 1.53 instead of 0.45 at rho = 0,
and power collapsed. The values are now centred at the weighted mean,
`w_i * (ybar_i - mu_w)`, identical when the weights are equal. This run is
after the fix. Test: "unequal sample counts alone carry no group signal".

**Open: the group term with `propensity_formula` is conservative (6 gates,
rejection 0.000-0.001).** The propensity model is fitted on the labels it
then resamples; the logistic MLE makes the observed labelling's covariance
with the covariates equal its expected value, and the pooled group statistic
is dominated by that covariate-driven component, so `T_obs` sits at the
centre of its null. Prototype (two groups, cross-sectional, 400 datasets):
residualising each peptide's subject values on the propensity covariates
before the group comparison gives 0.060 under the null and power 0.32 at an
effect where the current statistic has 0.003.

**Open: the interaction is not robust to unequal within-subject variance
(6.4: 0.109 / 0.025; 6.3a: 0.024; K = 3, T = 2: 0.068 / 0.063).** Its null is
about mean changes, and with unequal turnover or visit patterns between
groups the subjects are not exchangeable: the classic Behrens-Fisher pattern
(anti-conservative when the small group varies more). Prototype (6.4
scenario, 600 datasets): a Welch-type studentised statistic, recomputed per
draw, gives 0.052 / 0.060, and 0.043 in a balanced design (same as now).

**Chance: 6.2 null at 20/20 (0.074-0.080).** The four statistics share each
dataset, so they miss together. A re-run with fresh seeds gave 0.048 at
alpha = 0.05 (0.008 at 0.01, 0.102 at 0.10).

**Passing, as expected.** Time and interaction under a subject-level
confounder, with and without the propensity model; the residualisation of
6.3a (without it the interaction rejects in every run); dropout diagnostics
(flag 100% of informative-dropout runs, 0% otherwise; informative dropout
pushes the time term to 0.42, which the diagnostics are there to warn about);
rho recovery (|bias| < 0.002, clamped at 0 in half the runs at rho = 0);
6.6 / 6.7 tests.

**6.2 conclusion.** srlr's pooled null is offset by about 2.2 at 10/30 and
its power collapses in one direction (0.002 vs 0.63); the Wald t is
asymmetric the other way (0.26 vs 0.63 at 30/10); the score z and the
standardised srlr are symmetric and close to each other. Score stays the
default.

**6.5 decision input.** Power: "pooled" is best at true rho = 0.2 (0.64 vs
0.56 "one", 0.62 "zero") and matches "one" at rho = 0.8 (0.32 vs 0.33).
Per-stratum rho on the example data ranges 0.27-0.76 around a pooled 0.49,
rising with prevalence.
