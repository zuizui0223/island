# Chapter 1 all-data probability model upgrade — 2026-09-15

Status: **exploratory model upgrade; not yet a replacement for the frozen v10 submission surface**.

## Motivation

The current Chapter 1 database starts from 8,265 geographic islands and 4,505 islands with observed GBIF floras. The frozen v10 analysis then becomes much smaller because formal floristic-status strata (`all_native`, `native_nonendemic`, `endemic`) are fail-closed. This is scientifically conservative, but it mixes two separate restrictions:

1. trait-evidence quality;
2. floristic-status resolution.

The present upgrade separates them explicitly.

## Evidence hierarchy

Primary trait evidence:

- `all_analysis_eligible` = High + Medium + trait-specific Validated Low.

Sensitivity:

- `direct_only` = High + Medium only.

No arbitrary minimum number of trait-scored species per island is imposed. For each trait outcome, the denominator is the number of species on that island with a resolved value for that outcome.

## Flora hierarchy

Broad sampling-frame analysis:

- `all_observed` = every observed island-flora record with usable trait evidence, irrespective of unresolved/native/introduced status.

Status-resolved sensitivities:

- `all_native`;
- `native_nonendemic`.

The `all_observed` result is a statement about **observed island-flora composition**. It is not by itself evidence of native colonisation filtering, historical evolution, or in-situ adaptation. Those interpretations require persistence in the status-resolved sensitivities.

## Probability model

For island *i* and trait outcome *j*:

\[
Y_{ij} \sim \mathrm{BetaBinomial}(N_{ij}, p_{ij}, \kappa_j)
\]

where `N_ij` is the number of trait-resolved species and `Y_ij` is the number in the predeclared focal state. The mean model is

\[
\mathrm{logit}(p_{ij}) = \alpha_j + \beta_{j,d} z(\mathrm{distance})
+ \beta_{j,a}z(\mathrm{area}) + \sum_{k=1}^4 \beta_{j,k}z(\mathrm{climatePC}_k).
\]

Each response has its own beta-binomial concentration parameter `kappa_j`. This allows extra-binomial island-to-island variation instead of relying only on a grouped-binomial working variance.

Inference uses spatial-block cluster-robust sandwich covariance. A direct North–Tropical test fits context interactions and uses a joint Wald test of the distance × context coefficients.

## Primary response vector

Accessibility/generalisation components:

- generalized floral form;
- actinomorphic symmetry;
- shallow/open tube.

Reproductive-assurance components:

- self compatibility;
- selfing mating system;
- autonomous selfing capacity.

The probability analysis therefore operates on six explicit species-count outcomes rather than treating a continuous syndrome score as if it were a count probability.

## Fixed interpretation rules

- `all_analysis_eligible` is the primary evidence scope.
- `direct_only` is an evidence-quality sensitivity.
- `all_observed` is the broad sampling-frame composition result.
- `all_native` and `native_nonendemic` test whether the broad result persists after floristic-status restriction.
- Missing trait data are never coded as trait absence.
- No pollinator identity enters the model.
- A significant result in one context and a nonsignificant result in another is not a context difference; the interaction-vector test is required.
- This branch does not retroactively change the frozen v10 claim ceiling until the new analysis and its sensitivity hierarchy are reviewed.

## Reproducibility

Implementation:

- `src/island_v2/chapter1_all_data_probability.py`
- `config/chapter1_all_data_probability.yml`
- `.github/workflows/run-chapter1-all-data-probability.yml`

Frozen upstream input remains Chapter 1 progressive-analysis run `34232450884`, artifact `chapter1-progressive-analysis-34232450884`.
