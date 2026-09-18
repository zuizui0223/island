# Chapter 1 v14 preflight results — 2026-09-18

Status: **local reproduction from the exact frozen input artifacts used by the v14
workflow; CI artifact still pending**.

This note is not a result lock. It records the numerical preflight used to check
whether the H1/H2 reorganization changes the scientific conclusion before the PR
workflow is promoted.

## Frozen inputs

- progressive-analysis run: 34232450884
- all-data primary run: 35100991898
- H1 model: beta-binomial logit, spatial-block cluster-robust sandwich covariance
- H1 covariates: log island area + climate PC1–PC4 + standardized log distance
- H2 continuous models: equal-island OLS with spatial-block cluster-robust covariance
- H2 plain-colour model: beta-binomial logit with selfing_core + H1 covariates
- evidence scopes: all-analysis-eligible primary and Direct-only sensitivity
- flora scope reported here: all_observed

## H1 — seven-response island syndrome

Adding plain colour to the six v13 atomic responses does **not** remove the
multivariate isolation response. The seven-response joint Wald test remains supported
in all four predeclared geographic replication strata in both evidence scopes.

### All-analysis-eligible

| context | joint p | BH q |
|---|---:|---:|
| northern mid-latitude | 1.77e-6 | 2.36e-6 |
| northern high-latitude | 2.61e-6 | 2.61e-6 |
| tropical | 1.89e-7 | 3.78e-7 |
| southern extratropical | 4.01e-13 | 1.60e-12 |

### Direct-only

| context | joint p | BH q |
|---|---:|---:|
| northern mid-latitude | 0.00461 | 0.00461 |
| northern high-latitude | 4.99e-9 | 9.98e-9 |
| tropical | 4.42e-6 | 5.90e-6 |
| southern extratropical | 2.73e-24 | 1.09e-23 |

### Three biological domains

All-analysis-eligible standardized isolation-coefficient means:

| context | reproductive assurance | colour dulling | accessibility/generalization | equal-domain orientation |
|---|---:|---:|---:|---:|
| northern mid-latitude | +0.0190 | -0.0001 | +0.0200 | +0.0130 |
| northern high-latitude | +0.0730 | +0.0188 | +0.1783 | +0.0900 |
| tropical | +0.1838 | +0.0367 | +0.0810 | +0.1005 |
| southern extratropical | +0.1490 | +0.0808 | +0.0175 | +0.0824 |

Direct-only equal-domain orientation is also positive in all four contexts:
+0.0159, +0.0925, +0.0854, +0.0855.

The important qualification is colour. The beta-binomial plain-colour coefficient is
approximately zero in northern mid-latitudes and positive in the other three contexts.
Therefore v14 supports a recurrent **three-domain syndrome direction overall**, but it
does not claim a statistically uniform colour-dulling coefficient in every region.

## H2 — selfing versus pollinator-facing floral decomposition

### H2a reproductive-assurance route

The selfing_core isolation coefficient is positive in all four contexts under both
evidence scopes. In the all-analysis scope it is strongest in tropical and southern
extratropical islands; northern mid-latitude is positive but imprecise.

### H2b accessibility after conditioning on selfing_core

The generalized_accessible isolation coefficient remains positive in all four contexts
after conditioning on selfing_core.

All-analysis-eligible:

| context | conditional distance beta | p | primary-H2b q |
|---|---:|---:|---:|
| northern mid-latitude | +0.0214 | 0.0833 | 0.133 |
| northern high-latitude | +0.1248 | 0.000360 | 0.00144 |
| tropical | +0.0733 | 0.00185 | 0.00494 |
| southern extratropical | +0.0536 | 0.1085 | 0.145 |

Direct-only has the same positive sign in all four contexts; northern high-latitude is
FDR-supported, tropical is borderline after the eight-test H2b correction.

The secondary attraction_shift coefficient is also positive in all four contexts in
both evidence scopes, with the clearest support in northern high-latitude and tropical
islands.

### H2b colour after conditioning on selfing_core

Plain colour is more heterogeneous after selfing adjustment.

All-analysis-eligible:

| context | conditional distance beta | p | primary-H2b q |
|---|---:|---:|---:|
| northern mid-latitude | -0.0026 | 0.894 | 0.894 |
| northern high-latitude | -0.0061 | 0.749 | 0.856 |
| tropical | +0.0611 | 0.0251 | 0.0501 |
| southern extratropical | +0.0803 | 0.000115 | 0.000924 |

Direct-only retains a strong positive southern extratropical effect
(beta=+0.0937, p=3.75e-7, q=3.00e-6), while the other contexts are not supported.

## Interpretation

The reordering is empirically coherent:

1. **H1:** a recurrent island-syndrome response remains after colour is included.
2. **H2a:** reproductive assurance increases with isolation.
3. **H2b:** the structural accessibility/generalization response is not reducible to
   selfing_core; this is the strongest plant-side evidence for a parallel
   pollinator-facing route.
4. **Colour is part of H1 but not the strongest H2 mechanism marker.** It shows a
   broad positive tendency in three contexts, with a particularly strong
   selfing-independent signal in the southern extratropical stratum.
5. **H3 and H4 remain logically downstream:** H3 independently tests whether pollen
   limitation increases with isolation; H4 tests trait–pollen-limitation functional
   compatibility.

The conditional H2 models do not prove causal mediation or direct pollinator
selection. They reject the narrower explanation that the measured selfing core alone
accounts for the full floral-architecture response.
