# H1 traitwise reanalysis — 4 October 2026

The directional score and three domain scores are absent. Seven binary outcomes remain individually interpretable. Broad contemporary flora is primary; WCVP regional-native compatibility is the sole active floristic-origin sensitivity. Direct-only is an evidence-quality sensitivity. H2–H4 were not refitted or changed by this H1 redesign. Their previously defined two functional covariates are not the removed omnibus H1 score.

## Design and interpretation

This is a user-requested post-result redesign, not a prospective confirmation. Beta-binomial logit models adjust for standardized corrected isolation, island area and climate PC1–4 in each region. Spatial-block sandwich uncertainty uses t(G−1); intervals are pointwise 95%. Two-sided Holm adjustment controls the 28-test family within each flora/evidence scope without assuming independence. Failed tests would retain their family slot. All 112 fits converged. No cross-trait or cross-region composite is tested.

Positive = increasing prevalence of the specified state with isolation. Negative = decreasing prevalence. A nonsignificant coefficient is uncertainty, not evidence of no response. Differences in significance between regions are not a formal region-by-isolation interaction test. Standardization is within each trait/region/scope; coefficient magnitudes are not identical-distance contrasts.

WCVP retains existing source-native rows and upgrades unresolved rows only with accepted TDWG-L3 native compatibility; explicit introduced records are never overwritten. This is not WCVP-only evidence nor exact island nativity. Coverage before trait/covariate filtering: 513,320 island-species records, 2,372 islands. Broad flora uses 1,039,757 ledger rows before filtering.

## Supported trait-region associations

| Evidence | Flora | Region | Trait | Slope | 95% CI | Holm P |
|---|---|---|---|---:|---|---:|
| all | broad | northern_midlatitude | selfing_mating_system | 0.04634 | [0.01807, 0.07460] | 0.0420058 |
| all | broad | southern_extratropical | plain_colour | 0.07419 | [0.03763, 0.11074] | 0.00559448 |
| all | broad | southern_extratropical | shallow_open_tube | -0.21660 | [-0.34046, -0.09275] | 0.0302621 |
| all | wcvp | tropical | plain_colour | 0.11908 | [0.04961, 0.18855] | 0.0259535 |
| all | wcvp | southern_extratropical | self_compatibility | 0.39945 | [0.17301, 0.62588] | 0.0263245 |
| all | wcvp | southern_extratropical | plain_colour | 0.16214 | [0.09373, 0.23055] | 0.000725947 |
| direct | broad | northern_midlatitude | selfing_mating_system | 0.06573 | [0.02777, 0.10369] | 0.0226099 |
| direct | broad | tropical | autonomous_selfing | 0.15420 | [0.07790, 0.23050] | 0.00299271 |
| direct | broad | tropical | generalized_form | 0.18933 | [0.09708, 0.28157] | 0.00249917 |
| direct | broad | southern_extratropical | plain_colour | 0.09349 | [0.05965, 0.12732] | 5.39634e-05 |
| direct | wcvp | northern_midlatitude | selfing_mating_system | 0.06729 | [0.02697, 0.10762] | 0.03698 |
| direct | wcvp | southern_extratropical | plain_colour | 0.18041 | [0.10182, 0.25901] | 0.00119985 |

## Biological reading

The most consistent corrected signal is increasing plain colour in southern extratropical floras (all four evidence/flora scopes). Northern-midlatitude selfing increases in broad all-analysis and both Direct-only flora scopes. Southern shallow/open tubes decrease in broad all-analysis, but that contrast does not pass Holm in Direct-only or WCVP scopes. Tropical autonomous selfing and generalized form are supported in broad Direct-only; tropical plain colour is supported in WCVP all-analysis. These scope-specific results must not be combined into a universal checklist. No northern-high-latitude individual contrast passes the current 28-test correction; this does not overturn H2, which tests a different conditional estimand.

Changes in support reflect both removal of aggregation and replacement of one-sided score inference with two-sided traitwise Holm inference. They do not imply that the input data or fitted isolation slopes changed.

## Reproducibility

See manifest.json for exact input SHA256 and runtime; traitwise_results.csv contains all 112 estimates, uncertainty, support and convergence. Historical pooled results remain in results/h1_final_directional_20261003 but are superseded for H1. No causal claim about within-lineage evolution or named pollinators follows from these assemblage associations.

Replay instructions: [exact inputs and command](../../docs/REPLAY_H1_TRAITWISE_20261004.md). Figure: [PDF](traitwise_H1.pdf), [PNG](traitwise_H1.png).

## Verification

All 112 fits converged; 21 focused tests passed. An independent read-only review recomputed Holm adjustments and checked table cardinality and input/config receipts. Its two documentation findings (obsolete WCVP aggregate wording and missing replay command) were corrected. H2–H4 source results were not modified. This closes the requested analysis replay, not journal submission packaging or the poster revision.
