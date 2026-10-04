# Current submission pipeline — corrected geography

The **only active Chapter 1 submission selector** is
`config/chapter1_submission_current.json`.

Submission-facing prose lives in `submission/chapter1_current/`.
Corrected H2-H4 tables live in `results/geography_20260924/`; final H1 inference lives in `results/h1_traitwise_20261004/`.
The H1 replay entry point is `docs/REPLAY_H1_TRAITWISE_20261004.md`; geography/H2–H4 replay remains documented in `scripts/geography_correction/README.md`.

## Pipeline at a glance

```text
GSHHG 2.3.7 geography + GBIF island floras + trait evidence
                         |
                         v
             corrected 8,264-unit universe
                         |
                         v
               global island trait data
                         |
              +----------+----------+
              |                     |
              v                     v
        H1 pattern              H2 pathway
  seven separate traits  assurance-adjusted models
  + WCVP sensitivity      + raw colour/architecture
              |                     |
              +----------+----------+
                         |
                         | exact species scores
                         |                     independent GloPL
                         |                           |
                         |                           v
                         |                     H3 pressure
                         |                 weighted regression
                         |                           |
                         +-------------+-------------+
                                       |
                                       v
                                  H4 function
                         exact-species matched regression
```

## Geography correction

The corrected exposure uses source-matched GSHHG 2.3.7 continental coastlines and minimum minor-great-circle arc separation on a mean-radius sphere.

- 1,113 spurious H1 island zero distances were repaired;
- one Eurasian continental split component was removed;
- retained analysis universe = **8,264**;
- broad H1 union = **4,379**;
- 996 true continental GloPL zero-distance sites remain zero.

Source of truth:
- `config/chapter1_submission_current.json`;
- `docs/chapter1_corrected_submission_20260924.md`;
- `results/geography_20260924/`.

## H1 — Individual trait responses

H1 now tests seven previously defined binary outcomes separately: self-compatibility, selfing mating system, autonomous/delayed selfing, plain colour, generalized form, actinomorphy and shallow/open tubes. Reproduction, colour and structure are organizational domains only; neither domain scores nor an omnibus directional score are calculated. This user-requested revision on 4 October 2026 follows inspection of earlier results and is explicitly retrospective.

For each trait and region, counts of species expressing the state are modelled against their informative species denominator using a beta-binomial logit regression. Predictors are standardized corrected log isolation, log area and climate PC1–4. Spatial-block sandwich uncertainty uses a t reference with G−1 degrees of freedom; pointwise 95% intervals accompany every slope. Two-sided Holm adjustment covers all 28 region-by-trait tests within each flora/evidence scope, retaining failed tests in the family. Positive slopes indicate increasing prevalence, not an imposed direction of evolution. Cross-region differences in significance do not constitute a formal interaction test.

Broad contemporary flora is primary. The sole active floristic-origin sensitivity is WCVP regional-native compatibility; Direct-only is separately retained as an evidence-quality sensitivity. Historical pooled-score, strict-native, complementary-origin and high-dimensional omnibus analyses remain archived, not active evidence for H1. H2–H4 retain their distinct estimands and existing functional covariates.

All 112 individual-trait models converged. In broad all-analysis, three of 28 associations survived two-sided Holm correction: selfing mating system increased in northern mid-latitudes (slope 0.04634, adjusted P = 0.04201); plain colour increased in southern extratropical floras (0.07419, P = 0.00559); and shallow/open tubes decreased in the latter region (−0.21660, P = 0.03026). These are associations of assemblage composition with isolation.

The southern plain-colour increase was supported in every evidence/flora scope. Northern-midlatitude selfing was also supported in both Direct-only flora scopes. Other results were scope-specific: broad Direct-only supported tropical autonomous selfing and generalized form; WCVP all-analysis supported tropical plain colour and southern self-compatibility. The southern shallow/open-tube decrease did not survive correction in Direct-only or WCVP analyses. No northern-high-latitude individual contrast passed this 28-test correction. All coefficients and pointwise intervals, including unsupported and opposing estimates, are reported in results/h1_traitwise_20261004/traitwise_results.csv. These tests replace, rather than supplement, the historical omnibus directional score.

## H2 — Conditional decomposition

Data: reproductive assurance + floral responses.

Models:
- `reproductive assurance ~ isolation + area + climate`;
- `floral response ~ isolation + reproductive assurance + area + climate`.

Primary selfing-adjusted accessibility finite-cluster q-values:
- N mid: `0.1764`;
- N high: `0.00243`;
- Tropical: `0.01913`;
- S extra: `0.2977`.

Tropical Direct-only accessibility is not FDR-supported after finite-cluster correction (`q=0.1267`).

## H3 — Ecological pressure

Data: independent GloPL pollen-supplementation experiments.

Model:
`pollen limitation ~ isolation + geographic context + experimental design`

Weighted regression with publication-total weight one and SE clustered by publication.

Corrected primary estimate:
`beta=0.09191`, `SE=0.03806`, finite-publication `p=0.01594`.

## H4 — Functional compatibility

Data: exact species overlap between H2 scores and GloPL.

Model:
`pollen limitation ~ trait score + isolation + geographic context + experimental design`

Corrected primary estimates:
- reproductive assurance: `beta=-0.29830`, finite-publication `p=0.00417`;
- generalized accessibility: `beta=-0.29566`, finite-publication `p=0.02334`.

H4 is explicitly post-hoc functional triangulation.

## Submission package

Use:
- `submission/chapter1_current/MANUSCRIPT.md`;
- `submission/chapter1_current/FIGURE_CAPTIONS.md`;
- `submission/chapter1_current/COVER_LETTER_DRAFT.md`;
- `submission/chapter1_current/SUBMISSION_CHECKLIST.md`.

## Historical replay surfaces

### Superseded v14 surface
The following are preserved only for provenance and exact replay:
- `config/chapter1_v14_canonical_result_lock.json`;
- `docs/chapter1_manuscript_full_v14_reordered_hypotheses_20260918.md`;
- `docs/chapter1_submission_freeze_v14_20260918.md`;
- `.github/workflows/run-chapter1-v14-reordered-hypotheses.yml`.

### Frozen v13 parent paper pipeline
v13 remains parent provenance. It is not the current submission surface.

### Historical alpha1 database
The alpha1 bundle retains its original 8,265-unit contract for exact database replay. This does not restore the excluded continental component to the corrected analysis.
