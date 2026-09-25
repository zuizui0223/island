# Current submission pipeline — corrected geography

The **only active Chapter 1 submission selector** is
`config/chapter1_submission_current.json`.

Submission-facing prose lives in `submission/chapter1_current/`.
Corrected tables live in `results/geography_20260924/`.
The replay entry point is `scripts/geography_correction/README.md`.

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
  beta-binomial GLM       assurance-adjusted models
  + multivariate Wald      + raw colour/architecture
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

## H1 — Pattern

Data: island trait prevalences.

Method:
`trait prevalence ~ isolation + island area + climate PC1-PC4`

Model: beta-binomial GLM with spatial-block cluster-robust covariance; seven-response multivariate Wald test within four predeclared regions.

Corrected joint q-values:
- N mid: `3.216e-10`;
- N high: `2.433e-5`;
- Tropical: `3.498e-7`;
- S extra: `2.504e-18`.

Joint support does not imply that every atomic trait is positive.

## H2 — Conditional decomposition

Data: reproductive assurance + floral responses.

Models:
- `reproductive assurance ~ isolation + area + climate`;
- `floral response ~ isolation + reproductive assurance + area + climate`.

Primary selfing-adjusted accessibility q-values:
- N mid: `0.1703`;
- N high: `0.000847`;
- Tropical: `0.01616`;
- S extra: `0.2873`.

Tropical Direct-only accessibility is not FDR-supported after correction (`q=0.1196`).

## H3 — Ecological pressure

Data: independent GloPL pollen-supplementation experiments.

Model:
`pollen limitation ~ isolation + geographic context + experimental design`

Weighted regression with publication-total weight one and SE clustered by publication.

Corrected primary estimate:
`beta=0.09191`, `SE=0.03806`, `p=0.01575`.

## H4 — Functional compatibility

Data: exact species overlap between H2 scores and GloPL.

Model:
`pollen limitation ~ trait score + isolation + geographic context + experimental design`

Corrected primary estimates:
- reproductive assurance: `beta=-0.29830`, `p=0.00396`;
- generalized accessibility: `beta=-0.29566`, `p=0.02187`.

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
