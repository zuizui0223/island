# Current submission pipeline — corrected geography

The **only active Chapter 1 submission selector** is
`config/chapter1_submission_current.json`.

Submission-facing prose lives in `submission/chapter1_current/`.
Corrected H2-H4 tables live in `results/geography_20260924/`; final H1 inference lives in `results/h1_final_directional_20261003/`.
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
  directional 1-df score  assurance-adjusted models
  + H1a/H1b synthesis      + raw colour/architecture
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

## H1 — Directional tendency + heterogeneity

Confirmatory data: the seven frozen v14 indicators, all oriented so positive means the classic island-syndrome direction.

Primary score gives equal total weight to reproductive assurance, colour dulling and accessibility/generalization.

Regional models retain the corrected beta-binomial specification and spatial-block cluster-robust covariance, but final inference is one-dimensional:
- finite-cluster t with G−1 df;
- Rademacher wild-cluster sign-flip sensitivity;
- strict four-region recurrence by intersection-union test.

Paper-level synthesis separates:
- **H1a:** positive global-average direction using Paule-Mandel random effects and modified Hartung-Knapp inference;
- **H1b:** regional heterogeneity using Cochran Q / I2.

Current result:
- H1a supported: all-analysis mean 0.0691, one-sided p=0.0256; Direct-only mean 0.0635, p=0.0238;
- H1b supported: I2=0.819 / 0.685;
- strict four-region recurrence unsupported because northern mid-latitudes are weak (IUT p=0.134 / 0.106).

The three raw measurement axes remain a descriptive phenotype/provenance audit. Their high-dimensional direction-free omnibus tests are not confirmatory evidence for the classic island-syndrome direction and cannot rescue a failed H1 directional endpoint.
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
