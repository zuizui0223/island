# Current submission pipeline — corrected geography

The **only active Chapter 1 submission selector** is
`config/chapter1_submission_current.json`.

Submission-facing prose lives in `submission/chapter1_current/`.
Corrected H2-H4 tables live in `results/geography_20260924/`; final H1 inference lives in `results/h1_final_traitwise_t_20261004/`.
The H1 replay entry point is `docs/REPLAY_H1_FINAL_TRAITWISE_T_20261004.md`; geography/H2–H4 replay remains documented in `scripts/geography_correction/README.md`.

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

## H1 — Final separate-trait inference

### H1: seven separate trait responses

The primary analysis is broad contemporary island flora with All evidence. Each of seven binary traits is fitted separately in each of four regions by beta-binomial logit regression, adjusting for standardized corrected log isolation, log island area and climate PC1–4. Reproduction, colour and structure are labels, not aggregate scores. Spatial-block sandwich standard errors use a finite-cluster t reference with G−1 degrees of freedom and pointwise 95% intervals. All individual p values are two-sided and unadjusted; no Holm correction is applied. The final rule was user-approved after comparison of alternatives on 4 October 2026 and is retrospective, not preregistered. Individual significance does not establish family-wise significance or formal differences between regional slopes.

WCVP regional-native-compatible flora is the sole active origin sensitivity; Direct-only is a separate evidence-quality sensitivity. Both use the identical model and inference rule. WCVP compatibility retains existing source-native records and upgrades unresolved records only under accepted TDWG-L3 native-range compatibility; introduced records are not overwritten. This establishes regional compatibility rather than exact focal-island nativity. H2–H4 and their previously defined functional covariates are unchanged.

### H1: recurrent reproductive responses with regional floral differences

All 112 fits converged. Primary All evidence supports 17 of 28 individual associations at nominal two-sided P < 0.05, the same set of significant traits as the original poster's normal approximation. Self-compatibility increases with isolation in all four regions (P = 0.00569, 0.00671, 0.01812 and 0.01123 in northern mid-latitude, northern high-latitude, tropical and southern extratropical floras, respectively). Selfing mating system increases in the first three regions. Autonomous selfing increases in northern high latitudes and the tropics. Generalized form increases in northern mid-latitudes, the tropics and southern extratropics; actinomorphy increases in northern high latitudes and the tropics. Shallow/open tubes increase in northern high latitudes but decrease in southern extratropics. Plain colour increases in southern extratropics.

In WCVP All, 10 of 28 individual associations meet the same nominal threshold. Increasing self-compatibility remains supported in all four regions. Northern-high autonomous selfing, actinomorphy and shallow/open tubes remain positive and supported; tropical actinomorphy remains supported. Plain colour increases in both tropical and southern extratropical WCVP floras. Several broad-flora associations, including the southern shallow/open-tube decrease, no longer meet the threshold. This can reflect changed composition, precision and coverage, rather than proving an introduced-species mechanism. Direct-only yields 15/28 supported associations in broad flora and 11/28 in WCVP flora. Complete coefficients, intervals and unadjusted p values are reported in results/h1_final_traitwise_t_20261004/traitwise_results.csv.

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
