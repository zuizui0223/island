# Chapter 1 v14 submission freeze — 2026-09-18

Status: **canonical v14 reproduced analysis surface**.

## Canonical reproduction

The Chapter 1 v14 analysis was independently rerun from the frozen upstream artifacts
and passed the fail-closed reproduction gate.

- scientific freeze head: `7793c6b2a4568a609df2095cc796c367252ea68f`
- v14 workflow run: `35314955780`
- artifact: `10535020072`
- artifact digest:
  `sha256:4665e68341cea35bfb16c33afeb30ccf2705b5a47982f81e409e8507204be811`
- workflow result contract: `chapter1_v14_reordered_hypotheses_result_v1`
- canonical lock: `config/chapter1_v14_canonical_result_lock.json`

The reproduction verifier confirmed:

- H1 reproduced;
- H2 reproduced;
- the exact-H2-score H4 bridge reproduced;
- the atomic-reconstruction H4 sensitivity reproduced.

The v13 audit also passed on the same scientific freeze head. v13 remains immutable
parent provenance.

## Canonical H1–H4 order

### H1 — global floral/reproductive island syndrome

H1 contains seven identically oriented atomic responses across three biological domains:

- reproductive assurance: self-compatibility, selfing mating system, autonomous selfing;
- colour dulling: `plain_colour`;
- accessibility/generalization: generalized form, actinomorphy, shallow/open tube.

The seven-response beta-binomial multivariate isolation vector is supported in all four
predeclared geographic replication strata in both evidence scopes. The equal-weight
three-domain descriptive orientation is positive in all four strata.

The colour component is not uniform: `plain_colour` is approximately flat in northern
mid-latitudes and positive in northern high-latitude, tropical and southern
extratropical strata. H1 therefore establishes recurrent multivariate syndrome
orientation, not identical atomic responses everywhere.

### H2 — selfing versus pollinator-facing floral decomposition

H2 asks whether the floral part of H1 is reducible to measured reproductive assurance.

Primary conditional decomposition:

- `selfing_core ~ isolation + area + climate`;
- `generalized_accessible ~ isolation + selfing_core + area + climate`;
- `plain_colour ~ isolation + selfing_core + area + climate`.

`selfing_core` is positive in all four contexts. The selfing-adjusted
`generalized_accessible` coefficient is also positive in all four contexts in both
evidence scopes.

Pollination-syndrome concordance is evaluated from **raw reported phenotype states**,
not weighted guild scores:

1. five raw flower-colour states after `selfing_core` adjustment;
2. raw colour × raw form/tube joint prevalence;
3. `P(raw architecture | raw colour, architecture resolved)` after selfing adjustment.

The raw five-colour vector changes with isolation in northern mid-latitude, tropical
and southern extratropical floras. Colour-conditioned architecture shows replicated
loss of blue/purple specialized/deep architecture in northern high latitudes,
increased yellow/orange deep-tube coupling in tropical Direct evidence, and mixed
southern restructuring.

Weighted `large_bee_like / butterfly_like / bird_like` scores are historical
secondary summaries only. The raw-state analyses do not identify realized pollinator
identity or direct causal selection by a named pollinator.

### H3 — independent global pollen-limitation pressure

H3 is the frozen full-global GloPL distance analysis, formerly v13 H2. The statistical
model and estimate are unchanged:

- 2,969 experimental effects;
- 1,248 sites;
- 919 publications;
- standardized distance coefficient `+0.079368`;
- SE `0.037734`;
- two-sided `p=0.03543`;
- one-sided positive `p=0.01772`.

H3 establishes an independent isolation-associated pollen-limitation gradient. It does
not establish mediation of H1/H2.

### H4 — post-hoc functional bridge

H4 has separate inferential layers.

#### Atomic-trait discovery

The frozen v13 exact-species bridge remains post-hoc. Autonomous selfing is the
strongest atomic association with lower current pollen limitation
(`beta=-0.446724`, `p=2.60e-8`).

#### Exact H2 species-score discovery

The primary v14 family-aligned bridge uses the **literal Direct-only H2 species scores**
rather than reconstructing them from atomic H4 traits.

- `selfing_core`: 455 species / 409 publications / 734 cells;
  `beta=-0.29706`, two-sided `p=0.00411`.
- `generalized_accessible`: 143 species / 143 publications / 259 cells;
  `beta=-0.29601`, two-sided `p=0.0222`.

Thus the same species-level H2 response scores that increase with island isolation are
associated with lower current experimental pollen limitation in the post-hoc GloPL
overlap.

#### Atomic reconstruction sensitivity

Reconstructing comparable family scores from frozen Route A/B atomic states yields the
same direction. This is sensitivity evidence, not the primary H2-to-H4 bridge.

No global primary colour-to-pollen-limitation bridge is promoted because the H1/H2
colour response is context dependent and does not define one globally directional
functional prediction.

## Supplementary prospective validation audit

The post-2015 wild-plant and PolLimCrop attempts are retained as a **non-promoted
validation audit**, not as part of H4.

- both were prospectively frozen before outcome access;
- both stopped at predeclared outcome-blind support gates;
- no prospective pollen-limitation outcome was opened;
- the audit documents the missing confirmatory evidence and does not alter the H4 result.

Source of truth:
`config/chapter1_h4_prospective_validation_summary_lock.json`.

## Claim ceiling

v14 may claim:

- a recurrent seven-response floral/reproductive island-syndrome direction;
- partially separable reproductive-assurance and pollinator-facing floral responses;
- context-dependent raw colour and colour–architecture reorganization not reducible
  to measured reproductive assurance;
- an independent global increase in experimental pollen limitation with isolation;
- post-hoc functional compatibility of both literal H2 species scores with lower
  current pollen limitation.

v14 must not claim:

- causal mediation from pollen limitation to the global syndrome;
- that the H2 conditional models identify direct selection;
- a universal bee, butterfly or bird mechanism;
- global pollinator abundance or visitation decline;
- historical trait selection by pollen limitation;
- that either supplementary prospective audit is a biological null;
- within-lineage evolution rather than assemblage composition.

## Canonical files

Read in this order:

1. `docs/chapter1_manuscript_full_v14_reordered_hypotheses_20260918.md`;
2. `config/chapter1_v14_canonical_result_lock.json`;
3. `docs/chapter1_submission_freeze_v14_20260918.md`;
4. `config/chapter1_v14_hypothesis_architecture.yml`;
5. `config/chapter1_v14_h2_decomposition.yml`;
6. `config/chapter1_v14_h4_evidence_hierarchy.yml`;
7. `config/chapter1_v14_h4_exact_h2_score_bridge.yml`;
8. `docs/PAPER_PIPELINE.md`.

The v13 manuscript, result locks and submission freeze remain unchanged parent
provenance below this surface.
