# Chapter 1 P3 joint observation-bias checkpoint

## Decision

P3 is complete under the prospectively frozen joint observation-bias contract.

The parent analyses were reproduced first on the pinned PR142 artifact:

- V5 trait-resolution MNAR reproduction: passed;
- V6 species-list detection reproduction: passed;
- baseline numerical reconciliation against the frozen PR142 branching outputs: passed.

Only after those checks passed was the joint surface opened.

Canonical execution:

- run: `34949880409`;
- artifact: `10389197309`;
- artifact name: `chapter1-joint-observation-bias-34949880409`;
- digest: `sha256:b5cfe366686f83975dd989275904d2e0ddfbf43be27ffbd503dd05439a5264b2`.

## Joint contract

The finite joint surface combines two distinct observation processes without making them interchangeable:

1. V5 state-dependent trait resolution acts on `selfing_core` among recorded flora;
2. V6 distance-dependent flora-list completeness plus state-dependent species recording acts on `generalized_accessible`.

The joint finite grid is:

- V5 shared trait-resolution OR_R: 9 levels;
- V6 median completeness C0: 5 levels;
- V6 distance-completeness OR_C: 5 levels;
- V6 state-recording OR_D: 7 levels;
- total: `9 × 5 × 5 × 7 = 1,575` surfaces per evidence scope.

All hypothetical or completed species retain the original observed regression information weight. The maximum no-differential score discrepancy was `3.33e-16`. Grid fractions are sensitivity-domain summaries, not probabilities over the true missingness process.

The partial-identification envelope combines the six deterministic V5 bound assignments with eight prespecified V6 corner settings, yielding 48 surfaces per evidence scope. These are adversarial identification bounds, not fitted latent truth.

## Finite joint surface

### Formal P2 contrast

The formal post-baseline H2/P2 comparison remains `northern_midlatitude` versus `tropical` within `analysis_regime`.

For the primary direct-only × native-nonendemic profile:

- robust cells: `1541 / 1575`;
- finite-grid robustness fraction: `0.9784`.

The 34 failures are concentrated under strongly adversarial joint bias, especially the strongest distance-dependent incompleteness (`OR_C = 0.25`) combined with severe positive-state under-recording and low trait-resolution odds. Thus the formal vector contrast is highly stable over the fixed finite joint domain, but not invariant to all allowed assumptions.

### Palearctic accessibility

Palearctic remains the strongest observation-robust within-context branch.

Native non-endemics:

- all-analysis: `1575 / 1575` robust;
- direct-only: `1575 / 1575` robust.

The finite-grid distance-slope ranges are:

- all-analysis: `0.0402` to `0.1381`;
- direct-only: `0.0392` to `0.1428`.

Thus the positive Palearctic accessibility direction is not erased anywhere in the prospectively frozen finite joint surface for the native-nonendemic stratum.

### Tropical accessibility

The tropical accessibility axis is much more observation-sensitive.

Native non-endemics:

- all-analysis: `1161 / 1575` robust (`73.7%`);
- direct-only: `1269 / 1575` robust (`80.6%`).

The finite-grid slope ranges cross zero in both evidence scopes. This axis must therefore be described as observation-fragile rather than as equally defended with the Palearctic branch.

## Partial-identification envelope

The partial-identification result is intentionally stricter than finite-grid robustness.

### Palearctic accessibility

Native non-endemics:

- all-analysis estimate envelope: `[0.0402, 0.1381]`; expected positive sign identified; support identified across all 48 corners;
- direct-only estimate envelope: `[0.0392, 0.1428]`; expected positive sign identified, but FDR support is not retained in every corner.

This is the strongest observation result in P3. The direction itself is identified across both evidence scopes under the frozen corner envelope; full support is identified only in all-analysis.

### Formal North–Tropical vector contrast

For direct-only native non-endemics, finite-grid robustness is high (`1541/1575`), but formal support is not retained across every deterministic partial-identification corner.

The correct claim is therefore:

> finite-domain robust, but not identified under the full deterministic observation envelope.

### Tropical accessibility

Native non-endemics:

- all-analysis envelope: `[-0.0856, 0.0399]`;
- direct-only envelope: `[-0.0942, 0.0391]`.

Both cross zero and lose support in some corners. The tropical accessibility direction is not partially identified.

## Robust versus fragile regions

The P3 result creates three distinct observation categories:

1. **Robust core — Palearctic accessibility.** Positive direction survives the entire finite joint surface for native non-endemics, and the deterministic envelope preserves the positive sign in both evidence scopes.
2. **Finite-grid robust / partially identified — formal North–Tropical vector difference.** The direct-only native-nonendemic contrast survives 97.8% of the frozen finite surface, but the deterministic envelope includes unsupported corners.
3. **Observation-fragile — tropical accessibility.** A substantial part of the finite surface breaks and the deterministic envelope crosses the null.

These categories are inferential labels over the prespecified assumption domain. They do not estimate how likely any bias mechanism is in nature.

## Manuscript consequence

The previous wording that reported V5 and V6 separately is now incomplete. The manuscript should state that separate robustness does not imply arbitrary joint robustness and should report the P3 result explicitly.

The defensible hierarchy is:

- Palearctic within-context accessibility is the strongest observation-defended plant-side pattern;
- the formal H2 North–Tropical vector contrast remains strongly finite-grid robust but is not point-identified under arbitrary deterministic joint bounds;
- the tropical accessibility component is observation-sensitive and must not be presented as equally robust;
- none of these observation analyses identifies a pollination mechanism or a latent true flora.

## Figure 4 consequence

Figure 4 should replace the V6-only observation panel with:

- the formal P2 joint robust/fragile surface;
- Palearctic versus Tropical finite-grid robustness;
- partial-identification envelopes;
- the existing geometry and H5 mechanism-identifiability boundaries.

The figure must explicitly state that grid-cell fractions are not posterior probabilities.

## Claim ceiling

P3 supports saying that the strongest Palearctic accessibility result is robust to a prospectively fixed joint class of trait-resolution and species-recording biases, and that the formal North–Tropical vector difference is highly stable over that finite class. It does **not** support saying that arbitrary missingness has been eliminated, that the formal context contrast is identified over all deterministic bounds, that tropical accessibility is robust, that true species-list completeness has been estimated, or that pollinator loss caused any observed response.
