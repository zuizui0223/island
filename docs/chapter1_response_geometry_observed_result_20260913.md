# Chapter 1 observed response-geometry result — 2026-09-13

## Decision

The nonlinear response-geometry extension does **not** promote a global breakpoint, threshold, hinge or reversal claim for Chapter 1.

All 12 cross-scope cells are classified as `monotonic_or_unresolved` under the frozen observed-data gate. No cell is promoted to `identified_step` or `identified_reversal`.

This is the result of a staged design audit, not a post-hoc failure to find a convenient breakpoint.

## Frozen analysis sequence

1. `G0_flat`, `G1_cline`, `G2_step`, `G3_hinge` and `G4_reversal` were declared before opening any new observed nonlinear result.
2. V1 showed that naive AICc geometry selection was unsafe under the observed spatial-block structure, so observed geometry remained closed.
3. V2 used separate calibration and validation simulation seeds and retained the real island exposure distribution, trial counts, baseline covariates, spatial blocks and response-specific missingness.
4. V2 qualified 11/12 cross-scope design cells for `G2_step`, 0/12 for `G3_hinge` and 4/12 for `G4_reversal`.
5. Only after that qualification table was frozen was the observed geometry opened once.

## Frozen inputs and receipts

### Final Chapter 1 input

- workflow run: `34232450884`
- artifact: `chapter1-progressive-analysis-34232450884`
- artifact ID: `10058653212`
- digest: `sha256:b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`

### V2 recovery audit

- workflow run: `34764755666`
- artifact: `chapter1-response-geometry-v2-34764755666`
- artifact ID: `10320410517`
- digest: `sha256:1374fa264b9b94eb401a24d108a5a18d4beb6fede04d5622faf81569afc27838`

### Observed geometry

- workflow run: `34765070227`
- artifact: `chapter1-observed-response-geometry-34765070227`
- artifact ID: `10320630085`
- digest: `sha256:5ac34231b0cfa4acc1d536ab798a574a3cbcdd9423b6f0f2995195aab78a2f23`

## Cross-scope result

The 12 cells are three pre-existing broad measurement-domain contrasts (`plain_colour`, `generalized_form`, `self_compatibility`) × two biogeographic contexts (`northern_midlatitude`, `tropical`) × two floristic strata (`all_native`, `native_nonendemic`).

Every cell failed the frozen requirement that nonlinear evidence pass the independently calibrated gate in both `all_analysis_eligible` and `direct_only` evidence scopes with a V2-qualified common shape.

Two scope-specific signals crossed the calibrated nonlinear gate but were not promoted:

- tropical, all-native `generalized_form`: `direct_only` selected a step with `D=16.912 > critical D=14.239`, but the all-analysis scope did not pass;
- tropical, native-nonendemic `self_compatibility`: all-analysis selected a step with `D=18.460 > critical D=13.378`, but the direct-only scope did not pass.

These remain evidence-scope-sensitive descriptive signals only.

Several step coefficients also have cluster-robust 95% intervals excluding zero even where the calibrated nonlinear gate does not pass. For example, tropical `plain_colour` shows positive best-step contrasts in both scopes, but `D` remains far below the simulation-derived critical values. This demonstrates why a significant step coefficient or best raw AICc alone is insufficient evidence for a breakpoint in this design.

## Scientific interpretation

The current global Chapter 1 data do not support reframing the floral island response around a general threshold or non-monotonicity result. The main inferential spine therefore remains:

1. biogeographic context changes the **direction/composition** of the floral and reproductive response;
2. the response is expressed at different **hierarchical/taxonomic depths** under bounded conditions;
3. source and lineage assembly remain central to interpreting the observed assemblage pattern;
4. a universal pollinator-channel mechanism is not promoted by the failed prospective N1 chain;
5. nonlinear response geometry is retained as a transparent negative audit rather than rescued by alternative breakpoint choices.

The appropriate Chapter 1 framing is therefore not "isolation crosses a universal floral threshold". It is that isolation is associated with context-dependent response vectors and assembly depth, while the data do not support one universal response geometry.

## Connection to Chapter 2 (`izu-core`)

This negative global result strengthens, rather than removes, the role of the original Izu threshold programme. `izu-core` can ask a different and more identifiable question within a focal biological system: whether floral, outcrossing and autonomous-reproduction channels show different response geometries across directly observed pollination regimes, and whether a predeclared threshold co-localizes with a functional transition.

Thus the chapter division becomes:

- **Chapter 1 (`island`)**: where/whether, response direction, biogeographic context and assembly depth; nonlinear geometry was tested but not promoted globally.
- **Chapter 2 (`izu-core`)**: how/why within a focal system, including cline-versus-threshold response geometry under direct functional and reproductive evidence.

## Claim ceiling

This result does not identify pollinator loss, effective service, relaxed selection, selfing selection, source assembly, in-situ evolution, or a physical open-water dispersal threshold from response geometry alone.
