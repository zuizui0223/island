# Chapter 1 H5 — strict within-island identity filter

## Result

The strict within-island identity-preserving test is **not identifiable with the current exact-island occurrence evidence**.

This is not a null biological result. The test fails before an outcome coefficient can be interpreted because both the frozen support gate and the design-rank gate are incomplete.

Canonical provenance:

- design contract: `chapter1_h5_within_island_identity_filter_v1`
- design frozen before focal-outcome inspection at commit `576b15ba88eb4ceea2e14786a56b13220d79b06d`
- final workflow run: **35062545145**
- artifact: **10433495032**
- artifact digest: `sha256:0499c2141c59cee05c7a94220947b3cddcdf943f98ff75163a57fc316303d76d`
- final workflow head: `3d18156ea1b910fdda16336944038a21afb9f25e`
- validation: 4 tests passed; Ruff passed.

## Why this test was added

The earlier pooled five-channel bridge was testable but unsupported. The later identity-aware stacked diagnostic retained channel identity and produced a suggestive direct-only tropical reproductive-assurance association, but it still compared different islands and did not survive the complete inferential gate.

The new test removes that remaining island-level comparison. Its unit is **island × identity-matched pollinator channel**, and it asks whether, within the same island, the floral architecture matched to a strict disrupted channel is lower than the architecture matched to retained channels.

The fixed mapping is:

- `bombus` -> `large_bee_like`
- `lepidoptera` -> `butterfly_like`
- `flower_visiting_birds` -> `bird_like`

Non-Bombus bees and Diptera are excluded because no equally specific predeclared named floral template was frozen for them.

The model uses island fixed effects and channel × biogeographic-context fixed effects. Therefore all island-level covariates are absorbed, while static differences among channel identities and syndrome templates cannot generate the disruption coefficient.

## Frozen support gate

An island enters only if it has at least two strict evaluable matched channels and contains both a retained and a disrupted channel. A context needs at least 10 such mixed islands.

The exact same support is obtained for all-analysis and Direct evidence:

| context | mixed islands | channel rows | clusters | frozen minimum | support gate |
|---|---:|---:|---:|---:|---|
| northern mid-latitude | **5** | 13 | 5 | 10 | fail |
| tropical | **17** | 40 | 9 | 10 | pass |

Total mixed islands = **22**.

Thus the North–Tropical joint test already fails the frozen overlap gate because northern support is only half the predeclared minimum.

## The deeper problem: design-rank non-identifiability

Lowering the minimum from 10 to 5 would not solve the problem.

After island fixed effects and channel × context fixed effects are entered, the two context-specific disruption columns should add two independent dimensions to the design matrix. They add only one:

- nuisance rank = **26**
- full design rank = **27**
- target rank increment = **1**
- required increment = **2**

Context-specific estimability is therefore:

- northern mid-latitude: **not estimable**
- tropical: **estimable in principle**

The implementation now checks this explicitly and suppresses pseudo-inverse coefficients whenever the target does not add the required design rank.

### Why North is structurally non-estimable

Among the five northern mixed islands:

- Bombus disrupted: **5/5**; Bombus retained: 0
- Lepidoptera retained: **5/5**; Lepidoptera disrupted: 0
- flower-visiting birds retained: 3; disrupted: 0

So within the northern informative set, `disrupted` is completely tied to pollinator-channel identity. Once channel identity is controlled, no independent northern disruption contrast remains.

In the tropics:

- Bombus disrupted: 16; retained: 0
- Lepidoptera retained: 17; disrupted: 0
- flower-visiting birds retained: 6; disrupted: 1

That single bird-disrupted island breaks exact channel-state confounding and supplies one independent target dimension, which is why the tropical target is structurally estimable while the two-context joint target is not.

## Guardrails added during implementation

Two failures were found and repaired with explicit regression tests.

1. **Conflicting island context/block assignments.** The first implementation could silently keep the first duplicated covariate row. RED run **35061682393** reproduced the failure. The code now rejects any island with conflicting `analysis_regime` or `spatial_block` assignments before model construction.
2. **Pseudo-inverse masking non-identifiability.** The first implementation could numerically return coefficients even when disruption was in the nuisance design space. RED run **35062222553** established this failure. The code now requires the context-specific disruption target to increase design rank before any coefficient or p-value is emitted.

The final green run **35062545145** passes all four tests and the lint gate.

## Falsification rotations

Two fixed syndrome-label rotations were frozen with the design. They have the same support geometry and the same target-rank increment of one. Therefore neither rotation is testable under the complete two-context gate either.

This is expected: the present limitation is in the independent channel-state geometry, not the focal syndrome outcome.

## What this changes in H5

The evidence hierarchy is now sharper:

1. plant-side reproductive-assurance and pollination-associated floral responses are supported as partially separable response components;
2. named syndrome templates do not identify visitor identity because most architecture is shared;
3. global GloBI interaction breadth does not identify one robust universal mechanism;
4. individual strict channels lack balanced cross-context support;
5. pooling channels creates overlap but does not support a simple total-channel-attrition mechanism;
6. identity-aware stacking reveals a suggestive tropical pattern but fails FDR/sensitivity promotion;
7. **a stricter same-island, identity-preserving comparison cannot adjudicate the question because the exact-island channel states are structurally non-identifiable in the northern context.**

The seventh result is useful because it closes a tempting analytical escape route. Reweighting, adding island covariates, or relaxing the frozen `n=10` threshold cannot recover a northern identity-specific disruption effect from the current exact-island state geometry.

## Claim ceiling

The defensible conclusion is:

> **Current exact-island occurrence states do not contain enough identity-crossing variation to identify a strict within-island pollinator-channel effect across the northern-midlatitude and tropical contexts. The northern target is structurally confounded with channel identity, so failure to estimate it is an identifiability result rather than evidence for absence of pollinator effects.**

This analysis does not identify effective pollination service, historical pollinator loss, abundance decline, visitation decline, pollen delivery, or causal selection on floral traits.

The remaining biological identification target therefore still requires independent variation in pollinator identity/service: ideally abundance, visitation or effective pollen delivery measured on islands where the same functional channel can be observed in both retained/high-service and disrupted/low-service states within the same biogeographic context.
