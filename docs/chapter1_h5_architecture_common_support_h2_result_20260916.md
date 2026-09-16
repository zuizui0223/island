# Chapter 1 H5 — architecture common-support H2 adjudication

## Result

The six-atomic North–Tropical H2 contrast is **not robust across evidence scopes on the strict named-architecture support**.

Canonical provenance:

- design frozen before the common-support atomic outcome fit at commit `45a913c574ad251629a082d6c7d8cf17c2bb2594`;
- RED run: **35083140904**;
- final GREEN run: **35083454858**;
- artifact: **10440734711**;
- digest: `sha256:60f76fd6a09195943730b64bfea2a6188205444851f8c73da6632cdb7767b64b`;
- validation: **4 tests passed** and Ruff passed.

This is a post-hoc adjudication of the orthogonalized named-architecture analysis. It is not a replacement for the full-data H2 result.

## Why this gate was necessary

The orthogonalized named-architecture test uses a much narrower support than the primary six-atomic analysis. A species must have complete scores for all three sampled templates (`large_bee_like`, `butterfly_like`, `bird_like`), and an island must contain at least 50 such complete-case species.

A null result for shared or template-specific architecture can only be interpreted as informative about H2 representation if the original six-atomic H2 remains on that same support. Otherwise the support restriction itself may have removed or destabilized the H2 contrast before the architecture decomposition is considered.

The gate therefore re-fits exactly six atomic responses on the same complete-template species and strict-support islands:

- `generalized_form`;
- `actinomorphic_symmetry`;
- `shallow_open_tube`;
- `self_compatibility`;
- `selfing_mating_system`;
- `autonomous_selfing`.

The model remains grouped-binomial logit with island area and climate PC1–4 controls and spatial-block clustered covariance. The formal comparison remains `northern_midlatitude` versus `tropical`.

## All-analysis evidence

The complete-template pool contains **2,072 species**. The architecture support contains 155 North/Tropical islands before atomic-response complete-case loss (82 northern mid-latitude, 73 tropical). The six-atomic fit uses **154 islands / 57 spatial blocks** and retains all six outcomes.

The joint North–Tropical distance-response difference is:

- chi-square = **12.999**;
- df = **6**;
- p = **0.04305**.

Thus the common-support H2 remains nominally supported in the all-analysis evidence scope.

Descriptively, the largest tropical-minus-northern slope differences are `autonomous_selfing` (+0.253, p=.0051) and `generalized_form` (+0.225, p=.0065). These component p-values are descriptive here; the frozen gate is the six-dimensional joint test.

## Direct evidence

The Direct complete-template pool contains **1,453 species**. The strict architecture support contains 145 North/Tropical islands before atomic-response complete-case loss (75 northern mid-latitude, 70 tropical). The six-atomic fit uses **144 islands / 54 blocks** and again retains all six outcomes.

The joint North–Tropical distance-response difference is:

- chi-square = **7.563**;
- df = **6**;
- p = **0.27188**.

Therefore the H2 contrast does **not** survive the Direct common-support sensitivity.

`generalized_form` and `autonomous_selfing` retain positive descriptive slope differences, but the complete six-dimensional vector is not supported. A subset of favorable components cannot substitute for the frozen joint gate.

## Consequence for the orthogonalized architecture result

The preceding structural correction found:

- shared architecture North–Tropical interaction: p=.432 (all-analysis), p=.278 (Direct);
- rank-two named-template residual H2: p=.614 (all-analysis), p=.304 (Direct).

Those null results remain valid descriptions of the strict architecture-support data. What changes is their interpretation.

Because the primary six-atomic H2 itself is **not robustly retained on that same support**, the architecture null cannot be promoted to:

> “H2 is present, but it lies outside named pollination-syndrome geometry.”

The defensible statement is instead:

> **The strict complete-template support does not provide a robust bridge between the broad six-atomic H2 and the named-syndrome decomposition. The all-analysis H2 survives weakly, but the Direct H2 does not. Therefore the shared/residual architecture null is support-limited with respect to the primary H2 estimand.**

## H5 implication

This makes the role of pollination-syndrome concordance narrower, not stronger. Raw large-bee-, butterfly-, and bird-like scores remain useful biological compatibility descriptors, but they cannot carry the mechanism argument because:

1. about 86–87% of their covariance is one shared plant-architecture dimension;
2. neither the shared dimension nor the rank-two residual subspace shows a robust formal North–Tropical H2;
3. and the strict complete-template support does not robustly retain the original six-atomic H2 across evidence scopes.

This does not imply that pollinators are irrelevant. It means the present syndrome-based plant-trait layer cannot identify how much of the broad H2 response is specifically pollination-syndrome geometry.

Independent pollinator/service data remain necessary for the upstream mechanism: abundance, visitation, pollen delivery, and channel/service variation with sufficient within-context overlap.

## Claim ceiling

- The full-data H2 result is not erased by this post-hoc restricted-support analysis.
- Named pollination architecture is not excluded as a possible biological contributor.
- Pollinator identity, historical loss, abundance, visitation and effective service are not identified.
- The common-support gate is post hoc and cannot be relabelled confirmatory.
