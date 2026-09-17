# Chapter 1 v13 functional bridge result — 2026-09-17

Status: **reproduced and locked; post-hoc functional triangulation only.**

Canonical execution:

- workflow run: `35141624253`
- job: `104947270671`
- head: `1df459cbadf8118475124c077c0062cc4ac5ae53`
- artifact: `10465048981`
- digest: `sha256:50b6ca693ec989690854a9e9aadb5cd8a672b5fbcf3ff84e5aabbab06a6abc8e`
- result lock: `config/chapter1_v13_functional_bridge_result_lock.json`

The workflow first passed four unit tests and Ruff, then verified the frozen Route A and Route B artifact ZIP digests before reading their exact-species `MATCHED_EFFECT_ROWS.csv.gz` outputs. No trait state, species matching rule, support threshold, publication weighting rule, or GloPL effect definition was changed.

## Question

The historical frozen Route A/B analyses asked whether reproductive assurance or generally accessible floral architecture *moderated the distance slope* of experimental pollen limitation. Those moderation families were not supported and remain not supported.

The v13 post-hoc analysis asks a different functional question:

> after controlling for geographic distance, broad biogeographic context and GloPL measurement definitions, are the already frozen trait states associated with lower current experimental pollen limitation?

The primary model is

`PL ~ trait_state + z_distance + context_intercepts + measurement_fixed_effects`,

with publication total weight fixed to one and publication-cluster robust covariance. Sensitivities are supplemental-only, no-zero-constant, within-publication, and within-publication×site comparisons.

## Reproductive-assurance pathway

### Autonomous selfing

Autonomous or delayed selfing showed the strongest functional compatibility signal.

- primary: beta = **-0.446724**, SE = 0.080253, two-sided p = **2.60e-08**, directional negative p = **1.30e-08**;
- supplemental-only: beta = **-0.303724**, directional p = **2.69e-06**;
- no-zero-constant: beta = **-0.454798**, directional p = **9.18e-09**;
- within publication: beta = **-0.306034**, two-sided p = **0.0295**, 30 mixed publications;
- within publication×site: beta = **-0.186358**, directional p = **0.0417**, 36 mixed publication-site groups.

Thus the association is not restricted to between-publication geography. The result is consistent with autonomous selfing reducing the current reproductive consequence measured by pollen-supplementation experiments. It does not show that historical island pollen limitation caused autonomous selfing to evolve.

### Self-compatibility

Self-compatible species also had a negative point estimate, but the association was imprecise:

- primary beta = -0.117727, directional p = 0.0768;
- supplemental-only beta = -0.037990;
- no-zero-constant beta = -0.121010.

Within-publication designs were rank deficient. This trait is therefore directionally compatible but not promoted.

`selfing_mating_system` remains non-evaluable under its already frozen parent support gate.

## Floral-accessibility pathway

### Generalized floral form

Generalized/open floral form had lower current pollen limitation in the primary global adjusted model:

- primary beta = **-0.184104**, two-sided p = **0.0449**, directional p = **0.0225**;
- no-zero-constant beta = -0.165775, directional p = 0.0454;
- supplemental-only beta = -0.113540, directional p = 0.1240.

Only seven mixed publications and nine mixed publication-site groups were available, below the existing clustered-design requirements. The result is therefore compatible with the proposed functional pathway, but not robust within studies.

### Actinomorphic symmetry

Actinomorphic species had a strong negative global association:

- primary beta = **-0.381194**, two-sided p = **1.18e-05**, directional p = **5.89e-06**;
- supplemental-only beta = **-0.169098**, directional p = **0.0227**;
- no-zero-constant beta = **-0.377765**, directional p = **1.22e-05**.

Within-publication and within-site point estimates remained negative but were imprecise (directional p = 0.254 and 0.314). This is global functional concordance, not within-study causal evidence.

`shallow_open_tube` remains non-evaluable under the frozen parent support gate.

## Integrated interpretation

Three of four evaluable frozen trait contrasts had negative directional primary associations at p<=0.05, but their robustness differs sharply. Autonomous selfing is the strongest bridge because it is retained under both global sensitivities and within-publication comparisons. Actinomorphy is robust across global measurement sensitivities but not within studies. Generalized floral form is primary-compatible but less robust, and self-compatibility is directionally compatible but imprecise.

These results add a functional compatibility layer to the v13 synthesis:

`isolation -> higher experimental pollen limitation`

and independently

`isolation -> more reproductive assurance + more generally accessible floral composition`,

while frozen trait states expected to reduce dependence on successful pollinator service are associated with lower current experimental pollen limitation, most clearly for autonomous selfing.

## Claim boundary

This analysis was designed after the Route A/B moderation outcomes were known. It is therefore **post-hoc functional triangulation**, not a confirmatory mediation test. It does not establish historical selection, causal trait evolution, mediation of the isolation effect, global pollinator abundance decline, or a named pollinator mechanism. The failed frozen distance-by-trait moderation families remain failed and are not reclassified.
