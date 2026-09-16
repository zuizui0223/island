# Chapter 1 H5 — multichannel pollinator-disruption two-route bridge

Status: **completed; simple global two-route mechanism not supported.**

## Why this test was added

The earlier upstream bridge privileged Bombus. The broader biological hypothesis is not specifically "Bombus loss", but whether loss or disruption of pollination opportunity can affect island floras through two partially separable plant routes:

1. **reproductive assurance** — reduced pollination reliability favors self-compatibility / selfing / autonomous assurance;
2. **pollinator-facing floral architecture** — reduced or altered pollinator availability changes the value of floral accessibility, attraction and mechanical matching independently of measured selfing.

To test that broader hypothesis, this bridge pools five independently observed functional channels: Bombus, non-Bombus bees, Lepidoptera, flower-visiting birds and Diptera.

## Exposure definition

The exposure is deliberately conservative. A channel contributes only when it was source-available and its island observation state is evaluable under the frozen exact-island Search policy.

- `disrupted` = source-available + `adequate_non_detection`;
- `retained` = `detected` **and** the same background-effort threshold needed to call a non-detection is satisfied;
- insufficient-effort and unresolved states remain missing.

The symmetric effort requirement is important because detections can otherwise be confirmed at much lower effort than non-detections.

The primary island exposure is `any_channel_disrupted`, requiring at least **two effort-qualified evaluable channels**. Models adjust for which channel set is evaluable as well as distance, area, climate PC1-4 and spatial block.

This is a source-to-island **functional-channel disruption** metric. It is not temporal decline, abundance decline, historical extinction, visitation loss or effective pollination-service loss.

## Primary support

Direct High/Medium all-observed data, minimum two evaluable channels:

| context | islands | any disrupted | no documented disruption |
|---|---:|---:|---:|
| northern mid-latitude | 303 | 15 | 288 |
| tropical | 99 | 21 | 78 |

Unlike every single-channel test, the pooled exposure passes the frozen reference overlap requirement of at least 10 islands in both states in both primary contexts. The broader pollinator-disruption hypothesis is therefore testable rather than blocked by the Bombus-specific overlap failure.

## Route A — reproductive assurance

Prediction: documented functional-channel disruption should increase `selfing_core`.

Primary direct results:

- northern mid-latitude: beta = **-0.0113**, SE = 0.0338, p = **0.738**;
- tropical: beta = **+0.0588**, SE = 0.0948, p = **0.536**.

The predicted positive direction is therefore not reproduced across both contexts and neither coefficient is supported.

There is a threshold-sensitive tropical hint: when the minimum evaluable-channel requirement is relaxed to one channel, tropical `selfing_core` is +0.1895 (p=0.0257), and at minimum three channels it is +0.1799 (p=0.00213). This does not survive as the frozen primary min-two result and is not accompanied by the second pathway, so it is retained as suggestive sensitivity only.

## Route B — floral accessibility after reproductive assurance

Prediction: if pollinator disruption directly changes the value of pollinator-facing architecture rather than acting only through selfing, `generalized_accessible` should increase with disruption even after conditioning on `selfing_core`.

Primary direct results:

- northern mid-latitude: beta = **-0.00496**, p = **0.834**;
- tropical: beta = **-0.00863**, p = **0.897**.

Neither context supports the predicted positive accessibility shift. The secondary common named-architecture summary and the legacy attraction-shift diagnostic are also unsupported in the primary min-two analysis.

## Joint decision

The four primary route × context tests have FDR q approximately **0.897**. The overlap gate passes, but the directional two-route requirement and FDR requirement both fail.

Therefore:

> **Pooling Bombus, other bees, butterflies/moths, flower-visiting birds and flies solves the single-channel overlap problem, but it still does not support a simple global mechanism in which source-to-island pollinator-channel disruption drives both reproductive assurance and an independent floral-accessibility shift.**

This is stronger than the Bombus-only null because the broader exposure is actually estimable in both primary contexts.

## What this does and does not mean

It does **not** show that pollinators are unimportant. Several biologically different possibilities remain compatible with the data:

- effective pollination service can decline without complete loss of a functional channel;
- abundance and visitation frequency may matter more than presence/absence;
- different channels may compensate for one another;
- floral architecture may respond to channel identity/turnover rather than total channel attrition;
- native lineage assembly may absorb some interaction effects before a residual plant response is measured;
- local within-lineage selection is not identified by global cross-sectional occurrence states.

The current global data therefore support the two plant response components, but do not identify **overall pollinator-channel loss** as their common upstream cause.

## Provenance

- contract: `chapter1_h5_multichannel_two_route_v1`
- workflow run: **35051029608**
- artifact: **10429046082**
- digest: `sha256:ec56ad97899bc90774c05fd08c7384bebc853c703a4fed4746f2741aa2b0f065`
- result lock: `config/chapter1_h5_multichannel_two_route_v1_result_lock.json`
