# Chapter 1 H5 — identity-aware pollinator-disruption two-route diagnostic

Status: **completed; suggestive context-specific signal, no promoted mechanism.**

## Why this test was needed

The pooled five-channel bridge solved the Bombus-specific overlap problem but treated all functional channels as interchangeable. That is biologically restrictive: losing a bee channel, a bird channel or a lepidopteran channel need not have the same floral consequence.

This diagnostic therefore keeps the broader pollinator-reduction hypothesis but makes channel identity explicit.

- unit: island x pollinator channel;
- retained: detected **and** the same strict effort gate needed for a non-detection;
- disrupted: source-available + strict adequate non-detection;
- insufficient effort / unresolved: missing;
- island rows are reweighted so each island has equal total weight regardless of the number of evaluable channels;
- channel x biogeographic-context fixed effects absorb static differences among channel identities and regions;
- controls: source separation, island area, climate PC1-4;
- inference: spatial-block cluster-robust covariance.

This is post hoc and remains an occurrence-state analysis. It does not measure temporal abundance decline, visitation frequency, pollen delivery or effective pollination service.

## Route A — reproductive assurance across all five channels

Prediction: if reduced pollination opportunity favors reproductive assurance, strict channel disruption should be associated with higher `selfing_core` after channel identity and context are controlled.

Direct High/Medium results:

| context | channel rows | islands | disrupted rows | beta | p | FDR q |
|---|---:|---:|---:|---:|---:|---:|
| northern mid-latitude | 1,127 | 631 | 19 | +0.06735 | 0.3563 | 0.3563 |
| tropical | 435 | 277 | 34 | **+0.21983** | **0.01528** | **0.06113** |

The tropical coefficient is nominally supported and in the reproductive-assurance direction, but it narrowly misses the four-test FDR family. More importantly, it is not reproduced by the broader all-analysis evidence sensitivity:

- northern mid-latitude: beta = -0.00207, p = 0.9555;
- tropical: beta = +0.07396, p = 0.5416.

Therefore the tropical direct-only result is **suggestive, not promoted**.

## Route B — channel-matched floral architecture after selfing

A generic `generalized_accessible` response assumes all pollinator losses should move flowers in the same direction. The identity-aware alternative instead asks whether disruption of a channel is associated with lower concordance to the floral architecture predeclared for that channel.

Primary identity matches are:

- Bombus -> `large_bee_like`;
- Lepidoptera -> `butterfly_like`;
- flower-visiting birds -> `bird_like`.

Non-Bombus bees and Diptera are excluded from this identity-matched Route B because no equally specific named floral template was frozen for them.

The model conditions on `selfing_core`, so this is explicitly the floral-architecture component remaining beyond measured reproductive assurance.

Direct High/Medium results:

| context | channel rows | islands | disrupted rows | beta | p | FDR q |
|---|---:|---:|---:|---:|---:|---:|
| northern mid-latitude | 905 | 605 | 9 | +0.06469 | 0.1251 | 0.1668 |
| tropical | 356 | 261 | 29 | **-0.13820** | 0.1215 | 0.1668 |

The tropical sign is exactly the channel-matched prediction: disruption is associated with lower concordance to the corresponding pollination-associated architecture after selfing is controlled. But uncertainty is too large for support.

The all-analysis sensitivity keeps the tropical negative direction but remains unsupported:

- northern mid-latitude: +0.04687, p = 0.2976;
- tropical: -0.12066, p = 0.1624.

Thus Route B is **directionally compatible in the tropics but not identified**.

## Interpretation

The new result changes the mechanistic reading in a useful way.

The pooled five-channel null should not be interpreted as evidence that pollinator reduction has no relation to the two plant responses. Once channel identity is retained, the direct-only tropical data show a sizeable reproductive-assurance association and the matched floral response points in the expected negative direction. However, neither survives the full inferential gate.

The current evidence therefore supports this hierarchy:

```text
isolation-associated plant response
  ├─ reproductive assurance                      supported
  └─ pollination-associated floral architecture supported and partly selfing-independent

simple total channel attrition -> both routes    not supported

identity-aware channel disruption
  ├─ tropical reproductive assurance             nominal/direct-only signal
  └─ tropical matched floral architecture         expected direction, unsupported

causal pollinator mechanism                      not identified
```

The biologically plausible next target is no longer merely "how many pollinator channels were lost". It is **pollination-service limitation plus functional identity/turnover**: abundance, visitation and effective pollen delivery within particular visitor channels, with compensation among channels explicitly represented.

## Provenance

- contract: `chapter1_h5_identity_matched_two_route_v1`
- workflow run: **35053988220**
- artifact: **10430156411**
- digest: `sha256:9f845fbb1b61e862e5b9956c43d132322304f93c82ad47c538828afa7105ec30`
- result lock: `config/chapter1_h5_identity_matched_two_route_v1_result_lock.json`
