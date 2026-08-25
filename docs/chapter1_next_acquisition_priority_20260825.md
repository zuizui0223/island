# Chapter 1 next acquisition priority — 2026-08-25

## Decision after observation-robust WHEN/WHERE run

The northern-midlatitude and tropical WHEN/WHERE result is already observation-robust under the current locked inputs. Additional acquisition must not be used to strengthen, rescue, or tune that result.

Future acquisition is justified only to expand the **testable geographic boundary** into currently underpowered regions.

## Current bottleneck is not GBIF occurrence volume alone

The frozen full-island observation layer contains:

| context | all islands | flora-recorded islands | islands with any direct Chapter 1 trait | islands entering native-status trait surface |
| --- | ---: | ---: | ---: | ---: |
| northern mid-latitude | 3,119 | 2,206 | 241 | 240 |
| tropical | 3,013 | 1,558 | 138 | 136 |
| southern extratropical | 776 | 317 | 34 | 34 |
| northern high-latitude | 1,357 | 424 | 12 | 12 |

Thus southern and northern-high-latitude regions already contain hundreds of islands with realised-flora records, but only a small subset can enter the native-status direct-trait analysis.

The next bottleneck is therefore primarily:

```text
recorded flora
→ floristic status resolved
→ direct trait evidence resolved
→ enough islands per atomic response for vector testing
```

not simply “download more GBIF occurrences.”

## Acquisition order

### Priority 1 — floristic-status resolution on already recorded floras

Resolve native / introduced / unresolved status independently of trait outcomes for species on underpowered islands.

Highest-value regions:

1. **northern high-latitude** — 424 flora-recorded islands but only 12 currently contribute native-status trait surfaces;
2. **southern extratropical** — 317 flora-recorded islands but only 34 currently contribute native-status trait surfaces.

This step is outcome-blind: islands/species are prioritized by missing status evidence and geographic testability, never by the sign or significance of a trait effect.

### Priority 2 — direct floral-form evidence

Once status is resolvable, floral form is the highest-value trait domain because:

- it contributes strongly to the confirmed northern-vs-tropical vector difference;
- direct coverage is low in underpowered regions;
- all-native median direct coverage is only about 0.141 in southern extratropical and 0.097 in northern high-latitude floras.

### Priority 3 — direct SI/SC evidence

SI/SC remains scientifically important even though SC itself does not show a clear positive distance-gradient slope in the current sample.

The goal is not to force a Baker-rule result. The goal is to make reproductive-state response vectors testable in additional regions.

### Priority 4 — colour only after the above

Colour is lower priority because direct coverage is already substantially higher:

- southern extratropical median ≈ 0.727;
- northern high-latitude median ≈ 0.815.

Colour acquisition is useful only where it increases island-level response support, not as a general database-filling exercise.

## Confirmatory target

The declared confirmatory vector threshold is >=50 islands per retained atomic response.

This threshold is an analysis-support rule, not a biological boundary. Acquisition should aim to create enough independently represented islands for multiple atomic responses to meet this support level.

Because only 34 southern and 12 northern-high-latitude islands currently enter the native-status trait surface, **trait acquisition alone cannot close the confirmatory gap** unless floristic-status coverage is expanded first.

## What not to do

Do not:

- impose a mainland-distance cutoff;
- preferentially acquire species because they carry the expected trait state;
- preferentially acquire islands because preliminary slopes are large;
- interpret an unrecorded species as absent;
- interpret an island with no trait record as trait-zero;
- launch a Bombus campaign to rescue Chapter 1;
- continue indiscriminate global trait filling when it does not change regional testability.

## Acquisition priority score

A future automated queue should be based only on outcome-blind variables such as:

```text
priority =
    underpowered_region
  × missing_status_evidence
  × missing_direct_trait_domain
  × island_contribution_to_confirmatory_support
```

Observed trait values, effect directions, p-values, and pollinator hypotheses must not enter this score.

## Source-pool extension

A separate, stronger analysis can be attempted where a defensible candidate source pool can be declared independently of island trait outcomes:

```text
P(species occurs on island | trait, candidate source pool, geography)
```

This is secondary to the present opportunistic-flora WHEN/WHERE framework and should not replace the observation-bias-aware primary analysis.
