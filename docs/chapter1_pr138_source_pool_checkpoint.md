# PR138 source-pool and pathway checkpoint

## Status

This checkpoint is a sensitivity layer. It does not replace the observed assemblage primary estimand and it does not identify historical ancestry, in-situ evolution, pollinator loss, or a pollinator mechanism.

Source-pool assignment is outcome-blind. The source universe uses suitable GIFT 3.2 `Mainland` entities with strict native Tracheophyta flora. Four source definitions were frozen before reading source-adjusted ecological outcomes:

- `geo_k5`
- `geo_k10`
- `geo_k20`
- `geo50_climate10`

Run `32954909953` completed successfully. The GIFT source universe contained 748 mainland entities with native flora, 746 with coordinates and complete frozen climate PCs, and 1,666,186 unique native entity-species rows. The same frozen PR138 species-level syndrome scores were used for mainland and island scoring.

## Source-pool-standardized syndrome result

The Northern mid-latitude all-native classic projection remained positive and FDR-supported under all four source definitions:

- `geo_k5`: +0.0429, q=0.00268
- `geo_k10`: +0.0447, q=0.00320
- `geo_k20`: +0.0435, q=0.00162
- `geo50_climate10`: +0.0421, q=0.00177

The Northern–Tropical direct syndrome-vector difference also remained supported under all four source definitions.

Within Northern mid-latitude all-native flora, the source-adjusted attraction components remained stable: `large_bee_like` decreased and `generalized_accessible` increased. The broad `selfing_syndrome` component attenuated relative to the unadjusted analysis.

Palearctic source-adjusted classic projections were supported under all four source definitions. The Palearctic–Neotropical direct vector difference was source-definition-sensitive in all-native flora but supported under all four definitions in native non-endemics, consistent with a floristic/endemic-composition contribution layered on top of the source-adjusted contrast.

## Source-adjusted attraction versus strict selfing core

Run `32956197702` completed successfully.

### Northern mid-latitude

Across all four source definitions and both all-native and native-nonendemic strata:

- strict `selfing_core ~ distance` was not supported;
- `attraction_shift ~ distance` was positive and FDR-supported;
- the positive distance effect on attraction shift remained after conditioning on source-adjusted strict selfing core;
- the conditional selfing-core coefficient was not supported.

All-native attraction distance estimates after conditioning on selfing core were approximately +0.0475 to +0.0526 across source definitions. Strict selfing-core distance estimates were only about +0.0140 to +0.0188 and unsupported.

This is evidence that the Northern attraction/accessibility reorganization is not reducible to the measured reproductive selfing core or to simple mainland source-flora composition. It is not causal evidence for Bombus loss or relaxed pollinator selection.

### Tropical

Source-adjusted attraction shift was strongly negative under all four source definitions and remained negative after conditioning on selfing core. Native-nonendemic tropical strict selfing core increased with distance under all four source definitions.

Thus tropical flora do not behave as a simple absence of reproductive assurance. Instead, reproductive-assurance change and floral-access architecture can move in different directions.

### Formal realms

Palearctic showed both components independently after source adjustment:

- strict selfing-core distance slope: positive and supported under all four source definitions;
- attraction-shift distance slope: positive and supported under all four source definitions;
- attraction-shift distance slope remained supported after conditioning on selfing core.

Nearctic remained below the 30-island pilot gate for the strict pathway decomposition and was not promoted to a replication claim.

## Current interpretation ceiling

The strongest supported statement is now:

> Mainland distance is associated with a biogeographically contingent reorganization of floral/pollination-syndrome composition that is not explained away by simple external mainland source-flora composition. In Northern mid-latitude islands, the most robust component is a shift away from large-bee-compatible / restricted-access architecture toward generalized accessibility, and this shift is not reducible to the measured strict reproductive selfing core. Within the Palearctic, reproductive assurance also rises independently.

Remaining causal boundaries:

- no direct identification of Bombus loss or relaxed attraction selection;
- Nearctic strict-pathway replication remains below pilot support;
- source-pool sensitivity is a GIFT-based proxy, not a reconstruction of true historical colonization probabilities;
- assemblage turnover versus in-situ evolution remain distinct unresolved processes.
