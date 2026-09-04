# PR138 selfing-interaction checkpoint

Run: `32961139430`

This checkpoint tests whether the source-adjusted isolation-associated attraction/accessibility response becomes stronger as measured strict reproductive assurance increases.

The model is

`attraction_shift ~ distance + selfing_core + distance:selfing_core + island area + climate PC1`

with spatial-block cluster-robust inference. `attraction_shift = (-large_bee_like + generalized_accessible) / 2`. The interaction uses continuous `selfing_core`; no selfing threshold is introduced. This is a mechanistic-concordance decomposition, not causal mediation.

## Confirmatory primary-stratum result

Across four prespecified external-mainland source-pool modes and the two primary floristic strata (`all_native`, `native_nonendemic`):

### Northern mid-latitudes

- Distance effect on attraction/accessibility shift: positive and FDR-supported in **8/8** scenarios.
- Distance estimate range: **+0.0462 to +0.0527**.
- Distance q range: **2.37e-4 to 0.00637**.
- `distance × selfing_core` interaction: **0/8** FDR-supported.
- Interaction estimate range: **-0.00375 to +0.01569**.
- Interaction q range: **0.323 to 0.985**.
- Selfing-core main effect: **0/8** nominally supported.

Thus the northern isolation-associated attraction/accessibility reorganization persists after source-pool adjustment but is not detectably stronger on islands with greater measured reproductive assurance.

### Palearctic

- Distance effect: positive and FDR-supported in **8/8** scenarios.
- Distance estimate range: **+0.1013 to +0.1152**.
- Distance q range: **7.54e-6 to 0.00148**.
- `distance × selfing_core` interaction: **0/8** FDR-supported.
- Interaction estimate range: **-0.00801 to +0.02072**.
- Interaction q range: **0.270 to 0.989**.

The strong Palearctic attraction/accessibility distance response is therefore not detectably contingent on the continuous strict selfing core in this source-adjusted analysis.

### Tropical broad region

- Distance effect: negative and FDR-supported in **8/8** scenarios.
- Distance estimate range: **-0.1475 to -0.1103**.
- Distance q range: **1.64e-5 to 8.35e-5**.
- `distance × selfing_core` interaction: negative and FDR-supported in **8/8** scenarios.
- Interaction estimate range: **-0.1052 to -0.07775**.
- Interaction q range: **4.07e-5 to 0.00746**.
- Selfing-core main effect itself is not nominally supported in these eight fits.

This is a clear counter-pattern to a globally uniform selfing-syndrome explanation: in tropical islands, greater reproductive assurance does not couple to a more northern-classic attraction/accessibility response; instead, the distance slope becomes more negative as selfing core increases.

### Neotropical realm

The Neotropical result is heterogeneous across source definitions and floristic strata. It does not justify a single realm-level interaction claim. One all-native source mode supports a positive distance effect and two all-native source modes support a positive interaction, whereas native-nonendemic fits are mostly null.

## Mechanistic interpretation ceiling

The current evidence supports the following decomposition:

1. The northern / Palearctic attraction-accessibility response is **not reducible to the measured strict reproductive-assurance syndrome**.
2. A simple `selfing syndrome -> plainer / more generalized flowers` explanation is therefore insufficient for the northern pattern.
3. This is more consistent with an attraction-selection / pollination-function shift operating independently or in parallel with reproductive assurance than with selfing syndrome as the sole pathway.
4. The tropical negative interaction further argues against a globally fixed coupling between reproductive assurance and floral simplification.

This analysis still does **not** identify relaxed attraction selection, Bombus loss, or any pollinator guild as a causal mechanism. A species-level SI / obligate-outcrossing restriction would be a stronger discriminant, but the currently frozen analysis artifacts do not contain the island-by-species native-incidence table needed to reconstruct that restriction without a new incidence/status input.
