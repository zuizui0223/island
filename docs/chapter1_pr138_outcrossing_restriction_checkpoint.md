# PR138 SI / outcrossing-restricted attraction checkpoint

Run: `32962391754`

This is the strongest current decomposition of the two candidate explanations for the Palearctic floral response using the frozen PR138 data.

The existing literature-predeclared `large_bee_like` and `generalized_accessible` syndrome definitions are unchanged. `attraction_shift = (-large_bee_like + generalized_accessible) / 2` is recomputed only within reproductively restricted species subsets, then the external GIFT mainland source-pool expectation is recomputed within that same subset before fitting the distance model.

Restrictions are exact and fail-closed:

- `si_only`: `self_incompatibility == SI`
- `predominantly_outcrossing`: `mating_system == predominantly_outcrossing`
- `si_and_predominantly_outcrossing`: both exact states

Ambiguous reproductive states are excluded rather than coerced.

## Coverage

- SI species: 3,879 reproductively classified species; 2,088 with attraction-syndrome scores; 385 islands with source-adjusted scores.
- Predominantly outcrossing species: 7,517 classified; 3,995 with attraction-syndrome scores; 362 islands with source-adjusted scores.
- SI AND predominantly outcrossing: 730 classified; 334 with attraction-syndrome scores; 272 islands with source-adjusted scores.

## Palearctic result

Across four prespecified external-mainland source definitions and both primary floristic strata:

### SI-only

- all-native: positive in 4/4 and FDR-supported in 4/4; slope `+0.0714 to +0.0788`, q `~3.2e-8 to 7.1e-6`.
- native-nonendemic: positive in 4/4 and FDR-supported in 4/4; slope `+0.0643 to +0.0687`, q `~1.1e-5 to 8.1e-5`.

### Predominantly outcrossing

- all-native: positive and supported in 4/4; slope `+0.0418 to +0.0536`, q `~0.0054 to 0.0250`.
- native-nonendemic: positive and supported in 4/4; slope `+0.0488 to +0.0608`, q `~0.0009 to 0.0056`.

### SI AND predominantly outcrossing

- all-native: positive in 4/4 and supported in 3/4; slope `+0.0517 to +0.0585`; the `geo50_climate10` mode is borderline (`q≈0.070`).
- native-nonendemic: positive and supported in 4/4; slope `+0.0550 to +0.0639`, q `~0.0062 to 0.0280`.

Thus the Palearctic isolation-associated attraction/accessibility response persists even among species that cannot self via the measured SI axis, and also among predominantly outcrossing taxa. Measured selfing syndrome is therefore not a necessary explanation for this floral response.

## Broad northern mid-latitudes

All restricted-sample slopes are positive, but none is FDR-supported at the broad northern-midlatitude level. This is consistent with the accumulating evidence that the robust confirmatory floral branch is specifically Palearctic rather than universal across all northern mid-latitudes.

## Tropical counter-pattern

Among SI species, tropical source-adjusted attraction_shift is negative and FDR-supported in all 8 primary-stratum scenarios:

- all-native: `-0.1485 to -0.1358`
- native-nonendemic: `-0.1355 to -0.1173`

So even where reproductive selfing is unavailable by the measured SI axis, tropical isolation does not converge on the Palearctic generalization direction. This further rejects a globally fixed selfing-syndrome coupling.

## Mechanistic interpretation ceiling

This restriction directly rules out measured selfing capacity as a necessary condition for the Palearctic attraction/accessibility reorganization. Together with the continuous selfing-core conditional and interaction analyses, the strongest supported decomposition is:

1. Palearctic isolation is associated with a floral attraction/accessibility shift that persists after source-pool subtraction.
2. The shift persists after conditioning on reproductive assurance.
3. It is not detectably stronger where reproductive assurance is higher.
4. It persists in self-incompatible and predominantly outcrossing taxa.

Therefore a simple `reproductive assurance -> selfing syndrome -> floral simplification` pathway is insufficient for the Palearctic pattern.

This still does not identify relaxed attraction selection, Bombus loss, or any specific pollinator guild causally. Those remain mechanistic interpretations requiring independent pollinator evidence.
