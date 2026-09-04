# PR138 joint source-pool + observation-selection checkpoint

Status: robustness checkpoint only. The observed-assemblage primary estimand is unchanged.

Run: `32958663810` (`pr138-joint-source-selection-32958663810`).

## Frozen construction

The analysis combines two already-predeclared sensitivity layers without redefining the response:

- source-pool adjustment: observed island syndrome score minus the expectation from GIFT native Mainland floras;
- observation-selection adjustment: stabilized inverse-probability weights fitted only to whether a branch-axis score is available in the full 8,265-island universe.

The primary branch vector is pollinator-name-free:

- `accessibility_generalization` = source-adjusted `generalized_accessible`;
- `reproductive_assurance` = source-adjusted strict `selfing_core`.

Four source definitions are kept separate (`geo_k5`, `geo_k10`, `geo_k20`, `geo50_climate10`) and three selection modes are run for each (unweighted, stabilized IPW clipped at 0.05–0.95, stabilized IPW clipped at 0.025–0.975). No source definition is selected from the ecological result.

Neither source assignment nor propensity fitting uses island branch-score values, fitted effect sizes, or p-values.

## 1. North–Tropics heterogeneity survives both corrections

The direct two-axis distance-response vector contrast between northern mid-latitude and tropical islands is supported in **24/24** source-mode × selection-mode × stratum scenarios.

All-native direct-vector FDR q ranges:

- unweighted: `1.4e-5` to `1.64e-4`;
- stabilized IPW 0.05: `2.34e-4` to `1.28e-3`;
- stabilized IPW 0.025: `3.73e-4` to `1.26e-3`.

Native-nonendemic direct-vector FDR q ranges:

- unweighted: `2e-6` to `6.6e-5`;
- stabilized IPW 0.05: `3.31e-4` to `1.94e-3`;
- stabilized IPW 0.025: `8.25e-4` to `2.98e-3`.

Therefore the robust global result is **biogeographic heterogeneity of the plant-side response vector**, not a universal northern syndrome.

## 2. Broad northern branch is not observation-selection robust

After source adjustment but before IPW, northern mid-latitude all-native accessibility generalization remains positive in all four source modes (`+0.0426` to `+0.0516`, axis q `0.0063` to `0.0123`). Strict reproductive assurance remains unsupported (`+0.0140` to `+0.0188`, q `0.36` to `0.52`).

After either stabilized-IPW specification, however, the **northern two-axis context vector is unsupported in all 16 source-mode × IPW-mode × stratum scenarios**.

All-native source-adjusted accessibility slopes remain positive but weaken under IPW:

- IPW 0.05: `+0.0356` to `+0.0405`, axis q `0.141` to `0.196`;
- IPW 0.025: `+0.0286` to `+0.0337`, axis q `0.389` to `0.481`.

Thus the broad “Northern mid-latitude accessibility branch” cannot be claimed independently after measured observation-selection correction, even though its sign remains positive.

## 3. Palearctic branch is robust to both corrections

Palearctic is supported as `accessibility_generalization + reproductive_assurance` in **24/24** source-mode × selection-mode × stratum scenarios.

All-native ranges across the four source modes:

- unweighted accessibility: `+0.1153` to `+0.1238`; reproductive assurance: `+0.0572` to `+0.0644`;
- IPW 0.05 accessibility: `+0.1278` to `+0.1344`; reproductive assurance: `+0.0395` to `+0.0436`;
- IPW 0.025 accessibility: `+0.1176` to `+0.1235`; reproductive assurance: `+0.0367` to `+0.0403`.

Native-nonendemic gives the same branch in all scenarios.

The strict reproductive-assurance slope attenuates under IPW but remains supported, while accessibility remains strongly positive. This supports a Palearctic-specific two-response pattern rather than a generic northern-latitude effect.

## 4. Tropical and Neotropical responses remain different

Tropical native-nonendemic floras commonly show the combination:

- `accessibility_generalization < 0` (specialized-access maintenance/increase), and
- `reproductive_assurance > 0`.

The reproductive-assurance slope is positive and supported across all four source modes under both IPW variants. Thus tropical divergence from the Palearctic branch is not equivalent to “no selfing”.

Neotropical accessibility is not consistently resolved, while reproductive assurance becomes more evident under IPW. The specific branch label therefore remains region- and selection-sensitive even though the direct Palearctic–Neotropical vector contrast strengthens under IPW.

## 5. Palearctic–Neotropical direct contrast

With observation-selection weighting, the Palearctic–Neotropical two-axis vector difference is supported in **16/16** source-mode × IPW-mode × stratum scenarios.

Unweighted support is weaker (3/8 scenarios), showing that the direct realm contrast should not be summarized from unweighted significance alone.

## 6. Weight diagnostics and claim boundary

- maximum realized stabilized weight: `11.9546`, below the frozen cap of 20;
- minimum effective-sample fraction: `0.3883`;
- therefore the IPW layer is informative but still a sensitivity analysis, not a new population-representative primary model.

The result supports:

> Mainland distance is associated with biogeographically contingent plant-side floral/reproductive response vectors. After both external mainland source-pool standardization and measured observation-selection weighting, the broad northern-vs-tropical vector difference persists, while an independent broad northern branch does not. The clearest positive classic branch is Palearctic, where accessibility generalization and strict reproductive assurance both rise with distance.

It does **not** identify historical source ancestry, in-situ evolution, Bombus loss, attraction-selection relaxation, realized pollinator identity, or effective pollination service.
