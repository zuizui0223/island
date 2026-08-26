# PR138 regional lineage decomposition checkpoint

Status: **secondary assembly decomposition / claim-boundary checkpoint; not a causal pollination model**.

The exact-SI analysis showed that the Palearctic attraction/accessibility slope does not exceed a genus-composition-preserving source-pool null. This checkpoint asks the broader Chapter 1 question: how much of the two primary regional response axes is represented by genus composition, and what remains after subtracting the deterministic genus-mean expectation?

## Frozen inputs and reproduction gate

The decomposition uses only frozen PR138 inputs:

- observed source-adjusted island scores from the source-adjusted pathway artifact;
- species-level pathway syndrome scores, including `selfing_core`;
- frozen island floristic status;
- frozen GIFT native mainland floras;
- the four frozen outcome-blind source assignments;
- frozen geography / area / climate / spatial blocks;
- frozen observation-selection weights from joint source+selection run `32958663810`.

For each axis, every scored species is replaced by the mean score of scored species in its genus. Island membership, mainland membership, score missingness and source assignment remain fixed. Mainland source expectations are then recomputed from these genus-mean scores.

The observed branch is never reconstructed from genus scores. It is the frozen observed source-adjusted response.

Independent frozen-input reproduction gives:

- **192/192 focal observed context slopes reproduced**, maximum absolute coefficient difference `2.9e-16`;
- **48/48 focal two-axis joint regional contrasts reproduced**, maximum absolute Wald-statistic difference `1.7e-13`.

Thus this decomposition is on the existing PR138 estimand rather than a replacement analysis population.

## 1. Northern-midlatitude versus Tropical response geometry

Across 4 source definitions × 3 selection modes × 2 primary floristic strata = **24 scenarios**:

- observed two-axis regional vector: supported in the existing joint analysis in **24/24** scenarios;
- genus-expected two-axis vector: nominal `p < 0.05` in **24/24** scenarios;
- observed-minus-genus residual vector: FDR-supported in **24/24** scenarios.

The key point is that the two axes decompose very differently.

### Accessibility generalization

The interaction is Tropical minus Northern-midlatitude; therefore a negative value means the tropical isolation slope is more specialized-access / less generalized than the northern slope.

Across all 24 scenarios:

- observed interaction: `-0.2896` to `-0.1694`;
- genus-expected interaction: `-0.2709` to `-0.1555`;
- residual interaction: only `-0.01924` to `-0.01305`.

The median genus-expected / observed interaction ratio is **0.928** (range `0.905–0.935`).

Therefore roughly **91–94% of the measured North–Tropics accessibility contrast is represented by genus composition under this deterministic decomposition**. A small residual remains, but the floral-accessibility contrast is overwhelmingly lineage-associated.

### Reproductive assurance

For the same 24 scenarios:

- observed interaction: `+0.1068` to `+0.1514`;
- genus-expected interaction: `+0.0121` to `+0.0640`;
- residual interaction: `+0.0777` to `+0.1008`.

The median genus-expected / observed interaction ratio is only **0.254**.

Thus the residual two-axis regional vector is not evidence that the large floral-accessibility contrast survives lineage control. It is driven mainly by **reproductive-assurance heterogeneity that is not represented by genus means**.

## 2. Within-region decomposition

### Northern mid-latitudes

Accessibility generalization:

- observed slope: `+0.0223` to `+0.0516`;
- genus expected: `+0.0168` to `+0.0455`;
- residual: `+0.0051` to `+0.0078`.

Reproductive assurance:

- observed: `-0.0028` to `+0.0257`;
- genus expected: `-0.0003` to `+0.0215`;
- residual: `-0.0025` to `+0.0080`.

The broad northern branch therefore remains weaker and heterogeneous, as already established by the observation-selection analysis.

### Tropical

Accessibility generalization:

- observed: `-0.1361` to `-0.0698`;
- genus expected: `-0.1274` to `-0.0661`;
- residual: `-0.0094` to `-0.0037`.

The tropical specialized-access direction is therefore also largely represented by genus composition.

Reproductive assurance:

- observed: `+0.0746` to `+0.1232`;
- genus expected: `+0.0157` to `+0.0506`;
- residual: `+0.0551` to `+0.0780`.

All 24 tropical reproductive residual fits have nominal `p < 0.05`. The tropical reproductive-assurance increase is therefore the clearest component that remains beyond the deterministic genus-mean expectation.

## 3. Palearctic versus Neotropical decomposition

Across the 24 source × selection × stratum scenarios:

- observed realm vector has nominal `p < 0.05` in `19/24`;
- genus-expected vector has nominal `p < 0.05` in `16/24`;
- residual vector is FDR-supported in **24/24**.

Again, the residual is not a surviving large floral generalization effect.

### Palearctic within-realm slopes

Accessibility generalization:

- observed: `+0.1083` to `+0.1344`;
- genus expected: `+0.1008` to `+0.1299`;
- residual: only `+0.0038` to `+0.0079`.

The genus-expected / observed ratio has median **0.954**. The all-species Palearctic accessibility signal therefore has a small positive within-genus residual, but about 95% of its coefficient is represented by genus composition. This does not contradict the stricter exact-SI null, where no positive beyond-genus residual was recovered.

Reproductive assurance:

- observed: `+0.0367` to `+0.0693`;
- genus expected: `+0.0559` to `+0.0678`;
- residual: `-0.0199` to `+0.0015`;
- none of the 24 Palearctic reproductive residual fits has nominal `p < 0.05`.

Thus the Palearctic reproductive-assurance slope itself is also compatible with lineage composition rather than a positive beyond-genus residual.

### Neotropical within-realm slopes

Accessibility residuals are small and positive (`+0.0146` to `+0.0213`). Reproductive-assurance residuals are much larger (`+0.0439` to `+0.0734`), with nominal support in `19/24` fits.

The Palearctic–Neotropical residual vector is therefore primarily a difference in the residual reproductive component, with a smaller accessibility correction that partly counteracts the lineage-expected realm contrast.

## Scientific interpretation

The regional decomposition changes the simplest reading of the Chapter 1 pattern.

The data do **not** support treating the two-axis regional response as one mechanistic syndrome whose components all survive lineage control equally.

Instead, the current decomposition is:

1. **Floral accessibility geography is predominantly lineage-associated.**
   - Palearctic generalization and tropical specialized-access responses are largely reproduced by genus composition.
   - exact-SI Palearctic attraction/accessibility does not exceed the genus-fixed source-pool null.
   - all-species Palearctic accessibility retains only a small positive residual.

2. **Reproductive assurance has a different decomposition.**
   - the Palearctic positive reproductive slope is also largely represented by genus composition;
   - tropical / Neotropical reproductive-assurance responses retain a large positive residual beyond genus means;
   - this residual reproductive component drives much of the remaining regional vector heterogeneity after genus subtraction.

The strongest current Chapter 1 ordering is therefore more precise than a single `isolation -> syndrome` chain:

`isolation -> region-specific lineage assembly -> much of floral accessibility composition`

plus

`region-specific residual reproductive filtering / species sorting`, especially in tropical and Neotropical floras.

“Residual reproductive filtering” here is statistical shorthand for observed-minus-genus expectation. It does **not** identify Baker's law, pollinator limitation, in-situ evolution, or any other causal process.

## Claim ceiling

Supported:

- biogeographic contingency is real at the assemblage level;
- most of the floral-accessibility regional contrast is represented by genus turnover;
- the remaining two-axis regional heterogeneity after genus subtraction is mainly reproductive rather than floral-accessibility;
- Palearctic and tropical assemblages are therefore not simply different magnitudes of one universal island syndrome.

Not supported:

- a causal pollination-channel explanation for lineage turnover;
- a large lineage-independent Palearctic floral generalization effect;
- equating the residual reproductive component with selfing evolution or Baker's law;
- historical source ancestry or within-island evolutionary change.

A dedicated Actions workflow is staged to materialize this deterministic decomposition from the same frozen upstream artifacts and to gate on exact reproduction of the existing observed coefficients. Connector-authored pushes have not auto-triggered that new workflow, so this checkpoint remains **provisional until workflow materialization**, although the frozen-input independent calculation reproduces the existing observed estimand to numerical precision.
