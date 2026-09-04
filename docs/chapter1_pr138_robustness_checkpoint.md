# PR138 syndrome robustness checkpoint

Status: **provisional / not frozen**.

This checkpoint records robustness of the pattern-first Chapter 1 result after formal pollination/selfing syndrome analysis, source-pool adjustment, selection sensitivity, reproductive restriction, spatial leverage analysis, and the exact-SI genus-fixed source-pool null. It does not identify causal pollinators.

## Headline

The strongest global result is not a universal classic island syndrome. It is **biogeographic contingency of floral/reproductive assemblage responses to mainland isolation**.

The clearest positive floral branch is Palearctic: with increasing mainland distance, assemblages shift away from `large_bee_like` architecture and toward `generalized_accessible` architecture. Tropical assemblages follow a different response vector, and the current Nearctic pilot does not reproduce the Palearctic direction.

## 1. Evidence and template robustness

Main analysis uses nonduplicated all-analysis-eligible evidence:

- species-direct evidence first;
- validated genus-consensus evidence only fills species x trait gaps;
- no species x trait is counted twice.

Direct-only independently reproduces the main regional response direction. Genus-consensus-only is weaker and does not independently create the main result. Outcome-blind syndrome-template variants, information-weight alternatives, and raw / sqrt / log1p distance forms preserve the major North-vs-Tropical contrast.

## 2. Spatial and realm structure

The original leave-one-spatial-block analysis preserves the North-vs-Tropical vector difference in **81/81** deletions. `lat12_lon20` (Aegean/eastern Mediterranean) is influential for the scalar classic/selfing component but does not create the attraction/accessibility direction.

Formal RESOLVE realm analysis sharpens the geography:

- **Palearctic** strongly reproduces the positive generalized-accessibility branch;
- **Nearctic** currently reaches only pilot support and does not reproduce the Palearctic three-axis direction (`q≈0.306`);
- southern / Oceanian counter-patterns remain pilot and heterogeneous.

Therefore the supported scope is Palearctic, not a universal northern-temperate syndrome.

## 3. Source-pool and observation-selection sensitivity

Under joint external GIFT mainland source-pool subtraction and outcome-blind observation-selection IPW (`32958663810`):

- the broad North-vs-Tropical two-axis response-vector difference is supported in **24/24** source-mode x selection-mode x primary-stratum scenarios;
- the broad northern positive branch weakens under IPW;
- the Palearctic branch remains the clearest positive branch;
- tropical native-nonendemic floras commonly combine specialized-access direction with increased reproductive assurance rather than converging on the Palearctic response.

This establishes regional heterogeneity more strongly than a universal syndrome.

## 4. Selfing decomposition

Multiple tests show that a simple `reproductive assurance -> selfing syndrome -> floral simplification` sequence is insufficient as a complete description of the Palearctic assemblage pattern:

- conditioning on strict `selfing_core` leaves the attraction/accessibility distance association;
- `distance x selfing_core` is unsupported in Palearctic fits;
- exact-SI restriction retains a positive source-adjusted Palearctic attraction shift;
- exact-SI remains positive across all 77 Palearctic leave-one-block-out fits;
- with observation-selection IPW recomputed after each deletion, exact-SI remains positive everywhere and FDR-supported after 76/77 block deletions; the sole support-loss deletion is the preidentified `lat12_lon20` Aegean block.

Thus measured selfing capacity is not necessary for the observed Palearctic floral response.

## 5. Critical lineage decomposition: the SI slope does not exceed the genus-fixed null

The next stricter test changes the mechanistic boundary.

`docs/chapter1_pr138_si_genus_source_null_checkpoint.md` applies a genus-composition-preserving null directly to the exact-SI source-adjusted response:

- `large_bee_like` and `generalized_accessible` scores are permuted only among scored SI species within the same genus;
- island membership, floristic status, score missingness, GIFT mainland membership, source-region eligibility and frozen source assignments are fixed;
- island and mainland source scores are recomputed on every draw;
- 1,000 permutations use seed `20260827`;
- the observed all-analysis exact-SI coefficients reproduce the frozen restriction artifact to numerical precision (maximum coefficient difference `2.9e-16`).

### All-analysis-eligible evidence

Across the four frozen source definitions:

- all-native observed SI slope: `+0.0714` to `+0.0788`;
- all-native genus-null mean: `+0.0752` to `+0.0841`;
- all-native observed-minus-null residual: `-0.0065` to `-0.0038`;
- native-nonendemic observed SI slope: `+0.0643` to `+0.0687`;
- native-nonendemic genus-null mean: `+0.0708` to `+0.0744`;
- native-nonendemic residual: `-0.0065` to `-0.0042`.

Empirical one-sided probabilities that the observed slope exceeds the genus null are high (`0.739–0.845`), not small. No positive beyond-genus residual is recovered.

### Strict direct-only sensitivity

Using both direct SI status and direct-only floral syndrome evidence:

- observed slope remains positive (`+0.0543` to `+0.0646`);
- genus-null mean is almost identical (`+0.0556` to `+0.0650`);
- residual slope is only `-0.0015` to `+0.0001`;
- residual FDR q is at least `0.964` and empirical one-sided q at least `0.616`.

Therefore the lineage result is not a genus-consensus gap-fill artifact.

### Consequence

The SI restriction and lineage null answer different questions:

1. **Selfing is not required** for the Palearctic attraction/accessibility assemblage slope.
2. **But genus turnover is sufficient to reproduce that SI-restricted slope.**

The current data therefore do **not** identify an interaction-mediated floral filter acting beyond measured genus composition. A pollination filter may still operate by filtering entire genera/lineages, but this dataset does not separate that possibility from biogeographic lineage assembly.

## 6. Current scientific boundary

Strongest supported statement:

> Mainland isolation is associated with a **biogeographically contingent reorganization of floral and reproductive assemblage strategies**. The clearest positive floral branch is Palearctic and is not a simple measured-selfing by-product. However, the exact-SI source-adjusted Palearctic attraction/accessibility slope is fully reproduced by a genus-composition-preserving null, so the present data do not support a residual floral filtering effect beyond lineage composition. Tropical floras follow a different assemblage branch.

Recovered from the original hypothesis:

- strong support for biogeographic contingency;
- a strong Palearctic assemblage-level attraction/accessibility shift;
- evidence that the Palearctic shift is not simply generated by measured selfing capacity;
- evidence that between-genus turnover is sufficient to reproduce the SI-restricted floral slope.

Not recovered / unresolved:

- a universal Northern-temperate classic syndrome across Palearctic and Nearctic realms;
- residual within-genus / beyond-genus attraction filtering;
- causal attribution to Bombus loss or relaxed attraction selection;
- historical source ancestry or in-situ evolutionary change;
- a confirmed tropical butterfly or southern bird replacement mechanism;
- confirmatory southern-realm counter-patterns.

## 7. Interpretation for Chapter 1

The evidence now favors the ordering:

`isolation -> biogeographic / lineage assembly -> floral and reproductive assemblage composition`

rather than a demonstrated residual sequence of:

`isolation -> pollination-channel loss -> within-lineage floral filtering`.

This is still a substantive Chapter 1 result: the study identifies where island floral/reproductive assemblage responses differ and shows that a superficially mechanistic SI-restricted pattern can be absorbed by lineage turnover when tested directly. Pollination-channel causation remains a separate question requiring independent channel evidence / Chapter 2-style mechanism identification.

## Validation state

Core synthetic CI, formal-realm CI, source-pool sensitivity, joint source+IPW sensitivity, continuous selfing interaction, exact reproductive restriction, and the existing within-Palearctic restricted-IPW workflow are green. The exhaustive block-deletion and exact-SI genus-fixed source-pool modules/tests/workflows are committed. Connector-authored commits may not auto-trigger push workflows; the genus-null frozen-input deterministic run completed locally with 1,000 permutations and exact reproduction of the existing SI estimand. PR remains **draft / unmerged / not frozen**.
