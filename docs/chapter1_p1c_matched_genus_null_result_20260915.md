# P1c matched-complexity genus null — final result — 2026-09-15

## Question

Does the large attenuation of the Palearctic floral-island response after true genus adjustment reflect biological genus structure, or can a generic fine partition with the same within-family grouping complexity absorb a comparable amount of signal?

P1c is a post-baseline specificity test. It does not reclassify historical H3 and cannot identify a causal assembly process.

## Frozen design

The canonical P1 contract was frozen on branch `ch1/nee-island-first-design` before the full null was inspected.

Primary test:

- evidence scope: **direct-only**;
- context: **Palearctic**;
- profiles: four frozen source modes × two floristic strata = **8**;
- randomization unit: accepted species;
- randomization constrained within family;
- each family retains exactly the true number of genus groups and the exact multiset of genus group sizes;
- species scores, island memberships, source assignments, support thresholds, covariates and regression contract are unchanged;
- primary statistic: median conditional family-to-genus attenuation across the eight profiles;
- permutations: **2,000**;
- seed: **20260915**;
- one-sided test: true-genus statistic greater than or equal to the matched pseudo-genus null;
- pass threshold: `p <= 0.05`.

Qualification was outcome-closed with respect to pseudo-genus attenuation and passed before the full null was opened:

- run `34935831652`;
- artifact `10382904745`;
- digest `sha256:6972bb32da46ebf10c2124679e1c78b02fe5b2967280b8101fa99bb75764abbb`.

## Permutation execution and aggregate repair

The full permutation run was `34936193944` at head `5784beecbcfbfd68125bb0dade90cc1b30f4a6c9`.

All **40/40** shards completed successfully and generated the complete predeclared schedule of **2,000** permutations. The original aggregate job failed only because it required exact bitwise equality among independently recomputed copies of the same observed floating-point statistic.

No permutation was regenerated. An aggregate-only repair replaced exact float equality with a fail-closed numerical-equivalence criterion `max - min <= 1e-12`. The actual observed-statistic spread among shards was only `3.33e-16`.

Final aggregate:

- workflow run `34941827774`;
- artifact ID `10385820775`;
- artifact `chapter1-p1c-matched-genus-null-final-34941827774`;
- digest `sha256:861d18fe87b9190619f925a2446be5fd4d460b818825578930883257c0a6ed16`.

## Result

All **2,000/2,000** permutations were valid.

Observed true-genus primary statistic:

- median conditional family-to-genus attenuation = **0.7206615**.

Matched pseudo-genus null:

- 2.5% quantile = `-0.30056`;
- 25% quantile = `0.04496`;
- median = **0.22434**;
- 75% quantile = `0.40363`;
- 97.5% quantile = `0.72954`.

There were **57/2,000** matched partitions with a primary statistic at least as large as the observed value. With the frozen finite-randomization correction,

- one-sided randomization `p = (1 + 57) / (1 + 2000) = 0.02899`.

Frozen classification:

> **`true_genus_exceeds_matched_complexity_null`**

## Interpretation

This closes an important alternative explanation for H3. The large genus attenuation is not reproduced as readily by arbitrary fine grouping once family membership, genus-group count and the entire within-family genus-size distribution are preserved. True genus boundaries therefore contain biological/taxonomic structure relevant to the Palearctic response beyond generic partition complexity.

P1c should be read together with P1d. The paired spatial-block bootstrap showed that the exact *incremental magnitude* of the family-to-genus drop is imprecisely localized, even though total genus attenuation remains large. The combined inference is therefore:

> **The Palearctic floral-island response contains genus-specific taxonomic structure that exceeds matched arbitrary fine grouping, but the exact size of the incremental family-to-genus attenuation is not precisely estimated across spatial blocks.**

This is stronger than `fine taxonomic composition sensitivity`, but narrower than `a precisely measured assembly depth at one exact taxonomic boundary`.

## Claim ceiling

P1c does not show:

- that dispersal alone causes genus sorting;
- that within-lineage evolution is absent;
- that genus attenuation is causal mediation;
- that pollinator loss or service disruption is the upstream mechanism;
- that the family-to-genus incremental attenuation magnitude has a narrow uncertainty interval.

It supports **genus-specific taxonomic localization against a matched-complexity alternative**.
