# Three-level filtering versus evolution audit — 2026-10-05

Status: post-baseline mechanistic design / evidence audit. This does not replace the final Chapter 1 H1–H4 submission inference.

## Why two levels are not enough

An island-mainland difference within the same species does not by itself identify
post-colonization evolution. Island populations may have been founded by a non-random
subset of source genotypes that already carried the island-associated state.

The relevant hierarchy is therefore:

1. **between-lineage assemblage filtering** — which species or genera from the source pool
   occur and persist on islands;
2. **within-species founder / colonist-genotype filtering** — which source genotypes
   successfully found or dominate island populations;
3. **post-colonization evolution** — derived change after island establishment.

## Current evidence ledger

The curated ledger contains **10 system/taxon rows**:

- 1 current Chapter 1 global assemblage result;
- 2 repository population-pilot systems;
- 7 external same-species island-mainland validations.

Among the nine population-level rows:

- **3** are null or explicitly non-uniform with respect to a simple island selfing shift;
- **6** show directional same-species island-mainland reproductive divergence;
- **3** retain common-environment or genetic-marker support that weakens immediate
  environmental plasticity as the sole explanation;
- **0** separate founder-genotype filtering from post-colonization evolution strongly
  enough to identify a level-3 effect.

The ledger therefore rejects both extreme models:

- the global island pattern is not adequately described as **pure assemblage filtering
  only**, because persistent within-species island-mainland reproductive divergence exists
  in several systems;
- it is also not a **universal within-lineage island evolution** response, because direct
  population comparisons include clear null/non-uniform systems and because present-day
  divergence does not distinguish founder filtering from derived evolutionary change.

## Best current model

> Island isolation first filters lineages and genera. Within established lineages,
> founder-genotype filtering and/or post-colonization evolution can add species-specific
> reproductive responses.

Chapter 1 now directly supports level 1 in the source-evaluable contexts. Population
comparisons show that levels 2+3 jointly matter in some species, but the partition between
them remains unresolved.

## Evidence examples

Directional island-mainland reproductive divergence is retained for systems including
*Lycium carolinianum*, *Waltheria ovata*, *Oeceoclades maculata*, *Limonium lobatum*,
*Nicotiana glauca* and *Campanula punctata*. The ledger explicitly does not promote these
to post-colonization evolution where founder filtering remains an alternative.

Null/non-uniform cases include the repository *Arabidopsis lyrata* and *Phormium tenax*
pilots and the published *Metrosideros excelsa* comparison.

## What would identify level 3

The highest-value designs are:

- recent island colonizations with a documented mainland source and historical/source
  genotype frequencies;
- SI/SC systems with causal S-locus alleles, allowing source pre-existence of the island
  allele to be distinguished from an island-derived mutation;
- common-garden island/source populations paired by genome-wide ancestry;
- replicated independent island colonizations showing repeated derived change after
  establishment.

Conceptually, for source genotype (g),

```
total within-species shift
= founder-genotype filtering
+ post-colonization change

Delta zbar = Cov_S(w_g, z_g^S) / wbar + E_I(z_g^I - z_g^S)
```

The first term captures which source genotypes become represented; the second captures
change after establishment.

## Claim boundary

The current evidence supports:

- assemblage filtering concentrated at genus composition in the testable Chapter 1
  source-matched contexts;
- heterogeneous within-species island-mainland reproductive divergence;
- rejection of a universal rule that island populations always become more selfing.

It does **not** yet quantify how much of the global H1 response is produced by
founder-genotype filtering versus post-colonization evolution.

## Provenance

Validated workflow:
- workflow: Three-level filter evolution audit
- run: **37258885009**
- commit: **cc5e000718c17a70a089f14798c6fe54dae517c3**
- artifact: **11323872236**

Core files:
- `data/v2/curation/within_lineage_literature_triangulation_20261005.csv`
- `docs/chapter1_three_level_filter_evolution_design_20261005.md`
- `src/island_v2/chapter1_three_level_filter_evolution.py`
- `tests/test_chapter1_three_level_filter_evolution.py`
