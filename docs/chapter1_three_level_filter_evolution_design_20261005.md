# Chapter 1: three-level separation of filtering and evolution — 2026-10-05

## Why the original two-way contrast is insufficient

A difference between island and mainland populations of the same species does not by itself
identify post-colonization evolution. The island may have been founded by a non-random
subset of genotypes already present in the source population.

The biologically relevant decomposition is therefore:

1. **between-lineage assemblage filtering** — which species/genera from the mainland
   source pool are represented on islands;
2. **within-species founder / colonist-genotype filtering** — which genotypes or mating
   types within a species successfully colonize and persist on islands;
3. **post-colonization evolution** — change after establishment within the descended
   island lineage.

Chapter 1 now has direct evidence for level 1 in source-evaluable northern-midlatitude and
tropical assemblages. Population comparisons can observe levels 2+3 jointly, but most
published island-mainland contrasts do not separate them.

## What the current repository already identifies

The source-matched decomposition on native non-endemics gives the exact identity

```
raw island H1 mean
= source-species expectation
+ genus-structure enrichment
+ within-genus species-sorting enrichment
```

In northern mid-latitudes, the supported isolation response is concentrated in genus
structure and additional within-genus species sorting is not robust across all four
source definitions. In the tropics, the source-species expectation points toward H1 but
realized genus structure points in the opposite direction.

This is positive evidence for an assemblage-composition route. It does not identify
arrival, establishment, persistence or extinction.

## Population-level triangulation

The first repository pilot and selected external comparisons show that the within-species
response is not universally directional.

| System | Design | Result | What it identifies |
| --- | --- | --- | --- |
| *Arabidopsis lyrata* | 3 named Great Lakes island populations vs 15 other populations | island mean tm 0.638 vs 0.744; exact P=0.670 | no simple categorical island shift |
| *Phormium tenax* | two matched island-mainland pairs | no significant seedling-tm difference in either pair | no simple island shift |
| *Metrosideros excelsa* | two offshore-island vs three mainland populations | no significant island-mainland difference in outcrossing | no simple island shift |
| *Nicotiana glauca* | California mainland vs two recently colonized Channel Islands | island plants have greater self-pollination capacity, but no poorer pollinator service/current selection detected | favours a colonization-filter interpretation over ongoing pollinator-driven selection |
| *Lycium carolinianum* | two mainland vs two Hawaiian populations | Hawaii SC; mainland predominantly SI; S-RNase diversity reduced in Hawaii | persistent genetic differentiation, but founder filtering vs post-colonization loss of SI is not separated |
| *Limonium lobatum* | five Canary vs five Iberian populations, plants raised under common greenhouse conditions | selfed seed-set higher on islands; novel pollen-stigma combinations strongly enriched on islands | persistent within-species island divergence; plasticity reduced as explanation, founder filtering vs evolution still unresolved |
| *Campanula punctata* | mainland Honshu and Izu Island populations | outcrossing declines from mainland/Oshima to northern and southern islands; strong allozyme differentiation | pronounced island divergence, but the island group may represent an independently evolving taxon |

The corresponding curated evidence ledger is
`data/v2/curation/within_lineage_literature_triangulation_20261005.csv`.

## The key inference

The available evidence rejects both extreme models:

- **pure assemblage filtering only** is too strong, because some species show persistent
  island-mainland reproductive divergence;
- **universal within-lineage island evolution** is also too strong, because several
  direct population comparisons show no island-mainland mating-system difference and the
  direction varies strongly among species.

The best current model is hierarchical:

> island isolation first filters lineages and genera; within the lineages that establish,
> genotype filtering and/or post-colonization evolution can add species-specific responses.

The remaining unresolved question is not whether within-lineage divergence ever exists.
It does. The unresolved question is **how much of the global assemblage response is
founder filtering within species versus evolutionary change after establishment**.

## Design that can actually separate founder filtering from post-colonization evolution

A clean test needs source-population genotype information, not only present-day island and
mainland phenotypes.

For genotype or haplotype (g), let:

- (z_g^S) = reproductive phenotype in the source before/at colonization;
- (w_g) = relative representation of that source genotype among successful island
  founders/established descendants;
- (z_g^I) = present island phenotype for descendants of that genotype.

Then the within-species island-mainland shift can be written conceptually as

```
total within-species shift
= founder-genotype filtering
+ post-colonization change
```

with a Price-equation form

```
Δz̄ = Cov_S(w_g, z_g^S) / w̄ + E_I(z_g^I - z_g^S)
```

The first term is level 2; the second is level 3.

### Minimum empirical routes

**Route A — recent introductions with known source.**
Use a species with documented colonization date and source population. Genotype historical
source material and present island descendants. If the island trait state was already
enriched among founding/source genotypes, that is founder filtering; a derived island
change after establishment supports evolution.

**Route B — S-locus / mating-system causal alleles.**
For SI→SC transitions, identify causal S-locus loss-of-function alleles. If the same SC
allele exists in the source at low frequency and is enriched on islands, founder filtering
is supported. A derived island-specific loss-of-function mutation on an island
monophyletic background supports post-colonization evolution.

**Route C — common garden plus ancestry-resolved paired populations.**
Common-garden island and source populations to remove immediate environmental plasticity,
then use genomic ancestry to pair island haplotypes to their closest source haplotypes.
A persistent trait difference within ancestry-matched pairs is the minimum phenotype-level
evidence for post-colonization divergence.

**Route D — replicated independent colonizations.**
Fit the same ancestry-matched contrast across multiple independently colonized islands.
Repeated derived shifts after independent founding events are much harder to explain by a
single founder event.

## Analysis ladder

1. **Assemblage component:** retain the current source-matched genus/species decomposition.
2. **Population component:** build species-specific island-mainland contrasts only from
   direct population measurements.
3. **Founder-filter component:** test whether source genotypes carrying the island-like
   state are preferentially represented among island ancestry.
4. **Evolution component:** within ancestry-matched source-island pairs, estimate derived
   phenotype or causal-allele change after colonization.
5. **Cross-species synthesis:** meta-analyse the level-2 and level-3 effects separately;
   never pool them into one “within-lineage” coefficient.

## Claim ceiling

Current data can support:

- assemblage filtering concentrated at genus composition in the testable Chapter 1
  source-matched contexts;
- heterogeneous within-species island-mainland reproductive divergence across published
  systems;
- rejection of a universal rule that island populations always become more selfing.

Current data cannot yet quantify the global fraction of H1 produced by founder-genotype
filtering versus post-colonization evolution.

## Highest-value next empirical target

The most informative additional datasets are not arbitrary species with large mating-system
variance. They are systems with:

- documented island colonization history or source population;
- causal SI/SC loci or dense genome-wide ancestry markers;
- direct island and mainland/source reproductive measurements;
- preferably common-garden phenotyping;
- replicated island colonizations.

These data would close the only remaining mechanistic gap left after the assemblage
decomposition.
