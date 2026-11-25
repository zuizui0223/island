# Functional-trait turnover beyond species sharing: first verified source-cohort audit

**10 November 2026 | exploratory macroecological results | fd draft PR #13**

This report measures the **recorded plant-side trait composition**, not observed
flower-visitor performance, causal selection, ecological functional redundancy,
or a geologically independent island syndrome. The rights-protected original
species/trait/island rows are **not** redistributed.

## The distinct question

If two islands share similar floral access or reproductive-assurance trait
prevalence, how much of that similarity depends on **shared species**, versus
**shared genera**, rather than merely changing how many observed trait carriers
are present on each island?

The sampled fd G0 public cohort includes **7,618 overlapping island-pair
comparisons across four broad geographical regimes**, using **1,419 distinct
sampled islands**. This is *not* a set of 7,618 phylogenetically or spatially
independent contrasts. The evidence scopes **All and Direct are nested**.

Original source:
- island progressive source-recovery run `34232450884`, artifact
  `10058653212`, original ZIP SHA256
  `b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`;
- island source-scale trait run `34191508045`, artifact `10042377593`,
  ZIP SHA256 `5e9c887e093c42eb720f3abd2dc85701ce123fb7fe888c916dc4098ac898ccdd`;
- fd public G0 workflow `42769575992`, artifact `13741845313`,
  ZIP SHA256 `2efa7d2bb40adc60cc7d22b3c0c1043eb63c7c608e492fb0a422a14b3a563801`.

All original taxon-name, island, genus, and raw source data remain protected.
The released [aggregate JSON](functional_stagewise_real_aggregate.json)
contains only coarse regional/source-scope summary values and identifiers of
the archived source *artifacts*, not source taxa or islands.

## Seven source-frozen binary traits

Use the **existing Chapter 1 seven-outcome coding**, not new functional labels:
- reproductive assurance: `self_compatibility`, `selfing_mating_system`,
  `autonomous_selfing`;
- floral accessibility proxies: `generalized_form`,
  `actinomorphic_symmetry`, `shallow_open_tube`;
- secondary colour/signal proxy: `plain_colour`.

All measures are binary recorded species states; ambiguous descriptions and
unrecorded values are **missing**, never inferred as 0. Each trait/pair needs
**20 known binary carriers on each island**. The primary domain averages at
least **two eligible traits from its three frozen traits**, with equal
component weights. No post-hoc all-seven-domain syndrome score is used.

Most importantly, **shared genera are calculated from the full island species
lists, including names with missing trait data**, not just from SC-observed
species. This corrects a prior overestimate of eligible genus-excluded SC
pairs (6,524 -> **5,990 All**, and **5,400 Direct**).

## First actual eligibility gate

| Trait scope after shared-genus removal | All eligible pairs | Direct eligible pairs |
|---|---:|---:|
| Self compatibility alone | 5,990 | 5,400 |
| Reproductive assurance (>=2 of 3 traits) | 5,054 | 4,414 |
| Floral accessibility (>=2 of 3 traits) | 5,121 | 4,753 |

The original colour/form 606-island same-archipelago prediction cohort is NOT
interchangeable with this distinct 7,618-pair set; it has not yet been
cross-matched with the functional-trait evidence gates.

## Sequential equal-carrier thinning control

For each ORIGINAL pair and each independently eligible frozen trait:

1. Start with both islands' recorded taxon lists and exact trait states;
   remove taxon names shared between the two islands.
2. Independently identify genera shared between the **original whole
   floras** and remove their remaining carriers from both species-exclusive
   assemblages.
3. Retain only comparisons with >=20 binary trait carriers on both sides
   **after** these exclusions.
4. Construct **399 hypergeometric without-replacement random draws**
   directly from the already **species-exclusive** observed carrier pools
   down to the final genus-exclusive carrier counts (on each island).
5. Define `excess_gap = observed absolute prevalence difference after
   genus removal - mean random equal-carrier prevalence difference`.
6. Average within each domain using two or three traits available on the
   SAME pair. Report group-level means and medians, without pretending pair
   observations are independent hypothesis tests.

**Positive excess gap** means that removal of shared genera increases
floral-architecture prevalence dissimilarity beyond losing the same number
of randomly chosen recorded taxon carriers. That is **consistent with**
shared genera contributing to observed trait similarity. It does NOT prove
phylogenetic trait conservatism, historical colonization filtering,
functional replacement, evolutionary convergence, or geographic adaptation.
The outcome is plant-side **prevalence of categorical structural proxies**.

### Verified aggregate results

| After shared-genus exclusion, controlling sample size | All | Direct |
|---|---:|---:|
| Reproductive assurance mean excess gap | **+0.0103** (5,054 pairs) | **+0.0024** (4,414) |
| Floral accessibility mean excess gap | **+0.0182** (5,121 pairs) | **+0.0166** (4,753) |
| Reproductive share of positive pair excess | 57.9% | 50.4% |
| Accessibility share of positive pair excess | 66.6% | 66.0% |

Mean differences are **absolute prevalence-gap units**, on 0–1 scale.
For example +0.0182 corresponds to about +1.82 percentage points in the
absolute mean prevalence gap. They are descriptive means, not estimates of
selection strength or population fitness. All and Direct reuse taxa,
pairs, regions and donor species histories.

The accessibility excess is positive in all four broad regimes in the
two nested scopes, but with distinctly different magnitudes:

| Region | Accessibility All / Direct | Reproductive All / Direct |
|---|---:|---:|
| Northern high latitude | +0.0087 / +0.0210 | +0.0189 / +0.0124 |
| Northern midlatitude | +0.0161 / +0.0117 | −0.0006 / −0.0025 |
| Southern extratropical | +0.0275 / +0.0223 | +0.0055 / −0.0033 |
| Tropical | +0.0173 / +0.0145 | +0.0231 / +0.0109 |

The broad-region directions are post-hoc descriptive patterns. The northern
high-latitude sampled cohort uses only about 90 distinct islands and the
southern cohort only about 60, so **thousands of reused pair rows do not imply
thousands of independent islands**. Most sampled pairs are not near neighbors.

## Crucial interpretation boundary

These findings distinguish a reproducible, proxy-level **lineage-linked
assemblage signature** from the stronger claim that functional traits
converge independently of taxonomy. The result is, if anything,
**evidence against interpreting high nominal floral-profile sharing as
species/genus-independent ecological functional equivalence**.

The negative control adjusts **recorded binary carrier sample size only**.
It does not preserve genus frequencies under a phylogenetic null, regional
mainland source availability, real colonization history, native status,
trait-observation detection, population-specific phenotypes, or ecological
pollen-transfer service. The 7,618 pair records reuse islands and genera;
**no p-value, independently replicated origin effect or evolutionary causal
effect is justified by this analysis**.

The analogous original PR321 colour×form profile result (median nonshared-
taxon profile Jaccard about 0.212 over 2,000 descriptive sampled pairs)
measures a separate **nominal phenotype repertoire**, not proof of the
seven-trait access/reproduction function.

## Reproducible analysis and next identification gate

[Source-verifiable pair-statistic and sequential null module](../../src/fd/functional_support_null.py) accepts a rights-authorized private
pair×trait table with `scope, region, island_a, island_b, trait` and
`{all,species,genus}_{n,k}_{a,b}` counts. Example:

    python -m fd.functional_support_null \
      --private-counts /protected/fd_frozen_pair_counts.csv \
      --output work/functional_support_null_private_replay \
      --min-carriers 20 --draws 399 --seed 20261110

The CLI **never writes pair-, island-, genus-, species-, or row-level input**;
only coarse released aggregates after thresholding and a status JSON.
A second independent conservative [traitwise rarefaction module](../../src/fd/functional_rarefaction.py)
counts spatial blocks and suppresses under-supported regional aggregates.
Synthetic and public GBIF frequency verification are in the dedicated CI.
The source archives required to **reconstruct** the private pair-count input
are rights-protected and not published.

The NEXT gate is to compare **the same focal islands** with (a) nearest
within-archipelago donors and (b) **distance-matched other-archipelago donors**,
while ensuring that pairwise exclusion of shared species/genera retains enough
verified carriers. This has **not yet been run** on the seven functional
traits. If support collapses, the local-source functional hypothesis is
underpowered, not rejected. Independent spatial and source/genus nulls are
required before discussing functional maintenance or adaptation.

The active Chapter 1 still supplies the independent **GloPL experimental
pollen limitation** evidence, but this descriptive fd analysis has **not**
produced new realized pollinator-activity or pollen-limitation measurements.
