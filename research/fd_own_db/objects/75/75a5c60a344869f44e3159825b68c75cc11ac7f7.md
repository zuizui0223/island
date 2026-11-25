# 14 November 2026 — Same-genus endemic floral colour: geological-type comparison with real shared genera

**Original source data completed, two independently seeded numerical runs SUCCESS.** The question is whether the recorded floral phenotype of independently labelled *endemic species* differs across geological oceanic and continental islands even after keeping botanical genus identity shared. Original GIFT candidate mainland source genera are observational checklists, NOT true ancestors. All source matching uses high/medium original Direct accepted species trait observations, and neither pollinator-mediated selection nor within-species evolution is measured.

## Source and exact biological sample

- Protected original Actions [#43577486353](https://github.com/zuizui0223/island/actions/runs/43577486353) and [#43577798107](https://github.com/zuizui0223/island/actions/runs/43577798107) (20261114 and 20261115 seeds, each original data job SUCCESS), anonymous aggregate artifacts #14107201181 and #14107656150 from the final run.
- Source code [fd original taxon-deduplicated shared-genus analysis](../../src/fd/geological_endemic_shared_genus.py) and [synthetic source-/variation-guard tests](../../tests/test_geological_endemic_shared_genus.py), original [island source runner](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_endemic_shared_genus_geology.py).
- **434** original geologically typed islands (independently matched Roeble 2024/Zizka 2022), original GIFT-native five-nearest *candidate* mainland source checklists, **10** trait-specific original high/medium species_direct categories. Two independent original status columns must agree that each focal species is native AND endemic; source names are removed if they occur on the same focal island.
- **4,307** directly measured *island × endemic-species × trait* records have actual GIFT candidate same-genus source category data; these repeat species/genera and are **not** 4,307 independent endemic species/radiations. Each trait is treated independently; a species need not have all three core flower traits.
- For the main geological comparison, all exact accepted endemic plant species observed in **both** island geology types are removed, and duplicate island observations of the SAME accepted name are averaged **before** computing one mean per shared genus. Botanical genus IDs are shared across the two geological categories, and source-relative floral mismatch is calculated against each island's own distinct accepted GIFT source congeners.

## Main tropical comparison: flower primary colour

| Original source support | Oceanic | Continental |
| --- | ---: | ---: |
| Distinct endemic species with original Direct colour and original GIFT candidate source congener | **57** | **295** |
| Distinct botanical genera containing those endemic species | **28** | **127** |
| **Botanical genera shared by both geological island types** | **13** | **13** |
| Geographical islands carrying at least one of those shared genera | **17** | **6** |
| Geographical 10-degree blocks represented | **12** | **6** |
| Different first-nearest GIFT source signatures | **11** | **6** |

**9/13 shared botanical genera have different recorded endemic floral-colour category states across oceanic and continental original island floras**, whereas 4 genera carry the same directly recorded colour category. This makes colour a comparatively informative within-lineage variation candidate, unlike the symmetry endpoint.

The original source-relative, distinct-*genus*-equal **oceanic minus continental** mean colour-category discordance = **−0.0526699771**.

- Numerical seed 20261114 genus-bootstrap 95% interval **[−0.280750,+0.168061]**.
- Numerical seed 20261115 genus-bootstrap 95% interval **[−0.283890,+0.157344]**.
- Fully enumerated **2^13 = 8,192** paired-genus sign reassignments: two-sided descriptive **p=0.67578125**. This is NOT a true geological randomization test, since the sampled genera and island/source contexts are not randomly assigned across geology.
- **Leave one genus out**: effect ranges **−0.112615 to +0.026274**, and the initial negative sign is retained under **12/13** exclusions. Source-qualified original genus units share clade history and geographic backgrounds, and no phylogenetic or independent-archipelago causal claim follows.

**Inference:** No convincing directional oceanic-versus-continental floral-colour divergence effect *within the 13 actually shared genera* is established. A nonzero source-relative category difference is present but changes sign under removal of a genus and does not exceed its genus-label reallocation null.

## Morphological symmetry and other traits: structural support limits

For original tropical *floral symmetry* (different single-trait source cohort): oceanic **189 distinct endemic species / 78 genera**, continental **1,924 / 237 genera**, and **51 botanical genera are shared by island geology types**. In **zero of those 51 shared genera** do the recorded endemic plants have different symmetry category states by geological type. Consequently **symmetry has no observed focal within-genus state variation** with which to estimate a geological association; absence of a directional test is NOT proof of conserved phenotype under island adaptation. These 51 source-qualified genera differ from the 13 colour-qualified genera.

In northern-midlatitude primary colour: **14 shared genera**, of which **10** have different observed endemic colour categories, but the oceanic side comprises only **3** 10-degree geographic blocks and **2** first-GIFT source signatures; fails the pre-outcome geography/source support floor (5/5), so no regional effect coefficient is released.

Across ten primary/exploratory trait categories, **only tropical flower colour** meets all original support gates (>=5 shared genera, >=2 state-variable genera, >=5 spatial blocks and >=5 independent first-GIFT donor regions on EACH geological type). No other trait comparison receives unsupported pseudo-zero coefficient. Data for matching are evaluated before looking at flower outcomes; shared genus availability depends on observed recorded floras.

## Causal ceiling, FCP next opportunity and manuscripts

Two independent observations must be kept distinct:

1. **Assemblage sorting/colonization proxy**: original island flower trait distributions and candidate GIFT mainland source genus pools, with geological origin; no observation of propagule arrival or community establishment probability. Previously geographically disjoint geology comparison supports only **three oceanic/continental matched pairs per major regime**. Original source-repertoire effects are not universal after lineage/geographic/area/source adjustment.
2. **Potential post-colonization adaptation**: colour variation *among different endemic plant species within shared genera* is **NOT** a measured island vs mainland genetic/phenotypic divergence. A global one-colour-label-per-species Direct database cannot estimate within-species local floral response. Island-mainland population photographs, validated native endemism, credible sister species true-tip phylogeny and pollinator reproductive success data are required to demonstrate selective evolution.

A useful next DATA gate is to overlay independent original **FCP 4,475 species / 447,500 genuinely measured photos** with geological/native-endemic island GIFT floras, count original exact accepted-name coverage, and then retrieve the *locality-preserved individual photo measurements* for island and mainland populations rather than relying on a global species-wide D/W photographic summary. This should not be treated as independent genetic polymorphism or as proof of trait adaptation.

The geological origin × native/endemic × morph/colour project remains an observational Draft, not yet a confirmed pollinator-syndrome evolutionary paper.
