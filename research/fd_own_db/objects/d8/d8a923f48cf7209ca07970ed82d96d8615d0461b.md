# Genus-conditioned and GIFT mainland-source floral repertoire diagnostics

**Date:** 2026-11-08 JST  
**Research repository:** [zuizui0223/fd](https://github.com/zuizui0223/fd)  
**Status:** completed exploratory **private-checkpoint reanalysis**; identity-level data are NOT publicly distributed. This report and its CSVs contain only non-identifying aggregate results.

## 1. Why this matters

[Previous origin-matched analysis](../20261108_matched_flower_repertoire/README.md) found that 37 continental–oceanic island pairs (23 northern midlatitude, 14 tropical) differ in the **exact recorded colour×form repertoire** under paired label swaps (All/Direct BY q≈0.00042), even though the categorical-similarity-kernel difference was not BH-supported. Geological history could affect floral repertoire; however **which plant taxa colonized and survived**, and their regional sources, are competing explanations.

This follow-up tests (A) **genus- and geographically available taxon-conditioned profile assignment**, and (B) independently constructed GIFT **potential mainland-source floral repertoires**. The two references condition on DIFFERENT aspects of the assembly process; rejection under one and non-rejection under the other are not logically inconsistent.

## 2. Provenance: exact PR321 reconstruction

Authorized frozen `island` Actions artifacts:
- [Trait source run 34191508045](https://github.com/zuizui0223/island/actions/runs/34191508045), artifact ID 10042377593, ZIP SHA256 `5e9c887e093c42eb720f3abd2dc85701ce123fb7fe888c916dc4098ac898ccdd`.
- [Progressive flora run 34232450884](https://github.com/zuizui0223/island/actions/runs/34232450884), artifact ID 10058653212, ZIP SHA256 `b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`.

The original PR321 state-set and quality rules were reconstructed without multiplying multistate descriptions. Checks:
- 106,124 lexically eligible taxon-name units after name gate.
- 2,847 reviewed Direct cells restored.
- All **16,739 complete-name units / 496 profile types**; Direct **8,455 / 467**.
- **All 1,532 All and 1,442 Direct** eligible island `n_names` AND `n_profiles` values reproduced **exactly** from the private island×taxon presence checkpoint; no mismatches.
- The exact **37 matched physical island pairs** and the original `fd` 37-pair L2/kernel statistics were reproduced to numerical precision (All L2=0.00075638885, kernel=0.00069637962; Direct L2=0.00086057341, kernel=0.00059851945).

This validation establishes that the following conditional tests use the SAME source identities/trait definitions and match targets as the published aggregate exploration, not a new name/crosswalk approximation.

## 3. Primary conditional null — one coherent taxon profile, fixed genus and geographic availability

Randomize **whole recorded colour×form profile IDs** across different name units **within genus × observed four-regime island-flora occurrence signature**. Preserve every island×taxon-name incidence, each taxon's observed regional availability, genus distribution and trait-complete evidence support. One profile assignment per taxon is reused across every island. The 4 regions are broad observed areas, not historical mainland ancestral sources.

Under **1,999 randomizations**, compare the origin-pair mean normalized exact-profile difference and categorical similarity-kernel difference to their null distributions; family of **4 global comparisons** (two metrics × All/Direct).

| Global 37 pairs | Observed | Conditional null median | Approx. upper-tail p | BH q |
|---|---:|---:|---:|---:|
| All, exact-profile L2 | 0.00075639 | 0.00070242 | 0.1035 | 0.1847 |
| All, categorical-kernel | 0.00069638 | 0.00055197 | 0.1385 | 0.1847 |
| Direct, exact-profile L2 | 0.00086057 | 0.00078728 | 0.0675 | 0.1847 |
| Direct, categorical-kernel | 0.00059852 | 0.00049513 | 0.1860 | 0.1860 |

**No pooled-global test remains BH- or BY-supported.** The minimum BY q is approximately 0.385. This is evidence that the *previous paired label-swap association is no longer clearly exceptional under a genus- and observed geographic-availability-conditioned randomization*. It is **NOT proof that genus composition completely explains it**, nor a demonstration of equivalent phenotypes/pollination outcomes.

A second globally coherent **genus-only** null gives the same broad result: no global test is BH-supported (minimum BH q=0.209).

### Regional heterogeneity

The more restrictive genus×geographic-signature null yielded suggestive **unadjusted** residuals in the tropical 14 pairs (e.g., Direct exact-profile L2 upper-tail p=0.0070; All tropical categorical-kernel p=0.018). **None of the 8 regional comparisons remains BH-supported** (minimum BH q≈0.056) or BY-supported. A looser genus-only global null gives some tropical secondary BH/BY support, showing that regional availability constraints materially affect inference. The overall global claim must follow the **more restrictive and globally coherent** conditional reference.

## 4. Native-only sensitivity

From the SAME 37 original geography-matched pairs, retain only those with **≥20 explicitly native complete floral name units at each endpoint in BOTH All and Direct**: **28 pairs** (20 northern-midlatitude, 8 tropical). Other origin statuses are treated as unknown/excluded here, **not** declared alien.

Under **999 genus×geographic-signature randomizations**, no pooled-global test was supported:

| Native-only global 28 pairs | Upper-tail p |
|---|---:|
| All exact-profile L2 | 0.492 |
| All categorical kernel | 0.679 |
| Direct exact-profile L2 | 0.319 |
| Direct categorical kernel | 0.908 |

Minimum BH q≈0.905 for the four strict source-signature tests; genus-only sensitivity likewise yields no adjusted global effect. This subset cannot be interpreted as fully representative because many records have unresolved native/introduced status. The original complete-trait matched-island records were ~52% origin-status unresolved.

## 5. Genuine mainland source data: GIFT 3.2 native checklists (separate reference)

The **same frozen progressive artifact** includes the earlier [GIFT 3.2 source-pool sensitivity](https://github.com/zuizui0223/island/tree/main) from `island` PR138:
- **748** source mainland regions with native lists, **1,666,186** native source-region × taxon rows.
- Existing outcome-blind `geo_k5`, `geo_k10`, `geo_k20`, `geo50_climate10` mainland assignment modes; all **74 matched islands** have assignments under all four.
- Exact GIFT `work_species` name match to the frozen PR321 complete-floral name population; no invented synonyms or source citations.
- The 37-pair `geo_k5` source reference covers a median **91.8–92.3% of island-recorded profile types** depending on evidence/origin scope. This is *broad categorical overlap*, not evidence of actual island ancestry. Continental islands generally share a larger fraction of their recorded **genera** with nearest mainland lists than oceanic islands (All medians **81.6% vs 72.2%**), suggesting source/lineage access is worth explicit modelling; this difference is descriptive.

### Source-repertoire-subtracted paired comparison (NOT a genus-null test)

For each island, build a normalized unique-profile incidence vector from the union of five nearest GIFT mainland native flora lists with resolved colour/form descriptions. Subtract this **potential mainland source** vector from the observed island vector before constructing the oceanic-minus-continental matched contrast. Perform 9,999 paired sign swaps; BH/BY over the four global geo_k5 residual tests.

| Source-subtracted residual | All upper-tail p | Direct upper-tail p | Interpretation |
|---|---:|---:|---|
| Exact-profile L2 | **0.0002** | **0.0002** | BH q=0.0004, BY q≈0.00083 in both |
| Categorical similarity kernel | 0.0429 | 0.3039 | BH q=0.0572 (All), 0.3039 (Direct): **no adjusted support** |

The exact-profile difference **persists after simple GIFT source-repertoire subtraction** under this model. However, this *does not preserve island genus composition, filter by ancestral dispersal ability, equalize mainland and island species diversity, or control trait-resolution bias*. The source repertoire is typically more diverse than a sampled island flora, and subtracting normalized vectors is **not** a probabilistic colonist-draw null. Hence this result does **not** overturn the genus-conditioned global non-rejection and cannot establish geological selection. `geo_k10`, `geo_k20` and `geo50_climate10` give the same **unadjusted exact-L2 residual direction** as sensitivity, but are correlated reference definitions.

## 6. The updated independent paper's message

**Demonstrated exploratory composition association:** origin-matched island floras differ in the exact recorded colour/form combinations.

**More discriminating explanation:** after conditioning on the taxon/genus incidence and observed regional availability, the global contrast is not exceptional; native-only observations also show no clear additional departure. This makes **taxonomic/source assembly of flower-profile repertoires** a viable explanation; it has not been uniquely identified as causal.

**Unresolved:** whether finer historical source-pool identity and taxon arrival opportunities, phylogenetic relatedness within genus, native status/missingness, or pollinator selection explain any residual regional divergence. Broad flower-colour/form similarity does **not** establish interchangeability of pollination function.

A careful article thesis is:

> *Geological island classes differ in the recorded combinations of floral phenotypes, but disentangling which plants arrive and persist from within-lineage floral change is essential for interpreting that difference.*

This remains an **assemblage-level observational study**, not evidence of within-species floral evolution.

## 7. Reproducibility and public/private boundary

- `src/fd/genus_source_null.py` accepts privately reconstructed PR321 species-profile maps and private island-flora records (see [replay contract](../../docs/GENUS_SOURCE_CONDITIONED_PROTOCOL.md)), but emits **only aggregate statistics and source SHA256**.
- `src/fd/mainland_source_residual.py` accepts the original GIFT mainland assignments and list records privately, again emitting **only summary tables**.
- Committed source outputs: [global null table](global_genus_null_summary.csv), [GIFT source residual primary tests](gift_source_residual_primary.csv), [GIFT coverage primary audit](gift_source_coverage_primary.csv).
- **Execution provenance:** source member identities and ZIP digests verified locally; 1,999 global/genus draws, 999 native-only draws, 9,999 GIFT source subtractions; scripts run against authorized user-owned frozen inputs. **The private-data test itself did not run on public GitHub Actions**; CI tests use synthetic source fixtures and validate the implementation without handling private rows.
- Missing provenance: full species→trait/source lineage cannot be redistributed under the [island database rights audit](https://github.com/zuizui0223/island/blob/3bfa818ffbc7a0dea4319933c4b0a0eb41cfabb5/docs/CHAPTER1_DATABASE_RIGHTS_AUDIT.md). Do not upload original source archives, taxon×profile maps, island×taxon rows, source entity IDs, geographic source-ID crosswalk, or private taxon names.

The next scientific step is to test whether **explicitly modelled genus/source-region assembly** predicts the observed repertoire at the paired-island level, with uncertainty due to incomplete taxonomic and trait coverage. An exact source- and phylogeny-controlled evolutionary test would require additional authorized genealogical/biogeographic information, not merely more category permutations.
