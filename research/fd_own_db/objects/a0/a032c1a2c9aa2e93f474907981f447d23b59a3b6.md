# fd — functional differentiation of island floral assemblages

**Status:** independent exploratory analysis line, migrated from [island PR #321](https://github.com/zuizui0223/island/pull/321). This repository contains new analyses, **not** a new copy of the original plant-trait database. Its first runnable stage is implemented; no new global significance or evolutionary result is implied by creating the pipeline.

## Research question

**Can island floras occupy locally constrained combinations of recorded floral traits while remaining functionally differentiated from other islands?**

The target is *floral phenotype assemblage structure*, not demonstrated pollination effectiveness. The present two features are the recorded **flower-colour state set** and **floral-form state set**. A multi-state description stays one profile; it is never expanded into hypothetical colour × form flowers.

The main argument is tested in stages: (1) local spread, (2) between-island functional separation, (3) source/lineage/area-conditioning, then (4) ecological interpretation. Hypotheses and nonclaims are frozen in [docs/ANALYSIS_CONTRACT.md](docs/ANALYSIS_CONTRACT.md).

## What already exists upstream (not new fd findings)

From the frozen public [PR #321 report](https://github.com/zuizui0223/island/blob/d4a5dafab5f0737ea6c1c1c0263fc20decf6035c/results/floral_space_20261106/README.md):

- **496** distinct All profiles (467 Direct) from recorded colour × form state sets.
- Northern midlatitude **single-state** island profile richness is ~10.8% (All) and 11.4% (Direct) below genus-conditioned permutation medians (both BH q=0.016); this does **not** hold as an adjusted result for multistate-inclusive scopes.
- Genus-conditioned *shared-name-exclusive* between-island profile similarity has no BH/BY-supported effect in 16 comparisons.
- Whole-island **fixed-margins** profile Jaccard is lower than its approximate MCMC null in 16/16 comparisons, with BH 13/16 and BY 12/16 supported; this is a different endpoint **including shared taxa** and not conditioned on genus.
- An observed profile-turnover share of 0.690–0.902 (median across analyses) is *descriptive*, not above-null evidence.
- Gaussian KDE hypervolumes were **not** calculated. Two PCoA axes retain only 15.35% of the common positive inertia; 80% requires 17 axes. Do not call a 2-D KDE area total floral hypervolume.

## Local spatial predictability of floral profile repertoires (2026-11-09)

[44-source-archipelago heldout prediction report](results/20261109_archipelago_flower_transfer/README.md) · [Source/distance validation protocol](docs/ARCHIPELAGO_FLORAL_HOLDOUT.md) · [all 12 comparisons](results/20261109_archipelago_flower_transfer/archipelago_holdout_tests.csv)

**New evidence from outside the original 37-pair geological subset:** The original All/Direct common support includes 711 islands with named Weigelt archipelago labels, of which 606 from 44 named archipelagos can be evaluated by holding out each island and predicting its 496 recorded flower-colour/form profile distribution from three other islands in the same source group. Compared with the regional **outside-group** average, the within-archipelago predictor lowers exact-profile mean error from 0.006861 to 0.005077 (All) / 0.007351 to 0.005473 (Direct), and the original full categorical-kernel error from 0.010609 to 0.006203 (All) / 0.011137 to 0.006574 (Direct). Group-level sign flips give BH q=0.0003/BY q≈0.00093 across 12 correlated post-hoc test slots for these broad improvements.

**But this is almost entirely spatial, NOT a demonstrably special source-archipelago identity effect.** For **581/606 targets (95.9%)**, the three nearest physical islands are *exactly the same three islands* as the three nearest within the target archipelago. A **geography-only nearest-three** predictor, without using archipelago identity to choose donors, achieves virtually identical heldout profile error. The 44-group within-vs-nearest geographic improvement is not supported (All exact p=0.5913, Direct p=0.7033; original-kernel All p=0.8565, Direct p=0.9296). A stricter distance-matched different-archipelago control on 51 targets in 18 source groups likewise does not establish an additional named-source effect.

This supports **local predictability of recorded floral assemblage composition**, not geological trait filtering, evolutionary source ancestry, out-of-species transfer, confirmed functional pollination niches or true floristic absence. Nearby islands may share the exact same plant taxa. The next unresolved ecological question is whether profile similarity persists after **shared species/genera are excluded**, using authorized species-identity-level data only; the public aggregate alone cannot answer it.

## 2026-11-09 geographic source-disjoint rematching: no confirmed geological effect

[Independent-pair floral repertoire report](results/20261109_geography_independent/README.md) · [Source-disjoint matching contract](docs/GEOGRAPHIC_INDEPENDENCE.md) · [Exact test aggregates](results/20261109_geography_independent/independent_origin_repertoire_tests.csv)

The original 37 continental/oceanic island comparisons share only five geographically connected source groups. We therefore **re-matched the entire eligible geological flora population before examining flower traits**, maximizing different physical island comparisons under the original environment, area, isolation, name-support and coverage gates while prohibiting reuse of spatial blocks and/or source archipelago labels among pairs. This gives **13 block-disjoint pairs, 11 source-archipelago-disjoint pairs, 11 satisfying both, and 10 with both constraints plus verified archipelago labels**.

The strictest **10-pair** original categorical-profile exact L2 comparison is **All p=.1680 / Direct p=.2188** under complete `2^10` exact pair flips; across the **full post-hoc 16-test family none passes BH or BY** (smallest BH q=.2734). A nominal Direct p=.0408 in a **different 13-pair block-only subset** is not supported after multiple-testing correction. Original categorical-similarity-kernel comparisons likewise remain unsupported.

**This supersedes any claim of confirmed oceanic-island adaptive floral syndrome or source-independent structural divergence.** Observed category repertoires differ in particular samples, and known-genus/native-mainland opportunity predicts some recorded assemblage structure. But no currently available set of geographically disjoint pairs provides adjusted support for an independent geological-origin effect. More separate source archipelagos and identity-aware detection/phylogenetic controls are needed. Geographic labels are pragmatic proxies, not definitive phylogenetic independent histories.

## Critical 2026-11-09 spatial-independence audit — earlier pairwise significance is NOT confirmatory

[Exact conservative spatial-component analysis](results/20261109_coverage_rarefaction/README.md) · [Component graph and method](docs/COVERAGE_RAREFACTION_SENSITIVITY.md)

**The original 37 geological oceanic/continental matched island pairs are not 37 demonstrably independent geographical tests.** Under the frozen spatial-block graph, linking pairs that share either endpoint's geographic block gives only **5 connected components (23, 11, 1, 1, 1 pairs)**. All 23 northern-midlatitude pairs belong to ONE such connected component; a second coarse Weigelt archipelago grouping independently gives **5 components (18, 13, 4, 1, 1)**. Whole-component exact origin-label flips offer only 32 configurations, with minimum possible p=.0625 for the sign-invariant squared multivariate statistic.

**Actual original-profile exact-L2 results under conservative block grouping:** All exact p **0.2500**, Direct **0.3750** (37 matched pairs); in the earlier 21-pair GBIF-effort-balanced subset, All **0.2500**, Direct **0.3125**. Thus the previous pair-independent exploratory q-values in PR #4 and PR #8 **are not supported as spatially robust discoveries**. Raw category differences still exist, but cannot be claimed to exceed geographically dependent null expectations. Coarse blocks connected by graph transitivity are a deliberately severe stress test, not proof of fully dependent species histories.

The same 37 pairs do retain exact recorded colour×form composition differences after **equal-profile/complete-name rarefaction**, but this does NOT repair spatial pseudoreplication or the stronger **genus×regional-availability conditional null non-rejection**. The primary manuscript framing must shift to **source-lineage assembly and predictive, explicitly observational floral repertoires**, rather than claiming confirmed geological causal filtering or independent oceanic functional niche contraction.

## Sampling-support-standardized geological floral repertoire (2026-11-09)

[Validated 37-pair rarefaction report](results/20261109_coverage_rarefaction/README.md) · [Methods and bounds](docs/COVERAGE_RAREFACTION_SENSITIVITY.md)

**Exact recorded colour×flower-form combinations remain different after controlling the observed amount of information per island.** The frozen 37 paired continental–oceanic island floras all have ≥30 distinct recorded profiles and ≥40 complete taxon-name carriers. Under 99 repeat uniform observed-support draws and 1,999 paired sign swaps per draw:

- Exactly **30 distinct recorded profiles per island**: median exact-profile L2 p **0.007 All / 0.0035 Direct**; 83.8% / 93.9% of the repeated thinnings have raw p below .05.
- Exactly **40 sampled complete name carriers per island**: median p **0.007 All / 0.004 Direct**; 82.8% / 93.9% of thinnings below raw .05.
- Per-pair exact profile richness matched to the poorer side: median p **0.0005** in both nested evidence scopes.
- Restricting to ≥20 **single-state-only profile types** per pair retains **11 All / 10 Direct** pairs: exact-profile p 0.0517 / 0.093. This restricted subset has insufficient support for a comparable positive generality claim.
- The full original categorical **trait-similarity kernel** remains unsupported across these support screens. The distinction is largely exact recorded category-profile membership, **not demonstrated functional niche-space disjunction**.

**These 99 resamples are correlated thinnings of the same 74 original islands, not independent replications or inferred missing species.** This strengthens the evidence that simple observed sample size alone is insufficient, but it does **not** override prior genus×regional-source conditional null non-rejection, remove trait-dependent reporting bias or demonstrate geological adaptation/colonization. [CI #42554306259](https://github.com/zuizui0223/fd/actions/runs/42554306259) ran the full aggregate source pipeline with 50 tests passed.

## Origin-matched floral repertoire screen under independent GBIF observation-effort metadata (2026-11-08)

[Full effort-robustness results](results/20261108_gbif_effort_robustness/README.md) · [Pinned source/method contract](docs/GBIF_EFFORT_SENSITIVITY.md)

Rechecked the **same 37 original oceanic–continental island pairs** using the independent, earlier `island` v2 **Tracheophyta** exact-island GBIF effort metadata (record counts, dataset counts, coordinate uncertainty). All 37 pairs had older effort rows and ≥30 GBIF IDs/≥2 datasets on both islands. Requiring additional **≤4× record-count asymmetry and ≤10 km reported median coordinate uncertainty** retained **21 physical pairs** (12 north-midlat, 9 tropical).

The exact colour×form profile-repertoire separation remains supported in **All** (21-pair signed origin-label swap p=0.0025, 16-slot BH q=0.0057, BY q=0.0193), and more weakly Direct (p=0.0089, BH q=0.0178, **BY q=0.0602**). The original **categorical-similarity-kernel** separation is **not supported** in either scope (All p=0.474; Direct p=0.756).

**Strict interpretation:** this is a **post-hoc gross observation-effort screen** in a *different/earlier sampling campaign*, not a true per-genus detection model or independent replication. It does not establish biological absences, colonization failure, phylogenetic independence, or adaptive flower evolution. It also does **not** override the genus×geographic-availability conditional null non-rejection in [PR #5](results/20261108_genus_source_control/README.md).

## Mainland lineage opportunity × floral assemblage prediction and 17D projected KDE (2026-11-08)

[Corrected source / two-stage / 17D results](results/20261108_source_opportunity_hv17/README.md) · [Methods and privacy contract](docs/SOURCE_OPPORTUNITY_HYPERVOLUME17.md)

**Taxon-name leakage audit changes the interpretation.** Allowing GIFT mainland lists to include the **same species** present in tested island floras makes source–island overlap biologically descriptive, but not an independent out-of-species prediction. After excluding every observed taxon name from all 434 geology-classified islands globally from the mainland source lists, the earlier apparent nearby-source genus ranking improvement on continental islands **disappears/reverses**. The strict source-ranking support is 399 islands (174 continental / 225 oceanic) and 301 in explicitly native records (141/160). Its original 37 geographical pair differences are **not BH/BY-supported** over 12 tests (minimum BH q≈0.487). Unrecorded mainland candidate genera are **unknown, not true ecological absences**.

In the stricter independent two-stage flower-assemblage prediction, all names from the original 74 matched test islands are also excluded from GIFT **genus opportunity source features and floral trait donors**. The source-only genus mixture, without using observed island genus identities, lowers original categorical similarity-kernel **floral frequency** error from 0.020476 to 0.018107 (All, **11.6%**) and 0.024446 to 0.021582 (Direct, **11.7%**). Matched-pair bootstrap improvements remain positive; exact individual-taxon profile log-loss improvements do not consistently remain. The known-island-genus benchmark is an *oracle*, not prospective genus colonization prediction.

The new [17-dimensional KDE HDR implementation](src/fd/hypervolume17.py) goes beyond the previous two axes (15.35% inertia): 17 common PCoA axes retain **80.49%** of the 496-profile positive inertia. Runs use two bandwidths, 12,000 draws ×2 per island, 260 geological All and 243 geological Direct islands. Raw oceanic volume differences are largely accounted for by differences in the number of distinct recorded flower profiles, and no profile-richness-adjusted contrast has a 95% pair-bootstrap interval excluding zero. These are **projection-limited Gaussian-smoothed volumes, not total discrete floral functional hypervolumes or pollination niches**.

**Scientific boundary:** mainland source geography can constrain recorded plant lineage repertoires, but existing observation checklists cannot distinguish actual arrival, establishment, observation, native/introduced history or adaptation causally. Source-only flower-repertoire prediction is a modest assemblage-level association, not an ecological sorting mechanism proof. No restricted GIFT taxon identities/source rows or original island evidence are republished.

## Out-of-taxon mainland source and genus prediction (2026-11-08)

[Held-out predictive results](results/20261108_heldout_mainland_prediction/README.md) · [Private-input-only reproducibility contract](docs/HELDOUT_MAINLAND_PREDICTION.md)

The next test asks **how well island floral-profile identities can be predicted** when the island's observed genus composition is known. Every name recorded on any of the **74 original geology-matched island endpoints** is excluded from **all 748 native-GIFT mainland training regions**, eliminating exact taxon-name leakage. The independent predictors use either global mainland colour/form frequencies, nearest-mainland floral repertoires, known genus profiles from other mainland species, or nearest-source genus profiles.

**Main positive result:** Knowing genus reduces out-of-taxon exact categorical log loss by **0.709 All** (37-pair bootstrap range 0.644–0.774) and **0.338 Direct** (0.309–0.369) versus the global mainland distribution, alongside gains in original category-distance accuracy. **Nearest-mainland genus information improves categorical-distance losses** by about 0.014 All/0.011 Direct at fixed shrinkage 8, but does **not robustly improve exact-profile log loss**; its sign changes with smoothing. Global genus donor coverage is ~73% All/~66% Direct; five-nearest mainland genus coverage only ~38%/~28%.

Within 28 explicitly native-supported matched pairs, the genus predictive benefit remains. Residual exact-profile geological origin differences under the imperfect held-out predictive model do not contradict the separate, more rigorous genus×regional occurrence randomization **non-rejection** in PR #5. These tests **do not forecast taxon/genera colonization**, reconstruct true ancestral sources, or prove selection/pollination function.

## Genus and potential mainland-source conditioning (2026-11-08)

[Verified conditional analysis](results/20261108_genus_source_control/README.md) · [Rights-restricted replay contract](docs/GENUS_SOURCE_CONDITIONED_PROTOCOL.md)

The original **37-pair exact colour×form repertoire shift is observed**, but **is not clearly exceptional once the actual island×taxon incidence and genus×four-region geographic availability are held fixed** (1,999 conditionally shuffled whole floral profiles; pooled All exact L2 p=0.1035, Direct p=0.0675; none of the four global tests survives BH/BY, min BH q≈0.185). Genus-only global permutation likewise has no supported pooled contrast. The secondary tropical-genus residual is suggestive but not supported under the more restrictive regional-signature null's eight-comparison BH adjustment. An **explicitly native-only** sensitivity retains 28 original matched pairs and also shows no supported global excess.

The archived **GIFT 3.2** geographically near-native mainland source flora (748 source regions; 1,666,186 records) allows a separate test. Over 90% of the broad *recorded* profile types in each island are represented in the potential five-nearest mainland-source repertoire. Subtracting that normalized profile spectrum leaves the exact-profile origin contrast (All/Direct BY q≈0.00083), but not a consistent categorical-similarity-kernel contrast after FDR. **Source-spectrum subtraction is not the genus conditional null or an ancestral colonization model.**

Interpretation is now centered on **which taxa and floral repertoires become assembled** rather than a claim of evolutionarily smaller or functionally equivalent flowers on oceanic islands. The private trait source archives and taxon/island identity-level tables are never redistributed: only aggregate results, code and hashes are checked into fd.

## Matched floral repertoire results (2026-11-08)

[Matched repertoire report](results/20261108_matched_flower_repertoire/README.md) · [Protocol and inference limits](docs/MATCHED_REPERTOIRE_PROTOCOL.md)

The new analysis extends the geological-origin comparison to **which exact flower-colour × flower-form state-set profiles** are represented. In 37 outcome-blind 1:1 continental–oceanic pairs (23 northern-midlatitude, 14 tropical), the exact-profile repertoire composition differs under paired label swaps (**All/Direct BY q≈0.00042** across the fixed four global comparisons), but the **categorical-trait-similarity-aware** global test is not supported (All p=0.0828, Direct p=0.1609). No individual profile survives BH across all 496 tested categories. Colour/state presence is broadly shared across islands; the distinction lies in *combinations and their recorded prevalence*.

Geographically closer matched pairs (≤5,000 km, n=25; ≤3,000 km, n=14 northern-midlatitude) retain exploratory exact-profile differences. These are **post-hoc sensitivities**, not additional independent experiments. The result does **not** identify the source-pool/genus sorting process or establish geological selection, pollinator replacement or true functional niche equivalence; source identity remains the main unresolved gate.

## Geological origin follow-up (2026-11-08): continent-derived versus oceanic islands

The [geological-origin report](results/20261108_geological_origin/README.md) and [literature/novelty audit](docs/NOVELTY_LITERATURE_AUDIT.md) distinguish existing island trait-space research from the new question: **do geological histories change the *within-island* recorded floral dispersion and *between-island* functional separation relative to well-defined references?**

- Independently sourced Roeble (2024) and Zizka (2022) geological classifications mapped conservatively to **434/8,264** GSHHG islands, of which **260 All (130 continental / 130 oceanic)** and **243 Direct** meet the complete-floral support gate. Unmatched and mixed islands remain unresolved; continental distance is never used as a geological label.
- Separate Weigelt (2013) last-glacial mainland-connectivity data matched **1,334 All** support islands. A last-glacial unconnected island is **not** automatically geological oceanic.
- Oceanic All islands have lower observed equal-profile Rao Q (median **0.75705 vs 0.77447**). An exploratory area, distance, climatic, region, observed richness and trait-coverage adjusted coefficient is **−0.00589** (uncorrected spatial-block wild-bootstrap 95% interval [−0.00988, −0.00227]). Similar direction for Direct.
- **Crucial falsification:** For the 12 testable of 32 specified *origin-difference* tests under independent binary within-origin profile-margin nulls, **0 were BH/BY-supported**. Class-specific profile counts/prevalence can already reproduce apparent local Q origin differences; these observations do **not** establish enhanced packing, selection, within-lineage evolution, source-independent geological causation or pollinator functional replacement.
- The local and pairwise outcomes are **categorical colour/form phenotype proxies**. 2-axis KDE is a projection (15.35% positive inertia), not the full floral-functional hypervolume.

The paper's next identification gate is **profile-repertoire composition × source/genus/arrival constraints**, not treating any null rejection as a positive ecological mechanism.

## Implemented G1 and projected-hypervolume extension (2026-11-08)

The [verified analysis report](results/20261108_g1_hypervolume/README.md) records the successful 11-test CI and full results. Under the fixed binary profile margins, 13/32 exploratory comparisons pass BH adjustment, **none passes BY**. In northern midlatitudes, local equal-profile Rao Q is slightly lower while between-island nearest-profile dissimilarity is higher than the fixed-margin reference.

The new hypervolume implementation computes **finite categorical occupancy** and **90%-mass 2-D Gaussian KDE projected area** for 4,722 scope–island records, with globally fixed bandwidth/weighting sensitivities. Only 15.35% of the reference positive inertia lies in the first 2 PCoA axes, so these KDE areas are **not** the full functional floral-space hypervolume.

Definitions, nonclaims, and replay are in [G1 and hypervolume methods](docs/G1_HYPERVOLUME.md).

## Implemented fd stage: original-distance analysis

Original-distance metrics do not use PCoA; the supplementary Gaussian KDE is explicitly projection-limited and does not enter the functional-distance significance tests.

- **Local Rao's Q:** mean colour/form dissimilarity over pairs of recorded taxon-name units. Includes a finite-name version and an equal-profile-weight sensitivity.
- **Between-island functional nearest-profile dissimilarity:** symmetric, carrier-weighted average of the nearest recorded profile in the other island; zero for two identical profile supports even if carrier frequencies differ.
- **Cross-island mean profile distance:** carrier-weighted expected dissimilarity between two sampled name units, one from each island. This is **not** a pure between-island beta component.
- **Exact-profile Jaccard:** the original nominal overlap comparator.

All distances are computed from **equal-weight mean set-Jaccard** colour/form mismatches, not from PCoA coordinates. Pair comparisons reuse PR #321's already selected matched-pair IDs. They concern **whole-island** aggregate profiles, including any shared taxa. The raw taxon-excluded/within-genus tests remain in island.

## Replay from allowed public aggregates

Python 3.11+:

```bash
python -m pip install -e ".[test]"
pytest -q
python scripts/fetch_island_aggregates.py --output work/island_public
python -m fd.pipeline --input work/island_public --output results/run_fd_v1
```

The importer pins source commit `d4a5dafab5f0737ea6c1c1c0263fc20decf6035c` and fetches **only** public aggregate island-by-profile carrier tables, island metrics, the 496-profile reference, and the public matched-pair cohort. Actual file checksums go to `source_manifest.json`; the runner independently records input hashes in `manifest.json`.

Outputs: four `*_islands.csv` files, four `*_pairs.csv` files and `manifest.json`. All/Direct and single-state subsets overlap and are sensitivity scopes, not independent replications. No species-level traits, private source ledgers or restricted redistribution material enter this repository.

The [GitHub Actions workflow](.github/workflows/analysis.yml) runs tests and attempts the pinned public-data analysis, attaching its outputs as a workflow artifact rather than automatically committing them.

## Claims gate

The public-aggregate stage provides **descriptive measured distances only**. It cannot estimate genus/source-conditioned null significance, exclude shared taxon names, identify pollinators, prove ecological functional equivalence, establish evolutionary divergence, or attribute differences to selection.

Next inference stage: implement function-distance nulls while retaining the exact existing PR #321 selection and distinguishing a fixed-profile-margin test (public aggregates suffice) from a genus/source-conditioned test (requires rights-cleared private input). Then audit area, isolation and trait-coverage effects on the *same* cohorts before manuscript claims.

**Ownership:** `island` is the frozen data/Chapter 1 provenance and PR #321 exploratory origin. `fd` owns the new distance-based floral-functional diversity question, implementation, validation and any eventual separate manuscript.
