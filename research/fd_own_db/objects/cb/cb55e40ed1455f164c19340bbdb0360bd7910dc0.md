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
