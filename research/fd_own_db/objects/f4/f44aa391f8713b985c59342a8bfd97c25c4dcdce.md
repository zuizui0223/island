# GloBI provider-provenance gate and leave-one-botanical-family-out audit

**2026-11-11 · source-pinned fd PR #14 · exploratory.** This is about **documented GloBI visitor-channel reference composition**, not actual island-local pollinator abundance or pollination efficiency. The remaining provider-independent prediction has NOT been validated; the original per-record provider/reference matrix is not in the archived genus summary.

## 1. Completed real-data result: 88 separate family holdouts

The earlier strictly unambiguous candidate GBIF family mapping retains **484 genera from 88 families**, with original Direct high/medium strictly single-state flower traits and GloBI genus source-reference channel shares. Using the exact same original source ZIPs and fixed 5-channel labels, train a model on EVERY FAMILY EXCEPT ONE and predict all genera of that unseen family; repeat deterministically over all 88 families. Each genus is scored exactly once; there is no intra-family training leakage. Train-only scaling, fixed ridge alpha100, and projection onto five nonnegative channel proportions summing to 1.

Comparison to the same citation-reference-effort / represented-plant / measured-plant-genus coverage predictors:

| Additional morphology | Genus-equal five-channel MSE gain | Family-equal MSE gain | Families with positive mean gain |
| --- | ---: | ---: | ---: |
| Open-radial form fraction only | **+0.001147** | **+0.001090** | **54/88** |
| Actinomorphic fraction only | +0.000414 | +0.000355 | 47/88 |
| Open-form + actinomorphy | +0.001306 | +0.001167 | 48/88 |
| Open-form + actinomorphy + green/brown fraction | +0.001266 | +0.001205 | 47/88 |

The first-source-check family-equal conditional bootstrap over **88 distinct family means** (4,999 draws) for **open-form ONLY** is **[+0.000093,+0.002076]**; for **open-form+actinomorphy** it is **[−0.000145,+0.002458]** (crosses zero). The same directional reduction occurs in mean absolute share error and Jensen–Shannon compositional divergence. Fixed ridge alpha10 does not reverse the direction. The results are not 88 evolutionary replications: botanical family histories, source/provider usage and taxonomic coverage are correlated, and trait feature/model choice was post-hoc.

This extends the prior 5-family-fold test without treating its 5 fold variants as independent biological data.

## 2. Evidence opportunity sensitivity

To probe whether the form signal is made solely by poorly documented GloBI genera or single-photo/single-trait-species genus summaries, rerun the same leave-one-family-out comparison after filtering genera by **independent GloBI interaction references** (3, 5, 10, 20, 50, 100) and by number of Direct complete morphology species per genus (1, 2, 3, 5). Require >=60 test genera and >=15 source families in each subcohort; never impute unknown traits.

There are **22 supported threshold combinations**; the genus- and family-equal open-form MSE gain is positive in **22/22**, though only **17/22** conditional equal-family bootstrap intervals have positive lower bounds. A stricter subset:

| Threshold | Genus | Families | Family-equal gain | Conditional 95% family bootstrap |
| --- | ---: | ---: | ---: | ---: |
| >=10 independent references | 366 | 77 | +0.001008 | [+0.000113,+0.001856] |
| >=20 independent references | 302 | 70 | +0.000842 | [+0.000165,+0.001520] |
| >=3 Direct morphological species/genus | 140 | 44 | +0.001211 | [+0.000353,+0.002162] |
| >=10 references and >=3 morphological species/genus | 115 | 41 | +0.000899 | [+0.000230,+0.001714] |
| >=100 references | 152 | 52 | +0.000734 | [−0.000026,+0.001535] |

The thresholds were explored post-hoc and alter the biological composition of the test set, so positive signs here do not give a single confirmatory p-value or document broad phylogenetic universality.

## 3. Major incomplete gate: genuine original provider-heldout GloBI validation

The archived `globi_source_genus_breadth.csv` has **2,228 genus-level rows** with five aggregated `refs__...` channel counts, `n_independent_references`, and GloBI reference volume. It **does NOT contain** original per-record source DOI, source dataset, dataset title, reference DOI/citation or publisher identifier. The separately retrieved `channel_taxon_catalog*.csv` gives pollinator-side evidence summaries; it does NOT link individual plant genus×visitor channel reference units to original provider IDs. Therefore **a genuine leave-one-original-provider-out prediction has NOT been run**. Reusing 88 family folds or applying a coarser reference-count permutation cannot substitute.

The frozen original GloBI source-of-truth is the 2026-06-04 version 0.9 Interpreted Data Products at DOI `10.5281/zenodo.20546682`, `interactions.tsv.gz` SHA256 **`095e57398420d1e0746e61233c796beedcfb8e03bd0c6a732aa17ae62043bd17`**, about 2.7GB compressed. The original rows contain source and reference citations and source-provider identifiers; these are needed to reconstruct provider-independent target counts. The GloBI source itself warns of geographic/taxonomic sampling bias and asks derivative works to cite original source contributors (not merely GloBI).

**New fail-closed provider contract** in `src/fd/globi_provider_provenance_gate.py` requires **verified original source SHA**, curated `plant_genus`, `visitor_channel`, `provider_id`, `reference_id`, `argument_type=SUPPORTS` and `evidence_class=pollination|flower_visit`. It rejects aggregated tables with missing provider IDs, invalid/refuted claims, and unknown reference provenance. It counts cross-provider reused reference keys and drops those citations from both sides before any source split. A second function privately reconstructs remaining genus×channel reference units after deleting the largest outcome-independent provider. Synthetic tests cover collisions and refusal to misrepresent the present aggregate data as provider-independent evidence.

When version-verified per-reference source rows can be obtained, first determine which original providers contribute to each tested genus and channel, deduplicate reference overlaps, then score the same common-genus family-heldout outcome table before and after removal of the largest provider and/or test on provider-independent channels. **Do not quote positive ecological source-provider results before this is actually run.**

## Scientific conclusion and limits

**Current best supported proxy result:** among the Direct-measured and GloBI-documented botanical genera considered, the **fraction of simple open-radial flower forms** shows a small, comparatively stable association with the **relative global literature reference composition** of five visitor groups, including in entirely untrained botanical families. The source-family LOFO test does not establish actual presence/absence of local island pollinator guilds or visitation effectiveness. The study is post-hoc and uses a selected sample of described floral trait species; source-provider concentration and true phylogenetic history remain untested.

**Unsupported claims:** true pollination-syndrome classification, actual moth versus butterfly effects (original GloBI channel collapses all Lepidoptera), floral functional diversity in realized pollen transfer, island pollinator loss, source/provider-disjoint generality, or adaptive causal species selection.

## Public code and confidentiality

- `src/fd/globi_family_lofo.py` and `tests/test_globi_family_lofo.py`: leak-free all-family-heldout proper scores, multiple error metrics and family-equal uncertainty; result is anonymous.
- `src/fd/globi_provider_provenance_gate.py` and `tests/test_globi_provider_provenance_gate.py`: original-source-preserving citation/provider input validation, required evidence classes and duplicate-reference collision guards. The latter has *synthetic-only* tests until the original source rows are licensed and processed.
- Released anonymous aggregates: `lofo_88_family_scores.csv` and `lofo_source_effort_morphology_thresholds.csv`; no original family/genus/plant/provider IDs.
- Protected source members and taxonomy lookup are not committed. Public GitHub CI tests the algorithms on synthetic data and cannot independently replay the original 484-genus private input.
