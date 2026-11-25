# Held-out source/genus prediction — distinct from a genus-conditional null

**2026-11-08 JST.** Research repo `fd` | upstream [island PR #321](https://github.com/zuizui0223/island/pull/321), frozen commit `d4a5dafab5f0737ea6c1c1c0263fc20decf6035c`.

## Question and observation unit

Given a known island's set of **recorded complete-floral taxon-name units and their genera**, how well can a mainland reference predict the held-out exact flower-colour × floral-form state-set profile for each name? Is this prediction improved by genus identity and/or geographic proximity of the source mainland flora?

This is a **conditional taxon-trait prediction** task, NOT prediction of which genera arrive on an island, native lineage evolution, pollinator effectiveness, or geological adaptation.

Same existing **37 continental–oceanic geographically matched island pairs** (23 north-midlat, 14 tropical) and 74 islands are used for All and Direct. No new islands or favorable profile categories are selected by the outcome. Exact profile state sets remain intact; do not expand multistate labels to fictional flowers.

## Train/test isolation and source evidence

Use privately reconstructed 16,739 All or 8,455 Direct source-validated floral profiles. GIFT 3.2 native mainland lists (748 mainland regions, 1,666,186 source-region × taxon records) are the only training flora. **All taxon names occurring in ANY of the 74 test islands are removed from the entire mainland training pool** within the corresponding evidence scope, including when the same taxon appears in multiple native mainland regions. Test names do not contribute to the global profile prior, genus model, or local source donor pool.

Training has 6,321 unique mainland floral-profile donors (All) or 4,172 (Direct) after whole-cohort name holdout; tested 13,141 taxon×island observations / 3,316 unique All names, and 11,664 observations / 2,605 Direct unique names. The reference profile catalog remains 496 unique colour/form state-set types (Direct 467 observed, against the common All reference).

The training input includes no original island origin labels except the outcome-blind choice of nearest mainland candidate source entities. Re-use of generic trait evidence across sources and biological relatives is possible; **taxon-name holdout is not phylogenetically independent validation**.

## Four probabilistic competitors

Let `p_global` be the marginal GIFT-native mainland frequency of each recorded complete profile (a tiny fixed **0.02 pseudo-count per one of 496 reference states** ensures nonzero probability), and `N_g` be the donor count in a genus.

- `global_mainland`: identical global frequency `p_global` for every held-out name.
- `local_mainland`: frequency among native GIFT names belonging to the island's assigned five nearest mainland source entities, with **8 global-prior pseudo-counts**.
- `global_genus`: frequency among other mainland GIFT taxa from the **same genus**, with 8 global-prior pseudo-counts; if no external donor genus exists, fall back to `p_global`.
- `local_genus`: frequency among other same-genus native GIFT taxa **in those five source entities**, with 8 `global_genus` pseudo-counts; if none exists, fall back to `global_genus`.

Each prediction sums to 1 over the fixed 496 states. The known test name's recorded **genus** is an input; the actual test name's flower traits are excluded from training. Evaluate: (a) proper categorical log loss on each held-out taxon name, (b) expectation of original mean state-set Jaccard colour/form distance, (c) approximate normalized predicted island incidence, scored by L2 and colour/form kernel distance. Expected binary incidence is calculated as `1 − product_taxa(1 − predicted_profile_probability)`; normalizing this expectation is an approximation, not an exact `E[normalized incidence]`.

Report equal-island mean of per-name scores; paired 37-cohort-cluster nonparametric bootstrap (4,999 draws) gives conditional percentile intervals for error improvement. **Identical taxon names may recur across matched islands**, so pair-bootstrap uncertainty is conditional on the selected repeated population, not an independent-taxon evolutionary sample. Treat these as model comparison diagnostics.

## Origin residual test — separate estimand

For each model, compute observed minus **predicted normalized island profile-incidence vector** for each island; compare residual oceanic minus continental repertoire vectors within the 37 fixed pairs using 9,999 paired label sign flips. The exact-profile squared L2 and category-distance-aware kernel are different endpoints. For the stated nearest-source genus model (2 scopes × 2 metrics), BH/BY over a fixed four-test family; global/mainland/genus-only models are descriptive sensitivities.

Even a residual paired-swap rejection **is not a source/genus-conditioned taxon-assignment permutation test**: predicted source donor coverage is limited, expected incidence uses a nonlinear approximation, names repeat, and unmeasured species-level history can remain. In particular, do not read this test as refuting the earlier [genus×regional availability null](../results/20261108_genus_source_control/README.md), which used one globally coherent profile randomization per taxon and found no pooled global BH/BY departure.

## Donor coverage gate

Record both the proportion of names whose genus occurs anywhere in heldout GIFT training and in the chosen nearby mainland source. In the complete 74-island cohort (equal island mean):

| Genus available in donor pool | All continental/oceanic | Direct continental/oceanic |
|---|---|---|
| Any native GIFT mainland source | 73.7% / 72.7% | 66.3% / 65.0% |
| Five-nearest native mainland sources | 38.6% / 36.6% | 29.2% / 27.4% |

Thus a substantial fraction of `local_genus` predictions fall back to a broader reference. Never call a high overlap of **profile types** (~92% in prior GIFT analysis) proof of high coverage of **same-genus source taxa**.

## Predefined primary, post-hoc sensitivities

The initial descriptive score comparison uses `geo_k5`, 8 pseudo-counts, original 37 pairs, All/Direct. Post-hoc sensitivity additionally varies `geo_k10`, `geo_k20` and pseudo-count 2/20. Genus-global improvement over global-mainland is positive under these smoothings; extra nearest-source benefits in **exact-profile log loss** depend on smoothing and can change sign. Extra nearest-source benefits in **original categorical trait distance** are positive but smaller than adding genus information. Do not select the most favorable smoothing based on test-set scores and present it as external confirmation.

An explicitly native-only version keeps **28 of the original 37 pairs** (20 northern midlatitude, 8 tropical) after requiring at least 20 complete explicitly native names per island under both nested All and Direct scopes; donor names are held out again for that reduced sample. Native/invader status is missing or unresolved for many source rows, so this sensitivity cannot represent the global native flora without qualification.

## Actual-source replay

Supply privately rebuilt and verified PR321 scoped `species_profiles.csv.gz`, private `chapter1_status_flora.csv.gz`, `GIFT 3.2 gift_native_mainland_flora.csv.gz`, `island_source_assignments.csv.gz`, and the locally matched `island_origin_crosswalk.csv`. Original data must not be committed to the public `fd` repository.

```bash
python -m fd.heldout_mainland_prediction \
  --private-all /secure/island/floral_space_20261106/all/species_profiles.csv.gz \
  --private-direct /secure/island/floral_space_20261106/direct/species_profiles.csv.gz \
  --flora /secure/island/chapter1_status_flora.csv.gz \
  --source-assignments /secure/island/island_source_assignments.csv.gz \
  --mainland-native /secure/island/gift_native_mainland_flora.csv.gz \
  --crosswalk /secure/fd/island_origin_crosswalk.csv \
  --g0 work/fd_g0 \
  --reference work/island_public/floral_space_20261106/all/profile_coordinates.csv \
  --source-mode geo_k5 --shrink 8 --bootstrap-draws 4999 --permutations 9999 \
  --output work/private_fd_predictions
```

For restricted taxon-level data, **public CI runs only synthetic tests**. This runner emits only `heldout_model_scores.csv`, `paired_prediction_improvements.csv`, `mainland_genus_support.csv`, `origin_repertoire_residuals.csv`, `manifest.json`. It must not publish original taxon-to-profile mappings, source area IDs, named island pairs, taxon names, or GIFT provider lineages. Explicitly rerun on held-out native-only and alternative source modes/smoothing as sensitivity checks.

## Interpretation gate

This analysis tests whether conditional prediction is better with knowledge of the island's genus identities. It **does not** predict which taxa or genera successfully colonize. Demonstrating an ecological mechanism requires a separate model of **arrival, source colonist sampling, naturalization and within-lineage change**, ideally with independent phylogeny and pollination interactions. Positive genus prediction and non-significant conditional-genus permutation results can both coexist: they target different questions.
