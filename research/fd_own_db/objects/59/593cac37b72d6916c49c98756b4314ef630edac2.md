# Observation-process sensitivity of origin-matched floral repertoires

2026-11-08. Independent `fd` analysis, post-hoc to the [37-pair continental/oceanic floral repertoire result](../results/20261108_matched_flower_repertoire/README.md) and the [genus-source results](../results/20261108_source_opportunity_hv17/README.md).

## What this can and cannot test

The original 37 geographically matched island pairs compare recorded complete-colour/form categorical profile mixtures, and the strongest earlier result was **a difference in exact profile repertoire**, not a supported difference in the original categorical-similarity geometry. A key confounder is uneven floristic observation effort. A source-only mainland genus model also distinguishes **recorded** from **unrecorded** genera, not confirmed presence versus true absence.

The `island` v2 GBIF collector has an independent, public, exact-island **Tracheophyta** effort table from an **earlier** collection campaign. It records raw GBIF ID count, data-source count, specimen/human observations, coordinate uncertainty and year span. It is useful to screen gross recording asymmetry, but it is neither (a) a repeated name-specific survey with documented zero detections, (b) the later reviewed angiosperm floral database, nor (c) a detection probability model.

**No undocumented nonrecord is converted into biological absence**, even on a well-observed island.

## Pinned public source and population

- Public `island` repository: SHA `3bfa818ffbc7a0dea4319933c4b0a0eb41cfabb5`, `data/v2/staging/gbif/collected/island_observation_effort.csv`. Hash of the actual downloaded bytes is included in the runtime manifest.
- The source collection status at that checkpoint reports **103 succeeded exact-polygon GBIF blocks**, **1,039,757 island × taxon name pairs**, **4,549 islands with GBIF records**, and **115,328 accepted species names**. These counts describe the source campaign and should NOT be equated with later complete-floral coverage (All 1,532, Direct 1,442).
- `fd` geography and colour/form profiles remain pinned to the original [island PR #321](https://github.com/zuizui0223/island/pull/321). The same independent Roeble/Zizka geological labels and **exact 37 outcome-blind paired islands** are reconstructed.
- The independent full PR4 matched set has 23 northern-midlatitude and 14 tropical pairs. All and Direct share the physical pairs; scopes are nested sensitivities, not independent confirmations.

## Fixed screens

The same two original endpoints use equal recorded-profile weights across original 496 reference profiles, with original full colour/form Jaccard similarity `1−D`:

1. `all_original_pairs` — reproduce original 37-pair exact-profile L2 and categorical-kernel label-swap reference.
2. `both_effort_available` — both islands have an actual metadata row with finite GBIF ID count and dataset count; no metadata is **missing, not zero**.
3. `min30records_2datasets` — both islands have ≥30 unique GBIF records and ≥2 independent GBIF dataset keys in that **older, whole-Tracheophyta source**. These are predeclared exploratory screening cutoffs, not biological detection calibrations.
4. `min30_balanced4x_uncertainty10km` — screen 3 plus ratio of `log1p(record count)` within log 4 and **both reported** median coordinate uncertainties ≤10 km. If either uncertainty is missing, the pair is excluded rather than falsely called precise.

For each screen/scope, require ≥8 whole distinct physical pairs; otherwise report `insufficient_pairs` and do not calculate a p-value. Compute **9,999 paired sign swaps** for each of two multivariate statistics with deterministic seeds. Report effect size, null 95% range, raw approximate p and a **single 16-slot post-hoc BH/BY family** (four screens × two nested scopes × two endpoints; unsupported slots assigned p=1). These are correlated post-hoc selections, not independent replication; family adjustment does not solve outcome-conditioned cohort selection or imperfect geography.

Also report how many matched pairs possess both older effort records, balance of their counts and regional representation. For context, compare among all geological origin-labelled complete-floral islands the fraction with any older GBIF effort metadata and median raw effort/dataset counts.

## Interpretation

If a difference persists in an effort-balanced subset, the simplest **overall GBIF record-count/dataset imbalance** does not fully account for the descriptive association under that subset. It does **not** eliminate genus-specific detection variation, source region dissimilarity, trait evidence collection bias, true distribution turnover, origin label uncertainty or incomplete intraspecific flower phenotype records.

If support disappears or the quality subset is too small, the original PR4 difference is **not disproved**; its robustness to observation effort remains unidentified at present.

Original documented observation-bias sensitivity in `island` Chapter 1 ([V6](https://github.com/zuizui0223/island/blob/3bfa818ffbc7a0dea4319933c4b0a0eb41cfabb5/docs/chapter1_v6_species_detection_tipping_result_20260914.md)) is an **assumption-bounded sensitivity**, not fitted species occurrence detection. Do not copy its numerical sensitivity results to `fd` as though they were recalculated for this research question.

## Rights and reproducibility

Public source file is pinned; only the hash and **aggregate** matched support/functional statistics are published. The workflow rebuilds geological matching locally, uses public **island × profile aggregate** carrier counts, and never attaches row-level island-geology crosswalk, individual paired island IDs, taxon source lineages or private donor checklists.

```bash
python scripts/fetch_island_aggregates.py --output work/island_public
python -m fd.pipeline --input work/island_public --output work/fd_g0
python -m fd.origin_analysis --g0 work/fd_g0 --hypervolume work/fd_hv \
  --crosswalk work/fd_crosswalk --output results/origin_strata
python -m fd.effort_audit \
  --input work/island_public --g0 work/fd_g0 \
  --crosswalk work/fd_crosswalk/island_origin_crosswalk.csv \
  --output results/gbif_effort_robustness --draws 9999
```

This step **cannot estimate species arrival, persistence, establishment, native status or annual colonization events**. `island` Database 2.0 declares `first_record_year`, `last_record_year` and `occurrence_record_count` as potentially useful provenance fields; without repeated independent survey opportunities and complete zero-record protocols, first record ≠ actual first arrival.
