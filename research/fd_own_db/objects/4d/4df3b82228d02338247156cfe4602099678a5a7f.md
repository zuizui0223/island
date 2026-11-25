# Observation-effort robustness of geological island floral repertoire differentiation

> **2026-11-09 spatial-dependence correction:** After the effort screen, 21 matched pairs still share geographic blocks (five connected components, size 12/5/2/1/1). Exact component-level sign flips for the exact-category repertoire test yield **All p=0.2500** and **Direct p=0.3125**, versus the old independent-pair p=.0025/.0089. The BY-supported **All q=.0193** below is **not evidence of robustness to geographic pseudoreplication**; it is valid only under independent pair-label exchangeability. Do not claim confirmed geological-filtering significance. The effort-screened descriptive difference remains. [Spatial audit](../20261109_coverage_rarefaction/README.md).

**2026-11-08** | Independent `fd` research line | [method and source contract](../../docs/GBIF_EFFORT_SENSITIVITY.md) | [original 37 matched geological pairs](../20261108_matched_flower_repertoire/README.md)

**Result status:** completed exploratory observation-process sensitivity, with 47 synthetic tests passed in the [successful GitHub Actions run #42465727478](https://github.com/zuizui0223/fd/actions/runs/42465727478). The analysis does **not** identify ecological colonization, real biological absence, a calibrated detection process, or selection on flower traits.

## 1. Reproducible public effort source

The `island` project v2 collected Tracheophyta GBIF exact-island occurrences and separately saved a public 17-field, island-level effort table:

- Frozen commit: `3bfa818ffbc7a0dea4319933c4b0a0eb41cfabb5`.
- Path: [`data/v2/staging/gbif/collected/island_observation_effort.csv`](https://github.com/zuizui0223/island/blob/3bfa818ffbc7a0dea4319933c4b0a0eb41cfabb5/data/v2/staging/gbif/collected/island_observation_effort.csv).
- Source file SHA256 from the exact test run: `db044b73faab466b6f5f37e00471844b6cb55dc536ff378909ae214fb968e8d8`.
- **4,549** islands have source effort records from the earlier GBIF campaign. Across all of these, median 111 raw GBIF occurrence records per island, interquartile range 10–820, and median five distinct dataset keys.
- This prior campaign has **different taxonomic scope and time** from the later complete-floral reviewed snapshot. These fields are effort **proxies**, not independent species-level encounter histories.

Among All floral-support geological islands, **130/130 continental and 130/130 oceanic** had a previous effort row. Median unique GBIF IDs: **4,797 continental vs 3,870 oceanic**. Median dataset counts: **31.5 continental vs 29.5 oceanic**. Such broad effort summaries do not eliminate differential recording within genera.

## 2. Reconstructed matched cohort and quality screens

Reconstructed the exact same **37 outcome-blind matched geological continental/oceanic island pairs** used in [fd PR #4](https://github.com/zuizui0223/fd/pull/4), with 23 northern-midlatitude and 14 tropical pairs. All and Direct use the **same physical island pairs**, so these are nested evidence sensitivities, not independent studies.

Four declared **post-hoc** screenings:

| Screening | Pairs | North-mid | Tropical |
|---|---:|---:|---:|
| All original matched pairs | 37 | 23 | 14 |
| Both have finite GBIF IDs/dataset metadata | 37 | 23 | 14 |
| Both have ≥30 unique GBIF records and ≥2 dataset keys | 37 | 23 | 14 |
| Previous screen + record-count ratio ≤4 and both median spatial uncertainties ≤10 km | **21** | **12** | **9** |

The median log GBIF effort ratio (oceanic minus continental) across the 37 original pairs is **−0.6546**, i.e. oceanic has around **0.52×** the continental-side raw GBIF record count in the median matched pair. Effort balancing reduces this asymmetry at the cost of sample size. The result should not be interpreted as detection effort measured identically to the specific floral trait observation process.

## 3. Complete categorical repertoire association after effort screening

For each eligible island, construct the normalized binary incidence vector of the **same 496 recorded flower-colour state-set × floral-form state-set profiles**. Assess origin difference using both (a) exact-profile L2 and (b) the original full colour/form set-Jaccard similarity kernel. Fixed original matching is unchanged and the pair is the sign-swap unit; 9,999 Monte Carlo sign swaps per comparison. One 16-slot **post-hoc** multiple-test family (4 overlapping effort screens × 2 evidence scopes × 2 endpoints); BH/BY computed over all 16 slots.

| Selected subset and endpoint | n pairs | Direction vs paired sign-swap reference | p (upper-tail) | BH q (16) | BY q (16) |
|---|---:|---|---:|---:|---:|
| All original, exact-profile L2 | 37 | higher | 0.0001 | 0.00027 | 0.00090 |
| Direct original, exact-profile L2 | 37 | higher | 0.0001 | 0.00027 | 0.00090 |
| **Balanced strict, All, exact-profile L2** | **21** | **higher** | **0.0025** | **0.00571** | **0.01932** |
| Balanced strict, Direct, exact-profile L2 | **21** | higher | 0.0089 | **0.01780** | 0.06018 |
| Balanced strict, All, categorical similarity kernel | 21 | not supported | 0.4739 | 0.50549 | 1.000 |
| Balanced strict, Direct, categorical similarity kernel | 21 | not supported | 0.7555 | 0.75550 | 1.000 |

**Observed association survives a conservative effort symmetry screen as an exact recorded-category repertoire contrast**, most clearly in All. The full original `results/gbif_effort_robustness/effort_screen_repertoire_tests.csv` includes all 16 comparisons and their empirical null 95% ranges. **Full original dissimilarity-space separation is still not established.**

This post-hoc screen does **not** eliminate the possibility of source/geographic lineage composition, taxonomic detection differences, trait evidence imbalance or geological matching errors. The 21 selected pairs are not an independent random validation and differ in geographic balance from the original 37. No absent lineage or genus may be inferred from lack of a GBIF record.

## 4. Relation to the previous genus/source tests

These results complement, and do not refute, the [genus × region signature whole-flora conditional permutation](../20261108_genus_source_control/README.md), which showed **zero supported global excess** in its main four-test family (minimum BH q≈0.1847). The effort screen measures whether observed *composition association* remains under more comparable gross GBIF observation histories. The genus-conditioned null asks whether the association exceeds that expected from the fixed **taxon→island incidence and genus-specific profile repertoire**. They condition on different sources of confounding and have different null distributions.

The stronger developmental claim is **geographic and lineage assembly of recorded flower phenotypes** rather than a demonstrably source-independent adaptive floral niche shift.

## 5. Reproduction

- [Code](../../src/fd/effort_audit.py), [full limitations](../../docs/GBIF_EFFORT_SENSITIVITY.md).
- Executed in the existing [origin analysis workflow](../../.github/workflows/origin-analysis.yml) with source-pinned public aggregate floral data and independently sourced geological origin labels.
- Outputs: `effort_coverage_by_origin.csv`, `effort_pair_support_summary.csv`, `effort_screen_repertoire_tests.csv`, `manifest.json`, all saved as aggregated artifact only (no island IDs, taxon names, joined geological source entities or private species trait ledgers).
- No new annual species-by-island detection data were inferred from `year_min`/`year_max` alone; first observation cannot be treated as first biological arrival.

**Next identification gate:** Distinguish source-lineage assembly from trait-recording differences by obtaining independently defined historical source areas, phylogeny, and genus-specific survey/detection opportunity; absent these, true colonization failure cannot be identified from the GBIF name checklist.
