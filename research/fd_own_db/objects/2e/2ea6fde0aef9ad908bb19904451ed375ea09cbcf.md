# 13 November 2026 · source-pinned island-by-GBIF temporal observation pilot

**Completed original-data result — NO statistically identified historical pollinator decline or guild switching.** Source runner [island Actions #43343152605](https://github.com/zuizui0223/island/actions/runs/43343152605) **SUCCESS**, fully anonymous [artifact #13999959972](https://github.com/zuizui0223/island/actions/runs/43343152605/artifacts/13999959972). The related code lives in [fd PR #15](https://github.com/zuizui0223/fd/pull/15) and [island PR #325](https://github.com/zuizui0223/island/pull/325).

## Why this pilot was necessary

Original 928-island GBIF five-guild qualified states report potential partner-channel presence, adequately surveyed *nondetection*, or unresolved. They are not a time series. Original GloBI v0.9 is also event-date-unavailable. To avoid mislabelling historical channel disappearance, we directly queried the public GBIF Occurrence Search API for original exact-island **dated Apidae** source records and Bombus genus documentation within them. The source-qualified 928 Direct morphological flora eligibility was fixed **before** querying insect outcomes.

The pilot selected **20 islands per each of the four predeclared geographic regimes (80 total)** with deterministic hashed selection favoring unique physical 10-degree spatial blocks. Selection never read flower metric values or pollinator statuses. Of 928 source-linked Direct floral island IDs, 927 had *exact original GSHHG WKB-hash-identical* polygons (source ZIP SHA256 `8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41`). The one unresolved original island is **unassessable**, not presumed guild-absent.

## Frozen original-source temporal contract

- Early: **2010–2017** inclusive; recent: **2018–2025** inclusive.
- Apidae accepted taxonomic identity resolved via GBIF at runtime, original exact polygon WKT components; background Search pagination <=300 and **900 total raw records maximum per island×period across polygon components**. Truncation/API error means not exhaustive, even if a taxon was not returned.
- Reject rows with `ABSENT` status, captive/cultivated/known-introduced entries, invalid year or original-polygon coordinates, duplicates and missing dataset key. No genus- or study-outcome-based taxon search refinements.
- An observable **dataset×year** requires >=3 original qualified Apidae records. A persistent dataset needs >=2 observable years in **both periods** in a single island; island-level temporal documentation can be compared only where searches in both periods are complete. A regional estimate further requires >=20 comparable islands and >=8 distinct blocks. These thresholds were frozen **before querying GBIF**.
- Measured contrast, *if estimable*, is proportion of observable dataset-years with at least one Bombus genus occurrence in the **same dataset** between periods — opportunistic documentation only, not occupancy/abundance, visitation, effective service or flower change.
- The API is a 2026-11-13 **live Search snapshot**, not a citable GBIF download DOI; it would need a new pinned GBIF download for confirmatory reproduction.

## Original-data support results (not a biological treatment effect)

| Broad region | Original sampled islands | Both period searches exhaustive | Islands with actual same-dataset repeated years | Early raw records examined* | Recent raw records examined* |
| --- | ---: | ---: | ---: | ---: | ---: |
| Northern high latitude | 20 | 17 | **0** | 9 | 70 |
| Northern midlatitude | 20 | 15 | **2** | 1,034 | 2,323 |
| Southern extratropical | 20 | 10 | **0** | 27 | 146 |
| Tropical | 20 | 7 | **1** | 107 | 374 |
| **Total** | **80** | **49** | **3** | **1,177** | **2,913** |

*Records examined under capped source queries are **NOT** independently sampled bees, abundance, visitation events or population trends. Some islands/windows are censored by the hard 900-record cap. A single GBIF record may reflect a different contributor/data-processing stream than earlier surveys.

Only **3/80 islands** have sufficient persistent dataset×year evidence under the declared observation filters. **0/4 regions** reaches the 20-island/8-block support gate. Accordingly, the pipeline publishes only ascertainment/feasibility counts, **no apparent decline rate, no regression to flower functional diversity and no biological treatment effect**. The much larger recent raw record total chiefly warns that time-varying recording opportunity is entangled with any naïve count trend; it cannot be used to infer a biological increase or decrease.

## Relationship to the main fd ecological question

The floristic source-relative marginal and joint floral diversity contrasts are **cross-sectional** with one fixed species-level trait label per accepted plant species. They cannot be paired as actual floral time-series response to 2010–2025 pollinator history. This original pilot directly falsifies the tempting shortcut 'GBIF year count difference = pollinator decline' in the present time-matched island study cohort.

To test historical decline and guild-switch impact, prioritize **standardized repeated site/plant flower-visitor datasets** with documented visits, observation duration, detection effort, and an observed repeated response. Sampling-event data and survey protocols may allow paired site-year effort, unlike this pooled opportunistic GBIF Search stream. Existing source/GBIF values can still support an explicitly **observational plant-assemblage source-relative functional diversification** paper without asserting a visitor-loss mechanism.

## Reproducibility

- Source [scripts/run_fd_gbif_temporal_apidae_pilot.py](https://github.com/zuizui0223/island/blob/analysis/fd-functional-guild-impact-20261113/scripts/run_fd_gbif_temporal_apidae_pilot.py) and [workflow](https://github.com/zuizui0223/island/blob/analysis/fd-functional-guild-impact-20261113/.github/workflows/run-fd-temporal-apidae-pilot.yml).
- Original [fd analytical module](../../src/fd/gbif_temporal_opportunity.py), [synthetic tests](../../tests/test_gbif_temporal_opportunity.py), and [pre-outcome frozen thresholds](../../docs/GBIF_TEMPORAL_POLLINATOR_PILOT.md).
- Raw original islands, exact geometry coordinates, source plant taxa, GBIF IDs and dataset identifiers remain within the authorized source session; the original artifact is anonymized regional aggregate only.

**Final status: TEMPORAL_DOCUMENTATION_CHANGE_NOT_IDENTIFIABLE. No historical Bombus decline, bee-to-bird switch, pollination service effect or coevolution was measured.**
