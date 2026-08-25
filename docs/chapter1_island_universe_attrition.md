# Chapter 1 island-universe attrition audit

## Answer

The 4,370-island Chapter 1 table is not the original physical-island universe.
The frozen GSHHG universe contains 8,265 candidate islands after the operational
area rule. The largest reduction occurs because 3,716 candidates have no plant
record assigned inside their exact polygon in the completed GBIF snapshot. This
is observation/data availability, not evidence that those islands lack plants.

The previously reported label `islands_lost_at_covariate_join = 83` was wrong.
Every island with a strict species is present in the 8,265-island covariate
table. The 83 islands disappear when strict flora are joined to the configured
direct traits: none of their strict species has a resolved value on any of the
eight configured trait axes. The subsequent covariate join loses zero islands.

## Verified stage counts

The sets are nested and keyed by the same geometry-derived `island_id`.

| Sequential stage | Islands | Lost at stage | Interpretation |
|---|---:|---:|---|
| GSHHG candidates with frozen geography/covariates | 8,265 | — | physical candidate-island universe; area 5–7,000,000 km2 |
| At least one exact-polygon GBIF plant record | 4,549 | 3,716 | no exact record in the snapshot for the remainder; not biological absence |
| At least one accepted-species aggregate row | 4,505 | 44 | raw records exist, but no accepted-species row was retained |
| At least one species in the strict 106,295-species universe | 4,453 | 52 | accepted species exist, but none enters the strict analysis universe |
| At least one strict species with any resolved configured direct trait | 4,370 | 83 | strict flora exist, but no species has a resolved value on the eight configured traits |
| After distance/area/climate/region covariate join | 4,370 | 0 | no whole-island loss at this join |

Area, distance, region, and spatial block are complete for all 4,370 analysis
islands. Climate PCs are missing for 20, so individual complete-case models can
lose observations even though the island remains in the materialized analysis
table.

## Selection is geographically structured

The retained fraction differs substantially across the frozen candidate
universe.

| Region | Candidate islands | Analysis islands | Retained |
|---|---:|---:|---:|
| Northern | 3,119 | 2,167 | 69.5% |
| Northern high latitude | 1,357 | 414 | 30.5% |
| Tropical | 3,013 | 1,495 | 49.6% |
| Southern | 776 | 294 | 37.9% |

Size is a major observation filter: final retention rises from 37.7% in the
smallest area quartile to 77.9% in the largest. The no-exact-record islands have
median area 11.7 km2, compared with 26.1 km2 in the final analysis set.

Distance-related selection is heterogeneous rather than globally monotonic.
Overall retention is 51.3%, 60.9%, 53.3%, and 45.9% from the nearest to the
farthest distance quartile. Within Northern islands it falls to 55.0% in the
farthest quartile, whereas within Southern islands it rises to 61.9%. A distance
coefficient can therefore be biased differently by region if the observation
process is ignored.

## Consequence for Chapter 1

The 8,265-island universe remains the denominator for the observation-selection
audit. The 4,370 islands are the currently trait-observable analysis universe.
Missing islands must never be coded as zero trait presence or as biologically
empty floras.

The primary pattern workflow should therefore add a frozen selection analysis:

1. model entry into the 4,370-island set from region, log area, distance, and
   climate without using trait outcomes;
2. report distance support shared by retained and omitted islands;
3. repeat pattern estimates with prespecified inverse-probability weights and
   overlap trimming, alongside the unweighted analysis; and
4. retain region-specific sensitivity because the distance-selection direction
   is not common across regions.

This audit does not manufacture trait data for unobserved islands. It measures
how far the observed analysis set may depart from the physical-island universe.

The implemented outcome-blind ridge-logistic audit converged in six iterations.
Of 4,370 included islands, 3,580 lie inside the frozen propensity overlap range
0.1–0.9. Stabilized selection weights have median 0.77, 95th percentile 1.99,
99th percentile 3.69, and a prespecified cap of 10. These weights are now
available for floral-pattern sensitivity analyses; they do not add pollinator
variables.

## Reproducible inputs

Counts were recomputed from the following immutable snapshots rather than copied
from a prose summary:

| Input | SHA256 |
|---|---|
| `purpose_shortest_island_data.csv` | `9aa10f6a9364a42be0a0c6b644581044c57986bdcc06a0bd912ded886469982a` |
| `island_observation_effort.csv` | `27dce385d83046a68f47484597e39e4beced8392c05e00342a75ea04c037c3f4` |
| `island_species_occurrences.csv.gz` | `ae33ad190a3f029bc625525a547b50dee532a0ab4137f98efc3d8f84b6a1216d` |
| `strict_species_axis_coverage.csv.gz` | `254e8df2f7d27bc065c849e025eaac18a77c2c204abbf4447e34f19b8f07b12f` |
| `island_trait_composition.csv.gz` | `ce3501bc929670646a4c778488ea8cc372741dff4c504cf00c5cf7b50293b59c` |

The reusable implementation is `island_universe_selection.py`; its focused tests
verify the nested-set contract and fail closed when a downstream island appears
outside its parent universe.
