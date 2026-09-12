# Chapter 1 NEE credential-independent island Search observation

## Purpose

This is a prospective operational fallback for the N1 island-observation side of the double-filter test when authenticated GBIF bulk downloads are unavailable. It does not relax the frozen biological evidence rules.

## Eligibility

Only island × channel rows whose independently constructed primary `geo_k5` source proxy is `available` are queried. Source availability is constructed without island channel observations, focal plant traits, or N1 effects.

## Frozen acquisition boundary

- Provider: GBIF Occurrence Search API.
- Geometry: exact polygons from the frozen GSHHG island universe.
- MultiPolygon handling: unchanged exact component polygons, ordered by area descending then stable WKT.
- No geometry simplification, buffering, or representative-point substitution.
- Shared maximum: 900 raw records per island × channel across every component and every frozen background subquery.
- Background groups remain the already frozen groups: Apidae for Bombus, seven bee families for non-Bombus bees, Lepidoptera, Aves, and Diptera.
- Target detections are restricted to the canonical confirmatory functional-channel catalog.

## Critical completeness rule

A target detection remains a valid `detected` observation even if the broader background retrieval terminates early after that detection.

A zero-target result can become `adequate_non_detection` only when **all** frozen background subqueries for **all** exact geometry components reach GBIF `endOfRecords` before the shared 900-record budget is exhausted, and the already frozen effort gate also passes.

Therefore:

- truncated zero-target search -> `insufficient_effort`;
- API/geometry error with no target detection -> `unresolved`;
- complete zero-target search with inadequate effort -> `insufficient_effort`;
- complete zero-target search with adequate effort -> `adequate_non_detection`.

A truncated Search API result is never promoted to `adequate_non_detection` and therefore cannot manufacture a `disrupted` channel state.

## Validation

Validation run `34701517523` passed:

- 8 synthetic fail-closed tests;
- Ruff;
- prospective completeness assertions;
- live taxonomic resolution for all frozen background groups.

Production workflow `run-chapter1-nee-island-search-observations.yml` was also syntax/no-op validated in run `34703521134`: with no trigger file, input preparation completes as a no-op and all observation jobs are skipped.

The production workflow requires a completed source-proxy workflow run ID in `data/v2/control/run_nee_island_search_observations.trigger`. No island occurrence outcome is queried before that source-side gate exists.

## Claim boundary

`detected` means that a potential partner-channel taxon was observed on the exact island under the frozen catalog. It does not establish realized visitation, effective pollination service, or functional replacement.

`adequate_non_detection` means only that the target was not retrieved under an exhaustive bounded Search retrieval with the frozen effort thresholds. It is not proof of biological absence.
