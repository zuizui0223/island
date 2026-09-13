# Chapter 1 NEE island Search geometry transport audit — 2026-09-13

## Scope

This audit records a pre-outcome execution-layer geometry failure and its repair in the frozen Chapter 1 NEE exact-island Search acquisition. No island channel observation was opened before the repair, and no scientific source-state, channel, effort, Search cap, or observation-classification rule was changed.

## Canonical upstream state

The source-side positive scan and source-proxy projection were already frozen before island Search was launched:

- canonical source scan run: `34704844604`;
- canonical source-proxy run: `34729087812`;
- source-proxy artifact: `10308673426`;
- source-proxy artifact digest: `sha256:6288c60234a56de19201076394e21709503a7832cf01f8a6abeeec113189cc8d`;
- primary `geo_k5`: 8,265 islands × 5 channels = 41,325 rows.

## First island Search attempt — pre-outcome failure

The first sharded island Search attempt was workflow run `34731049099`. Its trigger pinned `source_proxy_run_id=34729087812`.

The prepare job completed source-proxy and GloBI catalog download successfully, but the workflow rebuilt the island geometry by calling `island-v2-gshhg build`. On that runner, live access to the legacy GSHHG distribution failed and the builder used its recorded Natural Earth 10m fallback. The resulting geometry contained 4,470 islands rather than the frozen 8,265-island analysis universe.

The prepare job then failed at the hard equality check between the rebuilt geometry and the frozen source-proxy island IDs. Consequently:

- no `scan-search-shard` job ran;
- no island-channel Search result was produced;
- no `detected`, `adequate_non_detection`, `insufficient_effort`, or `unresolved` outcome was opened;
- no N1 result was fitted.

The failure is therefore an execution/transport failure before outcome opening, not a biological result.

## Existing canonical island-universe artifact

The repository already contained `config/island_universe_artifact_lock.json`, which freezes the global island universe used by the full v2 analysis:

- workflow run: `28659417688`;
- artifact: `8066083419`;
- artifact name: `gbif-three-block-pilot-28659417688`;
- ZIP digest: `sha256:bee33e14672ec7ff4ed1f7acaea36b32cc6170a45c3d14b2caeb5f43bfaa623b`;
- source backend: GSHHG;
- GSHHG version: 2.3.7;
- resolution: high (`h`);
- minimum island area: 5 km²;
- source-file digest: `sha256:8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41`;
- frozen islands: 8,265.

Independent inspection of the artifact confirmed that its GeoPackage has 8,265 unique polygon island IDs in EPSG:4326 and that its island-ID set is exactly equal to the canonical `geo_k5` source-proxy island-ID set (zero IDs missing in either direction).

## Frozen repair

Before any island observation was opened, the Search execution contract was amended only at the execution layer:

1. canonical Search must download artifact `8066083419` rather than rebuild island geometry;
2. the artifact ZIP digest must equal the frozen digest above;
3. `source_policy.json` must identify the GSHHG backend, version 2.3.7, high resolution, minimum area 5 km², source-file digest above, and 8,265 islands;
4. the GeoPackage island IDs must exactly equal the 8,265 source-proxy island IDs;
5. Natural Earth or any other fallback is forbidden for canonical Search;
6. the deterministic partition remains sorted frozen island ID index modulo 20, before source-state filtering;
7. all Search and observation-classification rules remain unchanged.

The infrastructure retry boundary was also frozen before any island observation was opened: only a whole shard job that emitted no artifact may be retried once for runner timeout/infrastructure failure using the same code, inputs, channel, shard membership, and 900-record cap. Row-level retry of negative/error/unresolved islands and cross-run positive unions are forbidden.

## Corrected execution

Corrected Search run `34731405488` was triggered by workflow-repair commit `e505ae253b696c92b9375a632049b67fd59d4875`, still pinned to source-proxy run `34729087812`.

Its prepare job passed the frozen artifact digest check, GSHHG source-policy checks, exact 8,265-ID equality, and the deterministic 20-shard partition. The prepared input artifact is `10310100355`, digest `sha256:97628cc5dc270fa83b6a7cf28caa2230d27e75cd1ce5323fe822174c6a34f55c`, containing exactly five shards of 414 islands and fifteen shards of 413 islands, with no duplicate island IDs.

## Completed canonical Search

Run `34731405488` completed successfully: all 100 channel × shard jobs and all five channel aggregates succeeded without a row-level retry or cross-run union.

| channel | source available | detected | adequate non-detection | insufficient effort | unresolved | Search errors | aggregate artifact |
|---|---:|---:|---:|---:|---:|---:|---:|
| bombus | 7,154 | 824 | 23 | 5,983 | 324 | 324 | `10313926043` |
| non_bombus_bees | 8,186 | 1,088 | 11 | 6,697 | 390 | 390 | `10312514973` |
| lepidoptera | 8,265 | 2,113 | 0 | 5,758 | 394 | 394 | `10314040773` |
| flower_visiting_birds | 8,265 | 3,014 | 45 | 4,781 | 425 | 425 | `10313228947` |
| diptera | 8,265 | 1,261 | 7 | 6,605 | 392 | 392 | `10314045819` |

All five aggregate ZIP SHA-256 digests matched GitHub artifact metadata. Across the five completed tables there are zero duplicate island-channel rows, zero rows above the frozen 900-record cap, zero `adequate_non_detection` rows with incomplete Search, and zero non-detected Search-error rows promoted beyond `unresolved`.

The exact run, artifacts, digests, counts and fail-closed checks are frozen in `config/chapter1_nee_island_search_result_lock.json`.

## Transport bias limitation

The Search errors are not spatially random. Using the already frozen Chapter 1 covariate artifact `8270544465` (`sha256:695f35b97bae07e81b05deab537dd73fa687b2d99e9efdb1cb3babd2fa12dfb6`), overall Search-error rates are about 4.5–5.1% by channel, but errors are concentrated on large, geometrically complex islands. In the top decile of island area, channel-specific error rates are approximately 43.6–45.6%; in the bottom 80% they are below 0.5%.

This audit is reporting-only. Search-error islands remain `unresolved`; they are not re-queried, simplified, excluded to improve model fit, or unioned with another run. The transport audit cannot promote or demote a channel and cannot change the frozen N1 gate. Consequently any Search-based N1 inference must be stated as conditional on transport-evaluable exact island geometries.
