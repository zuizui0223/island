# Chapter 1 NEE source positive-scan transport audit — 2026-09-12

## Scope

This audit records an execution-layer failure and repair in the prospective source-side pollinator-channel availability scan. The scientific source-state contract was not changed: only a quality-filtered confirmatory GBIF occurrence can create `available`; no hit, geometry failure, or query failure remains `unresolved`; this scan never creates structural absence.

## Pinned inputs

- Frozen PR138 source-pool workflow run: `32954909953`
- Frozen source artifact: `9601703747` (`pr138-source-pool-sensitivity-32954909953`)
- Source artifact digest: `sha256:51cff9087a40b76655ff043d59aceabf25800a6029ed005a77b8db30e78345d2`
- Canonical GloBI SUPPORTS catalog workflow run: `34673860959`
- Catalog artifact: `10291643834`
- Catalog artifact digest: `sha256:c946e36cd04cb319fc63a5bde90534c1cf3c7f41480fc672638f85f5b2a89b02`
- Confirmatory catalog taxa: 6,928 across the five frozen channels.

## Failure observed in the first pilot

Pilot run `34697255116` used the outcome-blind source entities `345`, `346`, and `1384` across all five channels. Entity `345` returned `source_scan_error=RuntimeError` for every channel, while the same pipeline returned valid positive or unresolved states for the other pilot entities. The failure was therefore treated as a transport/infrastructure problem, not as biological absence.

Diagnostic run `34697492723` showed that GIFT entity `345` is a valid 72-part MultiPolygon. Its stable full WKT was 9,307 characters and produced a GBIF occurrence-search request URL of about 10,900 characters; all five channel queries returned HTTP 400.

Follow-up diagnostic run `34697572822` queried the unchanged component polygons separately. The full MultiPolygon still returned HTTP 400, while each of the first ten exact component polygons returned HTTP 200. These component WKTs were approximately 108–167 characters and generated request URLs of approximately 230–295 characters.

## Frozen repair

The transport strategy is now:

1. preserve the complete frozen GIFT 3.2 geometry with no simplification, buffer expansion, or representative-point replacement;
2. for a Polygon, query the Polygon as before;
3. for a MultiPolygon, query its exact component Polygons in deterministic order: descending component area, then stable WKT;
4. retain one shared maximum of 900 examined records per source entity × channel, divided deterministically across components;
5. stop after the first quality-filtered confirmatory target hit;
6. keep no-hit, component query failure, and exhausted fixed-budget cases as `unresolved`, never as structural absence.

Pure tests enforce the shared 900-record budget, deterministic component ordering, and fail-closed behavior after component errors. Live validation run `34700348137` passed the unit tests, Ruff, the frozen config assertions, and a real GBIF entity-345 component smoke query.

## Repaired pilot

Repaired pilot run `34700386669` passed all five channels with zero geometry failures and zero query-error entities. Entity `345`, previously a RuntimeError in 5/5 channels, was queried successfully and returned positive evidence in all five channels without changing the frozen biological eligibility rules.

Artifacts:

- Bombus: `10299588549`, digest `sha256:9904ec6331bd77891753541160f62a416fc1e2aaa33115c4e3955199ec3a8e90`
- non-Bombus bees: `10299733213`, digest `sha256:a7b826ece12ffd4dda852f2db192d055d0fda3d5f90f7c47127f5cd685865acf`
- Lepidoptera: `10299304369`, digest `sha256:c318d4a4a44e8bb30c6c779cdffb8266b13d9039fe6c3a32aee08a21303060c0`
- flower-visiting birds: `10300272304`, digest `sha256:3c7cb00a715da5a7f306cbd80d77089c7fd5472da202fce76d862d01f85bbdca`
- Diptera: `10299813551`, digest `sha256:464ae3ac1c56fd68037949ee4fc4782d125eac0a0b7e55b8d44a4ea44b605e21`

Entity `345` repaired positive hits:

- Bombus: `GBIF:6163123470`
- non-Bombus bees: `GBIF:5938501523`
- Lepidoptera: `GBIF:5938030199`
- flower-visiting birds: `GBIF:5938619449`
- Diptera: `GBIF:6133294529`

This change is interpreted only as repair of a query-transport failure. It is not a post-hoc biological source reassignment and it does not alter the source-state threshold, confirmatory taxon catalog, record budget, or channel definitions.

## Full scan

Full source positive-scan run `34700705549` was launched only after the repaired 5/5 pilot passed. Its frozen universe is 8,265 islands, 730 unique source entities across all four prespecified source modes, including 481 primary `geo_k5` entities. The five channel scans use the same exact-component strategy and fixed 900-record entity×channel budget.

The original full run used one serial 730-entity job per channel with `max-parallel=2`. Before any channel completed and before any full-scan source-state artifact existed, that execution layout was identified as the operational bottleneck. No source availability counts or channel outcomes had been observed from the full run at that point.

A deterministic execution-only sharding was therefore frozen before opening any full result:

- sharded workflow run: `34704844604`;
- source entity partition rule: sort the 730 frozen unique `entity_ID` values, then assign `index mod 10`;
- number of shards: 10;
- entities per shard: exactly 73 in every shard;
- all five channels use the same partition;
- each entity×channel is still processed by the unchanged `chapter1_nee_source_positive_scan_v1` scanner;
- the exact GIFT geometry strategy, confirmatory catalog, establishment filters and 900-record entity×channel cap are unchanged;
- shard membership uses no source occurrence outcome, island channel outcome, focal plant trait, N1 result or plant assemblage result;
- after scanning, the ten shard tables are reassembled and hard-checked to exactly 730 unique entity rows before a canonical channel artifact is emitted.

The sharded input artifact is `10300558673`, digest `sha256:9f434c2a3a20108df2305bed60403419930c98f0d53c01e9fc167cdc2fd87863`. Independent inspection confirmed 730 unique entities, zero duplicate entity IDs, and exactly 73 entities in each of the ten shards. Example previously diagnosed entities are separated by the outcome-blind rule (`345` in shard 8; `346` in shard 9).

## Canonical full-scan result lock

Canonical sharded run `34704844604` completed without a retry. All 50 scan shards and all five channel aggregation jobs completed successfully. Each full channel table contains exactly 730 rows, 730 unique source entities, zero duplicate entity×channel rows, zero structural absences, and zero GIFT geometry failures. Downloaded ZIP SHA-256 values matched the GitHub artifact digests exactly.

| channel | artifact | SHA-256 | available | unresolved | query-error entities | query-error + unresolved | multipart entities | max records examined |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Bombus | `10303941023` | `052db36a257ace02dfa4d4039eeb125477567fd1dad89e356f61433661ccddfe` | 492 | 238 | 39 | 34 | 240 | 48 |
| non-Bombus bees | `10303419534` | `47f0076e15fceed53b9d7a46c001a01e96bfb3d96d0b7861e0b0094ce4a1719f` | 676 | 54 | 38 | 19 | 240 | 900 |
| Lepidoptera | `10303921122` | `cc79318a7135cc8322d29606e36058ded1175b1fb1daca30ed10c7340da5ab45` | 712 | 18 | 39 | 13 | 240 | 179 |
| flower-visiting birds | `10303114973` | `afa31238543aa4bf9718da581612602061e8a66d292e2fa45330f39cbf7074e9` | 710 | 20 | 38 | 11 | 240 | 900 |
| Diptera | `10303509284` | `301d1c93b563d3fc570c67291114b8ef7f9eea504865a531ab4f569652e2ef11` | 668 | 62 | 38 | 19 | 240 | 900 |

`source_scan_error` is an audit flag, not a biological absence. If another exact component yielded a valid confirmatory positive, the entity may still be `available`; if no positive was established, the entity remains `unresolved`. No failed/error row was promoted to structural absence. The noncanonical serial run must not be unioned with these states and cannot be used to rescue any canonical unresolved entity.

These values are the fixed source-side inputs for the prespecified source-proxy projection. They are not a test of N1 by themselves, and no island retention/disruption outcome was used to choose or modify them.
