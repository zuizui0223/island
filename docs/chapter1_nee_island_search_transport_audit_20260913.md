# Chapter 1 NEE island Search transport audit — 2026-09-13

Status: reporting-only diagnostic. This document cannot change observation states, channel qualification, the N1 gate, the 900-record cap, island geometry, or canonical run membership.

## Canonical execution

- Search run: `34731405488`
- Search execution SHA: `e505ae253b696c92b9375a632049b67fd59d4875`
- Source-proxy run: `34729087812`
- Frozen island universe artifact: `8066083419`
- Frozen island universe digest: `sha256:bee33e14672ec7ff4ed1f7acaea36b32cc6170a45c3d14b2caeb5f43bfaa623b`
- Prepared Search input artifact: `10310100355`
- Prepared Search input digest: `sha256:97628cc5dc270fa83b6a7cf28caa2230d27e75cd1ce5323fe822174c6a34f55c`

The first pre-outcome attempt (`34731049099`) opened zero island-channel observations and failed because live GSHHG acquisition fell back to Natural Earth (4,470 islands). The corrected canonical execution uses the already frozen 8,265-island GSHHG artifact. See `config/chapter1_nee_island_search_execution.yml`.

## Bombus — complete 20-shard audit

All 20 Bombus shards completed successfully without retry. Independent reassembly of the shard CSVs produced exactly the frozen source-available set:

- rows: **7,154**
- unique islands: **7,154**
- duplicate island × channel rows: **0**
- missing expected source-available islands: **0**
- extra islands outside source-available support: **0**

Observation states:

- detected: **824**
- adequate non-detection: **23**
- insufficient effort: **5,983**
- unresolved: **324**

Search-completeness audit:

- complete searches: **6,005**
- truncated searches: **1**
- Search errors: **324**
- adequate non-detection with incomplete Search: **0**
- truncated non-detection incorrectly promoted above insufficient effort: **0**
- maximum raw records examined: **900**

The sum of all 20 shard receipts exactly matches the independently reassembled CSV counts above.

Under the already frozen channel qualification thresholds, Bombus alone has 847 evaluable source-available islands (= 824 retained + 23 disrupted) spanning 95 source-region IDs; this exceeds the per-channel confirmatory support minima. This does **not** open N1: primary N1 still requires at least three confirmatory channels.

## Bombus transport error structure

The 324 Search-error rows are all fail-closed `unresolved`. Error classes from the canonical shard outputs are:

- HTTP 400 Bad Request: **260**
- URL query too long after retries: **18**
- HTTP 429 / rate-limit failures after retries: **22**
- server disconnected after retries: **24**

The missingness is strongly geometry-dependent and is therefore not treated as random:

- median exact Polygon WKT length, non-error islands: **473 characters**
- median exact Polygon WKT length, error islands: **6,677 characters**
- mean canonical `log_island_area_km2`, non-error islands: **3.13**
- mean canonical `log_island_area_km2`, error islands: **7.73**
- mean canonical `log_distance_to_continent_km`, non-error islands: **3.77**
- mean canonical `log_distance_to_continent_km`, error islands: **4.61**
- area top decile: **312 / 716 = 43.6%** Search-error rate
- lower eight area deciles: **0** Search errors in the Bombus source-available subset

This is a material transport-evaluability limitation concentrated on large, geometrically complex islands. It does not justify row-level retry, geometry simplification, point/buffer substitution, a higher record cap, or union with another run. If N1 later passes its frozen gate, the claim must explicitly be bounded to the transport-evaluable island subset and this diagnostic must be reported.

## Remaining channels

non-Bombus bees, Lepidoptera, flower-visiting birds and Diptera remain in the same canonical run. Their transport diagnostics will be appended only after each channel completes all 20 fixed shards. No threshold will be selected from the observed error patterns.
