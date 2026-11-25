# fd prospective GBIF temporal opportunity pilot — 13 Nov 2026

**Outcome-blind design, not a historical pollinator-loss result.** Frozen original eligible universe: 928 islands with >=30 direct colour/form/symmetry species. Exactly 80 original polygon-verified islands are sampled without reading GBIF pollinator records: 20 per broad region, with spatial-block dispersion prioritized and SHA256 deterministic ranking used only to break ties. Original GSHHG WKB hashes must match the 2026 source checkpoint. No centroids, island buffers or simplified polygons.

## Two fixed windows

- Early: **2010–2017** inclusive.
- Recent: **2018–2025** inclusive.
- Use **Apidae** occurrences as the common observation background and the **Bombus genus** as a proxy target. GBIF taxonomy resolved once at run-time, by accepted name and family rank; do not choose another target or threshold based on outcomes.
- Fetch each exact original island polygon component and each window from the public GBIF occurrence Search API, pagination <=300 and MAXIMUM **900 raw records per island×window across components**. GBIF Search results >900 are *truncated*, never complete. Any failed query, unmatched polygon or missing year is **unassessable**, not absence. No claim of total Apidae sampling or biological absence when incomplete.
- Keep only original coordinate-matching, year-valid, noncaptive, non-living, non-ABSENT georeferenced occurrence records; deduplicate by GBIF ID; do not treat unspecified establishment as native.
- The source island metadata, detailed GBIF records, dataset IDs, coordinates, and private trait/guild records must never be committed or uploaded. Return aggregated statistical evidence only.

## What constitutes two-window within-sampling-frame comparison?

A **dataset×year cell** is observable only when >=3 qualified Apidae records are found that year in the exact original island; Bombus presence in that observable cell is an indicator of *documented occurrence opportunity*, NOT abundance or effective pollination. Require BOTH windows to be exhaustive. A dataset key must recur between windows. For a qualifying dataset, require >=2 observable years in **each window**. For a qualifying island, require >=1 such persistent dataset with >=2 observable years in both windows. A truly regional comparison additionally requires >=20 qualifying islands and >=8 distinct physical spatial blocks. Without the regional support gate, publish only coverage/feasibility counts. No MSE/ratio smoothing or "0 means extinct" fallback.

**Observed statistic** (if identifiable): within each qualified island × persistent GBIF dataset, the proportion of observable dataset-year cells with >=1 Bombus. Difference = recent − early. The period difference measures **GBIF documentation/recording pattern under stated effort filters**, not change in occupancy, abundance, flower visitation, or a before-after experimental treatment. Resample whole spatial blocks and datasets in sensitivity. A different-data-provider term is always a possible confound.

## Follow-on biology is explicitly separate

If >=20 islands and >=8 blocks have comparable repeated dataset histories AND the original flora/source matched Direct morphological outcomes, evaluate a frozen source-relative joint Rao Q and open-radial prevalence against historical documented Bombus probability **as a cross-sectional association**, controlling climate, area, regional/source fixed effects and evidence effort. Single global accepted-species morphology cannot establish within-species evolution or temporal changes in floras. To infer genuine pollinator shifts, need actual same-plant flower visits/visitor effectiveness and time or truly matched independent island–mainland populations. Never fit a causal FD-impact slope from this pilot when overlap fails.

Live GBIF API technical reference: https://techdocs.gbif.org/en/openapi/v1/occurrence (Search at most 300 per page; Search hard limit 100,000). The 900-record cap is predeclared and far below the limit, but inevitably censors heavily sampled islands.

## Output states

`source_pinned_temporal_effort_preflight`: complete/incomplete search fraction by region and period; persistent dataset/year coverage; explicit `not_identifiable` if support fails; no original identifiers. No claim of pollinator decline, pollinator shifts or effective ecosystem service.
