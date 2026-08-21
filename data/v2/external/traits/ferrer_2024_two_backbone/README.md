# Ferrer 2024 strict synonym recovery

This package records the 2026-08-21 audit of the 1,171 unique safe direct
Ferrer et al. (2024) names left unmatched after exact-master and strict GBIF
reconciliation. The source workbook is not redistributed; it remains
downloadable from `https://doi.org/10.1093/aob/mcae056` and is pinned by
SHA-256 `8ef2f5ac99780ca19a15b847442272457634064577f8a12fbb9710f5521e5899`.

`ferrer_two_backbone_responses_20260821.jsonl.gz` freezes the WFO June 2026
and GBIF responses for all 1,171 names. The derived mapping audit shows that
WFO and GBIF alone rescued no additional fixed-master species. Three WFO
accepted names that are in the fixed master were independently corroborated
by Kew POWO/WCVP species search results and retained in
`powo_wfo_selected_mappings_20260821.csv`. The accepted trait rows remain
Medium because Ferrer is a cited literature compilation rather than the
underlying experiment.

No fuzzy match, single-backbone match, family conflict, family inference,
global fallback, or genus-by-axis-only join is accepted.
