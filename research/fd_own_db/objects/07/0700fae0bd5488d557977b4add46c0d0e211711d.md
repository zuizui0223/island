# Prospective external validation: standardized flower-visitor observation, NOT GBIF occurrence abundance

**13 November 2026 · source discovery and acquisition gate, no original dataset analysis completed.**

## Why new evidence is needed

The existing original dated GBIF Search preflight [Actions #43343152605](https://github.com/zuizui0223/island/actions/runs/43343152605) admitted only **3/80** outcome-blind island flora candidates with same-dataset repeated Apidae sampling years, and no region met the source support gate. The original 928-island static GBIF channel labels cannot identify genuine historical guild loss; species-wide floral trait records cannot identify in-situ plant phenotype evolution. The original GloBI v0.9 source is also undated and substantially aggregated by source archive.

## Verified public candidate: UK PoMS flower-insect timed counts

Original publisher: [UK Pollinator Monitoring Scheme — datasets](https://ukpoms.org.uk/data) / UK Centre for Ecology & Hydrology Environmental Information Data Centre (EIDC). Latest published source at verification:

- **Flower-Insect Timed Count data, 2017–2023**, 2026 publication, DOI **10.5285/d706088c-9d48-4be5-b156-6ed6ad5e0618**, catalogue https://catalogue.ceh.ac.uk/documents/d706088c-9d48-4be5-b156-6ed6ad5e0618 .
- Counts visitors **landing on target plant flowers** inside a **50 cm × 50 cm quadrat during a fixed 10-minute observation**, grouped by insect guild. Separate *public volunteer* FIT counts are ad-hoc selected locations from the *1-km square* FIT surveys, which target up to ~95 UK squares and up to four summer visits/year (with gaps).
- Publisher lists Open Government Licence terms with dataset acknowledgement and DOI citation. Dates and flower-target identities exist in the study's survey design, but the actual source CSV schema, location exactness at downloaded record granularity, per-island support and temporal repeat measurements have **not** been verified in fd.
- Main source download entry: https://data-package.ceh.ac.uk/data/d706088c-9d48-4be5-b156-6ed6ad5e0618 ; direct file access through the current research environment was not available during this check. **Do not invent source row counts, treat UK-only evidence as globally general, or claim sampled islands until the original file has been inspected.**

## Precise feasibility gate after authorized acquisition

1. Pin original ZIP/CSV SHA256, original publication DOI and license, years, survey mode. Deduplicate by source event ID and validate target flower taxonomy and recorded insect group.
2. Map coordinates to original GSHHG **exact original WKB polygon** IDs with verified spatial precision; record unassessable localities if reported grid-cell precision overlaps sea or multiple islands; do not use centroid or radius heuristics.
3. Keep only genuinely comparable 10-minute flower visitor counts. Treat volunteer public-location events separately from randomized 1-km-square events, and include flower target, date/season, visit effort and square/site as repeated strata.
4. Predefine focal insect effect: change in **within-site, within-target-flower**, standardized observed visitor group composition, not GBIF occurrence record share. Zero observed insects from a complete survey is a valid survey zero; missing survey is not an absence.
5. Gate for within-island longitudinal comparison: ≥2 sampled years per geographic island/site and flower identity, ≥2 distinct repeated sites per island, ≥10 islands across ≥5 spatial blocks or source clusters, and both comparison guild types supported. Source geography, flower availability, climate/year effects and baseline land use included; do not cherry-pick a detected guild change after outcome inspection.
6. Floristic FD relation: retain island-scale **same-carrier Direct reproductive and morphology** outcomes only if temporally matched or explicitly describe association as cross-sectional. No inference of historical flower evolution from globally fixed accepted-species trait labels.

## Claim ceiling

This UK dataset would offer an independent **actual flower-visitation event** domain relative to static GBIF and global GloBI documentary guild references, not pollen deposition/seed set, long-term evolutionary divergence, global replicated geography or a ready-made island survey panel. First task is a **source-licence/schema and exact-island overlap count**, not an island pollinator-shift p-value.
