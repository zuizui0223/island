# GBIF collection readiness audit

## Scope

This note records the current v2 status before the first full collection run.
It does not change the frozen `v1-freeze` branch.

## Confirmed implementation

The active v2 collector already:

1. downloads archives only for ledger entries whose GBIF status is `succeeded`;
2. parses GBIF `SIMPLE_CSV` records and retains occurrence-process fields;
3. assigns coordinates to the original, unbuffered island polygons rather than
   the buffered GBIF request catchments;
4. emits island-by-species candidates, a taxa hand-off table, and an
   island-level observation-effort table; and
5. runs unit tests before collection.

## Required hardening before scaling collection

- Include exact-boundary points (`covered_by`), not only strictly interior
  points (`within`).
- Deduplicate `gbifID` values across overlapping regional download blocks before
  species and effort aggregation.
- Write a per-block collection manifest with archive SHA256, source-row counts,
  valid-coordinate counts, exact-island assignments, buffer-only records, and
  duplicate-record counts.
- Stream large `SIMPLE_CSV` archives in chunks rather than loading a complete
  regional archive into memory.
- Expand observation effort to include unique GBIF IDs, basis-of-record counts,
  coordinate-uncertainty summaries, and establishment-means summaries.

## Acquisition scheduling rule

Do not rely on a push chain from commits made with `GITHUB_TOKEN` to start the
next workflow. Submit and poll must each retain direct scheduled triggers; the
poll/submit push chain is only an optional latency optimisation.
