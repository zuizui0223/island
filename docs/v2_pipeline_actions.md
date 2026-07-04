# v2 pipeline — reproducible GitHub Actions

Every stage of the v2 data pipeline is driven by a committed workflow or a
versioned CLI contract, so the whole analysis can be reproduced without hidden
manual transformations. The active set (all others were retired pilots and have
been removed):

| Stage | Workflow / command | Trigger | What it does |
|-------|--------------------|---------|--------------|
| 0. Island universe | `acquire_gshhg.yml` | manual / push to gshhg source | Build the exact island polygons (GSHHG 2.3.7 ≥ 5 km², Natural Earth 10m fallback) → `islands_v2.gpkg` + manifest + source policy. |
| 1. Submit | `submit-gbif-frozen-full-acquisition.yml` | `*/15` schedule / manual | Regenerate the 103-block request manifest from the frozen universe and submit up to `batch_size` (default 3) pending blocks, capacity-aware, into the ledger. |
| 2. Poll | `poll-gbif-full-acquisition.yml` | `7,22,37,52` schedule / manual | Advance each active download `submitted → running → succeeded/failed`, record DOI. |
| (recovery) | `reconcile-gbif-full-acquisition.yml` | manual | One-off ledger reconciliation helper. |
| CI | `validate-v2.yml` | PR / push | `ruff` + ontology validation + `pytest`. |
| 3. Collect | `collect-gbif-full-acquisition.yml` | `41 */6` schedule / manual | Download succeeded SIMPLE_CSV archives, assign occurrences to the original exact island polygons, and write `island_species_occurrences.csv`, `island_observation_effort.csv`, and `island_taxa.csv`. |
| 3.5 Flora diagnostics | `island-v2-diagnostics run` | — | Consume `island_observation_effort.csv` → per-island observation-process flags + realised-coverage `analysis_included` gate (`config/observation_diagnostics.yml`). Inclusion follows data coverage, not island area. |
| 3.6 Bombus diagnostics | `island-v2-bombus-diagnostics run` | manual after pollinator records are assigned to the exact islands | Consume an island-assigned pollinator occurrence table → `detected` / `adequate_non_detection` / `insufficient_effort` / `unresolved` evidence under `config/bombus_observation_diagnostics.yml`. A missing record alone is never an absence claim. |
| 4. Traits | `island-v2-traits run` (manual) | — | Consume `island_taxa.csv` (`accepted_species, genus, family`) → LLM web-search trait candidates for review. |
| 4.5 Phase 1 attrition audit | `island-v2-attrition run` | manual after applicability, Bombus diagnostics, and accepted trait evidence exist | Build direct-evidence trait coverage by island and the sequential Core/M1/M2–M3 attrition table under `config/attrition_audit.yml`; no outcome models are fitted at this stage. |

Stages 1–2 share the `gbif-full-acquisition` concurrency group and rebase before
pushing the shared ledger (`config/gbif_full_acquisition_v2.json`), so scheduled
runs never race.

## Resuming as data arrives

The campaign fills in over ~a day, but the downstream flow does not wait for all
103 blocks. `collect` is **idempotent and cumulative**: every run reprocesses
whichever blocks are currently `succeeded` and regenerates the full island-level
outputs, so the `island_taxa.csv` species list simply grows as more blocks
complete. The handoff to trait acquisition is therefore stable — point
`island-v2-traits run --taxa-csv .../island_taxa.csv` at the latest collect
output and re-run it on the enlarged list whenever more blocks land. No stage
needs the campaign to be finished before it can start on the data in hand.

The pollinator and attrition commands are likewise re-runnable. They require
versioned input tables rather than silently mutating the flora outputs: (1) an
outcome-blind applicability table, (2) a separately assigned pollinator
occurrence table, and (3) reviewed trait evidence. This preserves the boundary
between data generation, measurement quality, and analysis eligibility.

## Reproducibility note

The submit workflow currently pins the frozen island universe to a GitHub
Actions **artifact** (`FROZEN_ARTIFACT_ID`), which expires 2026-10-01. Before
then, the universe should be rebuilt deterministically from the committed
`config/island_source_gshhg.yml` via `acquire_gshhg.yml` (or the prepared
`islands_v2.gpkg` committed to the repo) so the campaign stays reproducible past
artifact expiry. Tracked as follow-up.
