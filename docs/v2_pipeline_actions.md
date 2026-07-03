# v2 pipeline — reproducible GitHub Actions

Every stage of the v2 data pipeline is driven by a committed workflow, so the
whole thing is reproducible from the repository with no local steps. The active
set (all others were retired pilots and have been removed):

| Stage | Workflow | Trigger | What it does |
|-------|----------|---------|--------------|
| 0. Island universe | `acquire_gshhg.yml` | manual / push to gshhg source | Build the exact island polygons (GSHHG 2.3.7 ≥ 5 km², Natural Earth 10m fallback) → `islands_v2.gpkg` + manifest + source policy. |
| 1. Submit | `submit-gbif-frozen-full-acquisition.yml` | `*/15` schedule / manual | Regenerate the 103-block request manifest from the frozen universe and submit up to `batch_size` (default 3) pending blocks, capacity-aware, into the ledger. |
| 2. Poll | `poll-gbif-full-acquisition.yml` | `7,22,37,52` schedule / manual | Advance each active download `submitted → running → succeeded/failed`, record DOI. |
| (recovery) | `reconcile-gbif-full-acquisition.yml` | manual | One-off ledger reconciliation helper. |
| CI | `validate-v2.yml` | PR / push | `ruff` + ontology validation + `pytest`. |
| 3. Collect | `collect-gbif-full-acquisition.yml` | `41 */6` schedule / manual | Download succeeded SIMPLE_CSV archives, assign occurrences to the original exact island polygons, and write `island_species_occurrences.csv`, `island_observation_effort.csv`, and `island_taxa.csv`. |
| 4. Traits | `island-v2-traits run` (manual) | — | Consume `island_taxa.csv` (`accepted_species, genus, family`) → LLM web-search trait candidates for review. |

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

## Reproducibility note

The submit workflow currently pins the frozen island universe to a GitHub
Actions **artifact** (`FROZEN_ARTIFACT_ID`), which expires 2026-10-01. Before
then, the universe should be rebuilt deterministically from the committed
`config/island_source_gshhg.yml` via `acquire_gshhg.yml` (or the prepared
`islands_v2.gpkg` committed to the repo) so the campaign stays reproducible past
artifact expiry. Tracked as follow-up.
