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
| 3. Collect | _(next)_ | — | Download succeeded SIMPLE_CSV archives and assign occurrence coordinates to the original exact island polygons + build the observation-effort table. |

Stages 1–2 share the `gbif-full-acquisition` concurrency group and rebase before
pushing the shared ledger (`config/gbif_full_acquisition_v2.json`), so scheduled
runs never race.

## Reproducibility note

The submit workflow currently pins the frozen island universe to a GitHub
Actions **artifact** (`FROZEN_ARTIFACT_ID`), which expires 2026-10-01. Before
then, the universe should be rebuilt deterministically from the committed
`config/island_source_gshhg.yml` via `acquire_gshhg.yml` (or the prepared
`islands_v2.gpkg` committed to the repo) so the campaign stays reproducible past
artifact expiry. Tracked as follow-up.
