# v2 pipeline — reproducible GitHub Actions

Every stage of the v2 data pipeline is driven by a committed workflow or a
versioned CLI contract, so the whole analysis can be reproduced without hidden
manual transformations. The active set (all others were retired pilots and have
been removed):

| Stage | Workflow / command | Trigger | What it does |
|-------|--------------------|---------|--------------|
| 0. Island universe | `acquire_gshhg.yml` | manual / push to gshhg source | Build the exact island polygons (GSHHG 2.3.7 ≥ 5 km², Natural Earth 10m fallback) → `islands_v2.gpkg` + manifest + source policy. |
| 0.5 Bombus applicability freeze | `island-v2-bombus-applicability build` | manual, before outcome traits are inspected | Join the source-region registry, **source-region Bombus-status evidence**, island-to-source assignments, and the **separate island-assignment evidence** (geography/source-pool/dispersal) → every frozen island is explicitly `applicable`, `structurally_not_applicable`, or `unresolved`, plus Core/contingency/out-of-domain tables. Climate, Bombus records, floral traits, reproductive traits, and model results are rejected as inputs. Only rows a **human reviewer** has marked `accepted` are frozen; `agent_drafted_pending_review` rows stay `unresolved`. The current curation set is a **partial applicability registry** — the first Phase 0.5 review queue (3 source regions + 2 island assignments), not a complete freeze of the frozen universe; every other island stays `unresolved` by design. |
| 0.6 Core-pilot nomination | `island-v2-core-pilot-nomination nominate` | manual, after 0.5 and after collection | Nominate Core-pilot islands from **only** frozen applicability (`applicable` & not `regional_analysis_only`) **and** whether the current collection produced ≥ N **raw exact-island GBIF species labels**. This second gate is a **data-availability gate to START Core-pilot work, not analysis inclusion**: the count is raw collection availability, not reviewed flora richness, native establishment, or accepted taxonomic scope. Never uses traits, Bombus records, or model results. When none are eligible it reports the gap rather than inventing a Core island. A nominee still needs a flora anchor and taxon/establishment review before any analysis table. |
| 1. Plant submit | `submit-gbif-frozen-full-acquisition.yml` | `*/15` schedule / manual | Regenerate the 103-block plant request manifest from the frozen universe and submit up to `batch_size` pending blocks, capacity-aware, into the ledger. |
| 2. Plant poll | `poll-gbif-full-acquisition.yml` | `7,22,37,52` schedule / manual | Advance each active plant download `submitted → running → succeeded/failed`, record DOI. |
| (plant recovery) | `reconcile-gbif-full-acquisition.yml` | manual | One-off plant-ledger reconciliation helper. |
| CI | `validate-v2.yml` | PR / push | `ruff` + ontology validation + `pytest`. |
| 3. Plant collect | `collect-gbif-full-acquisition.yml` | `17,47 * * *` schedule / manual | Download one succeeded plant SIMPLE_CSV archive per run, assign occurrences to original exact island polygons, and update cumulative `island_species_occurrences.csv`, `island_observation_effort.csv`, and `island_taxa.csv`. |
| 3b. Selected-block collect | `collect-gbif-selected-block.yml` | manual only | Process one named succeeded, uncollected GBIF block through the identical exact-island collector. Input must include an outcome-blind geographic/source-region coverage rationale; traits, Bombus records/non-detection, reproductive outcomes, and model results are never inputs. The selected ID and rationale are written to `collection_status.json`. |
| 3.5 Flora diagnostics | `island-v2-diagnostics run` | — | Consume `island_observation_effort.csv` → per-island observation-process flags + realised-coverage `analysis_included` gate (`config/observation_diagnostics.yml`). Inclusion follows data coverage, not island area. |
| 3.6 Pollinator collect | `island-v2-gbif-pollinator-collect collect` | manual after a separately declared Apidae campaign has succeeded | Reuse GBIF block archives but retain every exact-island pollinator record with taxon, dataset, year, basis, coordinates, and exact island ID. It produces a raw, auditable input table rather than a biological absence table. |
| 3.7 Bombus diagnostics | `island-v2-bombus-diagnostics run` | manual after pollinator collect | Consume exact-island pollinator records → `detected` / `adequate_non_detection` / `insufficient_effort` / `unresolved` evidence under `config/bombus_observation_diagnostics.yml`. A missing record alone is never an absence claim. |
| 4a. Free trait source discovery | `discover-core-pilot-trait-sources.yml` | manual, only after a Core-pilot island passes nomination | For a single nominated island's staged taxa, collect free OpenAlex / Crossref / Unpaywall literature and public-access leads. The staged input must match the nominated island and preserve its release gate. Outputs are unreviewed M0 leads only. |
| 4b. Public PDF locator | `island-v2-trait-pdf-locator locate-pages` | after 4a | Download bounded, explicitly OA PDF candidates to memory, save only matched keyword page numbers, then discard PDF bytes. A page hit is a reading locator, never a biological confirmation. |
| 4.5 Phase 1 attrition audit | `island-v2-attrition run` | manual after applicability, Bombus diagnostics, and accepted trait evidence exist | Build direct-evidence trait coverage by island and sequential Core/M1/M2–M3 attrition table under `config/attrition_audit.yml`; no outcome models are fitted at this stage. |

The current scheduled GitHub Actions workflows concern the plant campaign only.
The separate pollinator campaign must have its own declared ledger and should be
run first for the predeclared Core subset, not automatically launched as a second
unbounded global download campaign.

## Pollinator acquisition contract

Prepare a declared Apidae acquisition manifest with the same frozen island
geometry and generic block builder, then submit and poll it through a separate
campaign ledger. The declared taxon resolution, catchment parameters, and block
membership must be retained alongside the downloads.

```text
island-v2-gbif-blocks prepare
  --islands-gpkg <frozen islands_v2.gpkg>
  --output-dir <apidae campaign directory>
  --taxon-name Apidae
  --taxon-rank FAMILY
  --target-scope apidae
  --grid-degrees <predeclared value>
  --max-wkt-chars <predeclared value>
  --max-islands-per-block <predeclared value>
  --query-buffer-m <predeclared value>
  --query-simplify-m <predeclared value>
```

The target-group choice is a measurement decision, not a claim that Apidae
exhaustively measures pollinator service. `Apidae` is the primary effort
background for the predeclared Bombus observation diagnostic; sensitivity
analyses can use an alternative target group only when a comparably auditable
record table exists.

## Resuming as data arrives

The plant campaign fills in over time, but the downstream flow does not wait for
all 103 blocks. `collect` is **idempotent and cumulative**: every run processes
one currently succeeded block and regenerates cumulative island-level outputs,
so the raw `island_taxa.csv` candidate list grows as more blocks land. It remains
an occurrence-derived candidate list, not reviewed flora richness.

The default scheduled collector selects the smallest uncollected succeeded block
to keep runtime bounded. A selected-block run is allowed only for a predeclared
geographic/source-region coverage reason. This changes processing order, not the
frozen island universe, GBIF query geometry, exact-polygon assignment, or any
trait/pollinator/model definition.

The applicability, pollinator, and attrition commands are likewise re-runnable.
They require versioned input tables rather than silently mutating the flora
outputs: (1) an outcome-blind source-region registry, (2) a separately assigned
pollinator occurrence table, and (3) reviewed trait evidence. This preserves
the boundary between data generation, measurement quality, and analysis
eligibility.

## Reproducibility note

The submit workflow currently pins the frozen island universe to a GitHub
Actions **artifact** (`FROZEN_ARTIFACT_ID`), which expires 2026-10-01. Before
then, the universe should be rebuilt deterministically from the committed
`config/island_source_gshhg.yml` via `acquire_gshhg.yml` (or the prepared
`islands_v2.gpkg` committed to the repo) so the campaign stays reproducible past
artifact expiry. Tracked as follow-up.
