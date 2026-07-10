# v2 pipeline — reproducible GitHub Actions

Every stage of the v2 data pipeline is driven by a committed workflow or a
versioned CLI contract, so the analysis can be reproduced without hidden manual
transformations.

| Stage | Workflow / command | Trigger | What it does |
|---|---|---|---|
| 0. Island universe | `acquire_gshhg.yml` | manual / source change | Build exact island polygons and the frozen manifest. |
| 0.5 Bombus applicability freeze | `island-v2-bombus-applicability build` | manual / reviewed evidence change | Join source-region Bombus evidence and separately reviewed island-to-source assignments. Outcome traits and Bombus observations are prohibited inputs. |
| 0.6 Core nomination | `island-v2-core-pilot-nomination nominate` | manual / applicability change | Identify islands where a later Bombus-channel analysis may begin. This is not the active trait-acquisition selector. |
| 1–3. Plant acquisition | submit / poll / collect workflows | scheduled | Acquire exact-island plant occurrences and maintain `island_species_occurrences.csv`, `island_observation_effort.csv`, and the global unique-species master `island_taxa.csv`. |
| 3.5 Flora diagnostics | `island-v2-diagnostics run` | rerunnable | Audit the plant observation process. Missing or sparse records are never converted directly into biological absence. |
| 3.6 Pollinator collect | `island-v2-gbif-pollinator-collect collect` | after a separately declared Apidae campaign | Preserve exact-island pollinator records for effort-aware Bombus diagnostics. |
| 3.7 Bombus diagnostics | `island-v2-bombus-diagnostics run` | after pollinator collection | Classify `detected`, `adequate_non_detection`, `insufficient_effort`, or `unresolved`. |
| 4. Global trait campaign | `drive-global-trait-campaign.yml` / `island-v2-global-trait-campaign run` | hourly / manual | Advance one family-balanced wave across the global unique-species master. The ordered tasks are reproductive and pollen-vector evidence, floral-access traits for direct biotic-vector candidates, then alternative-pollinator guild evidence. |
| 4a. Regional source discovery | `discover-core-pilot-trait-sources.yml` | manual only | Optional deep reading for a named island or regional mechanism case. It no longer drives global acquisition. |
| 4b. Public PDF locator | `island-v2-trait-pdf-locator locate-pages` | manual / after source discovery | Save matched public-PDF page locators without treating page hits as biological confirmation. |
| 4.5 Attrition audit | `island-v2-attrition run` | after reviewed evidence exists | Quantify trait coverage and island attrition before any outcome model. |
| CI | `validate-v2.yml` | PR / relevant push | Validate workflows, lint, ontology contracts, and run the full test suite. |

## Global-first trait acquisition

The active automated path begins with `data/v2/staging/gbif/collected/island_taxa.csv`.
It does not nominate a preferred island group. Each wave is selected with a
deterministic family-balanced round robin, using only species identity, family,
and occurrence-derived availability counts.

The campaign proceeds through four ordered source tasks:

1. reproductive and pollen-vector statements from Wikimedia;
2. reproductive and pollen-vector statements from OpenAlex;
3. floral symmetry, form, tube depth, size, and colour for species carrying a
   direct machine candidate for biotic pollen vector;
4. focused alternative-pollinator functional-guild evidence.

The ledger is resumable. New species added to the global master are appended
without resetting completed work. A successful lookup with no matching statement
is a zero-hit lookup, not evidence that the trait is absent. Failed lookups are
retried and remain explicit in the audit state. Full policy and output contracts
are documented in `docs/global_trait_campaign.md`.

The Izu and other regional registries remain useful for mechanism-focused checks,
flora-anchor work, and later island-level joins. They do not determine which
species enter the global machine campaign.

## Pollinator acquisition contract

Prepare a declared Apidae acquisition manifest with the frozen island geometry
and generic block builder, then submit and poll it through a separate campaign
ledger:

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

`Apidae` is an observation-effort background for the Bombus diagnostic, not a
claim that bees exhaustively represent pollinator service. Alternative target
groups may be used only as labelled sensitivity tracks with equally auditable
record tables.

## Reproducibility note

The plant submit workflow currently pins the frozen island universe to a GitHub
Actions artifact that expires on 2026-10-01. Before then, rebuild it
deterministically from `config/island_source_gshhg.yml` or commit the prepared
GeoPackage so the acquisition contract remains reproducible after artifact
expiry.
