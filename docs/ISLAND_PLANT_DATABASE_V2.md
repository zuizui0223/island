# Island Plant Database 2.0

## Goal

Database 2.0 is a reusable, public **island × plant** database. It is not a Chapter 1 analysis matrix and it is not defined by the current floral-trait hypotheses.

The core object is the statement that a taxon is associated with an island, together with explicit evidence and establishment status. Floral and reproductive traits are optional linked layers.

Database 1.0 remains immutable and citable as the Chapter 1 species-trait analysis database. Its Zenodo DOI refers to the rights-filtered public derivative only. Database 2.0 does not overwrite Database 1.0.

## Design principle

Do not make the public database depend on whether every trait source can be redistributed.

Separate five objects:

1. **islands** — the geographic island universe;
2. **taxa** — normalized plant names and taxonomy;
3. **island_taxa** — island × taxon membership/status;
4. **occurrence_evidence** — auditable evidence supporting an island-taxon statement;
5. **traits** — optional taxon-level or island-population-level trait observations.

This makes the island flora useful even when a particular trait is missing or legally non-redistributable.

## Canonical tables

### `islands.csv`

One row per frozen island polygon.

Required fields:

- `island_id` — stable project identifier derived from the frozen island source/version;
- `source_island_id` — source polygon identifier;
- `island_name` — nullable human-readable name;
- `archipelago` — nullable;
- `country_or_territory` — nullable descriptive field, not the spatial key;
- `area_km2`;
- `centroid_lat`, `centroid_lon`;
- `geometry_source`;
- `geometry_version`;
- `geometry_sha256`;
- `release_status`;
- `source_license`.

The existing v2 GSHHG pipeline is the initial source. Geometry definition must remain versioned because the island universe is itself a scientific choice.

### `taxa.csv`

One row per normalized accepted plant taxon.

Required fields:

- `taxon_id` — stable project taxon identifier;
- `accepted_name`;
- `authorship`;
- `taxonomic_rank`;
- `family`;
- `genus`;
- `backbone_name`;
- `backbone_key`;
- `backbone_version`;
- `taxonomic_status`;
- `release_status`;
- `source_license`.

Raw names from occurrence sources never become canonical taxa without a recorded normalization decision.

### `island_taxa.csv`

One row per `island_id × taxon_id`.

Required fields:

- `island_id`;
- `taxon_id`;
- `membership_status` — `candidate`, `supported`, `accepted`, `rejected`, `unresolved`;
- `establishment_status` — `native`, `introduced`, `cultivated`, `transient`, `unknown`;
- `endemic_status` — `island_endemic`, `archipelago_endemic`, `nonendemic`, `unknown`;
- `first_record_year` — nullable;
- `last_record_year` — nullable;
- `occurrence_record_count` — non-negative evidence summary, not absence evidence;
- `specimen_record_count` — nullable/non-negative;
- `evidence_count`;
- `review_status`;
- `release_status`;
- `provenance_id`.

A GBIF hit creates at most a `candidate`. It does not by itself establish native presence, establishment, or biological absence.

### `occurrence_evidence.csv`

One row per evidence unit used to support or reject an island-taxon statement.

Required fields:

- `evidence_id`;
- `island_id`;
- `taxon_id`;
- `source_type` — e.g. `gbif_download`, `specimen`, `checklist`, `flora`, `literature`;
- `source_record_id` or `source_url`;
- `dataset_key_or_doi` — nullable;
- `basis_of_record` — nullable;
- `event_date` — nullable;
- `coordinate_uncertainty_m` — nullable;
- `evidence_role` — `supports_presence`, `supports_establishment`, `supports_native`, `supports_rejection`, `context_only`;
- `source_license`;
- `rights_status` — `redistributable`, `reference_only`, `review_required`;
- `review_status`.

For restricted/reference-only sources, retain identifiers and review decisions without copying protected prose or other restricted source content.

### `traits.csv`

Long-format trait observations. Do not force every trait into a single species-wide value.

Required fields:

- `trait_record_id`;
- `taxon_id`;
- `island_id` — nullable; null means species-level/general evidence, non-null means island/population-specific observation;
- `trait_name`;
- `trait_value`;
- `trait_unit` — nullable;
- `trait_ontology_version`;
- `evidence_id`;
- `quality`;
- `source_license`;
- `rights_status`;
- `release_status`.

Database 1.0 floral/reproductive values can be imported only through a versioned mapping that preserves their original values, quality labels, lineage provenance, and redistribution gate.

## Release layers

### Layer A — public flora core

Target: make this as complete as defensibly possible.

Contains:

- `islands.csv`;
- `taxa.csv`;
- redistributable `island_taxa.csv` facts;
- redistributable/reference-safe provenance identifiers and summaries.

This is the primary Island Plant Database product. It must remain useful without any floral trait columns.

### Layer B — public trait extension

Contains only trait observations whose values and provenance pass redistribution review.

Missing rows here mean **trait unavailable for public redistribution**, not plant absence.

### Layer C — research provenance registry

Contains review state and reference-only provenance needed to reproduce decisions. Restricted source content itself is not redistributed.

### Layer D — analysis products

Chapter-specific matrices, SEM/RDA/INLA inputs, syndromes, source-region assignments, and model-derived variables belong here, not in the public database core.

## Why this fixes the Database 1.0 release problem

Database 1.0 made the scientific unit a species × trait-axis cell. If the trait evidence was not redistributable, that cell could not enter the public release.

Database 2.0 makes **island × plant membership** the central public object. A species can therefore be fully represented in the island flora even when one or all floral traits are unavailable for redistribution.

The target is not `100% trait completeness`. The target is:

- near-complete auditable island/taxon coverage;
- explicit establishment uncertainty;
- transparent observation effort;
- trait coverage reported separately by trait and rights status.

## Existing v2 components that already fit

The current v2 rebuild already supplies much of the acquisition spine:

- GSHHG public island polygons → `island_manifest.csv`;
- exact-polygon GBIF requests;
- `island_species_raw.csv` as occurrence-based candidate names;
- explicit separation of taxonomy, establishment audit, and trait acquisition;
- raw pollinator data kept separate from biological interpretation.

Database 2.0 should therefore build on v2 rather than restart from `legacy/v1` or reshape Database 1.0 into an island flora.

## First build milestone

Freeze **Island Plant Database 2.0-alpha1** with no Chapter 1 trait requirement:

1. build the declared GSHHG island universe;
2. collect exact-polygon vascular/angiosperm plant candidates through the existing v2 GBIF workflow;
3. normalize taxonomy into `taxa.csv`;
4. aggregate evidence into `island_taxa.csv` while keeping `candidate` distinct from `accepted`;
5. generate table-level and row-level rights ledgers before any Zenodo release;
6. publish the flora core independently of the trait-extension schedule.

Only after alpha1 is frozen should Database 1.0 traits be mapped into the optional `traits.csv` extension.

## Non-negotiable boundaries

- no absence claim from missing occurrence records;
- no native/introduced inference from occurrence alone;
- no silent taxonomic synonym collapse;
- no copied restricted prose in public evidence tables;
- no global umbrella licence used to erase upstream obligations;
- no Chapter 1 hypothesis field required for database membership;
- no mutation of immutable Database 1.0.
