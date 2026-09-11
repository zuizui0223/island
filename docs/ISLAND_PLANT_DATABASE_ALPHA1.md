# Island Plant Database 2.0-alpha1

## What alpha1 is

`2.0.0-alpha1` is the first materialized **island x plant taxon** database bundle.
It is intentionally broader than the Chapter 1 trait analysis and narrower than a
verified island flora checklist.

The frozen build target is:

- 8,265 GSHHG 2.3.7 island polygon units at the primary 5 km2 threshold;
- 115,328 provisional plant name identities recovered from the completed GBIF campaign;
- 1,039,757 candidate island x taxon rows;
- 1,039,757 aggregate evidence rows;
- 4,549 islands with at least one GBIF-supported candidate taxon;
- 103 GBIF occurrence downloads represented in the rights ledger.

## Scientific boundary

A GBIF record is evidence that a taxon name has been observed and spatially assigned to
an exact island polygon. It is **not**, by itself, evidence that the taxon is native,
established, endemic, or part of a complete checklist.

Therefore alpha1 freezes all GBIF-derived island x taxon membership as `candidate` and:

- does not infer absence from missing GBIF records;
- does not infer native or introduced status;
- does not infer endemicity;
- keeps the occurrence-derived name identity provisional until explicit taxonomy
  normalization is run;
- keeps trait coverage optional.

## Dateline provenance

GBIF query construction can split an island that crosses the antimeridian into multiple
query-side catchments. Those synthetic catchment IDs are never database island IDs.
The builder collapses them through `analysis_island_id` and records every GBIF download
that could contribute evidence to the original island.

## Rights boundary

The island geometry core is redistributable under the source-specific GSHHG or Natural
Earth terms recorded by the island builder.

GBIF occurrence downloads can contain constituent datasets licensed under different
GBIF-supported licences. Alpha1 therefore retains the candidate-flora rows for database
construction but marks their value redistribution `review_required` until
constituent-dataset licence lineage has been recovered. Download DOI/key provenance is
retained so that audit can be completed without repeating the occurrence campaign.

## Canonical alpha1 artifact

The build workflow emits:

- `islands.csv`
- `taxa.csv.gz`
- `island_taxa.csv.gz`
- `evidence.csv.gz`
- `RIGHTS_LEDGER.csv`
- `DATABASE_MANIFEST.json`
- `VALIDATION_REPORT.json`
- `BUILD_RECEIPT.json`

`BUILD_RECEIPT.json` records SHA-256 and byte size for every frozen bundle component.
The CI build fails rather than emitting the artifact if the exact island universe or the
frozen row counts change.

## Next database upgrades

The next upgrades are independent layers, not reasons to mutate alpha1 silently:

1. resolve provisional taxon names to explicit backbone identifiers and accepted/synonym
   status;
2. recover constituent GBIF dataset licences and release-safe evidence lineage;
3. add native / introduced / endemic status from auditable checklist or flora sources;
4. attach Database 1.0 floral and reproductive traits as a species-trait extension;
5. enrich island names, archipelagos, and country/territory membership without changing
   geometry identity.
