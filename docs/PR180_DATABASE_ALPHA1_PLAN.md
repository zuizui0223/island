# PR180 — Database 2.0 alpha1 candidate-flora materialization

This PR advances Island Plant Database 2.0 from schema-only infrastructure to a
materialized candidate-flora database.

## Frozen target

The PR must build and validate all of the following before merge:

- 8,265 GSHHG 2.3.7 islands;
- 115,328 provisional taxa;
- 1,039,757 island x taxon candidate rows;
- 1,039,757 aggregate GBIF evidence rows;
- 4,549 islands with exact-assigned GBIF occurrence rows;
- 4,505 islands with species-level candidate-flora rows;
- 44 exact-occurrence islands with no retained species-level candidate row;
- 103 GBIF download provenance objects.

The 44-island gap is an explicit unresolved coverage state and must not be interpreted as
plant absence.

## Non-claims

The build does not promote GBIF occurrence to verified flora membership, native status,
endemicity, or absence evidence. It does not alter Database 1.0, Chapter 1 H1-H5, or
legacy/v1.

## Merge gate

Merge only if the dedicated alpha1 workflow:

1. rebuilds the exact GSHHG island universe without fallback;
2. regenerates GBIF block membership and collapses antimeridian query parts through
   `analysis_island_id`;
3. materializes the candidate flora tables;
4. validates relational keys and provenance references;
5. matches all frozen counts exactly, including the 4,549 -> 4,505 -> 44 coverage split;
6. uploads a bundle with a SHA-256 build receipt.

The workflow also runs on the merged `main` commit. That merged-main artifact, rather
than the pull-request merge-ref artifact, is the canonical frozen alpha1 build.

After merge, the next database task is taxonomy normalization and constituent GBIF
dataset licence recovery. Those upgrades must be new versioned layers rather than silent
changes to the alpha1 artifact.
