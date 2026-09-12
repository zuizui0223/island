# Island Plant Database 2.0-alpha2 taxonomy crosswalk

This layer normalizes the 115,328 provisional alpha1 plant-name identities against the current Catalogue of Life Extended Release through GBIF's v2 species matcher.

It is a non-destructive crosswalk layer. Database 2.0-alpha1 remains frozen and unchanged.

The full build is sharded into 12 deterministic segments and combined only after source SHA, segment continuity, checklist identity, metadata identity, crosswalk SHA, row count, and unique submitted-name checks pass.

The combined artifact contains:

- `taxonomy_crosswalk.csv.gz`
- `taxonomy_review_queue.csv.gz`
- `gbif_colxr_metadata.json`
- `TAXONOMY_MANIFEST.json`

An automatic-resolution candidate is not equivalent to human curation or a flora-membership decision. Taxonomy normalization does not decide whether a plant is native, introduced, established, endemic, or absent from an island.
