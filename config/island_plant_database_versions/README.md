# Island Plant Database version receipts

This directory pins immutable build receipts for the reusable Island Plant Database.

Each version manifest records the scientific scope, rights boundary, canonical GitHub
Actions run/artifact receipt, and SHA-256 digest for every materialized table.

`current.yml` is the exact active receipt for the materialized **island x taxon core**.
It remains `2.0.0-alpha1` until a later database version actually replaces or extends
that core bundle.

Taxonomy normalization is a non-destructive layer over the frozen alpha1 core.
`current_taxonomy.yml` is therefore the exact active receipt for that separate taxonomy
crosswalk layer. It currently points to `2.0.0-alpha2-taxonomy` and does not imply that
alpha1 island-flora membership, establishment status, or trait values were mutated.

A version receipt does not imply that every table is publicly redistributable. Release
status is recorded separately from scientific inclusion. In particular,
`2.0.0-alpha1` retains all candidate GBIF-derived island x taxon rows for database
construction while keeping their value redistribution `review_required` until
constituent GBIF dataset licence lineage is recovered.

Published releases must add their persistent DOI/record URL to the versioned manifest
without changing the pinned scientific table hashes.
