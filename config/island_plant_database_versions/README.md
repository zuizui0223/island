# Island Plant Database version receipts

This directory pins immutable build receipts for the reusable Island Plant Database.

Each version manifest records the scientific scope, rights boundary, canonical GitHub
Actions run/artifact receipt, and SHA-256 digest for every materialized table. The
`current.yml` file is an exact copy of the version manifest currently designated as the
active database version.

A version receipt does not imply that every table is publicly redistributable. Release
status is recorded separately from scientific inclusion. In particular,
`2.0.0-alpha1` retains all candidate GBIF-derived island x taxon rows for database
construction while keeping their value redistribution `review_required` until
constituent GBIF dataset licence lineage is recovered.

Published releases must add their persistent DOI/record URL to the versioned manifest
without changing the pinned scientific table hashes.
