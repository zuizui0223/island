# Active v2 build status

## Current milestone

Rebuild GBIF island species-name acquisition before trait collection.

## Implemented

- exact-polygon standardization for a documented island vector source;
- per-island GBIF asynchronous `SPECIES_LIST` request manifests;
- no bounding-box substitution, no facet truncation, and no conversion of API failures to empty floras;
- download-key, status, DOI, and raw-name provenance retention;
- GitHub Actions validation plus a manually triggered GBIF workflow;
- tests covering multipart geometries, area filtering, and exact geometry preservation.

## External input still required before the first global GBIF submission

The original vector source used by v1 (`GlbIslands.gdb`, layer `BigIslands`) is not committed to this repository. It must be placed at `data/v2/external/islands/GlbIslands.gdb` or substituted with a documented equivalent source, then the manual workflow stage `prepare` can run.

The `submit` workflow stage also requires repository/environment secrets `GBIF_USERNAME` and `GBIF_PASSWORD`; `GBIF_EMAIL` is used for notification-address metadata.

See `docs/gbif_flora_v2.md` for the exact command sequence and interpretation boundary.
