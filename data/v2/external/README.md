# External data retention

Large external inputs and GBIF download archives are not silently committed by the pipeline.

- The island geometry source should be archived under a documented release or DOI and referenced in `island_manifest.csv`.
- GBIF request manifests, download keys, DOI values, and checksums are versioned as provenance.
- Raw downloads may be retained in a release/archive location when too large for the main Git repository.

The v2 workflow never overwrites a previous manifest or treats a failed GBIF request as zero species.
