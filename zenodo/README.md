# Zenodo release surfaces

Use two records rather than one mixed record.

## Dataset record

Deposit the output of `chapter1_database_public_subset` as a Zenodo **Dataset**. This record receives the database DOI cited by the paper.

The canonical scientific Database 1.0 remains immutable and may contain cells that are not authorized for redistribution. Do **not** upload the full `species_axis_coverage.csv.gz` while its full-database release gate is false. Instead, publish only the rights-filtered derivative produced from the canonical post-PR172 cell-rights audit. The public derivative must satisfy the pinned invariants in the release workflow and must include `PUBLIC_SUBSET_MANIFEST.json` documenting the immutable source database SHA and rights-audit receipt.

Database 1.0 public subset version 1.0.0 is distributed as a compilation under **CC BY-SA 4.0**, while row-level provenance and upstream source-specific attribution/licence obligations remain preserved.

A later Database 2.0 should be a new version of the dataset record so Database 1.0 remains permanently citable.

## Software record

The GitHub repository may separately be enabled in Zenodo's GitHub integration. Repository releases are the software citation surface. `CITATION.cff` supplies repository citation metadata.

The dataset DOI and software DOI should be cross-linked after both records exist.
