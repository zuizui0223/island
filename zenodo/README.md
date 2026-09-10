# Zenodo release surfaces

Use two records rather than one mixed record.

## Dataset record

Deposit the output of `chapter1_database_release` as a Zenodo **Dataset**. This record receives the database DOI cited by the paper. Do not publish while `RELEASE_MANIFEST.json` reports `release_ready_for_public_zenodo: false`.

A later Database 2.0 should be a new version of the dataset record so Database 1.0 remains permanently citable.

## Software record

The GitHub repository may separately be enabled in Zenodo's GitHub integration. Repository releases are the software citation surface. `CITATION.cff` supplies repository citation metadata.

The dataset DOI and software DOI should be cross-linked after both records exist.
