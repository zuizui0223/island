# Zenodo release surfaces

Use two records rather than one mixed record.

## Dataset record

The Chapter 1 Database 1.0 rights-filtered public derivative has been published as a Zenodo **Dataset**:

- DOI: `10.5281/zenodo.22704973`
- Record: `https://zenodo.org/record/22704973`
- Public derivative: 46,274 of 222,688 resolved species × axis cells
- Compilation licence: CC BY-SA 4.0

The canonical scientific Database 1.0 remains immutable and contains cells that are not authorized for redistribution. The Zenodo DOI therefore identifies the rights-filtered public derivative, not redistribution of the full `species_axis_coverage.csv.gz` analysis database. Omitted cells are omitted for rights reasons and must not be interpreted as biological missingness.

The published package is derived from the canonical merged-main public artifact and retains `PUBLIC_SUBSET_MANIFEST.json`, the immutable source-database SHA, the post-PR172 rights-audit receipt, public cell-level rights rows, source-rights policy, and release notes.

A later Database 2.0 should be a new version of the dataset record so Database 1.0 remains permanently citable.

## Software record

The GitHub repository may separately be enabled in Zenodo's GitHub integration. Repository releases are the software citation surface. `CITATION.cff` supplies repository citation metadata.

The dataset DOI and software DOI should be cross-linked after both records exist.
