# Chapter 1 database release and version-swap contract

This document separates three objects that must not be conflated:

1. the **versioned Chapter 1 trait database**;
2. the **frozen Chapter 1 analysis contract**;
3. the **paper outputs produced from one database version under that contract**.

The design goal is simple: a future Database 2.0 must be able to replace Database 1.0 without changing the scientific model definitions after outcome inspection.

## 1. Database 1.0 used by the paper

Canonical version manifest:

- `config/chapter1_database_versions/v1.0.0.yml`
- active pointer: `config/chapter1_database_versions/current.yml`

Frozen source:

- workflow run: `34191508045`;
- artifact: `source-scale-batch-integration-34191508045`;
- artifact ZIP digest: `sha256:5e9c887e093c42eb720f3abd2dc85701ce123fb7fe888c916dc4098ac898ccdd`;
- analysis ledger: `species_axis_coverage.csv.gz`;
- ledger SHA-256: `a6eef8d731b2730a99e388ddd0683d7ac9b2af30d43f1f9569c857f89f664b2a`.

Database dimensions:

- 106,295 species;
- 3 axes;
- 318,885 species-axis rows;
- 222,688 resolved cells;
- flower colour: 82,556;
- floral structural complexity: 91,635;
- reproductive assurance: 48,497.

The schema consumed by the paper is:

- `accepted_species`;
- `axis`;
- `trait_composition`;
- `trait_names`;
- `source_groups`;
- `source_lineages`;
- `quality`.

The database manifest locks both the content hash and the expected dimensions. A file with the right filename but different bytes is therefore rejected.

## 2. One entry point for Database 1.0, 2.0, ...

The analysis contract remains:

- `config/chapter1_progressive_analysis.yml`;
- contract ID: `chapter1_progressive_analysis_v1`.

The database version is selected independently through a small manifest. Validate the selected version with:

```bash
python -m island_v2.chapter1_database_manifest validate \
  --manifest config/chapter1_database_versions/current.yml
```

Dispatch the unchanged canonical analysis with:

```bash
python -m island_v2.chapter1_database_manifest dispatch \
  --manifest config/chapter1_database_versions/current.yml \
  --execute
```

The same operation is exposed in GitHub Actions as:

- `.github/workflows/dispatch-chapter1-from-database-version.yml`.

The dispatcher does not duplicate the scientific pipeline. It only translates the selected database manifest into the inputs of:

- `.github/workflows/run-chapter1-progressive-trait-analysis.yml`.

Therefore a database update changes the data pointer, not H1-H5, model order, support gates, or claim ceilings.

## 3. How to create Database 2.0

Do not edit the old manifest. Instead:

1. produce a new immutable source-scale integration artifact;
2. copy `v1.0.0.yml` to `v2.0.0.yml`;
3. change only the database version, run ID, artifact name, file hash, coverage counts and optional `previous_snapshot` pointer;
4. validate the new ledger;
5. move `current.yml` to the new manifest content;
6. dispatch the same analysis contract.

Conceptually:

```text
new trait sources
      |
      v
immutable species-axis ledger
      |
      v
Database v2.0 manifest
      |
      +------------------------------+
                                     v
                      chapter1_progressive_analysis_v1
                                     |
                                     v
                         regenerated paper outputs
```

A new database is allowed to change biological estimates. It is not allowed to silently change the analysis contract in response to those estimates.

## 4. Zenodo: dataset DOI and software DOI are separate

The trait database should be deposited as a Zenodo **Dataset** record. The repository code can separately be archived as a Zenodo **Software** record through the GitHub integration.

This separation gives the paper two citable objects:

- **data DOI** — the exact Database 1.0 snapshot used by the paper;
- **software DOI** — the analysis/data-construction code version.

Database 2.0 should become a new Zenodo dataset version rather than overwriting Database 1.0. After publication, add the version DOI / record URL to the corresponding database manifest; do not alter the content hash.

## 5. Build the Zenodo candidate

Download the source artifact and run:

```bash
python -m island_v2.chapter1_database_release \
  --source-dir /path/to/source-scale-batch-integration \
  --manifest config/chapter1_database_versions/v1.0.0.yml \
  --source-policy config/chapter1_database_release_source_policy.yml \
  --output-dir chapter1-database-v1.0.0
```

This produces a candidate containing the exact analysis ledger, direct ledger, integration summary, baseline recovery manifest, checksums in the release manifest, a data dictionary, and a source-license inventory.

For a public release, add `--public`. Public mode is **fail-closed**: every provenance token represented in resolved cells must be explicitly marked `redistributable` in the source-policy file. Unknown or mixed provenance blocks publication rather than being silently assigned a global license.

The same process is available as:

- `.github/workflows/build-chapter1-database-release.yml`.

## 6. Why the public gate is necessary

The final ledger combines many source groups and derived evidence layers. A global license cannot safely be inferred from the fact that a source was publicly accessible. The release policy therefore records redistribution status source by source.

The candidate bundle may be built before this audit is complete, but the `--public` gate must remain closed until the redistribution inventory is reviewed. This is a distribution/licensing gate only; it does not change which database was used in the scientific analysis.

## 7. Zenodo upload sequence after the gate passes

1. Run the release workflow with `public_release=true`.
2. Download the `chapter1-database-release-candidate-*` artifact.
3. Create a Zenodo record with upload type **Dataset**.
4. Upload the bundle files and set version `1.0.0`.
5. Add creators, description, keywords, funding/acknowledgements and the audited license(s).
6. Reserve or publish the DOI.
7. Record the DOI and record URL in `v1.0.0.yml` and `current.yml` without changing the database SHA-256.

GitHub-to-Zenodo integration can be enabled separately for repository releases so the analysis code receives its own software DOI.

## 8. Reproducibility rule

A paper result should always be reconstructable from the tuple:

```text
(database version manifest,
 database file SHA-256,
 analysis contract ID,
 code commit,
 analysis output artifact)
```

The database DOI is a citation surface; the manifest and hash are the execution identity.
