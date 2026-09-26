> **HISTORICAL DATABASE-PROVENANCE NOTE.** This document records the trait-database snapshot and release mechanics used during Chapter 1 development. Database provenance remains valid, but database manifests do **not** select the current scientific analysis. The current paper selector is `config/chapter1_submission_current.json`.

# Chapter 1 trait-database release and provenance

## 1. Database 1.0 snapshot

The immutable trait snapshot used in the current paper lineage is identified by:

- manifest: `config/chapter1_database_versions/v1.0.0.yml`;
- legacy database pointer: `config/chapter1_database_versions/current.yml`;
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

The manifest and hashes identify the database bytes. They do not define H1–H4 or choose a paper model.

## 2. Database manifest versus scientific selector

Two separate objects must remain distinct:

```text
trait database identity
    -> config/chapter1_database_versions/*.yml

current paper analysis identity
    -> config/chapter1_submission_current.json
```

The legacy database manifest can still be validated with:

```bash
python -m island_v2.chapter1_database_manifest validate \
  --manifest config/chapter1_database_versions/current.yml
```

Historical dispatcher code may remain for replay, but it must not be described as the current or canonical Chapter 1 analysis route. The former progressive H1–H5 contracts and their runners are superseded scientific provenance.

## 3. Future Database 2.0+

A future trait-database release should be versioned without rewriting old manifests:

1. produce a new immutable species-axis ledger;
2. create a new version manifest rather than editing `v1.0.0.yml`;
3. lock the new file hash and dimensions;
4. validate source rights and provenance;
5. evaluate any new scientific analysis under an explicitly declared analysis selector/contract.

A data update must never silently redefine the scientific hypotheses after inspecting the new outcomes.

## 4. Zenodo: dataset and software are separate objects

The trait database should be deposited as a Zenodo **Dataset** record. Repository code can separately be archived as a Zenodo **Software** record.

This gives two citable objects:

- **data DOI** — exact trait-database snapshot;
- **software DOI** — code version.

A later Database 2.0 should be a new dataset version, not an overwrite of Database 1.0.

## 5. Build a release candidate

The retained release builder is:

```bash
python -m island_v2.chapter1_database_release \
  --source-dir /path/to/source-scale-batch-integration \
  --manifest config/chapter1_database_versions/v1.0.0.yml \
  --source-policy config/chapter1_database_release_source_policy.yml \
  --output-dir chapter1-database-v1.0.0
```

Public mode remains fail-closed on redistribution rights. Every provenance token represented in resolved cells must have an explicit release decision.

## 6. Reproducibility identity

A published result should be reconstructable from:

```text
database manifest
+ database SHA-256
+ scientific analysis selector/contract
+ code commit
+ result artifact
```

For the present Chapter 1 submission, the scientific selector is `config/chapter1_submission_current.json`.
