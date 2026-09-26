> **Historical/provenance pointer.** These manifests identify immutable trait-database snapshots. They do not select the current Chapter 1 scientific analysis; use `config/chapter1_submission_current.json` for the corrected paper surface.

# Chapter 1 database versions

Each YAML file identifies one immutable species × trait-axis database snapshot.

- `v1.0.0.yml` identifies the frozen Database 1.0 snapshot.
- `current.yml` is retained as the legacy database-version pointer used by database tooling. It is **not** the current paper-analysis pointer.

Never rewrite an old version manifest. A future Database 2.0 should be added as a new version with its own hashes and dimensions.

Changing a database pointer must not silently reactivate or redefine the superseded progressive H1–H5 analysis contracts.

## Rights-aware provenance for Database 2.0+

New database versions should ship a SHA-locked `rights_registry` alongside the species-axis ledger.

The registry should record, for every exact `source_lineage`:

- normalized source family;
- provider;
- source URL;
- dataset DOI where available;
- license;
- redistribution status (`redistributable`, `review_required`, or `not_redistributable`);
- evidence supporting the rights decision.

A future manifest can declare:

```yaml
rights_registry:
  relative_path: source_lineage_registry.csv.gz
  sha256: <64-hex SHA256>
  required_for_public_release: true
```

Validation must fail if a lineage used by a resolved species-axis cell is missing from the registry. Rows marked `redistributable` require both an explicit license and rights evidence.

This is a data-provenance/release rule, not a scientific H1–H4 model definition.
