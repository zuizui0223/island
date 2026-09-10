# Chapter 1 database versions

Each YAML file identifies one immutable species × trait-axis database snapshot.

- `v1.0.0.yml` is the database used for the frozen Chapter 1 paper analysis.
- `current.yml` is the active execution pointer and initially has identical content.

When Database 2.0 is ready, add `v2.0.0.yml`; never rewrite `v1.0.0.yml`. Update `current.yml` only after the new ledger passes schema/hash validation.

Changing the database manifest must not change `config/chapter1_progressive_analysis.yml`.

## Rights-aware provenance for Database 2.0+

Database 1.0 predates the release-rights contract, so its source rights are being reconstructed conservatively after the fact. New database versions should instead ship a SHA-locked `rights_registry` alongside the species-axis ledger.

The registry contains one row per exact `source_lineage` and must record:

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

Validation then checks both the registry hash and that every source lineage used by every resolved species-axis cell is represented. A lineage missing from the registry is a hard failure. Rows marked `redistributable` must have both an explicit license and rights evidence.

This changes publication provenance, not the scientific H1-H5 analysis contract.
