# Chapter 1 database releases

Start with [`docs/CHAPTER1_DATABASE_RELEASE.md`](docs/CHAPTER1_DATABASE_RELEASE.md).

## Current role

The database layer preserves **data identity, validation, rights auditing and release packaging**. It does not select or dispatch the current scientific analysis.

- Database 1.0 manifest: `config/chapter1_database_versions/v1.0.0.yml`
- Legacy database pointer: `config/chapter1_database_versions/current.yml`
- Manifest validator: `python -m island_v2.chapter1_database_manifest validate`
- Historical database-dispatch workflow path: `.github/workflows/dispatch-chapter1-from-database-version.yml` — validation only; scientific analysis dispatch is retired
- Zenodo candidate builder: `.github/workflows/build-chapter1-database-release.yml`
- Redistribution-rights audit: `docs/CHAPTER1_DATABASE_RIGHTS_AUDIT.md`

The **current Chapter 1 scientific selector** is `config/chapter1_submission_current.json`.

The `analysis_contract: chapter1_progressive_analysis_v1` field retained inside Database 1.0 manifests is historical provenance. It must not be interpreted as an active model selector.
