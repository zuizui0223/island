# Chapter 1 database releases

Start with [`docs/CHAPTER1_DATABASE_RELEASE.md`](docs/CHAPTER1_DATABASE_RELEASE.md).

- Active database pointer: `config/chapter1_database_versions/current.yml`
- Versioned Database 1.0 manifest: `config/chapter1_database_versions/v1.0.0.yml`
- One-click analysis dispatcher: `.github/workflows/dispatch-chapter1-from-database-version.yml`
- Zenodo candidate builder: `.github/workflows/build-chapter1-database-release.yml`
- Redistribution-rights audit: `docs/CHAPTER1_DATABASE_RIGHTS_AUDIT.md`

The scientific analysis contract remains `chapter1_progressive_analysis_v1`; database versions are replaceable inputs, not new model definitions. Public redistribution is a separate, fail-closed rights gate.
