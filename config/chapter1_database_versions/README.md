# Chapter 1 database versions

Each YAML file identifies one immutable species × trait-axis database snapshot.

- `v1.0.0.yml` is the database used for the frozen Chapter 1 paper analysis.
- `current.yml` is the active execution pointer and initially has identical content.

When Database 2.0 is ready, add `v2.0.0.yml`; never rewrite `v1.0.0.yml`. Update `current.yml` only after the new ledger passes schema/hash validation.

Changing the database manifest must not change `config/chapter1_progressive_analysis.yml`.
