# Zenodo note: Database 1.0 public subset

The full immutable Chapter 1 Database 1.0 is the scientific analysis database but is not yet wholly redistributable. The rights-filtered public subset contains only cells authorized by the canonical post-PR172 cell-level rights audit.

Deposit the artifact produced by `.github/workflows/build-chapter1-database-public-subset.yml` as a Zenodo Dataset only after that workflow passes its pinned invariants.

The Zenodo description must state that:

1. this is a rights-filtered derivative of Database 1.0, not the canonical analysis database;
2. 46,274 of 222,688 resolved cells are included under the frozen post-PR172 audit;
3. omitted cells are omitted for redistribution-rights reasons and must not be interpreted as biological missingness;
4. the public compilation uses CC BY-SA 4.0 while retaining all upstream source-specific attribution and licence obligations;
5. the paper's analyses remain pinned to immutable Database 1.0, not to the public subset.
