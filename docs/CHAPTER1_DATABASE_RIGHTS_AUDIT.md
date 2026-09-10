# Chapter 1 database redistribution-rights audit

Public release of the exact Database 1.0 species × axis ledger is **fail-closed**. The scientific database identity is already frozen; this document governs only whether the bytes may be redistributed publicly.

## Why `source_groups` cannot be the licensing unit

`source_groups` records acquisition/integration routes such as `latest_public_web` or a validation wave. One such group can contain lineages from several unrelated providers with different reuse terms. It is therefore useful provenance but not a rights decision.

The release builder now audits the exact `source_lineages` column and normalizes row-level identifiers to a provider/source family. It writes both:

- `SOURCE_LICENSE_INVENTORY.csv` — family-level rights decisions and affected resolved cells;
- `SOURCE_LINEAGE_DETAILS.csv` — exact row-level lineage tokens mapped to those families.

A public Zenodo bundle is emitted only if every source family represented in resolved cells is explicitly marked `redistributable` and has a recorded source license.

## Seeded decisions

The policy starts with only decisions that have a clear authoritative basis.

| source family | current decision | recorded license | reason |
| --- | --- | --- | --- |
| `dataset:dryad` | redistributable | CC0-1.0 | Dryad publishes accepted research datasets under CC0 and permits reuse/redistribution. |
| `dataset:austraits` | redistributable | CC-BY-4.0 | AusTraits states that compiled releases are distributed under CC BY 4.0. |
| `database:pladias` | review required | — | Pladias states that provider-specific conditions apply and third-party provision can require provider/Governing Board permission. |
| `domain:en.wikipedia.org` | review required | CC-BY-SA-4.0 | Reuse is possible but release-level attribution/share-alike compatibility has not yet been implemented. |
| `derived:validated_low_without_direct_source_lineage` | review required | — | The stored derived lineage is not sufficient to assign upstream redistribution rights automatically. |
| `unresolved:*` | review required | — | Exact rights cannot be inferred safely from the stored lineage. |

Authoritative pages used to seed the audit:

- Dryad End User Terms: https://datadryad.org/terms
- Dryad reuse guide: https://datadryad.org/help/guides/reuse
- AusTraits project access page: https://austraits.org/
- Pladias data-use conditions: https://pladias.cz/en/download/features

These decisions affect redistribution only. They do not add, remove, or reclassify scientific trait evidence.

## What remains before public Zenodo release

The remaining source families must be audited by provider/domain or resolved to a more precise upstream source. In particular, generic DOI/citation hashes and derived Validated-Low lineages must not be promoted merely because their source was reachable on the public web.

If a source family cannot be redistributed, there are three honest release choices:

1. obtain explicit redistribution permission;
2. publish an open subset plus a machine-readable reconstruction/provenance manifest while keeping Database 1.0's full SHA identity documented;
3. publish a derived aggregate analysis surface that does not reproduce restricted source rows, clearly labelled as different from the exact species-axis Database 1.0.

Option 2 or 3 must never be described as the exact full Database 1.0.
