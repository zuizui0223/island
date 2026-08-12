# Palmpedia / Palmweb archive checkpoint (2026-08-13)

This package is a credential-free, source-backed Medium pilot for unresolved
Arecaceae traits. It does not use search snippets as evidence.

## Acquisition

- Discovery used the Internet Archive CDX API for archived
  `palmpedia.net/wiki/*` pages. The inventory contained 14,761 archived wiki
  URLs.
- The 500 highest-priority exact master-name matches were retrieved as archived
  full pages. Page-title identity accepted 498 pages and rejected two.
- Every evidence row records the original URL, archive URL and timestamp, CDX
  digest, retrieval time, response SHA-256, normalized content fingerprint,
  exact quote, and normalized-quote fingerprint.
- Full archived HTML is not committed. Exact evidence quotes and immutable
  archive/content identifiers are committed, so each decision can be replayed
  without redistributing page bodies.

## Review and production gate

All 243 extracted candidates were reviewed line by line. The reproducible page
audit contains all 154 candidate pages plus 46 no-candidate pages sampled with
seed `20260813`, for 200 manually reviewed pages total.

The common source-package gate is trait-specific:

| Trait | Reviewed | Correct | Precision | Production |
| --- | ---: | ---: | ---: | --- |
| `flower_primary_color` | 88 | 86 | 0.9773 | approved |
| `flower_size_class` | 99 | 99 | 1.0000 | approved |
| `inflorescence_display` | 46 | 46 | 1.0000 | approved |
| `floral_form` | 10 | 8 | 0.8000 | rejected |

The four rejected rows concern calyx shape, anther colour, or fruit-stage
persistent perianth rather than the claimed strict trait. Cultivar
contamination was zero in the reviewed sample. Approved rows remain
species-direct Medium. `pollen_vector_mode` and `reward_type` are not mapped to
the flower-structure axis.

Rebuild with:

```powershell
$env:PYTHONPATH = "src"
python analysis/prepare_palmpedia_archive_source_package_20260813.py
```

The manifest records input/output hashes and the exact review gate result.
