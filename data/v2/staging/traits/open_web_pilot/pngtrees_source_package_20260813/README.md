# PNGTreesKey flower-trait source package (2026-08-13)

This credential-free package materializes structured species descriptions from
the *Guide to Trees of Papua New Guinea* by Barry J. Conn and Kipiro Q. Damas.
The guide is hosted by the National Herbarium of New South Wales and Papua New
Guinea National Herbarium.

## Scope and evidence contract

- The acquisition started from the guide's complete tree-description index,
  not from search snippets.
- Only exact species-rank page headings already present in the fixed 106,295
  species master were considered.
- The family printed on each page had to agree exactly with the master. Eleven
  outdated or conflicting family assignments remain visible in the page
  manifest as `family_rejected` and contribute no evidence.
- The accepted values come only from controlled fields in the page's
  `Flowers:` paragraph: flower symmetry, the guide's explicit 10-mm diameter
  class, inner-perianth colour, and flower arrangement.
- Multiple stated values remain state sets. Ficus syconium shape and
  arrangement are not treated as flower symmetry or a solitary flower.
- All pages share one compilation lineage. They are not counted as independent
  sources in genus-rule lineage validation.
- Every accepted row is species-direct Medium. The package does not create
  genus, family, or global inference.

The credential-free acquisition made 313 HTTP requests (one index and 312
species pages), used zero search-API queries, and cost USD 0. The frozen result
contains 312 fetched exact-rank target pages, of which 301 passed identity and
family gates, and 1,143 species-trait evidence rows. A
stable hash-stratified audit covers 200 unique pages (50 per trait). The shared
source-package gate reports precision 1.0 and cultivar contamination 0.0 for
the four approved traits. These are extraction-audit statistics for this
structured source, not a claim that the source itself is biologically complete.

## Rebuild

From the repository root:

```powershell
$env:PYTHONPATH = "src"
py -3.13 analysis/acquire_pngtrees_source_package_20260813.py --skip-fetch
py -3.13 -m pytest -q tests/test_pngtrees_source_package.py
```

Omit `--skip-fetch` to reacquire the live index and pages. The page manifest
records each source URL, exact flower excerpt, retrieval time, response SHA-256,
identity decision, and rejection reason. The source-package manifest pins all
committed input and output hashes.
