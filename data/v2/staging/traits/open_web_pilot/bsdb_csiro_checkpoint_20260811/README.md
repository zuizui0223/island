# BSdb + CSIRO reviewed source checkpoint (2026-08-11)

This checkpoint combines two independently acquired, source-lineage-preserving
packages without treating the wrapper as an independent biological source:

- Zell et al. (2025) BSdb (`10.1111/nph.20234`), pinned to
  `dirtyplants/BSdb@9e87946d1e3121d39e657b702cf9f92ccc10936e`;
- CSIRO Australian Tropical Rainforest Plants, online edition 8, whose current
  species index was retrieved from
  `https://apps.lucidcentral.org/rainforest/text/entities/index.htm?v=20260811`.

The combined evidence file has 1,240 rows. BSdb contributes 574 explicit
species-level `SC`/`SI` rows retaining 284 original-study lineages. CSIRO
contributes 666 species-direct flower-section rows (527 species) after exact
accepted-name and family gates and exclusion of previously resolved direct
species-traits. No snippet, family inference, global fallback, or axis-level
trait substitution is used.

The deterministic audit contains 400 rows: 200 BSdb structured records and
200 CSIRO page excerpts. Every row records the exact value, quote or structured
record, URL, reviewer, time, and decision. All 400 passed the extraction,
identity, provenance, and cultivar gates (precision 1.0; cultivar contamination
0.0). This validates extraction from the cited records; it is not a new
biological remeasurement of the upstream studies.

CSIRO acquisition is resumable: downloaded treatments are stored in a local
gzip cache, and a rerun requested zero completed pages. The committed index,
treatment inventory, summary, package manifest, and SHA-256 manifest allow the
selection and outputs to be independently checked without committing raw page
copies.

Rebuild CSIRO candidates with:

```bash
python -m island_v2.lucid_rainforest_checkpoint \
  --master-csv data/v2/staging/gbif/collected/island_taxa.csv \
  --strict-coverage-csv strict_species_axis_coverage.csv.gz \
  --baseline-direct-csv direct_species_trait_ledger.csv.gz \
  --cache-dir csiro-rainforest-cache \
  --output-dir csiro-rainforest-output
```

Formal promotion must use `bsdb_csiro_evidence_20260811.csv.gz` together with
`bsdb_csiro_manual_audit_400_20260811.csv` through the common reviewed source
package gate and PR #131's `genus x trait_name` Validated Low rebuild.
