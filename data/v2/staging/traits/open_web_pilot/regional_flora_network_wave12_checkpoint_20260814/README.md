# Regional-flora network wave 12 checkpoint (2026-08-14)

This checkpoint records the first production acquisition from the shared
Zambia, Mozambique, Malawi and Botswana flora interface.  GitHub Actions Run
`31767694551`, artifact `9207019490`
(`regional-flora-network-acquisition-31767694551`, artifact digest
`sha256:ff556cb3744fad5a81e39d5858d5cb511759dbef05dcd4101f3ae585f7c0ae64`)
scanned family IDs `[1,320)` on all four sites.

The run retrieved 887 family indexes successfully, inventoried 21,110 source
records, matched 7,163 records exactly to the fixed 106,295-species strict
denominator, and fetched 735 prioritized species pages.  It emitted 37
species-direct candidates for 32 species.  All 37 were reviewed against their
stored exact excerpt.  Thirty-one were accepted and six were rejected:

- three component dimensions incorrectly proposed as whole-flower size;
- one fruit dimension and one between-flower spacing dimension proposed as
  flower size;
- one cyathial gland/appendage colour proposed as primary flower colour.

Thus the full extraction audit precision is `31/37 = 0.8378378`, below the
production-wide approval threshold.  No provider or trait scope is approved
wholesale.  Only the 31 individually reviewed rows enter the combined curated
ledger.  The review audit records every rejection and its reason; the shared
extractor was tightened so these error classes fail closed in later runs.

Accepted rows are all Medium species-direct evidence with exact accepted-name
and family matches, original URLs, exact excerpts, retrieval timestamps,
content SHA-256 values and provider-treatment source lineages.  They comprise
14 flower-colour, eight inflorescence-display, seven flower-size, one
floral-form and one floral-symmetry row.  `pollen_vector_mode` and
`reward_type` remain separate traits and are absent from this checkpoint.

Formal coverage is not inferred from this directory.  The common all-evidence
rebuild must apply High > Medium > Validated deterministic Low precedence and
recompute genus rules by `genus x trait_name`.  Family inference, global
fallback and `n=2` formal inference remain excluded.

## Rebuild

```bash
python -m island_v2.regional_flora_wave11_checkpoint \
  --source-snapshot-csv data/v2/staging/traits/open_web_pilot/regional_flora_network_wave12_checkpoint_20260814/regional_flora_network_wave12_source_rows_20260814.csv \
  --master-csv data/v2/staging/gbif/collected/island_taxa.csv \
  --prior-curated-evidence-csv data/v2/staging/traits/open_web_pilot/regional_flora_wave11_checkpoint_20260814/combined_curated_evidence_20260814.csv \
  --prior-curated-audit-csv data/v2/staging/traits/open_web_pilot/regional_flora_wave11_checkpoint_20260814/combined_curated_manual_audit_20260814.csv \
  --output-dir data/v2/staging/traits/open_web_pilot/regional_flora_network_wave12_checkpoint_20260814 \
  --source-group regional_flora_network_wave12_checkpoint_20260814 \
  --output-stem regional_flora_network_wave12 \
  --created-at 2026-08-14T04:05:00Z \
  --reviewer "Codex exact-excerpt regional-flora audit" \
  --contract regional_flora_network_wave12_checkpoint_v1 \
  --acquisition-manifest-json data/v2/staging/traits/open_web_pilot/regional_flora_network_wave12_checkpoint_20260814/source_acquisition_manifest_31767694551.json
```
