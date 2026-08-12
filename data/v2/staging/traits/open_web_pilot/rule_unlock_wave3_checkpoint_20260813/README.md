# Rule-unlock wave 3 checkpoint (2026-08-13)

This checkpoint records seven individually reviewed species-direct trait
statements selected from the formal `current_support == 2` acquisition queue,
plus one direct upgrade of a Validated Low cell. It does not emit genus
inference itself.

The source pages or documents were retrieved before review. Every row stores an
exact supporting excerpt, content SHA-256, accepted-name identity, source
lineage, retrieval time, and a one-to-one audit decision. The shared
all-evidence workflow later deduplicates lineages, resolves conflicts and tests
`genus x trait_name` rules.

Pre-formal queue targets are deliberately not reported as achieved coverage:

- `Bidens x autonomous_selfing_capacity`: 47 unresolved congeners
- `Coccothrinax x flower_size_class` and `inflorescence_display`: 50 distinct
  unresolved species-axis cells in total, not 100
- `Bambusa x inflorescence_display`: 38 unresolved congeners
- `Bulbophyllum tseanum x self_incompatibility`: one Low-to-direct candidate

Thus the maximum distinct strict coverage target is 135 species-axis cells.
The exact realized increment must come from the formal GitHub Actions artifact.

Rebuild from the preceding checkpoint:

```powershell
$env:PYTHONPATH = "src"
py -3.13 -m island_v2.rule_unlock_wave3_checkpoint `
  --master-csv data/v2/staging/gbif/collected/island_taxa.csv `
  --prior-curated-evidence-csv data/v2/staging/traits/open_web_pilot/rule_unlock_wave2_checkpoint_20260812/combined_curated_evidence_20260812.csv `
  --prior-curated-audit-csv data/v2/staging/traits/open_web_pilot/rule_unlock_wave2_checkpoint_20260812/combined_curated_manual_audit_20260812.csv `
  --output-dir data/v2/staging/traits/open_web_pilot/rule_unlock_wave3_checkpoint_20260813
```

No family inference, global fallback, formal `n=2` rule, search snippet, or
cross-trait substitution is present in this checkpoint.
