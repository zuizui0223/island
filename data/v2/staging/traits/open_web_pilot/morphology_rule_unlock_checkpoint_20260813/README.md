# Morphology support-2 rule-unlock checkpoint (2026-08-13)

This checkpoint contains nine individually reviewed species-direct records
selected from genus × trait rules that previously had two agreeing supporting
species. The original species page or official PDF was re-read for every row;
the exact quote, URL, retrieval hash, source lineage, quality tier, reviewer,
and review time are frozen in the evidence and audit tables.

The checkpoint excludes `Gonystylus macrophyllus × flower_size_class` because
the same PNGTrees treatment is already present in the PNGTrees source package.
It also records, but does not promote, the tempting `Sloanea dasycarpa ×
floral_symmetry` hit because the statement is a family description rather than
a species description.

The 106 queued cells are only a pre-rebuild theoretical potential. Formal
coverage is determined only after the common all-evidence implementation
rebuilds genus × trait rules, lineage leave-out validation, conflicts, and
quality precedence. No genus, family, or global inference is emitted here.

Rebuild from the repository root:

```powershell
$env:PYTHONPATH = "src"
py -3.13 -m island_v2.morphology_rule_unlock_checkpoint `
  --master-csv data/v2/staging/gbif/collected/island_taxa.csv `
  --prior-curated-evidence-csv data/v2/staging/traits/open_web_pilot/kudo_2022_alpine_reproduction_checkpoint_20260813/combined_curated_evidence_20260813.csv `
  --prior-curated-audit-csv data/v2/staging/traits/open_web_pilot/kudo_2022_alpine_reproduction_checkpoint_20260813/combined_curated_manual_audit_20260813.csv `
  --output-dir data/v2/staging/traits/open_web_pilot/morphology_rule_unlock_checkpoint_20260813
```
