# Reproductive rule wave 4 checkpoint (2026-08-14)

This checkpoint uses one original controlled-pollination article that compared
autogamous, geitonogamous and allogamous hand pollination and evaluated fruit,
seed viability and germination. Six target-master species are explicitly
reported as self-compatible and are retained as species-direct High evidence.

All six rows retain the same DOI source lineage. They are six species records,
but one publication lineage. `Orchis patens subsp. brevicornis` is not promoted
to the species-level master record because the result is infraspecific.
Self-compatibility is not converted into autonomous selfing or mating system.

The current `Ophrys × self_incompatibility` queue has two agreeing direct
species and 52 unresolved congeners. The new exact evidence includes two more
`Ophrys` species, but 52 is only a pre-rebuild ceiling. The shared all-evidence
workflow must still apply lineage deduplication, High > Medium > Validated Low
precedence, dominance and both leave-one-species and leave-one-lineage checks.

Rebuild with:

```text
python -m island_v2.reproductive_rule_wave4_checkpoint \
  --master-csv data/v2/staging/gbif/collected/island_taxa.csv \
  --prior-curated-evidence-csv data/v2/staging/traits/open_web_pilot/visible_morphology_rule_checkpoint_20260813/combined_curated_evidence_20260813.csv \
  --prior-curated-audit-csv data/v2/staging/traits/open_web_pilot/visible_morphology_rule_checkpoint_20260813/combined_curated_manual_audit_20260813.csv \
  --output-dir data/v2/staging/traits/open_web_pilot/reproductive_rule_wave4_checkpoint_20260814
```
