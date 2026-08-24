# Reproductive rule wave 5 checkpoint (2026-08-14)

This checkpoint adds seven exact species x reproductive-trait records selected
from support=2 genus-rule acquisition targets. Sources are primary studies,
official government treatments, or curated botanical database treatments. Each
row preserves the exact quote, URL, content fingerprint, and one canonical
source lineage.

The records keep self-incompatibility, autonomous selfing capacity, and mating
system separate. In particular, self-compatibility is not converted to
autonomous selfing, and rare successful selfing in *Fritillaria meleagris* is
retained alongside the authors' `predominantly_outcrossing` classification.

The queue estimated 463 potentially unlockable cells. That is a prioritization
ceiling, not coverage. Only the shared all-evidence workflow may rebuild
trait-specific genus rules and accept Low cells after dominance, masked
leave-one-species-out, and leave-one-source-lineage-out validation.

Rebuild with:

```text
python -m island_v2.reproductive_rule_wave5_checkpoint \
  --master-csv data/v2/staging/gbif/collected/island_taxa.csv \
  --prior-curated-evidence-csv data/v2/staging/traits/open_web_pilot/reproductive_rule_wave4_checkpoint_20260814/combined_curated_evidence_20260814.csv \
  --prior-curated-audit-csv data/v2/staging/traits/open_web_pilot/reproductive_rule_wave4_checkpoint_20260814/combined_curated_manual_audit_20260814.csv \
  --output-dir data/v2/staging/traits/open_web_pilot/reproductive_rule_wave5_checkpoint_20260814
```
