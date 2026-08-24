# Reproductive counterevidence checkpoint (2026-08-13)

This checkpoint follows the support=2 reproductive queue from integrated Run
`31699299351`. It adds nine individually reviewed species-direct records from
primary papers, an official agency assessment, a doctoral thesis and a primary
monograph. Positive candidates and counterexamples are frozen in the same
batch so that an attractive genus rule cannot be unlocked by omitting contrary
evidence.

The clean expected rule candidate is `Abutilon × autonomous_selfing_capacity`:
the two current agreeing species plus *A. theophrasti* may resolve up to 55
currently unresolved cells, but this is a queue ceiling rather than claimed
coverage. The Allium, Calanthe, Liparis and Grevillea records deliberately add
counterevidence; they should reduce or block unsafe genus-wide inference. The
shared all-evidence implementation must determine the actual result after
High > Medium > Validated Low precedence, direct-conflict resolution,
source-lineage deduplication, dominance, species leave-one-out and lineage
leave-one-out validation.

Each row keeps the exact `species × trait_name` key. Self-compatibility,
autonomous selfing, mating system and cleistogamy are not interchangeable.
This checkpoint emits no genus inference, family inference, global fallback or
`n=2` formal result.

Rebuild with:

```text
python -m island_v2.reproductive_counterevidence_checkpoint \
  --master-csv data/v2/staging/gbif/collected/island_taxa.csv \
  --prior-curated-evidence-csv data/v2/staging/traits/open_web_pilot/cached_evidence_recovery_checkpoint_20260813/combined_curated_evidence_20260813.csv \
  --prior-curated-audit-csv data/v2/staging/traits/open_web_pilot/cached_evidence_recovery_checkpoint_20260813/combined_curated_manual_audit_20260813.csv \
  --output-dir data/v2/staging/traits/open_web_pilot/reproductive_counterevidence_checkpoint_20260813
```
