# Current own-DB syndrome research workspace

**Owner: island.** Focal island flora is the project's GBIF-derived incidence roster, not official GIFT island checklists. The historical label `GIFT4505` is retained only for exact replay; official GIFT mainland donor data are a distinct auxiliary source.

## Verified input chain

- Target: protected island artifact 10058653212, `fixed/canonical/input/chapter1_status_flora.csv.gz` (4,505 islands); GBIF occurrence-derived roster with appended status.
- Traits: same archive, `fixed/source/current/evidence/direct_only_trait_ledger.csv.gz` (own integrated trait ledger).
- Donor: `fixed/source/source/gift_native_mainland_flora.csv.gz` (748 GIFT mainland candidate regions); not the target flora.
- Geography: pinned corrected island geography (8,264 units).
- Thus every Q1–Q4, pair atlas, origin-size and colour-size sensitivity using these inputs is island-owned, including mixed donor analyses. No independent GIFT replication is implied.

## Migration checkpoints

- Original fd PR #69 snapshot: `85769a811d151b05940b554b27dc1a662c28b1df` ([archive](../../snapshots/pr069.json)).
- Latest inspected fd PR #69 head: `3097c8553dc96fe17a9e2d60653e89d5bb6b34ff` ([archive](../../snapshots/pr069_3097c855.json)); three changed files copied byte-for-byte, including the anonymous origin-size summary and colour-taxon sensitivity script.
- Runner is preserved in `docs/legacy_workflows/`, not activated. Historical names and source bytes are kept, not silently re-labelled as new results. Future execution and development belong to island. No protected raw data copied, no new numerical result claimed.

The frozen GEB submission selector is unchanged.
