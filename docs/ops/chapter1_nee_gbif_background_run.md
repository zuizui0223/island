# Chapter 1 NEE GBIF background runbook

The four new background campaigns are acquisition-only. They do not classify N1 and they do not inspect plant outcomes.

Execution order:

1. Rebuild the frozen global island geometry.
2. Prepare bounded GBIF request blocks for the frozen acquisition taxon.
3. Submit and poll all blocks.
4. Reassign every returned coordinate to the original exact island polygon.
5. Remove cross-block duplicate GBIF IDs.
6. Apply the frozen analytical-background filter (required for Hymenoptera -> seven bee families).
7. Persist the GBIF download ledger, exact-island collection receipt and filtered-background receipt.

The resulting occurrence rows become input to `chapter1_nee_channel_observation.py`; target-catalog detection and effort classification happen only after this acquisition stage closes.
