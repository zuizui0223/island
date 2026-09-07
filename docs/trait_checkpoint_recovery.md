# Recoverable trait checkpoint

The historical private TRY + Wave55 Batch 5 value (222,759 cells) remains recorded, but its completed row-level ledger was not found in the inspected GitHub branches/artifacts or accessible local files on 2026-09-07. Do not use that number as proof of a reproducible current ledger, and do not infer that its 384-cell difference from the public checkpoint is false or lost biological evidence.

The public Wave53 checkpoint can be restored now:

- Source Run: 33605098044; artifact: `wave53-support-two-reproductive-unlock-v2-33605098044`; artifact ID: 9836760588.
- Universe: 106,295 species, 318,885 species-axis cells; no duplicate keys.
- Filled: 222,375 (69.7352%); unresolved: 96,510.
- Structural complexity: 91,613; colour: 82,502; reproductive assurance: 48,260.
- High / Medium / Low: 55,582 / 47,037 / 119,756.
- Species with 0 / 1 / 2 / 3 filled axes: 9,431 / 14,130 / 39,957 / 42,777.

These are restored historical results, not new acquisition. Hash verification proves identity with the pinned public files; it does not independently re-adjudicate each biological assertion or original-source independence. Preserve existing counterevidence holds, including the Wave58 source-independence review, when deciding new promotions.

## Restore without acquiring pages again

```sh
gh run download 33605098044 --repo zuizui0223/island --name wave53-support-two-reproductive-unlock-v2-33605098044 --dir source-wave53
python scripts/recover_verified_trait_checkpoint.py --source source-wave53 --output recovered-public-checkpoint
```

The output path must not exist. The CLI checks the pinned coverage SHA-256, source-file hashes, fixed denominator, key uniqueness and quality tiers, and saves all public source files alongside `recovery_manifest.json`. The `Recover verified public trait checkpoint` workflow performs the same operation and publishes one artifact. It performs no web acquisition and handles no private TRY rows.

## Private reconstruction boundary

The original TRY request was recovered locally and matched against the exact public species universe. It prepared 9,385 source candidate rows for 6,056 species, not 9,385 new filled cells. It remains private. Reconstruction must preserve original-source conflicts and must compare species-trait and species-axis keys before claiming gain. Wave55's public 14-row review queue is not the missing Batch 5 ledger.

The historical checkpoint, recovered public checkpoint, private diagnostic reconstruction, and unpromoted Wave56-58 candidates must remain separately labelled. Replace the historical working claim with a newer checkpoint only after the new row-level ledger, direct/Low split, source provenance, gains and losses are validated. Never silently fill the 384-cell difference or automatically rerun previously completed acquisition tasks.
