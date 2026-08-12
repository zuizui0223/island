# Combined reviewed source package through IOSPE (2026-08-12)

This package deterministically merges the previously reviewed Go Botany,
BASEFLOR, PROTEUS, and FlorML/size evidence with the IOSPE orchid flower-colour
checkpoint. It re-runs the repository-wide trait-specific review gate over all
1,200 audit rows.

The merged audit has 1,192 correct rows (precision 0.9933), zero cultivar
contamination, and nine approved traits. Of 2,513 candidate evidence rows, 2,505
remain after eight explicitly reviewed failures are removed. The formal workflow
then performs source-lineage deduplication, direct conflict resolution, and
`genus x trait_name` Validated Low reconstruction against the prior ledger.
