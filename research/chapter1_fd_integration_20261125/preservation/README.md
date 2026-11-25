# Chapter 1 source-preservation checkpoint

This records actual GitHub metadata at 26 November 2026 JST, not a newly preserved copy. Source restrictions and the original H1–H4 freeze remain intact.

| Source | Artifact ID | Expiry (UTC) | Role |
|---|---:|---|---|
| source-scale-batch-integration-34191508045 | 10042377593 | 2026-12-07T05:40:13Z | current lineage |
| chapter1-progressive-analysis-34232450884 | 10058653212 | 2026-12-07T13:31:06Z | current lineage |
| chapter1-database-v1-rights-audit-34482428564 | 10154261753 | 2026-12-09T13:24:19Z | current lineage |
| chapter1-database-v1-cell-rights-34566271397 | 10186293306 | 2026-12-10T05:31:10Z | current lineage |
| island-plant-database-2.0.0-alpha1-34665438815 | 10288721725 | 2026-12-11T01:41:36Z | developmental only |

The immutable v1.0 ledger expires first. Preserve the four current-lineage archives together; alpha1 is optional and does not replace the paper database. The original web-provenance, WCVP mapping/range and H4 parent artifacts listed in the manifest are not returned by the API. Missing listings are not proof of irrevocable loss. The original corrected H4 derived outputs and frozen result tables still exist in Git. Scan every surviving archive for exact matches to the two missing H1 input hashes; never substitute a nearby WCVP version.

The canonical release audit (run 34566271397) counted 46,274 redistributable, 176,413 rights-review-required and 1 semantic-review-required cells. Review-required is not equivalent to prohibited. 5,032 blocker families; unresolved citation lineages 505. These are the historical audit results, not new permissions or new science. The public subset DOI is not the full input. Do not publish restricted sources, source-level rows or the original ZIP.

## Actual custodian backup (authorized existing storage only)

From a checkout of this branch, with authenticated GitHub CLI and Python 3, run:

```sh
python research/chapter1_fd_integration_20261125/preservation/preserve_sources.py --destination /path/to/authorized/existing/storage --ack-authorized-storage --receipt /path/to/authorized/existing/storage/source_receipt.json
```

The destination must already exist, be access-controlled and authorized by applicable licences. Optional `--include-alpha1` retains the separate developmental database. The script verifies outer bytes, archive CRCs and exact nested pins, refuses mismatching existing files and does not upload or modify GitHub. No destination is inferred. A verify-only runner uses temporary storage and deletes the source bytes; its receipt cannot be represented as a backup.

[Metadata and exact pins](SOURCE_PRESERVATION_MANIFEST.json). The user/custodian still needs to run the backup command or specify an approved storage location accessible for that operation. No external account, publication or permanent archive has been authorized by this commit.

## First verification failure and strict diagnostic

Run 45660302290 reached the repository-pin check after the ZIP and inner-pin verification loop, then failed at the H1 YAML SHA. The file's Git blob `6173ab3432cd7fbf0d51b1fab32e11bab797ac2d` is identical in its introducing commit `a7b19298dc0bf5ff00aff651913037ac1dfaff81`, the frozen main, and this branch. The follow-up tests **only** exact LF-to-CRLF reconstruction against the historical SHA256, with the original Git blob additionally pinned. No new scientific settings or accepted hash substitution. An unresolved mismatch remains a failed gate; diagnostics are retained even on failure. No assertion of CRLF equivalence until the real check runs.
