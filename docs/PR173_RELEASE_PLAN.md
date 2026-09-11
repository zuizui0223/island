# PR173 release plan

Purpose: close the release-layer handoff after PR #172 without altering scientific Database 1.0.

Acceptance gates:

- canonical Database 1.0 SHA remains unchanged;
- canonical post-PR172 rights artifact is pinned by run id, artifact id, name, and digest;
- all 222,688 resolved cells are checked against the frozen rights ledger;
- exactly 46,274 `redistributable` cells are exported;
- exactly 176,414 blocked/review cells are excluded;
- the single semantic mismatch is excluded;
- public subset contains only high/medium/low resolved cells;
- no trait values, quality tiers, H1-H5 contract, or provenance tokens are rewritten;
- output manifest records both immutable scientific database SHA and rights-audit receipt;
- analysis remains pinned to Database 1.0 rather than the public derivative.

After merge, the next scientific step is a canonical H1-H5 rerun from immutable Database 1.0; further source-rights expansion becomes release maintenance rather than the Chapter 1 main line.
