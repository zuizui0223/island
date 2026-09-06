# Wave38 source-locked extraction audit

Wave38 revisits the already-frozen extraction packets under
`data/v2/staging/traits/llm_evidence_extracted`.  A model response is never
treated as a biological source.  Promotion requires all of the following:

- an exact quote present in the frozen source packet;
- a species-level identity match already recorded by the packet;
- fail-closed mapping to one specific trait ontology;
- source-lineage and normalized-excerpt deduplication;
- a deterministic, domain-by-trait audit sample;
- at least 90% `accepted_correct / reviewed` and at most 2% cultivar
  contamination for the promoted domain-by-trait scope.

`wave38_manual_review_decisions.csv` is the frozen audit decision table.  The
GitHub workflow rebuilds the candidate audit, reviewed direct ledger,
trait-specific `genus x trait_name` rules, and three-axis overlay from this
file and the source packets.  The formal Wave33 input is pinned to GitHub Run
`32932103226`, artifact `integrated-trait-coverage-32932103226`, artifact ID
`9593698966`, and its immutable digest.  Wave37 Run `33143109604` supplies the
immediately preceding lossless coverage ledger.

The strict overlay excludes pollen-vector and reward evidence from the flower
structure axis, does not use family inference or global fallback, and refuses
any High/Medium downgrade or previously filled species-axis loss.
