# IOSPE orchid flower-colour checkpoint (2026-08-12)

This checkpoint recovers species-direct flower-colour text from the Internet
Orchid Species Photo Encyclopedia (IOSPE). IOSPE is treated as a specialist
Web source and therefore as a Medium candidate, not as an official database or
an image-derived source.

The acquisition matched 5,644 exact accepted-species headings in the island
master, selected the 3,768 cells whose strict flower-colour axis was unresolved,
and fetched 3,655 source pages. The 113 failed pages remain unresolved. Only a
source-page description with an exact species heading and a locally attached
flower, petal, sepal, tepal, lip, spur, or pouch colour can enter the candidate
ledger. Common names, image captions, snippets, cultivar or hybrid text,
comparators, vegetative colours, family inference, and global fallback are
excluded.

The parser was developed against two disjoint 200-page samples. The production
gate uses a third deterministic, disjoint 200-page holdout. Direct review found
193 correct rows (precision 0.965 under `accepted_correct / all_reviewed`) and
zero cultivar contamination. Seven known failures are recorded in the audit and
excluded by the shared source-package gate. The approved output has 661 evidence
rows covering 631 species-trait cells before source-lineage and prior-ledger
deduplication.

Reproduction information:

- `iospe_index_fetch_manifest_20260812.csv` records the 25 source indexes and
  their content hashes.
- `iospe_exact_master_index_20260812.csv.gz` records every exact island-master
  match; `iospe_unresolved_colour_candidates_20260812.csv.gz` records the
  acquisition queue.
- `iospe_full_page_fetch_manifest_20260812.csv.gz` records every requested URL,
  HTTP outcome, byte count, and content hash. Cache paths are acquisition-time
  relative paths and are not evidence locators.
- The three reviewed audit files retain development, intermediate holdout, and
  final holdout decisions with exact quotes, reviewer, time, reason, lineage,
  and content fingerprint. `iospe_prior_reviewed_400_for_exclusion.csv`
  deterministically excludes the two development samples from final sampling.
- `iospe_orchid_colour_evidence_20260812.csv.gz` is the candidate source package;
  the shared review gate excludes the seven final-holdout failures before formal
  promotion.

The immutable source quote, canonical URL, retrieval date, page SHA-256,
normalized-excerpt fingerprint, and source lineage are retained for every row.
