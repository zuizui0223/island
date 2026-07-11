# LLM recovery pilot, 2026-07-12

This directory records the first real execution of the evidence-locked LLM recovery lane after PR #86.

The fixed-site collector produced 13 baseline candidate rows across 7 of 25 validation species. Three additional source-backed LLM claims passed the repository validator, increasing any-trait coverage to 9 of 25 species. The result demonstrates that the validator and provenance contract work, but also shows that source discovery is the dominant bottleneck.

The exact workflow run and artifact identifiers are stored in `llm_recovery_pilot_20260712_summary.json`. The retained CSV contains the three validated incremental candidate rows used in the benchmark.

These rows remain unreviewed candidates and are not promoted to curated biological truth.
