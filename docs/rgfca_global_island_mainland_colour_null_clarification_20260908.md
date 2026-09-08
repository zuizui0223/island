# Pre-outcome clarification: RGFCA × island null statistic

Frozen 2026-09-08 (JST), before any island-versus-mainland colour contrast was computed.

This note clarifies one scalar comparison left implicit in `rgfca_global_island_mainland_colour_protocol_20260908.md` without changing eligibility, outcomes, sampling depth, repetitions, seeds, directions, or sensitivities.

For each primary outcome, the observed scalar used against the 1,000-label-permutation null distribution is the **median of the 1,000 species-equal empirical resampling estimates** generated with seed `20260908`.

- whitening p-value: `(1 + count(null >= observed_median)) / 1001`;
- chromatic-dulling p-value: `(1 + count(null <= observed_median)) / 1001`.

The null realizations use seed `20260909` and follow the already frozen within-species 5-versus-5 label-permutation rule. No further inferential clarification may be made after outcome opening except to correct a demonstrated implementation error while preserving the frozen estimand.
