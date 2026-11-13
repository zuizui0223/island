# Fixed-margin follow-up design (distinct from genus-conditioned PR #321 test)

## Question
Do matched island floras share floral profiles unusually often **after fixing**
(a) the number of profiles observed on each island and (b) each profile's
occurrence frequency among the sampled islands in that region?

This is a **complementary** null model. PR #321 excludes shared taxon names and
preserves observed taxon incidence and genus-level trait frequencies; its
permutations do **not** preserve endpoint profile richness. The public aggregate
data cannot reconstruct every pair's taxon-excluded incidence profile set.
Consequently the new analysis scores complete island profile sets **including
shared taxon names**, and the margins are fixed for that estimand.

## Fixed inputs and statistical unit
Read the public four-scope island/profile incidence files
(results/floral_space_20261106/<scope>/island_profile_carriers.csv.gz),
validate against island_metrics.csv, and reuse the same pair cohorts from
results/floral_space_pairs_20261106/pairwise_observed_and_null.csv.gz.
All and Direct must have identical pair IDs within their encoding and regime.
The statistic is median Jaccard similarity over the **same sampled pair IDs**
(each pair is weighted equally); pairs are not independent replicates.

## Conditional reference
Apply the standard symmetric Curveball switch to a binary island × profile
matrix within each of four regimes and four scopes; an entire sampled island
is assigned one profile set reused across every pair. The trade exactly fixes
each island's profile richness (row sum) and each profile's number of occupied
sampled islands (column sum), and the Markov chain targets the uniform
fixed-margin binary-matrix ensemble. Run two chains with independent seeds,
12 trades-per-island burn-in, three trades-per-island between samples, and
200 draws per chain (400 total).

Report observed minus null median, empirical null 2.5–97.5% range, chain Rhat,
lag-1 autocorrelation and number of distinct sampled statistics. Report
exploratory approximate two-sided tail frequencies **only** when
Rhat ≤1.05, maximum lag-1 ≤0.30 and at least eight distinct null values.
Calculate BH/BY over the fixed family of 16, treating missing p as 1 for
adjustment while keeping its reported value undefined. MCMC draws are not
independent uniform random samples; report p as approximate, not exact.

## Interpretation
This design controls categorical compression opportunities via exact observed
richness and profile prevalence. It does *not* fix genus frequencies, taxon
composition/range size, native status, evolutionary history, or historical
source pools; it can reflect lineage and floristic dissimilarity. Furthermore,
because it retains shared taxa, its Jaccard result cannot be substituted for
PR #321's shared-name-exclusive endpoint. The zero-shared-name pair subset is
retained as a descriptive sensitivity, not a new family of significance tests.

Do not report either rejection or non-rejection as evidence for/against
evolutionary convergence or pollination service replacement. Use only the
aggregate data permitted by the database rights audit, and do not commit
species-to-profile rows. The study is outside the Chapter 1 submission.
