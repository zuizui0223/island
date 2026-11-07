# Conditional sharing of floral profiles between islands

## Intent and scope

Continue the separate floral-profile exploration from commit `506ea1f6792d9bf63cfdd3b6ec0daafddbff2f04` and draft PR #321. Test whether different taxa on different islands carry the same recorded floral profiles more or less often than a genus-conditioned regional reference predicts. The user's instruction to proceed, together with the project instruction to execute available next steps without repeated confirmation, authorizes implementation of this already identified next objective.

The prior within-island richness analysis remains intact. This extension concerns between-island sharing and does not replace, reinterpret, or erase the earlier results. Source inputs remain the same pinned complete-profile mappings, corrected geography, and island flora. Public outputs contain aggregate values only.

## Approaches and implementation decision

1. Compare observed taxon and profile Jaccard values directly. Already done; categorical compression prevents an inference about excess sharing.
2. Independently permute each island or each pair. Rejected because the same taxon could receive contradictory labels and overlapping pairs would be treated as independent observations.
3. Permute profile labels once per regional taxon within genus and reuse that assignment in every sampled pair. Selected by the implementer because it conditions on observed membership and retains dependence among pairs sharing islands or taxa.

Only pairs within the same broad latitude regime are tested. Cross-regime pairs require a different globally coherent conditioning design and are outside this extension.

## Frozen comparison design

- Four evidence scopes: `all`, `direct`, `all_single_state`, `direct_single_state`.
- Four regimes, with the exact keys in `floral_space_null.REGIMES`.
- Compare All and Direct on identical island pairs within each encoding: all recorded state sets versus single-state records.
- Both islands must have at least 50 complete taxon-name records in both scopes of the encoding.
- Remove shared taxon names before converting the remaining names to profiles. Each side must retain at least 20 complete names in both scopes. This is a declared sample-support rule, not a biologically special threshold; record its effect on pair eligibility.
- Sample at most 2,000 qualifying unordered pairs uniformly without replacement per regime and encoding. Use sorted identifiers and deterministic reservoir sampling. Record every candidate and exclusion count and every sampled pair.
- Retain the full regional pool of complete taxa, including taxa on islands below the eligibility threshold and taxa removed from particular pair comparisons. Shared taxa are removed from pair scoring, not from the regional permutation pool.
- Pair selection uses membership and record counts only, never observed profile overlap or the permutation result.
- Primary statistic: median Jaccard similarity of remaining profile sets over the fixed sampled pairs. The statistic weights sampled pairs equally; it is not an equal-island or abundance-weighted estimate.
- Save the observed and permuted medians of profile intersection, union, and profile richness at each pair-exclusive endpoint as descriptive diagnostics. No additional p-values are attached to these components. The null does not condition on each endpoint's observed number of distinct profiles, so concentration and alignment can both contribute to Jaccard changes.
- Generate 999 assignments per scope and regime. Within-genus profile multisets are fixed, including duplicate states. Identity and repeated assignments remain valid draws.
- One assignment per taxon is reused across all pairs in its regime. Regimes are separate conditional comparisons, not reconstructed global species phenotypes.
- Main multiplicity summary: BH across the fixed 16 new pair comparisons, retaining missing comparisons as p=1 in the family. Also report the BY adjustment as a dependency-conservative sensitivity. The prior richness comparisons belong to the explicitly separate exploratory endpoint and remain reported as such.
- Inclusive, finite-corrected doubled smaller-tail p-values reuse the tested implementation in `floral_space_null`. A 95% null range is not an effect confidence interval.
- Seed: 20261106. Pair selection uses `SeedSequence([seed, 17, encoding_index, regime_index])`; permutation streams use `SeedSequence([seed, 29, scope_index, regime_index])`.

## Outputs and interpretation

Write new results to `results/floral_space_pairs_20261106/`:

- `pair_cohort_audit.csv`: matched eligible islands, candidate pairs, exclusions, qualifying pairs, sampled pairs, and sampled endpoint support.
- `pairwise_observed_and_null.csv.gz`: pair identifiers, original/shared/remaining taxon counts, observed profile overlap, and each pair's permutation mean and standard deviation. No per-pair significance tests.
- `permutation_summary.csv`: all 16 comparisons, observed median, null median and quantiles, difference, ratio when defined, p-values, BH and BY.
- `permutation_null_draws.csv.gz`: all regional summary statistics, one row per permutation.
- `manifest.json`: input and code hashes, parent commit, software, parameters, fixed hypotheses, interpretation boundaries, and output hashes.
- A report, verified figure, and validation record documenting results of either direction.

An excess indicates conditional sharing of recorded categorical profiles by nonshared taxa under the stated regional/genus reference. A deficit indicates conditional differentiation. Neither alone identifies evolutionary convergence, selective mechanisms, historical source pools, valid species absences, or pollination-service interchangeability. Profile frequencies are fixed per regional taxon, not per island occurrence. Associations between profile identity and taxon range size can therefore contribute to a result. All/Direct and encoding subsets overlap and are not independent replication.

## Implementation boundaries

Reuse `RegionPermutationDesign`, `build_membership`, `prepare_region`, `draw_profile_assignment`, `finite_two_sided_p`, and `bh_adjust_fixed`. Do not change the prior modules or results. Use Python, NumPy, pandas and SciPy already installed; no new runtime dependency. Keep the original rights exclusions for species-to-profile tables and detailed source audits. Continue the existing dedicated feature branch and update the existing draft PR after verification.

## Acceptance criteria

1. A fully enumerable small example reproduces the exact Jaccard distribution after shared taxa are removed, including duplicated profile labels.
2. Every randomized taxon assignment is shared across the entire pair set and preserves each genus's profile multiset.
3. All/Direct have exactly the same sampled pair identifiers within each encoding and regime; selection is deterministic and independent of input row order.
4. Invalid/self/duplicate pairs, unknown endpoints, mismatched counts or regions, and insufficient remaining taxa fail explicitly. Empty comparisons remain visible without invented statistics.
5. Production p-values and null quantiles can be reconstructed from saved draws; adjusted p-values match an independent implementation; public output contains no taxon-to-trait rows.

## State before inference

This design is written before running or inspecting the new pair-overlap permutation results. An independent method review and a support-only audit are in progress. Any design revisions and their timing will be recorded before production inference.

The method review approved the conditional design and requested the descriptive intersection/union/richness diagnostics and explicit sampled-null degeneracy diagnostics; these were added before implementation. The support audit verified 46 hashes across 24 files and the intended nesting without inspecting profile-sharing outcomes. It counted 835,646 candidate pairs and 775,734 jointly supported pairs. The frozen cap retains 7,618 full-state and 4,562 single-state pairs. Single-state northern high-latitude support is 171 pairs from 40 islands; uneven partner counts and the changed target population are material limitations. The regional assignments are explicitly separate conditional references, despite some taxa occurring in multiple regimes.
