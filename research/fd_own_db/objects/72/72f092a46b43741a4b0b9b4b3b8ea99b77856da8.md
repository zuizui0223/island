# fd — functional differentiation of island floral assemblages

**Status:** independent exploratory analysis line, migrated from [island PR #321](https://github.com/zuizui0223/island/pull/321). This repository contains new analyses, **not** a new copy of the original plant-trait database. Its first runnable stage is implemented; no new global significance or evolutionary result is implied by creating the pipeline.

## Research question

**Can island floras occupy locally constrained combinations of recorded floral traits while remaining functionally differentiated from other islands?**

The target is *floral phenotype assemblage structure*, not demonstrated pollination effectiveness. The present two features are the recorded **flower-colour state set** and **floral-form state set**. A multi-state description stays one profile; it is never expanded into hypothetical colour × form flowers.

The main argument is tested in stages: (1) local spread, (2) between-island functional separation, (3) source/lineage/area-conditioning, then (4) ecological interpretation. Hypotheses and nonclaims are frozen in [docs/ANALYSIS_CONTRACT.md](docs/ANALYSIS_CONTRACT.md).

## What already exists upstream (not new fd findings)

From the frozen public [PR #321 report](https://github.com/zuizui0223/island/blob/d4a5dafab5f0737ea6c1c1c0263fc20decf6035c/results/floral_space_20261106/README.md):

- **496** distinct All profiles (467 Direct) from recorded colour × form state sets.
- Northern midlatitude **single-state** island profile richness is ~10.8% (All) and 11.4% (Direct) below genus-conditioned permutation medians (both BH q=0.016); this does **not** hold as an adjusted result for multistate-inclusive scopes.
- Genus-conditioned *shared-name-exclusive* between-island profile similarity has no BH/BY-supported effect in 16 comparisons.
- Whole-island **fixed-margins** profile Jaccard is lower than its approximate MCMC null in 16/16 comparisons, with BH 13/16 and BY 12/16 supported; this is a different endpoint **including shared taxa** and not conditioned on genus.
- An observed profile-turnover share of 0.690–0.902 (median across analyses) is *descriptive*, not above-null evidence.
- Gaussian KDE hypervolumes were **not** calculated. Two PCoA axes retain only 15.35% of the common positive inertia; 80% requires 17 axes. Do not call a 2-D KDE area total floral hypervolume.

## Fixed absolute-distance reference (N1d, 2026-11-09)

[Verified 1,000/2,000 km source-pinned N1d report](results/20261109_fixed_radius_ses/README.md) · [Methods and limits](docs/FIXED_RADIUS_DISTANCE_SES.md) · [Original-distance implementation](src/fd/fixed_radius_ses.py) · [Successful full CI #42587303351](https://github.com/zuizui0223/fd/actions/runs/42587303351)

**Equalizing absolute donor distance also attenuates, but does not erase, the oceanic flower-category proximity pattern.** In original-D nearest-*other*-profile SES after uniform 30-profile thinning (All), a **1,000 km** reference gives median SES **−0.037 continental (127 islands) / −0.227 oceanic (115)**; a **2,000 km** reference gives **+0.001 continental (129) / −0.370 oceanic (124)**. On exactly the same oceanic islands, the median change from their own N1b broad-regional reference to fixed-distance N1d is **+0.295 SES at 1,000 km** and **+0.227 at 2,000 km**. The remaining oceanic median observed-minus-null D departures are small: **−0.0050 / −0.0088** original-distance units. Direct source-evidence shows the same descriptive direction.

**Fixed radius is not fixed observed source density.** Within 1,000 km the median number of other eligible donor islands is **41 continental versus 9 oceanic**; within 2,000 km **123 versus 24.5**. Inadequate donor support removes some islands, so comparisons explicitly use same-island subsets; All and Direct overlap. These radii were chosen *after* seeing N1c's geographic mismatch and are exploratory, **not** new geological-significance tests or evidence of selection. Taxon/genus sharing and detection still remain uncontrolled. The next necessary independent-data gate is source- and genus-conditioned heldout N2 rather than adding more geographic reference variants.

## Nearby recorded-island reference (N1c): geographic attenuation, not geological proof

[Validated same-regime nearest-ten source report](results/20261109_local_geography_ses/README.md) · [Method/limits](docs/LOCAL_GEOGRAPHIC_NULL.md) · [Original-D implementation](src/fd/local_geography_ses.py)

**The public-source N1c run succeeded** ([CI #42584462875](https://github.com/zuizui0223/fd/actions/runs/42584462875)). With 30 distinct profiles per island, All oceanic median nearest-profile D SES changes from **−0.555** against broad-regional island-prevalence N1b to **−0.231** against the ten geographically nearest *other* eligible islands' profile-prevalence N1c; the continental medians are **+0.166 → −0.051** (130 islands per origin). Direct repeats the attenuation (oceanic **−0.543 → −0.258**). All/Direct are nested and these comparisons are **descriptive**: no new geological-origin p-value is inferred. The oceanic tenth-neighbor distance is **~1,149 km** median versus **~387 km** for continental islands, so fixed ten-neighbor geography has very different absolute reach. This baseline still allows neighboring islands to share the same taxa/genera and cannot replace donor-held-out N2. **Verified paired same-island audit:** N1c−N1b nearest-profile SES difference has a median **+0.286 across the same 130 All oceanic islands** (continental −0.165); Direct oceanic **+0.306** (n=119). In the oceanic group, the median nearest-10 reference mean itself is **0.00731 original-D units lower** than the regional reference for the same islands. This confirms that the weaker negative oceanic SES is not caused by switching which islands are analyzed, but does not attribute its cause.

## 2026-11-10: direct marginal-fixed test shows little geographic colour/form coupling prediction

[Full source-verified result and limitations](results/20261110_joint_flower_source/MARGINAL_FIXED.md) ·
[12-row aggregate proper scores and whole-block intervals](results/20261110_joint_flower_source/marginal_fixed_placebo_aggregate.csv) ·
[Five heldout geographic fold audit](results/20261110_joint_flower_source/marginal_fixed_five_fold_audit.csv) ·
[Reproducible IPF model](src/fd/marginal_fixed_coupling.py) ·
[Four-figure manuscript v0.5](docs/FD_MANUSCRIPT_DRAFT_V0_5.md)

**The separate floral traits explain most of the apparent joint-profile source signal.** In exactly the same **1,007 All islands (119 blocks) / 950 Direct (116 blocks)** and the same two independently donor-coverage-/pool-size-matched alternative source assignments, the original full-rank proper **joint colour×form kernel** favors real mainland source geography by **+0.002378 All (95% block interval +0.000878 to +0.003753)** and **+0.001097 Direct (−0.000153 to +0.002178)**. We then use **iterative proportional fitting** to force each actual/placebo model's **separate colour-set and form-set predicted marginals to equal the same heldout global-genus marginals** on the original 496 observed profiles. The **marginal-fixed joint source-placebo** gain drops to just **+0.000214 All (+0.000036 to +0.000370)** and **+0.000140 Direct (−0.000062 to +0.000331)**; matched full-profile Brier gives **+0.000205 All** and **+0.000137 Direct**, respectively. Original unadjusted additive/joint gains were reproduced on the same islands and there were **zero IPF convergence failures** (tolerance 2×10⁻⁷).

This ~**91% All / 87% Direct score attenuation** is a *reference-conditioning sensitivity*, **not a causal fraction of phenotype evolution explained**. All retains a tiny conditional association signal while nested Direct's conservative confidence interval still includes zero. The paper should emphasize strong known-genus prediction and weak, donor-evidence-sensitive mainland geographic source information that is predominantly **marginal trait-linked**. No adaptive joint floral syndrome, pollinator-mediated selection or mainland ancestral origin has been identified.

**Further conservative audit:** The five complete heldout geographic test-folds all show positive marginal-fixed joint-kernel point estimates in both scopes, but the **32 possible full-fold sign configurations give a two-sided tail of 0.0625** if folds were exchangeable; they are not fully independent because donor models and taxa overlap. Dropping any one fold leaves an All mean of +0.000183–+0.000249 and Direct +0.000101–+0.000175. The paired **form minus colour** source-gain difference also remains uncertain: All **+0.001738 [−0.000221,+0.003288]**, Direct **+0.000662 [−0.000808,+0.001967]**. Neither the tiny coupling residual nor apparent form dominance should be described as a universal selection mechanism. [Reproducible aggregate-only fivefold audit](src/fd/fold_coupling_sensitivity.py) · [nonidentifying output](results/20261110_joint_flower_source/fold_coarse_sign_sensitivity.csv) · [paired marginal comparison](results/20261110_joint_flower_source/form_vs_colour_paired_source_gap.csv).

## 2026-11-10 original floral traits: additive scores do not identify joint colour/form arrangements

[Marginal versus joint trait-kernel methods](docs/FLORAL_TRAIT_COUPLING.md) · [Executable and source-audited kernel algebra](src/fd/trait_kernels.py) · [Synthetic coupling falsification](tests/test_trait_kernels.py) · [Four-figure manuscript decision draft](docs/FD_MANUSCRIPT_DRAFT_V0_5.md)

**Important correction to the meaning of previous N2 results:** the old original-distance similarity score uses **Kadd=(Kcolour+Kform)/2**. Therefore the earlier similarity-aware local native mainland source gain is an exact average of **colour-only** and **form-only** score improvements; it does **not** establish prediction of their recorded *joint arrangements*. The frozen 496-profile public reference has **31 distinct colour-set profiles and 46 form-set profiles**. Direct pivoted-Cholesky audits (rank tolerance 1e-8) yield Kcolour rank **31**, Kform **46**, Kadd **76**, and the **joint Hadamard product Kcolour×Kform rank 496/496**. Thus on the frozen full profile reference the new Kjoint is capable of strictly proper full-profile probability scoring; unlike Kadd it detects different colour-form pairings with identical marginal proportions.

**Joint source prediction now tested on the restored, SHA-verified private original checkpoints.** [Measured trait-kernel source placebo](results/20261110_joint_flower_source/README.md) · [aggregate block-bootstrap table](results/20261110_joint_flower_source/coverage_matched_trait_kernels.csv). With identical source-support-matched controls, **Kjoint** yields **+0.002378 All (1,007 islands, 119 spatial blocks; 95% interval +0.000878 to +0.003753)** and **+0.001097 Direct (950 islands, 116 blocks; −0.000153 to +0.002178)**; the latter interval crosses zero. The original **Kadd** is independently reproduced exactly as the marginal mean (+0.001860/+0.000927), validating the replay. For All strict source controls, form-only gain (+0.002729) is larger than colour-only (+0.000991) as a point estimate. The joint product kernel is strictly proper here (minimum eigenvalue **+0.003895**), but a positive joint gain can still arise from marginal improvements; **it is not evidence of interaction-specific evolution or pollinator adaptation**. Restricted species/source identities remain local, and public CI validates only code plus the public 496-profile reference.

## 2026-11-10 fairness audit: source opportunity matching substantially reduces geographical predictive advantage

[Five-fold true-versus-matched-source report](results/20261110_n2_source_placebo/README.md) · [Aggregate 95% intervals](results/20261110_n2_source_placebo/coverage_matched_summary.csv) · [Reproducible private-input implementation](src/fd/n2_source_assignment_placebo.py)

**The decisive new post-hoc falsification:** In taxon-name-heldout, spatial-block-folded mainland-genus prediction, the earlier real-source-versus-random-source difference in the original distance-based proper K-score (+0.00396 All / +0.00352 Direct) partly reflects unequal independent **genus donor coverage** (+5.37/+3.52 percentage points for real assignments). After matching two same-regime, different-block, non-identical alternative source assignments to within **±2 percentage points of independent local-genus donor support**, and additionally restricting their **total independent donor-name pool size to 0.5–2×**, the actual source advantage is **+0.001860 All (95% spatial-block interval +0.000687 to +0.002877; 1,007 islands/119 blocks)** and **+0.000927 Direct (−0.000031 to +0.001741; 950 islands/116 blocks)**. Brier agrees in direction, but Direct is not clearly separated from zero. The five geographical fold mean K-score improvements are positive in each scope, though folds are not fully independent botanical lineage experiments.

The **much smaller**, Direct-uncertain residual is the current inference ceiling for geographically specific mainland-source information, and cannot be called evolutionary selection, geological floral niche contraction or verified historical provenance. Control source/caliper choices are post-hoc; stricter matching changes the island population and does not allow a causal percentage-explained claim. Genus and geographic opportunity remain predictive, but **source evidence availability materially affects the size and apparent certainty of its incremental advantage**.

## New 2026-11-10 macro result: local mainland-source identity predicts floral similarity, not exact categories

[Full heldout N2 proper scoring report](results/20261110_n2_proper_scoring/README.md) · [Aggregate score/CI table](results/20261110_n2_proper_scoring/island_block_bootstrap.csv) · [One-name-one-vote genus sensitivity](results/20261110_n2_proper_scoring/taxon_name_balanced_genus_bootstrap.csv) · [Code](src/fd/n2_crossfit_macro.py)

**Empirical, source-audited private-checkpoint reconstruction, NOT a geological adaptation test:** five disjoint geographic block folds excluded all evaluated island taxon names from every GIFT mainland donor entity. Across 1,533 All islands (170 spatial blocks) and 1,443 Direct islands (159 blocks), globally known genus profiles greatly improve exact heldout phenotype log loss (+0.641 All / +0.360 Direct relative to global no-genus prior). Adding the five nearest putative *native-mainland genus donors* yields almost **zero additional exact-profile log-loss gain** (+0.00011 All; +0.00019 Direct; both spatial-block 95% intervals cross zero), but improves a **proper original-distance similarity kernel score** by **+0.00829 All [0.00633,0.01004] / +0.00948 Direct [0.00655,0.01191]**. The full-category strictly proper Brier score improves **+0.00221 All [0.00095,0.00335]**, while Direct's interval includes zero. K = 1 - original set-Jaccard D is numerically PSD but rank-76, so its quadratic score is proper but **not strictly proper for all exact categories**.

**Species/name repetition checked:** 91.6% of All and 95.5% of Direct island×taxon rows involve names seen on at least 10 islands. When each different taxon name is counted ONCE across heldout folds, and entire genera are bootstrapped, the source K-score gain remains positive (**+0.01354 All [0.01027,0.01716]; +0.00905 Direct [0.00696,0.01109]**). This is taxon/lineage-sensitivity, not an independent phylogeny or independent mainland ancestry proof.

**Critical caveat:** independently heldout *local genus donor coverage* is only **24.3% All / 19.6% Direct** of name carriers. The **specific geographic assignment** has now been tested against same-regime, different-block placebo assignments with independent genus-donor coverage and total donor count matched. Its additive source advantage is small (**+0.001860 All, +0.000927 Direct**), and Direct's spatial-block interval includes zero; see the [strict source-placebo report](results/20261110_n2_source_placebo/README.md). The source lists describe potential mainland supply and do not identify historical arrivals, within-island selection or pollinator function. The separate strict 37-pair/five-component and 10 geography-disjoint geological-origin N2 comparisons are **not newly confirmed**.

## N2 donor-held-out genus/source distance reference (private-only, 2026-11-09)

[Source/lineage conditional protocol](docs/N2_SOURCE_DISTANCE.md) · [Private-input implementation](src/fd/n2_source_distance.py) · [Synthetic tests](tests/test_n2_source_distance.py)

**The separate fixed-profile N2 geological-origin SES test has not been estimated.** In contrast, heldout N2 *predictive* genus/source comparisons and their donor-matched joint probability scores have been evaluated in the authorized local input workspace. The public island–profile aggregates alone omit the accepted-name/genus identities required for either private model. The new code uses rights-authorized GIFT native mainland donors, removes every evaluated island taxon name from every source and donor region, and generates model-based reference distributions conditional on **observed island genus counts**. The fixed-n weighted support draw is an approximation to a conditional genus null, not an exact genus-fixed permutation. Only one four-slot post-hoc family is planned for 37 pairs (five spatial components) and 10 geographically disjoint verified pairs. Donor support, the model-predicted number of distinct profiles and **true held-out taxon-name prediction accuracy** (exact-profile log loss and original categorical-distance loss versus a global-genus benchmark) must all be audited before interpreting N2. Support coverage, donor holdout, expected richness and sensitivity to smoothing remain inferential gates. **No pollinator selection, adaptation, or geological causation is claimed.**

## 2026-11-09 original-distance floral packing / fixed-support SES

[Validated original-distance SES report](results/20261109_distance_ses/README.md) ·
[Distance-native protocol](docs/DISTANCE_SES.md) ·
[Executable module](src/fd/distance_ses.py) ·
[Synthetic tests](tests/test_distance_ses.py)

**Source-pinned CI completed successfully** ([N1b-complete run #42579672053](https://github.com/zuizui0223/fd/actions/runs/42579672053)). None of the original 32 spatially controlled SES origin contrasts supports a geological-specific effect (min exact p=0.0879; minimum BH q=0.4464).
The 17D Gaussian-KDE HDR result is now considered a **diagnostic appendix**:
when narrow 17D kernels are mostly disjoint, its log volume tracks the number
of distinct recorded profiles and its richness-adjusted paired confidence
intervals all contain zero. Gaussian smoothing is not a natural physical
measure of a finite colour/form category repertoire.

The replacement uses the original 496 × 496 colour/form set-Jaccard distance
matrix, with two distinct-profile measures: **mean off-diagonal pair distance**
and **mean within-island nearest-profile distance**. Smaller values indicate
closer recorded categories, not measured pollinator function. For each island
the metrics are compared against uniform fixed-n profile draws from (N0) the
global 496-profile public reference and (N1) recorded states in other eligible
islands of the same broad region, with the focal island removed. **N1b** is a
more stringent descriptive sensitivity that draws regional profile types in
proportion to their occurrence in OTHER eligible islands; the universal negative
SES appearance under uniform N1 does not survive this frequency weighting
(All equal30 nearest D median SES: continental N1 −1.45 → N1b +0.18;
oceanic N1 −2.21 → N1b −0.52). These are recorded-support nulls, not genus
or true mainland-source conditionals.
Uniform observed 30-profile thinning is an additional support sensitivity,
not independent biological replication.

Geological SES contrasts are re-tested only as **post-hoc spatial
sensitivities**: entire geographic components are sign-flipped for the old
37-pair sample, and pairwise flips are enumerated for the newer geography-
disjoint **10 verified pairs**. There are 32 explicitly enumerated post-hoc
comparisons with one BH/BY family (N1b is descriptive and excluded). Source/genus conditional N2 still requires
rights-restricted underlying donor data and has **not** been implemented in
this public pipeline. No new causal claim or population-level floral
contraction is implied by adding the code.

## Local spatial predictability of floral profile repertoires (2026-11-09)

[44-source-archipelago heldout prediction report](results/20261109_archipelago_flower_transfer/README.md) · [Source/distance validation protocol](docs/ARCHIPELAGO_FLORAL_HOLDOUT.md) · [all 12 comparisons](results/20261109_archipelago_flower_transfer/archipelago_holdout_tests.csv)

**New evidence from outside the original 37-pair geological subset:** The original All/Direct common support includes 711 islands with named Weigelt archipelago labels, of which 606 from 44 named archipelagos can be evaluated by holding out each island and predicting its 496 recorded flower-colour/form profile distribution from three other islands in the same source group. Compared with the regional **outside-group** average, the within-archipelago predictor lowers exact-profile mean error from 0.006861 to 0.005077 (All) / 0.007351 to 0.005473 (Direct), and the original full categorical-kernel error from 0.010609 to 0.006203 (All) / 0.011137 to 0.006574 (Direct). Group-level sign flips give BH q=0.0003/BY q≈0.00093 across 12 correlated post-hoc test slots for these broad improvements.

**But this is almost entirely spatial, NOT a demonstrably special source-archipelago identity effect.** For **581/606 targets (95.9%)**, the three nearest physical islands are *exactly the same three islands* as the three nearest within the target archipelago. A **geography-only nearest-three** predictor, without using archipelago identity to choose donors, achieves virtually identical heldout profile error. The 44-group within-vs-nearest geographic improvement is not supported (All exact p=0.5913, Direct p=0.7033; original-kernel All p=0.8565, Direct p=0.9296). A stricter distance-matched different-archipelago control on 51 targets in 18 source groups likewise does not establish an additional named-source effect.

This supports **local predictability of recorded floral assemblage composition**, not geological trait filtering, evolutionary source ancestry, out-of-species transfer, confirmed functional pollination niches or true floristic absence. Nearby islands may share the exact same plant taxa. The next unresolved ecological question is whether profile similarity persists after **shared species/genera are excluded**, using authorized species-identity-level data only; the public aggregate alone cannot answer it.

## 2026-11-09 geographic source-disjoint rematching: no confirmed geological effect

[Independent-pair floral repertoire report](results/20261109_geography_independent/README.md) · [Source-disjoint matching contract](docs/GEOGRAPHIC_INDEPENDENCE.md) · [Exact test aggregates](results/20261109_geography_independent/independent_origin_repertoire_tests.csv)

The original 37 continental/oceanic island comparisons share only five geographically connected source groups. We therefore **re-matched the entire eligible geological flora population before examining flower traits**, maximizing different physical island comparisons under the original environment, area, isolation, name-support and coverage gates while prohibiting reuse of spatial blocks and/or source archipelago labels among pairs. This gives **13 block-disjoint pairs, 11 source-archipelago-disjoint pairs, 11 satisfying both, and 10 with both constraints plus verified archipelago labels**.

The strictest **10-pair** original categorical-profile exact L2 comparison is **All p=.1680 / Direct p=.2188** under complete `2^10` exact pair flips; across the **full post-hoc 16-test family none passes BH or BY** (smallest BH q=.2734). A nominal Direct p=.0408 in a **different 13-pair block-only subset** is not supported after multiple-testing correction. Original categorical-similarity-kernel comparisons likewise remain unsupported.

**This supersedes any claim of confirmed oceanic-island adaptive floral syndrome or source-independent structural divergence.** Observed category repertoires differ in particular samples, and known-genus/native-mainland opportunity predicts some recorded assemblage structure. But no currently available set of geographically disjoint pairs provides adjusted support for an independent geological-origin effect. More separate source archipelagos and identity-aware detection/phylogenetic controls are needed. Geographic labels are pragmatic proxies, not definitive phylogenetic independent histories.

## Critical 2026-11-09 spatial-independence audit — earlier pairwise significance is NOT confirmatory

[Exact conservative spatial-component analysis](results/20261109_coverage_rarefaction/README.md) · [Component graph and method](docs/COVERAGE_RAREFACTION_SENSITIVITY.md)

**The original 37 geological oceanic/continental matched island pairs are not 37 demonstrably independent geographical tests.** Under the frozen spatial-block graph, linking pairs that share either endpoint's geographic block gives only **5 connected components (23, 11, 1, 1, 1 pairs)**. All 23 northern-midlatitude pairs belong to ONE such connected component; a second coarse Weigelt archipelago grouping independently gives **5 components (18, 13, 4, 1, 1)**. Whole-component exact origin-label flips offer only 32 configurations, with minimum possible p=.0625 for the sign-invariant squared multivariate statistic.

**Actual original-profile exact-L2 results under conservative block grouping:** All exact p **0.2500**, Direct **0.3750** (37 matched pairs); in the earlier 21-pair GBIF-effort-balanced subset, All **0.2500**, Direct **0.3125**. Thus the previous pair-independent exploratory q-values in PR #4 and PR #8 **are not supported as spatially robust discoveries**. Raw category differences still exist, but cannot be claimed to exceed geographically dependent null expectations. Coarse blocks connected by graph transitivity are a deliberately severe stress test, not proof of fully dependent species histories.

The same 37 pairs do retain exact recorded colour×form composition differences after **equal-profile/complete-name rarefaction**, but this does NOT repair spatial pseudoreplication or the stronger **genus×regional-availability conditional null non-rejection**. The primary manuscript framing must shift to **source-lineage assembly and predictive, explicitly observational floral repertoires**, rather than claiming confirmed geological causal filtering or independent oceanic functional niche contraction.

## Sampling-support-standardized geological floral repertoire (2026-11-09)

[Validated 37-pair rarefaction report](results/20261109_coverage_rarefaction/README.md) · [Methods and bounds](docs/COVERAGE_RAREFACTION_SENSITIVITY.md)

**Exact recorded colour×flower-form combinations remain different after controlling the observed amount of information per island.** The frozen 37 paired continental–oceanic island floras all have ≥30 distinct recorded profiles and ≥40 complete taxon-name carriers. Under 99 repeat uniform observed-support draws and 1,999 paired sign swaps per draw:

- Exactly **30 distinct recorded profiles per island**: median exact-profile L2 p **0.007 All / 0.0035 Direct**; 83.8% / 93.9% of the repeated thinnings have raw p below .05.
- Exactly **40 sampled complete name carriers per island**: median p **0.007 All / 0.004 Direct**; 82.8% / 93.9% of thinnings below raw .05.
- Per-pair exact profile richness matched to the poorer side: median p **0.0005** in both nested evidence scopes.
- Restricting to ≥20 **single-state-only profile types** per pair retains **11 All / 10 Direct** pairs: exact-profile p 0.0517 / 0.093. This restricted subset has insufficient support for a comparable positive generality claim.
- The full original categorical **trait-similarity kernel** remains unsupported across these support screens. The distinction is largely exact recorded category-profile membership, **not demonstrated functional niche-space disjunction**.

**These 99 resamples are correlated thinnings of the same 74 original islands, not independent replications or inferred missing species.** This strengthens the evidence that simple observed sample size alone is insufficient, but it does **not** override prior genus×regional-source conditional null non-rejection, remove trait-dependent reporting bias or demonstrate geological adaptation/colonization. [CI #42554306259](https://github.com/zuizui0223/fd/actions/runs/42554306259) ran the full aggregate source pipeline with 50 tests passed.

## Origin-matched floral repertoire screen under independent GBIF observation-effort metadata (2026-11-08)

[Full effort-robustness results](results/20261108_gbif_effort_robustness/README.md) · [Pinned source/method contract](docs/GBIF_EFFORT_SENSITIVITY.md)

Rechecked the **same 37 original oceanic–continental island pairs** using the independent, earlier `island` v2 **Tracheophyta** exact-island GBIF effort metadata (record counts, dataset counts, coordinate uncertainty). All 37 pairs had older effort rows and ≥30 GBIF IDs/≥2 datasets on both islands. Requiring additional **≤4× record-count asymmetry and ≤10 km reported median coordinate uncertainty** retained **21 physical pairs** (12 north-midlat, 9 tropical).

The exact colour×form profile-repertoire separation remains supported in **All** (21-pair signed origin-label swap p=0.0025, 16-slot BH q=0.0057, BY q=0.0193), and more weakly Direct (p=0.0089, BH q=0.0178, **BY q=0.0602**). The original **categorical-similarity-kernel** separation is **not supported** in either scope (All p=0.474; Direct p=0.756).

**Strict interpretation:** this is a **post-hoc gross observation-effort screen** in a *different/earlier sampling campaign*, not a true per-genus detection model or independent replication. It does not establish biological absences, colonization failure, phylogenetic independence, or adaptive flower evolution. It also does **not** override the genus×geographic-availability conditional null non-rejection in [PR #5](results/20261108_genus_source_control/README.md).

## Mainland lineage opportunity × floral assemblage prediction and 17D projected KDE (2026-11-08)

[Corrected source / two-stage / 17D results](results/20261108_source_opportunity_hv17/README.md) · [Methods and privacy contract](docs/SOURCE_OPPORTUNITY_HYPERVOLUME17.md)

**Taxon-name leakage audit changes the interpretation.** Allowing GIFT mainland lists to include the **same species** present in tested island floras makes source–island overlap biologically descriptive, but not an independent out-of-species prediction. After excluding every observed taxon name from all 434 geology-classified islands globally from the mainland source lists, the earlier apparent nearby-source genus ranking improvement on continental islands **disappears/reverses**. The strict source-ranking support is 399 islands (174 continental / 225 oceanic) and 301 in explicitly native records (141/160). Its original 37 geographical pair differences are **not BH/BY-supported** over 12 tests (minimum BH q≈0.487). Unrecorded mainland candidate genera are **unknown, not true ecological absences**.

In the stricter independent two-stage flower-assemblage prediction, all names from the original 74 matched test islands are also excluded from GIFT **genus opportunity source features and floral trait donors**. The source-only genus mixture, without using observed island genus identities, lowers original categorical similarity-kernel **floral frequency** error from 0.020476 to 0.018107 (All, **11.6%**) and 0.024446 to 0.021582 (Direct, **11.7%**). Matched-pair bootstrap improvements remain positive; exact individual-taxon profile log-loss improvements do not consistently remain. The known-island-genus benchmark is an *oracle*, not prospective genus colonization prediction.

The new [17-dimensional KDE HDR implementation](src/fd/hypervolume17.py) goes beyond the previous two axes (15.35% inertia): 17 common PCoA axes retain **80.49%** of the 496-profile positive inertia. Runs use two bandwidths, 12,000 draws ×2 per island, 260 geological All and 243 geological Direct islands. Raw oceanic volume differences are largely accounted for by differences in the number of distinct recorded flower profiles, and no profile-richness-adjusted contrast has a 95% pair-bootstrap interval excluding zero. These are **projection-limited Gaussian-smoothed volumes, not total discrete floral functional hypervolumes or pollination niches**.

**Scientific boundary:** mainland source geography can constrain recorded plant lineage repertoires, but existing observation checklists cannot distinguish actual arrival, establishment, observation, native/introduced history or adaptation causally. Source-only flower-repertoire prediction is a modest assemblage-level association, not an ecological sorting mechanism proof. No restricted GIFT taxon identities/source rows or original island evidence are republished.

## Out-of-taxon mainland source and genus prediction (2026-11-08)

[Held-out predictive results](results/20261108_heldout_mainland_prediction/README.md) · [Private-input-only reproducibility contract](docs/HELDOUT_MAINLAND_PREDICTION.md)

The next test asks **how well island floral-profile identities can be predicted** when the island's observed genus composition is known. Every name recorded on any of the **74 original geology-matched island endpoints** is excluded from **all 748 native-GIFT mainland training regions**, eliminating exact taxon-name leakage. The independent predictors use either global mainland colour/form frequencies, nearest-mainland floral repertoires, known genus profiles from other mainland species, or nearest-source genus profiles.

**Main positive result:** Knowing genus reduces out-of-taxon exact categorical log loss by **0.709 All** (37-pair bootstrap range 0.644–0.774) and **0.338 Direct** (0.309–0.369) versus the global mainland distribution, alongside gains in original category-distance accuracy. **Nearest-mainland genus information improves categorical-distance losses** by about 0.014 All/0.011 Direct at fixed shrinkage 8, but does **not robustly improve exact-profile log loss**; its sign changes with smoothing. Global genus donor coverage is ~73% All/~66% Direct; five-nearest mainland genus coverage only ~38%/~28%.

Within 28 explicitly native-supported matched pairs, the genus predictive benefit remains. Residual exact-profile geological origin differences under the imperfect held-out predictive model do not contradict the separate, more rigorous genus×regional occurrence randomization **non-rejection** in PR #5. These tests **do not forecast taxon/genera colonization**, reconstruct true ancestral sources, or prove selection/pollination function.

## Genus and potential mainland-source conditioning (2026-11-08)

[Verified conditional analysis](results/20261108_genus_source_control/README.md) · [Rights-restricted replay contract](docs/GENUS_SOURCE_CONDITIONED_PROTOCOL.md)

The original **37-pair exact colour×form repertoire shift is observed**, but **is not clearly exceptional once the actual island×taxon incidence and genus×four-region geographic availability are held fixed** (1,999 conditionally shuffled whole floral profiles; pooled All exact L2 p=0.1035, Direct p=0.0675; none of the four global tests survives BH/BY, min BH q≈0.185). Genus-only global permutation likewise has no supported pooled contrast. The secondary tropical-genus residual is suggestive but not supported under the more restrictive regional-signature null's eight-comparison BH adjustment. An **explicitly native-only** sensitivity retains 28 original matched pairs and also shows no supported global excess.

The archived **GIFT 3.2** geographically near-native mainland source flora (748 source regions; 1,666,186 records) allows a separate test. Over 90% of the broad *recorded* profile types in each island are represented in the potential five-nearest mainland-source repertoire. Subtracting that normalized profile spectrum leaves the exact-profile origin contrast (All/Direct BY q≈0.00083), but not a consistent categorical-similarity-kernel contrast after FDR. **Source-spectrum subtraction is not the genus conditional null or an ancestral colonization model.**

Interpretation is now centered on **which taxa and floral repertoires become assembled** rather than a claim of evolutionarily smaller or functionally equivalent flowers on oceanic islands. The private trait source archives and taxon/island identity-level tables are never redistributed: only aggregate results, code and hashes are checked into fd.

## Matched floral repertoire results (2026-11-08)

[Matched repertoire report](results/20261108_matched_flower_repertoire/README.md) · [Protocol and inference limits](docs/MATCHED_REPERTOIRE_PROTOCOL.md)

The new analysis extends the geological-origin comparison to **which exact flower-colour × flower-form state-set profiles** are represented. In 37 outcome-blind 1:1 continental–oceanic pairs (23 northern-midlatitude, 14 tropical), the exact-profile repertoire composition differs under paired label swaps (**All/Direct BY q≈0.00042** across the fixed four global comparisons), but the **categorical-trait-similarity-aware** global test is not supported (All p=0.0828, Direct p=0.1609). No individual profile survives BH across all 496 tested categories. Colour/state presence is broadly shared across islands; the distinction lies in *combinations and their recorded prevalence*.

Geographically closer matched pairs (≤5,000 km, n=25; ≤3,000 km, n=14 northern-midlatitude) retain exploratory exact-profile differences. These are **post-hoc sensitivities**, not additional independent experiments. The result does **not** identify the source-pool/genus sorting process or establish geological selection, pollinator replacement or true functional niche equivalence; source identity remains the main unresolved gate.

## Geological origin follow-up (2026-11-08): continent-derived versus oceanic islands

The [geological-origin report](results/20261108_geological_origin/README.md) and [literature/novelty audit](docs/NOVELTY_LITERATURE_AUDIT.md) distinguish existing island trait-space research from the new question: **do geological histories change the *within-island* recorded floral dispersion and *between-island* functional separation relative to well-defined references?**

- Independently sourced Roeble (2024) and Zizka (2022) geological classifications mapped conservatively to **434/8,264** GSHHG islands, of which **260 All (130 continental / 130 oceanic)** and **243 Direct** meet the complete-floral support gate. Unmatched and mixed islands remain unresolved; continental distance is never used as a geological label.
- Separate Weigelt (2013) last-glacial mainland-connectivity data matched **1,334 All** support islands. A last-glacial unconnected island is **not** automatically geological oceanic.
- Oceanic All islands have lower observed equal-profile Rao Q (median **0.75705 vs 0.77447**). An exploratory area, distance, climatic, region, observed richness and trait-coverage adjusted coefficient is **−0.00589** (uncorrected spatial-block wild-bootstrap 95% interval [−0.00988, −0.00227]). Similar direction for Direct.
- **Crucial falsification:** For the 12 testable of 32 specified *origin-difference* tests under independent binary within-origin profile-margin nulls, **0 were BH/BY-supported**. Class-specific profile counts/prevalence can already reproduce apparent local Q origin differences; these observations do **not** establish enhanced packing, selection, within-lineage evolution, source-independent geological causation or pollinator functional replacement.
- The local and pairwise outcomes are **categorical colour/form phenotype proxies**. 2-axis KDE is a projection (15.35% positive inertia), not the full floral-functional hypervolume.

The paper's next identification gate is **profile-repertoire composition × source/genus/arrival constraints**, not treating any null rejection as a positive ecological mechanism.

## Implemented G1 and projected-hypervolume extension (2026-11-08)

The [verified analysis report](results/20261108_g1_hypervolume/README.md) records the successful 11-test CI and full results. Under the fixed binary profile margins, 13/32 exploratory comparisons pass BH adjustment, **none passes BY**. In northern midlatitudes, local equal-profile Rao Q is slightly lower while between-island nearest-profile dissimilarity is higher than the fixed-margin reference.

The new hypervolume implementation computes **finite categorical occupancy** and **90%-mass 2-D Gaussian KDE projected area** for 4,722 scope–island records, with globally fixed bandwidth/weighting sensitivities. Only 15.35% of the reference positive inertia lies in the first 2 PCoA axes, so these KDE areas are **not** the full functional floral-space hypervolume.

Definitions, nonclaims, and replay are in [G1 and hypervolume methods](docs/G1_HYPERVOLUME.md).

## Implemented fd stage: original-distance analysis

Original-distance metrics do not use PCoA; the supplementary Gaussian KDE is explicitly projection-limited and does not enter the functional-distance significance tests.

- **Local Rao's Q:** mean colour/form dissimilarity over pairs of recorded taxon-name units. Includes a finite-name version and an equal-profile-weight sensitivity.
- **Between-island functional nearest-profile dissimilarity:** symmetric, carrier-weighted average of the nearest recorded profile in the other island; zero for two identical profile supports even if carrier frequencies differ.
- **Cross-island mean profile distance:** carrier-weighted expected dissimilarity between two sampled name units, one from each island. This is **not** a pure between-island beta component.
- **Exact-profile Jaccard:** the original nominal overlap comparator.

All distances are computed from **equal-weight mean set-Jaccard** colour/form mismatches, not from PCoA coordinates. Pair comparisons reuse PR #321's already selected matched-pair IDs. They concern **whole-island** aggregate profiles, including any shared taxa. The raw taxon-excluded/within-genus tests remain in island.

## Replay from allowed public aggregates

Python 3.11+:

```bash
python -m pip install -e ".[test]"
pytest -q
python scripts/fetch_island_aggregates.py --output work/island_public
python -m fd.pipeline --input work/island_public --output results/run_fd_v1
```

The importer pins source commit `d4a5dafab5f0737ea6c1c1c0263fc20decf6035c` and fetches **only** public aggregate island-by-profile carrier tables, island metrics, the 496-profile reference, and the public matched-pair cohort. Actual file checksums go to `source_manifest.json`; the runner independently records input hashes in `manifest.json`.

Outputs: four `*_islands.csv` files, four `*_pairs.csv` files and `manifest.json`. All/Direct and single-state subsets overlap and are sensitivity scopes, not independent replications. No species-level traits, private source ledgers or restricted redistribution material enter this repository.

The [GitHub Actions workflow](.github/workflows/analysis.yml) runs tests and attempts the pinned public-data analysis, attaching its outputs as a workflow artifact rather than automatically committing them.

## Claims gate

The public-aggregate stage provides **descriptive measured distances only**. It cannot estimate genus/source-conditioned null significance, exclude shared taxon names, identify pollinators, prove ecological functional equivalence, establish evolutionary divergence, or attribute differences to selection.

Next inference stage: implement function-distance nulls while retaining the exact existing PR #321 selection and distinguishing a fixed-profile-margin test (public aggregates suffice) from a genus/source-conditioned test (requires rights-cleared private input). Then audit area, isolation and trait-coverage effects on the *same* cohorts before manuscript claims.

**Ownership:** `island` is the frozen data/Chapter 1 provenance and PR #321 exploratory origin. `fd` owns the new distance-based floral-functional diversity question, implementation, validation and any eventual separate manuscript.
