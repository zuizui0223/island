# Source × pollinator-guild availability: prospective identification and falsification contract
Date: 2026-11-16. Status: DESIGN ONLY; no new empirical coefficients or original-source run. Linked to fd PR #16 (4,505 original GIFT island floras) and separate fd PR #15 (guild/FD support audit). This is a prospective exploratory protocol, not a preregistration and not confirmation of a guild mechanism.

## Central question and explicit rival
Do mainland candidate-source floral phenotype proportions interact with independently measured island pollinator-guild availability to explain island floral-sector representation, beyond finite source sampling, botanical lineage composition, geography and observation? Rival H0: source taxonomic composition, botanical ascertainment, geography and dispersal opportunity account for apparent losses without any incremental pollinator availability contribution.

## Population, evidence and non-imputation
- Retain ALL 4,505 botanical island floras in the frame, with trait-specific missingness and eligibility; do not restrict main inference to historic 928 triple-complete / 30-carrier subset.
- High/medium Direct species-specific morphology only. Single-axis traits analyzed separately: symmetry, form, true ordinal tube depth, size, and colour. Joint sectors use only genuinely observed SAME-accepted-species two-trait combinations; never multiply marginal traits to invent plant profiles. Three- or four-way joint analysis is sensitivity only when directly recorded.
- Separate reproductive axes (SC/SI, autonomous capacity, mating system); do not include them as flower-shape coordinates.
- Original species-wide labels are constant across islands: conclusions concern recorded floristic assemblage composition, NOT within-species evolution.
- Main candidate source definitions: geographically nearest five GIFT native mainland regions (geo5); independently defined climate-nearest ten among geographic fifty (climate10) as source-choice sensitivity. Neither is actual seed arrival or known ancestral mainland population. All 748 source entities are potential shared pools, NOT 748 independent replicates.
- Retain status labels native nonendemic / endemic / introduced / unresolved; never turn unresolved into native. Main all-recorded description and strict native-source estimands are explicitly different, nonrandom populations.

## Primary sectors (fixed ontology, not pollinator labels)
1. zygomorphic symmetry (vs actinomorphic and asymmetric separately);
2. tubular or bell-campanulate form; open-radial comparator;
3. truly deepest ordinal tube-depth category;
4. genuine SAME-species tubular-or-bell × zygomorphic joint sector; genuine open-radial × actinomorphic comparator;
5. primary colour and colour×form only where same-species original direct labels exist.
These are morphological proxies. Do not call a sector Bombus-, bird-, moth-, butterfly-, fly- or other-bee-pollinated without independent plant-specific visitor evidence. Define independent sector×guild affinity weights using external interactions with study/provider/reference provenance, and report their uncertainty; if unsupported, leave w[k,g] missing and do NOT fill it from morphology or focal outcome.

## Estimands and arithmetic
For island i, sector k, and source definition m, let s[i,k,m] be source Direct-documented sector frequency and y[i,k] island Direct-documented sector frequency with denominator recorded separately. Let y0[i,k,m] be the finite source-draw expected island sector share, conditional on original source listing breadth, carrier count, recorded-state documentation, and optionally source-distance strata, without using observed island sector outcomes to tune source choices. Define source-relative deviation r[i,k,m] = y0[i,k,m] - y[i,k] (positive means deficit relative to reference); report raw y and s alongside r. This is recorded-share departure, not historical loss.
For guild g, define mainland availability M[i,g] and island availability I[i,g] as independently effort-calibrated probabilities or documented three-state retained/adequately nondetected/unresolved indicators. A putative lost channel can only be defined for documented source-positive guilds with qualified island nondetection. For probabilistic data, L[i,k] = sum_g w[k,g]*M[i,g]*(1-I[i,g]), with weights pre-fixed independently of island flower outcomes. Unresolved island guilds stay missing/interval-censored, never I=0.
Primary model comparison on eligible island×sector units:
  H0: r ~ source-sector share + log(area) + isolation + climate + Direct coverage + carrier effort + spatial/source effects;
  H1: H0 + guild loss opportunity + source-sector share × guild loss opportunity.
Main test is incremental interaction, with pre-fixed direction beta_interaction > 0 for source-rich sectors, estimated only within common support and reported with uncertainty. If r is mathematically dependent on s, use conditional finite-source draws / posterior predictive calibration and a pre-fixed residual model; compare observed statistic to its full source-conditioned null, not naive significance from raw (s-y) ~ s.
Secondary within-genus fixed-label null removes genus turnover and asks whether any residual association persists within the already assembled botanical genera; it is NOT the primary total assembly estimand because genus turnover may mediate filtering. Report both, never call the difference causal mediation.
Functional responses: (i) per-sector actual category occupancy/absence relative to finite frequency and genus-preserving nulls, (ii) island-source compositional displacement in a fixed Gower/PCoA space, (iii) original-distance Rao Q for trait-specific available states, and (iv) exact SAME-species joint Rao Q/Hill2. All with source carrier rarefaction and source-block uncertainty; 2D PCoA display is NOT total hypervolume; historical 17D KDE is supplementary only.

## Mandatory identification gate BEFORE fitting H1
For each guild × region × source-signature stratum: count retained, adequately nondetected, unresolved island IDs, distinct 10-degree geographic blocks, independent nearest-source signatures and observation providers; quantify source-rich/source-poor overlap. Minimum exploratory comparison: >=10 retained AND >=10 adequately nondetected islands, each spanning >=10 spatial blocks, and >=3 source signatures with BOTH states, as in PR #15. An additional stricter same-physical-block comparison check requires >=5 original 10-degree blocks containing BOTH states. This particular within-block gate was introduced on 2026-11-16 AFTER the first source/guild support counts (3 tropical bird slots) were inspected, but BEFORE fitting any guild effect: it is a transparent post-support identification safeguard, not a preregistered criterion. Require finite and nondegenerate source-sector share variation in both states and >=5 independent Direct plant carriers per trait endpoint. Pre-specify effective source-block independence and leave-one-archipelago-out robustness. If any gate fails, publish a coverage table and NOT_IDENTIFIABLE; no coefficient, no sign, no p-value, no opportunistic pooling across regimes.
The 2026-11-13 original PR #15 audit failed ALL 20 guild×four-region contrasts. Examples: Bombus northern-midlatitude retained/adequately nondetected 480/5, tropical 2/16. This is the current identification bottleneck, not a null ecological effect.

## Geography, geology, status, and alternative mechanisms
- Main geography uses 3 mutually exclusive regions; >=60-degree high-latitude nested sensitivity is NOT a fourth independent replicate. Southern extratropical is heterogeneous; split South America / Australasia / southern Africa only if independent blocks and support allow, never invent a uniform southern prediction.
- Geological type is independently verified for only 425/4,505 floras (249 oceanic, 176 continental), selected subset. Compare origin after matching common source opportunity, geography, climate, effort and guild support; report overlap failures. Geographic isolation is not a guild-absence proxy.
- Source-genus composition, genus sharing and island-only/source-only taxa are separate taxonomic-assembly endpoints. Do not infer true arrival/establishment failure from unrecorded mainland candidates.
- Endemic vs nonendemic congeneric contrasts are descriptive distinct-species differences; evolution requires sister/ancestor phylogeny, local phenotype and independent effective visitation/reproduction. Introduced taxa are not a natural randomized dispersal experiment.
- Negative controls: non-affinity sectors, unrelated guilds, source reassignment/placebo with matching source richness and climatic geography, plant Direct documentation and GIFT status missingness, spatial-block/source-group holdout, within-genus source reference, alternate geo5/climate10 source definitions.
- Fixed multiple-testing family for primary guild×sector interactions before outcomes; control Holm or hierarchical shrinkage, and publish all tests including unsupported slots. No significance claims from 16 posthoc sectors or 199 reselects.

## Decision table
A. Source-sector and genus turnover patterns supported, guild overlap fails -> source-conditioned floral assemblage sorting only; guild mechanism NOT_IDENTIFIABLE.
B. Guild overlap passes but interaction CI includes zero / null-calibrated effect not supported -> no incremental guild association demonstrated.
C. Interaction survives source-only, spatial/source, observation, status and independent guild-affinity controls -> source-dependent guild-associated assembly, still NOT historical decline or adaptive pollinator shift.
D. Same-plant actual visitor effectiveness, temporal data and independent sister/ancestral phenotypes additionally support it -> only then test pollinator-driven evolutionary change.
Any source, guild, or flower status missingness must be explicitly tabulated rather than imputed as absence.

## Existing empirical anchors (not newly recomputed here)
- PR #16: original 4,505 flora source/genus finite-sampling tests, source-conditional tropical zygomorphic −5.37 pp and genus expected −5.49 pp; high documentation imbalance.
- PR #15: north-midlat 743 all-recorded islands source-relative marginal Q −0.02571, joint Q −0.01420; guild positivity 0/20, strict-native selected coverage.
- 2026-11-16 original geo5 vs climate10 source audit: same qualitative source-genus deficit under alternative pools, not pollinator identification.
Links:
https://github.com/zuizui0223/fd/pull/16
https://github.com/zuizui0223/fd/pull/15
https://github.com/zuizui0223/fd/blob/analysis/discrete-functional-hypervolume-20261113/results/20261115_global_island_functional_beta/DIRECT_SOURCE_GEO5_VS_CLIMATE10_ORIGINAL_4505.md

## Implementation order
1. Produce original protected-source island×guild×region×source signature eligibility/positivity aggregate table with no identifying taxa or coordinates.
2. Produce sector-specific source-null standardized observed-minus-expected and trait-documentation support for all 4,505 floras.
3. Fit source-only sector models and null calibration; do not fit interaction unless positivity passes.
4. If supported, fit H1 and validate via geographic/source-held-out prediction and outcome-blind placebo.
5. Publish negative results, attrition, and scientific claim ceiling. Keep PR DRAFT until source-backed run, synthetic regression tests and independent ecological evidence pass.
