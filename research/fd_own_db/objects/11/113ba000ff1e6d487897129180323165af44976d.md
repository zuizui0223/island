# 14 November 2026 — Global FD conditioned on taxonomic composition WITHOUT changing the observed species mask

Original protected-data execution **SUCCESS**: [island Actions #43600916175](https://github.com/zuizui0223/island/actions/runs/43600916175). Analysis code was pinned to fd commit `13dc85e154205bc60ab795947bb8fea43a4db612`; workflow commit `84fab6329f1865c5e1a81124f5fa2cdb4a840ed6`. Anonymous artifact **14115299795**, ZIP SHA256 `9282b2928268ba26bc80495d83bfce573049fd677dda46e7dee8a457b59a66d6`.

## Primary scope stays GLOBAL

All **4,505 original GIFT botanical floras**, ten independent directly documented traits and three genuinely observed two-trait combinations remain in the main frame. All-recorded and original-source-native flora scopes remain separate. There are **117,130 original island × scope × endpoint rows**, not 117,130 independent observations. Geology is available for **425** actual floras; the one flora without a matching geographic crosswalk remains in the global denominator. No complete-three-trait, common-genus, phylogenetic-tip or photo-availability requirement was added.

The independently archived prior raw-FD aggregate was reconciled on **all 26 endpoint/status combinations**. Every scored island count is identical, and the maximum absolute difference in original raw mean finite-five categorical Rao is **1.11e-16**. Raw, genus-conditioned and family-conditioned scores therefore do not use different selected islands or different phenotype records.

## Method

For observed category counts n_c among n recorded species, Q5 = 0.8 [1 − sum_c n_c(n_c−1)/(n(n−1))]. This is the exact expected empirical Rao of five distinct sampled species, NOT a multidimensional hypervolume.

For each trait, globally permute its observed labels only within original botanical genera, or separately within source-consistent candidate families. Each island's recorded species and taxonomic counts remain fixed. Exact hypergeometric moments give the expected Q without Monte Carlo label simulation. For a pair, permute only the first trait and hold the second trait observed; no third trait is required. The residual is observed minus conditional expected Q, not an adaptation coefficient.

The union of trait-recorded species contains **59,836 names**, of which **57,571 (96.21%)** pass the existing strict source-consistent CoL XR candidate-family gate, across 391 candidate families. The other **2,265 species are retained**, not dropped: their labels are fixed in species-specific singleton strata. Candidate family labels are not reviewed taxonomy, and genus/family conditioning is not a dated phylogenetic correction.

## Global recorded-flora results: same sample before and after conditioning

All-recorded flora, island-equal descriptive means. Trait-specific samples differ across rows, so do not rank ecological importance by comparing raw Q across traits.

| Endpoint | Islands with >=5 observed species: raw / genus-null / family-null | Observed mean Q5 | Genus-conditional expected Q5 | Mean residual |
| --- | ---: | ---: | ---: | ---: |
| Primary flower colour | 3132 / 3132 / 3132 | 0.575305 | 0.589826 | −0.014520 |
| Floral form | 2806 / 2806 / 2806 | 0.640721 | 0.649391 | −0.008670 |
| Floral symmetry | 3145 / 3145 / 3145 | 0.304382 | 0.305468 | −0.001086 |
| Flower size class | 2781 / 2781 / 2781 | 0.526215 | 0.535789 | −0.009574 |
| Inflorescence display | 3212 / 3212 / 3212 | 0.467411 | 0.456099 | +0.011312 |
| Autonomous selfing capacity | 2411 / 2411 / 2411 | 0.336025 | 0.347197 | −0.011172 |
| Colour × floral form joint Q | 2324 / 2324 / 2324 | 0.720544 | 0.720837 | −0.000293 |

These are conditional descriptive averages, not significance tests or evidence of a universal island simplification. The original source-native subsets remain separately reported (e.g. 397 colour, 327 colour × form islands). Trait invariance is explicitly audited: for symmetry, **380/3145** all-recorded scored floras have no shufflable genus stratum. A zero residual in such a flora is not evidence of evolutionary stasis.

Observed Q is also partitioned exactly into same-stratum and different-stratum species-pair contributions. The proportion contributed by different-genus pairs is NOT explained variance or a causal percentage of assemblage formation; it partly reflects the number of cross-genus species pairs.

## Tropical oceanic-minus-continental comparisons on IDENTICAL masks

Geographic-block-equal contrasts; brackets are exploratory 999-draw conditional block-bootstrap 95% intervals. Shared blocks receive the SAME resampling multiplicity in both origin classes. All columns of a row use the same islands and observed species.

| Endpoint / original flora scope | Oceanic / continental islands | Raw Q difference | Family-conditioned residual difference | Genus-conditioned residual difference |
| --- | ---: | ---: | ---: | ---: |
| Flower size / original native | 77 / 29 | +0.04018 [+0.01986,+0.06015] | +0.01915 [−0.00065,+0.04059] | −0.00012 [−0.01511,+0.01367] |
| Colour × form joint Q / original native | 63 / 26 | −0.03271 [−0.06470,−0.00830] | −0.02126 [−0.04356,−0.00553] | −0.00527 [−0.01198,+0.00160] |
| Autonomous-selfing category Q / all recorded | 119 / 37 | −0.04354 [−0.06978,−0.01397] | −0.03645 [−0.05883,−0.01135] | −0.01524 [−0.03485,+0.00555] |

The positive raw flower-size contrast in documented native floras becomes approximately zero relative to the genus-conditioned expectation. The negative colour/form joint-Q contrast becomes much smaller. Unlike earlier matched-subset comparisons, this is not caused by discarding different islands/species between columns. It makes taxonomic composition an explanatory candidate; it does not establish colonization, extinction or adaptive change causally.

Results are not uniform: all-recorded flower-colour genus residual contrasts are **−0.00958 [−0.01834,−0.00075]** on 139/43 islands, whereas original-native-only colour gives **−0.00219 [−0.01366,+0.01076]** on 95/29 islands. Both scope and taxonomic-reference choices matter. The output contains 546 regional/trait/scope/metric contrasts, 315 with descriptive intervals; they are correlated and NOT multiplicity-corrected discoveries. Climate, shared species among islands and full phylogenetic dependence remain uncontrolled.

## Missing-state uncertainty without pretending to observe missing traits

All recorded flora members are retained in denominators. Conditional on the observed global category codebook being exhaustive, putting every missing state in a modal category minimizes Q; integer water filling maximizes Q without changing an observed state. These are mathematical sensitivity bounds, not imputed traits or confidence intervals.

Across the **3,580** all-recorded floras with at least five listed plants (including those with insufficient actual trait records), median PER-ISLAND bound widths are **0.35792 for colour, 0.50208 for form, 0.56465 for autonomous selfing, and 0.68629 for colour × form joint Q**, against a Q5 upper limit of 0.8. Thus full-flora generalization still depends strongly on the missing-state mechanism.

Single-island extrema are sharp under the closed codebook. Pair bounds relax information from species observed for just one trait, and globally shared unknown species are not forced to take one consistent category across every island. Regional combinations are conservative outer ranges, not jointly sharp bounds. Unrecorded new categories could widen the assumptions further.

## Validation and interpretation

Eleven mathematical/retention/baseline tests passed locally; the protected analysis workflow also passed the original partial-trait regression tests. An additional **120 randomized exhaustive-permutation cases** agreed with the analytical expectation to **2.22e-16**. Rechecking **1,380 additive identities** in the actual anonymous output gave maximum error **3.33e-16**. This is validation of this module and source execution, not a claim that all repository-wide CI is green.

The scientific main line remains global floral assemblage composition, with oceanic/continental context and taxonomic conditional benchmarks. A residual can guide hypotheses about post-establishment differentiation but is NOT by itself an effective pollinator shift, ancestral reconstruction, local population phenotype, selection or adaptation. No private source taxon/island/family identifiers were exported.

Code: [global_taxonomic_mask_fd.py](../../src/fd/global_taxonomic_mask_fd.py). Tests: [exhaustive expectations and missing-data extrema](../../tests/test_global_taxonomic_mask_fd.py).

Methodological context: Penone et al. 2014, doi:10.1111/2041-210X.12232, on missing trait records; de Bello et al. 2011, doi:10.1111/j.2041-210X.2010.00071.x, on within/between trait variability. The finite categorical conditional expectations used here are explicitly derived and independently checked by enumeration. Keep manuscript inference exploratory.
