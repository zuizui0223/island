# 2026-11-16 — Analytic within-genus fixed-margin null for *real* symmetry×tube-depth floral trait coupling

**New original-source-protected 4,505-island analysis SUCCESS**: [island Actions #43861441242](https://github.com/zuizui0223/island/actions/runs/43861441242), 14 synthetic tests passed, hash-verified original Direct species-trait ledger and botanical flora, independently rebuilt island physical 10-degree blocks. Anonymous original [artifact #14224515054](https://github.com/zuizui0223/island/actions/runs/43861441242/artifacts/14224515054), extracted JSON SHA256 `1ead3ca06bb3a915f9db7f901b00a45ea756d077aa46331ae7524d9860ba997a` (244,413 bytes). Source implementation: [pure analytic genus null](../../src/fd/real_joint_within_genus_coupling.py), [tests](../../tests/test_real_joint_within_genus_coupling.py), [original protected runner](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_real_joint_within_genus_coupling.py).

All original island accepted species, genus labels, full source-region/GBIF identifiers and coordinates remain ONLY in protected original CI. No original taxon/island lists or locations were published. Complete original **4,505** island botanical flora framework (not n928 complete cases) retained; **45** analytical region×status×true-two-trait morphology panels, **24** with >=30 eligible islands and >=10 physical blocks.

## Question: Does source-relative floral trait *coupling* remain after preserving genus-level A and B trait margins?

The prior actual same-species `fd` results showed on tropical islands a source-relative excess of `actinomorphic × deep` combinations and a source-relative scarcity of `zygomorphic × deep` combinations beyond the product of each trait's ordinary marginal frequency (real joint trait co-occurrence `C=P(A∩B)-P(A)P(B)`). The independent exact FOUR-term genus accounting decomposed the *JOINT PREVALENCE* `J=P(A∩B)`, not the non-linear `C`. Calling the latter independent of genus without an additional model would be incorrect.

The NEW **analytically exact** conditional reference takes separately each focal island Direct species roster and each of its TWO candidate GIFT native mainland source rosters (geographical nearest-five; climate-nearest-ten among geographical fifty). It retains only accepted plant species with **BOTH relevant traits directly described for the SAME accepted species**. For each already documented botanical genus `g`, retain number of carriers `n_g`, count of trait A `a_g`, and count of trait B `b_g`. Under complete exchange of B between already Direct-coded plants **within that genus**:

- `E[J | observed genus margins] = (1/N) * sum_g (a_g*b_g/n_g)`.
- `E[C | genus] = E[J|genus] - (A/N)*(B/N)`.
- `R = C_observed - E[C|genus] = (J_observed - E[J|genus])`.
- Exact conditional finite-hypergeometric variance of joint count within genus is `Var(k_g)=a_g*(b_g/n_g)*(1-b_g/n_g)*(n_g-a_g)/(n_g-1)`, zero for singleton or fixed-margin genera; sum across genera, scale by `N^2` for coupling variance.

**The source-relative observed co-occurrence contrast decomposes exactly as `C_island-C_mainland = (E_genus C_island - E_genus C_mainland) + (R_island-R_mainland)`.** This is an in-sample genus/margin-preserving mathematical reference, not a phylogenetic causal model, observational time trend, independently measured ancestral state, visitor mechanism or experimental test.

All source mode choices are independently fixed original GIFT regional native botanical checklists; each candidate mainland roster excludes ALL focal island accepted species names. Original scopes: all-recorded; strict status-known native; strict status-known native nonendemic. Each endpoint requires >=5 jointly Direct-described island accepted plant species and >=10 candidate source-only accepted plants in EACH mode. Additional screens require <=10 percentage point difference in direct-trait annotation coverage under BOTH modes and independently at least one genus with BOTH A and B variable within the island AND each candidate source (at least three such Direct carrier plants per side). Distinct source modes and groups are NOT independent biological replication.

## MAIN tropical original result: genus reference reproduces the large apparent 'access architecture' coupling

Original tropical **all-recorded** fixed-eligible 666 island floras, **89 physical blocks**:

| True same-plant structure | Observed `C_island - C_geo5`, percentage points | Expected from fixed **within-genus A/B margins**, pp | Difference beyond genus null, pp [95% block-bootstrap] | Source-group residual, pp [95%] |
|---|---:|---:|---|---|
| **Actinomorphic × deep** | **+4.018** | **+4.057** | **−0.039 [−0.061,−0.018]** | −0.012 [−0.030,+0.008] |
| **Zygomorphic × deep** | **−3.425** | **−3.639** | **+0.214 [+0.166,+0.266]** | +0.178 [+0.125,+0.237] |

Alternate climate10 source mode: actinomorphic×deep residual **−0.048** pp, zygomorphic×deep **+0.173** pp. These are very small relative to the observed source-relative coupling contrasts, but some conditional group intervals exclude zero. They are posthoc correlated and source documentation dependent, not Holm-corrected evidence for independent within-lineage adaptation. In the original <=10pp coverage-balanced **435** tropical all-recorded islands, radial/deep observed C remains **+3.460** pp, genus expected **+3.502**, residual **−0.042 [−0.065,−0.020]**. This is a *descriptive recovery* of the source-relative joint phenotype difference, not true causal explanation by historical lineage filtering.

**CRITICAL opportunity limit, not a small p-value discovery:** Among those 666 tropical same-plant symmetry×deep qualified island floras:
- `actinomorphic × deep`: only **4/666** have even ONE original island genus with independent *within-genus variation in both A and B*, compared with **228/666** candidate geo5 and **336/666** climate10 source rosters. **ZERO** original islands have a sufficiently variable genus opportunity in ALL THREE rosters under the fixed support test.
- `zygomorphic × deep`: **7/666** with any island-side informative genus, versus **270/666** geo5 and **372/666** climate10; again **ZERO** fully source/island informative contrasts.

In the **96 strictly native-known** tropical same-plant pair floras, **0** island floras have informative genera for actinomorphic×deep or zygomorphic×deep. Nevertheless source-relative observed C for actinomorphic×deep **+7.082** vs expected **+7.141** pp (residual **−0.059 [−0.096,−0.027]**) and zygomorphic×deep **−6.596** vs expected **−6.794** (residual **+0.199 [+0.115,+0.285]**). These are *algebraically mostly constrained* by singletons, fixed genus margins and/or divergent mainland reference information. **A near-zero residual where within-genus states cannot reshuffle must NOT be described as independent falsification of genuine island floral evolution.**

## Falsification on a more genus-variable morphology pair

For tropical **tubular/bell × deep** all-recorded n**598** islands /88 physical blocks, observed coupling difference **−1.097** pp, genus expectation **−1.634**, residual **+0.537 [+0.427,+0.642]** geo5 (climate10 residual +0.398). However **only 35/598** island floras have enough independently variable within-genus A/B states on island AND both sources (and only 26 also have <=10pp direct-documentation gap). On the more informative **35 islands /21 physical blocks**, geo5 residual **−0.112 [−0.652,+0.402]** and source-group residual **−0.026 [−0.507,+0.459]**, both crossing zero and with a different sign. Thus the striking unconditional residual **does not survive the actual informative-genus opportunity screen**. For strictly native-known, only **3/93** such island floras have informative genera; no statistically supported mechanism can be asserted.

For northern midlatitude, all-recorded symmetry×deep shows only **1/1443** jointly informative island; tubular/bell×deep has **111/1424**. In high latitude or southern extratropical, true symmetry×deep generally has zero independently variable island/source combined-genus opportunities. This is a global data-structure identification limit, NOT evidence that plants do not vary biologically.

## Scientific decision

**Source-backed positive result:** True same accepted-species floral accessibility architecture, particularly tropical `actinomorphic × deep` versus `zygomorphic × deep`, differs from potential mainland source floras beyond simple overall trait-margin products. However this contrasts very closely with the **in-sample fixed-genus-and-per-genus-trait-margins expectation**. This strengthens the description of **taxonomically structured trait-pairing composition** rather than a universal floral-specialization loss.

**What remains unidentified:** independent evolutionary modification of floral pairing within lineages, pollinator-dependent selection, real bird/bee/butterfly/moth shifts, differential settlement after propagule arrival, or a unique diagnostic bird/butterfly flower syndrome. Most original island botanical genera have insufficient jointly documented within-genus species and/or variation for the proposed additional test. Current original global GBIF effort-qualified guild states fail most region×guild support comparisons; GloBI and EuPPollNet local exact plant/visitor coverage is insufficient for global effective pollination inference.

**Next data acquisition priority:** obtain *independent, locality-resolved directly measured symmetry and tube depth for MULTIPLE co-occurring congeneric plant species and their potential mainland relatives*, alongside standardized flower-visiting guild observations with provider/effort, to create genuine within-genus/visitor common support. A new genus-heldout phenotype donor or phylogenetic-branch heldout validation would also test whether the in-sample genus expectation predicts otherwise unseen island morphology. Do not convert observed morphological `actinomorphic×deep` into a bird- or butterfly-pollination state without independent visitors.

This finding is exploratory, conditional on original literature-coded species-wide Direct traits; it is NOT a new confirmatory guild-evolution result and [fd PR #16](https://github.com/zuizui0223/fd/pull/16) stays **DRAFT**.
