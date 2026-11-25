# Original GIFT 4,505-island genus-disjoint functional beta: verified first pass (2026-11-18)

**NEW ORIGINAL-DATA RESULT** — Source-protected [island Actions #44290748952](https://github.com/zuizui0223/island/actions/runs/44290748952) **SUCCESS** (including the two newly added synthetic genus-disjoint/unknown-name tests), after independent earlier [source run #44290656015](https://github.com/zuizui0223/island/actions/runs/44290656015) **SUCCESS**. The original source file SHA and accepted-species evidence reconstruction remained checked against the archival source. Only source-publication-anonymous aggregate means/counts were printed in the CI job logs; no original island, plant, genus, donor source identity or coordinates exported. The latest fd code is in [island_functional_beta.py](../../src/fd/island_functional_beta.py), the regression tests in [test_island_functional_beta.py](../../tests/test_island_functional_beta.py). This analysis adds a genuinely new endpoint but **does not** independently validate pollinator-mediated convergence.

## Second-stage verified falsification: matched source and disjoint blocks

The original SHA-locked matched source [#44293839934](https://github.com/zuizui0223/island/actions/runs/44293839934) and stronger independent-10°-block reselection [#44295181899](https://github.com/zuizui0223/island/actions/runs/44295181899) both **SUCCESS**. On 1,112 within-region donor/effort/distance-matched focal islands, the recorded full-flora genus-disjoint minus species-disjoint genus-sharing floral overlap is **−0.01826 focal-equal**, **−0.00646 geographic-focal-block-equal**, with exploratory focal-block 95% [−0.02143,+0.00893] that ignores reused donor blocks. The strict physically block-disjoint reselection admits only **19–26** triple-matched global units (tropical 8–13), with both positive and negative directions. **No robust extra similarity beyond shared-genus source opportunity was independently identified.** This is NOT a test of active guild replacement or local evolution. [Full matched source and independence report](MATCHED_KM_EFFORT_BLOCK_STRESS.md). The original study population stays 4,505 botanical island floras, and PR remains DRAFT.

## Design and original denominator

Parent source: all **4,505** accepted GIFT botanical island floras. On each island, use only true accepted SAME-SPECIES high/medium species-Direct colour × floral-form observations; ≥5 supported plant species per island, leaving **2,324** eligible island floras and all **2,699,326** unordered pairs. Functional overlap is **Σ_c min[p_i(c),p_j(c)]**, exact frequency total-variation overlap in [0,1] across actually observed joint flower categories; no imaginary cross-trait profiles, imputed missing states, or Gaussian 17D hypervolumes.

Define botanical genus identity from the first validated binomial of an original accepted plant name (not an invented phylogenetic graft). **Strict whole-flora genus-disjoint** means both islands have certified accepted-binomials for *all members in the selected original status scope*, and the **entire documented island floral species rosters** have zero shared genera; it does NOT merely check genera among Direct trait carriers. An unresolved accepted binomial disqualifies a prospective strict genus-disjoint pair, rather than implicitly certifying difference. Nested partial-Direct-genus-only estimates remain a different sensitivity.

**Crucial denominator qualification:** all 842,058 zero-shared-species original pairs are not partitioned by the 155,569 strict zero-shared-genus and 564,023 genus-shared counts. The remainder contains uncertified complete-genus rosters and must NOT be treated as genus-disjoint or genus-shared. Pairs and sequential numerical seeds are not independent biological lineages.

## Original verified main source totals

| Original all-recorded floral status; same real Direct colour×form island sample | Unordered pairs | Mean exact TV functional overlap | Original whole-genus-label-shuffle expected squared-frequency distance | Actual squared-frequency distance | Actual minus genus-null distance |
| --- | ---: | ---: | ---: | ---: | ---: |
| All eligible pairs | 2,699,326 | **0.386901** | 0.136539 | 0.147891 | +0.011352 |
| No shared accepted *species* (full flora) | 842,058 | **0.247183** | 0.207836 | 0.229663 | +0.021827 |
| No shared species, **shared genus** in certified full flora | 564,023 | **0.255535** | 0.198213 | 0.218174 | +0.019961 |
| **No shared genus** (entire certified flora, strict) | **155,569** | **0.214879** | **0.256538** | **0.280641** | **+0.024103** |
| No shared genus and different physical 10° blocks | **155,566** | **0.214880** | 0.256536 | 0.280639 | +0.024102 |
| No shared genus only among *Direct observed* plants (looser sensitivity) | 937,292 | 0.256808 | 0.206295 | 0.223890 | +0.017595 |

The **positive discovery is absolute category recurrence across separately documented botanical genera**: even with zero recorded genus overlap, the same colour×form categories have 0.2149 average compositional overlap. HOWEVER **actual squared-frequency distance is LARGER than the original genus-preserving category-shuffle expectation** (+0.0241), not smaller: this is **no evidence of excess selected morphological convergence**. This squared-distance null uses the global documented original species/trait pool and each fixed source island roster and measured species mask; it does not condition on geographic source/climate or ancestry.

## Regional original support: all-recorded status

| Original region; pairs within that region only | Strict genus-disjoint pairs | Unique eligible islands / original 10° geographic blocks | Mean observed colour×form overlap | Observed-minus-genus-shuffle squared distance |
| --- | ---: | --- | ---: | ---: |
| Northern midlatitude | **10,120** | 942 /75 | **0.227232** | **+0.008743** |
| Tropical | **321** | 139 /54 | **0.262854** | **−0.005682** |
| Southern extratropical | **56** | 42 /19 | **0.166336** | **+0.034324** |

All regional numbers are **descriptive**, do not compare geographic regions as if samples, source/evidence coverage and pair distances were equal. Tropical −0.00568 is simply a signed comparison of squared distance to a *particular* global within-genus reassignment expectation on the 321 qualifying pairs, **not** a significant or geographically source-independent excess convergence result. Other original geography regimes have no qualifying summary in this CI log and cannot be assumed to have no biological genus-disjoint species.

## Strict GIFT native-known original status sensitivity

Across **327** native-known colour×form Direct n≥5 islands there are **53,301** original native-roster island pairs. Source-certified strict **native-flora** genus-disjoint pairs: **2,106** spanning 235 original islands /67 independent physical 10° blocks; mean TV overlap **0.234340**, squared-frequency distance 0.313888 vs within-genus reference 0.302856 (excess **+0.011032**).

- North midlat native strict genus-disjoint: **265** pairs /71 islands /19 blocks, overlap **0.188778**, squared-distance excess **−0.011101** (exploratory descriptive).
- Tropical native strict genus-disjoint: **NO qualifying pair shown** in the original verified regional analysis, despite 91 original tropical native Direct eligible islands. Do not interpret as biological lack of shared-functional categories or proof of trait loss; source status support and taxonomic exclusion are too restrictive.
- The full **all-recorded** and **native-known** populations are *different subsets*, not effect estimates of introduction.

## Scientific conclusion (not overstated)

1. **Confirmed positive recorded observation:** sharing of broad flower-colour × floral-form category combinations is not mechanically dependent on sharing accepted plant **species or genera**. Original test across 155,569 full-roster genus-disjoint island pairs is a quantifiable, nonzero functional overlap.
2. **Important falsification:** observable nonzero overlap is **not** automatically an above-genus-null evolutionary convergence; squared functional dissimilarity is globally higher than genus-preserving shuffle expectation.
3. **Ecological mechanism unresolved:** genus-disjoint status is not phylogenetic independence; geographical region, different recorded genus richness, number of Direct species per island, floristic source, documentation effort, trait ontology and environmental availability are strongly different between same-genus and genus-disjoint island pairs. No effective pollination, pollen transport, selection or island-local phenotypic evolution measured.
4. **Immediate next control:** prospectively match G0 genus-disjoint and species-disjoint genus-shared comparisons **within same region and independent physical/source blocks** on inter-island km, minimum Direct plant support, full-flora genus richness and observation fraction. Recompute original equal-focal/cluster source uncertainty; preserve the existing sub-2,324-island Direct mask and no-missing-imputation. The frozen matching/null contract is now on [main](https://github.com/zuizui0223/fd/blob/main/docs/GENUS_DISJOINT_FUNCTIONAL_BETA_NEXT.md).

**Research positioning:** the robust positive finding is *cross-genus functional category redundancy amid botanical turnover*, embedded in heterogeneous regional floral repertoires. A mechanism of pollinator decline or replacement remains an explicitly separate identification gate. This first exploratory result belongs in draft PR #18 until outcome-blind covariate and geographic matching is completed, not promoted as an independently confirmed causal island syndrome effect.