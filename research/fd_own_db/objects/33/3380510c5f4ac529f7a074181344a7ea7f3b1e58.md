# 14 November 2026 — Does tropical island colour × architecture coupling survive taxonomic/source/native matching?

**ORIGINAL DATA analysis completed: both seeds SUCCESS.** Exploratory, conditional taxonomic-assembly test; no actual local bee/bird/moth pollination effectiveness, island-specific phenotype evolution, or evolutionary source identity observed.

- Original rights-protected [island workflow #43566068144](https://github.com/zuizui0223/island/actions/runs/43566068144) SUCCESS for both independently seeded 20261114 and 20261115 source executions; source artifacts [#14102630812](https://github.com/zuizui0223/island/actions/runs/43566068144/artifacts/14102630812) / [#14102741308](https://github.com/zuizui0223/island/actions/runs/43566068144/artifacts/14102741308).
- Reusable code: [fd native genus/family matched pairing](../../src/fd/native_taxon_matched_pairing.py), [synthetic guards](../../tests/test_native_taxon_matched_pairing.py); protected original [island runner](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_wcvp_taxon_matched_pairing.py) and [workflow](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/.github/workflows/run-fd-native-taxon-matched-pairing.yml).
- Original frozen Direct species morphology: **2,234** accepted names, three high/medium single-state colour/form/symmetry records. Original observed complete island floras **928**, GIFT strictly candidate native mainland donor source flora; ALL accepted plant names from the original focal island flora are omitted from its donor sources, irrespective of local nativeness. **2,190** Direct names have strict CoL XR original GBIF-family-agreed *candidate* taxonomy, which is still unreviewed.
- Original native-status eligibility of the >=30-Direct-species island population: **928 all recorded; 357 WCVP regionally native-compatible; 87 strictly GIFT-native**. WCVP TDWG-L3 regional distribution does not prove individual island nativeness; no originally introduced name is upgraded.
- **66 regional×origin×matching-rank×sample-n analytical panels**, **40** satisfy the outcome-blind minimum >=10 physical islands, >=5 geographic blocks and >=5 different first-nearest GIFT mainland source signatures. Unsupported cells do not receive a conditional effect.

## Question and design

The preceding finite-population analysis [WCVP native Rao/pairing](WCVP_NATIVE_RAO_PAIRING_928.md) reported original **source-relative tropical colour–architecture pairing departure** as positive, while source-relative marginal diversity was not consistently reduced:

- **106 tropical all-recorded islands**, exact finite *n=30* species E[Qjoint-observed minus colour-shuffled] source difference **+0.00651**, block CI **[+0.00258,+0.01077]**;
- **36 WCVP regional-native-compatible tropical islands**, same *n=30* outcome **+0.00881**, block CI **[+0.00269,+0.01605]**, source CI **[+0.00288,+0.01818]**.

This index is **Dpairing = (n−1)/n × (S_colour × S_(form,symmetry) − S_joint)**, the signed difference in expected categorical joint Rao between observed data and an analytical fixed-colour-margin whole-pool shuffle that leaves intact form+symmetry pairs. It is **not mutual information**, and even a positive *island-minus-source* value does not mean island absolute pairing excess is positive or that flower combinations are known pollination syndromes.

Now match taxonomic identity *before inspecting any triple-phenotype outcomes*:

1. **Same accepted genus, different accepted plant species**, one matched plant from each of n identical genus IDs on island and donor.
2. **Same strict CoL XR candidate family IDs**, one different accepted plant species per n matched family, potentially same genus.
3. **Same candidate family IDs but different genera inside every matched family**, one different accepted plant species per family.
4. Separate **count-only control** on the SAME source-supported island cohort, exactly n independent genus/family IDs per side but NOT required to have the same IDs; paired matched-minus-count-only uncertainty is reported. This is a different subset and estimand from n=30 full-flora source-relative Rao.
5. Fixed **n=5 and n=8** rank-identity units, **40 independently seeded rarefactions within each physical island** for each matching rank/sample n; one island average is the analysis row, and those 40 repeats are NOT independent island histories. Every source taxon name is excluded from the original complete focal island flora. Categorical marginal Rao, joint Rao and analytic pairing departure are computed with the n-specific (n−1)/n correction. Whole physical 10° spatial-block- and first-nearest-GIFT-source-signature-equal **1,999** conditional bootstrap draws are reported independently; those source/geographic groupings overlap.
6. CoL XR matching is strict exact accepted candidate name and original reported family/genus agreement, never a true tip tree, curated phylogenetic history or proven ancestral source.

No outcome-driven relaxation of the genus/family/source support minimum is permitted.

## Main result — tropical ALL recorded: broad +pairing is not confirmed after taxon matching

Each result below: **source-relative pairing departure on genus/family ID-matched n5 or n8 cohorts**, block-equal mean with conditional 95% interval. Seed #1 then #2; n-sensitivity cohorts are NOT identical to n=30 full-flora estimand.

| Taxonomic match | n species per side | Source-comparable islands / blocks / first GIFT donors | Seed 20261114 | Seed 20261115 |
| --- | ---: | --- | --- | --- |
| Genus identity | 5 | **81 / 33 / 30** | +0.00812 [−0.00150,+0.01781] | +0.00731 [−0.00143,+0.01665] |
| Genus identity | 8 | **50 / 21 / 19** | −0.00077 [−0.00988,+0.00957] | +0.00123 [−0.00702,+0.01028] |
| Family identity | 5 | **106 / 42 / 36** | −0.00147 [−0.00471,+0.00194] | −0.00016 [−0.00325,+0.00301] |
| Family identity | 8 | **102 / 41 / 35** | −0.00109 [−0.00471,+0.00287] | +0.00054 [−0.00328,+0.00477] |
| Family identity, different genus | 5 | **106 / 42 / 36** | −0.00031 [−0.00378,+0.00372] | −0.00166 [−0.00449,+0.00157] |
| Family identity, different genus | 8 | **99 / 40 / 33** | −0.00001 [−0.00344,+0.00379] | +0.00042 [−0.00327,+0.00440] |

**All six identity-match panels in each seed have conditional 95% block intervals spanning zero.** All corresponding first-GIFT-source-equal intervals span zero too. Genus n5 has a positive point estimate but broad CI, and genotype-family tests do NOT preserve a supported positive tropical departure. This is evidence against *promoting the original source-relative positive full-flora contrast as an already established lineage-independent tropical syndrome*, NOT proof the true within-genus or within-family effect is zero and NOT a causal mediation fraction.

## Stronger tropical WCVP regionally native-compatible subset (previous baseline 36 islands)

| Taxonomic match | n | Matchable / blocks / donor signatures | Seed 20261114 pairing | Seed 20261115 pairing |
| --- | ---: | --- | --- | --- |
| Genus | 5 | 30 / 16 / 13 | +0.00045 [−0.01362,+0.01640] | +0.00385 [−0.00986,+0.01898] |
| Genus | 8 | 14 / 9 / 6 | −0.00517 [−0.01440,+0.00436] | −0.00640 [−0.01633,+0.00332] |
| Family | 5 | 35 / 20 / 16 | −0.00221 [−0.00655,+0.00251] | −0.00307 [−0.00830,+0.00308] |
| Family | 8 | 35 / 20 / 16 | −0.00359 [−0.00893,+0.00177] | −0.00401 [−0.00990,+0.00238] |
| Family, different genus | 5 | 35 / 20 / 16 | −0.00263 [−0.00818,+0.00280] | −0.00104 [−0.00676,+0.00490] |
| Family, different genus | 8 | 34 / 19 / 15 | −0.00140 [−0.00740,+0.00467] | −0.00088 [−0.00664,+0.00527] |

**Every main matched-source pairing interval includes zero** under both seeds, whether equal-block or equal-first-source. Exact genus n8 is supported by only 14/36 originally region-native-compatible islands; it is not a global tropical within-genus test.

**Important count-only source-assembly signal (separate post-hoc estimand):** on the **SAME 35 physical islands**, equal **8 family units per side** without matching family IDs gives source-relative colour–architecture departure **+0.00349** (seed 1; block), while exact family identity matching yields **−0.00359**. The within-island change *family-ID match minus family-count-only* is **−0.00707 [−0.01270,−0.00129]** block; independently seeded **−0.00754 [−0.01318,−0.00176]** block. However first-GIFT-source-equal interval crosses zero in seed 1 (**−0.00636 [−0.01325,+0.00015]**) and excludes zero in seed 2 (**−0.00801 [−0.01382,−0.00195]**). Therefore family identity/composition is a plausible and testable contribution, but NOT established as a source-independent causal filter, and these 35 selected islands differ from the original 36 and 106 cohorts.

## Strict GIFT-native tropical subset: no decisive within-genus or within-family syndrome

Only **17** tropical floras have >=30 strictly GIFT-native Direct-triple species, yielding 10 geographic blocks/9 first source signatures. In exact-genus match n5 only **14 islands** have source support; n8 only **5 islands/3 blocks/3 sources**, below the minimum so NO n8 effect is reported. All available n5/n8 family-matched pairing source intervals span zero in both seeds, with 17 candidate native islands and provenance-selected names. More independent sampled native source populations would be needed to reject a modest within-lineage effect.

## Secondary northern midlatitude result: some family-with-different-genus residuals remain

Do not collapse all global regions into a positive or negative syndrome:

- **All-recorded** 736/735 source-comparable northern-midlatitude family-matched floras have *pairing* interval spanning zero after exact five/eight family ID matching, both seeds and both geographic/source groupings.
- **WCVP compatible** northern-midlatitude family-ID different-genus n8 284 islands / 32 spatial blocks / 39 candidate donors: **−0.00510 [−0.00998,−0.00078]**, seed 2 **−0.00460 [−0.00882,−0.00068]** block; source-signature seed 1 **−0.00451 [−0.00920,−0.00058]** but seed 2 **−0.00310 [−0.00713,+0.00039]** crosses zero.
- **Strict GIFT-native** northern-midlat, only **67 islands / 16 geographic blocks / 23 first donor signatures**, same family but DIFFERENT GENUS n8: **−0.01243 [−0.02028,−0.00590]** seed 1 block, **−0.01091 [−0.01657,−0.00612]** seed 2 block; source-equal **−0.00891 [−0.01521,−0.00348]** / **−0.00659 [−0.01180,−0.00155]**. Here one conditional descriptive negative residual survives both seed and weighting choices.
- Importantly, **exact same-genus n8 on 66 strict-native northern-midlat islands** is *not* directionally supported: block +0.00272 [−0.00302,+0.01004] / +0.00322 [−0.00305,+0.01099]. This contrast is consistent with genus-level turnover inside matched families, **not proof** of a within-genus adaptation or functional pollinator shift.
- None of the 66 nested/exploratory tests has a pre-registered independent multi-comparison family or phylogenetic origin-disjoint null. Do not promote the regional conditional n8 interval as a confirmatory global evolutionary result.

## Biological interpretation, claim ceiling, and priorities

### Supported observational description

1. The previously observed *source-relative tropical full-flora flower-pairing difference* does **NOT** have demonstrated genus- or family-identity-independent persistence in the current original source-matched n5/n8 controls.
2. The WCVP region-native tropical 35-island matched-family n8 vs family-count-only negative difference is compatible with **plant family composition/identity shaping recorded joint floral states**, but it has source and sampling sensitivity and cannot identify historical causal selection.
3. The northern midlatitude strictly GIFT-native family-but-not-genus matched residual shows **heterogeneity at botanical grouping depth**. Its interpretation is source-dependent and post-hoc; no named pollinator guild has been observed to cause it.
4. This strengthens an explicitly *lineage/provenance/source-conditioned floristic assembly* paper. It weakens a universal claim that island guild loss causes within-lineage floral syndrome convergence.

### Not established

Historical pollinator disappearance, a bee→bird/moth/fly replacement, actual island guild-service redundancy, pollen transfer/seed output, individual-island indigenous status from WCVP regional range, true dated phylogenetic independent tropical evolutionary contrasts, trait divergence within species, or selection as opposed to checklist/provider/source limitations. Non-rejection does not prove that taxonomic matching explains all floral differences. Two independently seeded finite simulated draws improve numerical stability but are NOT separate biological data sets.

### Next strongest falsification

Use **true dated plant tree branches**, matched island and candidate source **native species supported in the same geographic region**, and per-family or finer-clade outcome-blind matching with a demonstrably larger number of independent tropical donor species. Freeze a *single primary pairing effect* and spatial/source multiplicity plan before obtaining additional original pollinator data. Separately require standardized within-island/site/year visitor surveys and plant-specific pollen addition or seed set to infer ecological functions. Keep PR **Draft**; these conclusions are conditional floral phenotype assembly observations, not a proven adaptive pollination syndrome.
