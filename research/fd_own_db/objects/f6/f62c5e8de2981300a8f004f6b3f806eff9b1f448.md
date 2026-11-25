# 2026-11-17 — 4,505 GIFT island population: ORIGINAL common-island source-genus predictive asymmetry

**Original protected-source numerical workflow SUCCESS**: [island run #44024477079](https://github.com/zuizui0223/island/actions/runs/44024477079) and independent [replay with synthetic full-population tests #44024507798](https://github.com/zuizui0223/island/actions/runs/44024507798). Original anonymous artifact #14302457701 (first) and #14302267791 (tested replay). Replays match to floating-point rounding (~1e−17), not byte identity. Original GIFT accepted botanical **4,505** islands retained in source runner (including uncharacterized islands); **60** geography × nativeness × donor-gate × source-mode reporting panels, not 60 independent discoveries. Protected original accepted island and plant identities do NOT leave the CI artifact. This README summarizes actual original-source results; subsequent per-endpoint diagnostics are separate.

## Estimand and limits

For each island, flower structure endpoint t, source mode m, and documented source-only genus donor cohort, compute **gain(i,t,m) = Brier(source pooled model) − Brier(source genus-weighted model)**, where lower Brier means better probabilistic prediction of a real same-accepted-species flower joint state. Each source donor is a different accepted species; any name occurring in the entire focal island flora is excluded from both mainland source candidate modes. Genus weighting is **conditioned on actual observed island botanical genus composition**, not predicted botanical arrivals. Brier scores average across eligible documented original target species for each endpoint.

**New paired island outcome:** D(i,m) = gain(i,actinomorphic×deep,m) − gain(i,tubular/bell×deep,m). Compare endpoints on the **same physical island** (both endpoints require >=5 original directly documented target plants and >=10 matched mainland source donor species); group equal island means into original 10-degree geographical blocks or independent first source entities, then 999 resamples of entire groups. Strong support gate >=30 matched islands, >=10 blocks AND >=10 first source entities. **The two endpoints are not necessarily measured on the same set of plant species within a matched island.** This is NOT phenotype-provider holdout, whole-archipelago out-of-fold prediction, a bird/butterfly guild test, local flower evolution, or estimated ancestral mainland phenotype.

## ORIGINAL numerical results: matched island contrast D (Brier units)

| Original geography & scope | geo5 qualified common islands / blocks / first sources | geo5 block-equal D [95% conditional block-resampling] | climate10 qualified islands / blocks / first sources | climate10 block-equal D |
|---|---|---|---|---|
| Tropical, all recorded | **248 / 66 / 54** | **+0.157729 [+0.137795,+0.178638]** | **267 / 68 / 55** | **+0.157967 [+0.140072,+0.177370]** |
| Tropical, strict known native | **42 / 15 / 13** | **+0.100887 [+0.056608,+0.141904]** | **45 / 16 / 14** | **+0.111998 [+0.075477,+0.146886]** |
| Tropical, strict native/nonendemic | **41 / 15 / 12** | **+0.101945 [+0.053114,+0.146243]** | **44 / 16 / 13** | **+0.116048 [+0.072137,+0.157694]** |
| Northern midlatitude, all recorded | **1166 / 70 / 85** | **+0.049906 [+0.033545,+0.067441]** | **1189 / 72 / 87** | **+0.048473 [+0.033840,+0.063723]** |
| Northern midlatitude, strict known native | **139 / 25 / 29** | **+0.049816 [+0.027402,+0.076030]** | **143 / 27 / 31** | **+0.039769 [+0.018017,+0.061220]** |
| Northern midlatitude, strict native/nonendemic | **138 / 24 / 28** | **+0.034259 [+0.011913,+0.060418]** | **143 / 27 / 31** | **+0.023681 [+0.002507,+0.045565]** |

Tropical all recorded, >=2 distinct donor species per genus in BOTH mainland modes: geo5 **142 islands /47 blocks /40 first sources**, D **+0.146904 [0.116873,0.174007]**. Northern midlatitude same donor gate: geo5 **1002 islands /57 blocks /78 first sources**, D **+0.047328 [0.027611,0.068328]**. Thus same-island positive predictive asymmetry persists when strengthening mainland donor count, but regional sample and trait documentation change.

Equal **first-source group** sensitivity (not independent sample replication): tropical all-recorded geo5 D **+0.174000 [0.151612,0.196996]**, tropical strict native geo5 **+0.093550 [0.042524,0.138639]**, northern midlatitude all-recorded geo5 **+0.061881 [0.045132,0.077527]**.

Northern high latitude all-recorded has **89** paired geo5 floras, **17** geographical blocks but **only eight** first-source groups; southern extratropical **16** geo5 floras /six blocks/nine first sources. Both are explicitly **NOT_IDENTIFIABLE_COMMON_ISLAND_SUPPORT**, NOT zero morphological asymmetry. Geographic unresolved one flora has no paired support.


## NEW original same-island endpoint decomposition (2026-11-17)

Follow-up protected original source [run #44025078142](https://github.com/zuizui0223/island/actions/runs/44025078142) **SUCCESS** including synthetic aggregator tests; anonymous [artifact #14303140157](https://github.com/zuizui0223/island/actions/runs/44025078142), extracted `paired_anonymous.json` SHA256 `f517368ce36d83ff6d976fe137d205ac9c17903975f1057b3590194305c461db`. The source runner now independently exports both component Brier gains and their within-island difference, each with whole-geographical-block bootstrap intervals. This is a new original-data breakdown on the **same** geo5 matched-island cohorts; no additional new floral observations were collected.

| Same-island region/plant scope | Number of original supported islands | Gain for **actinomorphic×deep** genus vs pooled [95% block] | Gain for **tubular/bell×deep** genus vs pooled [95% block] | Their exact signed difference |
|---|---:|---|---|---:|
| Tropical all-recorded | **248** | **+0.115652 [+0.101069,+0.130290]** | **−0.042077 [−0.056185,−0.028884]** | **+0.157729** |
| Tropical strict native | **42** | **+0.087745 [+0.041415,+0.129707]** | −0.013141 [−0.030702,+0.001402] | +0.100887 |
| Tropical native/nonendemic | **41** | **+0.090416 [+0.043503,+0.132771]** | −0.011529 [−0.033350,+0.002642] | +0.101945 |
| Northern midlatitude all-recorded | **1166** | **+0.047145 [+0.031709,+0.064129]** | −0.002761 [−0.008887,+0.002921] | +0.049906 |
| Northern midlatitude strict native | **139** | **+0.037312 [+0.020111,+0.058001]** | **−0.012504 [−0.025733,−0.002233]** | +0.049816 |
| Northern midlatitude native/nonendemic | **138** | **+0.022280 [+0.004199,+0.042728]** | **−0.011980 [−0.025985,−0.001928]** | +0.034259 |

For original all-recorded tropical `climate10` same-island n267, positive actinomorphic×deep **+0.119269** and negative tubular/bell×deep **−0.038698**, difference **+0.157967**. For northern midlatitude `climate10` n1189, actinomorphic×deep **+0.047449**, tubular/bell×deep **−0.001024**, difference **+0.048473**. This demonstrates that the signed difference is **not merely an unexplained arithmetic subtraction of separate old cohort means**. The all-recorded tropical tubular/bell loss of source genus predictive skill is negative under both candidate source modes; it is less secure in the sparse strict native subset, where both signed intervals span zero.

**Claim:** real documented island flower morphologies show **phenotype-specific conditional genus predictability** relative to mainland candidate congener states. **Not claim:** shared target plant species across both endpoints, independent information-source validation, source provider out-of-fold generalization, adaptive selection, effective pollinator turnover, or cross-region contrast significance. Equal blocks and bootstrap intervals condition on original spatial/source group assignments; ecological historical processes and observational annotation biases are not thereby identified.

## Scientific interpretation and explicit claim ceiling

**Positive:** in independently validated original GIFT botanical source-name-excluded source-prediction comparisons, the association of real joint floral architecture with **observed island genus membership is phenotype-specific**, not a universally positive predictive help. The new same-island contrast eliminates *different island selection* as a complete explanation for the previous radial/deep versus tubular/deep sign difference, though endpoints can still have different sets of documented target species. Tropical effect is descriptively larger than northern-midlatitude, but a *formal interregional paired-source spatial test* was not performed, and this analysis was designed after seeing source-only pilot results. The tropical strict-native panel has 42 supported islands, not 4,505 fully eligible independent trait comparisons.

**Not established:** distinct pollinator guilds or effective services; pollen limitation; different island vs mainland population phenotype; evolution/selection after colonization; independent provider-holdout prediction; unbiased comparison under large original Direct trait documentation gaps. The old original tropical raw 331 vs 260 source-only endpoint comparisons are not a paired contrast and must not be reported as one. Even positive conditional intervals are not multiplicity-adjusted post hoc discovery p-values. Do not present this as an island-specific adaptation causal conclusion.

**Necessary next test:** re-evaluate on exactly common *plant species* with symmetry, form and depth all original Direct recorded and original source donors in both modes; then separately do provider/archipelago holdouts and documented native balance. Small triple-documentation sample is an auxiliary mechanism screen, not a replacement of the global 4,505 island macroecology population. The main article retains global functional beta and four-region real paired floral repertoire as primary analyses.

Sources: [original runner](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_original_paired_congener_gain.py), [protected workflow](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/.github/workflows/run-fd-original-paired-congener-gain.yml), [source-only precursor](../20261116_source_only_congener_real_joint_prediction/README.md). PR #16 remains DRAFT.
