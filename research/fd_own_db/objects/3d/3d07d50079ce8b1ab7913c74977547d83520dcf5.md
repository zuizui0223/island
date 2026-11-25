# 2026-11-16 — Original 4,505 GIFT island same-island geo5 vs climate10 floral source robustness

**Completed protected-original-source numerical execution**: [island Actions #43845535858](https://github.com/zuizui0223/island/actions/runs/43845535858), SUCCESS after original botanical/trait source hash checks, independently rebuilt geography and guild acquisition, and source-mode/guild synthetic tests. Original anonymous artifact [#14214953065](https://github.com/zuizui0223/island/actions/runs/43845535858/artifacts/14214953065), extracted JSON SHA256 `5f626d9ede4adb5b238967df9cc11a00d5214098ce34080c935ffffbfc70c52f` (434,609 bytes). Original island/species/GIFT entity/GBIF row IDs remain in protected CI, never in public output. Original [analysis runner](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_source_guild_source_mode_robustness.py), [fd statistics](../../src/fd/source_mode_robustness.py), [synthetic tests](../../tests/test_source_mode_robustness.py).

## Frozen comparison, population and interpretation

- Preserve all **4,505** original GIFT botanical island floras, including incomplete morphological, native-status and source-availability records in the population denominator.
- Compare TWO independently archived candidate native mainland region sets: original **geographical nearest five** (`geo_k5`, 41,325 GIFT assignment rows across the larger original island source frame) and **climate-nearest ten of geographic fifty** (`geo50_climate10`, 82,190 rows). Climate source membership is available to **4,485/4,505 botanical island floras**. They contain *different candidate accepted plant names* and are not known dispersal or actual ancestral populations.
- At every original focal island, **exclude every name from its complete botanical island flora from BOTH mainland candidate species pools**, calculate direct species-level floral sector frequency in each pool, and keep the exact SAME island plant members across modes. Three NESTED provenance scopes: all-recorded; strictly status-known native; strictly status-known native nonendemic. No missing origin status is imputed native.
- Five original high/medium truly Direct morphological sectors only: zygomorphic symmetry, tubular/bell form, open-radial form, separately coded **deepest tube**, and genuinely SAME species zygomorphic×tube/bell paired morphology. Requires >=5 island Direct carrier species and >=10 candidate-only mainland Direct carrier species under **BOTH modes**, separately for each sector. No requirement for three-trait complete cases or invalid 17-dimensional KDE.
- Outcome `island phenotype share − candidate-source Direct phenotype share`. Pairwise sensitivity `(island − climate10) − (island − geo5) = geo5 candidate share − climate10 candidate share` **exactly by algebra**, not an independent replication or climate causal effect. All confidence intervals below are *descriptive spatial-block bootstrap sensitivities*, not taxon/phylogeny/propagule or causal source selection adjustments. 999 whole original spatial-block resamples.
- An observational Direct-documentation sensitivity retains only SAME islands where original island proportion with trait annotation differs from **each** candidate mainland source mode by <=10 percentage points; no missing plant states are filled.

## Primary same-island direct morphology contrasts: all-recorded

Numbers are **percentage-point difference island − candidate source**, geography 10-degree-block equal conditional means and bootstrap 95% intervals. n counts islands in BOTH candidate-source mode direct-support scopes.

| Region, sector | Same islands | geo5 pp [95%] | climate10 pp [95%] | Same-island sign agreement | n island/source documentation balanced within 10 pp |
|---|---:|---|---|---:|---:|
| Tropical, zygomorphic | 945 | **−8.37 [−10.21,−6.63]** | **−8.58 [−10.21,−6.84]** | 93.8% | 106 |
| Tropical, tubular/bell | 773 | **+3.20 [+1.26,+5.10]** | **+3.73 [+1.62,+5.61]** | 92.4% | 304 |
| Tropical, deepest tube | 728 | **+9.91 [+7.96,+11.75]** | **+10.00 [+8.15,+11.75]** | 94.0% | 338 |
| Tropical, truly same-plant bilateral×tube/bell | 692 | +0.99 [−0.43,+2.58] | +1.77 [+0.37,+3.32] | 90.8% | 386 |
| N-midlatitude, zygomorphic | 1,750 | **−4.99 [−6.19,−3.91]** | **−7.12 [−8.34,−5.98]** | 92.9% | 213 |
| N-midlatitude, deepest tube | 1,537 | **+1.89 [+0.38,+3.43]** | −0.55 [−2.10,+1.04] | 91.9% | 1,091 |
| S-extratropical, tubular/bell | 157 | **−12.59 [−17.75,−7.12]** | **−12.58 [−18.05,−6.87]** | 92.4% | 31 |
| S-extratropical, open-radial | 157 | **+8.74 [+3.88,+13.19]** | **+10.76 [+5.99,+15.01]** | 90.4% | 31 |
| N-highlatitude, zygomorphic | 258 | **−15.17 [−18.58,−12.02]** | **−15.68 [−19.11,−12.47]** | 99.2% | 108 |

**Nontrivial falsification:** the northern-midlatitude deepest-tube positive source-relative direction under geo5 is **not stable** under climate10: the estimate changes from +1.89 pp to −0.55 pp (the latter interval spans zero). Thus not all structural axes are source-choice stable. Southern open-radial appears positive in all-recorded, but among the **same 31 islands satisfying both source documentary coverage gates**, the source-relative differences are **−0.13 [−5.14,+5.84] geo5** and **+1.57 [−2.79,+6.67] climate10**, so the broad positive cannot be interpreted as provenance/documentation-robust ecological selection. For S-extratropical bilateral symmetry, only 24 islands pass the coverage-gap gate, below the >=30 geographical support minimum.

**Positive descriptive regional pattern that survives the source choice:** tropical zygomorphic proportion is lower, while true deepest-tube and tubular/bell proportion are higher; a single generic "specialist flower loss" cannot capture those opposing directions. With <=10pp documentation-gap in tropical all-recorded, zygomorphic remains **−7.48 [−11.57,−3.85] vs −8.22 [−11.85,−4.87]** (106 islands); deepest tube **+7.84 [+5.50,+10.14] vs +7.60 [+5.36,+9.83]** (338); tubular/bell **+6.98 [+4.51,+9.29] vs +7.89 [+5.40,+10.14]** (304). Still observational literature-coded trait patterns, not lineage-independent adaptation.

## Native provenance: informative, nested, strongly incomplete

Same island and same morphologically scored species within each mode, but **all-recorded vs strict-native are NOT independent replicates**. Strict-native and native-nonendemic eligibility is a nonrandom original provenance/documentation subset.

| Strict plant scope, region, sector | Same islands / blocks | geo5 pp [95%] | climate10 pp [95%] | <=10pp original trait coverage-balanced islands |
|---|---|---|---|---:|
| Native, tropical zygomorphic | 125 / 42 | −6.43 [−8.94,−3.69] | −6.31 [−9.03,−3.59] | **12** (6 blocks: not estimable) |
| Native, tropical tubular/bell | 112 / 40 | +6.23 [+3.31,+9.68] | +6.99 [+3.90,+10.61] | 55 |
| Native, tropical deepest tube | 107 / 40 | +9.90 [+4.69,+15.26] | +9.40 [+4.45,+14.54] | 53 |
| Native, tropical real bilateral×tube/bell | 102 / 39 | +4.15 [+1.80,+6.50] | +4.78 [+2.53,+7.05] | 67 |
| Native nonendemic, tropical zygomorphic | 120 / 42 | −6.60 [−9.42,−3.26] | −6.49 [−9.65,−2.97] | **9** (5 blocks: not estimable) |
| Native, N-midlatitude deepest tube | 207 / 35 | −2.78 [−4.95,−0.54] | −4.98 [−6.95,−2.88] | 176 |
| Native nonendemic, N-midlatitude deepest tube | 207 / 35 | −3.58 [−5.72,−1.53] | −5.78 [−7.76,−3.88] | 167 |

**Interpretation:** Native-only northern-midlatitude deepest tube is lower than candidate source under both source definitions, whereas all-recorded geographic-nearest5 gives the opposite direction in a DIFFERENT island cohort (1,537 vs 207). This does NOT establish an introduction effect: missing nativeness changes both source islands and taxonomic composition. In the full strict-native high-latitude panel, zygomorphic has only 11 islands across 3 blocks; S-extratropical native zygomorphic only 27 islands. Their region-level source differences are **NOT_IDENTIFIABLE** under the frozen >=30-island/10-block gate. Highlatitude all-recorded -15pp cannot be interpreted as a demonstrated strictly native or phylogenetically repeated ecological loss.

## Same-island qualified guild effect falsification across source definitions

Original independent GBIF channel states `retained` and `disrupted` retain their strict meaning: documented available and source-available, sufficiently surveyed but NOT documented. `disrupted` does not indicate measured past extinction or historical guild change.

- Across **3 nested plant scopes × 5 regions × 5 guilds × 5 sectors = 375** support slots, only **3** pass the source/effort/physical-block gates, all in **TROPICAL × FLOWER-VISITING BIRDS × ALL-RECORDED**, same previous three correlated sectors. **372** fail. Strict native has **0** support-passing guild contrasts, not a null ecological effect.
- On the exact common 2-mode island cohort, zygomorphic original retained/adequately-nondetected **587/13** (10 mixed blocks), tubular/bell and open-radial **500/12** (9 mixed blocks). All are the same underlying bird guild, not independent visitor-system replications.
- Source-relative bird-retained minus bird-nondetected matched physical-block contrasts in percentage points (95% descriptive bootstrap):

| Tropical structural sector | geo5 pp [95%] | climate10 pp [95%] |
|---|---|---|
| zygomorphic | +2.14 [−8.11,+9.50] | +2.84 [−7.22,+10.70] |
| tubular/bell | −3.58 [−14.29,+6.11] | −3.25 [−12.95,+6.23] |
| open-radial | +1.34 [−5.23,+7.91] | +0.92 [−5.54,+7.75] |

All cross zero; comparing two similar GBIF proxies does not demonstrate true species-specific bird pollination, bee-to-bird shift, compensation, independent selection or fitness. Neither source mode identifies a guild mechanism. No new guild×source interaction coefficient was estimated.

## Scientific decision

- **Strongest new defensible description:** the tropical recorded flora has a directionally heterogeneous source-relative architectural composition (lower zygomorphy but higher deep tubes and tubular/bell forms), reproduced under two substantially different independently specified mainland candidate pools and partially retained in strictly native-known flora.
- **Still unresolved:** whether true taxonomic recruitment/colonization filtering beyond finite candidate-source genus composition and **differential floral Direct documentation** generates this composition. Earlier GIFT source-breadth finite-genus null explains much of tropical regional bilateral and tube/bell patterns; that rival remains intact, not eliminated by the present source-mode comparison.
- **Unsupported:** universal island floral simplification, invariant highlatitude flora, island-only evolution of morphology, and pollinator availability as an identified causal determinant.
- **Next falsification:** on EXACT geo5/climate10 same-island directly scored native-source cohorts, condition on original fixed source listing breadth and genus identity and audit positivity of true source-genus overlap; separately validate detailed plant–visitor interactions, preferably time-, study- and plant-identity-aligned. Do not repair missingness by guessing trait states or guild extinctions.

**Both fd PR #16 and island analysis PR remain Draft.** Analyses here are posthoc exploratory, not registered before seeing existing sector and source results.
