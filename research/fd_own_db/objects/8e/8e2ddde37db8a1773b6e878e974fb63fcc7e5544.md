# Full eight-motif GIFT×DoPI pollination-syndrome compatibility / selection audit — 17 November 2026

> **Follow-up taxonomic-identification audit:** The [source-verified equal-genus DoPI analysis](GENUS_EQUAL_WITHIN_STUDY_AUDIT.md) shows much of the visitor-guild versus floral profile association can result from different **weights of botanical genera within the same publication**; the previous shared-genus restriction was not genus-frequency standardization. The original eight motif sign test is not a causal or lineage-independent pollination syndrome result.

**Status: COMPLETED source-verified original anonymous aggregates; exploratory, NOT a validated guild shift or effective island pollination.** This audit **supersedes use of the three selected 3/3 direction concordances** as evidence of transferability. Original GIFT parent population remains 4,505 accepted botanical island floras; the animal evidence is European DoPI source-positive visitation of the same accepted botanical species, not visits observed locally on the island.

## Source verification and reproducibility

- **[Successful complete source-checked fd CI #44089006636](https://github.com/zuizui0223/fd/actions/runs/44089006636)** checked source-original anonymous GIFT true two-trait floral architecture artifact **[island CI #43857181963](https://github.com/zuizui0223/island/actions/runs/43857181963)** and external real DoPI v2 visitors **[island CI #44082152065](https://github.com/zuizui0223/island/actions/runs/44082152065)**. All sources were independently SHA256 checked from the public anonymous artifact, no original accepted island, species, study or locality identifiers exported. Final new anonymous artifact **14329574789**, source JSON digest **48b53133deac93a0808c2a2c122eb9b68320f49a7d26704a6ebde92a90a94e27** (120,776 bytes).
- Source original full botanical flora SHA256 `6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a`; GIFT Direct flora phenotype SHA256 `9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745`; DoPI v2 published CSV SHA256 `b95d3084399ed2c0978af2485801108a27d14e843b9fed2611934ea582b745c8`.
- [Reusable exact 8-profile audit module](../../src/fd/cross_source_syndrome_concordance.py), [noncausal null and source receipt regression tests](../../tests/test_cross_source_syndrome_concordance.py) and [source-verified original anonymous CI workflow](../../.github/workflows/fd-cross-source-syndrome-concordance.yml). Original 11 true Direct two-trait GIFT profiles vs DoPI 27 true 1/2/3/4 trait visitor profiles; **eight** original two-trait names are eligible with >=30 GIFT island floras / >=10 source geography blocks and >=10 European same-paper BOTH-guild positive morphology studies.
- Original 4 GIFT regions ×3 nested botanical provenance statuses ×2 frozen original potential mainland source definitions × eight jointly supported motifs = **128 distinct scenario×motif aggregate comparisons**, representing overlapping species and islands, not 128 independent evolutionary events. Detailed file includes 24 regional/status/source summaries.

## Original tropical full-recorded source-relative floral composition and European visitors

GIFT column = original tropical **island minus potential GIFT native mainland geo5 source** documented same-accepted-species true TWO-trait frequency difference. DoPI column = original EUROPEAN **butterfly-family positive visitor plants minus Bombus positive visitor plants**, within original papers positively documenting BOTH guilds; repeated studies are collapsed before study-equal means. Both values are in percentage points but have **DIFFERENT population, observation and functional estimands**; compare signs cautiously, NEVER subtract the columns as if they were same causal effects. All intervals below are original conditional 95% block/study bootstraps, post hoc and not multiplicity-adjusted.

| True species-level conjunction | GIFT tropical original island−source (pp), geographic 95% | DoPI Europe butterfly−Bombus (pp), publication 95% | Visitor study comparisons | Same sign? |
|---|---:|---:|---:|---|
| actinomorphic × directly deep tube | **+10.989** [+8.958,+13.313] | **+2.084** [+0.127,+4.554] | 70 | Yes |
| zygomorphic × directly deep tube | **−1.368** [−2.680,−0.050] | **−3.073** [−4.895,−1.097] | 70 | Yes |
| tubular/bell × directly deep tube | **+5.172** [+3.82,+6.61] | **+3.313** [+0.595,+6.439] | 67 | Yes |
| actinomorphic × tubular/bell | +4.483 [positive GIFT block interval] | +0.845 [−2.787,+4.328] | 89 | Yes (DoPI interval crosses zero) |
| zygomorphic × tubular/bell | +0.991 [GIFT block interval crosses zero] | +0.140 [−2.188,+2.604] | 89 | Yes (both weak) |
| red_pink × tubular/bell | −1.398 [negative GIFT block interval] | +1.435 [−0.416,+3.697] | 60 | **No** |
| **yellow_orange × tubular/bell** | **+3.714** [positive GIFT block interval] | **−4.282** [−7.639,−0.809] | 60 | **No; both intervals exclude zero** |
| white × tubular/bell | +0.413 [weak GIFT block interval] | −0.184 [−1.507,+1.371] | 60 | No (both weak) |

**Audit result: 5/8 signs agree, NOT a general 3/3 external validation.** The five structure-only conjunctions are descriptively 5/5 same sign; three **colour×form** conjunctions are 0/3 same sign. This contrast was noticed AFTER inspecting the data, motifs share accepted species and axes, and several component intervals cross zero: it is a potential **floral mechanical-access vs visual-signal ecological hypothesis**, not a confirmed syndrome-level ecological law or independent comparison of five versus three independent datasets.

**Most defensible evidence:** three structure/depth profiles each have both GIFT and DoPI original marginal uncertainty intervals of the same nonzero sign (actinomorphic/deep, zygomorphic/deep, tubular/deep). The yellow_orange/tubular profile has both original intervals NONZERO in OPPOSITE directions. Thus a strict count of individually supported marginal intervals is **three aligned versus one discordant**; the other four directions are uncertain. No family-wise correction.

## Full-eight category-selection / posthoc bias benchmark

Across the eight distinct eligible motifs: observed **Spearman rho +0.476**, Pearson r **+0.395**. Exact exhaustive **8! = 40,320 source-positive profile-label permutations** holding the two sets of eight values fixed gives:

- **P(at least five matched signs in eight labels) = 18/28 = 0.642857** under this *hypothetical shuffled-label benchmark*. In other words the observed **5/8** cannot establish a generic ecological alignment relative to this diagnostic.
- Absolute Spearman rank-correlation **tail fraction 0.243056** under the same 40,320 label shuffles. The original 2026 DoPI dataset, GIFT species and floral motifs were not independently sampled per trait.
- **CRITICAL:** neither fraction is a valid ecological/phylogenetic hypothesis p-value: flower characters are overlapping, the eight motifs were selected by data support and the access-vs-colour partition interpreted after results. Never advertise the visually striking 5/5 vs 0/3 partition as an independent chi-square/Fisher test.

## Geography, source & nativeness sensitivities

- The **tropical 5/5 structure versus 0/3 colour-form sign pattern** recurs as an *overlapping descriptive sensitivity* under BOTH GIFT geo5 geographic and climate10 candidate mainland pools and in all-recorded, strict-known-native and strict-native-nonendemic GIFT subsets. The native subsets are nonrandomly documented and have smaller and shared original island populations; source/trait biases remain.
- Northern-midlatitude **geo5 all-recorded**: 4/5 structure and 0/3 colour sign agreement; **climate10 all-recorded**: 4/5 structure and 2/3 colour. Therefore the apparent access/signal division is **not a universal worldwide rule**.
- Source-potential mainlands are not observed propagule source communities, and geographic 10-degree blocks are not independent historical shifts.
- DoPI within-shared-genus same-study morphology associations are greatly attenuated for **zygomorphic/deep** (−3.073 pp raw versus −0.255 pp genus-conditioned, CI overlaps zero) and **tubular/deep** (+3.313 pp raw versus +0.501 pp genus-conditioned, CI overlaps zero). Actinomorphic/deep remains +1.677 pp within shared genera, CI [+0.318,+3.187]. For yellow_orange×tubular/bell the raw visitor −4.282 pp becomes ~+0.014 pp within shared genera (interval zero). Real GIFT source-genus signed accounting and mainland species-donor LOO negative controls likewise show a substantial inheritance/lineage source component. Shared genus matching still does not establish independent site-level pollinator selection.

## What is and is NOT identifiable

**Yes:** species-wide directly documented flower **access architecture** can covary with which European plant taxa are positively visited by butterflies vs Bombus, and the same kinds of paired floral states have source-relative, regionally varying distributions among the original 4,505 global island floras.

**No:** bird vs butterfly/hawkmoth identity on tropical islands, effective pollen transfer/fruit/seed production, genuine local island visitor community replacement, source ancestry, within-species postcolonization flower evolution, causal guild sorting, or independent evolutionary replication. Original DoPI is **European**, source positive-only, literature & observation effort heterogeneous, and original GloBI can contain overlapping papers. GIFT actual accepted species-wide Direct colour **red_pink/yellow_orange** cannot resolve red from pink, yellow from orange, bee visual colour space, UV, nectar, scent or day/night flowering.

**Practical manuscript decision:** keep 4,505 global botanical island floristic/functional diversity and source-genus reassembly PRIMARY. Present external DoPI butterfly/Bombus as a **visitor compatibility calibration with a falsifying colour×form counterexample**, not as demonstrated island bee→butterfly shift or pollinator syndrome evolution. Follow-up requires common-site/period visitor comparisons and independent native morphology/source provider audits, ideally actual pollen transfer, but absence/underdocumentation must never be converted to missing guilds.

The older [three-profile direction note](TROPICAL_EXTERNAL_DIRECTIONAL_CONCORDANCE.md) is **superseded for inference by this source-complete all-eight audit**; its 3/3 subset remains mathematically correct but SELECTED.
