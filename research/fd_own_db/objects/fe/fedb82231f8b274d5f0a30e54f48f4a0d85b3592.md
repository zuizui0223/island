# DoPI v2 × original GIFT: within-publication **equal-genus** pollination-mode test (17 November 2026)

**Completed original source-protected test, exploratory and descriptive.** This is a critical correction to interpretations of the earlier `shared_genus_study` comparisons: merely retaining genera appearing in both visitor guilds DOES NOT fix their relative genus weights. Here every shared genus is weighted equally **within each source publication**, then each publication is weighted equally. A second, stricter test requires actual **within-genus variation in the flower motif**. Original botanical study population: all 4,505 GIFT islands. Source DoPI observations are EUROPEAN, not actual island visits.

## Original source verification

- [Protected original GIFT×DoPI v2 CI #44090688441](https://github.com/zuizui0223/island/actions/runs/44090688441) **SUCCESS**, with tested fd equal-genus model pinned at commit `9c80fce1193b4f6085faab0863d4f1bcb85a007a`. Anonymous artifact #**14332460249**; exact extracted `anonymous_original_dopi_gift_guild_motifs.json` SHA256 **0c75d35e9a2ec6bbaa788a0bcfe8724bf874efae363f3e1dfd1054d99ab2637c** (1,317,682 bytes).
- 4,505 GIFT accepted island floras retained; exact original GIFT SHA256 `6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a`, high/medium species-Direct floral evidence SHA256 `9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745`. Original DoPI v2 676,900 rows, exact CSV SHA256 `b95d3084399ed2c0978af2485801108a27d14e843b9fed2611934ea582b745c8`. Exact GIFT accepted species × DoPI matches **3,203** plants across 518 distinct original source-publication keys.
- 4 study evidence tiers × 5 pre-specified visitor guild pairs × 27 directly species-scored, nonexclusive 1/2/3/4-axis flower motifs = **540 panels**, nested and correlated, never 540 independent experiments. Only positive guild-source records grade 1–3 used for main contrasts. All source plant, citation, place and site identifiers remain inside protected original Actions. [Pure implementation](../../src/fd/dopi_study_genus_equal.py) and [synthetic source and Simpson-confounding controls](../../tests/test_dopi_study_genus_equal.py); [original runner](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_original_dopi_gift_guild_motifs.py).

## Three EXACT-matched comparisons

For every source paper reporting BOTH guilds with >=2 Direct-scored accepted plants per guild, and at least TWO botanical genera documented by BOTH guilds, compare exactly the same study population:

1. `raw_study_diff` = proportion of positively documented visited plants carrying motif in guild A minus guild B, over all Direct-characterized study plants.
2. `shared_genus_pooled_diff` = same, after limiting to genera observed with BOTH guilds, still weighting each genus according to its own number of plants in each guild. This was the previous genus "conditional" score; **not truly genus fixed**.
3. **`genus_equal_diff`** = for each mutually recorded genus in a study, compute the guild A−B motif prevalence within that genus; **average each genus equally within study**, and then average studies equally. This completely removes differing genus mixture weights on this *selected observed* cohort, not species sampling bias or locality/time confounding.

Additionally, `genus_equal_informative` restricts to those studies whose shared genera contain BOTH a positive-motif and negative-motif species, so within-genus species-level motif association is actually identifiable rather than mechanically zero. All estimates are conditional study-level 999-resample bootstrap intervals, not preregistered or family-wise corrected tests. Strong support: >=10 publications **and** >=10 studies with >=1 actually documented motif-positive plant in the paired-guild cohort. A null/zero phenotype support count is **NOT_IDENTIFIABLE**, never evidence of ecological absence.

## Main source-authenticated result: classic guild-flower signal is overwhelmingly in genus composition

Original DoPI (grades 1,2,3; actual visitors), same comparable source publications for all three columns. Scores are **percentage points (pp)** of morphological trait frequency difference, *not visit probability or effective pollen transfer*.

| Observed guild A − guild B, true Direct flower morphology | No genus-weight control: study equal | Shared genera only, different genus weights retained | **Within-study genus-EQUAL** [95% whole-study interval] | Papers with same-genus observations / actually informative within-genus floral state variation |
|---|---:|---:|---|---:|
| **Bombus − Syrphidae, zygomorphic** | **+15.007** | +3.248 | **+0.027 [−0.105,+0.213]** | **123 / 4** |
| **Bombus − Syrphidae, open radial** | **−16.571** | −2.827 | **−0.079 [−0.247,+0.072]** | **119 / 10** |
| **Butterflies − Bombus, zygomorphic** | **−5.006** | −0.901 | **+0.021 [0,+0.063]** | **91 / 1** |
| **Butterflies − Bombus, open radial** | **+4.320** | approximately 0 | **−0.144 [−0.480,+0.041]** | **88 / 4** |

**Interpretation ceiling:** changing genus frequencies of *documented visited plants* explain most of the raw study-level differences. But genus states can be nearly invariant, so the equal-genus score is then approximately 0 **by construction**: it is an *identification limitation*, NOT a demonstration that effective pollinators have no within-lineage role. The European DoPI source is source-positive only, has no true non-visitation controls and often the same plant appears in both guild-positive lists. Survey, study, site, year, floral abundance, host availability and publication heterogeneity remain uncontrolled. Nor do these positive European visitor data identify an island-local Bombus-to-butterfly shift.

### Source-denominator qualification: how many comparisons actually contain within-genus variation?

| Source evidence tier | Genus-equal publication-supported profiles / 135 | Profiles with >=10 original publications containing **within-genus variable flower states** /135 |
|---|---:|---:|
| Positive visitation or pollen grades 1–3 | **37** | **22** |
| Direct observation with both surveys coded `all` | **31** | **13** |
| Pollen-carriage/confirmed grade 1–2 | **8** | **0** |
| Confirmed effective pollination grade 1 | **0** | **0** |

Among those 22 positive-visit within-genus-variable supported profiles, **8 are single-axis and 14 two-axis; ZERO are three/four-axis nonexclusive profile combinations**. There is NO original source-supported within-genus variable *multi-axis bird/butterfly/hawkmoth-specific* floral syndrome evaluation.

### A few *exploratory*, selected positive within-genus examples—not confirmatory

- On the **within-genus informative study subset** (selected after original data screening), `Bombus − Syrphidae` **tubular/bell** is **−14.618 pp [−22.881,−5.951]**, 31 publication keys /57 genuine shared-genus study units. Direct observed `all`-survey sensitivity n18 gives −13.981 [−25.926,−3.796]. However conditioning on within-genus variable floral state changes the study sample and selecting among 540 nested panels inflates apparent evidence. The full genus-equal, fixed all-source-study estimate is just −0.442 [−0.747,−0.177] pp on 119 sources, dominated by genera with no recorded variation. Neither estimates an effective pollination service or causal guild preference.
- `other_bees − Bombus` **actinomorphic × tubular/bell** on n24 informative papers is **+25.434 pp [14.058,37.261]**, with overlapping direct-visit `all` n16 +26.736 [12.153,43.403]. In the all-genus-equal n160-paper subset it is only **+0.521 pp [0.152,1.000]**. This is a hypothesis-generating **within-bee-guild** difference in the positively *documented* species roster. It does not support a specific island syndromic transition, and is post hoc without 540-panel multiple testing correction.

## Most important biogeographic implication

The earlier original GIFT source-only genuine same-plant floral syndrome-composition changes are TRUE descriptive plant assemblage results. However external flower guild compatibility **largely collapses after genus fixed weighting** for the classical Bombus–Syrphidae and butterfly–Bombus architecture contrasts; many genus-specific flower phenotypes do not vary adequately to identify a pollinator-mediated difference. The result therefore favors describing **lineage-mediated functional repertoire reassembly** and **potential multiple-guild access**, not an independently validated island-adaptive bird/butterfly/bee syndrome shift.

The related earlier GIFT-vs-DoPI **[full eight-profile sign concordance audit](EIGHT_MOTIF_SELECTION_AUDIT.md)** found only 5/8 same directions, with 3 source/visitor same-sign intervals and 1 discordant (yellow_orange × tubular/bell), and no generic significance. The present strict equal-genus control makes those potential structural links even less diagnostic of active pollinator selection.

## Next real-data discriminating test

The actual original DoPI v2 columns include country, latitude, longitude, British grid fields, year, month and date alongside publication, methodology and survey scope ([DoPI source documentation](https://www.sussex.ac.uk/lifesci/ebe/dopi/about), [original Ecology data paper](https://doi.org/10.1002/ecy.3801)). Restrict visitor-guild pair comparisons further to the **same study × documented site × observation year/month**, without inventing site IDs for missing coordinates, and then match directly measured morphologies and genus. Report independent study/sample-site counts *before* claiming functional visitor differences. True absence, effective pollen transfer, and within-island vs ancestor morphology still require other sources/experiments.

**New source-code follow-up:** [Published-location × observation-year/month source qualification](../../src/fd/dopi_site_period_guild.py) with [synthetic tests](../../tests/test_dopi_site_period_guild.py) is now wired into the protected original [GIFT×DoPI workflow](https://github.com/zuizui0223/island/actions/workflows/run-fd-original-dopi-guild-motifs.yml), via source branch commit `70948978cd906d85b3368b6408748fec15962f70`. Treat its site/year/month counts as **unknown until the source run succeeds and the anonymous artifact is audited**. Site grid codes or matching printed lat/lon do not necessarily prove a single experimental plot; time metadata are observational.

**Paper decision:** GIFT 4,505 botanical floras, functional beta and relative regional real floral morphologies remain mainline; DoPI as external visitor *source-positive access compatibility* with substantial lineage confounding. Never use 17D uncalibrated KDE volume as an ecological syndrome estimate, and keep PR #16 DRAFT.
