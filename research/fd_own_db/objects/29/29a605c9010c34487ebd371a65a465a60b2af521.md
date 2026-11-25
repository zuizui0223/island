# DoPI v2 published-locality × year/month botanical visitor morphology audit — 2026-11-17

**Original protected GIFT/DoPI source CI SUCCESS.** This tests only positive **European** visitor associations of plant species which ALSO occur on the original GIFT **4,505** island floras; no original island-local flower visitation or pollinator absence is inferred. Even identical grid codes and published coordinates are not necessarily the same experimental flower patch.

## Provenance, reproducible sample and source access

- New source-protected [original island Actions #44092542811](https://github.com/zuizui0223/island/actions/runs/44092542811) **SUCCESS**; source executable pinned fd commit `38042cdcd5aafdb69220b7191a231e80f1f0d3b2` and original `island` source branch `analysis/fd-discrete-functional-volume-20261113`. Original anonymous artifact #**14333170566**, extracted `anonymous_original_dopi_gift_guild_motifs.json` SHA256 **774a7282b0d0f8df5cab81b84d60c463c14c72efb2917230d153e411aba296b1** (1,827,482 bytes). Source original GIFT Direct, GIFT botanical flora and published DoPI Version 2.0 raw 676,900-row CSV SHA guards verified inside runner.
- Original independent public data source: [DoPI Version 2.0 DOI 10.5281/zenodo.22876365](https://doi.org/10.5281/zenodo.22876365) and [Balfour et al. original Ecology 2022 database description](https://doi.org/10.1002/ecy.3801). Source exposes *publication, country, British National Grid or coordinates, year, month* when originally recorded.
- Scope check: **676,900** original DoPI records; **490,162** rows mentioning GIFT-island-listed accepted species with some directly described 1–4-axis flower phenotype; **112,754** of these lacked a usable combination of bibliographic, original site or year information and were not imputed; **320,563** source rows have supported geographic/year metadata, positive grade1–3 evidence and a supported insect visitor group. Original **301 publication keys** contribute site/year-evidenced source observations and **40,712** *positive site×guild×evidence×resolution units*. These 40,712 are NOT unique sites, studies or independent evolutionary events.
- [Pure source-matched site/year/month code](../../src/fd/dopi_site_period_guild.py), [synthetic exact-locality/year and pollen-evidence tests](../../tests/test_dopi_site_period_guild.py), [protected original runner](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_original_dopi_gift_guild_motifs.py) and [original workflow](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/.github/workflows/run-fd-original-dopi-guild-motifs.yml). No source study labels, original GIFT accepted plant, island, exact site or published locality are serialized into the public aggregate.

## What 'same location' and 'same time' actually mean

A source guild pair is considered co-recorded only when its positive flower interactions have the SAME original study identifier **AND** identical publication-provided British grid cell OR exact printed lat/lon coordinates with at least two decimal positions **AND** the SAME nonmissing observation year; stronger sensitivity also requires the SAME nonmissing month. Unknown site/time or a country-level centroid cannot manufacture a matched record. This is a stringent *published-coordinate equality criterion*, NOT proof that both guilds were surveyed simultaneously, that the exact same flower population was available, or that the original site metadata have sufficient precision to identify a field plot.

For every supported study-site-period and separately directly documented flower phenotype, both guild-positive plant sets must contain at least two species. Use study-year/month location-specific plant phenotype differences, average qualifying locations within their original publication, then bootstrap *whole publications* (999 draws, >=10 distinct source papers and >=10 genuinely motif-positive papers). A separate source-genus balancing sensitivity requires >=2 botanical genera with actual directly characterized flowers visited by BOTH guilds in that study-location-period. That stricter comparison typically has a DIFFERENT nested publication/plant cohort from the unrestricted comparison: the two means cannot directly be subtracted to estimate a matched within-study causal selection effect.

Across 5 visitor guild pairs ×27 real 1/2/3/4-axis Direct flower outcomes ×2 positive evidence tiers ×2 time resolutions = **540 reporting profiles**, heavily correlated, no family-wise multiple testing correction.

## How much evidence survives localization?

| Source site/time and original evidence scope | Direct profile/guild contrasts with sufficient publication + TRUE positive motif records (of 135) | Same-location profile contrasts with >=2 shared botanical genera (of 135) |
|---|---:|---:|
| Same publication × published grid/coordinates × **year**, positive visitor grades 1–3 | **36** | **30** |
| Same site × **year and month**, positive visitor grades 1–3 | **32** | **23** |
| Same site × year, **pollen carriage / confirmed grades 1–2** | **4** | **2** |
| Same site × year and month, **pollen carriage / confirmed grades 1–2** | **0** | **0** |

The site/year/month study comparability is a genuine improvement in observational control. Yet no locally matched grade1+2 pollen-evidence profile is sufficiently supported with both month and place; true effective pollen transport on island plant populations is NOT identified.

## External European visit pattern remains in source-matched same-site/year/month comparisons

Differences are **recorded floral morphology frequencies among positive guild-visited plants** (guild A minus guild B); NOT relative visitation preference after measuring floral supply/abundance, native island pollinator occurrence, extinction or seed/fruit success. Original 95% uncertainty is whole-publication conditional bootstrap, NOT multiplicity adjusted.

| Same documented place × year × month guild comparison | Original matched papers | Visitor-positive plant flower morphology frequency difference, pp [95% source-study interval] | Separate shared-genus-equal nested sample |
|---|---:|---|---|
| **Bombus − hoverflies**, zygomorphic | **64** | **+23.841 [+19.369,+28.011]** | 56 papers, **0.000 [0,0]** pp; genus mostly fixed |
| **Bombus − hoverflies**, open/radial | **64** | **−21.887 [−26.096,−17.646]** | 57 papers, +0.088 [0,+0.263] pp |
| **Butterflies − Bombus**, zygomorphic | **43** | **−8.655 [−14.300,−3.148]** | 36 papers, 0.000 [0,0] pp |
| **Butterflies − Bombus**, open/radial | 43 | **+6.426 [+3.088,+10.467]** | 37 papers, −0.450 [−1.351,0] pp |
| Other bees − Bombus, tubular/bell | 87 | **−5.179 [−7.401,−2.748]** | 76 papers, +0.500 [−0.030,+1.339] pp |

The older source-published study-only fully equal-genus test [here](GENUS_EQUAL_WITHIN_STUDY_AUDIT.md) demonstrates that same-genus flower morphology rarely varies enough for genuine within-lineage guild inference; here apparent exact 0.000 in fixed-genus outputs is mostly a **structural no-variation limit, NOT proof of ecological equivalence**. Because site-year/month matched-genera panels select a smaller cohort than the raw guild panels, do not calculate a causal attenuation percentage by subtracting these two columns.

**Same published location and year (month not always provided):** Bombus−Syrphidae zygomorphic +20.150 pp, 98 papers; open/radial −19.513 pp, 95 papers; butterflies−Bombus zygomorphic −6.755 pp, 65 papers. The direction is similar to month-restricted cohorts, not an independent replication.

## Main article interpretation and stop rules

1. **Source-supported positive finding:** recorded floral access structures (especially symmetry versus open radial architecture) associate with the distribution of documented positive visitor guilds in European literature even when comparing the same original study and its reported locality/year/month. This reinforces that floral functional traits can have pollinator-access relevance.
2. **Mechanism NOT identified:** genus-fixed morphology largely cannot vary in these observations; neither effective pollination on GIFT islands nor historical guild changes nor true pollinator absences were measured. GIFT flowers are species-wide morphological descriptions, not observed evolution on island populations.
3. **Competing hypothesis:** tropical source-relative GIFT radial/deep versus zygomorphic/deep flower-frequency shifts may reflect **sorting of ancestrally different botanical lineages and generalized multi-guild access** rather than actual bee-to-butterfly or bird transitions. Original eight-motif cross-source agreement was only 5/8 overall and its source-category label shuffle did not establish a generalizable syndrome correspondence. No channel-shift claim from these positive-only European visitor records.
4. **Next biological qualification:** need contemporaneous true island/native mainland visitor and floral community availability, independently observed adequate sampling effort, UV/visual traits, scent/anthesis and floral reward, plus pollen-transfer effectiveness. A verified source study×coordinate×year×month does not substitute for those measurements.

Do **not** shrink the original 4,505 GIFT global flora study to these positive externally documented European visitor subsets. Main paper remains global botanical taxonomic assembly + discrete functional beta + regional real paired flower architecture + exact signed source-genus accounting; DoPI is a carefully scoped external biological compatibility / identification ceiling.
