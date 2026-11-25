# Original GIFT × DoPI v2: same accepted plant visited by multiple guilds, within reported site/year/month — 2026-11-17

**Original source-protected quantitative analysis COMPLETED.** This tests whether true floral morphological profiles appear among plants with **POSITIVE documented visits from BOTH guilds** in the SAME original published study × site/grid × observation year/month, rather than inferring exclusive pollination syndromes from flower colour/shape. Original global GIFT parent study population is **all 4,505 botanical island floras**; the animal visits are independent European DoPI literature observations of species that are also listed on those islands, **not visits observed on the GIFT islands**.

## Reproduction and eligibility

- Original [island protected source CI #44098890581](https://github.com/zuizui0223/island/actions/runs/44098890581) **SUCCESS** on branch `analysis/fd-discrete-functional-volume-20261113`, workflow source commit `b0ce1d02a131e86d9e4803c5ea81be33500dc253`, reading tested `fd` model commit `49c83690d40ff3e1e815d6ddaf4fb56808eff651`. Source CI checks original GIFT accepted flora SHA, original real Direct trait ledger SHA, and immutable DoPI v2 data SHA. The four `fd` source workflow CI jobs at that model commit also **all SUCCESS**.
- Original anonymized [artifact #14334842586](https://github.com/zuizui0223/island/actions/runs/44098890581), extracted `anonymous_original_dopi_gift_guild_motifs.json` **SHA256 `18859a2dffa25741e65f6f6092c459e61e1676130a6b439a10df060cd34e23cf`**, 2,625,148 bytes. Original accepted species names, GIFT source IDs, published DoPI study identifiers and site coordinates are not in the exported JSON.
- Source original DoPI Version 2.0 DOI [10.5281/zenodo.22876365](https://doi.org/10.5281/zenodo.22876365), **676,900** original raw records; **301** distinct source-publication strings with relevant published location and period evidence, and **40,712** positive site×guild×evidence×time-resolution reporting units. Twenty-seven true original GIFT plant Direct flower morphological profiles (4 single, 11 same-accepted-species pairs, 9 triples and 3 quadruples). Five guild pairs ×2 site-time scopes (same year vs same year/month) ×2 evidence scopes (positive grade1–3 vs pollen carriage/confirmed grade1–2) ×27 flower profiles = **540 nested source panels**, NOT 540 independent biological tests.
- [Pure and reproducible analysis](../../src/fd/dopi_positive_coaccess.py); [synthetic ecological non-absence, source variation and finite-genus expectation tests](../../tests/test_dopi_positive_coaccess.py); [original source runner](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_original_dopi_gift_guild_motifs.py) and [SHA-pinned original workflow](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/.github/workflows/run-fd-original-dopi-guild-motifs.yml).

### True co-documentation estimand, not visitation preference or fitness

For each independently source-verified paper and its **identical published coordinate string (or British grid cell)**, observation year and optionally month, require positive evidence for guild A and guild B. Restrict to accepted species with the exact necessary Direct floral axes. Let `U` be the union of positively documented visited plant species for the two guilds at that source site/period, and `C` their SAME-species intersection positively documented by BOTH guilds. The raw source profile frequency contrast is:

`P(motif | C documented BOTH) - P(motif | U documented EITHER)`.

This denominator is **positive-recorded plant species**, not the full site's flower availability, and a plant with no documented link to guild B is not a tested ecological non-visit or absence. Require >=3 known positive source species in `U`, >=1 accepted species positively documented in `C`, >=10 independent publication keys and >=10 original distinct publications with at least one motif-positive plant in `U`. Average sites within each original publication, average publications equally, and bootstrap the entire original publication 999 times. CIs are exploratory and not family-wise adjusted.

**Genus-conditional finite-margins null:** Within each original study/site/period and each botanical genus, hold the number `n` of union positive plant species, `k` of directly motif-positive species and `c` of source co-documented guild species fixed. If true floral phenotype labels are exchangeable within that group, `E[O] = c k / n` for the number `O` of motif-positive co-documented species. Compute observed-minus-expectation, divide by the number of co-documented species in *informative* genera, and independently bootstrap source papers. A genus is informative only with `n>=2`, `0<k<n` and `0<c<n`; no genus variation is NOT a zero biological effect. This null conditions on botanical lineage and site source counts, not on true visit probability, plant abundance or source effort.

## ORIGINAL positive site/year/month evidence (study-equal)

Original literature-positive DoPI source, same published location and YEAR AND MONTH, visitation/pollen grade1–3. **All differences are percentage-points of genuine accepted plant category frequencies among positive recorded host species, not pollination fitness.** Original 95% intervals are study-cluster bootstrap, no multiple testing adjustment.

| Positive visitor pair A + B | Original real morphological motif | Independent source papers with qualifying positive plant co-documentation | Original site×month/year with same-accepted-plant both guilds recorded | Different accepted positive co-documented plant species | Observed co-documentation frequency contrast [95% interval], pp | Informative genus-null source papers |
|---|---|---:|---:|---:|---|---:|
| **Bombus + Syrphidae** | **open radial form** | **63** | **793** | **321** | **−6.841 [−11.347,−1.963]** | **1** — NOT_IDENTIFIABLE |
| Bombus + Syrphidae | open radial × shallow depth | **49** | **355** | **77** | **−2.538 [−4.644,−0.495]** | **0** — NOT_IDENTIFIABLE |
| Bombus + Syrphidae | zygomorphic symmetry | **63** | **850** | **351** | +4.978 [−0.486,+10.117] (uncertain) | **0** — NOT_IDENTIFIABLE |
| **Butterflies + Bombus** | zygomorphic × deep tube | **36** | **239** | **61** | **−2.385 [−4.523,−0.491]** | **0** — NOT_IDENTIFIABLE |
| Other non-Bombus bees + Bombus | open radial form | 94 | 1,582 | 461 | −4.131 [−7.638,−0.462] | 1 — NOT_IDENTIFIABLE |
| Other bees + Bombus | zygomorphic symmetry | 101 | 1,728 | 528 | +4.103 [+0.186,+8.065] | 1 — NOT_IDENTIFIABLE |

The positive open/radial co-documentation comparison is **negative**, thus a naive hypothesis that open/short-tube flowers necessarily have *more* dual-guild visitation records is NOT uniformly borne out in the observed source literature. However, all source papers differ in observed plant taxa/record frequency, and a high guild-overlap degree alone does NOT measure effective multispecies pollen service. These source results do **not** falsify broad access as a physical mechanism; they expose a missing population, plant abundance and independent source study control.

**Original source SUPPORT gate across all 135 guild×profile endpoints per site-time/evidence tier:**

| Source site/time × DoPI evidence grade | >=10 papers with eligible observed BOTH-guild plant and genuine motif support | Of which with >=10 genuinely informative genus-fixed source papers |
|---|---:|---:|
| Same published source site + year, positive visits grades 1–3 | **42/135** | **8/135** |
| **Same site + year + month**, positive grades 1–3 | **33/135** | **3/135** |
| Same site + year, pollen carriage/confirmed grades 1–2 | **4/135** | **0/135** |
| Same site + year + month, pollen grades 1–2 | **0/135** | **0/135** |

On the most stringent same site+year+month genus-fixed source subset, only **3/135** morphology×guild panels have enough source support for a numerical residual. The single panel whose unadjusted **95% study-cluster interval** excludes zero is *other bee + Bombus / actinomorphic×deep*, **−22.467 pp [−39.394,−4.281]** on **11** source publication keys and 31 site-periods where morphology varied within genuine genera. This is a highly selective, posthoc result among 540 nested panels. It is **not** family-wise corrected, not an island pollinator selection coefficient, and is not a confirmation of the old tropical bee-to-butterfly hypothesis. It is a candidate for separate replication, not a main claim.

Across the same site/year/month positive-visit source tier **7/135** raw motif intervals exclude zero without correction, but only 33/135 are support-qualified. A few apparently strong raw source contrasts are primarily fixed genus/motif effects with no within-genus test. The genus-fixed evidence at pollen grades 1–2 is NONE in both site-year and year/month tiers. Thus **functional pollination shifts are not independently identified.**

## Interpretation and main-paper placement

Three results must be kept separate:
1. **Global original GIFT botanical flora**: relative to *potential* mainland source species, floral functional repertoires are geographically heterogeneous; genus membership/lineage source composition makes large signed contributions. The same accepted species Direct traits do not establish island versus mainland within-species evolution.
2. **Independent European DoPI positive visitors**: real floral architecture predicts some source visitor-guild frequencies, and same accepted plant species can be documented by multiple guilds at the same published site and period; these are **visitor associations**, not an effective pollinator syndrome or unique guild membership. The observed co-documentation overlaps do not show universal enrichment for open flowers.
3. **The unmeasured mechanism**: actual local island bee/bird/butterfly/hawkmoth availability, pollen deposition, reward and phenology, original study effort, plant population floral morphology and genuine historical arrival/establishment remain unsupported. The old island-evolution/guild-shift inference would overstate evidence.

**Manuscript decision:** Keep 4,505 GIFT island functional beta / taxonomic assembly and original source-relative floral structure as main analyses. Present DoPI visitor association and this co-documentation test as qualified external syndrome-compatibility and positive-pollinator-portfolio evidence, not a guild shift. The evidence is consistent with *possible multiple-guild access* but does not demonstrate the realized number of effective pollinators.

A peer-reviewed oceanic island network example [Xia et al. 2020 (Yongxing Island)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7434763/) found more plant visitor species than traditional syndrome expectations; broader [Dellinger 2020](https://doi.org/10.1111/nph.16793) reviews the importance of secondary visitor guilds and why **effective pollen transfer** is a different estimand from simple positive visitor records. These motivate but do not validate extrapolating a European source network to the user's four global island biogeographic regions.
