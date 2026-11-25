# 14 November 2026 — Oceanic versus continental island floral assembly, source overlap and endemic evidence

> **Source evidence update, 18 Nov 2026.** The later [original GIFT 4,505-island genuine two-trait source-referenced comparison with independent geological origin and Weigelt LGM history](../20261118_original_geology_lgm_genus_flower_pair/README.md), protected [run #44245811804](https://github.com/zuizui0223/island/actions/runs/44245811804) SUCCESS, finds tropical oceanic−continental **radial×deep +6.69 pp** (95% block [3.56,10.00], 37/124 source-qualified original floras) and LGM-unconnected−connected **+6.90 pp** ([2.91,10.71], 257/347) in two-trait SAME accepted-plant floristic SOURCE-RELATIVE frequency. These use a broader independently trait-specific >=5-species sample than the older 928 triple-complete >=30 botanical specimen cohort and do not overturn this source's strict conclusion of ONLY THREE fully matched oceanic−continental independent tropical geological pairs. The new outcomes are observation-, original source-, and lineage-sensitive botanical composition, not two independent causal historical replications, confirmed native-only effects or effective pollinator turnover. Geology covers only 425/4,505 original botanical island floras and the two historical axes must not be conflated.

**Status: original-source real-data analysis complete, two independently seeded numerical executions SUCCESS.** Ecological and geological source effects remain observational and conditional. No confirmed oceanic/continental colonization success, in-situ adaptive evolution or flower-pollinator guild replacement.

- Original protected run [island Actions #43574635085](https://github.com/zuizui0223/island/actions/runs/43574635085) succeeded on **both** original numerical seeds (20261114 and 20261115), uploaded only anonymous regional aggregates (IDs #14106150932 and #14106175846). These are numeric/bootstrap repeatability checks on the SAME original plant/island observations, not independent biological studies.
- Reusable methods: [geological_establishment_adaptation.py](../../src/fd/geological_establishment_adaptation.py), [geological_disjoint_matched.py](../../src/fd/geological_disjoint_matched.py); source [island runner](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_geological_establishment_adaptation.py) and [source workflow](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/.github/workflows/run-fd-oceanic-continental-assembly.yml).
- Original reference comprises **2,234** high/medium source-Direct unambiguous species each with floral primary colour × floral form × symmetry, **928 original geographic islands** with at least 30 such recorded plants. Independent audited Roeble 2024/Zizka 2022 *geological* origin provides **176 continental / 258 oceanic** labels among the underlying full GSHHG universe; 8 mixed and remaining unverified are NOT classified as either geology type. **LGM mainland connection is a distinct geological-history exposure** and is not used to set origin class.
- Exactly **148 geological-label + original ≥30 Direct floral + five-nearest GIFT native mainland pool source-comparable** islands in the present study; northern midlatitude **50 continental/23 oceanic**, tropical **21 continental/43 oceanic**, with other geographic regions too small. A GIFT donor pool is a candidate botanical source checklist rather than observed propagule arrival, true ancestor or established local absence.

## Source-repertoire association and denominator falsification

A source-relative categorical floral profile response uses exact finite-population 30-species expected marginal Rao Q, joint actual-profile Rao Q and colour × intact(form,symmetry) pairing departure (NOT 17-dimensional Gaussian KDE). Across the 148 origin-qualified islands, the present source-normalized *marginal* and *joint* Rao origin contrasts for northern midlatitude and tropical broad floras have conditional 95% confidence intervals spanning zero under both geographic-block and candidate-source-signature weighting; these 30-species model results cannot show a universal oceanic floral functional-space contraction.

In northern midlatitude all-recorded flora the **colour–architecture pairing source-relative** oceanic-minus-continental block-equal point estimate is **+0.01241 [0.00205,0.02230]**, but nearest-GIFT-source-signature-equal CI **+0.00907 [−0.00166,+0.01973]** crosses zero. In tropical all-recorded floras, block **−0.00361 [−0.01342,+0.00585]**, source **−0.00297 [−0.01196,+0.00598]**. The posthoc geographic block pattern is NOT an independent geological discovery after source-level uncertainty and multiple overlapping endpoints. In WCVP regional-compatible smaller cohorts, strong origin contrasts are not supported by both geological origin class sample/block/source counts.

The initially striking asymmetric botanical genus overlap deficit on oceanic floras is **denominator sensitive**:

| Scope (all-recorded) | Oceanic − continental fraction of SOURCE genera shared with island [block 95%] | Oceanic − continental expected shared fraction under fixed 10 distinct island and 10 source genera [block 95%] |
| --- | --- | --- |
| Northern midlatitude | **−0.08552 [−0.14256,−0.02752]** | **−0.00153 [−0.01071,+0.00893]** |
| Tropical | **−0.07631 [−0.13115,−0.02535]** | **+0.00037 [−0.00643,+0.00781]** |

When re-denominated as the fraction of ISLAND genera also observed in source, or symmetrized with genus Jaccard/Sørensen, the oceanic-minus-continental 95% intervals include zero. Therefore this result is evidence of **recorded source and floristic richness/composition imbalance**, NOT an independently established establishment-probability or taxon-filter selection effect.

## Stringent geographic and archipelago independence — origin causality unidentifiable

Using morphology-outcome-blind calipers on original area, distance, Direct species evidence, GIFT donor genus opportunity and 4 climate PCs, then forbidding any reused physical geographic 10° block and any reused independently labelled Weigelt archipelago across different matched pairs:

| Original region and plant-origin tier | Oceanic–continental independent, covariate-compatible matched pairs |
| --- | ---: |
| Northern midlatitude / all recorded | **3** |
| Tropical / all recorded | **3** |
| Northern midlatitude / WCVP regionally native-compatible | 2 |
| Tropical / WCVP regionally native-compatible | 0 |
| Northern midlatitude / strict GIFT native | 1 |
| Tropical / strict GIFT native | 0 |

These counts were confirmed twice, exactly matching independent-seed outputs. With only 3 disjoint comparison pairs per broad region, no p-value or origin-specific evolutionary inference is presented. The previous broader 37-pair evaluation had pseudo-independent spatial components, and its strongest corrected alternative yielded only 10 independently matched pairs and no BH/BY FDR-supported geological repertoire effect.

## Critical correction: original ENDOMICITY is NOT native-origin status

Earlier preliminary geological code looked for an invalid `origin_status=="endemic"`, yielding an apparent **0 endemic** count. That result was **invalidated**: the independently archived source `chapter1_status_flora.csv.gz` has separate `origin_status` (native, introduced, unresolved), `endemic_status` (endemic, nonendemic, unresolved), and `floristic_status` (endemic, native_nonendemic, native_endemism_unresolved...). The new original runner now merges original status by EXACT island×accepted species after the standard 928 Direct flora selection, and fails closed on missing or contradicting source labels. Strict endemic requires three concordant original fields and is NEVER imputed from WCVP TDWG-L3 regionally native range.

Corrected original 148-island cohort, **island×accepted plant species rows**, not necessarily unique biological species nor distinct ancestral radiations:

| Broad regime | Geological origin | Strict original endemic plant×island Direct triple records | Endemic plant×island records with a distinct *same-genus* source-Direct botanical candidate |
| --- | --- | ---: | ---: |
| Northern midlatitude | Continental | **18** | **12** |
| Northern midlatitude | Oceanic | **4** | **0** |
| Tropical | Continental | **55** | **7** |
| Tropical | Oceanic | **4** | **0** |

Across ALL regions of the 148-island intersection, **21 original endemic plant×island source-congener comparable records** were found, but island-level paired endemic versus native-nonendemic source-congener comparisons have only **one qualifying island** in each major continental-region group and **zero** in both main oceanic groups. There is no sample support for a geological class-specific endemic floral divergence coefficient. **Zero comparable original source congeners is NOT zero real related mainland plants, nor proof of absent adaptation or dispersal.** The recorded source taxon opportunity and endemic primary traits are very sparse on oceanic islands.

## Inference boundary and concrete next acquisition

**Establishment / taxonomic sorting:** directly accessible only as candidate-source-to-island **floristic composition** conditional on source records, native status, genus, detection and original geography. Fixed 10-genus and geographical-independence tests fail to establish an oceanic-specific filtering beyond source and sampling opportunity. Biological arrival, colonization, population establishment and extinction have not been observed.

**Adaptation / in-situ evolution:** correct endemic labels now enable support counting, but **the oceanic island endemic same-genus source contrasts are unidentifiable in both the tropical and northern-midlatitude original 3-trait intersection**. A stronger test must first acquire more independent island-native/endemism-qualified floral traits on oceanic islands, record actual same-island populations with comparable mainland conspecific or phylogenetically justified sister taxa, and match pollinator function or pollen limitation/seed-set evidence. Original single state per species cannot identify population-level divergence.

High-leverage next gate: do **NOT** require triple colour+form+symmetry completeness before auditing original endemicity. First enumerate, by independently verified geological island and origin tier, available *single-trait* original Direct colour, form, symmetry, SC/SI, selfing and floral access records, and the original botanical source counterpart availability, before preregistering any cross-trait adaptive evolution coefficient. Add exact genus/family and published true-tip source limits, source/document provider independence, region and sampling effort guards. No original named taxa, source family/archipelago IDs, island IDs or exact positions may be released from protected Actions.

Scientific source and manuscript remain in **DRAFT** pending source, taxon, pollinator-effect and evolutionary inference tests.
