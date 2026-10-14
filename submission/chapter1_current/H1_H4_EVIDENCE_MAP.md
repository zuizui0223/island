# Chapter 1 — H1–H4 publication evidence map

> **Editorial navigation, not a new scientific lock or analysis.** Updated 14 October 2026.
> The active paper is selected only by [chapter1_submission_current.json](../../config/chapter1_submission_current.json). Its current H1 specification is the **4 October 2026 seven-trait, finite-spatial-cluster-\(t\), unadjusted-pointwise** analysis, not a pooled/domain syndrome score. Use [the current manuscript](MANUSCRIPT.md) for submitted prose and [the final H1 report](../../results/h1_final_traitwise_t_20261004/README.md) plus [the corrected H2–H4 evidence](../../results/geography_20260924/) for numbers.
> This map must not be used to replace frozen results or automatically promote exploratory findings.

## One publishable question

**Does island isolation recurrently favour reproductive assurance while allowing regional differences in floral display and architecture, and is that geographic pattern aligned with independently measured pollen limitation?**

The paper's distinctive contribution is the **joint evidence structure**—continuous geographic isolation across assembled island floras, separately modelled reproductive/floral responses, independent published pollen-addition evidence and exact-species functional compatibility—not a claim to have discovered Baker's law, the effect of selfing on pollen limitation, or a newly identified universal pollinator mechanism. See [the novelty statement](NOVELTY_STATEMENT.md) and [prior-art audit](../../docs/chapter1_prior_art_novelty_audit_20261008.md).

## Four hypotheses, four separate evidentiary roles

| Paper step | Estimand and comparison | Current defensible result | Status and key restriction | Primary figure |
|---|---|---|---|---|
| **H1 — Traitwise island responses** | Seven independent binary-trait slopes against corrected isolation, separately in four geographic strata; beta-binomial logit with area, climate and spatial-block inference | **Self-compatibility (SC) increases in all four regions** in broad All and WCVP regional-native-compatible analyses. The other six trait directions/support vary | The 17/28 broad All nominal associations are **two-sided unadjusted, pointwise**; no global pooled H1 score, no universal seven-trait effect or formal between-region slope test | 4 |
| **H2 — Conditional floral response** | \`generalized_accessible\` and colour responses against isolation **after \`selfing_core\` adjustment**; area/climate included | The accessibility-distance effect remains FDR-supported in northern high latitude and tropical **All**, but not northern midlatitude or southern extratropical; tropical Direct-only loses FDR support | **Conditional statistical decomposition**, not mediation, direct selection, or universal causal pathway independence | 5 |
| **H3 — Independent ecological pressure** | Experimental pollen-limitation log-response-ratio against corrected continental-isolation distance in published GloPL data | Corrected primary \`beta = +0.09191\`, \`SE = 0.03806\`, two-sided finite-publication \`p = 0.01594\`; offshore-only gradient also positive | Publication-weighted, adjusted **association**; supplemental-only measurement sensitivity unsupported; neither pollinator abundance loss nor historical cause identified | 6A |
| **H4 — Functional compatibility** | GloPL pollen limitation against the **frozen Direct-only H2 trait scores** on exact-matched species, adjusting for distance, context and study conditions | \`selfing_core\`: \`beta = -0.29830\`, \`p = 0.00417\`; \`generalized_accessible\`: \`beta = -0.29566\`, \`p = 0.02334\` | **Explicitly post-hoc functional triangulation**. Present-day trait/limitation covariance does not prove historical selection or cross-dataset mediation | 6B–C |

The causal bridge **from earlier pollen limitation to the current trait composition** remains unobserved. H1/H2 and H3 are associations from different data domains; H4 tests current functional compatibility on exact-species overlap, not historical time-ordering.

### Denominators and measurement boundaries

- Corrected geographic universe: **8,264** island units. The **broad H1 analysis union** includes **4,379** observed islands; the trait database includes **106,295** accepted angiosperm species.
- Plant trait cells resolved: **222,688 / 318,885 (69.83%)** across three raw evidence axes. The current H1 outcomes are **seven atomics**, not three aggregate response scores.
- Floristic-origin uncertainty is large: among **1,039,757** island-by-species records, **869,784 (83.65%)** have unresolved origin status. The primary All results concern **contemporary GBIF-observed floras**, not demonstrated native colonization or within-lineage evolution. WCVP TDWG-L3 regional-native compatibility does not establish exact focal-island nativeness.
- Independent GloPL: **2,969 experiments**, **1,408 measurement cells**, **1,248 sites**, **919 publications**; **996** sites are true continental zero-distance observations. The corrected coastline distance repaired **1,113** spurious H1 zero distances and excluded one continental fragment. This was a **post-hoc measurement repair**, not prospective replication.
- H4 exact-score support differs: \`selfing_core\` matches **455 species / 409 publications**, whereas \`generalized_accessible\` matches **143 species / 143 publications**. These are not matched island-plant manipulations.

## H1: what actually differs among the four regions

A plus/minus below means a **direction with nominal two-sided p < 0.05 in the primary broad All H1 model**, not a multiple-testing-adjusted discovery; omitted atomics lack such support. Source: [full final traitwise results](../../results/h1_final_traitwise_t_20261004/traitwise_results.csv) and [its readable summary](../../results/h1_final_traitwise_t_20261004/README.md).

| Region | Reproduction: primary broad All | Primary broad All floral structure / colour | Sensitivity and biological reading |
|---|---|---|---|
| **Northern mid-latitude** | SC **+**, selfing mating system **+** | Generalized form **+** | WCVP retains SC but does not support every broad floristic-form coefficient. There is no defensible generic “white-flower evolution” inference |
| **Northern high latitude** | SC **+**, selfing mating system **+**, autonomous selfing **+** | Actinomorphy **+**, shallow/open tube **+** | WCVP retains SC, autonomous selfing, actinomorphy and shallow/open tube; this is the clearest access-generalization combination |
| **Tropical** | SC **+**, selfing mating system **+**, autonomous selfing **+** | Generalized form **+**, actinomorphy **+** | WCVP retains SC and actinomorphy and supports plain colour, but other patterns change with floristic origin and trait-evidence scope |
| **Southern extratropical** | SC **+** | Plain colour **+**, generalized form **+**, shallow/open tube **−** | WCVP retains SC and plain colour; a single direction of floral simplification **does not** describe deeper tubes and generalized form together |

**Only the SC association is supported by the pointwise trait tests in both primary broad All and WCVP All in all four regions.** The other six atomic traits are not a jointly proven universal phenotype, and \`p < 0.05\` is not a family-wise significance guarantee.

Secondary reported raw flower-colour states (white, red/pink, yellow/orange, blue/purple, green/brown/inconspicuous) and **architecture conditional on focal colour** can describe display–form coupling. They must be presented as secondary, geographically context-dependent associations rather than identifying realized bird, butterfly, or large-bee pollinators. Use the **corrected** tables in [results/geography_20260924](../../results/geography_20260924/), not pre-geography-correction numeric summaries.

## H2: what the selfing-adjustment test can and cannot establish

\`selfing_core\` is constructed from self-compatibility/incompatibility, mating system and autonomous selfing; **it excludes colour and floral form**. The current conditional accessibility result is:

| Region | Broad All, finite-cluster FDR q | Interpretation |
|---|---:|---|
| Northern mid-latitude | 0.1764 | Not FDR-supported after measured reproductive assurance |
| Northern high latitude | **0.00243** | Retained conditional floral accessibility response |
| Tropical | **0.01913** | Retained in All; **Direct-only q = 0.1267**, not supported after correction |
| Southern extratropical | 0.2977 | Not FDR-supported after measured reproductive assurance |

This challenges the compulsory **\`isolation → selfing → all floral change\`** statistical account. It does **not** establish a direct pollinator-mediated selection route, because source pools, species sorting, unmeasured traits and measurement error remain alternatives. The within-species evolutionary origin of the differences is unresolved.

## H3 and H4: independent pressure, then a functional check

**H3.** Pollen addition measures current reproductive gain from added pollen, not the number of pollinators. With each publication contributing total weight one, the corrected primary isolation coefficient is **+0.09191**, finite-publication two-sided **p = 0.01594**. The post-hoc offshore-only slope is **+0.22031**, **p = 0.02459**; a binary mainland/offshore indicator alone was unsupported, while supplemental-only primary measurement sensitivity was also unsupported (**p = 0.31082**). Do not describe this as a controlled island treatment.

**H4.** The exact-species H2-to-GloPL check has two distinct sample sizes and corresponding negative slopes:
- Reproductive assurance: **455 species / 409 publications**, **−0.29830** (finite-publication **p = 0.00417**).
- Generalized accessibility: **143 species / 143 publications**, **−0.29566** (finite-publication **p = 0.02334**).

The retained atomic check is strongest for **autonomous selfing**; SC alone is weaker and imprecise. Historical Route A/B tests asking whether these states moderate the *isolation slope* of pollen limitation were unsupported. H4 was formulated after the parent analyses, so this is **post-hoc compatibility** rather than confirmatory causal evidence.

## Main paper structure and live locations

| Reader sequence | File(s) | Purpose |
|---|---|---|
| Submission gate / selected version | [current selector](../../config/chapter1_submission_current.json) | Final source of truth for active results, geography and H1 inference |
| Abstract, Introduction, Methods, Results, Discussion | [MANUSCRIPT.md](MANUSCRIPT.md) | Already structured as H1 → H2 → H3 → H4 |
| Figure 1 / 2 / 3 | [FIGURE_CAPTIONS.md](FIGURE_CAPTIONS.md) | Data/geography, pipeline, constraint–response triangle with **dashed unproven causal edge** |
| Figure 4 | [Final H1 traitwise](../../results/h1_final_traitwise_t_20261004/) | Four-region seven-atomic results and evidence/origin sensitivities |
| Figure 5 | [Corrected H2 tables](../../results/geography_20260924/) | Selfing-conditional floral form, colour and accessibility |
| Figure 6 | [Corrected H3/H4 tables](../../results/geography_20260924/) | GloPL isolation gradient and post-hoc functional bridge |
| Manuscript + supplemental delivery | [package manifest](PACKAGE_MANIFEST.md), [submission checklist](SUBMISSION_CHECKLIST.md) | Rights, metadata, DOI and journal-ready conversion still gated |
| Old H1–H5/v8–v14 scientific contracts | [historical index](../../docs/CHAPTER1_HISTORY.md) | Replay/provenance only; **not live H1–H4** |

### Historical or exploratory results that must not be promoted by tidy-up

- **Pre-correction** v13/v14 H1 composites, joint omnibus/directional tests and the historical H1–H5 ladder are *not* the final H1 inference.
- **5 October common-support and nonlinear diagnostics** in [common-support](../../results/h1_regional_common_support_20261005/), [region × isolation](../../results/pollination_syndrome_region_scale_interaction_20261005/) and [far-tail genus attribution](../../results/pollination_syndrome_far_tail_genus_20261005/) are **post-hoc** and outside the active submission. Regional distance distributions differ substantially; the tropical yellow/orange deep-tube signal is concentrated in very remote islands, is not a simple shared-range tropical-versus-northern contrast, and its native far-tail persistence could also reflect genus/dispersal filtering. Do not recast it as butterfly replacement or newly confirmed mechanism.
- Genus/species source assembly, pollination-channel labels, GloBI interactions and inferred bird/butterfly/large-bee concordance may motivate Discussion but cannot identify the historical source of the H1 traits. Do not treat absence of GBIF records as true pollinator absence.
- **H4 is not causal mediation**; the observed H2 traits, historical pollen limitation and subsequent selection have not been traced in the same population through time.

## Publication gate: what is done vs what blocks submission

**Ready as an evidence spine:** final seven-trait H1, corrected H2–H4 tables, claim-limited manuscript, current figure captions, and frozen provenance remain intact. This map only improves navigation and ensures the right results are cited; it changes no numerical model or selector.

**Blocking external/compliance tasks:** Ecology Letters data-rights policy decision, final author/affiliation and declaration fields, permanent software/results DOI, and compilation of journal-ready submission files. Track these in [SUBMISSION_CHECKLIST.md](SUBMISSION_CHECKLIST.md) and [SUBMISSION_READINESS.md](SUBMISSION_READINESS.md); do not tick them as completed without evidence.

**Recommended manuscript claim:** *Geographic isolation is associated with recurrent reproductive assurance, especially self-compatibility, but non-uniform floral change. Independent pollen-supplementation evidence reveals greater contemporary pollen limitation with isolation, while exact-species trait comparisons provide post-hoc functional consistency with lower limitation. The historical causal link remains unresolved.*
