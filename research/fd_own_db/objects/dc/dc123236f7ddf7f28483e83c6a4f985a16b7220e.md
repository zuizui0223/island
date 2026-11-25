# 13 November 2026 — Genus turnover versus within-shared-genus contribution to island rare flower-profile space

**Original-data macroecology, completed and independently code-reproduced.** Full 928-island Direct flora universe; **736 northern midlatitude**, **106 tropical**, **53 northern high latitude**, **24 southern extratropical** islands have ≥30 independent candidate GIFT mainland-source species after **excluding EVERY focal-island accepted species from the donor source lists**. Original protected source execution [island Actions #43385444722](https://github.com/zuizui0223/island/actions/runs/43385444722) **SUCCESS**, improved independent compiled-lookup replay [#43386044096](https://github.com/zuizui0223/island/actions/runs/43386044096) **SUCCESS**; anonymous [artifact #14019250011](https://github.com/zuizui0223/island/actions/runs/43386044096/artifacts/14019250011). The original aggregate JSON bytes have **exactly the same SHA256: 97e23194924a0e70011ff37bf3f7b6a1453fda9bee3a347985bbfb8a6d785057** in both independent workflows, despite radically different code execution speed. Original data rights and source IDs are protected; the public artifact contains region-level aggregates only. Scientific [fd PR #16](https://github.com/zuizui0223/fd/pull/16) and original-input [island PR #326](https://github.com/zuizui0223/island/pull/326) remain Draft.

## Principal macro question

When an island flora is relatively deficient in globally **rare recorded flower colour × floral form × floral symmetry phenotype categories** relative to candidate nearby mainland native flora, is that because of
1. **which genera are recorded at all** on each side (island-only versus donor-only observed genera, here called 'exclusive/turnover' genera), or
2. the differences between the *different species* recorded in **genera present on both sides**?

The distinction is about **documented species composition**, not loss of actual functional pollinators, within-species floral evolution, measured dispersal, historical colonization or extinction.

## Exact categorical functional breadth: original evidence and attribution

Original globally complete high/medium, unambiguous Direct trait evidence: **2,234 unique original accepted plant species**, 100 ACTUALLY observed three-trait categorical profiles. Source-defined rare reference = **29 profiles carried by just one or two Direct species globally**; common reference = **25 profiles represented by ≥26 species** (46 intermediate). The metric coverage neighborhood is Hamming distance **at most one of the three categorical traits**, evaluated only on the original actual global reference profiles.

For each island and its geographically assigned five-nearest candidate GIFT native mainland floras, sample **exactly 30 different original source species on each side, 40 replicates**, after globally removing all exact focal island species names from that island's source pool. For both sides, partition sampled species by **whether their entire ORIGINAL unthinned island or source flora contains that genus**. Categories may have multiple overlapping neighboring sampled species. We calculate *two-player Shapley contributions* to the observed union-of-neighborhood-coverage \`V(Shared ∪ Exclusive)\`:
\[
\phi_{\mathrm{shared}}
 = \tfrac12\left[V(\mathrm{shared})+V(\mathrm{all})-V(\mathrm{exclusive})\right],
\quad
\phi_{\mathrm{exclusive}}
 = \tfrac12\left[V(\mathrm{exclusive})+V(\mathrm{all})-V(\mathrm{shared})\right].
\]
These sum **exactly** to the original observed total coverage, splitting overlap equally rather than double-counting the same reference profile. This is a mathematical *attribution accounting*, not a causal Oaxaca regression or phylogenetic history.

Whole original 10° spatial blocks and the independently stored first-nearest GIFT donor signatures are two overlapping conditional equal-weight resampling levels (1,999 draws). Source identities are **potential geographically near source floras**, not directly confirmed phylogenetic ancestors.

## Original full island vs mainland *RARE-profile support* decomposition

Numbers below are **island-minus-mainland candidate source coverage, in global rare-reference FRACTION units** (multiply by 100 for percentage points), **block-equal mean and conditional 95% block interval**:

| Region, all recorded | Original comparable islands / blocks / first source signatures | Total rare support gap | **Shared-genera** Shapley gap | **Only one-sided observed genera** Shapley gap |
| --- | --- | --- | --- | --- |
| **Tropical** | **106 / 42 / 36** | **−0.02970 [−0.04803,−0.00990]** | +0.00945 [−0.00907,+0.02756] | **−0.03914 [−0.06488,−0.01214]** |
| **Northern midlatitude** | **736 / 51 / 64** | **−0.02715 [−0.04837,−0.00614]** | **+0.09867 [+0.07500,+0.12387]** | **−0.12582 [−0.14713,−0.10391]** |
| Northern high latitude | 53 / 7 / **3** | −0.09015 [−0.11756,−0.05548] | +0.09408 [+0.06268,+0.11928] | −0.18422 [−0.21908,−0.14615] |
| Southern extratropical | 24 / 11 / 9 | +0.01322 [−0.03345,+0.05790] | −0.01143 [−0.02636,+0.00443] | +0.02465 [−0.02999,+0.07348] |

The union decomposition sums exactly: in the tropical row, (+0.00945)+(-0.03914)=(-0.02970), rounded. Under equal first-nearest GIFT donor signature, tropical exclusive-genera rare gap **−0.03844 [−0.07095,−0.00683]**, versus shared-genera **+0.01398 [−0.00490,+0.03297]**. Northmid exclusive-genera rare gap **−0.12373 [−0.14142,−0.10468]** source-equal, shared-genera **+0.09741 [+0.07666,+0.11838]**. Both physically block-equal and candidate source-equal sensitivity show an excess deficit of rare-profile reference coverage attributable to **genera recorded on only one side**.

This decomposition is **not** a treatment experiment: a GIFT source may list many more observed/curated genera than an island, native/provenance resolution and source archive observation effort can create apparent one-sided groups. Source total sampled species counts were matched but the number of GIFT candidate source genera was **not fixed**. In this original flora cohort, the tropical island median is **45 distinct recorded genera**, candidate donor source **74**, with only **7 shared genera**; northmid medians are 54/151/23. These are sample-level descriptive medians, not ancestral lineage counters. Therefore the one-sided genus contribution can be structurally large even without real selection or missing pollinator service.

## The more stringent *same-five-distinct-genera* cross-island comparison

The earlier full 15-species / matching original genus composition test had 14/106 source-matchable tropical islands. To use a much wider *but distinct* original source-supported genus cohort, predeclare a different **n=5 distinct shared genera**, draw **one different accepted source plant species per genus** on each island/donor side (40 independent repeats) and estimate categorical **rare-minus-common one-trait coverage gap**; **identical GENUS composition at each draw by construction**. A second *n=8 distinct shared genera* independently checks target-scale sensitivity. No source species may be the exact accepted focal island species, and genus matching is predeclared from taxon names, not observed flowers.

Minimum report support: ≥10 distinct islands, ≥5 original physical spatial blocks and ≥5 distinct first-nearest source signatures. Results:

| Region and scope, same-five-genera controls | Islands / blocks / first source signatures | Rare-minus-common source-relative categorical support gap | 95% conditional block interval |
| --- | ---: | ---: | --- |
| **Tropical all-recorded** | **81 / 33 / 30** | **+0.00261** | **[−0.01579,+0.01982]** |
| **Northern midlatitude all-recorded** | **733 / 50 / 63** | **−0.00867** | [−0.02004,+0.00379] |
| Southern extratropical all-recorded | 17 / 11 / 7 | +0.03515 | [+0.00395,+0.06527] |
| Northern midlat strict native | 67 / 16 / 23 | −0.00065 | [−0.01741,+0.01556] |
| Tropical strict native | 14 / 8 / 7 | −0.01146 | [−0.05800,+0.05204] |

An equal-source-signature conditional interval for tropical n5 also crosses zero (numeric estimate in original aggregate). Source-excluding genus-matched **n=8** controls remain unsupported: **tropical 50/21/19**: **+0.00914 [−0.01475,+0.03532]**; northern midlat **714/46/57**: **−0.01169 [−0.02702,+0.00349]**. The northern high latitude n5/n8 cohorts have only **three** first-nearest GIFT source signatures, so their effects are suppressed under the original source-independence-support gate.

A positive separate southern 17-island genus-n5 comparison is an isolated small, posthoc conditional cohort, not global independent repeated history or a contradiction of the full southern original n30 source-relative rare coverage interval that includes zero. The new n5/n8 estimands *by definition* occupy much less original category reference than n30, and selected source-shared genera are biologically not a random sample of all island flora.

**Key synthesis:** Within the **same original genus identities**, there is no supported rare-profile neighborhood deficit in tropical and northmid source floras under the n5/n8 original reference scheme. By exact allocation, the all-recorded n30 rare reference shortfall is concentrated in *profile neighborhoods covered by different genus repertoires* between islands and candidate sources. **This is compatible with lineage-compositional filtering but does NOT prove ancestral genus exclusion caused the deficit, does NOT demonstrate within-species evolution, and does NOT show a loss of pollination services.**

## Caveats and next research decision

1. **GIFT source genus diversity is not held fixed** across island and mainland, so a sizeable "exclusive genus" effect can partly be a **source-repertoire richness mismatch** even though each biological phenotype draw has the same n=30 accepted plant species. Correct this with a conditional genus-richness/source-area/lineage-availability baseline before interpreting any environmental filtering.
2. The global 100-recorded-profile direct trait reference may be uneven in taxonomic/source evidence quality. Rare bins represent documented **1–2 Direct species**, not true global floral functional rarity. A globally species-prevalence-weighted original support measure previously failed to show robust source-relative tropical total contraction. Common-reference neighborhood coverage is near mathematical 100% in tropical flora, so rare-minus-common comparisons have a built-in ceiling sensitivity.
3. **Genus is not a phylogenetic covariance matrix**. Same-genus source plant names are not verified island ancestral lineages; all accepted species have one immutable direct flower phenotype label globally. Reanalysis cannot establish local species trait adaptation or trait plasticity.
4. 87 island floras only have enough strictly source-**native** Direct traits, including 17 tropical; the tropical strict-native Shapley rare support interval includes zero. This is a major provenance/sampling-selection limitation.
5. Original qualified five-pollinator guild states lack geographic positivity, original GBIF data lack 2-window comparable histories, and external same-species documented visitor comparisons are sparse. Do **not** attach pollinator decline/shift or reproductive fitness causality to this morphological composition pattern.

The next test that could identify an additional ecological mechanism is to acquire independent taxonomically matched, native-verified donor/recipient island floras with adequate source-region richness, then condition genus entry/survival against geographic area, isolation and documented donor availability. Data currently support a **descriptive macroecological assembly result**, not an evolutionary island-syndrome mechanism.

## Reproducibility, numerical invariance and source/privacy

- [Pure fd code](../../src/fd/genus_shapley_functional_volume.py) and [synthetic Shapley/additivity tests](../../tests/test_genus_shapley_functional_volume.py) operate on original Direct categorical phenotype arrays. Two-player Shapley is exactly additive in each 30-species sample, and source accepted taxon names are disjoint before sampling. A compiled per-species 100-bit original source reference lookup changes **only computation speed**, not the metric.
- Source-protected original runner: [island/scripts/run_fd_genus_shapley_functional_volume.py](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_genus_shapley_functional_volume.py); [original CI workflow](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/.github/workflows/run-fd-genus-shapley-functional-volume.yml).
- Original [source-run #43385444722](https://github.com/zuizui0223/island/actions/runs/43385444722) and [compiled-look-up rerun #43386044096](https://github.com/zuizui0223/island/actions/runs/43386044096) independently produced exactly **identical JSON bytes** with SHA256 **97e23194924a0e70011ff37bf3f7b6a1453fda9bee3a347985bbfb8a6d785057**; latest anonymous [artifact #14019250011](https://github.com/zuizui0223/island/actions/runs/43386044096/artifacts/14019250011). Both jobs SUCCESS. No original accepted plant species names, GIFT source entity IDs, genus names, geographic coordinates or island identifiers are published.

**Interpretation ceiling:** source-relative **floral categorical reference-space genus-composition accounting**, consistent with but not a causal demonstration of lineage sorting; no actual functional pollination niche, adaptation, partner loss or effective replacement measured.
