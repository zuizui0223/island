# Geological-origin test of global island floral phenotype-space structure

**8 November 2026 — independent exploratory analysis in** [`zuizui0223/fd`](https://github.com/zuizui0223/fd).  
**Scientific status:** repeatable descriptive and adjusted associations. **No confirmed geological-origin effect beyond fixed profile margins, no evolutionary-causal or pollination mechanism claim.**

## 1. Prior art and the testable novelty

Analysis of island plant functional trait space is **not novel by itself**. Existing work includes: Tenerife's functional-space assembly (Barajas Barbosa et al., 2023, [Nature](https://doi.org/10.1038/s41586-023-06305-z)); source-pool filtering of island flora traits (Schrader et al., 2024, [J. Biogeogr.](https://doi.org/10.1111/jbi.14703)); global colonization selection involving reproductive system and floral symmetry (Zell et al., [New Phytol.](https://doi.org/10.1111/nph.20234)); and the global covariation of 15 angiosperm taxonomic traits (>320,000 species, [Nature Plants 2026](https://doi.org/10.1038/s41477-026-02358-w)). Geological-origin differences in species–area responses are also known (Cáceres-Polgrossi et al. 2024, [Community Ecology](https://doi.org/10.1007/s42974-024-00195-5)). Roeble et al. (2024) explicitly compared 272 oceanic and 238 continental islands for Asteraceae richness/endemism ([Nature Communications](https://doi.org/10.1038/s41467-024-51556-7)).

**The potential new question is the joint, distance-aware test** of local **floral** phenotype dispersion (colour and form) and between-island functional separation **by independently established geological island type**, using matched support, profile-fixed conditional nulls, covariate/lineage controls, and missingness audits. This is an identifiable question, **not yet an established first discovery**: existing phenotype data do not directly measure pollination effectiveness and source/genus history is not conditioned out.

## 2. Geology and glacial connectivity are different variables

The island base is the 8,264-unit corrected GSHHG geocoded population at frozen [island PR#321](https://github.com/zuizui0223/island/pull/321), commit `d4a5dafab5f0737ea6c1c1c0263fc20decf6035c`. Independent geological labels were drawn from:
- Roeble et al. 2024 [Supplementary Data 3](https://www.nature.com/articles/s41467-024-51556-7#Sec28) containing 510 curated geological island labels (original 272 oceanic, 238 continental). Source table is downloaded at execution only, not committed.
- Zizka et al. 2022 [Zenodo input archive](https://doi.org/10.5281/zenodo.6325640), with 83 independently sourced additions to missing GIFT geology. Mixed records are not forced to either geological class.
- Weigelt et al. 2013 [world-island environments](https://doi.org/10.1073/pnas.1306309110) and [open table](https://zenodo.org/records/4971536) separately provide `GMMC`: **LGM mainland-connected / not mainland-connected**. An unconnected continental fragment is NOT geologically oceanic.

All joins require centroid proximity, maximum area ratio 1.8, and no ambiguous nearby or reused source identity. They are **conservative spatial crosswalks, not validated identity matches or geological proof**. Records are retained as `oceanic`, `continental`, `mixed`, or `unresolved`. The 12 oceanic-labelled records with LGM-connected matches illustrate that source definitions/timescales differ; these are audited, not overwritten.

| Evidence | Matched across 8,264 GSHHG islands | Complete All scope (≥50 name units) | Complete Direct scope |
| --- | ---: | ---: | ---: |
| Geological classes `continental` or `oceanic` | 434 (176/258) | **260 (130/130)** | **243 (123/120)** |
| Additional `mixed` geological class | 8 | 8 | inspected separately |
| LGM mainland connectivity, connected vs not | 6,727 (4,575/2,152) | **1,334 (987/347)** | **1,252** |
| Remaining geology unresolved | 7,822 | 1,264 | retained, not imputed |

The geological sample is a **nonrandom subset** of the recorded floras because of source coverage and matching. Only 17 complete All islands had strict source labels before integrating the 2024 510-island study.

## 3. Observed flora-space differences (descriptive)

| All-scope response | Continental (130 islands) | Oceanic (130 islands) |
| --- | ---: | ---: |
| Equal-profile-weight Rao Q | **0.77447** | **0.75705** |
| Name-carrier-weighted Rao Q | 0.75078 | 0.74459 |
| 90% 2-D PCoA KDE area at h=0.35 | 14.340 | 14.242 |
| Within-class NN across upstream sampled pairs | 0.05033 (111 pairs) | 0.05945 (240 pairs) |

Direct-only equal-profile Rao Q also differs in the same direction: continental **0.77362** vs oceanic **0.75617**. These are descriptive medians, not causal effects or type-specific G1 significance tests.

**Regional sampling:** Northern midlatitude All has 85 continental /48 oceanic; tropical All has 31/76; southern extratropical 13/6; northern high 1/0. The pair result is especially sensitive to unequal cohort composition: within tropical pairs median NN is 0.05387 continental vs 0.05554 oceanic, much smaller than the global difference. Pair counts reuse islands and are not independent observations.

## 4. Exploratory covariate-adjusted contrasts

Response `~ oceanic vs continental + log(area) + log1p(generic continent distance) + log(recorded name richness) + joint trait coverage + climate PCs 1–4 + region`. Spatial-block Rademacher **wild** bootstrap with 999 draws provides unadjusted exploratory 95% coefficient intervals. Class is not randomized; missingness, source pools, floral ancestry and colonization history are not fully controlled.

| Scope | Response | Oceanic − continental coefficient | Spatial-block 95% interval |
| --- | --- | ---: | ---: |
| All | Equal-profile Rao Q | **−0.00589** | **[−0.00988, −0.00227]** |
| Direct | Equal-profile Rao Q | **−0.00496** | **[−0.00910, −0.00097]** |
| All | Carrier-weighted Rao Q | −0.00305 | [−0.00667, +0.00060] |
| Direct | Carrier-weighted Rao Q | −0.00176 | [−0.00621, +0.00274] |
| All | Projected KDE area | −0.235 | [−0.480, −0.004] |
| Direct | Projected KDE area | −0.253 | [−0.498, +0.005] |

The equal-profile directional association survives the listed covariates in All and Direct. The name-carrier-weighted response does **not** show a clear corresponding separation. The KDE endpoint is explicitly a low-fidelity two-axis projection retaining only **15.35%** of PCoA positive inertia. Interval exclusion of zero is **not corrected for multiplicity**, and overlapping All/Direct scopes are **not independent validations**.

The separate **LGM-connectivity** adjusted All equal-profile effect was −0.00363 with interval [−0.00786, +0.00064], so a consistent broad binary mainland-connection effect was **not established**.

## 5. Direct fixed-profile-margin origin contrast: stringent exploratory gate

The additional G1 null uses two 200-draw Curveball chains within **each origin × geographic regime × evidence scope**, preserving each island's count of recorded profiles and each profile's incidence frequency within that stratum. This asks **directly** whether the observed oceanic-minus-continental difference is larger than the class-specific binary-margin expectation. The comparison is intentionally not a genus/source-conditioned or taxon-excluded null. Regions without enough matched islands or pairs are retained as unsupported.

**12 of 32 fixed contrast slots were testable; none supports an origin-specific excess after BH or BY adjustment (0/12).**

| Equal-profile Rao Q contrast | Observed oceanic − continental | Null-median difference | Approx. p |
| --- | ---: | ---: | ---: |
| Northern midlatitude, All | −0.009095 | −0.009364 | 0.913 |
| Northern midlatitude, Direct | −0.009014 | −0.009000 | 0.998 |
| Tropical, All | −0.013136 | −0.014458 | 0.494 |
| Tropical, Direct | −0.014252 | −0.015119 | 0.718 |

Tropical whole-island pair NN origin contrast (All) observed +0.00514 vs null +0.00216 (p=0.429), with no FDR support. More remote/limited island pairs, especially northern midlatitudes, do not meet the predeclared ≥20-pairs-per-origin gate.

**Interpretation:** A lower floral Rao Q on oceanic islands is repeatedly observed, but the stronger test indicates that its **class-specific profile prevalence and profile-count structure can reproduce the observed local contrast**, without invoking additional within-island packing beyond those margins. This does not prove why the profile prevalence differs, nor establish origin-caused filtering, ecological equivalence or pollinator loss. The upstream [island PR321](https://github.com/zuizui0223/island/pull/321) taxon-excluded genus-conditioned non-rejection remains unaffected.

## 6. What to do next scientifically

The strongest follow-up is **profile-repertoire assembly rather than an assumed universal volume reduction**: identify which flower-colour/form combinations differ in prevalence between independently sourced continental/oceanic floras while preserving genus/species incidence, source-region availability, name-level coverage and origin classification quality. Then compare that repertoire-filtering explanation to conditional within-island packing and direct pollinator/reproductive data.

As presently supported, the bounded conclusion is:

> *Recorded island floras differ in local floral phenotype dispersion between independently curated oceanic and continental classes, but the difference is not demonstrably greater than expected from class-specific categorical profile incidence margins. Whether the repertoire shift reflects biogeographic sources, lineage assembly, observation bias or ecological filtering remains unresolved.*

## Provenance and validation

- [Original island floral profile input and methods](https://github.com/zuizui0223/island/pull/321)
- [fd G0/G1 and KDE background report](../20261108_g1_hypervolume/README.md)
- [Original matched LGM and Roeble geography source contracts](../../docs/NOVELTY_LITERATURE_AUDIT.md)
- Original crosswalk script: `src/fd/origin.py`; response analyses: `src/fd/origin_analysis.py`, `src/fd/origin_null.py`.
- GitHub Actions source run [#42441658408](https://github.com/zuizui0223/fd/actions/runs/42441658408), then [#42441776287](https://github.com/zuizui0223/fd/actions/runs/42441776287) for All/Direct/single-state adjusted sensitivity; 15 tests passed. Results include `stratified_island_medians.csv`, `stratified_pair_medians.csv`, `adjusted_origin_contrasts.csv`, `origin_g1_groups.csv`, `origin_g1_contrasts.csv`.
- Expected hashes for the latest adjusted analysis: `adjusted_origin_contrasts.csv` SHA256 `112dd4d4b3bfe2ba6fe484820ee9fa9e42a5c6e91dac68012715487c0ead3e73`; `origin_g1_contrasts.csv` (the 12-test run) `bfc8b476250bf666c9ef4a26d773f3edc359a09dcf7f076c55ced2ad3194dabd`.
- Raw GIFT/Roeble original source tables and the island-level crosswalk are **not committed**: the repository retains aggregate analysis and source hashes to respect redistribution rights.
