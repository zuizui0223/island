> **HISTORICAL / SUPERSEDED — pre-corrected Chapter 1 surface.** Retained for provenance/replay only. The current submission is selected by `config/chapter1_submission_current.json`; use `submission/chapter1_current/MANUSCRIPT.md` and `docs/PAPER_PIPELINE.md` for current results.

# Chapter 1 v8 submission figure synchronization — 2026-09-15

This document is the manuscript-facing synchronization layer for the **actual rendered and locked main Figures 1–4**. It supersedes earlier figure-layout sketches when those sketches differ from the final rendered panel arrangement. It does not change any biological estimate.

## Canonical figure sequence

The main text must read in this order:

1. **Figure 1 — inference map:** what can generate an island floral syndrome, what the global data establish, what survives cross-examination, and why local mechanistic resolution is still needed.
2. **Figure 2 — biogeographic branching:** the primary H1/H2 response changes direction across contexts and the Palearctic branch replicates across evidence scopes/strata.
3. **Figure 3 — assembly depth:** the strongest Palearctic response attenuates most strongly from family to genus.
4. **Figure 4 — cross-examination:** observation bias, nonlinear geometry, independent pollination specificity and distributed-threshold identifiability define the final claim boundary.

The older eight-component response-fingerprint atlas remains a descriptive/Extended Data visualization and must not replace the primary H2 Figure 2.

---

## Figure 1 — From an island floral syndrome to a hierarchy-of-assembly test

Canonical lock:

- `config/chapter1_v8_figure1_result_lock.json`
- run `34862804485`
- artifact `10355018602`
- digest `sha256:551e5f9738ca0193083232c27ba1abbdccc92e718c9553f175d68c1cde4c5e66`

Actual panels:

- **A — rival generators:** within-lineage response, hierarchical lineage sorting, or both can produce an assemblage-level floral syndrome.
- **B — global inference:** universal response rejected; context branching followed by `4/4 -> 4/4 -> 0/4` taxonomic attenuation, with genus attenuation `78.8–85.9%` and conditional family→genus attenuation `70.6–79.1%`.
- **C — cross-examination wall:** V5/V6, H4, nonlinear geometry, N1, GloBI breadth, H5c and H5d are shown as inferential outcomes rather than software successes/failures.
- **D — scale bridge:** Chapter 1 identifies macroecological direction/context/assembly depth; `izu-core` resolves interaction state → effective service → reproductive outcome → phenotype and local response geometry.

### Final Figure 1 legend

**Figure 1 | From an island floral syndrome to a hierarchy-of-assembly test.** An isolation-associated shift in floral or reproductive composition can arise through repeated within-lineage response, differential assembly of source-available lineages, or both (A). The global analysis rejects one universal floral/reproductive response, identifies biogeographic branching, and then localizes the strongest Palearctic branch by taxonomic depth (B). The broad Palearctic response remains supported after family adjustment but fails the predeclared vector gate after source-matched genus adjustment (`4/4 -> 4/4 -> 0/4`); genus adjustment attenuates `78.8–85.9%` of the observed vector and most conditional attenuation occurs from family to genus. Independent robustness and falsification analyses then define the claim boundary rather than supplying a post hoc mechanism (C): the Palearctic branch survives strong specified missingness challenges, while area, a common global nonlinear shape, coarse pollinator-channel heterogeneity, sampled interaction breadth and independent biotic-versus-wind specificity do not yield a promoted global mechanism; distributed thresholds are not identifiable against heterogeneous smooth clines in the present design. Chapter 1 therefore identifies response direction, biogeographic context and assembly depth without identifying a pollinator-specific causal chain. A deeply resolved local system is required to distinguish interaction, effective-service, reproductive and phenotypic response geometries (D). Dashed connections indicate mechanistic links not identified by the global data.

---

## Figure 2 — Biogeographic branching of the primary plant response

Canonical lock:

- `config/chapter1_v8_figure2_result_lock.json`
- run `34810404509`
- artifact `10334412842`
- digest `sha256:c177e28ac75a8bccf00e9f0e1a467f2958907ba3c3eae87334a62e5924780675`

Actual panels:

- **A — all-native phase plot:** accessibility/generalization slope × reproductive-assurance slope for northern-midlatitude and tropical contexts, with all-analysis and direct-only estimates.
- **B — native-nonendemic phase plot:** the same primary vector in the floristic stratum that removes island endemics from the response.
- **C — Palearctic replication:** the two primary components across all/direct and all-native/native-nonendemic cells.
- **D — pathway decomposition:** Palearctic floral-architecture distance response after conditioning on `selfing_core` across four frozen source definitions.

### Final Figure 2 legend

**Figure 2 | Floral and reproductive responses branch among biogeographic contexts.** Primary response vectors are shown as the isolation-associated slopes of accessibility/generalization and reproductive assurance in all-native (A) and native-nonendemic (B) assemblages. Open and filled symbols distinguish all-analysis and direct-only evidence rather than separate hypotheses. The Palearctic branch is reproduced across evidence scopes and floristic strata (C), whereas tropical assemblages can combine increasing reproductive assurance with decreasing accessibility/generalization, so the tropical pattern is not a weaker scalar copy of the Palearctic response. The Palearctic floral-architecture association also remains positive after conditioning on the measured `selfing_core` across four frozen source definitions (D), showing component decoupling rather than causal mediation. The figure visualizes frozen H1/H2 estimates; direct between-context multivariate tests remain the inferential basis for biogeographic heterogeneity.

---

## Figure 3 — The strongest Palearctic syndrome has an assembly depth

Canonical lock:

- `config/chapter1_v8_figure3_submission_result_lock.json`
- run `34816668470`
- artifact `10336875587`
- digest `sha256:853fd1c590d86341d6302c1da8e998af93b0b5fa3b35c464dd3762378385bd4d`

Actual panels:

- **A — 16 taxonomic attenuation profiles:** normalized Palearctic response-vector magnitude through observed → family-adjusted → genus-adjusted stages across source modes, evidence scopes and floristic strata.
- **B — accessibility/generalization stage effects:** canonical source-mode estimates across observed, family and genus stages.
- **C — reproductive-assurance stage effects:** the corresponding `selfing_core` stage estimates.
- **D — vector gate matrix:** all four evidence-scope × floristic-stratum cells show `4/4 -> 4/4 -> 0/4`.

Frozen quantitative anchors:

- family attenuation `19.6–33.4%` of observed vector;
- genus attenuation `78.8–85.9%` of observed vector;
- conditional family→genus attenuation `70.6–79.1%` of remaining signal.

### Final Figure 3 legend

**Figure 3 | The strongest Palearctic floral/reproductive response is concentrated at the family-to-genus assembly transition.** Across 16 frozen source-mode × evidence-scope × floristic-stratum profiles, normalized response-vector magnitude declines modestly after family composition is accounted for and much more strongly after source-matched genus composition (A). The same progression is shown separately for accessibility/generalization (B) and reproductive assurance (C) under the canonical source comparison. The predeclared multivariate gate is retained at the observed and family-adjusted stages but not after genus adjustment in all four primary evidence-scope × floristic-stratum cells (`4/4 -> 4/4 -> 0/4`; D). Genus adjustment attenuates `78.8–85.9%` of the observed vector; after family adjustment, genus composition removes `70.6–79.1%` of the remaining vector. These quantities localize hierarchical expression and are not causal mediation: they do not prove dispersal alone, exclude within-lineage change, or identify why particular genera are sorted.

---

## Figure 4 — Cross-examination and claim boundaries

Canonical lock:

- `config/chapter1_v8_figure4_result_lock.json`
- run `34817349228`
- artifact `10337380746`
- digest `sha256:f41a90089ddda1d285186ca441d416d71cf1c276e742a7da917dbb572f39fbb0`

Actual panels:

- **A — V6 survival fraction:** Palearctic accessibility `99/100`, North–Tropical vector `70/75`, tropical accessibility `35/75` survive the complete frozen species-detection bias grid among baseline-supported surfaces.
- **B — calibrated nonlinear margins:** observed nonlinear evidence relative to each scope-specific critical value; no cell enters the all-analysis + direct-only promotion quadrant (`0/12`).
- **C — H5c specificity:** independently assigned biotic-vs-wind pollination mode gives `distance × biotic = +0.06495`, 95% CI `[-0.09030, 0.22020]`, `p=0.41221`; not promoted.
- **D — H5d identifiability:** distributed-threshold vs heterogeneous-cline classification fails the frozen gate in all `0/8` cells; false threshold selection under smooth clines remains approximately `19–25.5%`.

### Final Figure 4 legend

**Figure 4 | Cross-examination narrows the interpretation without erasing the plant-side result.** Species-list detection sensitivity is asymmetric across contexts (A): the Palearctic accessibility branch survives `99/100` baseline-supported bias surfaces, whereas the tropical accessibility branch survives `35/75`; the direct North–Tropical multivariate contrast is more robust (`70/75`). Calibrated response-geometry tests require nonlinear evidence to exceed a design-specific threshold in both all-analysis and direct-only scopes; none of the 12 observed cells meets that promotion rule (B). An independent biotic-versus-wind specificity test in the sole prospectively qualified cell yields a positive but imprecise interaction (`+0.06495`, 95% CI `-0.09030` to `0.22020`, `p=0.41221`) and does not identify a pollinator-specific global mechanism (C). Finally, distributed lineage thresholds cannot be distinguished reliably from heterogeneous smooth clines in the realized global design (`0/8` qualified; D). Negative or non-identified results are not inverted into evidence that pollinators, local thresholds or other mechanisms are irrelevant.

---

## Required manuscript cross-references

The canonical v8 manuscript should use the following placements when final prose is synchronized:

1. **Introduction, sequential inferential problem:** cite **Fig. 1** after the four-stage question sequence.
2. **Results H1/H2:** cite **Fig. 2A,B** after the statement that response direction differs among contexts; cite **Fig. 2C,D** after Palearctic replication and selfing-core decomposition.
3. **Results H3:** cite **Fig. 3A–D** after the quantitative attenuation paragraph (`19.6–33.4%`, `78.8–85.9%`, `70.6–79.1%`, `4/4 -> 4/4 -> 0/4`).
4. **Results response geometry:** cite **Fig. 4B** after `0/12` promoted nonlinear shapes.
5. **Results V5/V6:** cite **Fig. 4A** after the `99/100`, `35/75` and `70/75` species-detection sensitivity comparison.
6. **Results H5c/H5d:** cite **Fig. 4C,D** after the independent specificity and distributed-threshold non-identification results.
7. **Discussion opening:** cite **Figs. 2–4** after the three-result summary (context branching; family→genus attenuation; failed simple mechanism tests).
8. **Discussion Chapter 1 / Chapter 2 scale bridge:** cite **Fig. 1D** after the statement that global assembly depth and local response geometry answer different levels of the same problem.

## Claim ceiling retained across all four figures

- no `distance -> pollinator loss -> floral simplification` causal chain;
- genus attenuation does not prove dispersal alone or evolution absent;
- H5c `p=0.412` does not prove pollinators irrelevant;
- tropical accessibility is not presented as equally robustness-secure as the Palearctic branch;
- `0/12` global nonlinear promotion does not rule out local thresholds;
- `0/8` H5d identifiability keeps observed genus threshold distributions closed;
- endemicity is not used as a time axis.
