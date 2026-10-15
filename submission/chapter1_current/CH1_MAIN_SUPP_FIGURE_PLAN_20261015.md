# Chapter 1 Q1-only paper: Main/Supp and figure migration plan (15 October 2026)

**Editorial decision implemented in the working manuscript on this PR branch; not a new numerical analysis.**
The submitted paper is **Q1 only**: observational island floral and reproductive data plus a separately published pollen-supplementation synthesis. Ch2/Q2 evolution, selection and persistence simulations belong to **a separate `izu-core` paper**. Neither the simulation trajectories nor their Q1–Q4 results are presented as empirical validation of this Chapter 1.

## Research question and inferential ceiling

**Do pollen-supply constraints underlie floral island syndromes?**

Classic *proposed* route: isolation → reduced effective pollen/mate service → colonization filtering or selection → self-compatibility/autonomous selfing and possibly altered colour, form and attraction.

We directly measure **none of those historical arrows together in the same island species through time**. Current pollen limitation (PL) is the reproductive improvement after pollen addition, not direct pollinator numbers, disappearance or a pre-selection exposure; selfing and other extant traits can *already lower measured PL*. Our conclusions therefore concern three empirically distinct sets of associations rather than an estimated proportion of floral-change variance causally explained by PL.

**Stage 1 — Constraint (source H3):** published GloPL: isolation–PL; 2,969 experimental rows, 1,248 sites, 919 publication clusters. Main β=+0.09191, finite-publication p=0.01594, but supplemental-only p=0.31082; offshore-only β=+0.22031, p=0.02459 is post-hoc.

**Stage 2 — Trait response (sources H1 + H2, the principal new database contribution):** 106,295-species trait universe and 4,379 sampled broad-H1 islands (not 106,295 fully measured species). Fit seven reproductive, colour and structural outcomes separately in four regions. The observed four positive self-compatibility isolation slopes do not imply that mating system, autonomous selfing or flower traits universally move in the same direction. Test the selfing-syndrome explanation *inside this section*: the source-H2 isolation coefficient for floral accessibility/colour after `selfing_core` adjustment (not causal mediation). H2's `selfing_core` predictor coefficient is exploratory when read as R–F.

**Stage 3 — Functional compatibility (source H4):** independent trait-score–GloPL associations, post-hoc. `selfing_core`: 455 matched species / 409 publications, β=−0.29830, p=0.00417. `generalized_accessible`: 143 species / 143 publications, β=−0.29566, p=0.02334. The two match sets differ and are not paired tests of a buffering mechanism. Previously frozen distance × trait interaction tests did not establish moderation of the isolation–PL slope.

## Proposed Main display set: **five** figures (within six-item Letter cap)

The five Q1 Main figures have now been rendered and committed as PDF/SVG pairs under `figures/q1_main5/`, with source figures 1–6 preserved separately. Numeric H1–H4 source IDs remain immutable; the new Main figure numbers are 1–5. Automated file/manifest checks have passed, while manual visual QA is pending.

| Future main | Main visual / one question | Existing evidence for graphical reuse | What remains in Supp |
| --- | --- | --- | --- |
| **Fig. 1: Study design and distinct data sources** | One restrained map, data-coverage statistics, and the dashed *hypothesized* D→pollen-service→R/F causal path. Show island trait assemblages and GloPL sites as **separate** sample streams. | Source Figs 1+2+3 | Exact coastline repair, regional support, full database derivation / WCVP rights data: S1–S2 |
| **Fig. 2: Isolation and experimental PL** | Source H3 observed site/effect data and fitted association, with primary β and mandatory supplemental-only uncertainty. Do not depict % of H1 explained. | Source Fig 6A | Measurement-category, mainland-zero, offshore-only, leave-one-study-out details: S5 |
| **Fig. 3: Reproductive and floral responses** | **Largest data figure**: the 7 atomic outcomes × 4 regions separated into 3 reproductive and 4 floral response rows. Show effect direction and pointwise CI, not a pooled universal syndrome score. | Source Fig 4 | All × Direct × WCVP full estimates, raw trait evidence, origin uncertainty, multiplicity caveat: S2–S3 |
| **Fig. 4: Selfing syndrome versus additional floral association** | Source H2: island-level `selfing_core`, `generalized_accessible`, and `plain_colour`. Emphasize the **conditional D coefficient** after R adjustment and show the response definition in each panel. If an unadjusted-vs-adjusted contrast is desired, refit *the same* response on *identical islands* first; never compare H1 atomic outcomes with H2 aggregate scores as if identical. | Source Fig 5 | All/Direct q, complete-three-component sensitivity, raw colour × morphology, exploratory R coefficients, future area/climate ablation: S4 |
| **Fig. 5: Functional association with current PL** | Source H4 scores against PL as two distinctly labelled matched analyses: `selfing_core` and accessibility, with denominators 455 and 143 species; post-hoc label visible. | Source Fig 6B–C | Atomic-trait estimates, supplemental-only, no-zero-constant, failed buffering/mediation tests: S6 |

**Figure-number migration was executed in PR #319.** `FIGURE_CAPTIONS.md` now documents the five current Q1 figures and Main/SI cross-references use their new numbering. Old source-six captions are in `history/`. CI verifies file integrity and frozen-source mapping; manual figure legibility review and journal-format assembly remain incomplete.

## Main sections

- **Introduction:** classical loss-of-effective-pollination hypothesis; why PL is a better reproductive endpoint than pollinator-distribution inference but is not direct pollinator decline; alternative selfing, parallel floral and species-sorting routes.
- **Results 1 / ecological constraint:** short independent GloPL isolation–PL observation.
- **Results 2 / database centerpiece:** seven atomic isolation responses, then source-H2 conditional floral response in the same biological section.
- **Results 3 / functional consistency:** independent post-hoc matched GloPL score regressions.
- **Discussion:** PL is partly trait-dependent; test serial selfing syndrome versus parallel pollinator-facing associations versus source/colonization filtering; explicitly distinguish accessible flowers from dull or low-cost flowers; all historical directional arrows remain hypotheses.

## SI organization and limits

- **S0 source crosswalk** (added to the working SI) maps the three reader stages to original H3, H1+H2, H4; **six-edge D/R/F/P scheme becomes a descriptive relationship inventory, not six coequal confirmatory tests**.
- **S1–S2** coastline, island sample and trait evidence/missingness/redistribution; emphasize WCVP *regional* native compatibility only.
- **S3–S4** the most detailed supplementary material: seven atomic slopes and test rules; complete selfing adjustment; raw floral colour/structure; area/climate sensitivity **if computed and validated**.
- **S5–S6** GloPL geography and exact species bridge, unmatched denominators, experimental study controls, post-hoc limitations and unconfirmed moderators.
- **S7** actual hash-locked sources, software/data availability and submission rights gate. Do not imply the 46,274 public rights-cleared cells recreate all 222,688 scientific cells.

## Completion gates

1. Manuscript: three-stage text plus SI crosswalk updated; preserve frozen H1–H4 tables, H1 pointwise unadjusted tests, source H2 FDR and H4 post-hoc status.
2. Figures: physically regenerate five Main figures; move audit/extra panels to SI; update captions and **all numeric figure cross-references** without using synthetic data or unsupported predicted proportions.
3. Inference: validate original coefficients and sample denominators, within-5000-word Main / within-150-word abstract; new R–F claims require proper exploratory status and multiplicity declaration.
4. Submission: do not mark ready until figure assets, author metadata, code DOI and Ecology Letters third-party trait redistribution policy are actually resolved.

*The separate `izu-core` mechanism paper can cite this paper as motivation, but no Ch2 model result is to be imported into the Ch1 Main or Supp.*
