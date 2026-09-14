# Chapter 1 v8 figure architecture — 2026-09-14

## Design objective

The figures must tell the same inferential story as v8 without forcing the reader to learn the analysis names.

The main-text sequence is:

1. **What could generate the apparent island floral syndrome?**
2. **Does the response actually point in one direction globally?**
3. **At what taxonomic depth does the strongest response live?**
4. **Which alternative explanations survive cross-examination?**

This is deliberately different from organizing figures by analysis module or workflow version.

---

# Main Figure 1 — Adaptation or hierarchical assembly?

Canonical design spec:

- `docs/chapter1_figure1_hierarchical_syndrome_spec_20260914.md`

Purpose: establish the rival generators and evidence boundary before any detailed result.

Four panels:

A. original island-flower question and rival generators;
B. context branching and family→genus attenuation;
C. falsification wall;
D. Chapter 1 macroecology → Chapter 2 local mechanism scale bridge.

Figure 1 should contain only the minimum numerical anchors needed to orient the reader:

- `4/4 -> 4/4 -> 0/4`;
- genus attenuation `78.8–85.9%`;
- H5c `p=0.412`;
- H5d `0/8 identifiable`.

Do not turn Figure 1 into the full result display.

---

# Main Figure 2 — Biogeographic branching of the floral/reproductive response

## Biological question

**Does increasing source separation generate the same floral/reproductive trajectory across biogeographic contexts?**

This is the H1/H2 figure and should visually prove that the key result is *vector direction*, not one region being significant and another not.

## Panel A — response-vector phase plot

Axes:

- x = isolation slope for `accessibility_generalization`;
- y = isolation slope for `reproductive_assurance`.

Each context is one point/vector from the origin. Use direct between-context uncertainty if available; if confidence ellipses are not directly estimable under the frozen output, use separate horizontal/vertical confidence intervals rather than constructing unapproved ellipses.

Primary labels:

- Palearctic;
- tropical;
- other context points only if they pass the frozen support gate.

Visual quadrants:

```text
                    assurance ↑
                         |
 generalization ↓       |      generalization ↑
                         |
-------------------------+--------------------------
                         |
                         |
                    assurance ↓
```

The reader should immediately see that Palearctic and tropical responses are not scalar-strength variants of one vector.

## Panel B — evidence-scope / floristic-stratum replication

Compact forest or dumbbell plot for the two primary components across:

- all-analysis / all-native;
- direct-only / all-native;
- all-analysis / native-nonendemic;
- direct-only / native-nonendemic.

Do not show every sensitivity model. The aim is to show that the Palearctic direction is not created by low-confidence evidence or island endemics.

Key numerical anchors already frozen in v8:

Palearctic all-analysis all-native:

- accessibility/generalization `0.0795`, q=`0.00463`;
- reproductive assurance `0.0534`, q=`0.00463`.

Palearctic direct-only all-native:

- accessibility/generalization `0.0616`, q=`0.0389`;
- reproductive assurance `0.0966`, q=`8.48e-14`.

Tropical direct-only all-native:

- accessibility/generalization `-0.1005`, q=`0.00275`;
- reproductive assurance `0.1360`, q=`0.00428`.

## Panel C — component decoupling diagnostic

Show the attraction/access response before and after conditioning on `selfing_core` across the four frozen source definitions.

Primary message:

> the floral-architecture shift is not reducible to the measured reproductive-assurance component.

Annotate the all-native conditional distance range only once:

- approximately `0.091–0.101`;
- q=`0.0079–0.034`.

Avoid calling this mediation.

## What Figure 2 must not do

- do not put bee/butterfly/bird icons on the primary vectors;
- do not imply that tropical accessibility is equally robustness-secure as the Palearctic branch;
- do not use endemicity as the main replication axis;
- do not replace direct vector comparison with stars on separate regressions.

---

# Main Figure 3 — The syndrome has an assembly depth

## Biological question

**Where in the taxonomic hierarchy is the strongest Palearctic response expressed?**

This is the flagship explanatory figure.

## Panel A — attenuation trajectory

For each of the four primary scope × floristic-stratum cells, plot normalized vector magnitude through:

`observed -> family-adjusted -> genus-adjusted`.

Normalize each trajectory to observed = 1 only for visual comparison; also report the raw/absolute effect metric in the source table or supplement.

Expected visual:

- modest drop observed→family;
- large drop family→genus.

Directly annotate:

- family attenuation: `19.6–33.4%` of observed vector;
- genus attenuation: `78.8–85.9%` of observed vector;
- conditional family→genus attenuation: `70.6–79.1%` of remaining signal.

## Panel B — support-state ladder

A small matrix:

| scope × stratum | observed | family | genus |
|---|---:|---:|---:|
| all-analysis × all-native | pass | pass | fail |
| all-analysis × NNE | pass | pass | fail |
| direct-only × all-native | pass | pass | fail |
| direct-only × NNE | pass | pass | fail |

Headline annotation:

**`4/4 -> 4/4 -> 0/4`**

Do not make `fail` look like a software error. Label as `vector gate not retained` in the final artwork.

## Panel C — hierarchical depth differs by context

Use the independent hierarchical-depth audit only within the scope where its claim ceiling is valid.

Safest headline:

> Within widespread native non-endemic island assemblages, taxonomic attenuation of isolation-associated floral architecture differs between northern-midlatitude and tropical contexts.

Do not generalize this panel to all natives, endemics or evolution.

## Panel D — interpretation schematic

Small causal-neutral diagram:

```text
source-available lineages
        ↓
arrival / establishment / persistence / ecological filtering
        ↓
differential genus representation
        ↓
assemblage trait syndrome
```

Beside it, list mechanisms that remain unresolved:

`dispersal | habitat | demography | interaction dependence | history`

This makes clear that H3 localizes the response but does not name the cause.

---

# Main Figure 4 — Cross-examination: what survives and what does not identify mechanism?

## Purpose

Figure 4 is not a “robustness dump”. It should show why the final interpretation is narrower but stronger after adverse evidence is included.

Use four quantitative panels, each with one question.

## Panel A — species-list detection tipping (V6)

Show two matched heatmaps or tipping maps:

- Palearctic accessibility;
- tropical accessibility.

Axes:

- distance-dependent completeness decline (`OR_C` or the frozen completeness parameterization);
- state-dependent recording odds (`OR_D`).

Overlay the support/tipping boundary.

Headline:

- Palearctic: `99/100` baseline-supported surfaces survive; `80/80` remote-under-survey scenarios survive;
- Tropical: `40/75` baseline-supported surfaces tip.

This visually justifies asymmetric confidence across contexts.

## Panel B — calibrated response geometry

For each of the 12 broad response cells, plot observed nonlinear evidence statistic relative to its predeclared calibrated critical value.

Preferred scale:

`observed D - critical_D`

with zero as the promotion line.

All observed points should remain at or below zero under the frozen result.

Headline:

**0/12 promoted nonlinear response shapes.**

A small inset may show qualification performance for G2 step to prove that the negative result is not simply lack of detectability:

- G2 cross-scope qualification in 11/12 cells;
- the one miss was `0.796875` recovery against the `0.80` rule and was not rescued.

## Panel C — independent biotic vs wind specificity (H5c)

Single forest-style estimate:

`distance × biotic = +0.06495`

95% CI `[-0.09030, 0.22020]`, p=`0.41221`.

Side annotation:

- only 1/8 planned cells prospectively qualified;
- 5,771 unambiguous independent mode species;
- observed opening performed once.

This panel is intentionally small: it is a claim boundary, not a new headline.

## Panel D — distributed-threshold identifiability (H5d)

Plot, across the eight design cells:

- classification accuracy;
- false distributed-threshold selection under smooth clines.

Reference lines:

- accuracy requirement `0.80`;
- false-selection ceiling from the frozen contract.

Headline:

- `0/8` qualified;
- accuracy `0.733–0.778`;
- false threshold calls `0.19–0.255`.

Interpretation ribbon:

> global assemblage data cannot identify a distributed-threshold generator; local mechanistic resolution is required.

---

# Extended Data / Supplement figure allocation

## ED Figure 1 — trait evidence and missingness anatomy

- coverage by raw axis;
- evidence scope composition;
- taxonomic / geographic missingness summaries;
- V5 finite MNAR tipping surface.

Purpose: defend observation process without spending a main figure on data acquisition.

## ED Figure 2 — area/capacity moderation

- all 16 H4 cells;
- heteroskedastic-null comparison;
- measurement-sensitive modifier classification.

Headline: `0/16` promoted mechanism.

## ED Figure 3 — response-geometry qualification

- V1 naive selector failure under clustered clines;
- V2 calibrated critical-D distribution;
- held-out false nonlinear maximum <=0.088542;
- G2/G3/G4 recovery summary.

This is necessary to defend the credibility of the main-text `0/12` nonlinear result.

## ED Figure 4 — pollination-associated architecture factor

- shared factor loadings / explained variance;
- all-analysis `86.94%`;
- direct-only `86.44%`;
- template-specific residuals only where frozen support allows.

Purpose: demonstrate why named syndrome templates cannot be interpreted as visitor identities.

## ED Figure 5 — N1 independent channel audit

- descriptive channel slopes;
- joint test `p=0.65516`;
- power simulation showing limited sensitivity to modest channel heterogeneity;
- LOBO/source-deletion stability.

Purpose: show non-identification rather than burying the null.

## ED Figure 6 — source breadth / GloBI audit

Show effort qualification and 0/4 promotion result. Keep this out of the main story unless a reviewer specifically asks about source-side partner breadth.

---

# Main-text data economy

The paper should not attempt to put every successful or failed analysis in the main figures.

Main figures earn their place by answering four sequential questions:

1. **What are the rival generators?** — Figure 1
2. **Does the syndrome branch by context?** — Figure 2
3. **At what assembly depth is the strongest branch expressed?** — Figure 3
4. **Does the interpretation survive missingness and independent mechanistic cross-examination?** — Figure 4

Everything else supports one of those four questions in Extended Data or Supplement.

---

# Poster / talk reuse

The same architecture can be compressed for the November Island Biology poster:

- Figure 1A + Figure 2 phase plot = `WHERE / WHAT DIRECTION`;
- Figure 3 attenuation = `AT WHAT ASSEMBLY DEPTH`;
- Figure 4D + `izu-core` local chain = `WHY GLOBAL DATA STOP / WHY CHAPTER 2 IS NEEDED`.

For the poster, do not reproduce the full falsification matrix. Use one compact “cross-examination” strip and spend visual space on the Figure 2 vector contrast and Figure 3 family→genus drop.
