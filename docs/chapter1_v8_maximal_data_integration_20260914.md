# Chapter 1 v8 maximal-data integration — 2026-09-14

## Status

This note updates the v8 manuscript-writing target after PR #214. It does **not** change the frozen H1–H5 inferential contract. It specifies how to use more of the already frozen data without turning post-freeze descriptive synthesis into new confirmatory biology.

The manuscript should no longer treat H1/H2/H3, taxonomic depth, geometry, H4, N1 and GloBI as disconnected analyses. They are different dimensions of one isolation-response fingerprint.

## Revised one-sentence claim

> **Island isolation is associated with biogeographically distinct floral and reproductive response vectors whose component composition and hierarchical expression differ across contexts; the strongest Palearctic response is concentrated at the family-to-genus assembly transition, while several simple mechanistic shortcuts fail explicit promotion gates.**

This is stronger than “the island floral syndrome is not universal” but does not claim that the mechanism has been identified.

## Abstract replacement logic

Use five moves only.

### 1. Problem

Island floras are commonly summarized by a coherent floral island syndrome, but that framing conflates response direction, component coupling and the taxonomic level at which assemblage change is represented.

### 2. Scale and contract

Retain the fixed universe of 8,265 islands, 106,295 accepted angiosperms, 318,885 possible species × axis cells and 222,688 resolved cells. State that the progressive contract froze hypotheses, source-pool safeguards, model order and claim ceilings.

### 3. Direction and component composition

Retain the frozen H1/H2 headline:

- Palearctic separation: accessibility/generalisation ↑ and reproductive assurance ↑;
- Tropical direct-only: accessibility/generalisation ↓ while reproductive assurance ↑.

Then add one short descriptive sentence:

> Across the same eight atomic colour, structural and reproductive contrasts, northern-midlatitude and Tropical response vectors were nearly orthogonal to moderately opposed (95.3–116.5°), showing that regional differences are not well described as the same syndrome expressed at different amplitude.

The angle is descriptive and must not replace the frozen direct H1/H2 heterogeneity test.

### 4. Hierarchical expression

Replace vote-count-only language with effect-size attenuation:

> In the Palearctic, family adjustment reduced the primary two-axis response by 19.6–33.4%, whereas source-matched genus adjustment reduced the original response by 78.8–85.9%. Genus adjustment therefore removed 70.6–79.1% of the family-adjusted remainder, after which the vector no longer passed the predeclared post-genus gate.

Immediately add the bounded generality result:

> A matched context × taxonomic-stage audit showed that genus adjustment altered the same plant-architecture response differently between northern-midlatitude and Tropical contexts across all four source definitions in native non-endemics, although the result weakened when island endemics were included.

### 5. Explanation elimination

Compress the negative alternatives into one sentence rather than four separate null stories:

> Neither apparent area amplification, a globally reproducible step/reversal, prospectively tested pollinator-channel heterogeneity, nor an effort-matched source-genus GloBI interaction-breadth proxy earned promotion as an explanation for the global assembly pattern.

Then conclude:

> The result is therefore best interpreted as context-dependent island floral assembly with strong genus-level representation, not as one universal syndrome or an identified universal pollination mechanism.

## Results order

### Result 1 — response direction is context dependent

Keep H1/H2 first because it is the frozen primary family.

Report effect sizes and intervals for Palearctic and Tropical primary axes. Do not begin with the eight atomic outcomes.

### Result 2 — regional response composition is not merely an amplitude difference

Use the atomic fingerprint as a **descriptive decomposition of Result 1**, not a new test family.

Key examples:

- northern-midlatitude `actinomorphic_symmetry` is positive across both evidence scopes and both primary floristic strata;
- Tropical native-nonendemic `plain_colour` and `autonomous_selfing` are positive across both evidence scopes;
- no Tropical native-nonendemic structural-complexity contrast has a 95% interval excluding zero in both scopes.

Then report the eight-component vector angle: 95.27°, 100.21°, 116.49° and 99.53° across evidence-scope × stratum combinations.

Interpretation:

> regional differences involve re-weighting of trait components, not simply stronger versus weaker expression of one fixed syndrome vector.

### Result 3 — the main Palearctic response is concentrated at the genus transition

Lead with continuous attenuation rather than votes.

Across 16 frozen source-mode × evidence-scope × stratum profiles:

- family attenuation: 19.6–33.4%, median 26.2%;
- total genus attenuation: 78.8–85.9%, median 81.3%;
- conditional genus attenuation after family adjustment: 70.6–79.1%, median 75.7%.

Use `4/4 -> 4/4 -> 0/4` only as the robustness shorthand after the effect-size result.

Preferred wording:

> **Most attenuation occurs between the family-adjusted and genus-adjusted stages, localising the strongest Palearctic signal to a genus-level assembly transition rather than diffuse higher-taxonomic sorting.**

Do not call the percentages causal mediation.

### Result 4 — hierarchical expression itself is context dependent, but bounded

Keep the existing matched context × stage audit.

Headline:

> **The same plant-architecture response is altered by genus adjustment differently across northern-midlatitude and Tropical contexts in widespread native assemblages.**

Immediately state that the all-native result is weaker, which prevents a universal taxonomic-depth law.

### Result 5 — competing explanations fail promotion

Put four results in one compact falsification subsection:

1. area/capacity: 0/16 heteroskedastic-null promotions;
2. global nonlinear geometry: 0/12 step/reversal promotions after recovery calibration;
3. N1 channel heterogeneity: W=1.6187, df=3, p=0.65516, stopped before N2;
4. GloBI source breadth: 0/4 context × stratum cells promoted after source-prevalence, source-richness and reference-effort matching.

The GloBI predictor is especially useful as a transparency result because it was built outcome-blind from 24,577,183 raw interaction rows, retaining 714,785 explicit flower-interaction rows and 2,228 genus breadth estimates before any island outcome/distance join.

Preferred interpretation:

> **The genus-level assembly signal is strong, but it is not trivially recovered by island area, one universal threshold, coarse channel-specific pollinator isolation, or documented source-genus interaction breadth.**

### Result 6 — secondary pollination-associated architecture

Keep V4 after the main fingerprint. Use it to connect plant architecture to ecological interpretation, not to infer visitor identity.

## Figure architecture

### Figure 1 — conceptual hierarchy

Retain the v8 question:

> **Where, whether, and at what assembly level does an island floral response emerge?**

No temporal “when”. No direct causal arrow `distance -> Bombus loss -> floral trait`.

### Figure 2 — primary context branching

Effect-size plot for the frozen H1/H2 plant response. Palearctic and Tropical; all-analysis and direct-only side by side.

This remains the inferential anchor.

### Figure 3 — response-fingerprint atlas

This should become the main integrative figure.

Rows:

- northern-midlatitude;
- Tropical.

Columns grouped into three domains:

- colour;
- floral structure;
- reproductive assurance.

Within each domain show the atomic outcome estimates as points with 95% intervals, with all-analysis and direct-only visually paired.

Add two side panels:

**Panel B — cross-context orientation**

- display the 8-component vector angle for the four evidence/stratum combinations;
- label explicitly as descriptive geometry.

**Panel C — hierarchical attenuation**

- observed vector norm;
- after-family norm;
- after-genus norm;
- four source modes shown as thin repeated trajectories or uncertainty-like replicates;
- annotate family attenuation range (19.6–33.4%) and total genus attenuation range (78.8–85.9%).

The visual message should be: **what changes, and where does that response disappear in hierarchy?**

### Figure 4 — H3 biological assembly detail

Use genus entry/loading or the existing source-lineage representation output to show what genus-level assembly means biologically.

Do not duplicate Figure 3 attenuation.

### Figure 5 — explanation-elimination band

Use four compact blocks, not four large plots:

- area: 0/16 promoted;
- nonlinear geometry: 0/12 promoted;
- N1: p=0.655, stopped before N2;
- GloBI breadth: 0/4 promoted.

For GloBI include the outcome-blind pipeline mini-flow:

`24.6M raw -> 714,785 flower interactions -> 524,922 matched -> 2,228 genera -> effort-matched island test -> 0/4 promotion`

This figure demonstrates claim discipline rather than “many null results”.

## Discussion order

### 1. Response vectors, not one syndrome

Open with the primary context contrast and the near-orthogonal atomic fingerprints. The central point is that biogeographic context changes the mixture of floral/reproductive components, not simply the magnitude of one syndrome.

### 2. The hierarchy of assembly

Use the attenuation result to make H3 concrete. Family composition removes only a minority of the vector; genus composition removes most of the remainder. This makes “genus-level assembly” quantitatively interpretable.

### 3. Context changes hierarchical expression

Introduce the bounded context × stage interaction. Avoid “Tropical = evolution”; beyond-genus or stage-sensitive residuals can include within-genus sorting, unmeasured source structure and other filters.

### 4. Why simple mechanisms do not close the explanation

Integrate H4, geometry, N1 and GloBI breadth as deliberate attempts to explain the pattern.

The useful sentence is:

> **No single tested shortcut converted the hierarchical assembly result into a universal mechanism.**

This is stronger and cleaner than giving each null result an isolated paragraph.

### 5. Double-filter framework as the remaining hypothesis space

The plant/source side is demonstrated; the pollinator side remains biologically plausible but not identified. The failed N1 and GloBI breadth results should be used to justify why independent service data are needed, not to argue that pollinators are irrelevant.

### 6. Chapter 2 handoff

Chapter 1 establishes that the global pattern is context-dependent and hierarchy-dependent, while no universal response geometry or coarse pollinator mechanism is promoted.

Chapter 2 (`izu-core`) therefore asks a more identifiable within-system question:

> **Can directly observed changes in interaction/service produce trait-specific clines or thresholds within one biological system?**

This turns the global negative geometry result into a scale argument rather than a dead end.

## What not to add

Do not add new primary hypotheses to Chapter 1 after this point unless a genuine pre-outcome dataset exists.

In particular, do not:

- reinterpret vector attenuation as causal mediation;
- reinterpret the 8-component angle as a new confirmatory H2 test;
- treat GloBI missingness as specialization;
- rescue the GloBI result with another effort threshold or source mode;
- reopen N2;
- infer in-situ evolution from beyond-genus residuals;
- restore a global threshold headline after the calibrated geometry audit failed promotion.

## Submission implication

This maximizes the current data without inflating the claim ceiling. The EL-facing contribution is no longer simply “island floral syndromes vary among regions”. It is:

> **Ecological responses to the same geographic isolation gradient differ in component composition and in where they are expressed within assemblage hierarchy, while several intuitive mechanism shortcuts fail prospective or calibrated promotion tests.**

The island system supplies the global natural experiment; Chapter 2 supplies the direct functional mechanism test at local scale.
