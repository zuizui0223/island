# Chapter 1 response-fingerprint atlas — 2026-09-14

## Purpose

The final Chapter 1 dataset is most informative when its outputs are not treated as separate H1–H5 analyses, but as orthogonal dimensions of one isolation-response fingerprint.

The atlas does not create a new hypothesis family. It synthesizes frozen estimates, the post-freeze descriptive effect-size fingerprint, and the separately frozen GloBI source-breadth extension.

## Fingerprint dimensions

### 1. Response direction

Primary H1/H2 asks whether the pollinator-name-free plant response points in the same direction across biogeographic contexts.

- Palearctic: accessibility/generalisation increases and reproductive assurance increases with source separation.
- Tropical direct-only: reproductive assurance increases while accessibility/generalisation decreases.

This is the primary evidence that one serial universal syndrome is inadequate.

### 2. Component identity and orientation

The eight frozen atomic contrasts are retained as descriptive component fingerprints:

- flower colour: `plain_colour`;
- floral structure: `generalized_form`, `actinomorphic_symmetry`, `shallow_open_tube`, `small_flower`;
- reproductive assurance: `self_compatibility`, `selfing_mating_system`, `autonomous_selfing`.

The strongest cross-evidence atomic fingerprints are not identical between contexts.

- Northern-midlatitude: `actinomorphic_symmetry` is positive in both evidence scopes in both primary floristic strata.
- Tropical native-nonendemic: `plain_colour` and `autonomous_selfing` are positive in both evidence scopes, while none of the four structural-complexity contrasts has a 95% interval excluding zero in either scope.

Using the same ordered eight-component vector in both contexts, the descriptive northern-midlatitude–Tropical angle is:

- all-analysis all-native: **95.27°**;
- all-analysis native-nonendemic: **100.21°**;
- direct-only all-native: **116.49°**;
- direct-only native-nonendemic: **99.53°**.

Thus the atomic response is approximately orthogonal to moderately opposed across contexts rather than merely the same syndrome expressed at different strength. This geometry is descriptive and does not replace the frozen multivariate H2 family.

### 3. Hierarchical/taxonomic expression

Wave52 is canonical. The Palearctic primary vector is retained after family adjustment but not after source-matched genus adjustment.

The effect-size synthesis quantifies this attenuation rather than reporting only `4/4 -> 4/4 -> 0/4`.

Across all 16 frozen source-mode × evidence-scope × stratum profiles:

- family-stage attenuation from the observed vector is **19.6–33.4%** (median **26.2%**);
- total attenuation after genus adjustment is **78.8–85.9%** (median **81.3%**);
- conditional attenuation added by the genus stage after family adjustment is **70.6–79.1% of the family-adjusted remainder** (median **75.7%**).

By evidence scope and stratum, total genus attenuation is:

- all-analysis all-native: 85.8–85.9%;
- all-analysis native-nonendemic: 82.4–82.5%;
- direct-only all-native: 78.8–79.1%;
- direct-only native-nonendemic: 79.9–80.2%.

The main attenuation therefore occurs at the family-to-genus transition rather than at the family stage. These are descriptive attenuation fractions, not causal mediation proportions.

A matched context × taxonomic-stage audit additionally yields a bounded hierarchical-depth result: the effect of genus adjustment on the same plant-architecture response differs between northern-midlatitude and Tropical contexts robustly in native non-endemics, but weakens in all-native assemblages.

### 4. Response geometry

Candidate shapes were frozen and recovery-audited before observed nonlinear results were opened.

- G2 step was identifiable in 11/12 cross-scope design cells;
- G3 hinge in 0/12;
- G4 reversal in 4/12.

Observed result: 12/12 cells are `monotonic_or_unresolved`; no cross-scope step or reversal is promoted.

Thus global island data do not support one common nonlinear transition even where a step would have been detectable.

### 5. Competing-explanation elimination

Four attractive explanatory shortcuts were evaluated rather than assumed.

- Area/capacity: 0/16 primary classifications pass the heteroskedastic-null mechanism-promotion gate. Apparent small-island amplification remains measurement-sensitive.
- Universal nonlinear threshold: 0/12 observed geometry cells promote step/reversal.
- Channel-specific pollinator isolation response: prospective N1 gives W=1.6187, df=3, p=0.65516 and stops before N2 without rescue.
- Source-side sampled interaction breadth: after source prevalence, source richness and broad GloBI reference-effort matching, 0/4 context × floristic-stratum cells pass the frozen promotion rule.

These are not interchangeable nulls. Together they narrow the explanation space around the stronger positive result: biogeographically contingent response direction and hierarchical assembly structure.

### 6. Source-side interaction dependency

The final `chapter1_globi_source_breadth_v2` result is negative under the prespecified promotion rule.

The plant-side predictor was constructed before any island outcome or distance was loaded from the pinned GloBI SUPPORTS archive:

- raw GloBI archive rows: **24,577,183**;
- retained explicit flower-interaction rows: **714,785**;
- plant-matched evidence rows: **524,922**;
- genus breadth rows: **2,228**.

No GloBI record was treated as specialization. The primary breadth metric was effective functional-channel number based on independent reference × channel evidence. Island expectations were then matched on source prevalence, source species richness and a broad GloBI independent-reference effort class.

Primary result:

- northern-midlatitude / all-native: not promoted;
- northern-midlatitude / native-nonendemic: not promoted;
- Tropical / all-native: not promoted;
- Tropical / native-nonendemic: not promoted.

Representative primary slopes are near zero. For `geo50_climate10`, northern-midlatitude all-native is `0.0060` [−0.0092, 0.0213], q=0.971, and Tropical all-native is `0.0077` [−0.0052, 0.0207], q=0.361.

Therefore the H3 genus-assembly pattern is not simply recovered by this independently documented source-genus interaction-breadth proxy. This does not establish that pollination dependence is irrelevant; it removes this particular sampled-breadth shortcut as a promoted explanation.

## Figure-ready atlas layout

Rows should be biogeographic contexts and columns should be the biological domains/components. Each context/domain cell can carry four visual encodings:

1. effect direction and magnitude;
2. evidence-scope agreement;
3. hierarchical attenuation from observed -> family -> genus;
4. geometry label (`monotonic_or_unresolved`, or a promoted nonlinear label if one had passed).

A compact side annotation can show the eight-component vector angle between northern-midlatitude and Tropical contexts.

A separate lower band should show the four explanation gates:

- area/capacity: not promoted;
- nonlinear step/reversal: not promoted;
- prospective N1 pollinator-channel heterogeneity: not promoted;
- effort-matched GloBI source breadth: not promoted.

The key visual message is therefore not simply that different regions have different slopes. It is:

> **Isolation-response fingerprints differ in which components move and where in assemblage hierarchy the response is expressed, while several simple mechanistic shortcuts fail explicit promotion gates.**

## Manuscript use

The atlas should not become six separate Results sections. The manuscript should use it as an integrative figure after the primary H1/H2/H3 sequence:

1. direction/context branching;
2. effect-size attenuation across taxonomic depth;
3. component fingerprint and bounded context × depth interaction;
4. one compact falsification band for area, geometry, N1 and GloBI breadth.

This preserves the primary preregistered inference while extracting substantially more information from the same frozen dataset.

## Claim ceiling

The atlas is a synthesis of already frozen evidence. It does not convert descriptive component intervals into a new confirmatory family, vector attenuation into causal mediation, vector angle into a replacement H2 test, beyond-genus residuals into evolution, or source-side interaction breadth into pollinator-loss mechanism.
