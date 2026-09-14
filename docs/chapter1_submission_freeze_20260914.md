# Chapter 1 submission freeze — 2026-09-14

This file supersedes the 2026-09-09 narrative freeze. Earlier manuscript/framing files remain preserved in git history and may be consulted as provenance, but they are no longer the canonical submission surface.

## Canonical manuscript

The single canonical working manuscript is:

- `docs/chapter1_manuscript_full_v8_hierarchical_syndrome_20260914.md`

Canonical inference/figure framing:

- `docs/chapter1_court_evidence_and_theory_synthesis_20260914.md`
- `docs/chapter1_figure1_hierarchical_syndrome_spec_20260914.md`
- `docs/chapter1_literature_positioning_20260909.md`

## Frozen primary data and analysis

Trait snapshot:

- 106,295 analysis-applicable angiosperm species;
- 318,885 species-by-raw-axis denominator;
- 222,688 resolved cells = 69.83%;
- raw axes: flower colour, floral structural complexity, reproductive assurance.

Primary Chapter 1 execution:

- run `34232450884`;
- artifact `chapter1-progressive-analysis-34232450884`;
- artifact ID `10058653212`;
- digest `sha256:b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`.

No v8 narrative change alters the frozen H1-H4 estimates from that execution.

## Frozen central results

1. **Universal syndrome rejected as a global rule.** Floral/reproductive response vectors differ among biogeographic contexts.
2. **Component decoupling.** Palearctic accessibility/generalization and reproductive assurance increase with separation; tropical reproductive assurance can increase while accessibility/generalization decreases.
3. **Taxonomic depth is the strongest explanatory result.** Broad Palearctic support follows `observed 4/4 -> family 4/4 -> genus 0/4`. Genus adjustment attenuates approximately 78.8–85.9% of the observed vector; conditional family-to-genus attenuation is approximately 70.6–79.1%.
4. **Area is not promoted as a mechanism.** H4 promotion is 0/16 under the frozen safeguards.
5. **No common global nonlinear breakpoint.** After calibrated geometry qualification, observed broad atomic responses yield 0/12 promoted nonlinear shapes.

## Frozen robustness / adverse-evidence results

### V5 trait missingness

The broad Palearctic primary response survives the finite predeclared MNAR trait-resolution grid. Some reproductive details and context contrasts remain more sensitivity-bounded. Arbitrary MNAR is not ruled out.

### V6 species-list incompleteness

Canonical V6 run:

- run `34800498716`;
- artifact ID `10331282464`.

Interpretation:

- Palearctic accessibility survives 99/100 baseline-supported bias surfaces;
- all 80/80 remote-under-survey scenarios in the biologically concerning direction retain support;
- tropical accessibility is more sensitive: 40/75 baseline-supported surfaces tip;
- North–Tropical vector heterogeneity is more robust than the tropical single-axis result.

V6 does not estimate true island-flora completeness and does not protect against arbitrary taxon-dependent omission.

### N1 independent channel heterogeneity

- `W=1.6187`, df=3, p=0.65516;
- not promoted;
- limited power for modest channel heterogeneity means this is non-identification, not proof of equality.

### H5c independent biotic-vs-wind specificity

Qualification:

- run `34803493307`;
- 5,771 unambiguous independent `pollen_vector_mode` species;
- 1/8 cells fully qualified.

Observed opening:

- run `34803837463`;
- only direct-only x native-nonendemic x Palearctic opened;
- distance x biotic interaction = +0.06495;
- 95% CI [-0.09030, 0.22020];
- p=0.41221;
- classification: `no_pollination_mode_specificity_support`.

This result prevents promotion of independent pollination-mode specificity; it does not show that pollinators are irrelevant.

### H5d distributed-threshold identifiability

- run `34803545574`;
- 0/8 design cells qualified;
- classification accuracy approximately 0.733–0.778;
- false distributed-threshold selection under heterogeneous smooth clines 0.19–0.255;
- observed genus threshold distributions remain closed.

The threshold concept may be used as a scale bridge to a local system but not as an identified Chapter 1 mechanism.

## Canonical publication-facing claim

> **Geographic isolation does not impose one floral island syndrome. It reorganizes floral and reproductive composition through biogeographically contingent, hierarchically structured assemblage filtering; the strongest Palearctic response is concentrated at the family-to-genus transition, while several simple adaptive and interaction-based explanations fail explicit falsification tests.**

General conceptual statement:

> **Trait syndromes observed across environmental gradients need not be repeated organismal adaptations; they can emerge from hierarchical changes in community membership, and identifying the assembly depth of a syndrome is therefore a prerequisite for mechanistic interpretation.**

## Claim ceiling

The submission may claim:

- context-dependent floral/reproductive response vectors;
- component decoupling;
- strong Palearctic family-to-genus attenuation;
- robustness of the Palearctic core to the specified V5/V6 missingness challenges;
- failure to promote continuous area, a common global breakpoint, coarse channel heterogeneity, sampled source breadth, independent biotic-vs-wind specificity, or distributed-threshold generators.

The submission must not claim:

- pollinator loss caused the Palearctic response;
- H5c proves pollinators do not matter;
- genus attenuation proves dispersal alone or proves absence of within-lineage evolution;
- tropical accessibility is as observation-robust as the Palearctic branch;
- endemicity is a temporal axis or establishes in-situ evolution;
- smooth global responses imply that local ecological systems lack thresholds;
- observed genus-specific thresholds are identified from the current Chapter 1 design.

## Chapter 1 / Chapter 2 boundary

Chapter 1 establishes **where, which components, and at what assembly depth** the global syndrome is expressed.

Chapter 2 / `izu-core` is reserved for **how and why functionally** within one resolved biological system:

`interaction state -> effective pollination service -> reproductive outcome -> phenotype`.

The local system may test cline, threshold, shared-breakpoint, and channel-specific response geometries prospectively. It must not be used retrospectively to relabel the Chapter 1 H3 result as pollinator-caused.
