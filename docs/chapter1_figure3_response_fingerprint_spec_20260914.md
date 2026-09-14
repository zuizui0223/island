# Chapter 1 Figure 3 — response-fingerprint atlas specification

## Purpose

Figure 3 should be the paper's main integrative figure. It must show, without requiring the reader to follow analysis labels, that isolation-associated responses differ in **which components move** and **where in taxonomic hierarchy the response is expressed**.

It must not imply that descriptive atomic intervals are newly confirmatory, that attenuation is causal mediation, or that vector angle replaces the frozen H1/H2 heterogeneity test.

## Panel A — component response fingerprint

### Rows

Two focal biogeographic contexts:

1. northern-midlatitude;
2. Tropical.

### Columns / grouped domains

**Colour**
- `plain_colour`

**Floral structure**
- `generalized_form`
- `actinomorphic_symmetry`
- `shallow_open_tube`
- `small_flower`

**Reproductive assurance**
- `self_compatibility`
- `selfing_mating_system`
- `autonomous_selfing`

### Marks

For each component show the existing distance slope and 95% cluster-robust interval.

- all-analysis and direct-only should be paired horizontally or vertically;
- all-native and native-nonendemic can be separate subrows or use shape coding;
- zero line must be visually explicit;
- do not use stars or significance glyphs; interval crossing is enough for the descriptive panel.

### Required annotation

Use a narrow note rather than a legend headline:

> Atomic effects are descriptive decompositions of the frozen multivariate response; no new multiple-testing family was created.

### Biological callouts

Only annotate cross-evidence examples that are already clear:

- northern-midlatitude: `actinomorphic_symmetry ↑`;
- Tropical native-nonendemic: `plain_colour ↑`, `autonomous_selfing ↑`;
- Tropical structural simplification components: no cross-scope interval-excluding-zero pattern.

Do not annotate every point.

## Panel B — context orientation

Show one compact four-row display of the angle between the ordered eight-component northern-midlatitude and Tropical vectors.

| evidence | stratum | angle |
|---|---|---:|
| all-analysis | all-native | 95.27° |
| all-analysis | native-nonendemic | 100.21° |
| direct-only | all-native | 116.49° |
| direct-only | native-nonendemic | 99.53° |

Preferred graphic:

- two arrows from a common origin for each row, normalized to unit length;
- annotate the angle arc;
- alternatively use a four-row lollipop of angle degrees if arrow geometry becomes visually noisy.

Label:

> **Regional fingerprints differ in orientation, not only amplitude.**

Footnote:

> Vector angle is descriptive; direct H1/H2 tests remain the inferential basis for regional heterogeneity.

## Panel C — hierarchical attenuation

### Data

Palearctic primary two-axis response (`generalized_accessible` + `selfing_core`).

Stages:

1. observed;
2. after family expectation;
3. after source-matched genus expectation.

### Preferred graphic

A slopegraph / trajectory plot.

- x-axis: `observed -> family -> genus`;
- y-axis: Euclidean norm of the two-axis distance-response vector;
- 16 thin trajectories: 4 source modes × 2 evidence scopes × 2 primary floristic strata;
- optionally emphasize the median trajectory with a thicker line;
- do not average away the source-mode replicates entirely.

### Required quantitative labels

- family attenuation: **19.6–33.4%**, median **26.2%**;
- total genus attenuation: **78.8–85.9%**, median **81.3%**;
- conditional genus attenuation of family-adjusted remainder: **70.6–79.1%**, median **75.7%**.

Main callout:

> **Most attenuation occurs at the family-to-genus transition.**

Do not label this “80% explained by genus”.

## Panel D — explanation-elimination strip

Optional fourth panel if layout permits; otherwise move to Figure 5.

Four equal-width boxes:

1. **Area/capacity** — `0/16 promoted`
2. **Global nonlinear transition** — `0/12 promoted`
3. **Prospective pollinator-channel heterogeneity** — `W=1.619, p=0.655; N2 closed`
4. **Effort-matched GloBI source breadth** — `0/4 promoted`

The strip should visually read as a claim-ceiling audit, not as four failures.

Suggested heading:

> **Simple mechanism shortcuts do not close the assembly explanation.**

## Data sources / receipts

### Frozen Chapter 1 parent

- workflow run `34232450884`
- artifact `chapter1-progressive-analysis-34232450884`
- artifact ID `10058653212`
- digest `sha256:b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`

### Final effect-fingerprint synthesis

- workflow run `34796763611`
- artifact `chapter1-effect-fingerprint-34796763611`
- artifact ID `10330230959`
- digest `sha256:a2db683b6d46ca09d4eae6d9a3bfcaf85a4f172b5589227958bc388df8edba60`

Files:

- `atomic_response_fingerprint.csv`
- `atomic_cross_context_vector_geometry.csv`
- `taxonomic_vector_attenuation.csv`

### GloBI explanation audit

- workflow run `34790842133`
- artifact `chapter1-globi-source-breadth-v2-34790842133`
- artifact ID `10328696801`
- digest `sha256:475dda8b0eef4ce22ad1aa736fa82cc049042d801217c65d52641489288ca254`

## Caption skeleton

> **Figure 3 | Isolation-response fingerprints differ in component composition and hierarchical expression.** (A) Distance-associated atomic floral and reproductive contrasts in northern-midlatitude and Tropical island assemblages, shown as descriptive decompositions of the frozen multivariate response. (B) The ordered eight-component response vectors are approximately orthogonal to moderately opposed between contexts (95.3–116.5°), indicating differences in response composition rather than amplitude alone. (C) In the Palearctic primary response, family adjustment removes only part of the vector whereas source-matched genus adjustment attenuates 78.8–85.9% of the observed magnitude, with most attenuation occurring at the family-to-genus transition. (D) Prespecified or independently calibrated tests do not promote area/capacity, a universal nonlinear transition, channel-specific pollinator isolation, or source-genus sampled interaction breadth as simple explanations for the global assembly pattern. Atomic intervals, attenuation and vector angles are descriptive synthesis quantities and do not replace the frozen confirmatory H1–H3 tests.

## Design rule

The reader should understand the figure without knowing what H1, H2, H3, V2 or N1 mean. Those labels belong in Methods/Supplement, not in the primary visual narrative.
