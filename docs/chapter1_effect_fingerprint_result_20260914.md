# Chapter 1 effect-size fingerprint — 2026-09-14

## Decision

The frozen Chapter 1 results are more informative when expressed as effect-size fingerprints rather than only as support counts.

For the Palearctic primary two-axis response (`generalized_accessible` + `selfing_core`), source-matched genus adjustment reduces the Euclidean magnitude of the distance-response vector by **78.8–85.9%** across the four frozen source definitions, both evidence scopes and both primary floristic strata. Median attenuation across the 16 frozen profiles is **81.3%**.

The hierarchy is not simply “taxonomy matters”. Family adjustment alone reduces the vector by **19.6–33.4%** (median **26.2%**), whereas adding genus adjustment removes **70.6–79.1% of the family-adjusted remainder** (median **75.7%**). The dominant attenuation therefore occurs at the family-to-genus transition.

These are descriptive attenuation quantities, not causal mediation percentages. They complement, rather than replace, the frozen Wave52 classification `4/4 -> 4/4 -> 0/4`.

## Frozen input and receipt

- Chapter 1 workflow run: `34232450884`
- Chapter 1 artifact: `chapter1-progressive-analysis-34232450884`
- Chapter 1 artifact ID: `10058653212`
- Chapter 1 artifact digest: `sha256:b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`
- final effect-fingerprint workflow run: `34796763611`
- effect-fingerprint artifact: `chapter1-effect-fingerprint-34796763611`
- artifact ID: `10330230959`
- artifact digest: `sha256:a2db683b6d46ca09d4eae6d9a3bfcaf85a4f172b5589227958bc388df8edba60`

The synthesis generated no new biological p-values and did not change any H1–H5 gate.

## Palearctic taxonomic attenuation

| evidence | stratum | source modes | family attenuation | total genus attenuation |
|---|---|---|---:|---:|
| all-analysis | all native | 4/4 | 32.6–33.4% | 85.8–85.9% |
| all-analysis | native non-endemic | 4/4 | 19.6–20.2% | 82.4–82.5% |
| direct-only | all native | 4/4 | 26.6–28.8% | 78.8–79.1% |
| direct-only | native non-endemic | 4/4 | 25.6–25.8% | 79.9–80.2% |

The useful manuscript wording is therefore:

> **Source-matched genus composition reduces the magnitude of the broad Palearctic floral/reproductive isolation-response vector by roughly four-fifths across source definitions and evidence scopes, with most of that attenuation occurring between the family-adjusted and genus-adjusted stages.**

Do not write that genus composition “explains 80% causally” or that the residual is exactly zero.

## Component fingerprint

The frozen eight atomic response contrasts show why the assemblage response should not be described as one internally coherent simplification score.

### Northern mid-latitude

Across all-analysis all-native assemblages, `actinomorphic_symmetry` and `self_compatibility` have positive 95% intervals, whereas the other six atomic contrasts are uncertain. In native non-endemics, `plain_colour` also has a positive all-analysis interval. Direct-only evidence retains a positive `actinomorphic_symmetry` interval and, in all natives, a positive `selfing_mating_system` interval.

### Tropical

The strongest cross-scope native-nonendemic component signals are different from the northern set:

- `plain_colour` is positive in both all-analysis and direct-only evidence;
- `autonomous_selfing` is positive in both all-analysis and direct-only evidence;
- none of the four structural-complexity contrasts has a 95% interval excluding zero in either evidence scope.

Thus the broad tropical signal is not well summarized as a coordinated shift toward floral structural simplification. Colour and reproductive components can move while the individual structural contrasts remain weak or inconsistent.

These intervals are descriptive re-expressions of existing frozen estimates. They are not a new multiple-testing family and must not be promoted as newly confirmatory atomic discoveries.

## Cross-context vector orientation

Using the same ordered eight-component vector in northern-midlatitude and Tropical contexts, the descriptive vector angles are:

| evidence | stratum | angle |
|---|---|---:|
| all-analysis | all native | **95.27°** |
| all-analysis | native non-endemic | **100.21°** |
| direct-only | all native | **116.49°** |
| direct-only | native non-endemic | **99.53°** |

Cosine similarity is negative in all four comparisons. The component fingerprints are therefore approximately orthogonal to moderately opposed rather than the same vector expressed at different amplitude.

This is descriptive geometry only. The preregistered H1/H2 direct multivariate tests remain the inferential basis for biogeographic response heterogeneity.

## Figure use

The strongest use in the manuscript is an effect-size panel rather than another significance table:

1. observed -> family-adjusted -> genus-adjusted Palearctic two-axis vector attenuation;
2. colour / structure / reproductive-assurance component fingerprints for northern-midlatitude and Tropical contexts;
3. all-analysis and direct-only shown side by side rather than selecting the sharper evidence scope;
4. a compact annotation for the eight-component cross-context vector angle.

This gives direct visual answers to three different questions:

- **how much of the primary response remains after source/taxonomic adjustment?**
- **where in the family-to-genus hierarchy does the main attenuation occur?**
- **are northern and Tropical responses the same syndrome at different strength, or differently oriented component mixtures?**

## Claim boundary

The attenuation profile does not identify causal mediation, in-situ evolution, pollinator loss, effective service or the reason genera are differentially represented. The atomic fingerprint and vector angle do not replace the preregistered H1/H2 multivariate tests.
