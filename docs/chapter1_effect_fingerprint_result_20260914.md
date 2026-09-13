# Chapter 1 effect-size fingerprint — 2026-09-14

## Decision

The frozen Chapter 1 results are more informative when expressed as effect-size fingerprints rather than only as support counts.

For the Palearctic primary two-axis response (`generalized_accessible` + `selfing_core`), source-matched genus adjustment reduces the Euclidean magnitude of the distance-response vector by **78.8–85.9%** across the four frozen source definitions, both evidence scopes and both primary floristic strata. Median attenuation across the 16 frozen profiles is **81.3%**.

This is a descriptive attenuation quantity, not a causal mediation percentage. It complements, rather than replaces, the frozen Wave52 classification `4/4 -> 4/4 -> 0/4`.

## Frozen input and receipt

- Chapter 1 workflow run: `34232450884`
- artifact: `chapter1-progressive-analysis-34232450884`
- artifact ID: `10058653212`
- artifact digest: `sha256:b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`
- effect-fingerprint workflow run: `34790185690`
- artifact: `chapter1-effect-fingerprint-34790185690`
- artifact ID: `10328176316`
- artifact digest: `sha256:b29d24623b39260d4296f8e9609af36def5cf880e44ff8394f11c428aa9abbcc`

The synthesis generated no new biological p-values and did not change any H1–H5 gate.

## Palearctic taxonomic attenuation

| evidence | stratum | source modes | genus attenuation |
|---|---|---|---:|
| all-analysis | all native | 4/4 | 85.8–85.9% |
| all-analysis | native non-endemic | 4/4 | 82.4–82.5% |
| direct-only | all native | 4/4 | 78.8–79.1% |
| direct-only | native non-endemic | 4/4 | 79.9–80.2% |

The useful manuscript wording is therefore:

> **Source-matched genus composition reduces the magnitude of the broad Palearctic floral/reproductive isolation-response vector by roughly four-fifths across source definitions and evidence scopes, after which the vector no longer passes the predeclared post-genus gate.**

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

## Figure use

The strongest use in the manuscript is an effect-size panel rather than another significance table:

1. observed -> family-adjusted -> genus-adjusted Palearctic two-axis vector attenuation;
2. colour / structure / reproductive-assurance component fingerprints for northern-midlatitude and Tropical contexts;
3. all-analysis and direct-only shown side by side rather than selecting the sharper evidence scope.

This gives a direct visual answer to two different questions:

- **how much of the primary response remains after source/taxonomic adjustment?**
- **which biological components actually move together within each context?**

## Claim boundary

The attenuation profile does not identify causal mediation, in-situ evolution, pollinator loss, effective service or the reason genera are differentially represented. The atomic fingerprint does not replace the preregistered H1/H2 multivariate tests.