# Chapter 1 V3 area-support artifact checkpoint

## Decision

V3 is complete for the frozen Wave36 snapshot.  The Palearctic lineage-entry
interaction keeps the predicted direction—distance effects are stronger on smaller
islands—but it does not pass the full predeclared measurement-support falsification.

H4 must therefore remain:

> **area is a measurement-sensitive modifier, not established evidence for island
> capacity, founder filtering, or pollinator-population persistence.**

## Frozen input and design

- contract freeze commit: `718d38a7c1f015daa961b07ad01e8357f0751833`
- implementation commit: `97047994ddc4c9f7c30f8a1d0d7bea9e57cbcd18`
- minimum-response-support correction:
  `a69785d510b7eee80a47e92e773f6848e9a954b5`
- formal workflow run: `33145923109`
- formal artifact: `chapter1-progressive-analysis-33145923109`
- artifact ID: `9675896087`
- artifact digest:
  `sha256:4e31472005350e6e03f000340d044e470cc49a52e6e7ff83e5c1d4d5d93bd061`
- trait source run/artifact: `33137984367` /
  `wave36-trait-coverage-33137984367`

The primary target was Palearctic GIFT-source-matched genus `entry_enrichment` in the
all-native and native-nonendemic strata.  The broad all-analysis and newly constructed
broad direct-only lineage scopes used the same GIFT assignments,
prevalence-by-source-richness matching, and minimum five represented genera.

Three frozen sensitivity modes were fitted:

1. equal island weight;
2. `min(n_response_species, 50)` weights normalized to mean one within response;
3. equal-weight refits on an outcome-blind common area range shared by low- and
   high-support islands.

The common-support split used the median minimum response support in each analysis
cell, required at least 30 islands per half before trimming, and intersected each
half's area 5th--95th percentile range.  Response values never defined the split,
weights, or overlap.

The frozen 1,000-draw simulation fitted a no-interaction mean model and assigned
independent normal error variance proportional to mean support divided by island
support.  Species support was never inserted into the biological model as a causal
control.

## Support diagnostic

Palearctic lineage response support increased substantially with island area:

- all-analysis Spearman rho: `0.520--0.561` across strata/source modes;
- direct-only Spearman rho: `0.494--0.541`;
- equal-weight primary samples: 140 all-analysis and 137 direct-only islands;
- common-support samples: 87/93 all-analysis islands and 78/90 direct-only islands
  for all-native/native-nonendemic strata.

This is precisely the measurement-precision structure V3 was designed to attack.  It
does not show that richness is a nuisance variable: richness/support can itself be part
of island assembly.

## Primary result

Across every source mode and both evidence scopes, the point estimates retained a
positive distance slope at mean area and a negative distance-by-area coefficient.
Direct-only estimates therefore did not reverse the predicted small-island
amplification direction.

Formal support was not stable across all frozen gates:

- all three sensitivity modes formally supported the all-native result for `geo_k10`
  and `geo_k20`;
- common-support or capped-information axis support was lost in the remaining
  source-mode/stratum combinations;
- the heteroskedastic-null exceedance probability was `0.166--0.302` in all-analysis
  data and `0.190--0.325` in direct-only data, never at or below the frozen `0.05`
  threshold.

Consequently all 16 primary stratum-by-tier-by-source-mode classifications were
`retain_area_as_measurement_sensitive_modifier_only`.  No result was promoted to
`area_capacity_interpretation_strengthened`.

The formal Linux and local rerun classification tables matched exactly after CSV
parsing.  Their 16 null-exceedance probabilities and all pass/fail decisions also
matched; only non-headline floating-point text differed.

## Secondary plant-response branches

The two-axis plant family (`accessibility_generalization` and
`reproductive_assurance`) also failed to show one stable small-island amplification
across equal, capped-information, direct-only, and common-support modes.  A few
mode-specific Palearctic axes were supported, but capped-information fits supplied no
formally supported small-island plant-branch axis in either evidence scope.  These
isolated sensitivities are retained as non-robust rather than selected post hoc.

## Effect on Chapter 1

The baseline area interaction remains a useful descriptive clue because its direction
is consistent across source definitions and direct-only evidence.  V3 prevents it from
being used as proof of the proposed capacity pathway:

```text
island area -> pollinator establishment/persistence -> effective service
```

Chapter 1 does not observe the arrows in that pathway.  It also does not observe a
genetic founder effect.  Area can represent target size, habitat/resource diversity,
population persistence, sampling support, or several of these together.  Separating
those explanations requires independent service/population data or a design with
trait-measurement precision decoupled from island area.
