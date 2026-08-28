# Chapter 1 V1 climate-overlap and contingency checkpoint

## Decision

V1 is complete for the Wave36 snapshot as a post-baseline prospective
falsification. The observed regional pattern remains part of the PR #141
historical baseline, but the stronger interpretation that H2 is a categorical
biogeographic contingency **not reducible to measured climate is not supported**.

- North--Tropical island climates have poor common support under the frozen
  criterion, so their raw distance-response vectors are not directly
  transportable over a shared area/climate population.
- After adding distance by climate-PC interactions, the residual categorical
  North--Tropical vector difference is unsupported in both trait-evidence scopes.
- Palearctic--Neotropical area/climate overlap is adequate, but its direct vector
  difference is unsupported after overlap weighting and after the
  distance-by-climate sensitivity.

H2 must therefore be stated as **region-associated, plausibly
environment-dependent branching**, not as a measured-climate-independent
biogeographic effect. V1 does not erase supported within-region slopes or identify
which environmental process is causal.

## Question

Does the regional difference in plant-side distance-response vectors persist (1)
among islands with overlapping area and climate distributions and (2) after the
distance effect is allowed to vary continuously with frozen climate PCs?

## Estimand

The target remains the joint distance-by-context coefficient vector for the
predeclared pollinator-name-free axes:

```text
accessibility/generalization + reproductive assurance.
```

V1 evaluates that vector under three modes:

1. the reference model with climate main effects;
2. outcome-blind area/climate overlap weights;
3. the reference model plus `distance × climate_pc1..4`, testing the residual
   `distance × context` vector.

The fixed response-specific support gates remain `<30` not promoted, `30--49`
pilot, and `>=50` confirmatory count support.

## Alternative explanation being attacked

The H2 baseline includes climate PCs as main effects, but regional differences in
distance slopes can still arise if isolation responses vary continuously with
climate or if the compared regions occupy different environmental support. V1
attacks the alternative that the categorical context label is standing in for
those measured environmental differences.

## Prospective implementation

The following choices were committed before the formal Wave36 run:

- contrasts: North--Tropical and Palearctic--Neotropical;
- balance variables: log island area and climate PCs 1--4;
- context propensity: binomial logit on the complete covariate island universe;
- overlap weights: `1-p(context_b)` in context B and `p(context_b)` in context A,
  normalized to mean 1 within context;
- no trait availability, trait state, score, effect estimate or p-value enters
  the weight model;
- poor overlap: either context ESS fraction `<0.25` or post-weight absolute SMD
  `>0.10`;
- continuous sensitivity: add all four `distance × climate PC` interactions;
- multiplicity: BH family across the two predeclared contrasts within validation
  mode, axis set, stratum and support tier;
- all-analysis and direct-only trait evidence are run separately.

The weighted test is reported even when positivity fails, but it cannot be used as
a transportable contrast in that case.

## Formal input and immutable output

- code HEAD: `c743f547a49ba5623bf1dccd624c7927f093254e`;
- trait run/artifact: `33137984367` /
  `wave36-trait-coverage-33137984367`;
- formal analysis run: `33142355523`;
- formal artifact ID: `9674488076`;
- artifact name: `chapter1-progressive-analysis-33142355523`;
- artifact digest:
  `sha256:6328f870ce9e5ffe4d11b4d3b40ed957da107fe3481d2027cff519da08d9dad1`.

Selected internal SHA-256 values:

| Artifact member | SHA-256 |
|---|---|
| `progressive_analysis_manifest.json` | `6f7efd62855d0660489be23f78d9d3a5fde942b7c6e56b1ee656f96f9b93ddda` |
| `climate-overlap/all/v1_overlap_summary.csv` | `8e72c6af7ff66510285799377ed264beac408f94edc8dcfb7dc58fdf464c6a33` |
| `climate-overlap/all/v1_vector_tests.csv` | `3bd7463d2e992c47bb958115295712646da361ca6cbe3318a8c3494fbabc59e2` |
| `climate-overlap/direct/v1_vector_tests.csv` | `3cac302fb10176b3b3e34621a2242b4bf738eb1a93ee76ae2e88d023b778f786` |

## Output

### Outcome-blind common support

| Comparison | Complete covariate islands A/B | ESS fraction A/B | Maximum post-weight abs SMD | Frozen overlap verdict |
|---|---:|---:|---:|---|
| North--Tropical | 3,115 / 2,993 | 0.071 / 0.101 | `<1.8e-13` | poor overlap |
| Palearctic--Neotropical | 947 / 1,045 | 0.351 / 0.274 | `<4.2e-15` | passes |

The nearly zero weighted SMDs show that the algorithm balanced the declared
predictors. They do not rescue North--Tropical positivity: achieving that balance
leaves only a very small effective fraction of each region.

### Confirmatory direct vector tests

Values are V1-family-adjusted q-values. `AN` is all-native and `NNE` is native
non-endemic.

| Evidence | Contrast | Reference AN / NNE | Overlap-weighted AN / NNE | Distance×climate-adjusted AN / NNE |
|---|---|---:|---:|---:|
| all-analysis | North--Tropical | 0.910 / **0.011** | 0.199 / 0.896 | 0.778 / 0.881 |
| all-analysis | Palearctic--Neotropical | **0.022** / 0.117 | 0.132 / 0.896 | 0.076 / 0.881 |
| direct-only | North--Tropical | 0.103 / 0.089 | 0.187 / 0.454 | 0.124 / 0.177 |
| direct-only | Palearctic--Neotropical | 0.394 / 0.353 | 0.159 / 0.454 | 0.124 / 0.177 |

Bold values pass the V1-adjusted `q <= 0.05` gate. The reference values are a
side-by-side V1 refit under the new two-contrast family; they do not replace the
historical PR #141 FDR classification.

## Verdict

1. **North--Tropical is not a common-support comparison.** Its environmental
   separation is too strong for a population-level contrast over shared measured
   area/climate support under the frozen ESS rule.
2. **The native-nonendemic all-analysis North--Tropical difference is absorbed in
   the continuous climate sensitivity.** Its q changes from 0.011 in the V1
   reference refit to 0.881 after distance-by-climate terms.
3. **Palearctic--Neotropical does have usable measured common support**, but the
   reference all-native contrast (q = 0.022) becomes unsupported after overlap
   weighting (q = 0.132) and after climate interactions (q = 0.076).
4. **Direct-only evidence does not independently recover either categorical
   contrast** under the V1 family or either added sensitivity.
5. The data therefore do not justify the phrase “biogeographic branching beyond
   measured climate.” They remain compatible with region-associated response
   heterogeneity whose environmental basis is unresolved.

## Effect on the H1--H5 tree

| Node | V1 consequence |
|---|---|
| H1 universal-syndrome rival | unchanged: the analysis still does not recover one universal syndrome |
| H2 biogeographic branching | narrowed: regional pattern exists descriptively, but measured-climate-independent categorical branching is not recovered |
| H3 source / lineage assembly | unchanged: V1 neither removes nor localizes lineage assembly |
| H4 area / capacity moderation | unchanged: area enters overlap balancing, but the biological area interaction awaits V3 |
| Secondary pollination concordance | unchanged: V1 tests the pollinator-name-free primary vector only |
| H5 channel-gated mechanism | not promoted: climate dependence does not identify a missing or retained pollination channel |

## Claim ceiling

V1 supports saying that Chapter 1 detects regional plant-side response patterns,
but those contrasts are not demonstrated to be independent of measured
environmental composition or climate-dependent isolation responses. It does not
support naming a regional pollinator mechanism, interpreting poor overlap as
biological absence, or claiming that climate is the unique cause of the observed
pattern.
