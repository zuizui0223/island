# Chapter 1 V4 source-trained architecture decomposition checkpoint

## Status

V4 executed prospectively on the frozen Wave36 Chapter 1 snapshot and is retained as
a **secondary plant-trait architecture decomposition**.  It does not identify a
pollinator guild or effective pollination service.

- contract freeze commit: `e5fc228a6e933143842dea0ec2c6a1dc9dc63182`
- implementation commit: `cbcab569a4d8df86f288689f9224ec919588449c`
- formal workflow run: `33143552969`
- formal artifact: `chapter1-progressive-analysis-33143552969`
- artifact id: `9675018453`
- artifact digest:
  `sha256:ca0aafb5653f6b8feb8190d7aa34d9f18ce651a589fbc2402d20f1d9e93113f7`
- trait source run: `33137984367`
- trait source artifact: `wave36-trait-coverage-33137984367`

The first workflow dispatch accidentally used an obsolete trait-run id and was
cancelled before use.  Run `33143552969`, which names the trait run and artifact in the
frozen snapshot manifest, is the only admitted V4 formal result.

## Frozen estimand

The three sampled template scores (`large_bee_like`, `butterfly_like`, `bird_like`)
were decomposed into:

1. a shared architecture factor trained only on complete-case species matched to the
   frozen GIFT mainland native flora; and
2. one rank-one-reconstruction residual for each sampled template.

The source means, population standard deviations, and PCA1 loadings were fitted
separately in the all-analysis and direct-only evidence scopes.  Island status,
distance, area, climate, context, realm, and island outcomes were absent from the
factor fit.  Missing template scores were not imputed.

The four architecture components were then evaluated under the existing Chapter 1
distance, area, climate, spatial-cluster, stratum, support-tier, and FDR contract.
The same components were carried through all four frozen GIFT source-assignment modes
and a source-matched genus entry/loading/beyond-genus decomposition.  Source-mode
multiplicity is included in the V4 result family.

## Source-trained common factor

| evidence scope | complete scored species | GIFT source species | PCA1 variance | large-bee loading | butterfly loading | bird loading |
|---|---:|---:|---:|---:|---:|---:|
| all-analysis | 1,492 | 997 | 0.85294 | 0.54343 | 0.58532 | 0.60174 |
| direct-only | 1,049 | 829 | 0.84512 | 0.53816 | 0.58853 | 0.60334 |

All three source-factor correlations were positive and high (all-analysis
`0.869/0.936/0.963`; direct-only `0.857/0.937/0.961`).  The source factor therefore
captures a large shared component of the three sampled syndrome templates.  Raw
template slopes must not be read as three independent guild responses.

## Observed and source-adjusted WHEN/WHERE result

The observed all-analysis four-axis vector was supported in tropical islands and in
the all-native Palearctic and Neotropical strata, but the focal between-context vector
tests were not supported.  Direct-only evidence was sharper: both the
northern-midlatitude versus tropical and Palearctic versus Neotropical vector
differences were supported in all-native and native-nonendemic strata.

After GIFT source adjustment and correction across the four source modes, the following
confirmatory axis responses survived **all four** source definitions:

- tropical shared architecture increased with distance in both strata and both
  evidence scopes;
- in direct-only northern-midlatitude floras, the large-bee residual increased and the
  butterfly residual decreased with distance in both strata;
- in direct-only Palearctic floras, the bird residual increased with distance in both
  strata.

The source-adjusted northern-midlatitude versus tropical vector difference survived
all four source modes in both evidence scopes and both strata.  The Palearctic versus
Neotropical difference survived all four modes in direct-only evidence, but not in the
all-analysis scope.

These residual directions are contrasts among overlapping plant-trait templates.  For
example, `large_bee_like_residual > 0` means only that the large-bee template component
is high relative to the source-trained shared factor; it is not an observation of
large bees and is not evidence for Bombus gain, loss, or effectiveness.

## Lineage-depth result within V4

The direct-only genus decomposition retained source-mode-robust confirmatory
`loading_increment` responses, especially in northern-midlatitude and Palearctic
floras.  These include opposing template-residual changes and therefore show that the
architecture pattern is not represented by one shared common factor alone.

Tropical large-bee and butterfly residuals also showed beyond-genus responses in three
geographic source modes, but not uniformly in the climate-nearest mode.  No
beyond-genus component survived all four source definitions.  All-analysis genus
results were less stable and no lineage outcome/axis combination survived every source
mode.

Therefore V4 supplies a **candidate architecture residual** for later H5 mechanism
testing, but it does not promote H5 to supported.  The present ceiling remains:

> sampled syndrome residuals are plant-trait contrasts; functional replacement and
> realized effective service remain not evaluable in Chapter 1.

## Family gate and V2 handoff

V4 records eight family-component combinations as
`not_evaluable_pending_V2_taxonomy_contract`.  This does not mean that the repository
contains no family table.  It means that V4 did not declare and hash a comprehensive
species-to-family mapping as a formal progressive-workflow input before executing the
factor analysis.

The repository master contains 115,328 unique species, of which 115,145 currently have
a nonblank family.  V2 must first promote that table to a frozen, validated formal input
and then compare family composition, genus composition, and within-genus species
loading.  Family is a grouping variable only and must never be used to infer a missing
trait.

## Cross-platform audit

The formal Linux artifact and the independent Windows rerun had identical categorical
headline results.  Numerical differences were normally below `3.2e-12`.  One
non-headline pilot Palearctic-Oceanian source-adjusted comparison used covariance rank
3 on Linux and rank 4 on Windows because the matrix was nearly singular.  The formal
Linux artifact is authoritative; this pilot-only numerical portability issue does not
alter any focal confirmatory V4 conclusion and is retained for later numerical QA.

## Chapter 1 interpretation after V4

V4 strengthens three bounded statements:

1. the sampled guild templates mainly share one floral-architecture dimension;
2. the regional distance response is nevertheless not wholly reducible to that common
   dimension or to measured GIFT source composition; and
3. part of the direct-only pattern is organized through genus entry/species loading,
   while a smaller tropical residual can remain beyond genus under some source modes.

It does **not** establish that Bombus loss generated the Palearctic pattern, that birds
or butterflies replaced a lost service, or that any named pollinator was present,
absent, mobile, or effective.  Those remain discussion compatibilities and Chapter 2
mechanism targets.
