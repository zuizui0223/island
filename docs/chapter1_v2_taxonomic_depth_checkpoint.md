# Chapter 1 V2 family-to-genus taxonomic-depth checkpoint

## Decision

V2 is complete for the frozen Wave36 snapshot.  In the broad two-axis Chapter 1
estimand, the Palearctic distance-response vector remains after both family and genus
source-composition nulls.  The pattern therefore cannot be summarized as family or
genus sorting alone.

This is an assemblage-depth result.  It does not identify within-genus evolution,
pollinator loss, effective pollination service, colonisation history, extinction, or
diversification.

## Frozen input and design

- contract freeze commit: `0c0629a28ccadda58ace2ca496de1623fa7bbb90`
- implementation commit: `602daeb33598d52c76e157a2fd9620ce243b3aec`
- portability and progressive-workflow commit:
  `97047994ddc4c9f7c30f8a1d0d7bea9e57cbcd18`
- formal workflow HEAD: `a69785d510b7eee80a47e92e773f6848e9a954b5`
- formal workflow run: `33145923109`
- formal artifact: `chapter1-progressive-analysis-33145923109`
- artifact ID: `9675896087`
- artifact digest:
  `sha256:4e31472005350e6e03f000340d044e470cc49a52e6e7ff83e5c1d4d5d93bd061`
- trait source run/artifact: `33137984367` /
  `wave36-trait-coverage-33137984367`
- taxonomy: 115,328 unique species, 115,145 with nonblank family
- taxonomy newline-canonicalized SHA-256:
  `c0586264d9a877b88c26e866264f2c77b20423cd93380a8eebe269582d53352e`

Family and genus labels were grouping variables only; they did not infer or fill a
trait.  Family and genus nulls used the same observed island species.  Source group
positions required at least two scored GIFT mainland species per group, and group
availability was built from the outcome-blind GIFT taxonomy match.

The two response axes were `generalized_accessible` and `selfing_core`.  The formal
stages were observed score, residual after the source-family expectation, and residual
after the source-genus expectation.  Both all-analysis and direct-only evidence, both
floristic strata, both context layers, and all four frozen source modes were retained.

## Data support

| evidence scope | scored source species | source-position rows | island-response decomposition rows |
|---|---:|---:|---:|
| all-analysis | 6,620 | 927 | 3,867 |
| direct-only | 3,505 | 872 | 3,054 |

The outcome-blind taxonomy-to-GIFT bridge contained 916,984 memberships representing
49,681 unique matched GIFT taxonomy species in each evidence scope.

## Result

Only the Palearctic vector passed the rule requiring support under all four source
definitions.  The same result occurred in both evidence scopes and both floristic
strata:

| scope | stratum | observed modes | after-family modes | after-genus modes | classification |
|---|---|---:|---:|---:|---|
| all-analysis | all native | 4/4 | 4/4 | 4/4 | within-genus or non-taxonomic residual retained |
| all-analysis | native non-endemic | 4/4 | 4/4 | 4/4 | within-genus or non-taxonomic residual retained |
| direct-only | all native | 4/4 | 4/4 | 4/4 | within-genus or non-taxonomic residual retained |
| direct-only | native non-endemic | 4/4 | 4/4 | 4/4 | within-genus or non-taxonomic residual retained |

In all-analysis data, both after-genus axes were individually supported under all
source modes.  In direct-only data, the after-genus vector remained supported, but the
axiswise result was driven by `generalized_accessible`; the positive `selfing_core`
residual did not pass the V2 source-mode family (`q = 0.058--0.082`).

The other 22 context-by-stratum classifications in each evidence scope were retained
as `observed_vector_not_robust_across_source_modes`.  V2 therefore does not assign a
taxonomic depth to northern-midlatitude, tropical, Neotropical, or the other formal
realm contexts.

## Reconciliation with the earlier genus result

The earlier exact-self-incompatibility genus null did not recover a robust
beyond-genus response.  V2 does not overwrite that result.  The estimands differ:

- the earlier result asks about the much smaller exact-SI floral subset;
- V2 asks about the broad, source-matched two-axis plant-response vector on common
  family/genus support.

The defensible combined statement is that the broad Palearctic response is not
absorbed by measured family or genus composition, whereas a beyond-genus claim is not
recovered for the exact-SI restricted subset.

## Portability audit

Formal attempt `33144735467` stopped before producing a V2 result because the tracked
taxonomy CSV had CRLF in the Windows worktree and LF in GitHub Actions.  The contract
was changed from a platform-specific raw-byte hash to an explicitly declared
newline-canonicalized hash.  No taxonomy rows, model choices, outcomes, or thresholds
changed.  The post-fix local classification tables were byte-for-byte identical after
CSV parsing to the pre-fix local tables.  The formal Linux classification tables also
matched the post-fix local tables exactly after CSV parsing in both evidence scopes.

## Effect on Chapter 1

V2 strengthens H3 only to this bounded level:

> The source-matched Palearctic floral/reproductive assemblage gradient includes a
> residual not explained by measured family or genus composition.

It does not establish why that residual exists.  Plant source-pool filtering remains
partly controlled through GIFT group availability, but the residual can still combine
unmeasured source composition, within-genus species sorting, other ecological filters,
and individual trait change.  Named-pollinator and effective-service mechanisms remain
outside the Chapter 1 claim ceiling.
