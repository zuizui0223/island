# Chapter 1 Wave52 progressive reanalysis checkpoint

## Decision

Wave52 was reanalysed prospectively under the unchanged
`chapter1_progressive_analysis_v1` contract. The main Chapter 1 conclusion changes in
one important place: the broad Palearctic two-axis plant-response vector is no longer
retained after the source-matched genus expectation is removed.

The updated bounded conclusion is:

> A universal island floral syndrome is still not recovered. A broad Palearctic
> floral/reproductive distance gradient is retained and survives family-level source
> composition, but the primary two-axis response is now compatible with genus-level
> lineage assembly rather than a robust beyond-genus residual. Its measured-climate
> independence and area/capacity mechanism remain unestablished, and H5 remains not
> evaluable without independent pollination-channel data.

This is a biological conclusion change, not a CI failure.

## Formal inputs and run

- fixed scientific contract: `chapter1_progressive_analysis_v1`
- Wave52 trait run: `33470056127`
- Wave52 trait artifact: `wave52-support-two-reproductive-unlock-33470056127`
- Wave52 trait artifact ID: `9786227115`
- Wave52 trait artifact digest:
  `sha256:a6dcf159e79a678a35d643e6f72e15afffe979b801d411005d28c286bf079c89`
- previous lossless comparison snapshot: Wave51 run `33467582475`
- reanalysis run: `33587356935`
- reanalysis artifact: `chapter1-progressive-wave52-33587356935`
- reanalysis artifact ID: `9830619059`
- reanalysis artifact digest:
  `sha256:9720586dafc6a46f1a4b1a2c9ad3e28ee5dd12da91e66a4c9c52fe718d5e9ee0`

The one-shot replay executed the `analyse` run blocks from
`.github/workflows/run-chapter1-progressive-trait-analysis.yml` in their tracked order;
model definitions, non-trait inputs, support gates, FDR families and claim ceilings
were not changed.

## Snapshot transition audit

The fail-closed Wave51 -> Wave52 transition initially stopped on one audited value
revision. The transition was inspected before rerunning with explicit revision
permission:

- unchanged: 318,622 species-axis cells;
- newly resolved: 262;
- value revisions: 1;
- quality downgrades: 0;
- became unresolved: 0.

The single revision is `Helichrysum italicum` reproductive assurance. Existing High
`mating_system=["predominantly_outcrossing"]` was retained and Medium direct
`self_incompatibility=["SC"]` was added from an explicit sometimes-autogamous species
account. The acquisition review explicitly prohibited conversion to autonomous
selfing, so this is an evidence enrichment rather than a contradictory replacement.

## Analysis-usable coverage

Wave52 analysis-usable coverage is 184,917 / 318,885 species-axis cells = 57.99%.
Global fill remains descriptive only.

| axis | Wave36 usable | Wave52 usable | Wave52 resolved species |
|---|---:|---:|---:|
| flower colour | 52.59% | 52.80% | 56,122 |
| floral structural complexity | 86.03% | 86.19% | 91,613 |
| reproductive assurance | 32.81% | 34.98% | 37,182 |
| all cells | 57.14% | 57.99% | 184,917 cells |

The largest material increase is therefore reproductive assurance, which is exactly
where the progressive contract expected later evidence to be most consequential.

## H1-H5 transition

| node | Wave52 transition | decision |
|---|---|---|
| H1 universal-syndrome rival | stable | one universal island syndrome is still not recovered |
| H2 biogeographic/environmental branching | stable with evidence-scope sharpening | regional response heterogeneity remains, but measured-climate-independent categorical branching is still not established |
| H3 source/lineage assembly | **materially changed** | primary Palearctic vector now survives family adjustment but not genus adjustment; genus-level assembly becomes the strongest supported depth |
| H4 area/capacity moderation | stable | area remains a measurement-sensitive modifier, not established capacity/founder/pollinator-persistence evidence |
| H5 channel-gated residual mechanism | stable, not evaluable | no independent channel exposure, disruption, visitation or effective-service data enter Chapter 1 |

## V5: stable MNAR-bounded result

The Palearctic within-context universal plant-response vector remains baseline-supported
throughout every finite predeclared MNAR scenario in both all-analysis and direct-only
evidence, for all-native and native-nonendemic strata.

The lower-coverage reproductive branch still cannot be declared immune to arbitrary
MNAR. V5 therefore remains `stable`: broad Palearctic pattern retained, reproductive
detail bounded, H5 not promoted.

## V1: conclusion stable, direct-only evidence sharper

North--Tropical still fails the outcome-blind common-support positivity gate. Its
climate-adjusted categorical contrast therefore remains non-transportable under the
frozen criterion.

Palearctic--Neotropical has adequate common support. Wave52 is less one-sided than
Wave36 because some direct-only tests now pass:

| scope / stratum | overlap-weighted q | distance x climate q |
|---|---:|---:|
| all-analysis / all-native | 0.309 | 0.814 |
| all-analysis / native-nonendemic | 0.981 | 0.593 |
| direct-only / all-native | **0.0117** | **0.0141** |
| direct-only / native-nonendemic | 0.139 | **0.0403** |

Because the result does not replicate across both evidence scopes and both strata, the
Chapter 1 claim does not advance to measured-climate-independent biogeographic
branching. V1 is therefore `stable` at the claim level, with a direct-only precision
strengthening that should be retained rather than selected as the headline.

## V4: secondary architecture result strengthens

The source-trained shared factor now explains:

- 86.85% of all-analysis variance across the three sampled templates; and
- 86.44% of direct-only variance.

Wave36 values were 85.29% and 84.51%, respectively, so the conclusion that the named
guild templates mostly share one plant-architecture dimension is strengthened.

A secondary result does change. Wave36 had no beyond-genus architecture component that
survived all four source definitions. Wave52 now contains source-mode-robust
beyond-genus template residuals. The cleanest cross-evidence result is tropical
`large_bee_like_residual > 0`, which survives all four source modes in both all-analysis
and direct-only evidence and in both all-native and native-nonendemic strata.

This does **not** restore a Palearctic beyond-genus primary claim. V4 residuals are
contrasts among overlapping sampled plant-trait templates; they do not identify large
bees, butterflies, birds, loss, replacement or effective service.

## V2: material conclusion reversal at taxonomic depth

This is the decisive Wave52 change.

Wave36 Palearctic primary two-axis classification was:

| scope | stratum | observed | after family | after genus |
|---|---|---:|---:|---:|
| all-analysis | all-native | 4/4 | 4/4 | 4/4 |
| all-analysis | native-nonendemic | 4/4 | 4/4 | 4/4 |
| direct-only | all-native | 4/4 | 4/4 | 4/4 |
| direct-only | native-nonendemic | 4/4 | 4/4 | 4/4 |

Wave52 is:

| scope | stratum | observed | after family | after genus | classification |
|---|---|---:|---:|---:|---|
| all-analysis | all-native | 4/4 | 4/4 | **0/4** | compatible with genus-level assembly beyond family |
| all-analysis | native-nonendemic | 4/4 | 4/4 | **0/4** | compatible with genus-level assembly beyond family |
| direct-only | all-native | 4/4 | 4/4 | **0/4** | compatible with genus-level assembly beyond family |
| direct-only | native-nonendemic | 4/4 | 4/4 | **0/4** | compatible with genus-level assembly beyond family |

At the same time, the number of scored GIFT source species used by V2 increased from
6,620 to 7,697 in all-analysis and from 3,505 to 4,549 in direct-only evidence. The
changed depth classification is therefore exactly the kind of prospective update the
progressive workflow was designed to expose.

The old sentence

> the broad Palearctic response is not absorbed by measured family or genus
> composition

is no longer supported by the current snapshot and must not remain the current
headline.

The replacement is:

> the broad Palearctic response persists beyond measured family composition but is
> compatible with source-matched genus composition.

This strengthens H3 as a lineage-assembly explanation and weakens any interpretation
requiring a robust within-genus/non-taxonomic residual for the primary response.

## V3: stable precision-sensitive area result

All 16 frozen primary classifications remain
`retain_area_as_measurement_sensitive_modifier_only`.

- heteroskedastic-null pass: 0 / 16;
- direct-only contradiction: 0 / 16;
- 10 / 16 cells formally support small-island amplification through all three fitted
  sensitivity modes, but the frozen heteroskedastic-null gate prevents mechanism
  promotion.

Thus H4 does not change.

## Secondary sampled-guild concordance

The broad architecture contrast remains, but exact supported templates changed with
new evidence.

Wave52 confirmatory results:

- Tropical: bird-like, butterfly-like and large-bee-like concordance all increase with
  distance in both evidence scopes and both floristic strata.
- Neotropical: none of the three sampled guild axes is supported.
- Palearctic: butterfly-like and large-bee-like concordance decline; bird-like decline
  is no longer FDR-supported in either evidence scope.
- Northern midlatitude: butterfly-like decline is supported in both evidence scopes;
  large-bee-like decline is supported in direct-only but not all-analysis evidence;
  bird-like remains unsupported.

Therefore the ecological statement remains one of region-specific floral architecture,
not pollinator identity. Exact named-template membership is trait-snapshot-sensitive
and should remain secondary.

## Current Chapter 1 ceiling after Wave52

The strongest current statement is:

> Island isolation is associated with region-specific floral/reproductive assemblage
> structure rather than one universal syndrome. The broad Palearctic gradient is
> robust to finite reproductive-trait MNAR sensitivity and persists beyond measured
> family composition, but the primary two-axis response is now compatible with
> source-matched genus-level lineage assembly. Some secondary tropical floral-
> architecture residuals remain beyond genus, but they do not identify pollinator
> identity or effective replacement. Measured-climate independence and an area/capacity
> pathway are not established.

H5 remains not evaluable in Chapter 1. Independent outcome-blind pollination-channel
exposure plus effort-aware retention/disruption or effective-service evidence is still
required to test the proposed channel gate directly.
