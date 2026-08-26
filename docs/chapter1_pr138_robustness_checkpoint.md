# PR138 syndrome robustness checkpoint

Status: **provisional / not frozen**.

This checkpoint records robustness of the pattern-first Chapter 1 result after adding formal pollination/selfing syndrome analysis. It does not identify causal pollinators.

## Headline being stress-tested

Northern mid-latitude island floras shift with mainland distance toward a predeclared classic island-syndrome direction:

`(- large_bee_like + generalized_accessible + selfing_syndrome) / 3`

while tropical island floras do not show the same projection and their complete syndrome response vector differs from the northern vector.

## 1. Evidence-layer and information-weight robustness

Main analysis uses nonduplicated all-analysis-eligible evidence:

- species-direct evidence first;
- validated genus-consensus evidence only fills species x trait gaps;
- 132,288 direct rows expand to 247,450 nonduplicated species x trait rows.

Direct-only independently reproduces the northern classic projection. Validated genus-consensus-only is weaker and does not independently recover the northern directional projection.

Changing information weights while holding island trait shares fixed (canonical, cap 100/50/20, equal-island) keeps the northern classic direction positive; equal-island strengthens the north-vs-tropical difference.

## 2. Syndrome-template sensitivity

Run: `32927708306`
Artifact: `pr138-syndrome-template-sensitivity-32927708306`
Digest: `sha256:060c94b1131d9b6a1655ed117d36afb70c4f243041b4d9d0ea2dc5606bc97a01`

Thirteen outcome-blind template variants were prespecified:

- canonical;
- equal weights;
- no colour;
- pollination morphology only;
- drop each syndrome trait in turn;
- require >=3 informative traits.

All-native:

- northern classic projection supported: **13/13 variants**;
- north-vs-tropical vector difference supported: **13/13**;
- `large_bee_like < 0` and `generalized_accessible > 0`: **12/13**.

Native non-endemic:

- northern classic projection supported: **11/13**;
- north-vs-tropical vector difference supported: **12/13**;
- attraction-direction signs supported: **12/13**.

The two native-nonendemic variants that lose the classic-projection significance are:

- `drop_floral_form`: projection +0.0302, q=0.0529; regional vector difference remains supported (q=0.000955);
- `minimum_three_traits`: projection +0.0391, q=0.0689; regional vector difference q=0.0511.

Thus the result is not colour-driven or dependent on a single reproductive trait, but floral form and the stricter >=3-trait requirement are identifiable sensitivity points.

## 3. Distance functional-form sensitivity

Run: `32927708299`
Artifact: `pr138-syndrome-distance-sensitivity-32927708299`
Digest: `sha256:4eba14ceec8ed338f3d8c06194a62ce4abf0b86b2ebab4e4ecf4e24eb2ccdca1`

The island set is unchanged and no distance threshold is used.

Across raw, square-root and log1p mainland distance, in both all-native and native-nonendemic strata (**6/6 scenarios**):

- northern classic projection is positive and supported;
- tropical classic projection is not supported;
- north-vs-tropical syndrome vector difference is supported;
- northern `large_bee_like` slope is negative;
- northern `generalized_accessible` slope is positive.

Example all-native classic projection:

- raw: +0.0626, q=4.35e-11;
- sqrt: +0.0628, q=3.15e-5;
- log1p: +0.0443, q=0.00641.

Therefore the headline is not specific to the log-distance transformation.

## 4. Leave-one-spatial-block-out sensitivity

Run: `32927708223`
Artifact: `pr138-syndrome-block-deletion-32927708223`
Digest: `sha256:3f94d1b47a5003e3302be6792d55e19db2681733cdbcd866774dbdcbad577b90`

There are 81 eligible spatial blocks. Every block is deleted in turn and the confirmatory FDR family is refit.

For both all-native and native-nonendemic:

- headline is testable: 81/81 deletions;
- north-vs-tropical vector difference retained: **81/81**;
- northern attraction direction (`large_bee_like < 0`, `generalized_accessible > 0`) retained: **81/81**;
- northern classic projection retained: **80/81**.

The single influential deletion is `lat12_lon20`, an Aegean/eastern-Mediterranean block (~35-40 N, 20-30 E; 83 islands in the full geography table). When this block is removed:

- all-native classic projection +0.0244, q=0.2145;
- native-nonendemic +0.0221, q=0.2303.

However the attraction-direction signs remain correct and the north-vs-tropical vector difference remains strongly supported (q=3.1e-5 all-native; q=5e-6 native-nonendemic).

Therefore no single spatial block is necessary for the regional syndrome-vector difference, but the strength/significance of the scalar classic projection is materially leveraged by the Aegean/eastern-Mediterranean block and must be reported.

## 5. Formal biogeographic-realm sensitivity

Run: `32927060210`
Artifact: `pr138-realm-sensitivity-32927060210`
Digest: `sha256:6556a7938a52547ddf7b630a822d743b2ba5ca479589716919441b96dcfa0521`

RESOLVE Ecoregions 2017 realms are assigned by point-in-polygon only, with no nearest-realm imputation. 5,583/8,265 islands receive a formal realm assignment.

Palearctic strongly reproduces the northern classic syndrome:

- all-native projection +0.1027, q=3.46e-6;
- native-nonendemic +0.1017, q=9.35e-7.

Palearctic all-native slopes include:

- `large_bee_like` -0.1166;
- `generalized_accessible` +0.1086;
- `selfing_core` +0.0784;
- `selfing_syndrome` +0.0828;
- `bird_like` -0.1376;
- `butterfly_like` -0.0825.

Neotropical all-native shows a weaker classic-direction projection (+0.0325, q=0.0194), but this disappears in native non-endemics (+0.0175, q=0.0724). Palearctic-vs-Neotropical syndrome vectors differ directly in native non-endemics (q=0.00401), whereas the all-native vector difference is not supported (q=0.192). This indicates a material endemic/floristic-composition contribution to the Neotropical all-native signal.

The source-backed Nearctic status pilot was then rerun after an outcome-blind 15-island source inventory ([`chapter1_pr138_nearctic_wave1_source_inventory.md`](chapter1_pr138_nearctic_wave1_source_inventory.md)) and a direct audit of Le Hors's Saint-Pierre flora ([`chapter1_pr138_saint_pierre_source_audit.md`](chapter1_pr138_saint_pierre_source_audit.md)).

- Run: `32942359223`
- Artifact: `pr138-nearctic-status-pilot-32942359223`
- Digest: `sha256:26d4653bdbb4b2b59f51d86745e596ca2620acfd911a4bc3604d8b7968aa37f1`
- San Nicolas: 277/355 status rows resolved from the author-maintained CCH2 checklist (136 native under the source convention, 141 explicit introduced, 78 unresolved).
- Block Island: 165/286 status rows resolved from the explicitly defined Appendix A convention (105 unstarred native, 60 starred naturalized/not native, 121 unresolved).
- Cedros: none of 338 rows enters the binary origin analysis because the source's starred taxa are only *presumably* introduced and it provides no native residual convention.
- Saint-Pierre: only six exact-name rows enter from Le Hors passages that directly couple origin with target-island occurrence (five native, one introduced); the remaining 305 rows stay unresolved. The archipelago-wide total of 506 indigenous taxa is not used as a residual-island rule.

This addition raises Nearctic all-native support to the predeclared 30-island pilot gate for three axes: `butterfly_like`, `generalized_accessible`, and `large_bee_like`. The joint three-axis vector is now testable but unsupported (31 islands, 10 spatial clusters; chi-square = 3.616, df = 3, p = q = 0.306). The fitted distance slopes are `butterfly_like` +0.0021 (p = 0.985), `generalized_accessible` -0.0236 (p = 0.205), and `large_bee_like` +0.0483 (p = 0.261). None reproduces the Palearctic classic direction with statistical support. This is a pilot-scale multivariate non-replication, not a confirmatory falsification: the 50-island confirmatory gate is unmet, the classic scalar projection is incomplete because the reproductive-assurance axes remain below threshold, and native-nonendemic Nearctic remains untestable (27 islands for butterfly/generalized/large-bee and 26 for bird/selfing).

## 6. Attraction versus strict selfing-core decomposition

Run: `32922758001`
Artifact: `pr138-attraction-selfing-32922758001`
Digest: `sha256:5c573fdd75ebab4b85dd97cd11a064d9764eb3772c8fa37e124236eb841e6770`

Definitions:

- `attraction_shift = (-large_bee_like + generalized_accessible) / 2`;
- `selfing_core = SC + mating system + autonomous selfing`, excluding flower size.

Northern all-native:

- `selfing_core ~ distance`: +0.0224, p=0.347, unsupported;
- `attraction_shift ~ distance`: +0.0466, p=0.0133;
- after selfing-core adjustment, distance remains +0.0471, p=0.0149;
- conditional selfing-core coefficient on attraction shift: p=0.720.

Native non-endemic gives the same qualitative result (conditional distance +0.0435, p=0.0168).

Thus the northern attraction/accessibility shift is statistically separable from the measured strict reproductive selfing core. This is conditional decomposition, not causal mediation.

## Current scientific boundary

The strongest supported statement is:

> Mainland distance is associated with a biogeographically contingent reorganization of floral/pollination-syndrome composition. In northern mid-latitude, especially Palearctic, island floras shift away from large-bee-compatible architecture and toward generalized accessibility; this floral shift is not reducible to the measured strict selfing core. Tropical floras follow a different syndrome-response vector.

Still unresolved:

- direct causal attribution to Bombus loss or attraction-selection relaxation;
- confirmatory Nearctic replication and native-nonendemic Nearctic replication; the all-native pilot vector is now testable but does not reproduce the Palearctic direction;
- confirmatory southern-hemisphere bird-like pattern;
- why the Aegean/eastern-Mediterranean block contributes strongly to the scalar classic projection.

These are the next targets; the PR remains draft and the scientific result is not frozen.
