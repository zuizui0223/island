# Chapter 1 existing-data result: what is resolved about when and where

## Canonical verdict

After floristic status and genus composition are represented, the existing data
do **not** recover a broad context-specific floral syndrome. They do recover a
strong Northern isolation-associated floristic/endemicity gradient and one small
Northern colour-category residual. Thus the current Chapter 1 result is best
described as **context-dependent floristic restructuring with limited residual
floral filtering**, not as a Northern Bombus-driven syndrome.

The result is presented in that order: establish the compositional pattern
first, then evaluate pollination as a candidate explanation. The current data do
not permit the candidate explanation to be promoted to the result.

## Why 8,265 candidate islands become 4,370 analysis islands

The reduction is now set-audited rather than inferred from row counts.

| Stage | Islands | Loss |
|---|---:|---:|
| frozen GSHHG/covariate candidates | 8,265 | — |
| at least one exact-polygon GBIF plant record | 4,549 | 3,716 |
| at least one accepted-species row | 4,505 | 44 |
| at least one strict-universe species | 4,453 | 52 |
| at least one strict species with a resolved configured direct trait | 4,370 | 83 |
| after covariate join | 4,370 | **0** |

The old `83 lost at covariate join` label is false: all 4,453 strict-flora
islands occur in the 8,265-island covariate table. Those 83 islands lack a
resolved configured trait, not covariates. The main loss—3,716 islands—is absence
of an exact assigned GBIF plant record and cannot be interpreted biologically.

Final retention is geographically structured: Northern 69.5%, Northern high
latitude 30.5%, Tropical 49.6%, and Southern 37.9%. It also rises from 37.7% in
the smallest area quartile to 77.9% in the largest. The 4,370-island pattern is
therefore conditional on observation availability. The detailed audit and
prespecified correction strategy are in
[`chapter1_island_universe_attrition.md`](chapter1_island_universe_attrition.md).

## Primary evidence: status and lineage workflow

The canonical status/lineage results are from the PR #133 workflow verified at
commit `f9d070440ce42a213507b26f42db43830e19e1b8`.

| Question | Result | Interpretation |
|---|---|---|
| Does isolation restructure floristic status? | Northern endemicity: beta `+0.902`, robust SE `0.191`, p `2.4e-6`, 222 islands; at at least 95% resolved status, beta `+0.514`, p `0.0016`, 139 islands | robust Northern floristic/endemic turnover |
| Do broad native-nonendemic traits remain after genus control? | Northern and Tropical plain colour, generalized form, and SC are all null; Southern estimates are also null but supported by only 30–33 islands | no broad residual floral/reproductive syndrome |
| Does fine-category decomposition recover hidden responses? | Northern all-native red/pink residual declines: beta `-0.00251`, p `0.00369`, q `0.0147`, 238 islands | small, coverage-robust colour turnover; not a multivariate syndrome |
| Is there a robust tropical counter-syndrome? | Tropical native-nonendemic red/pink declines at baseline but collapses in the at least 40% direct-coverage subset | measurement-sensitive; not a manuscript-level regional syndrome |
| Do floral forms retain a signal? | no floral-form category with at least 50 islands survives FDR | no confirmatory form residual |
| Did the 30-island endemic-form pilot recover a tropical residual? | tropical endemic generalized form: beta `-0.04319`, robust SE `0.02790`, p `0.1216`, 29 complete-case islands/19 blocks | unsupported; acquisition remains stopped at 30 rather than expanding to 50 |
| Does Bombus compatibility add the expected pathway? | generalized form and SC do not improve out of sample; no distance-by-compatibility interaction is nominally supported | predicted Bombus proxy pathway not recovered |

The broad native-nonendemic coefficients are explicitly small and unsupported:

| Region | Plain colour p | Generalized form p | SC p |
|---|---:|---:|---:|
| Northern | 0.238 | 0.432 | 0.811 |
| Tropical | 0.809 | 0.688 | 0.323 |
| Southern | 0.856 | 0.726 | 0.878 |

Southern support remains below the 50-island confirmatory target, so its nulls
are support-limited rather than strong equivalence evidence.

The newer PR #135 post-pilot artifact also retains nominal Northern endemic
coefficients for generalized form (beta `+0.0804`, p `0.0155`) and SC (beta
`-0.1647`, p `0.0160`), but these use only 20 islands/12 blocks and 12 islands/7
blocks. They remain explicitly below pilot support and are not biological
findings. Tropical endemic SC is also unsupported (beta `-0.0967`, p `0.1405`,
27 complete-case islands).

## What is resolved about where

The strongest defensible geographic result is Northern midlatitude **floristic
status turnover**. Isolation is associated with greater regional endemicity,
and that result survives a strict status-resolution filter.

For residual floral traits:

- Northern has one small all-native red/pink decline after genus control;
- Tropical has no coverage-robust counter-syndrome;
- Southern has insufficient status-stratified support for a strong regional
  conclusion; and
- no region currently shows a supported multivariate residual trait vector that
  warrants the word `syndrome`.

Thus Chapter 1 currently answers where structural turnover is strongest more
clearly than where a floral syndrome emerges.

## What is resolved about when

The current meaning of “when” is a statistical/ecological condition, not
historical time:

```text
when region, floristic status, and lineage composition are represented,
which isolation-associated trait residuals remain?
```

The answer is: very few. This narrows the condition for a floral-filtering claim
to signals that remain after status, genus, coverage, spatial, and multiplicity
controls. The current data do not determine when a pollinator was lost, when a
trait evolved, or whether filtering occurred at dispersal, establishment,
persistence, or post-colonization evolution.

## Secondary evidence: pre-status 4,370-island analysis

The replaceable Chapter 1 trait pipeline found a stable Northern pattern toward
more `very_small` and fewer `very_large` flowers across 4,370 islands, with
different slopes in Tropical and Southern comparisons. That result remains
useful for trait prioritization and for deciding which new axes should enter the
status/lineage workflow.

It is not yet the canonical syndrome result because that analysis adjusts for
area and climate but does not separate floristic status or subtract a
genus-fixed expectation. It cannot overrule the primary status/lineage result.
Flower size, tube depth, symmetry, autonomous selfing, and mating system should
be promoted only after they run through the same structural controls.

## Pollination interpretation after the pattern

The candidate causal reading is deliberately conditional:

- if a coherent Northern simplification vector survives, loss of a dominant
  Bombus/large long-tongued function without effective replacement is a plausible
  formation mechanism;
- if that vector is weak or absent elsewhere, persistence or island colonization
  by alternative guilds is a plausible non-formation mechanism; and
- butterflies, hawkmoths, and birds imply different predictions—chromatic
  signalling, long-tube/outcrossing maintenance, and effective vertebrate
  replacement, respectively.

The literature supports examples of specialized moths repeatedly colonizing
remote islands, butterfly learning of floral colour signals, and colonist birds
partially compensating for lost endemic pollinators. It does not support the
blanket statements that all butterflies or birds have higher dispersal than
Bombus, or that Tropical/Southern islands must evolve showier flowers. The safe
counterprediction is weakening, absence, or heterogeneity of a Northern vector.

This is the intended use of pollination in Chapter 1. The next primary analysis
is not a Bombus-deficit model; it is the status/lineage-aware joint analysis of
the already collected size, depth, symmetry, form, colour, compatibility,
autonomous-selfing, and mating-system axes. Pollinator distributions remain
Discussion-level causal candidates.

## Bombus and hypervolume verdict

The climatic hypervolume is retained as a falsified supporting measurement model.
In the status/lineage residual checkpoint:

- adding compatibility worsened or failed to improve spatial prediction for
  generalized form and SC;
- plain colour improved by only about 1.1–1.7%, and the compatibility coefficient
  was opposite to the simple Bombus-deficit prediction; and
- none of the distance-by-compatibility interactions was nominally supported.

Therefore Bombus/hypervolume does not explain the current Chapter 1 result and
is not a direct Chapter 1 variable. It remains available for Discussion-level
context and provenance; no occurrence rescue campaign is justified.

## Hypothesis decisions

| Hypothesis | Verdict |
|---|---|
| H1 universal floral island syndrome | **not supported** after status and lineage control |
| H2 biogeographic/floristic contingency | **supported for floristic/endemic turnover; weak for residual floral traits** |
| H3 fine-category residuals | **one small Northern colour result; no multivariate syndrome** |
| Dissertation-level channel-gated assembly | **not directly tested by Chapter 1** |

This is a publishable boundary result: much of the apparent island floral
pattern is structurally tied to which floras and lineages occur, while little
residual trait filtering remains on current direct evidence.

## What remains unresolved, and how it enters the programme

The first and last rows are Chapter 1 completion work. The intervening causal
rows set the interpretation ceiling and are not prerequisites for this chapter.

| Unresolved problem | Next identifying evidence or sensitivity |
|---|---|
| Is there a multivariate floral syndrome? | move pending size/depth/symmetry/selfing axes through status, genus, FDR, and coverage gates |
| Is a Northern pattern caused by Bombus-channel loss? | Discussion hypothesis now; direct source exposure and realized functional deficit would be later causal evidence |
| Why might the pattern fail elsewhere? | Discussion hypothesis now; visitation times per-visit effectiveness would be later causal evidence |
| What part of distance is source pool/connectivity? | retain distance as the Chapter 1 total gradient; use status/genus controls now and explicit source-pool reconstruction only for mechanism attribution |
| Assembly filtering or island evolution? | status-aware phylogenetic/population contrasts and measured within-lineage floral divergence |
| Does incomplete island observation alter where/when? | model 8,265-to-4,370 selection and run frozen IPW/overlap sensitivities by region |

## Reproducible evidence

- PR #133: `https://github.com/zuizui0223/island/pull/133`
- Verified head: `f9d070440ce42a213507b26f42db43830e19e1b8`
- Current-main post-pilot port: PR #135,
  `https://github.com/zuizui0223/island/pull/135`
- Verified post-pilot head/run/artifact:
  `db1863f4057277b7a1d4f712918b22c8bdddc05d` / `32721017664` /
  `endemic-pilot-lineage-32721017664`
- Endemism checkpoint: run `32559322028`
- Status-stratified lineage checkpoint: run `32613129024`, artifact
  `status-stratified-lineage-32613129024`
- Category and coverage checkpoint: run `32614027696`, artifact
  `status-category-decomposition-32614027696`
- Northern Bombus residual checkpoint: run `32613425082`, artifact
  `northern-bombus-lineage-residual-32613425082`
- Secondary pre-status evaluation:
  `analysis/chapter1_channel_gate_evaluation_20260825/`
