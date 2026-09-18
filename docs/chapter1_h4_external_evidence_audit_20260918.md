# Chapter 1 H4 external evidence audit — 2026-09-18

## Purpose

This note separates four evidence roles that must not be conflated:

1. the locked v13 post-hoc exact-species functional bridge;
2. the prospectively specified post-2015 wild-plant temporal replication;
3. independent published analyses of trait–pollen-limitation relationships;
4. the separately frozen PolLimCrop crop-domain transportability test.

The goal is to strengthen triangulation without retroactively calling the v13 H4 confirmatory.

## 1. Locked v13 discovery remains post-hoc

The current Chapter 1 v13 functional bridge was designed after the earlier
`distance × trait` moderation results had been inspected. Its inferential role remains
`posthoc_functional_triangulation`.

No later analysis changes that label.

## 2. Prospective wild-plant temporal replication: stopped at support, before outcomes

The prospective protocol was frozen before post-2015 validation outcomes were extracted.
It retained the two co-primary predictions and their thresholds:

- H4a: autonomous selfing predicts lower current pollen limitation;
- H4b: an equal-weight generalized-form + actinomorphy + shallow/open-tube score predicts
  lower current pollen limitation;
- one-sided alpha = 0.025 for each;
- no threshold relaxation and no component dropping.

The final outcome-blind support decision is locked in
`config/chapter1_h4_prospective_temporal_support_decision_lock.json`.

Even the most permissive predeclared diagnostic ceiling reached only:

- H4a: 9 matched species across 6 publications;
- H4b: 9 matched species across 4 publications.

Both are below the frozen minimum of 30 species and 10 publications. Therefore the
post-2015 wild-plant cohort is **not evaluable**, and its pollen-limitation outcomes remain
unopened for the primary replication.

This is a support/design result, not a biological null.

## 3. Published external evidence: classify by outcome-data independence

### 3.1 Same or overlapping GloPL outcome base — corroboration, not independent replication

Burns et al. (2019, *New Phytologist*, DOI 10.1111/nph.15935) analysed global pollen
limitation and plant traits and treated autonomous seed set / autofertility as a key
moderator of pollen limitation.

A later global analysis (Nature Communications 2025,
DOI 10.1038/s41467-025-61032-5) integrated pollen-supplementation data with reproductive
and floral traits. It reported that autofertile plants were less pollen limited and that
specialized floral structures were more pollen limited.

Because these analyses use GloPL or substantially the same pollen-supplementation evidence
base as v13 H2/H4, they provide **published biological corroboration**, not an independent
outcome replication of v13.

### 3.2 Older regional/comparative analyses — prior external biological support

Wolowski et al. (2014, PLoS ONE, DOI 10.1371/journal.pone.0091238) analysed 126 Atlantic
Forest species and found higher pollen limitation in phenotypic/ecological specialists;
actinomorphic and more generalized plants had lower pollen limitation than specialized
counterparts.

Larson & Barrett (2000, Biological Journal of the Linnean Society,
DOI 10.1111/j.1095-8312.2000.tb01221.x) compared 224 animal-pollinated species across
breeding-system and floral-specialization classes.

These studies predate the current H4 formulation and support its biological direction,
but they were not prospectively selected as the present validation cohort and may overlap
historically with studies later incorporated into GloPL. They therefore remain
**prior/external concordance**, not the new confirmatory replication.

### 3.3 PolLimCrop — independent outcome domain

PolLimCrop (Siopa et al. 2023, Scientific Data,
DOI 10.1038/s41597-023-02797-6) contains a distinct agricultural
pollen-supplementation dataset: 1,169 experiments from 294 studies and 108 crops.

Its published effect size is a log response ratio of hand pollen supplementation to
natural pollination, with positive values indicating stronger pollen limitation.

A 2026 PNAS analysis of the crop dataset
(DOI 10.1073/pnas.2533418123) independently reports that autogamy capacity is associated
with lower pollen limitation. This is strong external support for the H4a direction.

However, because that aggregate published result became visible after the H4 predictions
were frozen but before the row-level crop analysis is completed, the repo does **not**
present it as the result of our prospective test. The predeclared PolLimCrop analysis
remains governed only by the frozen support gate and response-mapping contracts.

## 4. Evidence-role table

| Evidence layer | Outcome data independent of v13 GloPL? | Prospectively frozen here? | Allowed interpretation |
|---|---:|---:|---|
| v13 exact-species bridge | No | No | post-hoc functional triangulation |
| post-2015 wild temporal cohort | Yes | Yes | support insufficient; outcomes remain unopened |
| Burns 2019 / Nat Commun 2025 | No / substantially overlapping | No | published corroboration |
| Wolowski 2014 / Larson & Barrett 2000 | partly external, possible historical overlap | No | prior biological concordance |
| PolLimCrop row-level test | **Yes** | **Yes, before row-level outcome read** | secondary cross-domain transportability if support gate passes |
| PNAS 2026 PolLimCrop aggregate result | Yes | No (published externally) | independent published corroboration of H4a only |

## 5. Manuscript claim ceiling

The strongest safe wording before a successful support-admitted PolLimCrop row-level test is:

> The post-hoc v13 trait-state associations align with a substantial prior literature in
> which reproductive assurance and generalized floral access are associated with lower
> pollen limitation. A newly preregistered post-2015 wild-plant temporal replication was
> stopped before outcome unblinding because its frozen support thresholds were not met.

Do **not** write:

- that the wild temporal replication failed biologically;
- that H4 is now confirmatory;
- that historical pollen limitation caused the island syndrome;
- that the published crop result validates the wild-island mechanism;
- that specialized-pollinator loss was identified.

If the frozen PolLimCrop row-level analysis is support-admitted and passes, add only:

> The same predeclared trait-to-current-pollen-limitation direction transported to an
> independent agricultural pollen-supplementation domain.

That still does not repair the support-limited wild temporal replication or identify
historical mediation.
