# Chapter 1 H4 prospective validation summary — 2026-09-18

Status: **both prospective validation routes stopped at predeclared outcome-blind support gates**.

This document is an evidence-boundary audit. It does not change the post-hoc H4
discovery estimates and it does not treat support insufficiency as a biological null.

## 1. Post-2015 wild-plant temporal replication

The wild temporal protocol was frozen before validation outcomes. Under the most
permissive predeclared support ceiling:

- H4a autonomous selfing: 9 matched species across 6 publications
  (state 0: 7 species; state 1: 2);
- H4b accessibility/generalization: 9 matched species across 4 publications
  (low score: 2 species; high score: 2).

Both were below the frozen minimum of 30 matched species and 10 publications, with
the state/score-balance requirements also unmet. Outcome extraction was therefore not
authorized and the primary post-2015 pollen-limitation outcomes remained unopened.

Canonical decision:
`config/chapter1_h4_prospective_temporal_support_decision_lock.json`.

## 2. Independent PolLimCrop transportability audit

The crop-domain protocol was also frozen before row-level pollen-limitation outcomes.

The currently resolved Figshare v1 primary CSV is
`1PUB_PolLimCrop_dataset.csv` (article 24625299, file 43269153). Outcome-blind source
inspection found 1168 rows, 293 nonblank article codes, 106 crop species and 62
countries. The published Data Descriptor reports 1169 experiments, 294 studies,
108 crops and 62 countries. The exact locked source-byte scope was retained and the
discrepancy recorded rather than forcing the published aggregate counts.

Support under the frozen crop gates:

- H4a reproductive assurance: 36 matched species / 109 publications, but only
  **9 state-0 species** versus the required minimum of 10;
- H4b accessibility/generalization: 42 matched species / 141 publications and score
  SD 0.353, but only **6 low-score species** versus the required minimum of 8
  (high-score species = 17).

Thus neither co-primary hypothesis was admitted. The support lock records:

- `admitted_hypotheses=[]`;
- `outcome_extraction_authorized=false`;
- `outcomes_read=false`;
- no threshold relaxation;
- no trait-definition change.

The follow-up transportability workflow verified that support decision and skipped the
trait-download, outcome materialization, outcome-model and result-freeze steps.

Canonical crop support decision:
`config/chapter1_h4_pollimcrop_transportability_preflight_result_lock.json`.

## 3. Integrated interpretation

The two prospective routes now have the same evidential status for different reasons:

1. the post-2015 wild cohort is too sparse overall;
2. the crop domain has adequate total overlap but fails the prespecified balance/tail
   requirements.

Neither failure tests whether the biological H4 prediction is true or false. No
prospective pollen-limitation outcome was unblinded in either route.

Therefore the current H4 evidence hierarchy is:

- **post-hoc discovery:** retained;
- **exact H2-score functional bridge:** retained as post-hoc triangulation;
- **atomic reconstruction sensitivity:** retained;
- **prospective wild validation:** not evaluable because of frozen support limits;
- **independent crop transportability validation:** not evaluable because of frozen
  support limits.

The strongest allowed statement is that the post-hoc H4 functional bridge remains
biologically coherent but has not yet received an evaluable prospective replication.
The failed admission gates strengthen the transparency of that boundary; they do not
convert the discovery into a null result.

Machine-readable integrated lock:
`config/chapter1_h4_prospective_validation_summary_lock.json`.
