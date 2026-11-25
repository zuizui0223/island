# How much does mainland source geography add after donor opportunity is matched?

**10 November 2026 — post-hoc matched-control macroecological prediction audit.**
This is a **source-matched probability-prediction sensitivity**, not a
geological origin test, evolutionary selection evidence or actual dispersal
source reconstruction.

[Reproducible authorized-data-only implementation](../../src/fd/n2_source_assignment_placebo.py) ·
[synthetic control/rights tests](../../tests/test_n2_source_assignment_placebo.py) ·
[aggregate source-score results](coverage_matched_summary.csv) ·
[five nonidentifying geographical fold summaries](fold_direction_support.csv).

## Question

Earlier N2 geographically blocked taxon-name holdout showed positive source
assignment gains in the original categorical-trait similarity kernel K-score.
The true nearby mainland assignment scored better than two randomly borrowed
same-regime/different-spatial-block mainland-source assignments:

- All: mean advantage **+0.003961** (spatial-block bootstrap 95%
  interval **+0.002189 to +0.005701**, n=1,521 island targets)
- Direct: **+0.003518** (**+0.001359 to +0.005670**, n=1,383)

**But the true assignment had more independently trainable same-genus donor
coverage:** +5.37 percentage points in All and +3.52 in Direct. Therefore
this test alone does not identify source geography as a novel ecological
driver. Is an advantage retained if two alternative mainland source
assignments are selected on **source-support opportunity rather than
flower outcomes**?

## Fixed matching requirements

- Original 496 intact recorded floral colour-state-set × form-state-set
  profiles and full original set-Jaccard D; K=1-D is PSD (rank 76), so the
  quadratic K-score is a **proper but not strictly proper** soft-category
  score. Brier is the strictly proper probability-score check.
- Use five **spatial-block-disjoint** test folds, with every tested taxon
  name globally removed from EVERY native GIFT source/flower-trait donor
  before each fold's model fitting.
- Keep observed island genus/counts known, global-genus baseline, source
  trait profiles, and shrinkage = **8** fixed. Source geography predicts
  a candidate native-mainland genus trait mixture, NOT genus colonization.
- For each test island choose candidate donor assignments of **other test
  islands in the same broad biogeographic regime and heldout fold**, but
  from a **different geographic source block**. Their GIFT source entity
  set must not be identical to the focal island's true entity set.
- Select exactly **two** control assignments using NO flower trait
  observations or prediction losses. Condition A: the fraction of focal
  genus-name carriers with independent local-genus trait donors is within
  **±2 percentage points** of the focal true assignment. Among eligible
  control candidates, take the two with the nearest coverage, then nearest
  total donor-count ratio (stable tie order).
- Condition B applies A **plus** each selected alternative source
  assignment's total *independent complete native donor taxon-name count*
  being between **0.5× and 2×** that of the true assignment.
- If a focal has fewer than two qualifying alternatives, do not impute a
  nonexistent reference or claim its null effect. The selection gate is
  explicitly **post hoc** and changes the analysis population.

All outputs below average K-score or Brier loss **across target islands**
with each island given one vote (even if it contains more taxon-name
records), and confidence bounds resample **whole geographic spatial
blocks** 4,999 times. All/Direct are overlapping evidence scopes.

## Results (source assignment advantage, positive is better)

| Focal comparison | All | Direct |
|---|---:|---:|
| Original two random same-regime source controls | +0.003961 | +0.003518 |
| Matched on same-genus donor coverage ±2pp | **+0.002439** | **+0.001794** |
| Number of islands with coverage-matched controls | 1,137 (127 blocks) | 1,077 (129) |
| Coverage + independent donor-name pool within 0.5–2× | **+0.001860** | **+0.000927** |
| **95% block-bootstrap range, strict match** | **+0.000687 to +0.002877** | **−0.000031 to +0.001741** |
| Number of strictly matched islands | **1,007 (119 blocks)** | **950 (116)** |
| Strict-match **Brier-score** advantage | **+0.002011** | **+0.000811** |
| Strict-match Brier 95% interval | **+0.000788 to +0.003295** | **−0.000375 to +0.001908** |

The strict control opportunity match leaves mean focal-minus-control
same-genus donor coverage differences only **0.000503 All / 0.000603
Direct** as proportions (roughly **0.05–0.06 percentage points**).
The median control/focal available donor-name count ratio is **1.032
All / 1.029 Direct**. By design the control source entities are
not identical; their mean source-entity Jaccard overlap with true
assignments is about **0.144 / 0.138**.

Under the strict joint caliper, the K-score advantage is **positive in
5/5 All** folds and **5/5 Direct** folds. However, the independent
whole-spatial-block 95% interval still touches/includes zero in Direct;
the two scopes cannot be counted as independent confirmations. Leaving
out any ONE spatial block gives All mean K advantage in
**[0.001486,0.002009]**, Direct **[0.000614,0.001062]**. This is an
influence check, not a p-value or independent replication.

The raw-to-strict point-estimate reduction is roughly **53% All /
74% Direct**, but the cohorts differ due to support restrictions.
Therefore those percentages **must not be described as causal proportions
explained by donor availability**.

## Interpretation

The distinctive result is **a much smaller and scope-sensitive
geographical source identity signal once a fairer mainland donor
comparison is made**. All retains a modest positive probability-kernel
gain and Brier gain on support-matched comparisons. Direct's
conservative geographical-block bootstrap intervals overlap zero.
This weakens any claim that particular nearby mainland source regions
systematically organize recorded island floral traits beyond genus
availability.

Even strict matching does **not** equalize within-genus taxon composition,
genus-specific sample sizes, phylogenetic relationships, plant
observation probabilities, source-country habitat/climate, or true
ancestral source identity. The original sharing of accepted names
across islands persists (although test names were removed from
training). Source placebo controls are assigned from other islands and
do not constitute naturally randomized island colonization histories.

**Appropriate manuscript claim:** Recorded floral phenotypes have
transferable genus- and geography-associated *similarity* structure,
but the apparent extra skill from naming the nearby native mainland
source is small and partially contingent on the independent
source/donor evidence available. The robustly established component
remains genus-linked predictive structure; any residual source
identity effect needs richer external donor trait observations.

**Not supported:** geological origin adaptive floral convergence,
pollinator adaptation, source-independent niche-space contraction,
or confirmed mainland-colonist lineage history. No new geological
p-values or causal null rejections were generated.

## Data rights and reproduction

The source input consists of licensed PR321 reconstructed complete
taxon×floral-profile maps, original island×taxon records, archived
GIFT native mainland taxon×source entity lists, and geo_k5
source assignment tables. These private identities are **never**
written to this repository or user-facing Actions artifacts.
The implementation uses these tables only transiently and writes
grouped estimates and a SHA source audit.

The public CI validates synthetic matching rules, reproducible
geographic clustering, and privacy gates. The empirical summary was
independently computed in a restricted local workspace using a
low-rank factorization of the frozen PSD similarity kernel; the
computational shortcut preserves the same original K-score algebra.
The source-specific full CLI replay is available, but no new claim
of a completed full private-input GitHub Action is made.
