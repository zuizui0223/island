# Does mainland geography predict floral colour–form pairing beyond either marginal trait?

**10 November 2026 — empirical, exploratory, rights-restricted source replay.**
[Nonidentifying 12-row full result](marginal_fixed_placebo_aggregate.csv) ·
[Five-fold grouped audit](marginal_fixed_five_fold_audit.csv) ·
[IPF method and reproducible code](../../src/fd/marginal_fixed_coupling.py) ·
[Original fully joint source-placebo comparison](README.md).

## Why the full-rank joint kernel was not enough

The original 496-profile colour × form product-Jaccard kernel is full-rank
and its quadratic probability score is strictly proper; nevertheless,
a better *joint* score can come from better prediction of flower **colour
alone**, flower **form alone**, or their **cross-trait pairing**. It is
NOT sufficient to label a positive joint score a "floral trait
combination" effect.

We therefore perform a more restrictive same-island/source negative
control: **force competing models to have exactly the same predicted
colour-state-set and form-state-set marginal probability distributions**
while retaining differences in the *association* between the two
traits, and compare strictly proper scores for the original full
colour/form profiles.

## Fixed source and conditional-association design

- Exactly the same strict original-source-versus-two-placebo cohort:
  **1,007 All islands across 119 geographic blocks** and **950 nested
  Direct islands across 116 blocks**. No island/sample swapping or
  endpoint-driven placebo choice.
- Five heldout geographic folds, all tested accepted taxon names excluded
  globally from every GIFT native-mainland source before donor training,
  known evaluated island genus, source shrinkage **8**.
- Two alternative mainland-source assignments per target, from
  distinct geographic blocks in the same broad regime. They are
  selected WITHOUT flower outcomes and match independently available
  same-genus source donor coverage to ±2 percentage points and total
  complete-trait donor-name pool within a factor 0.5–2.
- For each heldout genus, obtain global-genus probability **g** and
  actual/alternative geographic-source probabilities **p** over the same
  496 intact observed flower-colour-set/form-set categories.
- Use iterative proportional fitting (IPF) to find the
  Kullback–Leibler-minimal multiplicative row/column rescaling of **p**
  for which *both* colour-set and form-set category marginals equal
  those predicted by **g**. The global-genus distribution itself is
  unchanged. On the structurally supported category matrix,
  the original source model's conditional association odds ratios are
  preserved by row/column scaling. Structural-zero categories are
  **never invented** from a colour × form Cartesian product.
- Hold the same target island/genus/name observations and two placebo
  assignments fixed. Compare mean placebo score minus mean actual
  score for (a) **full-rank product-Jaccard Kjoint** and (b)
  **full-profile strictly proper Brier**. Positive gain means actual
  geography adds joint-association prediction AFTER marginal colour
  and form probability distributions have been held fixed.
- IPF convergence tolerance **2×10⁻⁷**, 0 failures across all scored
  prediction contexts. Both older unprojected additive and joint
  source gains were reproduced to within 10⁻⁶, an essential
  same-model/sample identity check.
- Geographic uncertainty: **4,999 whole-spatial-block bootstrap
  resamples** of the strictly matched island population. These are
  *post-hoc* conditional intervals, not a preregistered randomized
  source experiment.

## Main results

Values below are the true geographic source **predictive LOSS advantage**
over two otherwise matched source assignments. Positive favours
the actual geographic source.

| Source-predictive score | All (n=1,007) | Direct (n=950) |
|---|---:|---:|
| Colour-only proper Jaccard kernel | +0.000991 [0.000315,0.001704] | +0.000596 [−0.000100,0.001258] |
| Form-only proper Jaccard kernel | +0.002729 [0.000678,0.004405] | +0.001258 [−0.000308,0.002582] |
| Original additive colour/form kernel | +0.001860 [0.000687,0.002877] | +0.000927 [−0.000031,0.001741] |
| **Joint product kernel; unadjusted marginals** | **+0.002378 [0.000878,0.003753]** | +0.001097 [−0.000153,0.002178] |
| **Joint product kernel; both marginals fixed** | **+0.000214 [0.000036,0.000370]** | +0.000140 [−0.000062,0.000331] |
| **Full Brier; both marginals fixed** | **+0.000205 [0.000022,0.000379]** | +0.000137 [−0.000087,0.000360] |

All five heldout geographic-fold averages show a small positive
marginal-fixed Kjoint gain in All and in Direct; one Direct fold is
**nearly zero (+0.00000184)**. The full Direct block-bootstrap
interval includes zero. **All/Direct overlap**, so do not multiply
the apparent number of independent confirmations.

The magnitude of the reported joint-product K-score gain decreases
by about **91.0% in All and 87.2% in Direct** when both trait
marginals are equated, leaving ~9% and ~13% of the original point
estimate respectively. These figures quantify **sensitivity of a
particular scoring rule to conditioning**, not biological proportions
explained, statistical mediation or causal trait selection.

## Conservative five-fold and paired trait-scale checks

The published spatial-block 95% intervals above treat 119 All / 116
Direct spatial blocks as resampling units, **conditional on fitted donor
models and the source-assignment matching design**. This does not make
the five cross-fitting model fits or shared donor taxa independent.

A second deliberately conservative check uses the **five geographic
test-fold means** as the five sign units, weighted by island counts.
With only 2^5=32 hypothetical sign configurations, the smallest
possible conventional two-sided sign tail for five concordant
nonzero folds is 2/32=**0.0625**:

| Strict source-placebo endpoint | All positive folds | Direct positive folds | Five-fold two-sided tail, if exchangeable |
|---|---:|---:|---:|
| Joint K, unrestricted marginals | 5/5 | 5/5 | 0.0625 / 0.0625 |
| **Joint K, both trait marginals fixed** | **5/5** | **5/5** | **0.0625 / 0.0625** |
| Full Brier, both trait marginals fixed | 5/5 | 4/5 | 0.0625 / 0.1875 |

The **fold exchangeability assumption is not established**: folds
share taxon traits and overlapping native-mainland donor training
opportunity. These sign tails are diagnostic stress statistics,
**not confirmatory permutation p-values** or proof of five independent
evolutionary replications.

Removing any one of the five test folds leaves positive
island-weighted marginal-fixed K-score means:
**All +0.000183 to +0.000249**, Direct **+0.000101 to +0.000175**.
This is directional stability, not identification of a causal effect.

The form-only K-score gain has a larger point estimate than the
colour-only gain, but **the paired difference on identical islands**
is not consistently positive under spatial-block uncertainty:
All form-minus-colour **+0.001738**, 95% block interval
**[-0.000221,+0.003288]**; Direct **+0.000662**,
**[-0.000808,+0.001967]**. The two marginal kernels also have
different baseline state-space structure, so raw source-skill
differences must not be interpreted as relative strengths of
natural selection on flower form and flower colour.

[Machine-readable five-fold stress table](fold_coarse_sign_sensitivity.csv) ·
[Paired form-minus-colour uncertainty](form_vs_colour_paired_source_gap.csv).
The coarse fold calculation can be reproduced directly from the
already-released anonymous five-fold table without accessing any
restricted source data:

    python -m fd.fold_coupling_sensitivity \
      --summary results/20261110_joint_flower_source/marginal_fixed_five_fold_audit.csv \
      --output results/coupling_fivefold_stress

All these checks are **post-hoc robustness analyses**, not separate
pre-registered hypothesis families. A full sensitivity of marginal-fixed
scores to donor-model shrinkage strength (e.g. 2/8/20) remains
**unevaluated** and is a further limitation before strong ecological
claims.

## Heterogeneity: does the small average gain characterize most islands?

The marginal-fixed source result is a **signed difference in heldout
prediction scores**, and meaningful island-level heterogeneity can be
hidden by a positive group mean. Using the **same strict 1,007 All /
950 Direct island cohorts**, and resampling **whole spatial blocks**
rather than raw islands, we computed medians, 10%-trimmed means,
10%-winsorized means, and the proportion of islands with positive
actual-source versus placebo advantages.

| Post-hoc marginal-fixed joint K-score diagnostic | All | Direct |
|---|---:|---:|
| Positive-score island fraction | **60.8%** | **54.4%** |
| Island-level median gain | +0.000104 | +0.000049 |
| 10%-trimmed mean gain | **+0.000165** | +0.000113 |
| 95% whole-block interval for trimmed mean | **+0.000012 to +0.000313** | −0.000079 to +0.000305 |
| 10%-winsorized mean gain | +0.000192 | +0.000129 |
| 95% whole-block interval for winsorized mean | +0.000018 to +0.000359 | −0.000070 to +0.000323 |

Large positive island deviations are offset by negative ones.
The top 10% of islands with the greatest gains contribute about
**48% of all positive island-level gains** in each scope. As a
share of the *net signed sum*, that group contributes **~82% All**
and **~120% Direct**; the latter can exceed 100% because other islands
have negative source effects. This concentration statistic is
descriptive and inherently unstable when the net gain approaches zero.

**Conclusion:** A small All conditional association advantage
is not driven *only* by single extreme outliers, because the
symmetrically trimmed average remains positive under the fitted-model
spatial-block bootstrap. But many recorded islands do not benefit
from the true source assignment, and Direct remains uncertain.
Therefore a claim of broadly shared joint floral-phenotype transfer
or universal source filtering would be too strong.

[Anonymous two-row results](coupling_heterogeneity_aggregate.csv) ·
[reproducible private-input summary calculation](../../src/fd/coupling_heterogeneity.py).
The per-island/block input tables are NOT published, and all model
and shrinkage limitations remain in force.

## Mechanistic interpretation and inferential stop

**Supported:** Geographic native-mainland donor lists contain
additional information on recorded flower colour and form prediction
beyond globally known genus; some information can be detected
in the conditional colour/form association under a very constrained
All comparison. The conditional-association gain is **an order of
magnitude smaller** than the original unadjusted joint-source gain,
and is not robustly separated from zero in Direct.

**Not established:** a globally transferable evolutionarily selected
flower colour/form syndrome, causal colonization-source identity,
within-lineage evolutionary change, pollinator-mediated selection,
or a geology-dependent adaptive collapse of floral niches. A
fully ranked strictly proper scoring kernel and valid IPF
marginal constraints cannot by themselves repair shared taxon
identities, GIFT source underreporting, geographic donor opportunity
or lack of independent phylogenetic ancestry. Results concern
recorded *complete categorical state sets* only.

**Manuscript decision:** Make the strong lineage/genus prediction
and the modest, observation-sensitive source geography effect the
primary macroecological findings. Present the small conditional
coupling result as **secondary**, with Direct uncertainty explicitly
shown. Do not recast this as demonstrated adaptive "combinatorial
floral evolution".

## Reproduce without redistributing licensed input

The authorized private CLI in
[src/fd/n2_source_assignment_placebo.py](../../src/fd/n2_source_assignment_placebo.py)
accepts the same source-pinned private inputs and adds
**--evaluate-marginal-fixed** to retain these grouped scores.
IPF and proper Brier/joint K-score evaluation live in
[src/fd/marginal_fixed_coupling.py](../../src/fd/marginal_fixed_coupling.py),
with synthetic tests proving that a *pure marginal change with
no association change* is removed and that changed association
with identical marginals remains detectable.

Only aggregate 12-row results and ten fold summaries are published.
No named plant, genus, island, source entity, specific geographic block
or private floristic identity is committed to GitHub or attached to
public Actions. Public CI validates the algorithm/synthetic null, but
does not run the rights-restricted source inference.
