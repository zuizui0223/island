# Genus-conditional native-mainland source signals in joint flower colour × form profiles

**10 November 2026 JST — completed strict-support source-placebo replay, observational.**

[Aggregated 20-row proper-score results](coverage_matched_trait_kernels.csv) ·
[Marginal/joint kernel code](../../src/fd/trait_kernels.py) ·
[Five-fold source-placebo implementation](../../src/fd/n2_source_assignment_placebo.py) ·
[Mathematical and inferential contract](../../docs/FLORAL_TRAIT_COUPLING.md).

## Central question

The original source-conditional floral distance score is based on

    Kadd(a,b) = (Jaccard(colour_sets_a,colour_sets_b)
               + Jaccard(form_sets_a,form_sets_b))/2.

Thus source-predictive Kadd gain is *exactly* the mean of separate colour-only
and form-only source gains. By itself it says nothing about the **pairing**
of colour and form. We now ask whether a geographic mainland-source assignment
provides additional **proper joint colour/form probability-prediction skill**
over two nonlocal source assignments matched to the same independently
available genus/floral donor opportunity.

Use Kjoint(a,b)=Kcolour(a,b)*Kform(a,b), the Hadamard product of the original
recorded-state-set Jaccard kernels. The frozen 496-profile public reference
contains 31 unique colour-state sets and 46 form-state sets. Its numerical
kernel ranks are **31 colour, 46 form, 76 additive, and 496 joint**;
Kjoint has positive minimum eigenvalue **+0.003895479**. Hence on the frozen
496 full flower-profile categories, **Kjoint's quadratic score is strictly
proper** for the predicted probability vector. It distinguishes different
colour/form *couplings* even when both separate trait marginals are identical.
This is a mathematical property, not a result of flowering evolution.

## Reproducibility/rights audit

Independently reconstructed the frozen [island PR321](https://github.com/zuizui0223/island/pull/321)
upstream raw Direct ledger and progressive All/Direct trait snapshots.
Full SHA256 digests matched the archived source/progressive checkpoint files:

- source checkpoint: \`5e9c887e093c42eb720f3abd2dc85701ce123fb7fe888c916dc4098ac898ccdd\`;
- progressive checkpoint: \`b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999\`.

The reconstructed **496 state-set profile ID → colour/form category** tuples
match the public pinned reference across all 496 IDs, independent FNV32
\`c995f500\`; the counts of recorded complete-profile names per ID
match independent checksum \`26e99151\`. Complete names remain **16,739
All / 8,455 Direct**, 496/467 distinct known profiles.

Native GIFT entity×taxon source records remain rights-restricted; **nothing
below publishes species/taxon names, individual island identities, genus
labels, GIFT source IDs or geographic block names.** Only grouped scores and
block-bootstrap intervals are retained in this repository. The private
name-level replay happened locally; public GitHub Actions test the algorithms
and public kernel reference, NOT a private GIFT full-source CI run.

## Source placebo comparison

- Use the existing **five disjoint geographic test-block folds**, removing
  every target-fold taxon name from EVERY GIFT native mainland donor entity
  before training. All/Direct are nested evidence scopes.
- Keep each evaluated island's **observed genus/carrier counts** known, and
  the same global genus probabilities, local source donor profiles and
  hierarchical shrinkage **8** used in earlier results.
- For each focal island, choose two **alternative source assignments**
  from other islands in the SAME fold and broad biogeographic regime but
  DIFFERENT spatial block and a genuinely different GIFT source-entity set.
  Select on donor evidence, **not flower traits or resulting scores**.
- The STRICTEST control holds potential same-genus independent trait donor
  coverage within ±2 percentage points and total donor name pool sizes
  within a factor 0.5–2 of the real assigned native mainland donor pool.
  The looser sensitivity matches only genus donor coverage.
- Give each focal island one score vote; estimate 95% uncertainty by 4,999
  whole-spatial-block bootstrap resamples. Positive means the original
  nearby mainland assignment predicts the heldout profiles better than
  the two support-matched pseudo-source assignments.
- All four trait kernels are scored on the **same exact focal islands and
  source-placebo controls**. Additive gains must agree with the previously
  frozen original distance-K score to numerical precision.

## Strict control: independently rechecked results

| Proper quadratic prediction score | All: 1,007 islands, 119 spatial blocks | Direct: 950 islands, 116 blocks |
|---|---:|---:|
| **Kcolour** marginal-source gain | +0.000991 [0.000315,0.001704] | +0.000596 [−0.000100,0.001258] |
| **Kform** marginal-source gain | +0.002729 [0.000678,0.004405] | +0.001258 [−0.000308,0.002582] |
| Kadd: exactly the marginal mean | +0.001860 [0.000687,0.002877] | +0.000927 [−0.000031,0.001741] |
| **Kjoint: full-profile colour×form** | **+0.002378 [0.000878,0.003753]** | **+0.001097 [−0.000153,0.002178]** |
| Original categorical identity Brier | +0.002011 [0.000788,0.003295] | +0.000811 [−0.000375,0.001908] |

The old **additive** source-placebo scores are reproduced from independent
original checkpoint join and original source assignment calculations
to better than six decimal places: **+0.001859982 versus published
+0.0018600 All**, **+0.000927052 versus +0.0009271 Direct**.
That identity and the original full-profile carrier-count checks
are the replay's primary error-detection gates.

**Main observation:** On strictly donor-support-matched island comparisons,
actual geographical native-mainland source labels provide **a small
additional probability signal that remains visible under the strictly
proper Kjoint scoring rule in All**, despite strong attenuation from
donor coverage controls. Direct has the same sign but its block-bootstrap
interval contains zero; this is **not** two independent confirmations.
The marginal form gain is larger than the marginal colour gain
on the All point estimates.

Across the five disjoint geographic folds, STRICT Kjoint source
gain was positive **5/5 All and 5/5 Direct**, while Kcolour source gain
was positive in only **4/5 All and 4/5 Direct**. Fold means share
trait donor population and genera; a positive sign count alone is
not confirmatory independence.

## Looser support control sensitivity

With only genus-donor-coverage matching, **n=1,137 All / 1,077 Direct**:

| Score | All, block 95% | Direct, block 95% |
|---|---:|---:|
| Kcolour | +0.001227 [0.000583,0.001876] | +0.000803 [0.000034,0.001513] |
| Kform | +0.003650 [0.001303,0.005707] | +0.002785 [0.000108,0.005324] |
| Kadd | +0.002439 [0.001111,0.003642] | +0.001794 [0.000204,0.003319] |
| **Kjoint** | **+0.002945 [0.001366,0.004359]** | **+0.002087 [0.000177,0.003840]** |

The strict additional donor-pool-size match reduces point estimates
and removes unambiguous Direct interval separation, further emphasizing
**evidence availability as a major limitation** rather than claiming
a strong source-specific floral syndromic mechanism.

## New marginal-fixed follow-up (completed, 10 November 2026)

[Same-island, donor-matched, **marginal-fixed** joint source analysis](MARGINAL_FIXED.md) · [Aggregate score/95% interval table](marginal_fixed_placebo_aggregate.csv) · [All five geographic folds](marginal_fixed_five_fold_audit.csv).

**Important falsification:** The positive *unconstrained* Kjoint advantage above is **mostly marginal**. After each true/placebo local-genus full-profile probability distribution is IPF-rescaled so that its **separate colour-set and form-set marginals BOTH exactly match the global-genus benchmark** on the 496 observed categories, the strict All joint K-score gain falls from **+0.002378** to **+0.000214 (95% spatial-block interval +0.000036 to +0.000370)**. Direct falls from **+0.001097** to **+0.000140 (−0.000062 to +0.000331)**. Strict full-profile Brier yields **+0.000205 (0.000022 to 0.000379)** All and **+0.000137 (−0.000087 to +0.000360)** Direct. The local source score's residual conditional coupling signal is thus **very small and Direct-uncertain**; there were no IPF convergence failures. All/Direct remain nested scopes and these post-hoc scores are not evolutionary causal effects.

## Scientific ceiling and next test

**Supported observational result:** Geography-linked native-mainland
donor assignments offer some additional predictive information for
recorded *joint* floral phenotype categories in the All evidence scope
after strict donor opportunity control. The effective sizes are **small**
and Direct's conservative interval includes zero. This extends the
previous additive/soft-trait result to a **full-rank, strictly proper
joint-profile score**, not merely the marginal mean of colour and form.

**New observational qualification:** The marginal-preserving
counterfactual has now been computed: the source assignment retains
a **very small positive All coupling-sensitive score** under a strict
donor-matched reference, but Direct's 95% interval includes zero.
This does NOT identify an evolved flower colour/form syndrome or
causal ecological selection. The residual may still reflect
taxonomy, observed trait reporting, private GIFT donor sampling,
shared historical geography and incomplete source-region coverage.
A row/column probability projection over only the 496 originally
recorded categories does not generate hypothetical colour/form
combinations. We have not measured pollinator efficacy, evolutionary
selection, actual mainland ancestors, dispersal success, biological
absence or independent geological adaptive shrinkage.

There is no justification to turn the current small source-score
improvement into a causal historical mainland-island claim. Report
the All association and Direct uncertainty together, and center the
paper on **genus-linked and geographic transferability of recorded
floral phenotypes**, not adaptive floral simplification.
