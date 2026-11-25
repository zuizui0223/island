# fd manuscript decision draft — geographical transfer of floral trait similarity

**Working version:** v0.5, 10 November 2026
**Status:** independent source-donor placebo and new **marginal-preserving
conditional colour–form association prediction** have now been calculated
on the same strict 1,007 All / 950 Direct island cohorts. The residual
association gain is small and Direct's 95% interval crosses zero.
Interaction-specific evolutionary selection is not identified.
The intent is one genuinely macroecological, source-backed paper, not many
weakly connected null-model results.

## Recommended main research question

**To what extent can the recorded floral phenotypes of island floras be
predicted from plant lineage composition and the geographic distribution
of potential native-mainland donors, and does this transfer apply to the
separate floral traits or to their joint combinations?**

This is a **prediction/assembly question**, not an adaptive-selection
claim. It has a potentially transferable mechanism (regional lineage
availability) without pretending that GIFT native-source assignments
are known historical ancestors.

### Working title

*Plant lineage composition dominates the geographic predictability of
island floral phenotypes*

**Alternative conditional title (ONLY if a completed donor-matched
product-kernel source audit succeeds):**

*Mainland lineage supply shapes the transfer of floral trait
similarity to island floras*

Do NOT use "joint floral syndrome", "adaptive floral compression",
"pollinator syndrome evolution", or "source-driven convergence"
in the title on the current evidence.

## Three core findings that ARE supported

1. **Most recorded floral prediction skill is lineage-associated.**
   Five geographically disjoint source-fold tests use 1,533 All and
   1,443 nested Direct eligible islands. Against a no-genus global
   donor baseline, known global genus reduces heldout exact colour/form
   profile log loss by **0.641 (block bootstrap 95% range 0.587–0.691)
   All** and **0.360 (0.297–0.413) Direct**. This is conditional on
   KNOWN genus of each evaluated island name and does not predict
   which genera successfully colonize an island.

2. **Geographically selected native-mainland donors add a small,
   similarity-aware predictive increment, but do not improve exact
   flower-profile probabilities consistently.** Global genus ->
   locally assigned genus changes exact-profile log loss by
   **+0.00011 (−0.00451,+0.00379) All** and
   **+0.00019 (−0.00606,+0.00517) Direct** (positive = better).
   Yet it reduces proper *additive original-D* K-score by
   **+0.00829 (0.00633,0.01004) All** and
   **+0.00948 (0.00655,0.01191) Direct**. The full-category
   strictly proper Brier increment is **+0.00221
   (0.00095,0.00335) All**, but Direct **+0.00178
   (−0.00022,0.00349)** remains uncertain. The additive K-score is
   NOT evidence of predicting colour/form coupling.

3. **Source geography must be distinguished from source-donor
   data availability.** The first true-source vs random-source
   comparison was **+0.00396 All/+0.00352 Direct** in proper
   additive-K score, but the real assignment had **+5.37/+3.52
   percentage points** more complete independent same-genus donors.
   Against two alternative source assignments matched on same-genus
   donor coverage (±2 pp) AND total donor name pool (0.5–2x),
   the same comparison is only **+0.001860
   (0.000687,0.002877) All** across 1,007 islands/119 blocks, and
   **+0.000927 (−0.000031,0.001741) Direct** across 950/116.
   Scope overlap and post-hoc matching preclude causal percentages
   explained. Within-source independent local genus donor coverage
   itself is sparse: **~24.3% All/~19.6% Direct**.

**Taxon duplication sensitivity:** In the All 16,386 distinct names
(2,262 genera) and Direct 8,244 names (2,214 genera), scoring each
different taxon name only once and resampling entire genera preserves
positive additive-K source gains (**+0.01354 [0.01027,0.01716]
All; +0.00905 [0.00696,0.01109] Direct**). But note that this
robustness is for local genus vs GLOBAL genus on the all-eligible
population, not the stricter 1,007/950 source-matched-placebo cohorts.

## Falsifications that define the claim ceiling

- The 17-dimensional Gaussian KDE "volume" tracks recorded distinct
  profile counts near-deterministically. KDE goes to Methods/Appendix
  only, not the main ecological estimand.
- All/Direct nested scope != external replication.
- Original 37 geological-origin pairs form only five conservative
  geographic connected components; the strict source/block-disjoint
  rematch retains ten geographically verified pairs. Multiple-testing-
  corrected geological-origin functional shrinkage **is not confirmed**.
- Uniform profile reference produces exaggerated negative local
  packing SES, which weakens when other-island profile incidence,
  neighboring islands and fixed absolute geographic radii are used.
- Source priors strongly depend on independent donor genus evidence
  and smoothing. Full all-island simultaneous taxon-name holdout
  retains only 228/9,017 All and 172/6,562 Direct known-profile donors;
  geography-fold crossfitting is necessary for predictive evaluation.
- An additive colour/form Jaccard K-score tells us nothing about
  colour/form *pairing beyond the separate marginals*. The frozen
  496-profile rank audit shows **Kcolour=31, Kform=46, Kadd=76,
  Kjoint=496**. Thus the new product kernel is **strictly proper**
  on the full recorded profile reference, not just a coarse marginal
  summary. New heldout donor-matched replay supports a **small
  joint-sensitive source prediction increment** in All
  **+0.002378 [0.000878,0.003753]** (1,007 islands, 119 blocks)
  and uncertain Direct **+0.001097 [−0.000153,0.002178]**
  (950 islands, 116 blocks). Mathematical full rank and positive
  prediction increment do NOT uniquely identify a colour×form
  evolutionary interaction, adaptation or ancestral donor origin.

## New decisive falsification: force colour and form marginals to agree

A strictly proper **joint** product kernel alone cannot distinguish
better colour/form marginals from prediction of their association.
We therefore constrain the competing actual and alternative native-mainland
source full-profile distributions to **exactly the same 31 recorded
colour-state-set and 46 form-state-set marginals**, equal to the heldout
**global-genus** predictor. Iterative proportional fitting (IPF)
preserves the source model's association odds on the structurally
supported original 496 colour/form profiles, without inventing
Cartesian-product flowers. All empirical IPF contexts converged
(tolerance 2×10⁻⁷; zero failures).

Using the **same strict, outcome-blind, donor-coverage-and-size matched
source assignment controls** and exactly the same 1,007 All/950
Direct islands:

| Full-profile proper score | All | Direct |
|---|---:|---:|
| Joint product K before fixing marginals | **+0.002378** [0.000878,0.003753] | +0.001097 [−0.000153,0.002178] |
| Joint product K with BOTH marginals fixed | **+0.000214** [0.000036,0.000370] | +0.000140 [−0.000062,0.000331] |
| Full-profile Brier with BOTH marginals fixed | **+0.000205** [0.000022,0.000379] | +0.000137 [−0.000087,0.000360] |

The *measured joint score advantage* is attenuated ~**91% All /
87% Direct** when the colour and form margins are equated.
This is a **scoring-rule/reference-conditioning sensitivity**, not a
causal partition of adaptation or a "percentage of flower evolution
explained". All five geographic folds retain very small positive
joint-association point estimates in both nested evidence scopes,
but Direct's block-bootstrap interval crosses zero.

**Conclusion:** Geographic source prediction primarily changes the
*separate* colour/form trait probabilities. The remaining
colour–form association signal is weak, and is not enough to
claim a globally generalizable evolved joint floral syndrome.
A small All conditional source signal can be described as a
post-hoc **recorded phenotype association**, subject to donor
sampling, shared taxonomy, detection and unknown historical
colonization.

[Full empirical marginal-fixed audit](../results/20261110_joint_flower_source/MARGINAL_FIXED.md) ·
[Aggregate scores and intervals](../results/20261110_joint_flower_source/marginal_fixed_placebo_aggregate.csv) ·
[Reproducible IPF implementation](../src/fd/marginal_fixed_coupling.py).

## Stronger inference audit: what remains uncertain

The marginal-fixed joint-source signal above was reevaluated using the
five entire cross-fitting fold means as coarse sign units rather than
119/116 resampled spatial blocks. The exact 32-configuration two-sided
stress tail is **0.0625** for the All/Direct joint-source K scores
(5/5 positive folds), and **0.1875** for Direct marginal-fixed
Brier (4/5 positive folds). This is not a valid independent ecological
randomization test because the training donor pools and botanical
lineages overlap between folds. It demonstrates that apparently
narrow conditional spatial-block bootstrap intervals do **not**
alone establish generality across truly independent histories.

Leave-one-test-fold-out means stay positive
(**All joint +0.000183..+0.000249; Direct +0.000101..+0.000175**).
This supports stability of the *direction* in this observational
cross-validation design, not a geographically universal effect.

The question whether flower form is more geographically predictable
than flower colour also requires a direct paired contrast. The
matched-source form-minus-colour gain is **+0.001738** All
(95% block interval **−0.000221..+0.003288**) and **+0.000662**
Direct (**−0.000808..+0.001967**). Both intervals include zero;
neither stronger form predictability nor more intense ecological
selection on floral architecture is supported.

Finally, the IPF conditional coupling result uses donor shrinkage
**8**. Variation across shrinkage 2/8/20 was previously important for
exact-category source prediction, and the **marginal-fixed coupling
sensitivity to alternative shrinkage has not been run**. Therefore
the tiny All residual cannot be presented as a robust independent
colour×form evolution discovery.

**One-sentence conclusion locked for the current evidence:**
*Geographically plausible native-mainland donors weakly improve
heldout predictions of recorded floral phenotypes beyond plant genus,
largely via separate colour and form distributions; their additional
conditional colour–form pairing signal is small, observational and
uncertain across source-evidence scopes.*

[Coarse five-fold and paired marginal audit](../results/20261110_joint_flower_source/MARGINAL_FIXED.md).

## Distribution of residual association: a minority-driven mean, not universal transfer

A further *island-level* heterogeneity audit on the SAME 1,007 All /
950 Direct tightly donor-matched source-placebo targets finds
marginal-fixed joint-K gains above zero in **60.8% All / 54.4%
Direct**. The median gains are **+0.000104 / +0.000049**, both
smaller than their respective means. The **10%-trimmed** average
is **+0.000165 (whole-spatial-block 95% interval +0.000012 to
+0.000313) All**, but **+0.000113 (−0.000079 to +0.000305)
Direct**. The top 10% largest positive-gain islands account for
about **48% of all POSITIVE** gains in each scope; because negative
gains offset some islands, they account for about **82% of the
NET SIGNED All gain**.

The All effect is not merely a single extreme outlier, but is
heterogeneous; nested Direct and coarse five-fold tests do not
give strong evidence for ubiquitous positive association gain.
The paper should say **a weak, geographically uneven conditional
predictive signal**, not universal flower colour/form coupling,
synchronous evolution, or a globally consistent ecological driver.

[Private-data-safe heterogeneity audit](../results/20261110_joint_flower_source/coupling_heterogeneity_aggregate.csv) ·
[reproducible code](../src/fd/coupling_heterogeneity.py).

## Main-text figure structure (maximum four figures)

**Figure 1: Geography, available lineage opportunity, and analysis
population.** Map all eligible islands (All 1,533; Direct nested 1,443),
five nonoverlapping geographic training/test blocks, distributions of
candidate native-mainland source region donor coverage. Show
source-donor data coverage instead of raw GIFT "millions of data"
as an independence claim. Aggregate or simplify where licence restricts.

**Figure 2: Predictive trait information has a strong genus component.**
Show island-weighted heldout exact log-loss and proper source K/Brier
improvement under no-genus, known-global-genus and local-source-genus.
Clearly plot improving LOSS, confidence intervals clustered by
spatial block, and one-name-one-vote genus bootstrap as a second row.
Explain that exact category identity does not improve meaningfully
with local source.

**Figure 3: Negative controls remove much of the source-assignment
advantage.** One axis: random same-region pseudo-sources;
another: support-caliper-matched;
third: support + donor-pool-size-caliper. Distinguish changing matched
sample n, spatial block count, and All/Direct separate but nested
estimates. Show actual minus placebo proper additive-K score on all
comparisons; publish Brier companion. Do not call relative change a
proportion of causal effect explained.

**Figure 4: Where does geographic source predictive skill reside?**
Plot colour-only, form-only, additive, unrestricted joint and
**marginal-fixed joint** source-placebo gains with spatial-block
bootstrap interval bars on the SAME 1,007 All/950 Direct source-matched
island populations. Highlight the ~91%/87% joint score attenuation
after IPF, the tiny All positive conditional association,
and the Direct interval crossing zero. A toy 2-colour×2-form
case shows why additive similarity can miss real pairing variation.
Do not equate a conditional prediction signal with trait
co-evolution or a functional pollinator syndrome.

**Supplement:** 17D KDE richness regression; N0/N1/N1b/N1c/N1d
support screens; G1 fixed margins; old 37-pair spatial component
audits; uncertainty of GIFT donor completeness; genus-classification
and native-only sensitivity; source shrinkage 2/8/20.

## Draft results prose (not for submission before final provenance review)

*The categorical colour and floral-form repertoires of recorded island
plants showed substantial geographic predictability, but the source
of this predictability was not a simple geological island syndrome.
Under whole-geographic-block test folds with all focal taxa held out
from native-mainland floral donors, information on a plant's genus
substantially improved exact-profile predictions compared with a
global trait prior. Conditioning additionally on nearby mainland
source-region donors improved similarity-aware predictions of
floral phenotypes, yet gave negligible mean improvement in the exact
profile log-likelihood. Distinct taxon-name equal-weight sensitivity
retained the similarity score improvement. However, after matching
alternative mainland source assignments on independently observed
same-genus donor coverage and overall donor pool size, the
incremental advantage of the true geographic assignment was much
smaller and remained uncertain in the stricter Direct data subset.
Matching predicted colour and form marginals across the competing
mainland-source models reduced the already modest joint-profile
source gain by approximately an order of magnitude, leaving a small
All conditional-pairing association and an uncertain Direct analogue.
These results indicate transferable genus-associated floral trait
structure and a weak, observation-sensitive geographic source signal
that is primarily marginal-trait-linked, not demonstrated floral
co-evolution, pollinator adaptation or trait simplification.*

## Decision rules before drafting the final manuscript

1. **Stop adding KDE/SES nulls.** They have served their falsification
   purpose. They do not support independent geology-driven adaptation.
2. **Joint predictive replay complete:** all 496 original profile
   identities and source checkpoint SHA hashes were independently
   verified; the old additive source gains were reproduced to six
   decimal places. Report strict All positive and nested Direct
   uncertainty as a single correlated evidence stream.
3. **Marginal-preserving comparison is now computed.** It shrinks
   the joint source-placebo advantage to +0.000214 All and +0.000140
   Direct; Direct is uncertain. This identifies a conditional
   *predictive association* given heldout known genus, NOT an evolved
   colour–form syndrome. Phylogenetic and detection safeguards
   remain an unresolved gate.
4. **Calibrate inference under repeated species and shared genera.**
   The earlier taxon-balanced improvement is not a proper substitute
   for name-balanced STRICT source-placebo comparisons.
5. **Retain descriptive, not causal language**, even after a positive
   product-kernel score. True source, arrival, selection, pollinator
   biology and trait observation are not separately identified.
6. End the paper with a broadly applicable question: how source
   opportunity and lineage structure constrain phenotypic repertoire
   *predictability*, including failed predictions, rather than
   forcing an adaptive island-syndrome explanation.



## Newly completed strictly proper joint-score result

The same source-name-heldout predictive experiment has now been evaluated
with a **full-rank, strictly proper joint colour-state-set ×
flower-form-state-set Jaccard product kernel** on **the same
geographically disjoint test folds** and **the same real-vs-placebo
source candidate donor opportunity controls**. Independent local
reconstruction from the two pinned trait checkpoints reproduces
all 496 categorical reference IDs and their source carrier counts.

In the strictly coverage- and total-donor-size-matched group of
**1,007 All islands from 119 spatial blocks**, the real nearest
native-mainland source assignment improves the complete joint
colour/form K-score over the mean of two geographic placebo
assignments by **+0.002378 [0.000878,0.003753]**. The nested
**950 Direct islands in 116 blocks** show **+0.001097
[−0.000153,0.002178]**. The additive kernel improves by
+0.001860/+0.000927 respectively, **exactly the average
of the colour and form marginal improvements**, and matches the
previously published source placebo to six decimals.
The separate form marginal gain is larger than colour in All,
but this is a descriptive difference between highly
correlated proper scores.

The key ecological interpretation is a small conditional signal
of geographic mainland lineage source information in probability
distributions over intact recorded flower profiles, not a
large adaptive island syndrome. The source-matched comparison
is still post-hoc and constrained by incomplete mainland donors,
repeated taxa across islands, known test genus identities,
and nonrandom observational availability.

[Empirical joint-score audit](../results/20261110_joint_flower_source/README.md)
[Strict matched-score CSV](../results/20261110_joint_flower_source/coverage_matched_trait_kernels.csv)

## Files and provenance

- [Heldout source/lineage proper-score results](../results/20261110_n2_proper_scoring/README.md)
- [Support- and donor-matched geographic-placebo results](../results/20261110_n2_source_placebo/README.md)
- [Marginal versus joint kernel contract](FLORAL_TRAIT_COUPLING.md)
- [Private N2 reproducibility and scope contract](N2_SOURCE_DISTANCE.md)
- [Pinned original source PR321](https://github.com/zuizui0223/island/pull/321)

All data rights remain unchanged; species/flower assignments and native
GIFT source records must never be copied into public fd outputs.
