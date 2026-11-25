# fd manuscript decision draft — geographical transfer of floral trait similarity

**Working version:** v0.4, 10 November 2026
**Status:** empirical findings through **independently replayed joint-colour/form
source-donor matched controls** verified in All; Direct direction agrees
but its conservative whole-spatial-block uncertainty includes zero.
Interaction-specific evolution is not identified.
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

**Figure 4: Marginal versus joint floral trait source transfer
(EMPIRICALLY FILLED).** On the same support-matched source placebo
cohorts, show colour-only, form-only, additive identity and full-rank
joint-product K-score improvements with spatial-block interval bars.
For All 1,007 islands: colour **+0.000991**,
form **+0.002729**, additive **+0.001860**, joint
**+0.002378 (95% [0.000878,0.003753])**.
For nested Direct 950: colour **+0.000596**,
form **+0.001258**, additive **+0.000927**, joint
**+0.001097 (95% [−0.000153,0.002178])**.
This joint result is a *strictly proper* distributional similarity
score in the 496-category reference, but it can improve via better
marginal predictions. The accompanying toy coupling case illustrates
the algebraic gain of testing full-category arrangements. There is
NO purified interaction mechanism or within-lineage adaptation test.

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
These results indicate transferable genus-associated floral trait
structure and a weak, observation-sensitive geographic source
signal, rather than a demonstrated continent-to-island trajectory
of floral functional convergence or adaptive trait loss.*

## Decision rules before drafting the final manuscript

1. **Stop adding KDE/SES nulls.** They have served their falsification
   purpose. They do not support independent geology-driven adaptation.
2. **Joint predictive replay complete:** all 496 original profile
   identities and source checkpoint SHA hashes were independently
   verified; the old additive source gains were reproduced to six
   decimal places. Report strict All positive and nested Direct
   uncertainty as a single correlated evidence stream.
3. **Do not infer a pure colour–form association mechanism** from
   joint-product improvement alone. To identify interaction-specific
   assembly, a valid marginal-preserving trait coupling counterfactual
   and independent phylogenetic/detection safeguards would be needed.
   They are beyond the current paper's confirmed result.
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
