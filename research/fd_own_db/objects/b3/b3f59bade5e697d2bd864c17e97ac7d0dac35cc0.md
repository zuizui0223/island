# fd manuscript decision draft — geographical transfer of floral trait similarity

**Working version:** v0.3, 10 November 2026
**Status:** empirical findings through source-donor matched controls verified;
joint-colour/form prediction check is **code only**, with **NO empirical joint result yet**.
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
  summary. Its real source-placebo prediction result is **NOT YET
  CALCULATED**. Mathematical identifiability != ecological evidence.

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

**Figure 4: Marginal versus joint floral trait transfer (UNFILLED
until real validation).** Colour-only K-score, form-only K-score,
additive identity and colour×form product joint K-score, measured on
IDENTICAL donor-matched island/control cohorts and heldout test names.
A toy coupling case illustrates why additive-K misses composition.
If joint source score does not survive, this figure is a falsification,
not a missing success; the narrative remains useful. If it does,
describe it as joint-sensitive predictive association, not
interacting trait selection.

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
2. **Run and inspect the NEW colour/form marginal and product-joint
   source-placebo scores on exactly the same geographic folds and
   support-matched strict cohorts.** Check that additive source gains
   reproduce the previously published +0.001860/+0.000927.
   A discrepancy is an implementation failure, not a new effect.
3. **Do not describe a joint colour–form signal from the existing
   additive-K results.** Joint transfer stays open until the product
   kernel check and its source assignment negative controls.
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

## Files and provenance

- [Heldout source/lineage proper-score results](../results/20261110_n2_proper_scoring/README.md)
- [Support- and donor-matched geographic-placebo results](../results/20261110_n2_source_placebo/README.md)
- [Marginal versus joint kernel contract](FLORAL_TRAIT_COUPLING.md)
- [Private N2 reproducibility and scope contract](N2_SOURCE_DISTANCE.md)
- [Pinned original source PR321](https://github.com/zuizui0223/island/pull/321)

All data rights remain unchanged; species/flower assignments and native
GIFT source records must never be copied into public fd outputs.
