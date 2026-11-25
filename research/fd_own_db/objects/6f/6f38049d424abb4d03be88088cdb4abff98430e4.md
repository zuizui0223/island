# Additive Gower functional diversity is blind to true multi-trait floral syndromes

Mathematical methods note for the GIFT global island functional-diversity
project. This is an exact identity, **not** a new measured historical
pollinator shift, plant within-species evolution, or mainland ancestor trait
sample. Source implementation and synthetic tests:
src/fd/paired_nominal_trait_diversity.py and
tests/test_paired_nominal_trait_diversity.py.

## Two genuinely paired species-Direct floral attributes

Suppose each plant species has a genuinely recorded SAME-species two-axis
state (A,B), e.g. colour×form or floral symmetry×autonomous-selfing capacity,
and the observed plant distribution on a particular island is p(a,b).
Do NOT invent pair states by multiplying independent observations from
different botanical species.

The usual additive nominal-Gower distance for two completely documented
axes is

\[
D((a,b),(a',b'))
=w_A\mathbf 1[a\ne a']+w_B\mathbf 1[b\ne b'],\quad
w_A+w_B=1.
\]

The usual double-sum Rao quadratic entropy of that distance is EXACTLY

\[
Q_D(p)=\sum_{ab,a'b'}p_{ab}p_{a'b'}D((a,b),(a',b'))
=w_A\bigl(1-\sum_a p_A(a)^2\bigr)
+w_B\bigl(1-\sum_b p_B(b)^2\bigr).
\]

**It depends only on the two marginal trait frequencies.** It cannot see
whether a particular flower colour is preferentially paired with a
particular flower form. Two island floras with exactly the same marginal
colour/form frequencies but fundamentally different species-level trait
pairing have the SAME additive nominal-Gower Rao Q. This is independent of
geography, botanical genus modelling and trait documentation effort.

## Concrete four-category contrast with equal marginals

Let colour be c0/c1 and form f0/f1. Consider four genuine possible same-
species category states, ordered (c0,f0),(c0,f1),(c1,f0),(c1,f1).

| Different floristic observed joint frequencies | 4 source categories | Marginal colour Q | Marginal form Q | Additive Gower Rao Q (equal weights) | Joint-category nominal Rao Q |
| --- | --- | ---: | ---: | ---: | ---: |
| Independent colour and form | (.25,.25,.25,.25) | .5 | .5 | **.5** | **.75** |
| Perfect positive coupling | (.5,0,0,.5) | .5 | .5 | **.5** | **.5** |
| Perfect negative coupling | (0,.5,.5,0) | .5 | .5 | **.5** | **.5** |

The additive-Gower score is always **0.5**, and therefore cannot identify
the radically different observed category pairings. Joint nominal Rao Q
detects that the number/evenness of *pair states* changes, but by itself
does not reveal WHICH coupling pattern occurs or WHY it exists.

A signed descriptive joint-state diversity excess against a same-marginal
independence counterfactual is

\[
S_{\rm joint}
=Q_{\rm joint}(p)-Q_{\rm joint}(p_Ap_B)
=\left(\sum_a p_A(a)^2\right)
 \left(\sum_b p_B(b)^2\right)-\sum_{ab}p_{ab}^2.
\]

The value is zero under a fully independent joint category frequency,
negative for common matched-pair coupling, and can be positive for some
different off-diagonal distributions; it is NOT a unique association
statistic or evidence of island directional selection. For full joint
behavioural interaction, test the entire categorical P(A,B) against the
product of its own marginals and, crucially, a genus-conditioned null
SUM_g P(g)P(A|g)P(B|g) to avoid labelling between-genus assembly as an
independent within-genus syndrome effect.

## Why this matters for the real source GIFT analysis

The source-original 4,505 GIFT island flora cross-regional tests:
- 13 high/medium genuine species-Direct floral/reproductive endpoints
  showed only 9/39 all-recorded source region×endpoint cells better
  than a strong local-region genus-blind baseline, with only 2/39
  additionally beating two harder alternative source predictors;
- six genuine flower architecture (form/symmetry/tube) ×
  self-incompatibility/autonomous-selfing Direct pairs, evaluated in
  [original source run #44544366942](https://github.com/zuizui0223/island/actions/runs/44544366942),
  had **0/11** independently source-group-qualified original all-recorded
  joint-predictor gains over a genus-conditioned independence null
  (7 losses, 4 ambiguous), and 19/36 full status×region cells unsupported.

This exact mathematical identity explains why an additive mixed-trait
distance functional diversity metric is structurally incapable of
**separately identifying multi-trait coupling**. It is perfectly valid
for marginal functional differentiation under its distance definition.
Use BOTH categories/association-sensitive tests (actual same-species
pair states, joint source predictive proper score and genus-conditioned
independence) AND a clearly defined marginal distance-based functional
diversity readout if the research question asks about **syndrome-like
combinations**. Do not claim any of these detect actual animal guild
replacement, pollination fitness, local evolutionary change or distinct
mainland ancestors without independent observations.

## Scope and limits

The code demonstrates discrete NOMINAL two-axis Gower, not a full
continuous morphology trait matrix. Ordinal tube-depth class distances,
measured geometric distances, UV/colour distances, flower size and
self-incompatibility require explicitly justified, consistently scaled
empirical dissimilarities and correct missing-state treatment; do not
invent biological distances. Incomplete GIFT Direct flora records are
documentation selected, not iid plant species observations, so these
comparisons need source support, native-only holdout and original
spatial/source-cluster uncertainty assessments.

The exact two-axis identity is general, not a 4,505-island inferential
effect. It does not raise the number of independent islands, phylogenetic
replicates or documented pollinator transitions in any study.
