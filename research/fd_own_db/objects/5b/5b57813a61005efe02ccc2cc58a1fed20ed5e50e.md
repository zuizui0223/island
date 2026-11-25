# Why accurate island floral-category prediction need not recover functional diversity

A general algebraic result, independently of whether an island botanical source
dataset can identify colonization, pollinator-mediated selection or phenotypic
adaptation. Code and reproducible synthetic tests:
- src/fd/nominal_diversity_identity.py
- tests/test_nominal_diversity_identity.py

## One empirical island plant cohort: exact category-score decomposition

Let q be the empirical K-category frequency of TRUE same-species Direct floral
traits among n botanical plant species on one scored GIFT island; p is the
category probability predicted by a no-focal-genus-oracle model. Both
distributions use exactly the SAME category ontology and target species.

The multiclass Brier score averaged over individual observed one-hot plant
categories is, by exact expansion:

\[
\overline{\operatorname{Brier}}(p)
=1+\sum_kp_k^2-2\sum_kp_kq_k
=\underbrace{1-\sum_kq_k^2}_{Q_{\rm plugin}(q)}
+\underbrace{\sum_k(p_k-q_k)^2}_{\text{categorical prediction error}}.
\]

Therefore a matched Brier improvement has a straightforward meaning: it
reduces squared categorical frequency error against that empirical island
composition q. It **does not** directly minimize the error in the nominal
categorical Rao/Gini–Simpson diversity quantity.

Define Q(p)=1−Σ_kp_k², the Rao quadratic form with zero distance within
a category and unit distance across different categories. Put e=p−q.
A second exact identity is:

\[
Q(p)-Q(q)=-2\sum_kq_ke_k-\sum_ke_k^2.
\]

Brier controls only the magnitude of the final squared-error term, whereas
the first, SIGNED projection of the category forecast error onto observed
frequencies also matters. In particular, an overprediction of the commonest
floral categories can reduce Q despite making the entire predicted categorical
frequency vector closer to the true one.

## Explicit numeric counterexample, same ten observed plant species

There are ten Direct species on one imaginary island: eight in category A
and two in category B. q=(0.8,0.2). Compare the blind categorical prior
b=(0.6,0.4) and a model p=(0.95,0.05). Both are positive probability
simplexes, and the model is closer to the actual observed category
frequencies in squared error:

| Exactly same ten-plant source cohort | Observed q | Blind source b | Improved category model p |
|---|---:|---:|---:|
| Category A probability | 0.8 | 0.6 | 0.95 |
| Category B probability | 0.2 | 0.4 | 0.05 |
| Empirical mean proper Brier | n/a | **0.400** | **0.365** (better by 0.035) |
| Nominal predicted Q | 0.320 (plugin) | 0.480 | 0.095 |
| Absolute error versus observed plugin Q=0.320 | 0 | **0.160** | **0.225** (worse by 0.065) |
| Absolute error versus iid-unbiased Q=0.320×10/9≈0.35556 | 0 | **0.12444** | **0.26056** (also worse) |

Thus a model can have reproducibly positive proper-score gain and negative
diversity-magnitude absolute-error gain on IDENTICAL plants. This is not a
paradox or evidence that the biological floral syndrome signal changed.

The sample-size adjustment Q_unbiased=n/(n−1) Q_plugin is an unbiased
population nominal Simpson estimator only under an appropriate iid
categorical sampling model. GIFT incomplete recorded species-Direct floral
traits, checklist effort and genus-filtered source eligibility need not
satisfy iid sampling. The n/(n−1) factor cannot repair taxon-dependent
observation bias.

## Exact extension to a real flower-trait distance matrix

For a SHARED fixed symmetric nonnegative zero-diagonal category-to-category
dissimilarity matrix D (a valid candidate *if* botanical flower morphology
distances are independently measured and aligned), define general Rao
Q_D(p)=p^T D p. For e=p−q:

\[
Q_D(p)-Q_D(q)=2e^T Dq+e^T De.
\]

This identity is implemented in rao_distance_decomposition; it does not
presume D is positive-semidefinite, and the quadratic term e^T D e is NOT a
squared Euclidean frequency error in general. A generalized Gower/mixed-trait
functional distance, including ordinal flowering phenotypes and reproductive
mating state, would require empirically justified pairwise distances,
taxonomic source coverage and an explicit missing-state treatment. A 17-D KDE
hypervolume, categorical nominal Rao Q and distance-weighted Rao Q are three
different estimands and should never be silently substituted.

## What this says about the original GIFT results

Previously SOURCE VERIFIED original GIFT 4,505 botanical island floras:
- all 13 true Direct category endpoints showed some regional predictive gains,
  with 9/39 all-recorded macroregion×trait panels clearing a strong LOCAL
  region-blind category Brier/log and double group interval gate;
- of those, only tropical floral SYMMETRY and Southern extratropical flower
  COLOUR cleared the posthoc three-comparator stricter category gate;
- original two Direct phenotype endpoints (form and genuine same-species
  colour×form) had **0/6** unbiased nominal-Q all-recorded wins over strong
  local-regional Q reference, despite some categorical successes.

The last check already establishes that nominal diversity prediction is a
separate task. The more comprehensive original 13-endpoint nominal-Q source
runner is DRAFT island PR #337 and its corrected source run is still
being checked separately; the source preliminary logs are NOT the
independently audited final biological result.

The equations above explain why Brier and Rao Q gains need not match as a
mathematical matter; they do NOT identify which of the projection,
incomplete trait documentation, or geographic taxonomic assembly actually
causes particular original GIFT region×trait failures. Claiming a biological
mechanism requires independently audited original source-category vectors or
an explicitly predeclared source-safe aggregate decomposition and additional
observations.

## Operational manuscript gate

Keep: original full GIFT 4,505 botanical parent, all original 13 trait
endpoints (do not require full trait completeness across all taxa), exact
same-accepted-species Direct pair traits, no focal genus oracle, source name
and physical/first-source embargo, dual proper losses and physical/source
group uncertainty, original native-status support failures, unbiased/plugin
Q distinction and honest negative results. Require independent local
pollinator/fitness or ancestry data before claiming mechanism.

This is a general identity with verified synthetic counterexample tests,
not a novel measured ecological effect, empirical Gower FD analysis or
independent biological replication.
