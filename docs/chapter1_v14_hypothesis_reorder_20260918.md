# Chapter 1 v14 — reordered H1–H4 analysis

Status: **new reanalysis layer; v13 remains immutable provenance**.

## Why the order changes

The v13 numbering placed the independent GloPL pollen-limitation test before the
plant-side pathway decomposition. The biological argument is easier to inspect in
the opposite order:

1. establish the island syndrome;
2. ask how its floral component is generated;
3. independently test whether the proposed ecological pressure increases with isolation;
4. test whether current pollen limitation is functionally associated with the syndrome traits.

No v13 result lock is rewritten. H3 and H4 reuse the already frozen GloPL and
functional-bridge results; H1 and H2 are refit/reorganized under the definitions below.

## H1 — Global floral/reproductive island syndrome

H1 now has three explicit biological domains.

| domain | atomic responses | positive island-syndrome direction |
|---|---|---|
| reproductive assurance | self-compatibility; selfing mating system; autonomous selfing | more reproductive assurance |
| colour dulling | plain colour | more white / green-brown-inconspicuous relative to vivid colour states |
| accessibility/generalization | generalized form; actinomorphy; shallow/open tube | easier floral access |

The formal H1 test is the **seven-response multivariate Wald test** within each of
the four predeclared geographic replication strata, using the same beta-binomial
logit models and spatial-block cluster-robust covariance as the existing all-data
probability analysis.

The three domain means are descriptive. If one scalar orientation is displayed,
the three domain means receive equal weight. This prevents the 3+1+3 number of
atomic traits from silently determining biological weights.

Plain colour is a colour-composition contrast. It is not a direct measure of
animal-visible conspicuousness, pigment investment or pollinator identity.

## H2 — Why does the floral part change?

H2 contrasts two non-exclusive routes.

### H2a — selfing-syndrome / reproductive-assurance route

Model: selfing_core ~ isolation + area + climate.

selfing_core contains SI/SC, mating system and autonomous-selfing evidence and
deliberately excludes floral attraction/display traits.

### H2b — pollinator-facing floral route

Primary conditional tests:

- plain_colour ~ isolation + selfing_core + area + climate
- generalized_accessible ~ isolation + selfing_core + area + climate

Pollination-syndrome concordance is evaluated from the observed trait states rather
than a weighted guild score. The retained hierarchy is:

1. five raw colour states after selfing_core adjustment;
2. raw colour × raw form/tube joint prevalence;
3. raw architecture conditional on raw colour and selfing_core.

Predeclared bee/butterfly/bird associations are labels for particular raw form/tube
combinations used to interpret concordance. The historical weighted
large_bee_like/butterfly_like/bird_like scores do not define v14 H2 evidence.

If the isolation coefficient persists after conditioning on selfing_core, the
floral response is **not reducible to the measured selfing core**. This is not
causal mediation and does not prove that pollinator decline directly selected the
trait.

## H3 — Independent pollen-limitation pressure

H3 is the former v13 H2, renumbered without changing its statistical model.

The frozen GloPL analysis tests whether experimental pollen limitation increases
with geographic separation using publication-normalized weights and
publication-cluster-robust covariance.

This establishes an independent ecological-pressure pattern. It does not establish
that pollen limitation mediated H1 or H2.

## H4 — Functional bridge

H4 is a **post-hoc functional bridge**. Prospective validation is not part of the H4 result.

The v13 exact-species atomic-trait comparison remains the original post-hoc discovery
layer. Autonomous selfing remains its strongest individual bridge.

v14 adds a more direct post-hoc H2-to-H4 test by using the **literal Direct-only
species-level H2 scores** already used by Chapter 1, with no score reconstruction:

- exact `selfing_core`: `beta=-0.29706`, two-sided `p=0.00411`;
- exact `generalized_accessible`: `beta=-0.29601`, two-sided `p=0.0222`.

Both H2 scores that increase with island isolation therefore point toward lower current
GloPL pollen limitation in the primary post-hoc bridge. The accessibility result also
remains supported in the supplemental-only subset and retains the predicted negative
direction in the no-zero sensitivity.

A separate atomic-family reconstruction is retained only as sensitivity because its
species-level scores are concordant with, but not numerically identical to, the literal
H2 scores. The reconstruction produces the same biological direction, including an
equal-weight accessibility definition that matches the prospectively frozen H4b
prediction.

No global primary colour bridge is imposed because the H1/H2 colour response is
context dependent and does not define one globally directional functional prediction.

## Supplementary validation audit — not part of H4

PR #236 prospectively froze a post-2015 wild-plant route and an independent PolLimCrop
route. Both stopped at outcome-blind support gates before pollen-limitation outcomes
were opened. These audits are retained only to document that confirmatory evidence is
still unavailable. They do not enter, upgrade, or downgrade H4.

## Analysis flow

    H1  WHAT?
    isolation -> reproductive assurance + plain colour + accessible/generalized form
                             |
                             v
    H2  HOW DOES THE FLORAL SHIFT ARISE?
              +--------------+----------------+
              |                               |
     selfing-syndrome route       selfing-adjusted floral route
              |                               |
              +--------------+----------------+
                             |
                             v
    H3  IS THE PROPOSED PRESSURE PRESENT?
              isolation -> pollen limitation
                             |
                             v
    H4  DO THE TWO SIDES FUNCTIONALLY CONNECT?
       literal H2 score -> current pollen limitation
                     +
       atomic reconstruction sensitivity

## Claim ceiling

The v14 reordering does not identify historical trait evolution, mediation by pollen
limitation, temporal pollinator decline, direct selection by a named pollinator, or
realized pollinator identity from floral syndrome scores.
