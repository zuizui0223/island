# Chapter 1 v13 raw flower-colour audit

Status: exploratory secondary audit. This does **not** change the canonical v13 H1–H4 submission surface unless the frozen reanalysis is reproduced and explicitly promoted later.

## Why this audit exists

The final Chapter 1 database contains a large flower-colour axis, but the v13 primary plant response deliberately uses reproductive assurance and floral accessibility/generalization. Collapsing colour to `plain_colour`, or embedding it inside named pollination-syndrome scores, discards biologically useful information.

This audit therefore returns to the reported colour states themselves:

- `white`
- `red_pink`
- `yellow_orange`
- `blue_purple`
- `green_brown_inconspicuous`

No colour is treated as a pollinator identity. Raw colours can be compared later with predeclared pollination-associated architectures, but a red, yellow/orange or blue/purple flower is not by itself evidence for a butterfly, bird or bee visitor.

## Estimand

For each island and evidence scope, the response for colour `c` is:

> share of focal-colour-resolved plant species explicitly reported to contain colour `c`.

Species with multiple reported colours remain multistate and contribute a success to every focal colour they explicitly contain. Species whose only colour state is non-focal (for example `other_described`) do not enter the focal-colour denominator.

This is intentionally not a mutually exclusive multinomial recoding.

## Models

Each colour is fitted separately under the existing Chapter 1 grouped-binomial + spatial-block cluster-robust framework.

### Unconditional colour response

`colour_c ~ isolation + island area + climate PC1-4`

### Conditional on reproductive assurance

`colour_c ~ isolation + selfing_core + island area + climate PC1-4`

The conditional model asks whether a colour-composition shift remains associated with isolation after accounting for the measured reproductive selfing core.

Persistence of the isolation coefficient is compatible with a pollinator-facing floral response that is not reducible to measured selfing. It is **not causal mediation** and does not establish an attraction mechanism by itself.

## Joint colour-composition test

Within each geographic replication stratum, all supported raw-colour isolation coefficients are also tested jointly with a spatial-block cluster-robust Wald test.

This asks whether the raw flower-colour vector changes with isolation even when no single colour alone is treated as the syndrome.

## Evidence scopes

- `all`: all analysis-eligible colour evidence in the frozen final species-axis snapshot;
- `direct`: High/Medium colour evidence only.

The primary floristic population for this v13 audit is `all_observed`, matching the current global-only paper surface. Native/status-resolved versions can remain later sensitivities without redefining the current estimand.

## Biological interpretation

The useful comparison is not `colour -> pollinator identity`. It is:

1. establish the raw colour response;
2. test whether it persists after `selfing_core` adjustment;
3. compare the *combination* of raw colour with independently measured floral form, tube depth and symmetry;
4. only then discuss concordance with predeclared bee-, butterfly- or bird-associated architectures.

Thus colour can strengthen the pollination interpretation without becoming a circular pollinator classifier.

## Claim ceiling

This audit can support statements such as:

- isolation is associated with a change in the representation of a specific reported flower colour;
- the five-colour composition vector changes with isolation;
- a colour shift persists after conditioning on measured reproductive assurance;
- the colour direction is concordant with a broader pollination-associated floral architecture when shape/tube/symmetry evidence agrees.

It cannot by itself support:

- colour intensity or visual conspicuousness to a particular animal visual system;
- reduced attraction investment as a measured physiological quantity;
- realized pollinator identity;
- historical pollinator loss;
- causal selection on flower colour;
- causal mediation from pollen limitation through colour evolution.
