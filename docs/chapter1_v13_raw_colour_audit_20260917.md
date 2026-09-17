# Chapter 1 v13 raw flower-colour audit

Status: exploratory secondary audit. This does **not** change the canonical v13 H1–H4 submission surface unless explicitly promoted later.

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

Within each geographic replication stratum, all five raw-colour isolation coefficients are also tested jointly with a spatial-block cluster-robust Wald test.

This asks whether the raw flower-colour vector changes with isolation even when no single colour alone is treated as the syndrome.

## Evidence scopes

- `all`: all analysis-eligible colour evidence in the frozen final species-axis snapshot;
- `direct`: High/Medium colour evidence only.

The primary floristic population for this v13 audit is `all_observed`, matching the current global-only paper surface. Native/status-resolved versions can remain later sensitivities without redefining the current estimand.

## Reproduced execution

Canonical PR run: `35210281545`  
Artifact: `10491539130`  
Digest: `sha256:d30ebbda7f1613373a39a9df9577e1b34cc32cddcd5c39b87719cd512972e263`

The dedicated tests and Ruff checks passed before the frozen-data analysis. Both evidence scopes completed and uploaded row-level counts, per-colour models and five-colour joint Wald tests.

### Raw-colour isolation response after conditioning on `selfing_core`

The strongest current result is that **flower-colour composition retains an isolation response after measured reproductive assurance is conditioned out**.

#### Northern mid-latitude

All-analysis:

- `red_pink`: `beta=-0.03209`, `p=0.00221`, `q=0.0111`;
- `green_brown_inconspicuous`: `beta=+0.02310`, `p=0.01385`, `q=0.0346`;
- five-colour joint Wald: `chi2=28.33`, `df=5`, `p=3.14e-5`.

Direct High/Medium:

- `red_pink`: `beta=-0.02526`, `p=0.00624`, `q=0.0312`;
- `yellow_orange`: `beta=-0.03598`, `p=0.0276`, `q=0.0690` after five-colour FDR;
- five-colour joint Wald: `chi2=38.24`, `df=5`, `p=3.38e-7`.

Thus the northern colour signal is not simply a generic increase in white. The clearest replicated component is a decline in `red_pink`; the broader all-analysis layer also shows an increase in green/brown inconspicuous states.

#### Tropical

All-analysis:

- `white`: `beta=+0.06704`, `p=0.000613`, `q=0.00307`;
- `green_brown_inconspicuous`: `beta=+0.04541`, `p=0.00277`, `q=0.00694`;
- `blue_purple`: `beta=-0.06016`, `p=0.0236`, `q=0.0393`;
- five-colour joint Wald: `chi2=33.73`, `df=5`, `p=2.70e-6`.

Direct High/Medium:

- `yellow_orange`: `beta=+0.04087`, `p=0.0150`, `q=0.0375`;
- `green_brown_inconspicuous`: `beta=+0.04643`, `p=0.000801`, `q=0.00401`;
- five-colour joint Wald: `chi2=22.93`, `df=5`, `p=3.49e-4`.

The tropical result is therefore **not** a simple mirror of the northern result. Raw colour shifts occur, but the component pattern changes with evidence scope. This is exactly why colour should be retained as a multicomponent response rather than compressed to one `plain_colour` contrast.

#### Southern extratropical

All-analysis five-colour joint Wald: `chi2=20.08`, `df=5`, `p=0.00121`.

Direct High/Medium:

- `white`: `beta=+0.06435`, `p=3.23e-6`, `q=1.62e-5`;
- five-colour joint Wald: `chi2=50.26`, `df=5`, `p=1.23e-9`.

#### Northern high-latitude

The five-colour conditional vector was not supported in either evidence scope (`p=0.111` all-analysis; `p=0.289` direct-only).

## Biological interpretation

The result supports keeping flower colour in Chapter 1 as a genuine secondary response family.

The useful comparison is not `colour -> pollinator identity`. It is:

1. establish the raw colour response;
2. test whether it persists after `selfing_core` adjustment;
3. compare the *combination* of raw colour with independently measured floral form, tube depth and symmetry;
4. only then discuss concordance with predeclared bee-, butterfly- or bird-associated architectures.

The current raw-colour result already clears steps 1 and 2. The next useful analysis is therefore step 3: test prespecified **raw colour × raw architecture combinations** without converting them into a weighted pollination-syndrome score.

Examples include:

- `red_pink` with tubular / salver / spurred / funnel forms;
- `yellow_orange` with tubular / funnel / bell forms;
- `blue_purple` or `yellow_orange` with bilabiate / papilionaceous / bell / funnel forms;
- the same colours crossed with deep/intermediate versus open/shallow tube states.

These combinations can be described as concordant with butterfly-, bird- or large-bee-associated floral architectures only when the structural component agrees. The analysis must not infer the realized visitor from colour alone.

## Claim ceiling

This audit supports statements such as:

- isolation is associated with a change in the representation of specific reported flower colours;
- the five-colour composition vector changes with isolation in northern-midlatitude, tropical and southern-extratropical all-observed floras;
- the colour-vector response persists after conditioning on measured reproductive assurance;
- northern-midlatitude `red_pink` decline is reproduced in both evidence scopes;
- the colour direction can be compared with broader pollination-associated floral architecture when shape/tube/symmetry evidence agrees.

It does not by itself support:

- colour intensity or visual conspicuousness to a particular animal visual system;
- reduced attraction investment as a measured physiological quantity;
- realized pollinator identity;
- historical pollinator loss;
- causal selection on flower colour;
- causal mediation from pollen limitation through colour evolution.
