# Chapter 1 v13 raw flower-colour audit

Status: reproduced secondary audit. This does **not** replace the canonical v13 H1–H4 submission surface unless explicitly promoted later.

## Why this audit exists

The final Chapter 1 database contains a large flower-colour axis, but the v13 primary plant response deliberately uses reproductive assurance and floral accessibility/generalization. Collapsing colour to `plain_colour`, or embedding it inside named pollination-syndrome scores, discards biologically useful information.

This audit therefore returns to the reported colour states themselves:

- `white`
- `red_pink`
- `yellow_orange`
- `blue_purple`
- `green_brown_inconspicuous`

No colour is treated as a pollinator identity. Raw colours are compared with raw floral form and tube-depth architecture only after their own isolation response is established.

## Estimand

For each island and evidence scope, the response for colour `c` is:

> share of focal-colour-resolved plant species explicitly reported to contain colour `c`.

Species with multiple reported colours remain multistate and contribute a success to every focal colour they explicitly contain. Species whose only colour state is non-focal do not enter the focal-colour denominator.

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

The primary floristic population for this v13 audit is `all_observed`, matching the current global-only paper surface.

## Reproduced execution

Raw five-colour run: `35210281545`  
Artifact: `10491539130`  
Digest: `sha256:d30ebbda7f1613373a39a9df9577e1b34cc32cddcd5c39b87719cd512972e263`

Full raw-colour + architecture + coupling run: `35216672430`  
Artifact: `10495522252`  
Digest: `sha256:8200c982d9a275c3d7db5fc25fa762f08428d51fa637a42c8312b067fdcbe377`

The dedicated tests, Ruff checks and repository-wide v13 audit passed before the frozen-data result was accepted.

## Raw-colour isolation response after conditioning on `selfing_core`

Flower-colour composition retains an isolation response after measured reproductive assurance is conditioned out.

### Northern mid-latitude

All-analysis:

- `red_pink`: `beta=-0.03209`, `p=0.00221`, `q=0.0111`;
- `green_brown_inconspicuous`: `beta=+0.02310`, `p=0.01385`, `q=0.0346`;
- five-colour joint Wald: `chi2=28.33`, `df=5`, `p=3.14e-5`.

Direct High/Medium:

- `red_pink`: `beta=-0.02526`, `p=0.00624`, `q=0.0312`;
- `yellow_orange`: `beta=-0.03598`, `p=0.0276`, `q=0.0690` after five-colour FDR;
- five-colour joint Wald: `chi2=38.24`, `df=5`, `p=3.38e-7`.

The clearest replicated component is therefore a decline in `red_pink`, not simply a generic increase in white.

### Tropical

All-analysis:

- `white`: `beta=+0.06704`, `p=0.000613`, `q=0.00307`;
- `green_brown_inconspicuous`: `beta=+0.04541`, `p=0.00277`, `q=0.00694`;
- `blue_purple`: `beta=-0.06016`, `p=0.0236`, `q=0.0393`;
- five-colour joint Wald: `chi2=33.73`, `df=5`, `p=2.70e-6`.

Direct High/Medium:

- `yellow_orange`: `beta=+0.04087`, `p=0.0150`, `q=0.0375`;
- `green_brown_inconspicuous`: `beta=+0.04643`, `p=0.000801`, `q=0.00401`;
- five-colour joint Wald: `chi2=22.93`, `df=5`, `p=3.49e-4`.

The tropical result is not a simple mirror of the northern result. Raw colour shifts occur, but the component pattern changes with evidence scope.

### Southern extratropical

All-analysis five-colour joint Wald: `chi2=20.08`, `df=5`, `p=0.00121`.

Direct High/Medium:

- `white`: `beta=+0.06435`, `p=3.23e-6`, `q=1.62e-5`;
- five-colour joint Wald: `chi2=50.26`, `df=5`, `p=1.23e-9`.

### Northern high latitude

The five-colour conditional vector was not supported in either evidence scope (`p=0.111` all-analysis; `p=0.289` direct-only).

## Raw colour × raw architecture

Two additional analyses keep both colour and architecture as raw states rather than converting them into a syndrome score.

### Joint prevalence

The first asks whether a raw colour × raw floral-form combination becomes more or less common among species with both trait families resolved. This layer detects combined phenotype shifts but can still reflect marginal changes in colour frequency.

### Colour-conditioned coupling

The stronger anti-confounding test asks:

> among species carrying a given raw colour, what fraction also carries the specified raw floral form or tube-depth state?

This removes the marginal frequency of the colour from the response.

The reproduced coupling results are documented in `docs/chapter1_v13_raw_colour_coupling_audit_20260917.md`.

Key outcomes are:

- northern high latitude: blue/purple coupling to butterfly-associated form and to intermediate/deep large-bee-associated tube architecture declines strongly with isolation in both evidence scopes;
- tropical direct evidence: yellow/orange coupling to deep butterfly-associated tube architecture increases with isolation;
- southern extratropical all-analysis: yellow/orange bird-associated form declines while intermediate/deep tube representation increases, giving a mixed rather than coherent named-syndrome response;
- northern mid-latitude: no colour-conditioned named architecture survives the FDR family, despite the robust `red_pink` raw-colour decline.

## Biological interpretation

The result supports keeping flower colour in Chapter 1 as a genuine secondary response family.

The evidence is now hierarchical:

1. raw colour composition changes with isolation;
2. the colour response persists after `selfing_core` adjustment;
3. in several contexts, raw architecture conditional on raw colour also changes with isolation.

That third result matters because it shows that the signal is not always just “there are fewer red flowers” or “there are more yellow flowers.” In some regions, the association between colour and specialized floral architecture itself is restructured with isolation.

However, the named pollination interpretation is context dependent. The northern-high-latitude result is compatible with erosion of specialized pollinator-facing architecture. The tropical direct result is compatible with retention or strengthening of a deep-tube component among yellow/orange flowers. Southern results are internally mixed. Northern-midlatitude colour shifts cannot currently be assigned to a particular guild after FDR correction.

Thus flower colour strengthens the paper most as a **raw visual/display phenotype family plus colour–architecture coupling analysis**, not as a one-colour-one-pollinator classifier.

## Relation to the attraction hypothesis

The data now support an isolation-associated change in pollinator-facing visual/display composition that is statistically separable from measured reproductive assurance. They also support context-specific changes in colour–architecture coupling.

They do **not** directly measure attraction intensity, pigment concentration, visual contrast in an animal visual system, UV signal, visitation rate or energetic investment. Therefore the strongest defensible wording is “pollinator-facing display composition and architecture” rather than “reduced attraction investment” as a directly measured quantity.

## Claim ceiling

Supported:

- isolation is associated with changes in specific reported flower colours;
- the five-colour composition vector changes with isolation in northern-midlatitude, tropical and southern-extratropical all-observed floras;
- the colour response persists after conditioning on measured reproductive assurance;
- northern-midlatitude `red_pink` decline is reproduced in both evidence scopes;
- some raw colour–architecture couplings change with isolation even after conditioning on the colour itself and on `selfing_core`;
- the raw colour results can be discussed for concordance with predeclared pollination-associated architecture when the structural evidence agrees.

Not supported:

- colour intensity or visual conspicuousness to a particular animal visual system;
- a raw colour uniquely identifies a realized pollinator;
- reduced attraction investment as a measured physiological quantity;
- historical pollinator loss;
- causal selection on flower colour;
- causal mediation from pollen limitation through colour evolution.
