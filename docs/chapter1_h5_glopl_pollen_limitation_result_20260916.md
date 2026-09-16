# Chapter 1 H5 — exact-island GloPL pollen-limitation result

## Decision

The outcome-blind geographic preflight passed, so the separately frozen GloPL effect-size model was evaluated. The proposed isolation-associated pollen-limitation gradient is **not supported** under that frozen model.

This is more direct mechanistic evidence than the exact-island pollinator-channel occurrence proxies because GloPL uses pollen-supplementation experiments and measures reproductive pollen limitation. It is still cross-sectional across islands and does not measure pollinator abundance, visitation, visitor identity, historical loss, or a selection gradient on floral traits.

## Provenance

- GloPL repository: `idiv-biodiversity/pollen-limitation-data-descriptor`
- pinned source commit: `abbe193770981b2c61713e75a0f01e004bdded8e`
- source file: `Output/GloPL.csv`
- source SHA-256: `264ca8c2237f126b3c209fcba039e19b291324ddf6e7a96992dc8235a55984de`
- outcome-blind preflight run: **35085378171**
- preflight artifact: **10442126584**
- model freeze commit: `0f7479c3607f86333c8d904a4786c0b938912b0b`
- final analysis run: **35086524129**
- analysis artifact: **10442910503**
- artifact digest: `sha256:5d39d3658100edabcbb9a7f6b5e6d7130ffc7d32e643e7900fcebfbd137efaa0`

## Outcome-blind admission gate

The preflight used only latitude, longitude, DOI, author, year and accepted species name. Pollen-limitation effect columns were forbidden at that stage.

Exact `within` matching to the frozen 8,265-island GSHHG universe found:

- 485 GloPL rows on 37 frozen islands;
- 197 unique coordinate sites;
- 103 study/publication keys;
- northern mid-latitude: **15 islands / 55 studies**;
- tropical: **14 islands / 23 studies**.

The frozen minimum was 10 islands and 5 studies in each primary context, so the effect-size stage was admitted without changing the threshold.

## Frozen effect-size estimand

`PL_Effect_Size` is the published GloPL master log response ratio comparing hand-pollen treatment with natural reproductive output. Positive values mean a supplementation advantage and therefore stronger experimental pollen limitation at the measured reproductive-output endpoint.

The primary unit was publication × island. Multiple GloPL rows from a publication on one island were averaged first, then each publication received total analysis weight one. The model used the North/Tropical context, standardized mainland distance, distance × tropical, island area, climate PC1–4 and spatial-block robust covariance.

The frozen directional predictions were:

1. pollen limitation increases with mainland isolation in northern mid-latitude islands;
2. the isolation slope is weaker in tropical islands.

## Primary result

Support after finite effects and frozen exact-island matching:

- **262** effect rows;
- **80** publication × island cells;
- 56 northern mid-latitude cells;
- 24 tropical cells;
- **28 unique islands**;
- **77 publications**;
- **25 spatial blocks**.

Northern distance slope:

- estimate = **+0.3762**;
- SE = 0.3215;
- approximate 95% normal CI = **[-0.254, 1.006]**;
- frozen one-sided positive p = **0.1209**.

The estimate is in the predicted direction but does not pass the frozen directional gate.

Tropical-minus-North interaction:

- estimate = **-0.0311**;
- SE = 0.2933;
- approximate 95% normal CI = **[-0.606, 0.544]**;
- one-sided negative p = **0.4577**;
- two-sided p = 0.9155.

The primary interaction is effectively near zero relative to its uncertainty. The tropical distance slope is +0.3451 (two-sided p = 0.4774).

## Frozen sensitivities

### Supplemental-pollen treatment only

- 50 publication × island cells;
- North slope = **+0.4318**, one-sided p = 0.1498;
- Tropical-minus-North interaction = **+0.3006**.

The North slope remains positive but unsupported, while the interaction changes to the opposite of the frozen prediction.

### Excluding GloPL zero-constant cases

- 77 cells;
- North slope = **+0.3231**, one-sided p = 0.1653;
- interaction = **-0.0893**, one-sided p = 0.3868.

### Equal-island guardrail

- 28 islands, exactly 14 North and 14 tropical;
- North slope = **+0.7473**, one-sided p = 0.1216;
- interaction = **-0.2284**, one-sided p = 0.3440.

Thus all frozen versions retain a positive northern point estimate, but none supports the northern directional estimand. The predicted North–Tropical divergence is also not robust: it is unsupported in the primary model and reverses sign in the supplemental-only sensitivity.

## H5 interpretation

The result is classified as:

> **`pollen_limitation_gradient_not_supported`**

This adds an important independent layer to H5. The earlier pollinator-channel analyses could fail because occurrence is a coarse proxy for service. GloPL bypasses that problem by directly manipulating pollen supply, yet the available exact-island experiments still do not recover the proposed northern isolation gradient or a North–Tropical difference under the frozen model.

The appropriate synthesis is therefore not “pollinators are irrelevant.” It is:

> **The current global data do not identify increasing pollination-service limitation with island isolation as the common upstream generator of the Chapter 1 plant-trait response.**

The positive northern point estimates are compatible with such a gradient but remain too imprecise to support it, and the context-specific prediction is not reproduced.

## Claim ceiling

Do not infer from this result that:

- pollen limitation is absent on islands;
- pollinators are unimportant to island floras;
- pollinator abundance or visitation does not decline with isolation;
- there has been no historical pollinator loss;
- pollination cannot affect floral evolution locally;
- the plant-side H2 response is unrelated to pollination in every system.

The result addresses only the frozen cross-sectional association between experimental pollen limitation and island isolation in the available exact-island GloPL subset.
