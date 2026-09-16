# Chapter 1 H5 — GloPL reproductive-assurance moderation

## Question

The full-scale GloPL analysis established a positive association between distance from the six seeded major continental landmasses and experimental pollen limitation. This test asked a narrower functional question: **do species carrying predeclared reproductive-assurance states show a weaker distance-associated increase in pollen limitation?**

This is a moderation test, not a mediation or evolutionary-acquisition test.

## Frozen design

Design contract: `config/chapter1_h5_glopl_reproductive_assurance_moderation_v1.yml`  
Design freeze commit: `baf5b7994171e31c01d91e9a4c936408bc4db6cf`

The three traits were frozen separately before trait-overlap and GloPL effect inspection:

- self compatibility: `SC` versus `SI`;
- selfing mating system: `predominantly_selfing` or `obligate_selfing` versus `predominantly_outcrossing`;
- autonomous selfing: `autonomous` or `delayed` versus `absent`.

Ambiguous/mixed states were excluded. Matching was exact after formatting-only normalization; no fuzzy, synonym or genus fallback was allowed.

The predeclared buffering prediction was a **negative distance × reproductive-assurance interaction**. Each publication retained total analysis weight 1 and uncertainty was publication-cluster robust. The full-scale GloPL distance standardization and measurement controls were inherited from the frozen parent analysis.

## Outcome-blind support audit

Canonical workflow: **35093274622**  
Job: **104784426051**  
Artifact: **10444749163**  
Artifact digest: `sha256:317221b0315137d2cfc0591d9f6f146c6956106d50d8a5b1dc782096bab55da9`

Before any GloPL effect columns were read, the exact species join yielded **642 GloPL species** with at least one frozen reproductive-assurance assignment.

### Self compatibility

Global support passed:

- 499 matched species;
- SI: 220 species, 267 sites, 199 publications, 75 offshore sites;
- SC: 279 species, 348 sites, 250 publications, 106 offshore sites.

The North–Tropical diagnostic gate also passed.

### Selfing mating system

Global support failed:

- 106 matched species total;
- predominantly outcrossing: 95 species, 172 sites, 125 publications, 40 offshore sites;
- predominantly/obligately selfing: **11 species, 14 sites, 11 publications, only 1 offshore site**.

This contrast is therefore **not evaluable**, not negative biological evidence.

### Autonomous selfing

Global support passed:

- 558 matched species;
- absent: 418 species, 504 sites, 379 publications, 130 offshore sites;
- autonomous/delayed: 140 species, 176 sites, 120 publications, 77 offshore sites.

The North–Tropical diagnostic gate also passed.

## Primary results

### Self compatibility

On 769 cells, 429 publications and 591 sites:

- SI distance slope: **+0.0362**;
- SC distance slope: **+0.0712**;
- distance × SC interaction: **+0.0351**;
- SE = 0.0754;
- two-sided p = **0.6420**;
- predeclared one-sided buffering p = **0.6790**.

Both frozen sensitivities also had **positive**, rather than buffering, interactions:

- supplemental-only: +0.1050;
- no-zero-constant: +0.0389.

Thus self compatibility does **not** support the functional-buffering prediction.

### Autonomous selfing

On 865 cells, 469 publications and 640 sites:

- no-autonomous-selfing distance slope: **+0.1175**;
- autonomous/delayed distance slope: **+0.0439**;
- distance × autonomous-selfing interaction: **−0.0736**;
- SE = 0.0723;
- two-sided p = **0.3086**;
- predeclared one-sided buffering p = **0.1543**.

Both frozen sensitivities retain the expected negative direction:

- supplemental-only: −0.0268;
- no-zero-constant: −0.0681.

This is directionally compatible with buffering but does **not** pass the frozen statistical gate.

## Secondary North–Tropical diagnostics

These were predeclared as descriptive diagnostics that cannot rescue the global route.

Self compatibility:

- northern distance × trait = +0.0342;
- tropical = −0.0250;
- tropical-minus-northern = −0.0592, p = 0.7693.

Autonomous selfing:

- northern distance × trait = −0.1013;
- tropical = −0.0557;
- tropical-minus-northern = +0.0456, p = 0.8034.

There is no supported North–Tropical difference in buffering for either evaluable trait.

## Decision

Frozen family rule: at least **2 of 3** reproductive-assurance traits must satisfy the global buffering rule.

Observed:

- self compatibility: **not supported**;
- selfing mating system: **not evaluable**;
- autonomous selfing: **directionally compatible but not supported**.

Therefore:

`reproductive_assurance_buffering_not_supported`

The supported broad global isolation–pollen-limitation gradient does **not** acquire a demonstrated functional bridge to the Chapter 1 reproductive-assurance response through these matched traits.

## Claim ceiling

This result does not imply that reproductive assurance is irrelevant to island isolation. In particular, it does not test whether island colonization or chronic pollen limitation caused the evolutionary acquisition of self compatibility, selfing or autonomous selfing. It also does not show that the under-supported selfing-mating-system route is absent.

What it does show is narrower: **within the current exact species-matched GloPL data, species carrying the predeclared reproductive-assurance states do not provide frozen-rule evidence that the isolation-associated increase in experimental pollen limitation is functionally buffered.**
