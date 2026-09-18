# Chapter 1 v14 preflight results — 2026-09-18

Status: **local reproduction from the exact frozen input artifacts used by the v14
workflow; CI artifact still pending**.

This note is not a result lock. It records the numerical preflight used to check
whether the H1/H2 reorganization changes the scientific conclusion before the PR
workflow is promoted.

## Frozen inputs

- progressive-analysis run: 34232450884
- all-data primary run: 35100991898
- H1 model: beta-binomial logit, spatial-block cluster-robust sandwich covariance
- H1 covariates: log island area + climate PC1–PC4 + standardized log distance
- H2 continuous models: equal-island OLS with spatial-block cluster-robust covariance
- H2 plain-colour model: beta-binomial logit with selfing_core + H1 covariates
- evidence scopes: all-analysis-eligible primary and Direct-only sensitivity
- flora scope reported here: all_observed

## H1 — seven-response island syndrome

Adding plain colour to the six v13 atomic responses does **not** remove the
multivariate isolation response. The seven-response joint Wald test remains supported
in all four predeclared geographic replication strata in both evidence scopes.

### All-analysis-eligible

| context | joint p | BH q |
|---|---:|---:|
| northern mid-latitude | 1.77e-6 | 2.36e-6 |
| northern high-latitude | 2.61e-6 | 2.61e-6 |
| tropical | 1.89e-7 | 3.78e-7 |
| southern extratropical | 4.01e-13 | 1.60e-12 |

### Direct-only

| context | joint p | BH q |
|---|---:|---:|
| northern mid-latitude | 0.00461 | 0.00461 |
| northern high-latitude | 4.99e-9 | 9.98e-9 |
| tropical | 4.42e-6 | 5.90e-6 |
| southern extratropical | 2.73e-24 | 1.09e-23 |

### Three biological domains

All-analysis-eligible standardized isolation-coefficient means:

| context | reproductive assurance | colour dulling | accessibility/generalization | equal-domain orientation |
|---|---:|---:|---:|---:|
| northern mid-latitude | +0.0190 | -0.0001 | +0.0200 | +0.0130 |
| northern high-latitude | +0.0730 | +0.0188 | +0.1783 | +0.0900 |
| tropical | +0.1838 | +0.0367 | +0.0810 | +0.1005 |
| southern extratropical | +0.1490 | +0.0808 | +0.0175 | +0.0824 |

Direct-only equal-domain orientation is also positive in all four contexts:
+0.0159, +0.0925, +0.0854, +0.0855.

The important qualification is colour. The beta-binomial plain-colour coefficient is
approximately zero in northern mid-latitudes and positive in the other three contexts.
Therefore v14 supports a recurrent **three-domain syndrome direction overall**, but it
does not claim a statistically uniform colour-dulling coefficient in every region.

## H2 — selfing versus pollinator-facing floral decomposition

### H2a reproductive-assurance route

The selfing_core isolation coefficient is positive in all four contexts under both
evidence scopes. In the all-analysis scope it is strongest in tropical and southern
extratropical islands; northern mid-latitude is positive but imprecise.

### H2b accessibility after conditioning on selfing_core

The generalized_accessible isolation coefficient remains positive in all four contexts
after conditioning on selfing_core.

All-analysis-eligible:

| context | conditional distance beta | p | primary-H2b q |
|---|---:|---:|---:|
| northern mid-latitude | +0.0214 | 0.0833 | 0.133 |
| northern high-latitude | +0.1248 | 0.000360 | 0.00144 |
| tropical | +0.0733 | 0.00185 | 0.00494 |
| southern extratropical | +0.0536 | 0.1085 | 0.145 |

Direct-only has the same positive sign in all four contexts; northern high-latitude is
FDR-supported, tropical is borderline after the eight-test H2b correction.

The secondary attraction_shift coefficient is also positive in all four contexts in
both evidence scopes, with the clearest support in northern high-latitude and tropical
islands.

### Raw colour × raw architecture concordance after conditioning on selfing_core

The v14 pollination-syndrome concordance layer does **not** use the weighted
`large_bee_like / butterfly_like / bird_like` scores as its primary evidence. It
reuses the reproduced v13 raw-state analyses.

Raw five-colour composition remains associated with isolation after conditioning on
`selfing_core` in three of four contexts:

| context | all-analysis joint p | Direct joint p |
|---|---:|---:|
| northern mid-latitude | 3.14e-5 | 3.38e-7 |
| northern high latitude | 0.111 | 0.289 |
| tropical | 2.70e-6 | 3.49e-4 |
| southern extratropical | 0.00121 | 1.23e-9 |

The strongest replicated raw-colour component in northern mid-latitudes is a decline
in `red_pink`, not simply an increase in white. Tropical and southern colour
composition also changes, but the component pattern differs among contexts and evidence
scopes.

The stronger syndrome-concordance test conditions architecture on the raw colour itself:
`P(raw architecture | raw colour, architecture resolved)`. This removes changes in
the marginal frequency of the focal colour.

Replicated/supporting results include:

- northern high latitude: among `blue_purple` species,
  `butterfly_form_given_colour` declines in both scopes
  (all `beta=-0.20186, q=0.00144`; Direct `beta=-0.19967, q=0.000681`);
- northern high latitude: `blue_purple × intermediate/deep large-bee-associated tube`
  declines in both scopes (all `beta=-0.25505, q=0.00144`; Direct
  `beta=-0.24724, q=0.000266`);
- tropical Direct evidence: `yellow_orange × deep-tube butterfly-associated
  architecture` increases (`beta=+0.11703, q=0.0355`);
- southern extratropical all-analysis: yellow/orange bird-associated form declines
  while intermediate/deep tube representation increases (both `q=0.0129`), a mixed
  restructuring rather than one coherent named syndrome;
- northern mid-latitude: no colour-conditioned named architecture survives the FDR
  family, despite the replicated raw `red_pink` decline.

The correct H2 interpretation is therefore **context-specific raw display–architecture
reorganization**, not a weighted pollinator-score effect and not realized pollinator
identity.

## H3/H4 status under v14

### H3 — unchanged analysis, renumbered role

H3 is the frozen full-global GloPL distance analysis that was v13 H2. The statistical
model and result do not change under v14: standardized isolation coefficient
`beta=+0.07937`, SE `0.03773`, two-sided `p=0.03543`, one-sided positive
`p=0.01772` across 2,969 effects, 1,248 sites and 919 publications.

The change is conceptual placement: H3 now follows the plant-side H2 decomposition and
serves as independent evidence that the proposed ecological pressure itself increases
with isolation.

### H4 — evidence hierarchy updated, not simply unchanged

The v13 exact-species bridge remains the discovery layer and remains explicitly
post-hoc. Its strongest association is autonomous selfing
(`beta=-0.446724, p=2.60e-8`), with weaker global architecture concordance for
actinomorphy and generalized form.

PR #236 adds two validation layers without relabelling that discovery:

1. a prospectively frozen post-2015 wild-plant temporal replication; the outcome-blind
   support gate reached only 9 matched species for each co-primary family, below the
   frozen minimum of 30 species and 10 publications, so the test stopped before outcome
   unblinding and is **not evaluable**, not a biological null;
2. an independent PolLimCrop agricultural-domain transportability test with frozen
   response mapping and support gates; this is still separate/pending and cannot rescue
   the wild temporal support failure or make v13 H4 confirmatory for wild floras.

The v14 discovery layer now also aligns directly to the two H2 family scores.

- **H2-equivalent reproductive assurance:** at least two of self-compatibility,
  selfing mating system and autonomous selfing; 470 species, 417 publications,
  739 cells. Family-score coefficient `beta=-0.31255`, two-sided `p=0.00222`.
  The no-zero sensitivity remains supported (`beta=-0.31726, p=0.00229`), while
  supplemental-only remains negative but not supported (`beta=-0.13138, p=0.148`).
- **H2-exact accessibility/generalization:** generalized form weight 1.0,
  actinomorphy 0.75 and shallow/open tube 1.0, requiring at least two components;
  228 species, 232 publications, 403 cells. `beta=-0.32932, p=0.00515`;
  supplemental-only `beta=-0.24989, p=0.0367`; no-zero
  `beta=-0.30663, p=0.0134`.
- The equal-weight accessibility definition subsequently frozen for prospective H4b
  gives almost the same post-hoc discovery result (`beta=-0.32509, p=0.00636`).

Thus the **same two broad H2 response families that increase with isolation are
associated with lower current pollen limitation in the post-hoc GloPL discovery layer**.
Accessibility/generalization is the more measurement-robust family bridge.

The prospective validation predictions remain separately frozen:

- H4a: autonomous selfing -> lower current pollen limitation;
- H4b: higher equal-weight accessibility/generalization
  (generalized form + actinomorphy + shallow/open tube) -> lower current pollen limitation.

A global primary colour bridge is **not** added because H1/H2 colour responses are
context dependent and do not define one globally predicted functional direction.

## Interpretation

The reordering is empirically coherent:

1. **H1:** a recurrent island-syndrome response remains after colour is included.
2. **H2a:** reproductive assurance increases with isolation.
3. **H2b:** the structural accessibility/generalization response is not reducible to
   selfing_core; this is the strongest plant-side evidence for a parallel
   pollinator-facing route.
4. **Colour is part of H1 but not the strongest H2 mechanism marker.** It shows a
   broad positive tendency in three contexts, with a particularly strong
   selfing-independent signal in the southern extratropical stratum.
5. **H3 and H4 remain logically downstream:** H3 independently tests whether pollen
   limitation increases with isolation; H4 tests trait–pollen-limitation functional
   compatibility.

The conditional H2 models do not prove causal mediation or direct pollinator
selection. They reject the narrower explanation that the measured selfing core alone
accounts for the full floral-architecture response.
