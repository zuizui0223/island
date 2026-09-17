# Chapter 1 v13 raw flower-colour coupling audit

Status: reproduced secondary audit attached to PR #233.

This audit asks whether raw flower-colour states become more or less tightly coupled to raw floral architecture with isolation, without converting either trait family into a weighted pollination-syndrome score.

The key estimand is conditional architecture frequency among species carrying a focal colour. For example, `red_pink__butterfly_form_given_colour` is the share of red/pink species on an island whose raw floral form is salverform, tubular, spurred or funnel/trumpet. This removes changes in the marginal abundance of red/pink itself from the coupling response.

The same logic is applied to raw tube-depth combinations. These are concordance tests only; they do not identify realized pollinator identity or causal selection.

## Frozen execution

Canonical coupling run: `35216672430`  
Artifact: `10495522252`  
Digest: `sha256:8200c982d9a275c3d7db5fc25fa762f08428d51fa637a42c8312b067fdcbe377`

The dedicated raw-colour tests, colour-by-architecture tests, coupling tests, Ruff checks and the repository-wide v13 audit all passed. Results below are from the conditional-on-`selfing_core`, confirmatory-support models.

## Reproduced coupling results

### Northern mid-latitude

No raw colour-conditioned architecture combination survives the within-context FDR family at `q < 0.05` in either evidence scope.

The clearest near-signal is `yellow_orange__large_bee_form_given_colour`:

- all-analysis: `beta=-0.05962`, `p=0.00449`, `q=0.0629`;
- direct High/Medium: `beta=-0.04257`, `p=0.0357`, `q=0.2498`.

`red_pink__bird_deep_tube_given_colour` is also negative nominally in both scopes but does not survive FDR.

Therefore the robust northern-midlatitude result remains the raw colour shift itself (`red_pink` declines with isolation in both scopes), not a demonstrated shift toward or away from one named pollination-associated architecture.

### Northern high latitude

This context shows the clearest replicated breakdown of specialized colour–architecture coupling with isolation.

Replicated in both evidence scopes:

- `blue_purple__butterfly_form_given_colour`
  - all-analysis: `beta=-0.20186`, `p=0.000205`, `q=0.00144`;
  - direct: `beta=-0.19967`, `p=9.7e-5`, `q=0.000681`.
- `blue_purple__large_bee_deep_tube_given_colour`
  - all-analysis: `beta=-0.25505`, `p=0.000193`, `q=0.00144`;
  - direct: `beta=-0.24724`, `p=1.9e-5`, `q=0.000266`.

Direct-only additionally supports:

- `yellow_orange__butterfly_form_given_colour`: `beta=-0.18773`, `p=0.00587`, `q=0.0274`.

The important interpretation is not that blue/purple identifies butterflies or large bees. Rather, among species retaining blue/purple colour, the representation of long-tongued/specialized floral architectures declines strongly with isolation. This is compatible with a weakening of specialized pollinator-facing architecture in remote northern-high-latitude floras.

### Tropical

Direct High/Medium evidence supports one positive coupling:

- `yellow_orange__butterfly_deep_tube_given_colour`: `beta=+0.11703`, `p=0.00253`, `q=0.0355`.

The all-analysis estimate for the broader yellow/orange butterfly-form coupling is negative only nominally and does not survive FDR.

Thus the strongest tropical raw coupling signal is not generalized floral simplification. Among yellow/orange species, deep-tube architecture becomes more represented with isolation in the direct evidence layer. This is concordant with retention or strengthening of a long-tongued-insect-associated structural component, but does not identify butterflies as the realized visitors.

### Southern extratropical

All-analysis evidence supports a mixed yellow/orange pattern:

- `yellow_orange__bird_form_given_colour`: `beta=-0.12495`, `p=0.00102`, `q=0.0129`;
- `yellow_orange__bird_deep_tube_given_colour`: `beta=+0.19238`, `p=0.00276`, `q=0.0129`;
- the same positive tube-depth estimate is shared by `yellow_orange__large_bee_deep_tube_given_colour` because the frozen bird and large-bee tube templates both classify intermediate/deep states as concordant.

These results do not form one coherent named bird syndrome: bird-associated form declines while intermediate/deep tube representation increases, and the effects do not survive FDR in the direct-only layer. The appropriate conclusion is a context-specific restructuring of yellow/orange floral architecture rather than evidence for a single pollinator guild.

## What the colour data now add

The raw-colour analyses establish three distinct levels of evidence:

1. **colour composition changes with isolation**: strongly supported in northern-midlatitude, tropical and southern-extratropical contexts after conditioning on `selfing_core`;
2. **raw colour × raw architecture joint prevalence changes**: supported for selected combinations, but this estimand can still be driven by marginal colour frequency;
3. **architecture conditional on raw colour changes**: supported in northern-high-latitude, tropical and southern-extratropical contexts, demonstrating that some colour–structure associations themselves change with isolation.

This justifies retaining flower colour as a real secondary phenotype family in Chapter 1 rather than reducing it to `plain_colour` or hiding it inside a weighted syndrome score.

## Relation to the attraction / display hypothesis

The present result supports an isolation-associated change in **pollinator-facing floral display composition and colour–architecture coupling** that is not reducible to measured reproductive assurance.

It does not directly measure attraction intensity, visual contrast, pigment concentration, UV signal, visitation rate or energetic investment. Therefore `reduced attraction investment` remains an interpretation rather than a measured estimand. The strongest defensible wording is that isolation is associated with shifts in visual/display states and, in some contexts, with changes in how those states are coupled to specialized floral architecture.

## Claim ceiling

Supported:

- raw flower-colour composition changes with isolation independently of measured `selfing_core`;
- some raw colour–architecture couplings change with isolation after conditioning on colour frequency and `selfing_core`;
- northern-high-latitude isolation is associated with reduced blue/purple coupling to specialized/deep floral architecture;
- tropical direct evidence supports increased yellow/orange coupling to deep-tube architecture;
- these results can be compared with predeclared pollination-associated architectures.

Not supported:

- a raw colour uniquely identifies a realized pollinator;
- northern-midlatitude `red_pink` decline is specifically caused by loss of butterflies, birds or large bees;
- colour proves reduced attraction investment in a physiological or behavioral sense;
- historical pollinator loss caused the observed colour shift;
- the coupling tests are causal mediation analyses.
