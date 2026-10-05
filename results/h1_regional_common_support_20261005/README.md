# H1 regional common-support diagnostic — 2026-10-05

Status: **post-hoc diagnostic; not part of the frozen Chapter 1 submission inference.**

This diagnostic asks whether apparent regional differences in Chapter 1 depend on the four regions being observed over different mainland-distance distributions.

## Why this was needed

The final H1 regions occupy very different isolation ranges. In the primary All H1 union:

| Region | Islands | median distance (km) | IQR (km) | 5–95% (km) |
|---|---:|---:|---:|---:|
| Northern mid-latitude | 2,173 | 12.21 | 1.53–127.56 | 0.092–783.02 |
| Northern high latitude | 411 | 79.85 | 2.60–887.58 | 0.251–1,413.02 |
| Tropical | 1,493 | 449.79 | 10.84–1,257.94 | 0.227–3,388.14 |
| Southern extratropical | 302 | 62.98 | 3.42–618.23 | 0.107–2,201.71 |

Thus northern mid-latitude islands are strongly concentrated near continental land, whereas the tropical set contains a much longer remote-island tail.

Two non-tuned common-support windows were defined from the primary All H1 union:

- **common 5–95% support:** 0.2515–783.0171 km, the intersection of all four regional 5th–95th percentile ranges;
- **common IQR support:** 10.8445–127.5559 km, the intersection of all four regional interquartile ranges.

The common-IQR window is too sparse for strong inference, especially in the southern region (68 union islands / 8 spatial blocks; several traits fall below the frozen 50-island threshold). It is therefore retained only as a stress test. The 5–95% intersection is the useful common-support diagnostic.

## Seven-trait H1 under common 5–95% support

The frozen beta-binomial H1 model was re-fit without changing covariates or trait definitions.

Primary All evidence retained:

- northern mid-latitude: 1,812 islands / 72 blocks (83.4% of full union);
- northern high latitude: 262 / 47 (63.7%);
- tropical: 817 / 68 (54.7%);
- southern extratropical: 203 / 21 (67.2%).

For northern mid-latitudes, the three full-support nominal H1 associations remain nominally supported under common support:

- self-compatibility;
- selfing mating system;
- generalized form.

The two sign changes among the seven All traits are plain colour and actinomorphy, both of which are unsupported in the full model. Thus the principal northern-midlatitude H1 signal is not an artefact of having many near-mainland islands.

The tropical All model falls from five nominally supported traits to two under common support, although all seven trait slopes retain their full-support signs. This indicates that part of the tropical precision/strength is supplied by the far-isolation tail.

The direct-evidence sensitivity shows the same broad pattern: northern-midlatitude selfing mating system remains supported, whereas tropical supported traits fall from five to one.

## Pollination-syndrome concordance under common support

A second diagnostic re-fit the frozen raw **architecture conditional on colour** models using the corrected geography, selfing-core adjustment, and the same 0.2515–783.0171 km common-support window. The full corrected results were first replayed and matched the committed corrected result tables.

### Northern mid-latitude

The corrected full model contains a strong negative
`yellow_orange__large_bee_form_given_colour` association:

- All full: beta = -0.08038, q = 1.08e-5;
- All common support: beta = -0.06904, q = 0.00663;
- Direct full: beta = -0.06324, q = 0.000903;
- Direct common support: beta = -0.05131, nominal p = 0.00620, q = 0.0869.

The direction and magnitude therefore persist under like-for-like distance support. Loss of Direct FDR support is a precision/multiplicity issue after restriction, not a reversal. The earlier statement that northern mid-latitudes have no FDR-supported conditional architecture shift came from the pre-geography-correction surface and is superseded.

### Northern high latitude

The strongest loss of specialized colour–architecture coupling is robust to common support. For example:

- blue/purple × butterfly-associated form remains negative and FDR-supported in All and Direct;
- blue/purple × large-bee-associated intermediate/deep tube remains negative and FDR-supported in All and Direct.

This is the clearest region in which isolation is associated with weakening of specialized/long-tongued-insect-associated floral architecture.

### Tropics

The previously highlighted Direct
`yellow_orange__butterfly_deep_tube_given_colour` association is **not common-support robust**:

- Direct full: beta = +0.11584, q = 0.03949;
- Direct common support: beta = +0.00844, q = 0.92348.

The All layer is null in both versions. The tropical positive butterfly-like/deep-tube signal therefore comes from information in the broader tropical isolation range, especially the remote-island tail, rather than a like-for-like response over the shared regional distance range.

This does not make the full tropical result invalid: remote tropical islands are real parts of the tropical system. It does mean that the result cannot be cleanly interpreted as a biome/pollinator difference from northern regions without also acknowledging the different isolation regime and possible nonlinearity.

### Southern extratropical

The full All result is internally mixed: yellow/orange bird-associated form decreases while bird-associated intermediate/deep tubes increase. Under common support both effects retain their signs and substantial magnitudes but lose FDR support with only ~19 spatial blocks in the focal models. This remains evidence of structural reorganization, not a coherent bird syndrome.

## Main interpretation

The initial concern was that northern-midlatitude biology might be hard to distinguish because this region contains many islands close to continents. The common-support diagnostics do **not** support that as the main explanation.

Instead:

1. the northern-midlatitude reproductive/generalized H1 signal persists on common distance support;
2. the corrected northern-midlatitude large-bee-associated architecture decline also persists;
3. northern-high-latitude specialized architecture decline is especially robust;
4. the positive tropical butterfly-like deep-tube signal is the result most sensitive to unequal regional isolation support.

The safe conceptual conclusion is therefore not “northern mid-latitude is too near the mainland to detect a pollinator effect.” A better statement is:

> **regional syndrome concordance is partly a function of which isolation regime is represented. Northern specialization loss is visible even under common support, whereas the tropical positive long-tongued-insect signal is specific to the broader tropical isolation range.**

## Submission boundary

Do not promote these 2026-10-05 common-support numerical results into the current Chapter 1 submission. They are post-hoc diagnostics.

They should be used to:

- prevent over-interpretation of the tropical full-range syndrome signal as a like-for-like regional contrast;
- correct discussion text to the already-existing 2026-09-24 corrected raw-coupling results;
- motivate a future explicit nonlinear/common-support analysis of regional island syndromes.

## Reproducibility

- H1 common-support workflow run: **37281904325**
- H1 artifact: **11332262722**
- Pollination-syndrome common-support workflow run: **37282896938**
- Pollination-syndrome artifact: **11332724753**
- H1 diagnostic: `scripts/audit_h1_regional_common_support.py`
- syndrome diagnostic: `scripts/audit_pollination_syndrome_common_support.py`
