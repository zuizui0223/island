# Chapter 1 axis mapping before response-geometry analysis

## Decision

Chapter 1 has three distinct layers that must not be treated as interchangeable:

1. **measurement axes** in the frozen trait snapshot;
2. **broad scalar contrasts** already defined before the response-geometry extension;
3. **derived multivariate response axes** used by H2/H3 and secondary architecture analyses.

The response-geometry extension must preserve these distinctions. It does not redefine the trait database after inspecting nonlinear fits.

## Frozen measurement axes

The canonical species-by-axis snapshot contains exactly three axes:

| measurement axis | component traits in `chapter1_trait_snapshot.py` | role |
|---|---|---|
| `flower_colour` | `flower_primary_color` | floral colour state |
| `floral_structural_complexity` | `floral_form`, `floral_symmetry`, `tube_depth_class`, `flower_size_class`, `inflorescence_display` | floral access / structural architecture |
| `reproductive_assurance` | `self_incompatibility`, `autonomous_selfing_capacity`, `mating_system`, `cleistogamy` | reproductive-assurance evidence |

These are database-coverage axes. They are **not** three fitted continuous outcomes.

## Existing broad scalar contrasts

The frozen PR138 broad-outcome layer already supplies one directional scalar contrast from each measurement domain:

| measurement domain | frozen broad contrast | positive direction | interpretation ceiling |
|---|---|---|---|
| colour | `plain_colour` | white / green-brown-inconspicuous relative to yellow-orange / red-pink / blue-purple | colour-composition contrast, not attraction intensity |
| structure | `generalized_form` | open-radial / brush-puff / composite-head relative to restricted forms | form-only accessibility contrast, not the full structural axis |
| reproduction | `self_compatibility` | SC relative to SI | one reproductive-assurance component, not autonomous selfing or realized mating system |

These three contrasts are used for the **first response-geometry identifiability audit** because all three have the same observed island-level count/share representation and the same pre-existing genus-fixed null infrastructure. This choice is made before fitting any nonlinear response shape.

They must not be relabelled as complete measurements of the three database axes.

## Current H2/H3 derived axes

The current primary plant-response vector is two-dimensional:

- `accessibility_generalization` -> `generalized_accessible`;
- `reproductive_assurance` -> `selfing_core`.

`generalized_accessible` combines floral form, symmetry and tube depth. `selfing_core` combines self-incompatibility, mating system and autonomous-selfing evidence while excluding floral attraction traits.

**Flower colour is not a component of the current primary H2 `universal_plant_response` vector.** Colour enters the broad PR138 contrasts and the secondary sampled architecture templates, where it is deliberately downweighted relative to shape/tube traits.

Therefore any manuscript or Figure 1 that presents H2 as a direct three-axis test of colour + structure + reproductive assurance is incorrect.

## Secondary architecture layer

The sampled `large_bee_like`, `butterfly_like`, and `bird_like` templates combine colour and structural traits. V4 decomposes them into a shared architecture factor plus template-specific residuals. These are plant-architecture concordances, not observed pollinator identities.

The geometry extension may later be replicated on this V4 response family, but it must not substitute guild labels for the three raw measurement axes.

## Geometry interpretation

Response geometry adds a new descriptive dimension:

- `G0_flat`: no systematic response along the exposure;
- `G1_cline`: smooth monotonic change;
- `G2_step`: level shift at one breakpoint;
- `G3_hinge`: continuous slope change at one breakpoint;
- `G4_reversal`: an interior turning point (hump or U shape).

A selected geometry is **not a mechanism label**. In particular:

- a step does not prove pollinator loss;
- a cline does not prove reproductive-assurance selection;
- a reversal does not prove functional replacement;
- a genus-adjusted residual does not prove within-lineage evolution.

Mechanism compatibility is assessed only after the geometry and taxonomic-depth results are known.

## Chapter 1 / Chapter 2 handoff

Chapter 1 asks whether the **direction, geometry and taxonomic depth** of floral/reproductive responses vary among biogeographic contexts.

Chapter 2 (`izu-core`) asks whether a focal interaction system can discriminate the functional processes generating those response shapes. This restores the original Izu `cline versus threshold` logic without importing its mechanism labels into the global Chapter 1 analysis.
