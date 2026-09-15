# Chapter 1 P2 — component non-concordance result

Status: post-baseline robustness result. This does not replace the frozen H2 estimand or retroactively alter H2 thresholds.

## Question

Does the classic floral-island response behave as one coupled response vector, or does the North–Tropical response difference persist after forcing the two primary components onto the same observation support?

The formal direct comparison is **northern_midlatitude versus tropical within the `analysis_regime` layer**. `Palearctic` belongs to the separate `biogeographic_realm` layer and is not interchangeable with northern_midlatitude in a direct between-context test.

## Frozen execution

- source PR142 run: `34232450884`
- source artifact: `10058653212`
- source digest: `sha256:b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`
- P2 run: `34945548775`
- P2 artifact: `10387261614`
- P2 digest: `sha256:f2038c37459986c5b4abc57265fb69026abb34d1c46eb4733acc504588bf2491`
- paired spatial-block bootstrap: 2,000 draws

The first two attempts stopped before outcome opening at test/lint. The successful run used the same frozen scientific contract; repairs were only floating-point test tolerance and removal of one unused import.

## P2a — common-island support

Both primary axes were restricted to exactly the same islands within each profile before refitting.

Primary direct-only × native-nonendemic result:

- common islands: **348**
- spatial blocks: **78**
- direct joint North–Tropical vector-difference `p = 0.0009473`
- northern-midlatitude vector: `(accessibility, assurance) = (0.02544, 0.03743)`
- tropical vector: `(-0.10394, 0.17466)`

The direct joint vector difference therefore survives forcing both axes onto a common island denominator.

The same direction of result appears in the other three profiles:

- all-analysis/all-native: `p = 0.1767`
- all-analysis/NNE: `p = 0.000229`
- direct/all-native: `p = 0.01724`
- direct/NNE: `p = 0.000947`

The direct-only primary profiles retain the joint difference; the all-analysis all-native cell is weaker and is not used to upgrade the claim.

## P2a geometry — stronger non-collinearity claim does not pass

For the primary direct/NNE profile:

- determinant = **0.007458**
- paired block-bootstrap 95% CI = **[-0.003193, 0.015991]**
- angle = **67.14°**
- angle 95% interval = **[12.27°, 169.18°]**

All four determinant intervals include zero. Therefore the data do **not** establish the stronger geometric statement that the two context vectors are demonstrably non-collinear rather than noisy scaled/rotated alternatives.

Safe conclusion:

> the same-layer North–Tropical joint response vectors differ, but their exact geometric non-concordance is imprecisely estimated.

## P2b — same co-observed species denominator

A stronger sensitivity analysis rebuilt both primary axes from only species with finite scores for both `generalized_accessible` and `selfing_core`.

Direct-only:

- co-observed species: **853**
- NNE islands: **280**
- direct North–Tropical vector-difference `p = 0.005167`
- northern vector: `(0.03668, -0.01827)`
- tropical vector: `(-0.06683, 0.20732)`

Thus differing species denominators are not sufficient to explain away the direct North–Tropical joint vector difference.

However, this common-species restriction changes the estimand and individual within-context component slopes need not reproduce the original H2 point estimates. In particular, it should not be used to claim a universally positive northern reproductive-assurance component.

## P2c — realm boundary

Frozen Palearctic–Neotropical direct vector tests remain unsupported:

- all-analysis/all-native `p = 0.0784`
- all-analysis/NNE `p = 0.1190`
- direct/all-native `p = 0.2838`
- direct/NNE `p = 0.3958`

Therefore manuscript language must distinguish:

1. **Palearctic** — a strong within-context branch and the focal H3/P1 genus-structure result;
2. **North–Tropical** — the formal same-layer direct H2 contrast;
3. **Palearctic–Neotropical** — not supported as a formal direct contrast under the frozen result.

## Integrated decision

P2 strengthens H2 in one important respect and narrows it in another.

Strengthened:

> The North–Tropical joint response-vector difference is not a by-product of using different islands for the two components, and it persists under a much stricter common-species denominator.

Narrowed:

> The stronger claim that context vectors are demonstrably non-collinear is not established, and Palearctic-versus-tropical wording must not be presented as the formal direct H2 test.

This is the P2 claim ceiling for the submission surface.
