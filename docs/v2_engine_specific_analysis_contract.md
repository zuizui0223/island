# v2 engine-specific analysis contract

## Central hypothesis

Island isolation does not necessarily generate a universal floral syndrome. Its consequences depend on whether the dominant regional pollination channel is disrupted. In northern-midlatitude systems, isolation may generate a Bombus deficit that changes floral composition through two biologically distinct paths:

1. a direct path from Bombus deficit to floral composition, representing altered pollinator-mediated selection; and
2. an indirect path from Bombus deficit through self-compatibility to floral composition, representing reproductive assurance.

Tropical and southern-extratropical regimes are falsification domains for the universality of this northern pathway. Birds, butterflies, moths, and other bees remain discussion-level functional-replacement explanations until comparable island-level channel measurements are available.

## The engines answer different questions

The Bayesian and INLA analyses are not duplicate implementations and are not combined into one model-selection contest.

### Bayesian pathway inference

The brms workflow is restricted to broad binary outcomes:

- self-compatible versus self-incompatible;
- plain versus non-plain colour composition;
- open/generalized versus all other floral forms.

Its purpose is to estimate:

1. Bombus deficit to self-compatibility;
2. direct Bombus-deficit effects on binary floral responses;
3. indirect Bombus-deficit effects through self-compatibility;
4. total direct plus indirect effects;
5. regime-specific isolation slopes used to falsify a universal island syndrome.

Bayesian output is interpreted at the syndrome/pathway level only. It is not used to identify which individual category changed.

### INLA category decomposition

INLA owns category-preserving inference for all observed categories retained by PR #115:

- colour: plain, yellow/orange, red/pink, blue/purple;
- form: open/generalized, tubular/trumpet, zygomorphic/specialized, composite/brush.

INLA determines which categories increase, decrease, oppose one another, or cancel under binary aggregation. It is not treated as a second estimate of the same binary pathway.

## Falsification design

The northern Bombus pathway is fitted only in the northern-midlatitude focal regime. Tropical and southern-extratropical analyses estimate isolation slopes without imposing Bombus as their causal mechanism. Weakening, disappearance, or reversal outside the north falsifies a universal island syndrome.

## Support and filtering

- Islands with `distance_to_continent_km <= 0` are audited and excluded before log transformation.
- Each equation uses its maximum available island support.
- Missing form evidence does not remove an otherwise usable island from SC or colour equations.
- The current exploratory implementation uses `sensitivity_all`; evidence-tier robustness remains separate.

## Interpretation rules

- Bayesian results answer whether broad direct and SC-mediated pathways are supported.
- INLA results answer which categories produce or contradict the aggregate pattern.
- Agreement is not counted as independent replication because both use the same underlying data.
- A binary effect without coherent category decomposition is reported as an aggregate effect.
- Category turnover hidden by binary cancellation is reported explicitly.
- Alternative-pollinator mobility remains a biological explanation rather than a measured causal path.
