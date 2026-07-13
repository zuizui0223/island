# v2 engine-specific analysis contract

## Central hypothesis

Island isolation does not necessarily generate a universal floral syndrome. Its consequences depend on whether the dominant regional pollination channel is disrupted. In northern-midlatitude systems, isolation may generate a Bombus deficit that changes floral composition through two biologically distinct paths:

1. a direct path from Bombus deficit to floral composition, representing altered pollinator-mediated selection; and
2. an indirect path from Bombus deficit through self-compatibility to floral composition, representing reproductive assurance.

Tropical and southern-extratropical regimes are falsification domains for the universality of this northern pathway. Birds, butterflies, moths, and other bees remain discussion-level functional-replacement explanations until comparable island-level channel measurements are available.

## The engines answer different questions

The Bayesian and INLA analyses are not duplicate implementations and are not combined into a single model-selection contest.

### Bayesian pathway inference

The brms workflow is restricted to broad binary outcomes:

- self-compatible versus self-incompatible;
- plain versus non-plain colour composition;
- open/generalized versus other floral forms.

Its purpose is to estimate the pathway structure:

1. Bombus deficit to self-compatibility;
2. direct Bombus-deficit effects on binary floral responses;
3. indirect Bombus-deficit effects through self-compatibility;
4. total direct plus indirect effects;
5. regime-specific isolation slopes used to falsify a universal island syndrome.

Bayesian output is interpreted at the syndrome/pathway level only. It is not used to decide which individual colour or floral-form category changed.

### INLA category decomposition

INLA owns category-preserving inference. It estimates retained flower categories separately and determines which categories increase, decrease, oppose one another, or cancel under binary aggregation.

Retained colour categories:

- plain;
- yellow/orange;
- red/pink;
- blue/purple.

Retained floral-form categories:

- open/generalized;
- tubular/trumpet.

INLA output is interpreted at the category-composition level. It is not treated as a second estimate of the same binary pathway.

## Support and filtering

- Islands with `distance_to_continent_km <= 0` are audited and excluded before logarithmic transformation.
- Each response equation uses its maximum available island support.
- Missing floral-form evidence does not remove an otherwise usable island from SC or colour equations.
- The current exploratory implementation uses `sensitivity_all`; evidence-tier robustness remains a separate question.

## Falsification design

The northern-midlatitude Bombus pathway is fitted only in the Bombus-primary regime. Outside that regime, the analysis tests whether isolation-associated shifts in SC, plain colour, and generalized form weaken, disappear, or reverse.

This is a falsification test of a universal isolation syndrome, not an attempt to impose a Bombus causal mechanism on tropical or southern-extratropical islands.

## Interpretation rules

- Bayesian results answer whether the broad direct and SC-mediated pathways are supported.
- INLA results answer which retained categories produce or contradict the broad binary pattern.
- Agreement is not counted as two independent confirmations because both engines use the same underlying island-trait data.
- A Bayesian binary effect without coherent category-level decomposition is reported as an aggregate effect.
- An INLA category effect that cancels in the binary analysis is reported as category turnover hidden by aggregation.
- Alternative pollinator mobility is a biological explanation, not a measured causal path in the current models.