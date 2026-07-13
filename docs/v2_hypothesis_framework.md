# v2 hypothesis framework: island pollination-function and floral phenotype reassembly

## Core question

How does island isolation reorganize pollination function, reproductive assurance, and the joint composition of flower colour and floral form, and when is that reorganization specifically associated with the northern Bombus channel?

The framework does not assume a universal island selfing syndrome or a universal shift toward plain or open flowers. Those are testable special cases within a broader compositional-reassembly hypothesis.

## H1. Pollination-function filter

Island isolation reorganizes whole-flora pollination-mode composition. Animal-pollinated and wind/mixed components are modelled as whole-flora outcomes and may respond in opposite directions.

Primary test:

`pollination-mode composition ~ isolation + island area + climate + spatial structure`

Interpretation target: whether isolation changes the functional composition of the flora before any Bombus-specific mechanism is invoked.

## H2. Northern Bombus channel

Within the northern-midlatitude domain, isolation is expected to increase Bombus deficit after accounting for island area, climate, and spatial structure.

Primary test:

`Bombus deficit ~ isolation + island area + climate + spatial structure`

Bombus is a domain-specific mechanism, not a global universal driver.

## H3. Reproductive-assurance partition

Self-compatibility may respond independently to Bombus deficit and to whole-flora pollination-mode composition. A Bombus-mediated selfing pathway is therefore tested rather than assumed.

Primary test:

`SC composition ~ Bombus deficit + wind/mixed composition + isolation + island area + climate + spatial structure`

The indirect pathway `Bombus deficit -> SC -> floral phenotype composition` is retained as a piecewise mediation test. Failure of the Bombus-to-SC path does not invalidate the wider island reassembly hypothesis.

## H4. Joint floral-phenotype reallocation

Isolation and pollination-channel change reorganize the joint composition of flower colour and floral form. The primary response is therefore a colour-by-form composition, not separate independent binary traits.

Primary fine-grained promoted colours:

- white
- blue
- yellow
- red
- green
- cream
- rare_other

Primary fine-grained promoted forms:

- tubular
- composite_head
- star
- bell_campanulate
- open_radial
- trumpet_funnel
- rare_other

The cross-classification yields at most 49 colour-by-form cells. Only globally observed cells are fitted, while every resolved low-support state is retained through `rare_other` so that island-level joint counts still exhaust the joint trial total.

Primary model ladder:

- FJ0: geography and climate
- FJ1: FJ0 + Bombus deficit
- FJ2: FJ0 + self-compatibility + wind/mixed composition
- FJ3: FJ0 + Bombus deficit + self-compatibility + wind/mixed composition

The 4-colour by 4-form joint model remains a coarse sensitivity analysis. The individual one-vs-rest binomial models remain diagnostic only.

## H5. Regime contingency

Isolation-associated reproductive and floral responses are not assumed to be globally universal. Northern midlatitudes, northern high latitudes, tropics, and southern extratropics are compared explicitly through isolation-by-regime interactions.

Primary test:

`outcome ~ isolation * regime + island area + climate + spatial structure`

The key result may be weakening, disappearance, or reversal of slopes outside the northern-midlatitude domain.

## Evidence hierarchy

1. Current exploratory analysis: `sensitivity_all` evidence tier.
2. Confirmatory analysis: `primary` evidence tier.
3. Primary flower-phenotype inference: promoted fine joint colour-by-form composition.
4. Coarse 4x4 joint composition: sensitivity analysis.
5. Separate category binomial models: diagnostic and interpretive sensitivity analyses only.

## Interpretation guardrails

- A conditional Bombus coefficient is an island-level association, not proof of direct evolutionary selection.
- `open_radial` is a morphological accessibility category, not a direct measurement of pollinator-network generalization.
- Self-compatibility is not equivalent to autonomous selfing or realized selfing rate.
- Alternative pollinator guilds are descriptive global outcomes and biological interpretation aids, not primary causal covariates.
- Species-level conflicting trait records remain unresolved rather than being chosen arbitrarily.
