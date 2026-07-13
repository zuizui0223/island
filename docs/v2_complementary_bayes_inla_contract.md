# v2 complementary Bayesian and INLA analysis contract

## Scientific hypothesis

Island isolation does not necessarily generate a universal floral syndrome. Its consequences depend on whether the dominant regional pollination channel is disrupted. In northern-midlatitude systems, isolation may generate a Bombus deficit that affects floral composition directly and indirectly through self-compatibility. Tropical and southern-extratropical regimes are used as falsification domains for the universality of this northern pathway.

Birds, butterflies, moths, and other bees are retained as biological explanations for functional replacement outside the northern regime. They are not included as primary causal predictors in the current models because the available evidence does not support consistent island-level measurement of those channels.

## Complementary engine roles

### Bayesian binary pathway analysis

The brms workflow estimates broad syndrome-level responses:

- self-compatible versus self-incompatible;
- plain versus non-plain colour composition;
- open/generalized versus other floral forms.

It estimates:

1. Bombus deficit to self-compatibility;
2. the direct Bombus-deficit effect on each binary floral response;
3. the indirect Bombus-deficit effect through self-compatibility;
4. the total direct plus indirect effect;
5. isolation-by-regime slopes as a global falsification test.

Each equation uses its maximum available island support. Islands with non-positive distance to the continent are audited and excluded before logarithmic transformation.

### INLA category-level decomposition

The INLA workflow retains the available flower categories rather than collapsing them into a single binary syndrome. It determines which colour and floral-form categories account for any broad pattern observed in the Bayesian analysis.

The retained colour categories are:

- plain;
- yellow/orange;
- red/pink;
- blue/purple.

The retained floral-form categories are:

- open/generalized;
- tubular/trumpet.

Category availability is checked explicitly before fitting. Each outcome uses its own maximum support.

## Interpretation rules

The two engines are complementary, not competing primary and replication analyses.

- Bayesian results answer whether a broad binary island syndrome and its direct/SC-mediated pathways are supported.
- INLA results answer which retained categories drive, oppose, or complicate that syndrome.
- Agreement between the two engines is not counted as two independent confirmations because both use the same underlying island-trait data.
- A binary Bayesian effect without a coherent category-level INLA decomposition must be described as an aggregate result rather than a universal category shift.
- A category-level INLA effect that cancels in the binary analysis is interpreted as category turnover hidden by aggregation.

## Falsification boundary

The northern Bombus pathway is not imposed globally. Tropical and southern-extratropical islands test whether isolation slopes weaken, disappear, or reverse outside the Bombus-primary regime. Functional replacement by highly mobile pollinators remains a discussion-level explanation unless direct comparable channel data become available.
