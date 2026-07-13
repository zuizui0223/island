# v2 analysis workflow map

## Canonical current analysis

The canonical v2 route is:

- workflow: `.github/workflows/run-v2-inla-category-preserving-north.yml`
- engine: INLA
- current evidence tier: `sensitivity_all`
- planned confirmatory tier: `primary`

The northern-midlatitude mechanism is evaluated as:

1. isolation and geography/climate predict Bombus channel deficit;
2. Bombus deficit predicts self-compatibility;
3. flower colour and floral form respond through a direct Bombus-linked path and an indirect self-compatibility-mediated path;
4. direct, indirect, and total effects are reported separately.

Tropical and southern-extratropical islands are falsification domains. The global comparison tests whether isolation-associated slopes weaken, disappear, or reverse outside the northern-midlatitude domain. Alternative pollinator guilds are not required primary-model covariates; birds, butterflies, moths, bats, other bees, and generalist insects remain biological interpretations of counter-patterns.

## Data support rules

- Preserve the frozen 8,265-island audit universe.
- Exclude `distance_to_continent_km <= 0` from fitted model support.
- Use the maximum available support for each equation rather than one colour/form/SC complete-case intersection.
- Preserve flower-colour and floral-form category resolution in INLA.
- Treat unknown evidence as unknown, not as a negative observation.

## Engine roles

INLA owns category-preserving colour/form inference, the northern direct/indirect decomposition, and the tropical/southern falsification analysis.

brms is noncanonical replication only. It may be used for compact scalar path checks, but it must not become a second competing category-analysis route.

## Evidence tiers

The current maximum-data run uses `sensitivity_all` to inspect support and model behaviour. The planned confirmatory run must use the `primary` species-direct tier. Neither tier may silently replace the other.

## Legacy workflows

Older all-data M0-M4 and Bayesian category routes remain only for provenance and method-development history unless explicitly promoted in both `config/v2_workflow_registry.yml` and `config/v2_analysis_contract.yml`.
