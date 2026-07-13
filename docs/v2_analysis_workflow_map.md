# v2 analysis workflow map

This document defines the scientific role of the v2 analysis workflows. The registry in `config/v2_workflow_registry.yml` is the machine-readable source of truth.

## Canonical manuscript path

```text
frozen island universe
        ↓
evidence-tiered species trait master (primary = species-direct)
        ↓
island trait evidence aggregation
        ↓
Bombus channel components + locked geographic/climate covariates
        ↓
validated analysis-ready input
        ↓
PRIMARY: INLA M0–M4
```

The only canonical primary-inference workflow is:

- `.github/workflows/run-v2-inla-m0-m4-main.yml`

Primary inference uses the `primary` evidence tier. `broad` and `sensitivity_all` are robustness tiers and must not silently replace the primary tier.

## Model ladder

- **M0**: geographic filter — isolation × island area + climate.
- **M1**: adds the Bombus channel to floral-trait components.
- **M2**: adds reproductive assurance.
- **M3**: combines Bombus and reproductive-assurance paths, including indirect effects.
- **M4**: category-preserving replacement/regime analysis.

The executable model contract is `config/v2_analysis_contract.yml`.

## Robustness analyses

Robustness workflows test whether the primary conclusion survives changes in evidence depth, spatial validation and trait uncertainty. They are not alternate definitions of the primary analysis.

- Evidence depth: `broad`, `sensitivity_all`
- Spatial robustness: spatial CV and common-support CV
- Trait uncertainty: probabilistic trait uncertainty and inference validation
- Alternative model engine: brms replication

## Replication

`.github/workflows/run-v2-bayesian-m0-m4-main.yml` is retained as a Bayesian engine replication using brms. It should use the same primary evidence contract when used for manuscript replication and should not be described as a second independent primary analysis.

## Legacy and development workflows

Older M0–M4 workflows are retained for provenance. Their existence records the method-development history; it does not make them part of the manuscript's canonical inferential path. Promotion of any legacy workflow requires an explicit edit to `config/v2_workflow_registry.yml`.

## Scientific guardrails

1. Unknown evidence is not a negative observation.
2. `primary` means species-direct evidence only.
3. `sensitivity_all` includes fallback-derived values and is sensitivity-only.
4. Bombus environmental compatibility must not be described as observed absence.
5. Input artifacts must be locked and digest-verified.
6. Producer and consumer schemas must be validated before a full analysis run.
7. The frozen island universe is defined by its manifest, not by an unexplained hard-coded row count.

## Intended repository shape

```text
Data production → frozen inputs → one primary analysis
                                ├─ evidence-tier sensitivity
                                ├─ spatial robustness
                                ├─ trait-uncertainty robustness
                                └─ alternative-engine replication
```

New workflows should be added only when they represent a genuinely new scientific robustness dimension. Small implementation experiments should be tests or local scripts rather than new manuscript-level workflows.
