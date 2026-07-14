# Manuscript submission contract

## Status

This document defines the repository surface that may support the manuscript. It is intentionally narrower than the full development history.

The repository contains many historical pilots, acquisition experiments, and alternative analyses. Their existence in git history does not make them part of the manuscript analysis. A result is manuscript-canonical only when it is listed here, uses a locked input artifact, and produces an auditable output manifest.

## Canonical scientific question

The manuscript tests whether island isolation is associated with changes in floral and reproductive composition, with a northern-midlatitude Bombus-channel hypothesis evaluated only where that channel is biologically interpretable.

The study is a global comparative analysis. It does not identify a historical pollinator-loss event and does not establish causal floral evolution from cross-sectional island data.

## Canonical analysis hierarchy

1. **Global composition analysis**
   - estimates regime-specific associations between isolation and retained floral/reproductive categories;
   - does not impose a Bombus mechanism outside its declared applicability domain;
   - preserves category composition whenever possible.

2. **Northern-midlatitude Bombus analysis**
   - uses only an analysis-ready Bombus predictor with explicit provenance;
   - a global all-species environmental maximum is diagnostic only and cannot be treated as Bombus absence or deficit;
   - environmental compatibility is not realized occurrence, abundance, visitation, or pollination service.

3. **Reproductive-assurance pathway analysis**
   - is an associational pathway analysis unless a formal causal mediation estimand, assumptions, common analysis support, and response-scale posterior contrast are implemented;
   - products of coefficients from separate nonlinear models must not be labelled causal indirect effects or total effects.

4. **Sensitivity analyses**
   - expanded genus/family inference and global fallback may be used only as explicitly labelled sensitivity analyses;
   - results based on `sensitivity_all` are not confirmatory evidence because global fallback values are model-imputed priors rather than species-level observations.

## Evidence-tier rule

The manuscript must report trait evidence resolution rather than hide it behind complete fill.

- **Confirmatory:** direct source-backed species evidence, with a separately declared broader direct-evidence track if justified.
- **Secondary robustness:** taxonomic inference at genus/family level.
- **Sensitivity only:** global fallback.

Every primary figure/table must state the evidence tier and the number of islands, species, and trait-resolved trials contributing to the estimate.

## Required model safeguards before submission

### Common support for pathway claims

Any coefficient-product or mediated-path comparison must be estimated on the same island set for all component equations. Equation-specific maximum support may be reported for separate descriptive models, but it cannot be used to construct a single indirect pathway estimand.

### Nonlinear-model mediation

For beta-binomial/logit models, multiplying logit-scale coefficients is not a natural indirect effect on the response scale. A manuscript mediation claim therefore requires posterior simulation or g-computation under declared exposure/mediator contrasts. Otherwise the output must be described as an exploratory coefficient-product pathway index.

### Spatial and lineage structure

A manuscript-canonical global or regional inference must include, or explicitly test sensitivity to:

- spatial dependence among islands; and
- inherited floristic/lineage composition or source-pool structure.

A model containing only island geography and climate does not by itself satisfy the source-pool/lineage safeguard defined in `docs/v2_pollination_regime_framework.md`.

### Isolation functional form

The isolation response must be predeclared and checked for leverage. Raw distance plus a quadratic term can be dominated by a small number of remote islands. The manuscript must report the observed distance distribution and compare the chosen functional form with a prespecified robust alternative such as log-distance or a smooth with constrained complexity.

### Zero-distance islands

`distance_to_continent_km <= 0` must be audited as a data/geometry class. Dropping these islands before modelling is acceptable only after reporting their count, geography, and reason for exclusion. They must not disappear silently from the analysis universe.

## Bombus environmental-niche contract

The retained estimator is the PR #112-style winsorized, standardized, ridge-regularized Mahalanobis ellipsoidal environmental niche model with boundary calibration and extrapolation diagnostics.

Before manuscript use, the real-data artifact must distinguish two things:

1. the frozen species universe used for reproducibility; and
2. the occurrence records actually used after the declared quality filters and spatial thinning.

Restoring an older canonical artifact is reproducible, but it is not evidence that newly added coordinate-uncertainty, year, or basis-of-record filters were applied. The final manuscript artifact must record counts before and after each filter and a checksum for the exact occurrence table used to fit the model.

## Reproducibility and archival rule

GitHub Actions artifacts with finite retention are not a permanent supplementary archive. Before submission, all manuscript-critical inputs and outputs must be deposited in a durable release/archive (for example a versioned GitHub release plus an external DOI-bearing repository where appropriate), with checksums recorded in the repository.

A fresh reader should be able to identify from the repository root:

- the exact manuscript release/tag;
- the locked input artifacts and checksums;
- one canonical command/workflow per main analysis;
- the software environment;
- the evidence-tier definition;
- the attrition table from the frozen island universe to each fitted model; and
- the files used for each manuscript figure and table.

## Noncanonical material

The following are development history rather than manuscript methods unless explicitly promoted by a later contract change:

- trait-source scouting and free-bulk-source pilots;
- validation-pilot and core-pilot acquisition workflows;
- EOL/TraitBank inventory experiments;
- superseded M0–M4 workflow variants;
- old v1/v2 bridge analyses;
- alternative Bayesian replications superseded by the engine-specific design.

Historical code may remain in git history, but it should not appear beside canonical workflows in the active submission surface.
