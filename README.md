# Island floral syndrome — v2

This repository supports **Chapter 1 / macroecological WHEN–WHERE analysis** of floral and reproductive trait filtering across island floras.

## Current question

> **When and where do floral/reproductive trait probabilities change along a mainland-distance/source-pool-accessibility gradient, and where do multivariate response vectors differ?**

GBIF occurrence data are treated as an **opportunistic, incompletely observed sample of realised island floras**, not a census. No-record islands are retained in a separate observation-process layer and are never treated as biological trait zeros.

Mainland distance is a composite geographic axis that may represent both dispersal limitation and changing accessibility to mainland/source species pools. It is not interpreted as a pure causal isolation treatment.

## Frozen result

Formal observation-robust workflow run: `32845980788`.

> **Within opportunistically observed island floras, floral and reproductive trait probabilities change systematically along the geographic/source-pool gradient in both northern mid-latitude and tropical regions. Their multivariate response vectors differ, and the result persists in native non-endemic assemblages. Current data remain confirmatorily unresolved for northern high-latitude and southern-extratropical floras.**

The headline survives the major observation-bias checks:

- **10/10** headline replications across canonical, capped-100, capped-50, capped-20 and equal-island information weighting × two status strata;
- **2/2** headline replications after response-specific direct-trait coverage adjustment;
- **6/6** across log1p, square-root and raw distance representations while retaining the same full island universe;
- **84/84** leave-one-spatial-block deletions in all-native flora and **84/84** in native-nonendemic flora.

The checkpoint therefore records `observation_robust_headline = true`.

See:

- [`docs/chapter1_when_where_frozen_result_20260825.md`](docs/chapter1_when_where_frozen_result_20260825.md) — canonical observation-robust result;
- [`docs/chapter1_next_acquisition_priority_20260825.md`](docs/chapter1_next_acquisition_priority_20260825.md) — outcome-blind next data-acquisition priorities;
- [`THESIS_CHAPTER_POSITIONING.md`](THESIS_CHAPTER_POSITIONING.md) — dissertation role and claim boundary;
- [`docs/chapter1_pollination_syndrome_concordance_20260825.md`](docs/chapter1_pollination_syndrome_concordance_20260825.md) — post-freeze Discussion-only pollination-syndrome audit.

## Three analysis layers

```text
Layer O — observation process
all 8,265 islands
→ where is flora recorded?
→ where is direct trait evidence available?

Layer T — trait-centric ecological surface
trait-resolved observed flora
→ P(trait state | geography, area, climate, context)

Layer W — WHEN / WHERE
→ within-context multivariate response-vector test
→ between-context response-vector comparison
→ floristic-status persistence
→ lineage / genus-composition safeguard
```

### Observation layer

All **8,265 islands** remain represented:

- no flora record: 3,760;
- flora recorded: 4,505;
- at least one direct Chapter 1 trait represented: 425;
- all three core direct-trait domains represented: 405.

Observation and trait-resolution probabilities are geographically structured, so they are modelled explicitly rather than assumed ignorable.

### Trait-centric layer

The ecological estimand is:

```text
P(trait state | directly trait-resolved observed realised flora, geography)
```

For example, current direct evidence does **not** support a simple “farther islands have more SC species” result: SC distance slopes are weak in both the northern-midlatitude and tropical all-native samples.

### WHEN / WHERE layer

Confirmatory response-vector tests support geographic trait filtering in:

- **northern mid-latitude** floras;
- **tropical** floras.

The same result persists in `native_nonendemic` assemblages, and the northern-midlatitude and tropical vectors differ directly in a formal multivariate contrast.

Northern high-latitude and southern-extratropical floras are **unresolved confirmatorily**, not demonstrated null regions. Southern extratropical floras have a pilot signal.

## Additional acquisition

Additional data must not be collected to rescue or strengthen the already-supported northern/tropical result.

The next acquisition objective is **regional testability**:

1. expand floristic-status resolution on already flora-recorded **northern high-latitude** and **southern extratropical** islands;
2. then prioritize direct **floral-form** evidence;
3. then direct **SI/SC** evidence;
4. colour is lower priority because current direct coverage is already comparatively high.

This order matters: southern extratropical has 317 flora-recorded islands but only 34 currently entering the native-status trait surface; northern high-latitude has 424 flora-recorded islands but only 12. Trait filling alone cannot close the >=50-island confirmatory gap without first widening status-resolved island coverage.

Acquisition priority must be outcome-blind. Preliminary trait values, effect directions and p-values must not enter the queue.

## Bombus and pollination-syndrome boundary

Bombus, bird, butterfly, moth, hawkmoth and other pollinator labels do **not** enter the primary Chapter 1 design.

They may be compared with the frozen response vectors for Discussion-level concordance or mismatch only. Trait vectors do not identify causal pollinator guilds.

## Thesis handoff

```text
Chapter 1 — island
WHEN / WHERE do trait probabilities respond to the geographic/source-pool gradient?
WHERE do multivariate responses differ?
        ↓
Chapter 2 — izu-core
WHY do those contexts generate different response architectures?
        ↓
Chapter 3 — shimahotarubukuro
WHAT phenotype axes actually diverge within one focal lineage?
```

## Repository layout

- `src/island_v2/` — reusable v2 data and analysis utilities
- `analysis/v2/` — statistical analysis scripts and execution notes
- `config/` — frozen contracts, ontology, and artifact locks
- `data/v2/` — external/staging/curated/template data layers
- `docs/` — scientific design, data policy, methods, and reproducibility notes
- `.github/workflows/` — active validation and canonical analysis workflows
- `legacy/v1/` — frozen v1 provenance only

## Reproducibility rule

GitHub Actions artifacts are temporary. Before submission, manuscript-critical inputs and outputs must be archived durably with checksums, with one canonical workflow identified for each main analysis and attrition reported from the full 8,265-island universe to every fitted ecological surface.
