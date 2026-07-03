# v2: global island floral-channel architecture

## Status and boundary

`v1-freeze` is retained as a historical exploratory baseline. v2 is a new analysis program and must not import v1 scripts, derived trait labels, or v1 model outputs as analysis inputs.

The v2 objective is to test how island geography and the functional availability of bumblebees are associated with **floral-signal**, **floral-architecture**, and **reproductive-assurance** trait composition, while separating competing biological channels and the observational data process.

This architecture adapts the identification discipline of `campanula-channel-identification` to a global comparative setting. It does **not** claim that global flora composition alone identifies evolutionary selection.

## The cross-scale outcome model

For a species trait state `z` in an island regime:

```text
W(z) = F(z) × E(z)
```

- `F(z)`: local reproductive contribution, including outcrossed and selfed viable seed production and paternal contribution where measured.
- `E(z)`: arrival, establishment, persistence, and observation conditional on viable seed output.
- `W(z)`: observed island-flora membership or trait composition.

In the global database, `W` is observed only through island-by-species occurrence evidence. It is therefore not valid to interpret a trait association with `W` as direct selection on `F` unless intermediate observations or defensible external calibration are available.

## Predeclared competing channels

### C0: geography and source-pool channel

```text
source pool + island area + isolation + climate + biogeographic region
    -> island flora trait composition
```

This is the baseline explanation. Every biological channel is evaluated against it.

### C1: direct bumblebee-functional channel

```text
bumblebee functional deficit
    -> floral signal / floral architecture composition
```

This represents trait sorting or selection consistent with changed attraction, handling, pollen placement, or pollen export. It requires that the relevant floral trait has a declared functional link to bumblebee use; visitor counts alone are not sufficient evidence.

### C2: reproductive-assurance / selfing channel

```text
bumblebee functional deficit
    -> reduced reliable outcross service
    -> reproductive assurance traits
    -> floral signal / floral architecture composition
```

Reproductive-assurance traits include self-incompatibility, autonomous selfing, mating system, herkogamy, dichogamy, cleistogamy, and sex system. Self-compatibility alone must not be interpreted as autonomous selfing.

### C3: joint channel

```text
bumblebee functional deficit
    -> direct floral channel
    -> reproductive-assurance channel
    -> floral trait composition
```

C3 is supported only if the direct association remains after conditioning on the measured reproductive-assurance domain, while uncertainty and proxy status are preserved.

### C4: establishment / reachability channel

```text
traits, dispersal, source pool, island accessibility, habitat and persistence
    -> E(z)
    -> island-flora composition
```

C4 prevents an apparent pollinator explanation from absorbing colonisation, establishment, extinction, or taxonomic source-pool effects.

### C5: observational channel

```text
survey effort + GBIF/data-source structure + spatial coverage + taxonomic coverage
    -> observed W(z)
```

C5 is not biology. It is modelled and reported separately so that an apparent island trait pattern is not an artifact of records, accessible coastlines, specimen concentration, or data-source composition.

## What the global dataset can and cannot claim

The global analysis can estimate whether island trait composition is more consistent with C0, C1, C2, C3, C4, or a combination, conditional on observed covariates and data quality.

It cannot, by itself, establish that bumblebee absence historically selected a floral phenotype. Strong evolutionary claims require independent evidence such as repeated island-mainland contrasts in matched lineages, functional observations, paternity or pollen-transfer data, and recruitment measurements. The Campanula/Izu field program provides that complementary scale.

## Trait domains retained without premature binarisation

### Floral signal

- primary colour: white, green/brown/inconspicuous, yellow/orange, red/pink, blue/purple, multicoloured/variable, unresolved
- colour variability and patterning where evidence is available
- nectar guides / contrast only when the source specifies them

### Floral architecture and access

- symmetry
- floral form and openness
- tube presence and tube-depth class
- inflorescence display class
- flower-size measurements or ordinal size class
- composite head / brush-puff / open radial / tubular / zygomorphic forms
- reward location and nectar accessibility when directly supported

### Reproductive-assurance domain

- self-incompatibility status
- autonomous-selfing capacity
- mating system
- herkogamy
- dichogamy
- cleistogamy
- sex system
- pollen limitation, seed set, paternity, and inbreeding-depression evidence where available

Binary `plain/conspicuous`, `generalized/specialized`, and `SC/SI` contrasts are retained only as transparent secondary summaries for continuity with v1 and sensitivity analyses.

## Functional availability of Bombus

The v2 predictor is not a binary GBIF-absence variable. It has separate components:

1. `environmental_compatibility`: whether island conditions overlap the environmental envelope of candidate Bombus species;
2. `occurrence_support`: quality-filtered, effort-aware evidence that Bombus occurs in or near the island system;
3. `functional_service_proxy`: evidence relevant to activity season, richness, abundance, or documented functional use when available;
4. `proxy_quality`: whether the preceding components are direct observations, calibrated proxies, or unresolved.

The old PDI is treated only as a candidate environmental proxy for the first component, not as direct proof of pollinator service.

## Analysis sequence

1. Build a polygon-exact, provenance-tracked island-by-species occurrence database and a separate observation-effort table.
2. Build a source-cited, LLM-assisted trait evidence database.
3. Estimate geographic/source-pool and observational baselines before introducing bumblebee predictors.
4. Model floral signal, architecture, and reproductive-assurance composition as separate multistate or continuous domains.
5. Compare C0--C4 using predeclared model comparisons and report C5 diagnostics separately.
6. Test whether direct bumblebee-channel associations remain after reproductive-assurance traits are represented.
7. Use Campanula/Izu data to evaluate the intermediate measurements required for stronger channel claims.

## Required outputs

- a channel claim table stating each pathway, its required intermediate observable, its current evidence tier, and its falsification target;
- data-process diagnostics per island;
- trait evidence coverage maps and source-type summaries;
- domain-specific composition models before any composite syndrome model;
- a final cross-scale interpretation separating global compatibility from field-based causal evidence.
