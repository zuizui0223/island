# v2: global island floral-syndrome analysis — channel-audit architecture

## Status and boundary

`v1-freeze` preserves the exploratory draft and its first complete results. v2 is rebuilt from new data products and new scripts; it must not import v1 scripts, derived trait labels, or v1 model outputs as analysis inputs.

**The scientific aim does not change from v1.** v2 asks whether island isolation and the functional deficit of bumblebees help explain island-level variation in floral and reproductive traits, and whether the association is better represented by:

```text
isolation / geography
    -> bumblebee environmental or functional deficit
    -> reproductive assurance / selfing
    -> floral-trait composition
```

with a possible additional direct route from bumblebee deficit to floral traits.

v2 does not adopt the Campanula/Izu project as a second empirical question, and it does not turn the global project into a full life-history experiment. It borrows only its **channel-identification discipline**: specify what each dataset can identify, distinguish direct measurements from proxies, compare predeclared competing explanations, and avoid treating a terminal pattern as proof of a particular mechanism.

## What stays from v1; what changes in v2

### Retained scientific core

1. Island geography is the starting context.
2. Bumblebee limitation or environmental mismatch is the focal pollinator hypothesis.
3. Reproductive assurance and selfing are the focal indirect mechanism.
4. Direct bumblebee-to-flower associations and selfing-mediated associations are compared explicitly.
5. Northern-temperate islands remain the expected strongest test region; tropical and southern regions test regional contingency rather than universal syndrome.

### Rebuilt implementation

1. Island floras are rebuilt using polygon-exact occurrence evidence, source-pool information, and observation-process diagnostics.
2. The trait database is rebuilt as a source-cited, LLM-assisted evidence system rather than a filled categorical table.
3. Floral signal, floral architecture, and reproductive traits remain multistate or continuous wherever possible; v1 binary contrasts become secondary summaries.
4. Geography, establishment/reachability, and observation bias are treated as required alternative explanations before biological interpretation.
5. The old PDI is treated as one candidate environmental proxy, not as direct evidence of bumblebee service.

## Channel-audit lens borrowed from Campanula/Izu

For a species trait state `z` in an island regime, the Campanula project writes:

```text
W(z) = F(z) × E(z)
```

For the global analysis this is **not a new fitted model**. It is an interpretation guardrail:

- `F(z)`: local reproductive contribution, including outcrossed and selfed viable seed production and paternal contribution where directly measured;
- `E(z)`: arrival, establishment, persistence, and reachability conditional on viable seed output;
- `W(z)`: observed island-flora membership or island trait composition.

The global dataset mainly observes `W`, through island-by-species occurrence evidence. Therefore a global association between bumblebee deficit and a trait must not automatically be called direct pollinator-mediated selection on `F`. It may be consistent with direct floral filtering, selfing-mediated change, source-pool structure, establishment, or observation bias. The purpose of the audit is to state those alternatives and ask which remain compatible with the data.

## Predeclared explanations for the v1 question

### M0: geography / source-pool baseline

```text
source pool + island area + isolation + climate + biogeographic region
    -> island floral and reproductive trait composition
```

This is the necessary baseline, not a side analysis.

### M1: direct bumblebee-functional route

```text
bumblebee functional deficit
    -> floral-signal and floral-architecture composition
```

This tests the v1 possibility that loss of an efficient functional pollinator channel is associated with floral-trait change independently of reproductive assurance. It is interpreted as trait sorting or a pattern consistent with direct selection, unless field or lineage evidence identifies the intermediate process.

### M2: reproductive-assurance / selfing-mediated route

```text
bumblebee functional deficit
    -> reduced reliable outcross service
    -> reproductive-assurance traits
    -> floral-signal and floral-architecture composition
```

Self-compatibility alone is not autonomous selfing. The reproductive-assurance domain must preserve uncertainty about the distinction.

### M3: joint direct-plus-mediated route

```text
bumblebee functional deficit
    -> floral traits
    -> reproductive assurance
    -> floral traits
```

M3 is supported only when the direct association remains after the measured reproductive-assurance domain is represented, while proxy status and uncertainty are retained.

### Required alternative explanations, not new biological aims

```text
establishment / reachability:
traits + dispersal + source pool + habitat + persistence
    -> island-flora composition

observation process:
survey effort + GBIF/data-source structure + spatial coverage + taxonomic coverage
    -> observed island-flora composition
```

These are included so that a direct or selfing-mediated bumblebee explanation is not an artifact of colonisation history, extinction, accessible-coast sampling, specimen concentration, or database structure.

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

Binary `plain/conspicuous`, `generalized/specialized`, and `SC/SI` contrasts remain transparent secondary summaries for continuity with v1 and sensitivity analyses; they are not the sole v2 outcomes.

## Functional availability of Bombus

The v2 predictor is not a binary GBIF-absence variable. It has separate components:

1. `environmental_compatibility`: whether island conditions overlap the environmental envelope of candidate Bombus species;
2. `occurrence_support`: quality-filtered, effort-aware evidence that Bombus occurs in or near the island system;
3. `functional_service_proxy`: evidence relevant to activity season, richness, abundance, or documented functional use when available;
4. `proxy_quality`: whether the preceding components are direct observations, calibrated proxies, or unresolved.

The former PDI is one candidate proxy for component 1, not direct proof of pollinator service.

## Analysis sequence

1. Build a polygon-exact, provenance-tracked island-by-species occurrence database and a separate observation-effort table.
2. Build a source-cited, LLM-assisted trait evidence database.
3. Estimate M0 and observation-process diagnostics before introducing bumblebee predictors.
4. Model floral signal, floral architecture, and reproductive assurance as separate multistate or continuous domains.
5. Compare M0--M3, including direct versus selfing-mediated routes, with predeclared diagnostics and sensitivity analyses.
6. Ask whether the Northern-temperate association is stronger than tropical and southern alternatives.
7. Use Campanula/Izu measurements only as complementary evidence for intermediate links that the global dataset cannot itself identify.

## Required outputs

- a v1-to-v2 claim table listing each original claim, its revised measurement, its competing explanations, and its allowed interpretation;
- data-process diagnostics per island;
- trait evidence coverage maps and source-type summaries;
- domain-specific composition models before any composite syndrome model;
- a final interpretation that distinguishes global comparative compatibility from stronger field-based causal evidence.
