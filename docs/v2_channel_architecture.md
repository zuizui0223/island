# v2: conditional island floral syndrome — channel-audit architecture

## Status and boundary

`v1-freeze` preserves the exploratory draft and its first complete results. v2
is rebuilt from new data products and new scripts; it must not import v1
scripts, derived trait labels, or v1 model outputs as analysis inputs.

The detailed scientific and operational policy is
[`v2_pollination_regime_framework.md`](v2_pollination_regime_framework.md).
This architecture is its compact channel-audit companion.

v2 does **not** ask whether islands globally share one universal
Bombus-driven white-flower or inconspicuous-flower syndrome. It asks whether,
where a Bombus-mediated pollination channel is biologically applicable,
variation in that channel is conditionally associated with pollen-vector mode,
reproductive assurance, and floral phenotype after source-pool, geography,
establishment, and observation processes are represented.

```text
isolation / geography / source pool
    -> arrivals, establishment, persistence, and recorded flora
    -> pollination-function regime
    -> pollen-vector mode and reproductive assurance
    -> floral-trait composition
```

The Bombus route is one predeclared, regional mechanism within this broader
pollination-function framework. It is not assumed to apply where Bombus is
non-native or biologically irrelevant to the assigned source-region context.

v2 does not adopt the Campanula/Izu project as a second empirical question, and
it does not turn the global project into a full life-history experiment. It
borrows only its **channel-identification discipline**: specify what each
dataset can identify, distinguish direct measurements from proxies, compare
predeclared competing explanations, and avoid treating a terminal pattern as
proof of a particular mechanism.

## What stays from v1; what changes in v2

### Retained scientific core

1. Island geography is the starting context.
2. Bumblebee limitation or environmental mismatch remains a focal pollinator
   hypothesis where it is ecologically applicable.
3. Reproductive assurance remains a focal indirect mechanism, distinguished
   from self-compatibility alone.
4. Direct Bombus-channel associations and reproductive-assurance-mediated
   associations are compared explicitly where data allow.
5. Northern-temperate islands remain the expected strongest confirmatory test
   region; other regions identify functional replacement, contingency, and the
   limits of the Bombus hypothesis rather than serving as failed projections.

### Rebuilt implementation

1. Island floras are rebuilt using polygon-exact occurrence evidence,
   source-pool information, and observation-process diagnostics.
2. Whole-flora pollen-vector composition is modelled separately from floral
   phenotype within animal-pollinated taxa.
3. The trait database is rebuilt as a source-cited, LLM-assisted evidence
   system rather than a filled categorical table.
4. Floral signal, floral architecture, pollen-vector mode, and reproductive
   traits remain multistate or continuous wherever possible; v1 binary
   contrasts become secondary summaries.
5. Geography, establishment/reachability, observation bias, and lineage
   composition are required alternatives before biological interpretation.
6. The former PDI is one candidate environmental proxy, not direct evidence of
   pollinator service.

## Channel-audit lens borrowed from Campanula/Izu

For a species trait state `z` in an island regime, the Campanula project writes:

```text
W(z) = F(z) × E(z)
```

For the global analysis this is **not a new fitted model**. It is an
interpretation guardrail:

- `F(z)`: local reproductive contribution, including outcrossed and selfed
  viable seed production and paternal contribution where directly measured;
- `E(z)`: arrival, establishment, persistence, and reachability conditional on
  viable seed output;
- `W(z)`: observed island-flora membership or island trait composition.

The global dataset mainly observes `W`, through island-by-species occurrence
evidence. Therefore an association between a pollination-channel indicator and
a trait must not automatically be called direct pollinator-mediated selection
on `F`. It may be consistent with direct floral filtering, selfing-mediated
change, source-pool structure, establishment, or observation bias.

## Island strata and applicability

The confirmatory Bombus analysis uses a three-way applicability classification
plus an unresolved state, defined and frozen in
[`config/bombus_applicability.yml`](../config/bombus_applicability.yml).

- **Core confirmatory set**: `applicability = applicable`, with adequate flora
  coverage and interpretable Bombus observation diagnostics.
- **Regional contingency set**: applicable islands used only to investigate
  heterogeneity, functional replacement, or transportability.
- **Out-of-domain set**: `structurally_not_applicable`; these islands are not
  recorded as Bombus-deficit islands.
- **Unresolved set**: incomplete source-region or evidence information; retained
  in diagnostics but excluded from confirmatory Bombus inference.

Applicability is defined from source-region native Bombus status and
biogeographic context, not island climate. Environmental compatibility is a
subsequent predictor component; using it to define applicability would be
self-referential.

## Predeclared explanations

### M0: geography / source-pool / establishment / observation baseline

```text
source pool + lineage composition + island area + isolation + climate
+ biogeographic region + establishment/reachability diagnostics
+ observation-process diagnostics
    -> pollen-vector composition
    -> reproductive-assurance composition
    -> floral phenotype composition
```

M0 is the necessary baseline, not a side analysis. A strong M0 is informative:
it may show that an apparent pollinator pattern is indistinguishable from
geography, source-pool structure, establishment, or data generation.

### M1: conditional Bombus-functional route

For Core confirmatory islands and animal-pollinated taxa:

```text
M0 covariates
+ Bombus environmental compatibility
+ Bombus occurrence evidence
+ Bombus functional-service proxy
+ alternative functional replacement where available
    -> floral-signal and floral-architecture composition
```

M1 evaluates whether the Bombus-channel components add explanatory or
out-of-sample predictive value beyond a climate-flexible M0. It is interpreted
as trait sorting or as a pattern compatible with pollination-function filtering;
it is not causal proof without field or lineage evidence.

### M2: reproductive-assurance-mediated route

```text
Bombus-channel state
    -> reduced reliable outcross service
    -> reproductive-assurance traits
    -> floral-signal and floral-architecture composition
```

M2 is evaluated only in an eligible reproductive-assurance subset.
Self-compatibility alone is not autonomous selfing, and insufficient coverage
means that the route is not identified, not that it is absent.

### M3: joint direct-plus-mediated route

```text
Bombus-channel state
    -> floral traits
Bombus-channel state
    -> reproductive assurance
    -> floral traits
```

M3 is evaluated only when M2 is feasible. It is supported only when the direct
association remains after the measured reproductive-assurance domain is
represented, while proxy status and uncertainty are retained.

### Required alternative explanations

```text
establishment / reachability:
traits + dispersal + source pool + habitat + persistence
    -> island-flora composition

observation process:
survey effort + GBIF/data-source structure + spatial coverage + taxonomic coverage
    -> observed island-flora composition
```

These are included so that a Bombus-channel explanation is not an artefact of
colonisation history, extinction, accessible-coast sampling, specimen
concentration, or database structure.

## Trait domains retained without premature binarisation

### Pollen-vector mode

- biotic
- abiotic wind
- mixed
- unresolved

Wind pollination, self-compatibility, autonomous selfing, and pollinator
independence are not interchangeable states and remain separate domains.

### Floral signal

- primary colour: white, green/brown/inconspicuous, yellow/orange, red/pink,
  blue/purple, multicoloured/variable, unresolved
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
- pollen limitation, seed set, paternity, and inbreeding-depression evidence
  where available

Binary `plain/conspicuous`, `generalized/specialized`, and `SC/SI` contrasts
remain transparent secondary summaries for continuity with v1 and sensitivity
analyses; they are not sole outcomes and are not direct measures of evolutionary
reduction.

## Functional availability of Bombus

Within applicable islands, Bombus state is not a binary GBIF-absence variable.
It has three biological components:

1. `environmental_compatibility`: whether island conditions overlap the
   environmental envelope of candidate Bombus species from the assigned source
   region;
2. `occurrence_evidence`: quality-filtered, effort-aware evidence classified as
   detected, adequate non-detection, insufficient effort, or unresolved;
3. `functional_service_proxy`: evidence relevant to activity season, richness,
   abundance, visitation, or documented functional use when available.

`proxy_quality` is a measurement layer, not a fourth biological component. It
records whether the preceding values are direct observations, calibrated
proxies, or unresolved, and is used for gating and sensitivity analysis.

The former PDI is one candidate proxy for environmental compatibility, not
direct proof of pollinator service.

The non-detection rules are frozen in
[`config/bombus_observation_diagnostics.yml`](../config/bombus_observation_diagnostics.yml).
A missing Bombus record is never alone treated as evidence of Bombus absence.

## Analysis sequence

1. Freeze applicability, observation rules, trait contrasts, M0 covariates,
   lineage safeguards, and continuation criteria before outcomes are examined.
2. Build a polygon-exact, provenance-tracked island-by-species occurrence
   database, a flora observation-effort table, and a parallel Bombus
   observation-diagnostic table.
3. Build a source-cited, LLM-assisted trait evidence database.
4. Produce the coverage / attrition / identifiability audit before fitting
   outcome models.
5. Model whole-flora pollen-vector and reproductive-assurance domains
   separately from floral phenotype within animal-pollinated taxa.
6. Compare M0 and M1 in the Core confirmatory set before attempting M2/M3.
7. Run M2/M3 only in the reproductive-assurance-eligible subset.
8. Use regional contingency and out-of-domain sets to interpret functional
   replacement, regional heterogeneity, and the boundary of the mechanism.
9. Use Campanula/Izu measurements only as complementary evidence for
   intermediate links that the global dataset cannot itself identify.

## Required outputs

- a v1-to-v2 claim table listing each original claim, revised measurement,
  competing explanations, and allowed interpretation;
- island-level applicability, observation-process, and attrition diagnostics;
- trait evidence coverage maps and source-type summaries;
- whole-flora pollen-vector models before animal-pollinated floral-trait models;
- domain-specific composition models before any composite syndrome model;
- lineage/source-pool constrained sensitivity analyses; and
- a final interpretation that distinguishes comparative compatibility from
  stronger field-based causal evidence.
