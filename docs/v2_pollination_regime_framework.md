# v2: conditional island floral syndrome  -  pollination-regime framework

## Status and purpose

This document refines the v2 scientific scope before new outcome analyses are
run. It does **not** posit a universal global `Bombus syndrome`, nor does it
treat absence of a Bombus record as evidence of pollinator loss.

The central question is instead:

> Across island floras, are shifts in pollen-vector mode, reproductive
> assurance, and floral phenotype conditionally associated with loss of a
> biologically applicable pollination-function channel, after source-pool,
> geography, establishment, and observation processes are represented?

The Bombus channel is one predeclared, confirmatory case within that broader
framework. It is expected to be most interpretable in Northern-temperate
systems where Bombus is a biologically meaningful component of the regional
source-pool context. Other regions are used to identify applicability,
functional replacement, contingency, and data limits; they are not forced into
the same Bombus mechanism.

This is a global comparative synthesis. It should not be described as a formal
study-level meta-analysis unless a separate database of published effect sizes
and sampling variances is constructed.

## Claims allowed by this design

The design can evaluate whether comparative patterns are compatible with a
predeclared pollination-function hypothesis. It cannot by itself prove that
selection through a particular pollinator caused an evolved floral state.

Allowed interpretation:

```text
conditional island-level association
+ robustness to source-pool, geography, and observation-process alternatives
+ evidence-quality diagnostics
= comparative support compatible with a channel hypothesis
```

Not allowed without independent field, experimental, or lineage evidence:

```text
Bombus loss directly caused floral evolution on an island.
```

## Why the syndrome is conditional

Island isolation can alter several distinct processes:

```text
isolation / geography / source pool
    -> arrivals, establishment, persistence, and recorded flora
    -> pollination-function regime
    -> pollen-vector mode and reproductive assurance
    -> floral phenotype composition
```

A lack of one pollinator group does not have one inevitable floral outcome.
Depending on the regional species pool and the island's functional replacement
options, an island may show:

1. lower reliable animal-pollination service with greater representation of
   reproductive-assurance traits;
2. continued animal pollination through alternative bees, birds, flies,
   butterflies, hawkmoths, or other generalist visitors;
3. greater representation of wind-pollinated or mixed-vector taxa; or
4. no detectable floral-composition shift once source-pool and establishment
   filters are accounted for.

Consequently, white or less visually conspicuous flowers are not treated as a
universal island endpoint. They are one possible phenotype among several that
may be associated with particular pollination-function regimes.

## Analysis domains

### Domain A: whole-flora pollen-vector and reproductive-assurance composition

This domain is not restricted to animal-pollinated species. It asks whether
island floras differ in the representation of:

```text
pollen-vector mode:
  biotic
  abiotic_wind
  mixed
  unresolved

reproductive-assurance state:
  self_incompatible
  self_compatible
  autonomous_selfing_capacity
  delayed_selfing
  cleistogamy
  unresolved
```

This domain addresses an upstream question omitted by v1: whether isolated
islands differ not only within animal-pollinated taxa, but also in the relative
representation of animal-dependent, wind-pollinated, mixed-vector, and
reproductively self-reliant taxa.

Wind pollination, self-compatibility, autonomous selfing, and pollinator
independence are not interchangeable traits. They are recorded and modelled as
separate domains.

### Domain B: floral phenotype within animal-pollinated taxa

This domain retains the v1 continuity question but restricts its confirmatory
Bombus analysis to islands where a Bombus channel is biologically applicable.
Primary outcomes remain multistate or continuous:

```text
flower colour
floral symmetry
tube-depth class
floral openness
flower size
inflorescence display
reward accessibility where directly evidenced
```

Predeclared binary or ordinal contrasts are secondary summaries for continuity
with v1 and sensitivity analyses. They must be described as
`Bombus-associated` versus `lower-Bombus-association` contrasts, not as direct
measures of evolutionary reduction. In particular, white flowers are not
presumed to indicate low animal-pollinator dependence.

## Island strata

### Core confirmatory set

Islands with `applicability = applicable` under
[`config/bombus_applicability.yml`](../config/bombus_applicability.yml), adequate
flora coverage, and interpretable Bombus observation diagnostics.

This is the only set used for the primary confirmatory comparison of a
Bombus-channel model against geography/source-pool alternatives.

### Regional contingency set

Islands where a Bombus channel is applicable but regional evidence,
functional-service information, or coverage is insufficient for the primary
confirmatory analysis, or where the analysis is explicitly designed as a
regional comparison.

Results are exploratory and are used to assess heterogeneity and transportability,
not to declare a global syndrome.

### Out-of-domain set

Islands with `applicability = structurally_not_applicable`.

These islands are not coded as Bombus-deficit islands. They provide evidence
about alternative pollination regimes, possible functional replacement, and the
limits of the Bombus hypothesis.

### Unresolved set

Islands with unresolved source-region assignment, native Bombus status,
observation process, or data quality. They remain in diagnostic outputs but
are excluded from confirmatory Bombus-channel inference.

## Bombus channel: definition and measurement layers

Within applicable islands, the biological state is represented through three
separate components:

1. `environmental_compatibility`  -  overlap between island conditions and the
   environmental envelope of Bombus taxa that are candidates from the assigned
   source region;
2. `occurrence_evidence`  -  detected, adequate non-detection, insufficient
   effort, or unresolved evidence after target-group effort diagnostics;
3. `functional_service_proxy`  -  evidence about activity season, richness,
   abundance, visitation, or documented functional use where available.

`proxy_quality` is not a fourth biological component. It is a measurement layer
that reports whether the preceding components are direct observations,
calibrated proxies, or unresolved. It is used for gating and sensitivity
analysis, not interpreted as pollination biology.

The channel is described as a proxy for the availability of a
Bombus-mediated pollination-function channel relative to the relevant source
region. It is not assumed that all Bombus are long-tongued, buzz-pollinating,
or the principal pollinator of every focal plant species. It is also not
assumed that island non-occurrence proves a historical loss event.

### Applicability precedes climate compatibility

`applicability` is assigned from source-region native Bombus status and
biogeographic context only. Island climate, environmental compatibility,
occurrence records, floral traits, reproductive traits, and model outcomes are
prohibited inputs to this decision.

This separation prevents the environmental-compatibility predictor from
partly defining its own analysis domain.

### Effort-aware non-detection

The primary classification follows
[`config/bombus_observation_diagnostics.yml`](../config/bombus_observation_diagnostics.yml):

```text
detected:
  >= 1 quality-filtered Bombus record

adequate_non_detection:
  Bombus records = 0
  AND target-group records >= 50
  AND minimum spatial, temporal, and dataset dispersion is met

insufficient_effort:
  Bombus records = 0
  AND one or more effort criteria is not met
```

A simple ratio such as `1 - Bombus records / target-group records` is not used
as the primary evidence score because it conflates abundance, collector focus,
taxonomic practice, and data-source structure. Occupancy models are optional
sensitivity analyses only when genuinely repeated observation units can be
constructed.

## Predeclared models

### M0: geography, source-pool, establishment, and observation baseline

```text
source pool + lineage composition + island area + isolation + climate
+ biogeographic region + establishment/reachability diagnostics
+ observation-process diagnostics
    -> pollen-vector composition
    -> reproductive-assurance composition
    -> floral phenotype composition
```

M0 is required. A strong M0 is scientifically informative rather than a failed
analysis: it may show that an apparent pollinator signal is indistinguishable
from geography, source-pool structure, establishment, or data generation.

### M1: conditional Bombus-channel association

For the Core confirmatory set and animal-pollinated taxa only:

```text
M0 covariates
+ environmental compatibility
+ occurrence evidence
+ functional-service proxy
+ alternative-function evidence where available
    -> floral phenotype composition
```

The key question is whether Bombus-channel information improves explanatory
or out-of-sample predictive performance beyond a climate-flexible M0. Climate
must be represented flexibly enough in M0 that a Bombus environmental proxy is
not merely a nonlinear rewrite of climate.

### M2: reproductive-assurance-mediated route

```text
Bombus-channel state
    -> reproductive-assurance domain
    -> floral phenotype composition
```

M2 is evaluated only in the subset with adequate reproductive-assurance
coverage. Self-compatibility is not treated as autonomous selfing.

### M3: joint direct-plus-mediated route

```text
Bombus-channel state
    -> floral phenotype composition
Bombus-channel state
    -> reproductive-assurance domain
    -> floral phenotype composition
```

M3 is evaluated only when M2 is feasible. A missing or weak M2/M3 result under
insufficient reproductive-trait coverage means that the mediation route is not
identified by the current data; it does not mean the route is absent.

## Alternative functional replacement

Where data permit, a predeclared `alternative_functional_replacement` layer is
recorded separately from Bombus state. It can include evidence for other bees,
birds, flies, butterflies, hawkmoths, or broad generalist visitor capacity.

This layer is not used to invent a universal cross-region pollinator index.
It is used to test the conditional prediction that the association between a
Bombus-channel deficit and floral/reproductive outcomes may be weakened,
reversed, or absent where another relevant pollination-function channel is
available.

A principal interaction for relevant analyses is:

```text
Bombus-channel state x alternative-functional-replacement
```

## Phased implementation

### Phase 0: outcome-blind preregistration freeze

Before trait outcomes are examined, freeze:

- source-region assignment and applicability categories;
- target group, record filters, spatial rule, and effort thresholds;
- primary and sensitivity trait contrasts;
- Core / Regional contingency / Out-of-domain rules;
- M0 covariates and climate representation;
- lineage/source-pool sensitivity design; and
- continuation thresholds.

### Phase 1: coverage, attrition, and identifiability audit

The first main output is a diagnostic attrition table, not an outcome model:

```text
all islands
-> flora-eligible islands
-> applicability = applicable islands
-> Bombus observation-evaluable islands
-> floral-signal-eligible islands
-> floral-architecture-eligible islands
-> reproductive-assurance-eligible islands
-> M1-eligible islands
-> M2/M3-eligible islands
```

The audit reports direct, broad-web, genus-distribution, and family-distribution
trait coverage separately. It also reports the source and lineage composition
of every retained island set.

Before outcomes are inspected, use the Phase 1 covariate structure in
outcome-blind simulations to determine whether M0 and M1 can be separated, and
whether M2/M3 coefficients are identifiable at plausible effect sizes.

### Phase 2: M0 versus M1

Run the primary M0-versus-M1 comparison in the Core confirmatory set. Analyse
multistate and continuous floral domains before any composite `syndrome` score.
Report the incremental predictive contribution of Bombus-channel components,
the effect of flexible climate controls, and sensitivity to observation and
trait-evidence tracks.

### Phase 3: M2/M3 only in an eligible reproductive-assurance subset

Proceed only if the trait-coverage and outcome-blind identifiability criteria
are met. Otherwise report the reproductive-assurance route as unresolved by the
current comparative data.

## Lineage and source-pool safeguards

Source-pool lineage composition is a required M0 component, but it is not the
only phylogenetic safeguard. Sensitivity analyses must include at least one of:

- source-pool constrained null draws that preserve genus or family composition;
- expected trait composition conditional on source-pool lineage composition;
- clade-stratified reanalysis; or
- a species-level phylogenetic random-effect analysis where a suitable tree and
  coverage exist.

The purpose is to test whether an island-level association remains beyond the
lineage composition inherited from its source pool.

## Trait evidence and review policy

All trait analyses are reported as parallel evidence tracks:

1. `direct-conservative`;
2. `direct-broad-web`;
3. `expanded-genus-distribution`; and
4. `expanded-family-distribution`.

Human review is allocated to high-leverage records first: taxa frequent in the
Core set, taxa shared among many islands, taxa that can change an island-level
proportion materially, and records that influence both pollination-vector and
reproductive-assurance domains.

Sensitivity analyses exclude correlated source clusters, for example a single
flora, monograph, or web-series, rather than assuming all species-level errors
are independent.

## Continuation rule

The following values are provisional operational minima and must be checked
against Phase 1 outcome-blind simulations before the Phase 0 freeze:

```text
M1 pilot feasibility:
  >= 30 M1-eligible islands

M1 confirmatory target:
  >= 50 M1-eligible islands
  AND representation across multiple source-region groups

M2/M3 exploratory feasibility:
  >= 30 reproductive-assurance-eligible islands

M2/M3 confirmatory target:
  approximately >= 50 reproductive-assurance-eligible islands
```

Decision outcomes:

```text
A. M1 is identifiable and meets the confirmatory target:
   proceed with M0 versus M1 as the primary analysis;
   run M2/M3 only where reproductive coverage is adequate.

B. Floral outcomes are feasible but reproductive-assurance coverage is below
   the M2/M3 target:
   retain M1 as the primary analysis;
   report M2/M3 as exploratory or not identified.

C. M1 has fewer than 30 evaluable islands, or Bombus occurrence evidence is
   not evaluable for fewer than 30 applicable islands:
   do not claim a new global Bombus-channel test.
   Reposition v2 as a rigorous reanalysis of the v1-supported regional subset,
   with improved treatment of data generation, source pools, and uncertainty.
```

## Final positioning

The intended contribution is not a claim of a universal white-flower or
inconspicuous-flower island syndrome. It is a conditional test of whether a
biologically applicable Bombus-channel deficit is associated with pollen-vector
mode, reproductive-assurance state, and floral phenotype after alternative
explanations are made explicit.

The global framework is valuable even if it narrows the reliable inference to
a regional subset: identifying the boundary of a mechanism is a scientific
result, not a failure of scale.
