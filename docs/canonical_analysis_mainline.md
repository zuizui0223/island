# Canonical Chapter 1 analysis mainline

## Outcome first

Chapter 1 directly tests **biogeographically contingent island floral
filtering**, not a Bombus-loss model. Its job is to determine where isolation is
associated with floristic, floral, or reproductive composition after floristic
status, genus composition, area, climate, spatial dependence, and measurement
support are represented.

The current evidence does not recover a coherent residual floral syndrome. It
recovers robust floristic/endemic turnover, a small Northern colour-category
residual, and broad null results after status and lineage control. The word
`syndrome` therefore remains a question, not the result.

The order is `pattern -> structural tests -> mechanistic concordance`. Pollinator
mechanisms do not enter until the pattern has survived the primary gates.

## Direct hypotheses

1. **H1 — Universal syndrome test.** Test whether a single pooled isolation
   gradient describes plain/generalized/SC or other floral-reproductive axes.
2. **H2 — Context contingency.** Test whether isolation associations differ by
   biogeographic region and floristic status after the observed genus
   composition is fixed.
3. **H3 — Category decomposition.** Test whether fine trait categories retain
   residual isolation associations after support, multiplicity, and coverage
   controls.

Bombus, birds, butterflies, hawkmoths, other bees, and flies are absent from
these direct hypotheses. They enter only as symmetric Discussion-level
candidate mechanisms if a multivariate trait vector survives the primary
analysis.

Direct pollinator deficit is not a missing Chapter 1 predictor. The primary data
are the collected floral and reproductive traits; pollination syndromes are the
Discussion vocabulary used to interpret, but not relabel, the resulting vector.

## Analysis sequence

### 0. Audit entry from the physical-island universe

The frozen denominator is 8,265 GSHHG candidates; the current direct-trait table
contains 4,370. Model entry into the observed analysis set from region, area,
distance, and climate before fitting trait outcomes. Missing islands are not
trait zeros. Report unweighted, prespecified inverse-probability-weighted, and
overlap-trimmed sensitivity estimates. See
[`chapter1_island_universe_attrition.md`](chapter1_island_universe_attrition.md).

### 1. Floristic status is a response

Every `island x species` record is assigned to `native_nonendemic`, `endemic`,
`introduced`, or `unresolved` from source-backed status evidence. Conflicts fail
closed to unresolved.

First estimate:

```text
endemic share among resolved native species
    ~ standardized log distance + standardized log area + climate PC1-4
```

Endemicity is not a nuisance covariate because it is itself a primary route by
which isolation may restructure an island flora.

### 2. Separate floristic strata

Each trait response is summarized for:

- all confirmed native species;
- native non-endemics;
- corroborated endemics where support permits; and
- introduced species as a secondary negative-control stratum.

A native-nonendemic association is compatible with colonization or
establishment filtering. An endemic-concentrated association is compatible with
endemic lineage turnover or in-situ differentiation. Neither proves the process.

### 3. Hold genus composition fixed

The primary trait response is:

```text
observed direct trait share - genus-fixed null mean
```

The null preserves each island's observed species membership and genus
composition while reassigning direct trait states only within genera. It asks
whether an apparent isolation association exceeds what lineage turnover already
predicts; it is not imputation.

### 4. Decompose supported categories

Broad binary outcomes can hide opposing categories. Fine colour, form,
architectural, and reproductive categories are tested only after status and
lineage control. Confirmatory category results require at least 50 islands,
spatial-block robust uncertainty, BH FDR within the declared family, and
direct-coverage sensitivity.

### 5. Require multivariate coherence before naming a syndrome

A single red/pink, size, or form association is not a pollination syndrome. The
term is reserved for a prespecified combination of multiple supported trait
axes with a coherent direction. Literature comparisons are trait-vector
concordance, not pollinator identification.

## Model and exposure

The baseline trait model is:

```text
genus-fixed trait residual
    ~ standardized log distance to continent
    + standardized log island area
    + climate PC1 + PC2 + PC3 + PC4
```

Models use direct species trials as weights and frozen spatial-block clustered
uncertainty. Observation support determines eligibility, weighting, and
sensitivity; it is not treated as a biological mechanism.

Distance to continent is a composite gradient of isolation, connectivity,
dispersal limitation, source-region turnover, and assembly history. An explicit
source-pool reconstruction is not required for the Chapter 1 total association,
but no coefficient may be described as a pure isolation or source-pool effect.
Region, floristic status, and the genus-fixed null represent the currently
observable structural context.

This is how the source-pool problem enters the mainline: distance can absorb part
of the mainland-connectivity and source-supply gradient in the total association,
while status and genus composition expose two important assembly pathways. It
cannot tell whether the remaining coefficient is caused by pollinator filtering.

## Trait scope and replaceability

The status/lineage workflow is currently verified for broad plain colour,
generalized form, and SC, and for fine colour and floral-form categories.

Flower size, tube depth, symmetry, autonomous selfing, and mating system remain
important candidate axes from the replaceable trait pipeline. Their current
4,370-island broad models are discovery/supporting analyses because those
responses have not yet passed through the same floristic-status and genus-fixed
workflow. A new trait source can replace an old ledger input, but an axis becomes
canonical only after the same status, lineage, spatial, multiplicity, and
coverage gates are rerun.

This distinction prevents the pre-status Northern flower-size gradient from
overriding the more diagnostic status/lineage evidence.

### Multivariate trait target

After pending axes enter the status/lineage workflow, report the residual slope
vector jointly rather than selecting a pollinator label from one coefficient.
Before inspecting that vector, freeze literature-based contrast directions for
Northern simplification/reproductive assurance and for bird-, butterfly-, and
hawkmoth/specialized-moth retention.

Run an omnibus multivariate region-by-distance test first, then report vector
concordance with those frozen contrasts. The labels belong in Discussion even if
the statistical contrast is supported. Direct pollinator data are not required
for this Chapter 1 endpoint.

## Current evidence hierarchy

| Level | Result | Canonical interpretation |
|---|---|---|
| Floristic status | Northern endemicity increases robustly with isolation | primary supported result |
| Broad residual traits | plain/generalized/SC slopes are approximately null after status and genus control | no broad residual syndrome recovered |
| Fine categories | small Northern all-native red/pink decline survives FDR and moderate coverage sensitivity | limited colour turnover, not a syndrome |
| Tropical category | non-endemic red/pink decline collapses at higher direct coverage | measurement-sensitive, not a main result |
| Floral form | no category with at least 50 islands survives FDR | no confirmatory form residual |
| Endemic form post-pilot | tropical generalized form: beta `-0.0432`, p `0.122`, 29 complete-case islands | unsupported at the frozen 30-island acquisition boundary |
| Bombus compatibility | no predicted incremental pathway after status and lineage | supporting falsification, not Chapter 1 mechanism |

## Bombus/hypervolume placement

`environmental_compatibility_max` is an ellipsoidal climatic-compatibility score.
It is not occurrence, absence, abundance, functional diversity, visitation,
service, deficit, or replacement. After status and genus composition were
represented, it added no supported predicted pathway for generalized form or SC;
the only tiny predictive improvement involved plain colour in the opposite
direction to the simple deficit prediction.

Therefore:

- it is not a canonical Chapter 1 predictor;
- it remains in provenance as an exploratory climatic-opportunity diagnostic;
- its unsupported result is retained rather than rescued or retuned; and
- no Bombus occurrence campaign is launched for Chapter 1.

## Discussion and chapter hand-off

If a future Northern multivariate residual vector survives, the formation
candidate is loss of a functionally important Bombus/large long-tongued channel
without effective replacement. If the same vector is weak or absent in Tropical
or Southern systems, the non-formation candidates are persistence, recolonization,
or replacement by birds, butterflies, hawkmoths, or other guilds.

Those guilds must not be collapsed into one story. Specialized moth dispersal
supports maintenance of long-tube/outcrossing systems; butterfly colour learning
supports chromatic signalling; bird persistence or recolonization can sustain
effective vertebrate pollination, and some bird-pollinated lineages converge in
avian-visible long-wavelength colour space. The predicted contrast is weakened,
absent, or heterogeneous Northern-type simplification—not a mandatory tropical
or Southern reversal. These comparisons do not identify a causal channel.

This is the intended Chapter 1 stopping point: a floral-trait pattern plus a
bounded syndrome-level Discussion. Island-level Bombus deficit or replacement
data would address causation later; their absence does not force Chapter 1 to
become a pollinator-distribution study.

Chapter 2 asks why patterns differ using direct pollinator identity, functional
diversity, effective service, replacement, and reproductive-assurance evidence.
Chapter 3 asks which floral components diverge within a lineage. Only the
dissertation synthesis can integrate those results to evaluate the broader
`Channel-gated island assembly` hypothesis.

## Verified implementation boundary

PR #133 at `f9d070440ce42a213507b26f42db43830e19e1b8` contains the full
status/lineage evidence workflow through:

- `flora_endemism_analysis.py`;
- `status_stratified_lineage_analysis.py`;
- `status_category_decomposition.py`;
- `status_category_fdr.py`; and
- `status_category_coverage_sensitivity.py`.

`northern_bombus_residual_analysis.py` is supporting rather than canonical.
The six focused test files for these modules pass 19 tests. Ruff also reports
pre-existing style findings on that PR head; they do not alter the verified
scientific outputs.

Because PR #133 now conflicts with current `main`, PR #135 ports the post-pilot
endemic form/SC subset onto current `main`. Its run `32721017664` passed and the
downloaded artifact `endemic-pilot-lineage-32721017664` verifies the tropical
endemic generalized-form null at the predeclared stop boundary. PR #135 does not
replace the broader PR #133 evidence; it updates the currently executable
post-pilot endpoint.

## Claim ceiling

Chapter 1 may claim context-dependent floristic and residual trait associations,
including null results. It may not claim historical Bombus loss, observed
pollinator deficit, functional replacement, effective pollination service,
causal floral evolution, single-trait pollinator identity, or a separation of
assembly from in-situ evolution.

## Remaining identification problems

The two remaining Chapter 1 targets are multivariate syndrome coherence and
selection caused by the geographically structured 8,265-to-4,370 observation
filter. Realized Northern functional deficit, effective alternative-guild
replacement, source-pool versus pollinator causal attribution, and assembly
versus in-situ evolution are claim boundaries for Discussion and later work,
not requirements for finishing the floral-pattern chapter.

The machine-readable contract is `config/canonical_analysis_mainline.yml`; the
current empirical verdict is fixed in
[`chapter1_when_where_existing_data_results.md`](chapter1_when_where_existing_data_results.md).
