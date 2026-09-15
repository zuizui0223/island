# Existing-data design for an NEE-targeted manuscript — island-first version

Status: development plan, not a preregistration of already inspected results.
Created 2026-09-15. Parent interpretation reference supplied for this design: `e4a73796a`.
The current Database 1.0, frozen empirical results, thresholds, failed gates and canonical v8 submission surface remain unchanged.

## Starting biological question

The manuscript must begin from the island problem rather than from a general theory imposed on islands:

> **Why are island flowers expected to become duller, more generalized and more self-reliant as isolation increases?**

The existing results force that classical question to change in a sequence:

1. Is there actually one floral island syndrome globally?
2. Do floral accessibility and reproductive assurance behave as one coupled syndrome?
3. If a strong island syndrome appears, is it a robust response beyond genus composition, or is it represented mainly through which lineages assemble?
4. Could the apparent island pattern be produced by trait/species observation processes?
5. Which named mechanism, if any, survives independent evidence rather than floral-template inference?
6. Only after those island questions are answered: does the island case reveal a more general property of macroecological trait syndromes?

The intended conceptual ascent is therefore bottom-up:

`island-flower puzzle -> regional/component contradiction -> taxonomic assembly depth -> falsification/claim boundary -> general trait-syndrome inference`

The general concept **assembly depth** is an inference emerging from the island problem, not the premise used to select the island system.

## Working thesis and claim ceiling

Working thesis:

> An apparent floral island syndrome can be represented by context-dependent lineage composition rather than one universally coupled organismal response; the strongest Palearctic branch is currently localized near the family-to-genus transition.

Composition adjustment localizes where the association is represented. It does **not** identify a causal assembly process, prove dispersal alone, exclude within-lineage evolution, or identify pollinator loss.

Nature Ecology & Evolution is an aspirational target, not a completion criterion. The paper is strengthened by narrowing unsupported claims, not by adding more nominally significant tests.

## Fixed evidence and rival explanations

Retain Database 1.0, the fixed 8,265-island universe, primary endpoints, source definitions, all-analysis/direct-only separation and historical result locks.

Do not reopen or rescue:

- N1 remains unsupported and N2 remains unopened after its gate;
- GloBI source breadth remains not promoted;
- H5c remains bounded/non-supporting;
- observed global geometry remains unpromoted;
- H5d remains non-identifiable;
- new simulation results cannot retroactively qualify those failed gates;
- thresholds must not be moved after outcome inspection.

Rival explanations to distinguish where the existing data permit:

1. one universally coupled floral-island response;
2. context-dependent assemblage composition / lineage sorting;
3. observation selection in island floras and trait resolution;
4. a robust beyond-genus response compatible with within-lineage change;
5. pollinator-mediated filtering, which requires evidence independent of floral templates.

## Ordered work packages and stop rules

### P0 — immutable claim-to-result reconciliation

Island question:

> **What exactly do the existing data already establish about the floral island syndrome?**

Build one inventory linking every headline statement to the exact workflow run, artifact, digest, output table, estimand, denominator, evidence scope, floristic stratum, source mode and uncertainty quantity.

Rules:

- verify artifact ID and digest before any recalculation;
- distinguish empirical output from descriptive re-expression and from interpretation;
- a missing/expired artifact is unavailable evidence, not a negative result;
- manuscript/result discrepancies are repaired before new modeling;
- no new result is allowed to alter the historical classification of an older gate.

Deliverables:

- `docs/chapter1_p0_claim_result_ledger_20260915.md`;
- reproduced figure-source table inventory;
- explicit list of any unresolved provenance mismatch.

Stop: P1 does not open if a central H2/H3 claim cannot be traced to an immutable source.

### P1 — is the observed floral island syndrome genuinely localized at taxonomic assembly depth?

Island question:

> **When the strong Palearctic island syndrome is observed, is it primarily expressed through taxonomic assembly, or does a robust beyond-genus response remain?**

This is the decisive NEE-facing gate.

Existing H3 already uses family and genus as grouping variables rather than trait imputers, and family/genus residuals are computed on the same observed island species after a common eligibility filter. The new work therefore targets three remaining alternative explanations: sample loss, uncertainty in attenuation, and generic model/grouping flexibility.

#### P1a — exact paired-support audit

For each `evidence_scope x stratum x source_mode x response_axis` cell, verify and hash the exact island-species membership used to form:

- `observed_score`;
- `after_family_residual`;
- `after_genus_residual`.

Required identity:

- same observed island species at all three stages;
- same island rows and model covariates;
- same information weights;
- family and genus both source-eligible before any stage-specific score is created.

Any cell violating exact support identity is excluded from the new P1 defense and reported as non-comparable; it is not repaired by changing thresholds.

Direct-only is the primary safeguard. All-analysis is complementary and cannot be counted as independent replication.

#### P1b — paired spatial-block uncertainty for attenuation

Estimate uncertainty in the attenuation itself rather than comparing separate significance labels.

Primary two-axis vector:

- `generalized_accessible`;
- `selfing_core`.

For every comparable profile define:

- `A_family = 1 - ||beta_family|| / ||beta_observed||`;
- `A_genus = 1 - ||beta_genus|| / ||beta_observed||`;
- `A_family_to_genus = (||beta_family|| - ||beta_genus||) / ||beta_observed||`;
- conditional genus attenuation `1 - ||beta_genus|| / ||beta_family||`.

Bootstrap contract:

- resampling unit: frozen `spatial_block`;
- all stages are refit/recomputed within the same bootstrap draw;
- draws: 2,000;
- random seed: `20260915`;
- interval: paired percentile 95% interval;
- no independent bootstraps for observed/family/genus stages.

Primary interpretive safeguard: the family-to-genus extra attenuation must retain a positive paired interval in the direct-only profiles used for the headline interpretation. If it does not, narrow the assembly-depth claim.

#### P1c — matched-complexity pseudo-genus null

Question:

> **Would an arbitrary fine grouping with genus-like complexity absorb the response just as strongly as biological genus membership?**

Primary null generation is outcome-blind with respect to island distance and island response values.

For each family:

- preserve the number of genera;
- preserve the empirical genus group-size multiset;
- repartition species among pseudo-genera **within family**;
- recompute source-group positions, source availability, island expectations and the full observed/family/pseudo-genus decomposition exactly as for real genera;
- retain the existing minimum source-scored-species and island-support rules;
- never use the observed island outcome to choose a partition.

Null repetitions: `999`.
Random seed: `20260915`.

Primary inferential statistic:

- direct-only median conditional family-to-genus attenuation across the frozen 4 source modes x 2 floristic strata;
- one-sided empirical tail probability with `(1 + exceedances)/(1 + 999)`.

Consistency report (not a second promotion test):

- cell-level real-vs-null percentile for all eight direct-only profiles;
- number of source-mode/stratum profiles in which real genus attenuation exceeds the null median;
- all-analysis values shown only as complementary sensitivity.

If real genus attenuation is not distinguishable from the matched-complexity null, replace “assembly depth at genus” with the narrower statement that the response is sensitive to fine taxonomic composition.

#### P1d — model-flexibility and spatial safeguards

Because the response model is the same across observed/family/genus stages, report explicitly that model formula, covariates and clustering are unchanged and only the response decomposition differs.

A leave-one-spatial-block-out summary may be added as a post-baseline robustness extension only if its rule is frozen before inspection. It is not required to open P2 if P1a-c are decisive.

### P2 — does the classical floral island syndrome behave as one coupled syndrome?

Island question:

> **Do floral accessibility and reproductive assurance actually move together as the classical floral island syndrome implies?**

Retain the disjoint primary accessibility and reproductive-core definitions. Use direct between-context vector contrasts and joint uncertainty; never infer context differences from separate significance tests.

Audit:

- common-island support between contexts;
- shared denominators and evidence scopes;
- covariance of the two primary components;
- paired effect visualization in the same estimand space.

Preserve regions where components decouple. Colour and named pollinator templates remain secondary; overlapping floral traits are not counted as independent syndrome confirmation. Species-level selfing evolution and temporal asynchrony remain outside this test.

### P3 — could observation processes manufacture the apparent island syndrome?

Island question:

> **Could the isolation-associated island pattern be an artefact of which species and traits are recorded?**

Reproduce V5 and V6 separately before adding any joint sensitivity extension.

A joint extension, if opened, must:

- specify flora-list missingness and trait-state selection simultaneously;
- preserve observed information weights;
- display robust and fragile parts of parameter space;
- treat the grid as assumptions, not a posterior over true completeness;
- avoid unconstrained occupancy/detection fitting to uncalibrated occurrence records.

If conclusions become prior-dominated or non-identified, report partial-identification bounds and stop.

### P4 — can a smooth global island gradient conceal local transitions? (optional)

Island question:

> **Could different local island systems undergo thresholds or offsets even though the global assemblage response is smooth?**

Do not open observed cross-component breakpoint fits first.

Feasibility must compare, under shared spatial dependence and unequal missingness:

- synchronous steps;
- offset steps;
- heterogeneous smooth curves;
- mixed shapes;
- nulls.

Use separate calibration and held-out simulation seeds. Predefine the location scale, minimum biologically relevant offset, multiplicity family, false-offset ceiling, interval coverage and recovery targets before simulation results are inspected.

If recovery fails, stop at feasibility and hand the biological question to `izu-core`. A future feasible design cannot reverse the existing geometry/H5d classifications. Distance offsets are not temporal lags or measured ocean-crossing limits.

P4 must never delay a defensible manuscript.

## Four-figure narrative — island problem first

1. **Why should island flowers change?** Rival generators, island sampling footprint and observation limits.
2. **Is there one floral island syndrome?** Joint accessibility/reproductive responses with direct between-context uncertainty.
3. **Where does the strongest island syndrome live in the lineage hierarchy?** Assembly depth, direct-only exact-pair safeguards, paired attenuation uncertainty and matched-complexity null.
4. **What can and cannot explain it?** Observation tipping surfaces plus identifiable versus non-identifiable geometry/mechanism routes.

The general macroecological claim belongs at the end of Discussion, not at the beginning of Introduction:

> The island case reveals that a cross-sectional trait syndrome can have an assembly depth; before assigning repeated adaptation or mechanism, one must ask at what lineage level the syndrome is represented.

## Manuscript order

### Introduction

1. classical island-flower expectation: generalized, less conspicuous, reproductively assured;
2. empirical tension: regions and components need not agree;
3. island-specific inferential problem: within-lineage response versus differential lineage assembly can generate the same assemblage pattern;
4. questions: one syndrome? which components? where in the taxonomic hierarchy? which alternatives survive?
5. only one closing sentence hints at the broader trait-syndrome implication.

### Results

1. one global floral island syndrome is not recovered;
2. the classical syndrome decouples by biogeographic context/component;
3. the strongest Palearctic island syndrome localizes at family-to-genus depth;
4. observation and grouping-complexity defenses;
5. independent global pollination/geometry routes remain unpromoted or unidentified.

### Discussion

1. what became of the floral island syndrome?;
2. island syndromes can be assembled;
3. why the family-to-genus transition is or is not biologically specific after P1;
4. what is still unknown: within-lineage evolution, pollinator mechanism, local thresholds;
5. Chapter 2 resolves local mechanism/geometry;
6. **only then**: the island case exposes the wider macroecological problem of assembled trait syndromes.

## Completion criteria

- every headline claim traces to a verified immutable result;
- P1 uses exact paired support and reports uncertainty in attenuation itself;
- P1 explicitly tests whether arbitrary genus-like grouping can reproduce the attenuation;
- the central inference is narrowed if any defense fails;
- every new analysis is labelled post-baseline and has a stop decision;
- software tests and empirical-output verification are separate;
- main text, supplement and figures use the same estimands and claim ceiling;
- public release respects source rights and distinguishes full analytical inputs from the redistributable Database 1.0 subset.

First executable gate: **P0, then P1**. P4 remains optional.
