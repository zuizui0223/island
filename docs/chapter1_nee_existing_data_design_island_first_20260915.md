# Existing-data design for an island-first NEE-targeted manuscript

Status: development plan, not a replacement for the current frozen submission surface or any preregistered result.
Created: 2026-09-15.

The supplied design referenced parent interpretation commit `e4a73796a`; that identifier is not resolvable in the current `zuizui0223/island` repository, so it is retained only as an external provenance note, not as verified Git ancestry.

## Framing principle

The manuscript must begin with the island-biological problem and only generalize after the island evidence has forced a broader inference problem.

Primary biological question:

> Why should island flowers become duller, more generalized or more reproductively self-reliant with increasing isolation, and is the apparent floral island syndrome actually one shared response?

The analysis then asks, in order:

1. Is there one global floral island syndrome?
2. Do floral accessibility and reproductive assurance remain coupled across biogeographic contexts?
3. If a strong island response exists, is it represented mainly by within-lineage change or by which lineages assemble on islands?
4. At what taxonomic depth does that response attenuate?
5. Can observation bias or simple mechanistic alternatives explain the pattern?
6. Only after those island-specific questions are answered: what does this imply for trait syndromes more generally?

The general concept `assembly depth` is therefore an outcome of the island analysis, not the starting premise imposed on the island system.

## Objective and central thesis

Build the strongest defensible manuscript from the fixed island data, not the largest number of significant tests. Acceptance at Nature Ecology & Evolution is an aspiration, not a completion criterion.

Working island thesis:

> The classical floral island syndrome is not one universally coupled adaptive trajectory. Its direction and component coupling vary among biogeographic contexts, and the strongest Palearctic response is represented largely through taxonomic composition at the family-to-genus transition.

Generalization, reserved for the late Discussion:

> The island case reveals that a macroecological trait syndrome can be expressed through hierarchical changes in community membership rather than repeated organismal change alone; identifying the taxonomic depth of a syndrome is therefore necessary before assigning mechanism.

Composition adjustment localizes association. It does not identify a causal assembly process, prove dispersal alone, or exclude within-lineage evolution.

## Fixed evidence and rival explanations

Retain Database 1.0, the fixed island universe, primary endpoint definitions, source definitions, direct/all-analysis separation and historical results. Keep N1 unsupported, N2 unopened after its gate, H5c bounded, observed geometry unpromoted and H5d non-identifiable. New simulations cannot retrospectively qualify those failed gates. Thresholds must not be changed to obtain support.

Island-specific rival explanations are:

- one universal coupled floral-island response;
- context-dependent lineage composition and sorting;
- observational selection in island records and trait resolution;
- within-lineage response that is not identified by endemic status alone;
- pollinator-mediated filtering, which requires evidence independent of floral phenotype.

## Ordered work packages and stop rules

### P0 — immutable claim-to-result reconciliation

Purpose: establish exactly what the present island evidence already proves before any new analysis is added.

Build one inventory linking each central island statement to the exact formal run, artifact, table, estimand, denominator, evidence scope and uncertainty interval. Verify hashes before recalculation. Reconcile any manuscript/result discrepancy before fitting a new model. A missing artifact is unavailable evidence, not a negative result.

Deliverables:

- claim ledger;
- reproduced figure-source tables;
- one canonical map from each island-biological statement to its frozen evidence.

Stop rule: no P1 extension is opened until the H2/H3 headline values reproduce exactly from immutable inputs.

### P1 — determine whether the observed floral island syndrome is genuinely expressed at genus-level assembly depth

This is the decisive work package.

Primary question:

> Is the strong Palearctic floral island response truly localized near the family-to-genus transition, or can the observed attenuation be produced by sample change, low-confidence evidence, model flexibility, spatial structure or arbitrary fine grouping?

Required safeguards:

1. **Exact paired support.** Compare observed, family-adjusted and genus-adjusted stages on the same species, islands, weights, source definition and response definitions wherever mathematically possible. Any support loss caused by changing the analysed sample must be separated from composition attenuation.
2. **Direct-only primary safeguard.** The direct-only result is the primary protection against genus-imputed Low creating artificial genus structure. Low-inclusive analyses are complementary, not independent replication.
3. **Paired uncertainty for attenuation.** Report uncertainty on the attenuation estimand itself, not only significance before and after adjustment. Preferred estimand is the paired change in response-vector magnitude, with spatial-block-aware resampling or an equivalent frozen-support uncertainty procedure.
4. **Matched-complexity null.** Test whether genus attenuation is larger than attenuation generated by non-taxonomic partitions with comparable group number and group-size distribution. This is a post-baseline robustness extension and its randomization contract, seeds and acceptance rule must be fixed before its result is inspected.
5. **Spatial robustness.** Use existing frozen spatial nulls first. Any new holdout or block-resampling extension must be labelled post-baseline.
6. **Model-flexibility audit.** Quantify whether adding genus structure mechanically absorbs comparable signal from outcomes where no taxonomic-depth localization is expected.

Optional only if available without derailing P1: compare genus with an equivalent phylogenetic-depth control. Do not make a new global phylogeny a prerequisite for a defensible manuscript.

Stop rule: if attenuation depends materially on sample loss, Low evidence or arbitrary matched-complexity partitions, narrow the claim from `assembly depth` to `composition sensitivity`. Do not switch to another favourable endpoint.

### P2 — test whether the classical floral island syndrome behaves as one coupled syndrome

Primary island question:

> Do accessibility/generalization and reproductive assurance move as one coherent island syndrome, or can they decouple across regions?

Retain the disjoint primary accessibility and reproductive-core definitions. Use direct between-context vector contrasts and joint uncertainty; do not infer regional differences from separate significance tests. Audit common-island support, shared denominators and covariance between components.

Primary robustness requirement: a common-support version should show that the directional contrast is not created merely by different data quality or island support between Palearctic and Tropical analyses.

Colour and named pollination-guild templates remain secondary. Shared floral traits must not be counted twice as independent confirmation of a syndrome.

Claim ceiling:

> broad biogeographic branching may be defended even if the tropical single accessibility axis remains observation-sensitive; the two regional branches should not be portrayed as equally robust.

### P3 — ask whether the apparent island syndrome could be generated by observation processes

Reproduce V5 and V6 separately before adding any joint extension.

The scientific question is not whether GBIF or compiled island floras are complete. It is:

> What observation processes would be required to erase the island-biological conclusions that survive P1 and P2?

Any joint extension must specify flora-list incompleteness and trait-state selection together, preserve original observed information weights, and display both robust and fragile regions. The sensitivity grid is an assumption set, not a posterior distribution over true completeness.

Do not estimate unconstrained occupancy and detection from opportunistic occurrence data. Bayesian fitting does not resolve non-identifiability by itself.

Stop rule: if the result becomes dominated by unconstrained assumptions, report partial-identification bounds instead of a fitted latent truth.

### P4 — ask whether a smooth global island gradient could conceal offset local transitions

Optional; not the manuscript spine and must not delay a defensible paper.

Do not open observed cross-component breakpoint fits. First test whether the existing design can distinguish a difference between two transition locations under:

- synchronous steps;
- offset steps;
- heterogeneous smooth clines;
- mixed response shapes;
- null models;
- shared spatial dependence;
- unequal missingness.

Predefine the location scale, minimum scientifically relevant offset, multiplicity family, false-offset ceiling, interval coverage and recovery criteria before simulation outcomes are inspected.

If recovery fails, stop at feasibility. Do not widen effect sizes or grids repeatedly until observed data become admissible. A successful future offset-feasibility design cannot reverse the frozen geometry or H5d conclusions. Distance offsets are not temporal lags or measured ocean-crossing thresholds.

## Island-first manuscript sequence

### Introduction

The Introduction should climb, not descend, in abstraction.

1. Start with the classic island question: why should island flowers become less conspicuous, more generalized or more reproductively assured as isolation increases?
2. Introduce the empirical complication: different regions and floral/reproductive components need not move together.
3. Only then introduce the rival generators of an assemblage syndrome: repeated within-lineage response, differential lineage assembly, or both.
4. End with the study question: where in the lineage hierarchy is the floral island syndrome expressed, and which simple explanations survive explicit tests?

Do not open with `trait syndromes have assembly depth`. That is the synthesis earned by the Results.

### Results

Preferred biological headings:

1. **There is no single global floral island syndrome.**
2. **Isolation reorganizes floral and reproductive components differently across biogeographic contexts.**
3. **The strongest Palearctic island syndrome is concentrated at the family-to-genus transition.**
4. **Observation bias and simple alternative explanations do not trivially account for the Palearctic assembly signal.**
5. **The upstream pollination mechanism and local threshold generator remain unresolved.**

### Discussion

The Discussion should remain island-centred for most of its length:

1. What remains of the classical floral island syndrome?
2. How can an island syndrome be assembled rather than repeated identically within lineages?
3. Why the family-to-genus transition is the key result, conditional on P1 surviving.
4. What is still unresolved: why those genera are sorted, the role of pollination, and within-lineage evolution.
5. Why global smooth gradients cannot resolve local threshold geometry, creating the Chapter 2 / `izu-core` handoff.
6. Only in the final synthesis: the island case reveals a broader problem for macroecological trait syndromes.

General examples such as urban, alpine, drought or fragmentation syndromes belong here, not in the opening motivation.

## Four-figure narrative

1. **Figure 1 — the island problem and rival generators.** Begin with the floral island syndrome; show repeated within-lineage response versus differential lineage assembly versus both, plus the observation boundary and Chapter 2 bridge.
2. **Figure 2 — does the island syndrome move as one vector?** Show joint accessibility/reproductive responses and uncertainty across contexts using common-support/direct safeguards.
3. **Figure 3 — where does the island syndrome live?** Make assembly depth the hero quantitative figure: paired observed → family → genus attenuation, direct-only protection, matched-complexity null and paired uncertainty if P1 passes.
4. **Figure 4 — cross-examination.** Observation tipping surfaces and the explicit boundaries of geometry and mechanism identification.

Selfing and pollination syndromes should appear as biological predictions and alternative mechanisms, not hard visitor labels. Distinguish arrival, establishment, persistence and effective service. Distance, area, climate and source pools are proxies; the current models do not separately measure each process.

## NEE-level decision rule

The paper should be promoted conceptually only if the island evidence earns it.

If P1 shows that the family-to-genus attenuation persists under exact paired support, direct-only evidence, paired uncertainty and matched-complexity nulls, the manuscript can make the stronger inference:

> The floral island syndrome has a measurable assembly depth.

If P1 does not survive those tests, the paper should retain the island result but narrow the conclusion to context-dependent floral/reproductive responses with strong composition sensitivity.

The journal target must not determine the inferential threshold.

## Completion criteria

- Every headline island claim is traceable to a verified immutable result.
- P1 determines whether `assembly depth` remains a measured result or is narrowed to composition sensitivity.
- Direct-only and observation-process safeguards are treated as defenses of the same primary island inference, not as independent discoveries.
- Any new analysis is explicitly labelled post-baseline and has a stop decision.
- Main text, supplement and figures use the same estimands and claim ceiling.
- Public release respects source rights and distinguishes the full analytical input from the redistributable database subset.

## First executable gate

**P0 first, then P1.**

P1 is now defined as the critical test of the island-biological interpretation, not as an abstract test of a pre-existing general theory. P4 remains optional and cannot delay a defensible manuscript.

Journal positioning reference retained from the supplied plan: Nature Ecology & Evolution aims page, accessed 2026-09-15. Broad significance must emerge from what the island evidence forces us to learn, not from selecting a general theory in advance.
