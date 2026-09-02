# Chapter 1 progressive analysis contract

## Decision

Chapter 1 is **not frozen because the current trait matrix is incomplete**.  PR #141 is retained as a reproducible baseline snapshot, not as the final data freeze.

What is frozen now:

- the scientific hypotheses;
- the island universe and geographic/status layers;
- trait ontology and evidence precedence;
- model order and covariates;
- support thresholds and multiplicity rules;
- observation/source/lineage/area robustness gates;
- the causal claim ceiling.

What is allowed to change:

- species-level trait evidence;
- response-specific island support created by that evidence;
- estimates, uncertainty, significance and biological conclusions.

Every new trait wave is converted to the same `chapter1_trait_snapshot_v1` schema and the same analysis ladder is re-estimated from the beginning.

## Central hypothesis

### Source-conditioned, channel-gated island assembly

> Island isolation does not impose one universal floral island syndrome. Isolation changes which source-available floral and reproductive strategies are represented, and a classic simplification/reproductive-assurance syndrome should be strongest only where a regionally available pollination channel is disrupted and not effectively replaced.

The hypothesis has two deliberately separated levels.

### Level A — directly testable with the global plant data

```text
geography / source accessibility
    -> biogeographically contingent trait-response vectors
    -> source-available lineage entry / loading
    -> area-conditioned assemblage filtering
```

This level can be tested every time trait coverage improves.

### Level B — mechanistic extension

```text
regionally available pollination channel
    -> retained / disrupted / structurally absent
    -> residual floral/reproductive filtering beyond source + lineage + observation controls
```

This level is **not** established by floral phenotype.  It is promoted only after independent, outcome-blind channel evidence is available.  Bombus remains the first natural-experiment candidate, not a globally imposed cause.

## Predeclared hypothesis ladder

### H1 — universal-syndrome rival

A common response vector appears across contexts.

This is weakened if direct between-context vector tests show heterogeneity.  A significant result in one region and a nonsignificant result in another is never used as evidence of difference.

### H2 — biogeographic branching

Isolation/source accessibility produces different multivariate floral/reproductive responses in different prespecified contexts.

Required evidence:

1. within-context multivariate vector test;
2. direct between-context vector-difference test;
3. persistence in native non-endemics where support permits;
4. observation and source-pool sensitivities.

### H3 — source/lineage assembly

A supported assemblage response is decomposed into:

- which source-available genera enter/are represented;
- additional species loading within represented genera;
- residual trait structure not represented by genus composition.

Lineage composition is part of the biological assembly process, not a nuisance to be discarded after the fact.

### H4 — area/capacity moderation

Distance effects may be stronger on smaller islands if target size, habitat capacity or persistence opportunities strengthen assembly constraints.  Area stays continuous; no small/large cut point is introduced after seeing results.

### H5 — channel-gated residual mechanism

Only after H2-H4 are accounted for, test whether an independently measured pollination-channel deficit predicts residual trait filtering.

Promotion requires all of:

- source-region channel exposure defined without floral outcomes;
- effort-aware occurrence/abundance or direct functional interaction evidence;
- retained vs disrupted vs structurally absent contrast;
- residual explanatory value beyond source, lineage, area, climate, observation and spatial controls.

Climate suitability alone cannot define realized loss, and structural absence is not a deficit.

## Fixed secondary pollination-syndrome concordance

After the pollinator-name-free primary H2 pattern is estimated, every trait snapshot also reruns one prespecified **secondary** architecture-concordance family:

- `large_bee_like`
- `butterfly_like`
- `bird_like`

The templates are fixed in `config/chapter1_pr138_pollination_syndromes.yml`; the test family is fixed as `sampled_guild_concordance` in `config/chapter1_global_branching.yml`. It is evaluated in both broad analysis regimes and formal biogeographic realms, with BH-FDR within the predeclared `context_layer x axis_set x stratum x support_tier` family.

This secondary analysis asks whether the independently estimated plant assemblage response is concordant with sampled floral architectures described for candidate pollination-service channels. It is deliberately non-exhaustive: moths, flies, bats, beetles, wasps and other channels are unmodelled, not absent.

The same family is regenerated on direct-only evidence and is carried through the source-adjusted/joint-source-selection sensitivity. It may strengthen, weaken, reverse or become newly testable as trait coverage improves; the templates and multiplicity family do not change in response to those results.

A positive or negative guild-labelled score is **not** a visitor assignment and cannot establish pollinator abundance, visitation, effective pollen transfer, historical loss, mobility or functional replacement. Any H5 mechanism promotion still requires independent channel evidence.

The PR #141 baseline at HEAD `660bc02bedca8e15548bb913bc69b6b904fb8d89` is documented in `docs/chapter1_progressive_pollination_concordance.md`.

## Trait coverage is not a global gate

The present matrix being below 60% for some domains is a reason **not to freeze the final result**, but `60%` is not itself an inferential threshold.

A global fill percentage can be misleading because the actual model denominator is response- and region-specific. Therefore:

- `<30` supported islands for a response: not promoted;
- `30-49`: pilot;
- `>=50`: count component of confirmatory support;
- pairwise regional comparisons require the same response to meet the same support tier in both regions.

A new trait wave may therefore do any of four useful things even if global coverage remains below 60%:

1. make a previously untestable response enter the pilot tier;
2. promote a pilot response to confirmatory support;
3. narrow uncertainty for an already testable response;
4. change a previous conclusion, which is scientifically allowed and must be recorded rather than treated as a CI failure.

## Versioned trait snapshot

Every acquisition wave must be materialized as one exact 106,295-species x 3-axis table with:

- `accepted_species`
- `axis`
- `trait_composition`
- `trait_names`
- `source_groups`
- `source_lineages`
- `quality`

A non-empty source `quality` label is not by itself an analysis-ready trait. A cell is
materializable only when `trait_composition` and its source provenance are present. The
snapshot manifest reports both source-reported fill and analysis-usable fill; reported
Low cells lacking a composition remain unresolved for analysis and are never converted to
trait zeros or silently inferred.

The materializer produces:

- `chapter1_trait_snapshot_species_axis.csv.gz`
- `chapter1_trait_ledger_all_analysis.csv.gz`
- `chapter1_trait_ledger_direct_only.csv.gz`
- `chapter1_trait_snapshot_coverage.csv`
- `chapter1_trait_snapshot_manifest.json`
- optional `chapter1_trait_snapshot_transition.csv.gz`

The transition audit explicitly labels new resolution, quality upgrades, evidence loss, quality downgrade and value revision. Evidence loss/revision fails closed unless explicitly acknowledged for a correction run.

## Progressive analysis order

The reusable workflow is `.github/workflows/run-chapter1-progressive-trait-analysis.yml`.

For each new trait snapshot it executes, in this order:

1. **Snapshot/QC** — exact denominator, ontology, provenance and direct-vs-Low separation.
2. **Atomic WHEN-WHERE** — status-stratified trait-state response vectors, all-analysis and direct-only.
3. **Plant response synthesis** — pollinator-name-free accessibility/generalization and reproductive-assurance scores.
4. **Secondary sampled-guild concordance** — fixed large-bee, butterfly and bird architecture templates under the separate secondary family.
5. **Biogeographic branching** — broad region and formal realm comparisons for the primary and secondary axis sets.
6. **Observation-selection sensitivity** — same response definitions under the frozen observation contract.
7. **Source-adjusted response** — fixed GIFT/source-assignment rules, recalculated with the new trait scores.
8. **Joint source + observation sensitivity** — including the fixed secondary concordance axes where source expectations are evaluable.
9. **Lineage entry/loading** — source-matched genus representation and within-genus species loading using the new functional positions.
10. **Continuous distance x area moderation**.
11. **Snapshot manifest** — coverage and result-table shapes are recorded without asserting that biological coefficients must equal the previous wave.

The workflow deliberately does **not** contain numerical assertions such as “Palearctic must remain significant”.  Such assertions belong only to an immutable historical snapshot validator.  In the progressive workflow, a result changing after better evidence is a result, not a software failure.

## How later trait waves upgrade the analysis

A new acquisition workflow only needs to produce a valid species x axis coverage file. Dispatch the progressive workflow with:

- source workflow run ID;
- source artifact name;
- relative path to the species x axis file;
- optionally the previous snapshot artifact for transition auditing.

No model code, hypothesis text, region definition or p-value threshold is changed when a new wave is added.

Recommended update labels for the scientific comparison between snapshots:

- `stable`
- `newly_evaluable`
- `no_longer_supported`
- `direction_changed`
- `precision_changed_only`
- `unresolved`

These labels describe evidence evolution; they are not reasons to discard a later snapshot.

## Role of PR #141

PR #141 remains valuable because it provides a fully reproducible baseline under the current evidence.  It is no longer treated as “the final Chapter 1 trait freeze”.

Its role is:

> **baseline snapshot against which later, higher-coverage snapshots are prospectively compared under an unchanged analysis contract.**

Final manuscript freezing should occur only after the trait campaign reaches a defensible stopping point and the major response families show acceptable support/stability.  The stopping decision must be based on predeclared coverage/recoverability and scientific information gain, not on whether the newest p-values are convenient.

## Claim ceiling

Even at high trait coverage, Chapter 1 alone does not identify:

- historical pollinator loss;
- effective pollination service;
- bird/butterfly/Bombus replacement from phenotype;
- historical source ancestry;
- in-situ evolution versus assembly;
- Baker's law from self-compatibility alone.

Those mechanisms require independent evidence or Chapter 2/field-system tests.
