# Chapter 1 submission freeze — 2026-09-09

## Status

`submission_freeze_candidate`

This file fixes the current Chapter 1 paper, analysis provenance, claim ceiling, and merge boundary. It does not claim that all possible trait data have been exhausted, and it does not promote the prospective pollinator-side H5 mechanism into the Chapter 1 results.

## Canonical manuscript

The single canonical Chapter 1 manuscript is:

- `docs/chapter1_manuscript_full_v7_submission_order_20260909.md`

Earlier v2–v6 and intermediate manuscript drafts are superseded and intentionally excluded from the current submission diff. Their development history remains available in git history.

## Canonical scientific contract

- `config/chapter1_progressive_analysis.yml`
- contract: `chapter1_progressive_analysis_v1`
- canonical workflow: `.github/workflows/run-chapter1-progressive-trait-analysis.yml`

The fitted PR142 H1–H5 analysis is not redefined by later paper-level interpretation.

## Final trait snapshot

- source run: `34191508045`
- source artifact: `source-scale-batch-integration-34191508045`
- fixed species-axis denominator: `318,885`
- resolved species-axis cells: `222,688` (`69.83%`)
- reproductive-assurance resolved cells: `48,497 / 106,295` (`45.63%`)

## Final PR142 analysis

- workflow run: `34232450884`
- artifact: `chapter1-progressive-analysis-34232450884`
- artifact ID: `10058653212`
- digest: `sha256:b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`

## Frozen Chapter 1 interpretation

### H1

One universal floral/reproductive island syndrome is not recovered.

### H2

Source-separation responses branch among biogeographic contexts. Palearctic assemblages show increasing accessibility/generalization together with increasing reproductive assurance; tropical assemblages can show increasing reproductive assurance while specialized/attractive floral architecture is maintained or strengthened. Measured-climate-independent categorical realm causation is not established.

### H3

The broad Palearctic primary response is `observed 4/4 -> after family 4/4 -> after source-matched genus 0/4` across all-analysis/direct-only and all-native/native-nonendemic combinations. The strongest supported plant-side depth is therefore genus-level lineage assembly beyond family composition.

H3 identifies where the plant assemblage response is represented; it does not identify the ultimate cause of differential genus representation.

### H4

Area remains a measurement-sensitive modifier. All 16 frozen primary V3 classifications remain `retain_area_as_measurement_sensitive_modifier_only`; zero heteroskedastic-null promotion gates pass.

### H5

The fitted PR142 H5 remains not evaluable from the Chapter 1 plant database because no complete independent source-channel -> retention/disruption -> visitation -> single-visit effectiveness -> effective-service chain enters the primary analysis.

The paper may discuss two prospective causal routes without treating them as fitted results:

- H5a: pollination-channel change may act upstream by changing lineage establishment/persistence and therefore generate H3 genus assembly;
- H5b: the original frozen PR142 residual route tests pollinator-side effects that remain after source/genus composition is fixed.

## Pollination-associated floral architecture

The fixed `large_bee_like`, `butterfly_like`, and `bird_like` templates remain a secondary interpretation bridge only. Approximately 87% of their variance is shared by one source-trained plant-architecture factor. They do not identify pollinator identity, mobility, retention/loss, visitation, effectiveness, or replacement.

## Distance interpretation

The formal exposure remains `log1p_distance_to_continent_km`.

Distance is interpreted as a composite source-separation / connectivity / accessibility gradient, not as a direct causal force on floral phenotype and not as the exact historical source distance for every island. The conceptual `d -> 0` state is a high-accessibility boundary, not an additional fitted mainland datapoint.

## Progressive-wave stopping rule

Wave52 contained `184,917 / 318,885` analysis-usable cells. The final snapshot contains `222,688`, a gain of `37,771` usable cells including `+11,315` reproductive-assurance cells. The main H1–H5 claim structure remained stable after this increase.

Strict new-source acquisition also entered strong diminishing returns: one high-yield Orchidaceae source added +192 reproductive cells, later reviewed packets generally added 1–3, and the committed residual audit of 1,187 staging files / >2.25 million rows found zero immediately strict-ready reproductive records under the fixed evidence rules.

The freeze is therefore based on inference stability plus declining recoverability, not on an arbitrary 69.83% completeness threshold.

## Canonical supporting documents

Keep in the submission branch:

- `docs/chapter1_submission_hypothesis_framework_20260909.md`
- `docs/chapter1_figure1_hypothesis_tree_spec_20260909.md`
- `docs/chapter1_literature_positioning_20260909.md`
- `docs/chapter1_h3_h5_causal_hierarchy_20260909.md`
- `docs/chapter1_progressive_analysis_contract.md`
- V1–V5 checkpoint documents
- latest trait and progressive-analysis provenance/checkpoint documents

## Explicit exclusions from the Chapter 1 submission merge

The prospective H5 evidence-chain implementation is not part of the Chapter 1 submission freeze. Its commits remain in repository history for later mechanistic work, but the current merge candidate excludes its workflow, config, code, tests, evidence projection, and implementation contract.

Superseded manuscript drafts are likewise excluded from the current merge candidate; v7 is the only canonical manuscript.

## Chapter 2 handoff

Chapter 1 closes the global WHEN/WHERE and plant-side taxonomic-depth problem. Chapter 2 (`izu-core`) addresses HOW and proximal WHY through partner arrival/loss, realized community, functional matching, effective service, dependency/assurance, and response branching. Izu remains the prospective depth system for directly linked visitor -> SVD -> effective service -> reproductive-dependency measurements.

This handoff is a discussion-level consequence of the Chapter 1 identification ceiling, not retroactive validation of the Chapter 1 regional patterns.

## Merge decision rule

PR #142 is merge-ready only after:

1. this freeze boundary is reflected in the final changed-file set;
2. canonical Chapter 1 tests/workflows are green on the cleaned head;
3. no prospective H5 implementation or superseded manuscript draft remains in the submission diff;
4. PR #142 remains mergeable against `ch1-submission-freeze`.

When these conditions hold, merging PR #142 into `ch1-submission-freeze` is scientifically appropriate. A regular merge is preferred over squash so the progressive-wave provenance remains visible.
