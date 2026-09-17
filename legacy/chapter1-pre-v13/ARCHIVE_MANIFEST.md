# Chapter 1 pre-v13 archive manifest

Date: 2026-09-18
Cleanup PR: #235
Pre-cleanup canonical commit: `ea60041ed31d1b159f2722074c6fcca26ebfeb4a`

## Purpose

Chapter 1 publication science is fixed to the v13 surface established by PRs #229–#234. This directory is the historical boundary for retired Chapter 1 storylines and implementations that are no longer part of the active package, CLI, workflow, or publication-facing documentation surface.

The cleanup is repository hygiene only. It does not change a v13 numerical result, trait snapshot, GloPL estimate, raw-colour result, or claim ceiling.

## Active files deliberately retained outside this archive

The active root keeps the canonical v13 manuscript/hypothesis/figure/result-lock surface, raw-colour and colour–architecture analyses, shared trait/database/geography infrastructure, and the frozen parent result locks that the v13 result lock validates directly:

- `config/chapter1_all_data_route_result_lock.json`
- `config/chapter1_v12_two_panel_result_lock.json`
- `config/chapter1_v12_h5_glopl_extension_result_lock.json`

These result locks remain active provenance dependencies; their retired implementation branches do not.

## Physically retained historical sentinels

For the main retired storylines, exact blobs are retained under their former relative paths below this archive:

- `src/island_v2/bombus_regime.py`
- `src/island_v2/chapter1_h5_bombus_upstream_bridge.py`
- `src/island_v2/chapter1_pr138_palearctic_restricted_ipw.py`
- `docs/chapter1_submission_freeze_20260915.md`
- `docs/manuscript_submission_contract.md`
- `config/chapter1_h5_bombus_upstream_bridge_v1.yml`
- `tests/test_chapter1_h5_bombus_upstream_bridge.py`
- `.github/workflows/run-chapter1-h5-bombus-upstream-bridge.yml`

## Retired from the active surface in PR #235

### Bombus / old pollinator-channel implementation

- `.github/workflows/prepare-bombus-niche-real-data.yml`
- `.github/workflows/run-bombus-occurrence-evidence.yml`
- `.github/workflows/run-chapter1-h5-bombus-upstream-bridge.yml`
- `config/chapter1_h5_bombus_upstream_bridge_result_lock.json`
- `config/chapter1_h5_bombus_upstream_bridge_v1.yml`
- `docs/chapter1_h5_bombus_upstream_bridge_result_20260916.md`
- `src/island_v2/bombus_absence_evidence.py`
- `src/island_v2/bombus_applicability.py`
- `src/island_v2/bombus_channel_score.py`
- `src/island_v2/bombus_diagnostics.py`
- `src/island_v2/bombus_island_analysis.py`
- `src/island_v2/bombus_niche_hypervolume.py`
- `src/island_v2/bombus_niche_inputs.py`
- `src/island_v2/bombus_niche_inputs_production.py`
- `src/island_v2/bombus_regime.py`
- `src/island_v2/bombus_spatial_thin.py`
- `src/island_v2/chapter1_h5_bombus_upstream_bridge.py`
- `src/island_v2/m1_m3_bombus_channel_input.py`
- `src/island_v2/trait_bombus_analysis.py`

### PR138 / Palearctic-centred implementation

- `src/island_v2/chapter1_pr138_biogeographic_pattern.py`
- `src/island_v2/chapter1_pr138_block_leverage_diagnostic.py`
- `src/island_v2/chapter1_pr138_direct_ledger_pattern.py`
- `src/island_v2/chapter1_pr138_equal_island_support.py`
- `src/island_v2/chapter1_pr138_evidence_ladder.py`
- `src/island_v2/chapter1_pr138_gift_source_pool.py`
- `src/island_v2/chapter1_pr138_information_weight.py`
- `src/island_v2/chapter1_pr138_lineage_representation_bridge.py`
- `src/island_v2/chapter1_pr138_outcrossing_restriction.py`
- `src/island_v2/chapter1_pr138_outcrossing_selection_stress.py`
- `src/island_v2/chapter1_pr138_palearctic_restricted_block_deletion.py`
- `src/island_v2/chapter1_pr138_palearctic_restricted_ipw.py`
- `src/island_v2/chapter1_pr138_pathway_decomposition.py`
- `src/island_v2/chapter1_pr138_realm_replication_queue.py`
- `src/island_v2/chapter1_pr138_realm_sensitivity.py`
- `src/island_v2/chapter1_pr138_regional_lineage_decomposition.py`
- `src/island_v2/chapter1_pr138_san_nicolas_status.py`
- `src/island_v2/chapter1_pr138_selfing_interaction.py`
- `src/island_v2/chapter1_pr138_si_genus_influence.py`
- `src/island_v2/chapter1_pr138_si_genus_source_null.py`
- `src/island_v2/chapter1_pr138_source_adjusted_branch_selection.py`
- `src/island_v2/chapter1_pr138_source_adjusted_pathway.py`
- `src/island_v2/chapter1_pr138_source_pool_sensitivity.py`
- `src/island_v2/chapter1_pr138_syndrome_analysis.py`
- `src/island_v2/chapter1_pr138_syndrome_block_deletion.py`
- `src/island_v2/chapter1_pr138_syndrome_distance_sensitivity.py`
- `src/island_v2/chapter1_pr138_syndrome_template_sensitivity.py`

### Retired pre-v13 publication/submission documents and scripts

- `docs/chapter1_p0_claim_ledger_20260915.md`
- `docs/chapter1_p1_final_decision_20260915.md`
- `docs/chapter1_submission_freeze_20260914.md`
- `docs/chapter1_submission_freeze_20260915.md`
- `docs/manuscript_submission_contract.md`
- `scripts/build_chapter1_v10_p2.py`
- `scripts/promote_chapter1_v11_p3.py`
- `scripts/promote_chapter1_v11_submission_surface.py`

### Retired tests

The Bombus, old H5, and PR138 tests corresponding to the modules above were removed from the active test suite. They remain recoverable at the pre-cleanup commit.

## Recovery rule

Historical files not duplicated physically in this archive are immutable and exactly recoverable from Git at:

`ea60041ed31d1b159f2722074c6fcca26ebfeb4a:<former/path>`

This avoids keeping a second importable copy of a large retired subsystem while preserving complete scientific and implementation provenance.

## Boundary going forward

Active Chapter 1 development must target v13 or later. A retired Bombus-, PR138-, Palearctic-centred, old H5, or pre-v13 submission branch must not be restored under active `src/`, `tests/`, `config/`, `docs/`, `scripts/`, or `.github/workflows/` without explicitly reopening the publication architecture.
