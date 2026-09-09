# Chapter 1 analysis database — construction and provenance

This document is the database-side entry point for the Chapter 1 paper. It distinguishes the **analysis database actually consumed by the final pipeline** from historical acquisition experiments and broad inventories.

## Canonical dimensions

The Chapter 1 trait database has a fixed denominator:

- 106,295 analysis-applicable species;
- 3 strict axes;
- 318,885 species-axis cells.

Axes:

1. `flower_colour`
2. `floral_structural_complexity`
3. `reproductive_assurance`

The final paper snapshot contains **222,688 resolved cells (69.83%)**, including **48,497 reproductive-assurance cells (45.63%)**.

The final snapshot is identified by:

- workflow run `34191508045`;
- artifact `source-scale-batch-integration-34191508045`.

## Data layers

### A. Island universe and geography

The island universe is fixed at **8,265 islands**. Geography is built from the GSHHG-based island geometry/covariate infrastructure. The primary exposure used by the paper is `log1p_distance_to_continent_km`, interpreted as a composite source-separation / connectivity / accessibility gradient rather than a literal randomized isolation treatment.

The geography/regime artifact is consumed as a frozen upstream input by the canonical Chapter 1 workflow.

### B. Realised island floras

GBIF occurrence records provide the observed island-flora layer. They are treated as opportunistic observations, not a census and not absence data.

The canonical workflow tracks the frozen island-taxon table at:

- `data/v2/staging/gbif/collected/island_taxa.csv`

No-record islands are not converted to trait zeros.

### C. Floristic status and source-pool infrastructure

The paper uses frozen source/status infrastructure for:

- native versus introduced/unresolved status;
- `native_nonendemic` and endemic strata;
- source-pool expectations;
- family/genus composition and source-matched lineage representation.

These layers are used as assembly controls/decompositions, not as proof of exact historical ancestry.

### D. Trait evidence

Every accepted trait row must retain source identity/provenance. The canonical evidence order is:

1. species-direct High/Medium;
2. trait-specific Validated Low only if species-direct evidence is absent;
3. unresolved otherwise.

The paper contract prohibits:

- family-level inference as a fill rule;
- global fallback;
- Low evidence overwriting direct evidence;
- silent downgrade of existing High/Medium evidence.

The detailed promotion/reporting rules are in:

- `docs/chapter1_trait_coverage_contract.md`

The extraction prompt retained for evidence acquisition is:

- `prompts/trait_evidence_extraction_v2.md`

## Database build path

The current reproducible line is:

```text
raw / external source
        |
        v
source-specific extraction + audit
        |
        v
reviewed species-direct packet
        |
        v
data/v2/staging/traits/source_batches/*.json
        |
        v
cumulative source-scale integration
        |
        v
106,295 x 3 species-axis ledger
        |
        v
Chapter 1 trait snapshot
```

### 1. Recover a complete public checkpoint

PR #151 established a hash-verified, row-level recoverable checkpoint and restored missing historical Validated-Low values/provenance from the pinned sidecar without relabelling them as new direct evidence.

Relevant entry points after the PR #141 merge into `main`:

- `.github/workflows/recover-verified-trait-checkpoint.yml`
- `scripts/recover_verified_trait_checkpoint.py`
- `docs/trait_checkpoint_recovery.md`

### 2. Acquire by source, not by species micro-batch

The final acquisition strategy is source-scale: complete tables, databases, floras, monographs, or reviews are extracted and audited as one source before integration.

Each promoted packet is represented by a manifest under:

- `data/v2/staging/traits/source_batches/`

The manifests are the shortest path from a final database cell back to the reviewed source packet.

### 3. Integrate cumulatively

Completed source packets are replayed cumulatively from the recovered baseline. This prevents a later packet from accidentally dropping an earlier gain.

Primary workflow after PR #141 merges into `main`:

- `.github/workflows/integrate-source-scale-batch.yml`

The final paper snapshot came from Run `34191508045`.

## What is not the canonical database

The following must not be confused with the final Chapter 1 database:

- the separate broad/all-master acquisition inventory;
- candidate/page-hit counts;
- unreviewed machine candidates;
- old species/page micro-batch queues;
- historical nominal checkpoints that cannot be reproduced row-by-row;
- pollinator occurrence or floral-syndrome labels treated as trait truth.

The paper analysis uses the frozen species-axis snapshot identified above, not a percentage copied from a historical README.

## Reproducibility boundary

Database construction and biological inference are deliberately separated:

- database scripts decide whether a species-axis cell is supported and preserve provenance;
- `config/chapter1_progressive_analysis.yml` defines the scientific hypotheses and inference order;
- `.github/workflows/run-chapter1-progressive-trait-analysis.yml` consumes a frozen database snapshot without rewriting the evidence rules after seeing outcomes.

This separation is required for the progressive-wave design: evidence can improve, but the scientific contract does not change post hoc.
