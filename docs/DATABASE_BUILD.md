# Chapter 1 analysis database — construction and provenance

This document describes the database consumed by the **current corrected Chapter 1 submission** and separates it from historical acquisition/replay bundles.

## Current analysis dimensions

Current geographic analysis universe:
- **8,264 island units** after excluding one verified Eurasian continental split component.

Trait denominator:
- **106,295 analysis-applicable species**;
- **3 strict raw axes**;
- **318,885 species-axis cells**;
- **222,688 resolved cells = 69.83%**;
- **48,497 reproductive-assurance cells = 45.63%**.

Raw axes:
1. `flower_colour`
2. `floral_structural_complexity`
3. `reproductive_assurance`

Final trait snapshot:
- workflow run `34191508045`;
- artifact `source-scale-batch-integration-34191508045`.

The geography correction changed the exposure and removed one continental fragment; it did **not** rebuild trait evidence or change the species-axis denominator.

## A. Geography

The current paper uses source-matched GSHHG 2.3.7 island/continental geometry.

Primary exposure:
- minimum coastline separation on a mean-radius sphere;
- artificial dateline split edges omitted;
- continental sibling polygons recombined;
- one Eurasian split component excluded;
- 1,113 formerly spurious island zero distances repaired.

Current selector:
- `config/chapter1_submission_current.json`.

Corrected tables:
- `results/geography_20260924/`.

Methods and audit:
- `docs/chapter1_corrected_submission_20260924.md`;
- `scripts/geography_correction/README.md`.

Distance remains a composite source-separation/connectivity/accessibility gradient, not a randomized treatment.

## B. Realized island floras

GBIF occurrence records provide the observed island-flora layer. They are opportunistic observations, not a census and not absence data.

Frozen candidate flora source:
- `data/v2/staging/gbif/collected/island_taxa.csv`.

No-record islands are never converted to trait zeros.

## C. Floristic status and source-pool infrastructure

Frozen support layers include:
- native / introduced / unresolved status;
- native-nonendemic and endemic strata;
- source-pool expectations;
- family/genus composition and source-matched lineage representation.

These are assembly controls and decompositions, not proof of exact historical ancestry.

## D. Trait evidence rules

Every accepted trait row retains source provenance.

Evidence priority:
1. species-direct High/Medium;
2. trait-specific Validated Low only when direct evidence is absent;
3. unresolved otherwise.

Forbidden as primary fill rules:
- family-level inference;
- global fallback;
- Low overwriting direct evidence;
- silent downgrade of High/Medium evidence.

Trait contract:
- `docs/chapter1_trait_coverage_contract.md`.

Extraction prompt:
- `prompts/trait_evidence_extraction_v2.md`.

## Database build path

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
source batch manifests
        |
        v
cumulative source-scale integration
        |
        v
106,295 × 3 species-axis ledger
        |
        v
frozen Chapter 1 trait snapshot
        |
        + corrected 8,264-unit geography
        |
        v
current Chapter 1 analysis
```

## Historical 8,265-unit alpha1 contract

The alpha1 database and several older audits retain the original **8,265-unit** universe because their purpose is exact historical replay. They should not be rewritten to 8,264.

The current paper excludes the verified continental split component at the analysis layer. Therefore:

- **8,264 = current corrected Chapter 1 analysis universe**;
- **8,265 = historical alpha1 / frozen provenance universe**.

This distinction is intentional.

## What is not current

Do not use as the current paper database/exposure:
- historical uncorrected v14 geography;
- v13 publication locks;
- pre-v13 mechanism branches;
- broad acquisition inventories or page-hit counts;
- unreviewed machine candidates;
- historical alpha1 geography as the corrected exposure.

The current analysis is selected only by `config/chapter1_submission_current.json`.
