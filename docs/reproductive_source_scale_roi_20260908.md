# Reproductive source-scale ROI — 2026-09-08

This table replaces page/species hit counting with **realized strict reproductive cells per completed source packet**.  The comparison basis is the recoverable public restart checkpoint; it is not a reconstruction of the historical private TRY + Wave55 state.

| Source | Source-scale status | Exact unresolved overlap before semantic gate | Realized / gated strict direct opportunity | Decision |
|---|---|---:|---:|---|
| Ackerman et al. Orchidaceae database, Zenodo 10.5281/zenodo.14601785 | complete workbook extracted and source methods reviewed | 192 reference-backed SI/SC species after strict mapping | **192** (188 SC, 4 SI) | **integrate as next formal source batch** |
| Freyman & Höhna Onagraceae dataset | complete source extracted | 20 | **2** | demote; retain two direct candidates, do not spend more acquisition time here |
| Delaney & Igić Fabaceae block in Meyer review workbook | complete Meyer block screened; original-source evidence-type validation not yet completed | 15 | **3 non-conflicting explicit SC/SI rows in the Meyer compilation** | lower priority than Orchid; source article/supplement must prove experimental basis before promotion |
| Ramírez et al. 2022 block in Meyer review workbook | complete Meyer block screened | 27 | **0 strict values in the compilation block** (`mating_system=review` only) | do not extract species by species from this block; only revisit if the original structured source exposes explicit breeding states |
| Razanajatovo block in Meyer review workbook | complete Meyer block screened | 11 | not a new source: dedicated Razanajatovo adapter already exists in the direct ledger | do not reacquire merely because it appears in a later compilation |
| Moeller et al. 2017 Dryad | previously source-scale extracted | current overlap now small | already represented in direct evidence | no routine reacquisition |
| Ferrer et al. 2024 SI database | dedicated full-source adapter already present | already represented | already represented | use only for conflict/novelty audit, not duplicate acquisition |
| Meyer et al. 2026 Dryad | dedicated full-source adapter already present | already represented | already represented | source for provenance/ROI discovery, not duplicate acquisition |
| Zell et al. BSdb | dedicated source-scale adapter already present | largely overlaps prior direct ledger | already represented where novel | use for exact current-gap audit before any new ingestion |

## Orchid semantic gate

The Ackerman paper explicitly defines SI/SC from controlled hand-pollination comparisons.  The completed 2024 Zenodo workbook retains species-level literature references, so `SI`/`SC` rows with exact fixed-universe names and non-empty references form one valid source-scale strict packet.

Three tempting fields are **not** crossed into strict traits:

- `Mixed mating` is a mixed pollination-system label (chasmogamy + autonomous reproduction), not a population-genetic `mating_system` estimate;
- `autonomous selfing/agagamospermy` conflates autonomous sexual selfing with apomixis;
- `evidence for selfing` is an evidence-strength field, not a trait state.

Therefore the batch contains only `self_incompatibility` and cannot train genus rules in the restart lane.

## Validated-Low lane

The last hash-verified all-evidence rule audit remains the Wave53 `current_min3` trait-specific audit. Restart/source packets are marked `genus_rule_training_allowed=false`, so they do not silently create new rules.  A batch-level Low headroom audit is run against the post-batch coverage using exactly those current thresholds; min2/relaxed rules, family inference, global fallback and genus×axis shortcuts remain excluded.

## Formal integration cadence

The Orchid packet is the first post-pivot packet large enough to justify a formal run.  It is integrated once through a committed batch manifest. Evidence-row edits do not trigger integration. Future formal runs require another completed multi-species source packet and a new manifest.
