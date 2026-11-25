# Original GloBI v0.9 source-provider coverage — confirmed real-data finding

**11 November 2026, source-audit workflow #43008185034, 2026-06-04 frozen GloBI v0.9.** This is an *original-source evidence concentration* diagnostic, **NOT** a family-heldout morphology-to-pollinator effect.

## Verified original scan

The source workflow [island Actions #43008185034](https://github.com/zuizui0223/island/actions/runs/43008185034) **completed successfully** and wrote an anonymous receipt [artifact #13850931901](https://github.com/zuizui0223/island/actions/runs/43008185034/artifacts/13850931901). It downloaded **the exact frozen original `interactions.tsv.gz`** (SHA256 `095e57398420d1e0746e61233c796beedcfb8e03bd0c6a732aa17ae62043bd17`; version DOI `10.5281/zenodo.20546682`) rather than a live or newer GloBI release. The focal plant-genera universe is the outcome-independent **2,228 genera** pinned in original `island` artifact #10328696801.

Original source-derived receipt, *not predicted or inferred values*:

| Evidence diagnostic | Original v0.9 observed quantity |
| --- | ---: |
| All raw GloBI relationship rows streamed | **24,577,183** |
| Matched, citation-attributable flower interaction source rows | **686,908** |
| Plant genera with any qualified cited flower relationship | **2,221** |
| Genus×visitor-channel×provider×unique-reference units | **369,659** |
| Distinct *archive-preferred provider-label* groups | **62** |
| Fraction of qualified reference units in the biggest archive-preferred provider label | **95.9584%** |
| Reference IDs seen across several provider labels, excluded | 32 |
| Genera with ≥3 source-reference units before removing biggest group | **1,674** |
| Genera with ≥3 units after removing biggest group | **846** |
| Matched genera with ≥3 units on both sides | **846** |
| Median change in documented five-channel composition after deletion (total variation distance) | **0.260904** |

## What is and is not established

The largest source label accounts for **95.96%** of the citation-backed *reference units*, NOT necessarily 95.96% of unique human research studies or of all raw GloBI interactions. The provisional provider grouping picks `sourceArchiveURI` first and falls back to source DOI or source namespace+citation. Thus an original archive URI may be an **aggregator containing multiple underlying databases**. The 62 groups are *not proven 62 independent scientific study providers*.

The fact that about half of the qualified genus source information disappears under this operation, and that median five-channel documentation shares shift by 0.261, is a new empirically grounded **source-dependence and identifiability risk**. It does NOT imply floral morphology has no biologically meaningful relation to the actual pollinator community. Conversely, the previously found open-form/GloBI bird-documentation positive prediction can no longer be called source-provider independent.

## New strict follow-up audit

New `island` [PR #323](https://github.com/zuizui0223/island/pull/323) extends the scanner to report **within-largest-archive counts of distinct source namespaces, source DOIs and source citations** (aggregate only). This tests whether the dominant archive is a single underlying dataset or a *coarse source-URI aggregator*.

In the same-source rerun, it also adds a new **source-before vs largest-source-after matched genus cohort** prediction using:
- original species-level High/Medium single-state Direct floral traits from existing pinned progressive source artifact #10058653212;
- original GBIF genus–family crosswalk, dropping genus labels with conflicting families;
- frozen original 2,228-genus eligibility, ≥3 independent source references both before and after source omission;
- **each botanical family fully excluded** from training of its genus test predictions, refitting documentation effort baselines after source removal;
- identical assessed genera and family universe on both sides, with separate *genus-equal* and *family-equal* out-of-sample morphology incremental MSE.

This second predictive endpoint has **not yet returned a completed numerical result** at publication of this progress note. Source hierarchy diagnostics and family CV must be judged on their own completed artifacts, not guessed from the successful coverage-only job.

## Claim ceiling

GloBI **global documented interaction-channel reference fractions** do not represent local island pollinator guild abundances, visitation effectiveness, pollen delivery, fruit set or selection on island flowers. Original pollination-labeled rows and flower-visit-only rows remain distinct. Lepidoptera still combines moths and butterflies. The source and direct-flower datasets reuse botanical lineages, and label diversity alone does not demonstrate historical source independence.

**New scientific decision:** No island-floral-syndrome mechanism can be promoted until the source group hierarchy is validated and the **actual provider-omitted family-heldout model** establishes whether any small morphology predictive gain survives on a common, independently documented genus set.
