# Original GloBI source-provider audit: reproduction route and status

**11 Nov 2026 — fd PR #14 — source-provenance acquisition, NOT an empirical provider-independent flower/pollinator effect.**

## Main unresolved inferential issue

The original 484-genus / 88-unambiguous-family Direct flower morphology analysis predicts **GloBI DOCUMENTED global reference-share composition**, not real local island pollinator abundances or pollen transfer. Only the genus summary was archived in the matched-score analysis, so its original source providers and bibliography keys are absent. No robust claim of independence from one high-volume source publication is currently justified.

## Source selection has been fixed and original processing started

The source-side GloBI study is now isolated as [island Draft PR #323](https://github.com/zuizui0223/island/pull/323), branch `analysis/fd-globi-provider-v09-20261111`, with GitHub Actions [run #43008185034](https://github.com/zuizui0223/island/actions/runs/43008185034) triggered. The script:

1. Downloads the already frozen **pre-island-outcome** source predictor, 2,228 original plant genera, from `island` workflow run `34790842133` / artifact `10328696801`. It does not select GloBI genera based on fd morphology or the previously favorable source classification.
2. Downloads **exact GloBI v0.9 (2026-06-04) `interactions.tsv.gz`**, DOI `10.5281/zenodo.20546682`, and hard-checks SHA256 `095e57398420d1e0746e61233c796beedcfb8e03bd0c6a732aa17ae62043bd17`. A newer GloBI v0.11 exists but would be a **different training source** requiring its own frozen analysis and cannot silently replace v0.9.
3. Retains the four prospectively frozen RO pollination/flower-visitation orientations (visitor can occur on source or target side), real Plantae and the previously defined five visitor taxonomy groups.
4. Requires a **source provider ID and a separate bibliographic citation ID**. No invented provider, study, citation, or source-independent count is permitted when those original fields are missing.
5. Deduplicates genus×channel×provider×reference units, detects the same bibliographic reference appearing across multiple providers and removes such cross-provider duplicates from both providers.
6. On the *same genus and channel evidence*, reports only an anonymous original-source coverage receipt: number of providers, fraction of reference units belonging to the largest provider, number of genera with >=3 citation-backed units before and after excluding that provider, and median change in global documented five-channel composition.

**Current status (at source report creation): the input/implementation and synthetic checks are complete; the full 2.7GB original source scan has not yet reported a source-provider result.** The job status must be checked before interpreting coverage metrics. The original provider/reference IDs, genus identities and original rows remain outside the fd public repository and outside released artifacts.

## Scientific decision gate after the source audit

A later *scientific* provider-independent test would require the complete same-subset **484 genus × 5 reference-channel composition** before/after deleting the largest actual source provider, on identical documented plant genera with >=3 qualifying references in both, the same Direct morphology predictors, and the same 88-family whole-family leaveout training design. It must report:

- Matched genus/family retention and group-wise dropout. A huge reduction in support may make the test **unidentifiable**, not 'negative'.
- The source-reference-effort-only baseline **retrained after provider omission**, as well as the open-radial-form incremental MSE gain; never reuse the old baseline on changed targets.
- Genus-weighted versus equal-family mean gains separately, and conditional uncertainty that acknowledges family and documentation source overlap.
- Separate visitor-channel results only if the corresponding source-support criteria are met. The GloBI Lepidoptera source still does **not** split butterflies and moths.
- Distinguish `flower_visit` and `pollination` evidence, retaining a primary visitor-reference share model and stricter source-evidence sensitivities without mislabeling either as effective seed-set evidence.

If the source scan fails a verified download, missing field or provenance check, stop and publish the blocker. Do not replace the original pinned source with a live API or newer v0.11 to obtain a favorable result. No independent-source claim may be added until a real original provider-independent family CV is complete.

## Code / audit links

- [Source side PR #323](https://github.com/zuizui0223/island/pull/323) — `src/island_v2/fd_globi_provider_v09.py`, CLI script, five synthetic tests and manual/push-action workflow.
- [Source-side GitHub Actions #43008185034](https://github.com/zuizui0223/island/actions/runs/43008185034) — pinned GloBI download/sha verification and receipt job.
- Current fd `src/fd/globi_provider_provenance_gate.py` — refuses original-provider independent claims if provider/ref fields are absent.
- [GloBI v0.9 Zenodo original](https://zenodo.org/records/20546682) — documentation requires attribution to original contributing datasets, with taxonomic/geographical/source biases noted.

**The narrow fd finding remains unchanged:** coarse open-radial floral form weakly predicts GloBI's **global literature-reference** visitor-channel composition across selected unseen botanical families. This has not been shown to represent actual island-local pollination guilds, successful pollen delivery or independent original publication providers.
