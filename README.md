# Island — Chapter 1 corrected submission baseline

> **Current scientific surface: corrected geography baseline (24 September 2026).**
> Use `config/chapter1_submission_current.json` as the machine-readable selector.
> The uncorrected v14 and v13 surfaces are preserved only as provenance.

## Start here

1. [Current submission manuscript](submission/chapter1_current/MANUSCRIPT.md)
2. [Current submission contract](config/chapter1_submission_current.json)
3. [Corrected methods and results](docs/chapter1_corrected_submission_20260924.md)
4. [Corrected result tables](results/geography_20260924/)
5. [Replay instructions](scripts/geography_correction/README.md)
6. [Compact paper pipeline](docs/PAPER_PIPELINE.md)
7. [Historical / superseded surface index](docs/CHAPTER1_HISTORY.md)

## Current baseline

- analysis universe: **8,264 island units**;
- broad H1 union: **4,379 islands**;
- plant species: **106,295**;
- resolved trait cells: **222,688 / 318,885 = 69.83%**;
- GloPL: **2,969 experiments / 1,248 sites / 919 publications**;
- geographic exposure: minimum source-matched GSHHG 2.3.7 coastline separation on a mean-radius sphere;
- **1,113** spurious island zero distances repaired;
- **996** true continental GloPL zero-distance sites retained.

The geography repair is a post-hoc measurement correction selected as the primary exposure baseline. It is not a new prospective confirmation.

## Current H1–H4 result spine

### H1 — Pattern
A seven-response floral/reproductive isolation vector is jointly supported in all four geographic strata in both evidence scopes. Individual traits are not uniformly positive: southern shallow/open tube is negative in the corrected analysis.

### H2 — Conditional decomposition
Reproductive assurance is separated from additional floral responses. Primary selfing-adjusted accessibility is FDR-supported in northern high latitudes and the tropics. Tropical Direct-only accessibility is nominally positive but **not** FDR-supported (`q=0.1196`). Colour and colour × architecture responses remain region dependent.

### H3 — Ecological pressure
Experimental pollen limitation increases with corrected geographic isolation:
`beta=0.09191`, `SE=0.03806`, two-sided `p=0.01575`.

### H4 — Functional compatibility
Exact-species post-hoc functional triangulation links both H2 trait families to lower current pollen limitation:

- reproductive assurance: `beta=-0.29830`, `p=0.00396`;
- generalized accessibility: `beta=-0.29566`, `p=0.02187`.

H4 remains an association, not causal mediation.

## Database

The trait database denominator remains **106,295 species × 3 raw axes = 318,885 species-axis cells**:

1. flower colour;
2. floral structural complexity;
3. reproductive assurance.

Trait construction and provenance are documented in [docs/DATABASE_BUILD.md](docs/DATABASE_BUILD.md).

The historical alpha1 database bundle retains its original **8,265-unit** contract for exact replay. That historical acquisition universe is not the corrected Chapter 1 analysis universe.

## Active repository map

```text
README.md
config/chapter1_submission_current.json
submission/chapter1_current/
  MANUSCRIPT.md
  FIGURE_CAPTIONS.md
  COVER_LETTER_DRAFT.md
  PACKAGE_MANIFEST.md
  SUBMISSION_CHECKLIST.md
docs/chapter1_corrected_submission_20260924.md
docs/PAPER_PIPELINE.md
docs/DATABASE_BUILD.md
results/geography_20260924/
scripts/geography_correction/
src/island_v2/
tests/
legacy/
```

## Superseded v14 surface (provenance)

The reproduced v14 surface is retained for audit and historical replay only:

- `config/chapter1_v14_canonical_result_lock.json`;
- `docs/chapter1_manuscript_full_v14_reordered_hypotheses_20260918.md`;
- `docs/chapter1_submission_freeze_v14_20260918.md`;
- `.github/workflows/run-chapter1-v14-reordered-hypotheses.yml`.

It must not be used as the current submission result surface.

Superseded v8–v14 publication/promotion/render workflows are retained as **manual-only historical replay**; `tests/test_chapter1_historical_workflows_manual_only.py` prevents automatic push/PR/schedule triggers from being reintroduced.

## Frozen v13 parent paper surface

v13 remains immutable parent provenance beneath v14. Pre-v13 publication branches are archived under `legacy/chapter1-pre-v13/`.

## Claim ceiling

The current submission may claim recurrent multivariate response, conditional decomposition, an independent pollen-limitation gradient, and post-hoc functional compatibility.

It must not claim:

- historical causal mediation from pollen limitation to trait evolution;
- global pollinator abundance or visitation decline;
- a universal named pollinator mechanism;
- uniform positive change in every H1 trait;
- tropical Direct-only H2 accessibility as FDR-supported after correction;
- within-lineage evolution rather than assemblage composition;
- corrected geography as prospective confirmation.
