# Minimum endemic-colour acquisition queue — closed 2026-08-24

## Decision

The targeted endemic-colour campaign is complete and must not be expanded.

Latest locked support snapshot: run `32695319042`.

| regime | endemic-status ceiling | direct-supported colour islands | gap to 50 | decision |
|---|---:|---:|---:|---|
| northern mid-latitude | 59 | **50** | 0 | **model ready; stop acquisition** |
| tropical | 70 | **50** | 0 | **model ready; stop acquisition** |
| southern extratropical | 13 | 10 | 40 | **status limited** |

The historical 23-species queue was a bounded stopping queue, not a standing acquisition target. Successive reviewed direct-evidence waves closed the northern and tropical gaps. Empty minimum-queue outputs in run `32695319042` are therefore the expected stopping state.

## Historical trajectory

The campaign started from:

- northern endemic colour: 31 direct-supported islands, gap 19;
- tropical endemic colour: 46 direct-supported islands, gap 4.

It was intentionally restricted to species that could add currently unsupported islands rather than maximize global species-axis coverage. Reviewed evidence was promoted only through the existing species-direct provenance and ontology gates in PR #132.

## Frozen stop rule

1. Do not acquire additional northern or tropical endemic colour solely to increase coverage.
2. Do not continue down the old ranked queue after 50 supported islands has been reached.
3. Do not treat southern endemic colour as a trait-acquisition problem; the current ceiling is floristic-status coverage.
4. Any future colour acquisition must be justified by a new, explicit analysis-support gap rather than by global missingness.

## Next analysis boundary

The remaining acquisition question concerns endemic **floral form** and **SI/SC**. These outcomes are handled separately in `docs/endemic_trait_recoverability_pilot_20260824.md` and may only be acquired to the 30-island pilot boundary before reanalysis.
