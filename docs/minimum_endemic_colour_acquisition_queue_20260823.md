# Minimum endemic-colour acquisition queue — 2026-08-23

## Purpose

This checkpoint converts the model-support rule into an explicit stopping queue.
It does **not** ask how many unresolved species can be filled. It asks which
species need a new direct flower-colour observation to bring the corrected
regional-endemic stratum to the 50-island continuation target.

Real-data run:

- workflow run: `32614335237`;
- artifact: `minimum-endemic-acquisition-queue-32614335237`;
- status: WCVP-corroborated run `32559322028`;
- direct trait evidence: PR #132 run `32551742699`;
- covariates: run `29228212586`.

Upstream ranking is fixed as:

1. unsupported endemic-island gain;
2. whether the island extends the current supported isolation 5–95% envelope;
3. GBIF island-record count as a coarse recoverability proxy;
4. stable species-name tie break.

The minimum queue then stops as soon as the ranked candidates cover the island
gap. A queue entry is a **search target only** and never implies a trait value.

## Northern mid-latitude endemic colour

Current direct support: 31 islands. Status ceiling: 59. Gap to the 50-island
target: **19 islands**.

The greedy stop rule selects exactly **19 species / 19 new islands**:

| order | species | GBIF records | extends current isolation envelope |
|---:|---|---:|---|
| 1 | *Rostraria azorica* | 142 | yes |
| 2 | *Adinandra ryukyuensis* | 66 | yes |
| 3 | *Odontarrhena lesbiaca* | 38 | yes |
| 4 | *Cirsium tanegashimense* | 34 | yes |
| 5 | *Sasa jotanii* | 22 | yes |
| 6 | *Mammillaria evermanniana* | 18 | yes |
| 7 | *Centaurea thasia* | 14 | yes |
| 8 | *Campanula sciathia* | 11 | yes |
| 9 | *Centaurea xylobasis* | 9 | yes |
| 10 | *Limonium aphroditae* | 7 | yes |
| 11 | *Astragalus kawakamii* | 1 | yes |
| 12 | *Oxytropis kunashiriensis* | 1 | yes |
| 13 | *Fritillaria rhodia* | 86 | no |
| 14 | *Cirsium ishizuchiense* | 75 | no |
| 15 | *Erysimum naxense* | 56 | no |
| 16 | *Hedysarum sachalinense* | 41 | no |
| 17 | *Cirsium umezawanum* | 38 | no |
| 18 | *Limonium retusum* | 29 | no |
| 19 | *Viola cephalonica* | 20 | no |

The first 12 targets preferentially add islands outside the current direct
colour isolation envelope. This is important because the original failure mode
was coverage declining along the isolation/endemicity gradient; the queue is not
allowed to improve only the already well-covered centre.

No additional northern endemic-colour species should be searched merely to
increase global coverage before these 19 targets are resolved or explicitly
shown unrecoverable.

## Tropical endemic colour

Current binary-colour direct support: 46 islands in the targeting contract.
Status ceiling: 70. Gap to 50: **4 islands**.

The minimum queue is only **4 species / 4 new islands**:

| order | species | GBIF records | extends current isolation envelope |
|---:|---|---:|---|
| 1 | *Dracaena fernaldii* | 60 | yes |
| 2 | *Ixora mooreensis* | 28 | yes |
| 3 | *Myrcia stenocarpa* | 28 | yes |
| 4 | *Myrsine fasciculata* | 25 | yes |

All four selected islands extend the current 5–95% direct-support isolation
envelope under the upstream ranking. This is therefore a much smaller and more
valuable acquisition problem than a global tropical-colour campaign.

## What is deliberately not queued

### Endemic floral form and SI/SC

Northern and tropical endemic form / SI-SC remain below the 30-island pilot
support threshold. They are **recoverability problems first**, not automatic
50-island acquisition campaigns. The next step is to estimate how many currently
unsupported islands can actually be unlocked by recoverable direct evidence.

### Southern endemic traits

The corrected endemic-status ceiling is only 13 southern extratropical islands.
Trait searches cannot solve this ceiling, so southern endemic trait acquisition
is not queued.

### Bombus occurrence

No Bombus-occurrence campaign is queued to rescue the environmental-hypervolume
mechanism. The northern lineage-residual checkpoint failed to recover the
predicted Bombus-deficit pathway.

## Stopping rule

After every accepted direct colour record:

1. rebuild the direct-evidence ledger;
2. recompute which targeted endemic island has become colour-supported;
3. remove all remaining species for that island from the active queue;
4. stop the regional colour campaign immediately when 50 supported islands are
   reached;
5. rerun status-aware category and broad-colour analyses before any further
   acquisition.

If a queued species cannot be recovered from the declared evidence sources, mark
that target as `unrecoverable_current_sources` and move to the next ranked species
from the full target artifact. Do not weaken evidence rules to force the queue to
complete.
