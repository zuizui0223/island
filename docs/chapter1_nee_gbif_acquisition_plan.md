# Chapter 1 NEE — prospective GBIF background acquisition

Status: **pre-occurrence acquisition plan**. This document does not contain N1 outcomes.

## Analytical backgrounds already frozen

The observation policy remains authoritative:

- `bombus`: Apidae, using the existing Bombus observation policy;
- `non_bombus_bees`: Andrenidae, Apidae, Colletidae, Halictidae, Megachilidae, Melittidae, Stenotritidae;
- `lepidoptera`: Lepidoptera;
- `flower_visiting_birds`: Aves;
- `diptera`: Diptera.

## Acquisition proxies

The GBIF acquisition taxon may be broader than the analytical background only when a deterministic predeclared post-filter restores the frozen background before effort is calculated.

| N1 channel | GBIF acquisition taxon | rank | post-filter before effort |
|---|---|---|---|
| non_bombus_bees | Hymenoptera | ORDER | keep only the frozen seven bee families |
| lepidoptera | Lepidoptera | ORDER | none; acquisition taxon equals frozen background |
| flower_visiting_birds | Aves | CLASS | none; acquisition taxon equals frozen background |
| diptera | Diptera | ORDER | none; acquisition taxon equals frozen background |

Bombus is not folded into this new campaign. Existing Apidae acquisition/diagnostics remain the primary Bombus observation route.

## Spatial acquisition

All campaigns reuse `island_v2.gbif_blocks` and the frozen global island geometry. Query catchments may be buffered for acquisition only. Every occurrence must be reassigned to the original exact island polygon before it can enter an observation state. Buffer-only mainland/coastal records are excluded.

## Frozen output boundary

The acquisition stage may emit exact-island occurrence rows and audit receipts. It must not:

- inspect the focal plant trait response;
- alter the functional target catalog;
- classify source availability;
- fit N1;
- change the effort thresholds;
- use target detections to alter the broad background definition.

Target-channel detection is applied only later by the already-frozen channel-observation classifier using the independently frozen confirmatory functional catalog.
