# Chapter 1 NEE challenge — channel observation contract

Status: `prospective_pre_channel_occurrence`

This contract prevents the non-Bombus channels from becoming raw presence/absence proxies.

## Target taxa

Target taxa come only from the independently frozen functional channel catalog. A taxon must first have explicit flower-interaction evidence before its island occurrence can count as a pollination-channel detection.

Primary N1 uses only `confirmatory` catalog taxa. Sensitivity-only taxa cannot promote N1.

## Broad background groups

Observation effort is measured using a broad taxonomic group, not the target catalog itself:

- Bombus: Apidae (existing frozen Bombus policy);
- non-Bombus bees: Andrenidae, Apidae, Colletidae, Halictidae, Megachilidae, Melittidae and Stenotritidae;
- Lepidoptera: Lepidoptera;
- flower-visiting birds: Aves;
- Diptera: Diptera.

Thus target non-detection is interpretable only when the island demonstrably has enough observation of the relevant broad group.

## Primary effort gate

The same gate is used for all newly added non-Bombus channels:

```text
>= 50 background records
>= 3 spatial units
>= 2 temporal units
>= 2 datasets
latest background record <= 20 years old relative to 2026
```

These thresholds are frozen before channel outcomes are opened. They are not optimized separately for birds, butterflies, flies or bees.

Prespecified sensitivity gates are 25-record/liberal and 100-record/strict variants. Sensitivity gates can diagnose dependence on observation policy but cannot rescue a failed primary N1 result.

## Observation states

`detected` requires at least one quality-filtered exact-island occurrence of a confirmatory catalog taxon. A detection remains a detection even when background effort is sparse.

`adequate_non_detection` requires zero target records **and every primary effort criterion**.

`insufficient_effort` is used whenever a zero-target island fails one or more effort criteria.

Missing or unresolved records are never silently converted to absence.

## Establishment status

Known introduced/alien target records do not create a primary natural-retention detection when establishment status is available. Unknown establishment status is retained with an explicit audit flag. A sensitivity analysis may include all non-captive records.

## Source availability

Positive source-region evidence is sufficient to classify a channel as `available`.

`structurally_absent` is much stricter: it requires accepted curated biogeographic evidence explicitly supporting absence of the channel from the source region. Zero GBIF records alone can never create structural absence.

Without positive source evidence or explicit structural-absence evidence, the source state remains `unresolved`.

## Claim ceiling

An island `detected` state means:

> a taxon independently known to use flowers and assigned to the frozen functional channel is geographically present on the island.

It does **not** establish realized flower visitation, per-visit effectiveness, effective service or functional replacement on that island.

Likewise, `disrupted` means source-available but adequately non-detected under the frozen observation policy. It does not prove historical extinction or zero abundance.

This is enough for N1 partner-side geographic filtering. N2 still requires independent lineage dependency before a plant-assembly mechanism can be claimed.
