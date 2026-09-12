# Chapter 1 NEE challenge — source-set proxy contract

Status: `prospective_pre_pollinator_occurrence`

N1 reuses the outcome-blind mainland source assignments already frozen for the Chapter 1 source-pool sensitivity. No new source region is chosen from pollinator retention or plant-trait outcomes.

## Fixed source assignment artifact

The source assignment is pinned to workflow run `32954909953`, artifact `pr138-source-pool-sensitivity-32954909953`, digest `sha256:51cff9087a40b76655ff043d59aceabf25800a6029ed005a77b8db30e78345d2`, file `result/island_source_assignments.csv.gz`.

That artifact predates the NEE channel challenge and explicitly records `source_assignment_uses_island_traits = false`.

## Primary source proxy

Primary N1 uses `geo_k5`: the five geographically nearest eligible GIFT mainland source-region representative points.

The five selected source entities are treated as a candidate source set, not as five independent observations.

For each island × channel:

- `available` if at least one of the five source entities is independently classified `available`;
- `structurally_absent` only if all five source entities have explicit accepted structural-absence evidence;
- `unresolved` when there is no positive evidence and at least one of the five source entities remains unresolved.

Thus partial absence does not become source absence.

All ranks 1–5 must be present exactly once. An incomplete or unexpected rank set fails closed.

## Source context used in the N1 model

The predeclared source-region fixed-effect label is the rank-1 entity in the frozen `geo_k5` set. This is chosen by geographic rank only, before any channel state is inspected. A channel detected only in rank 3, for example, does not change the source-context label.

The complete five-entity source set is retained in an audit table for every island × channel.

## Prespecified source-definition sensitivities

The following already-frozen PR138 source modes are retained as source-definition sensitivities:

- `geo_k10`;
- `geo_k20`;
- `geo50_climate10`.

They must be reported if N1 is run, but they cannot replace a failed primary `geo_k5` result after outcome inspection.

## Claim boundary

These source sets are accessibility proxies, not reconstructions of historical colonization routes. A positive N1 result may show that source-available pollination channels differ in their isolation-retention curves. It does not identify the historical mainland source of a plant or pollinator lineage.
