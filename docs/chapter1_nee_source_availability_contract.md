# Chapter 1 NEE challenge — source availability contract

Status: `prospective_pre_source_occurrence`

N1 distinguishes a channel that was available in a source region and then not retained on an island from a channel that was never part of the source-region pollination pool. This distinction must be made without using island outcomes.

## Positive source evidence

A source region is `available` for channel `g` when at least one **confirmatory functional-channel catalog taxon** has accepted positive source-region evidence. Accepted positive evidence is restricted to:

- curated native distribution;
- curated regional checklist;
- quality-filtered source-region occurrence.

A single accepted positive record is sufficient to establish presence. Presence does not require an observation-effort absence gate.

Sensitivity-only channel taxa cannot establish primary N1 source availability.

## Structural absence

`structurally_absent` is intentionally much harder to assign. It requires accepted evidence explicitly supporting biogeographic exclusion of the channel from the source region, such as an authoritative native-range exclusion.

The following can never establish structural absence:

- zero GBIF records;
- island non-detection;
- climatic unsuitability alone;
- floral phenotype or pollination-syndrome inference.

If no accepted positive evidence and no explicit structural-absence evidence exists, the state is `unresolved`, not absent.

## Evidence conflict

If accepted positive evidence and accepted structural-absence evidence coexist for the same source-region × channel pair, neither receives precedence. The pair is marked `unresolved`, and the conflict is emitted as an audit table. It cannot enter confirmatory N1 until reconciled on evidence grounds.

## Island-to-source assignment

A decisive source state requires a pre-existing accepted island-to-source-region assignment. Missing or pending assignments make all channels for that island `unresolved`.

The assignment itself cannot use island pollinator occurrence or Chapter 1 plant outcomes.

## Output

The source-region state is expanded to exactly one row per island × frozen N0 channel, in the schema already required by `chapter1_nee_channel_inputs_v1`:

```text
island_id
source_region_id
channel_id
source_state
evidence_id
evidence_type
source_citation
source_url
review_status
```

The builder never reads focal plant traits, island channel observations or N1 model outcomes.
