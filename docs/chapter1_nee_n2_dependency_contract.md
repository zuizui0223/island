# Chapter 1 NEE challenge — N2 lineage-dependency contract

Status: prospective, frozen before any N2 genus-entry outcome is inspected.

## Why this layer exists

N2 asks whether pollination-channel retention/disruption sorts source-available plant lineages according to independently measured functional dependence on that channel. The dependency variable cannot be inferred from the floral response that N2 is meant to explain.

The primary model family is therefore:

`source-matched genus entry/persistence ~ channel_retention × lineage_dependency + plant/source covariates + Baker rival + cross-classified island/genus effects`

The identifying contrast is within-island differential sorting: lineages experiencing the same island context should differ in entry/persistence according to independently measured channel dependence when that channel is retained versus disrupted.

## What counts as primary dependency

Primary D1/D2 evidence must quantify functional reproductive contribution of a named pollination channel. Accepted routes are channel-exclusion/reproductive-output experiments, visitor-specific effective-service partitioning, or repeated quantitative field measurements that directly attribute pollen transfer or reproductive contribution to the channel.

Dependency is represented on a 0–1 scale with uncertainty. Species-level evidence is propagated to source-specific genus × channel dependency hierarchically or through posterior/imputation draws. A single evidenced species may inform a genus only with explicit uncertainty; copying its point estimate to the entire genus as an exact value is prohibited.

## What does not count

GloBI edges and visitation-only networks are D3 association evidence. They can support sensitivity analyses and evidence discovery but cannot promote primary N2. Existing `pollination_guild`, pollination-syndrome templates, flower colour, flower shape, floral-architecture scores, or free-text pollination notes are not dependency measures.

Self-compatibility, mating system and autonomous selfing are also not channel dependency. They define the Baker/colonization-assurance rival and must remain a separate mechanism comparison.

## Qualification before N2 is fitted

The primary N2 fit remains closed unless all of the following are true without reading genus-entry outcomes:

- at least 3 confirmatory channels have usable dependency support;
- each eligible channel has at least 30 dependency-resolved genera;
- eligible channels jointly contain at least 150 unique dependency-resolved genera;
- each eligible channel has at least 20 islands classified retained and at least 20 classified disrupted;
- D3 evidence cannot satisfy these counts.

If this gate fails, the thresholds are not relaxed after looking at outcomes. N2 stops and the frozen H3 genus-assembly result remains the Chapter 1 mechanism ceiling.

## Robustness required if qualification passes

Primary source proxy is `geo_k5`; `geo_k10`, `geo_k20`, and `geo50_climate10` are prespecified sensitivities. Dependency uncertainty must be propagated. Leave-one-spatial-block and leave-one-source-context checks are required, as are single-species-dependency exclusion, D1-only sensitivity when supported, and explicit comparison against Baker/colonization-assurance terms.

A successful N2 supports channel-dependent lineage filtering. It does not establish realized pollination service on every island, within-lineage evolution after colonization, or adaptation of a floral trait to a named pollinator.
