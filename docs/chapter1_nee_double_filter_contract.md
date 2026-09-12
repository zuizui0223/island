# Chapter 1 NEE double-filter challenge

Status: `prospective_preoutcome_qualification`

Issue: #184

## Purpose

This challenge tests one additional causal arrow without reopening the frozen Chapter 1 H1-H4 results.

Frozen Chapter 1 currently identifies:

```text
source separation / isolation
        -> biogeographic response branching
        -> genus-level assemblage representation
        -> floral / reproductive response
```

For the broad Palearctic response, the plant-side signal is represented at genus assembly depth. The frozen analysis does not identify why source-available genera are differentially represented.

The NEE challenge asks whether independent pollinator-side geography closes that missing arrow:

```text
                         geographic isolation
                                  |
                 +----------------+----------------+
                 |                                 |
                 v                                 v
        plant-source filter              pollinator-channel filter
                                                    |
                                       x lineage dependency
                 |                                 |
                 +----------------+----------------+
                                  |
                                  v
                             genus assembly
                                  |
                                  v
                    emergent floral/reproductive
                         island response
```

Failure of this challenge leaves the frozen Chapter 1 claim intact. It does not trigger a refit, relabelling, alternative-channel search or renewed trait-acquisition campaign.

## N0 — qualification before outcomes

N0 is not a biological Result. It is the anti-overfitting gate.

The machine-readable contract is `config/chapter1_nee_double_filter.yml`; `src/island_v2/chapter1_nee_contract.py` validates the frozen decisions.

### Fixed channel ontology

The first challenge uses five taxonomic-functional channels:

- `bombus`
- `non_bombus_bees`
- `lepidoptera`
- `flower_visiting_birds`
- `diptera`

These are partner-side channels, not floral syndromes. `large_bee_like`, `butterfly_like` and `bird_like` remain plant-architecture concordance scores only and cannot qualify N1 or N2.

A channel may be excluded only by a predeclared data-quality failure before N1 outcomes are inspected. Exclusion remains visible in the qualification receipt. Channels cannot be added after N1 outcome inspection.

### Source availability and island state

Source-region state is defined independently of focal plant traits:

- `available`
- `structurally_absent`
- `unresolved`

Island channel state is:

- `retained`
- `disrupted`
- `unresolved`
- `structurally_absent`

`disrupted` is only available when the channel exists in the source region and observation effort is adequate. Raw absence, raw GBIF richness and climatic compatibility alone are insufficient.

### Lineage dependency

The N2 dependency variable is channel-specific and outcome-independent. Preferred evidence is, in order:

1. single-visit effectiveness by channel;
2. channel-attributed pollen transfer or seed set;
3. replicated exclusion/access experiments that resolve the channel;
4. curated effective-pollinator interactions;
5. repeated direct visitation or explicit literature statements as lower-tier evidence.

Flower colour, form, tube depth, symmetry and Chapter 1 pollination-architecture scores are prohibited inputs to dependency. Unresolved dependency is not coded as zero.

Self-compatibility or reproductive assurance is not used to define channel dependency. It remains a required Baker/colonization-assurance rival.

## N1 — pollinator-side geographic filter

Question:

> Do independently measured pollination channels show different source-conditioned retention/disruption curves with isolation?

Conceptual model:

```text
channel retention
  ~ isolation x channel
  + source availability
  + island area
  + climate / habitat
  + observation effort
  + regional / spatial structure
```

N1 passes only if channel retention is non-exchangeable after source availability and observation effort are handled. A result based only on raw presence/absence or richness cannot pass.

N1 failure stops the challenge before N2.

## N2 — channel-dependent lineage filtering

Question:

> Does channel geography explain H3 genus assembly?

Primary estimand:

```text
source-matched genus entry / persistence
  ~ channel retention x lineage dependency
  + plant/source covariates
```

The preferred design is cross-classified by island and lineage so that the key information is within-island differential loss: among lineages that could have arrived from the source pool, are those more dependent on a disrupted channel selectively under-represented?

Required rivals are:

- plant-only geographic/source filtering;
- Baker/colonization assurance.

A distance x SC result alone is not N2. N2 requires stable added information from `channel retention x lineage dependency` beyond the plant-only model.

N2 failure terminates H5a. No alternative-channel rescue is allowed.

## N3 — held-out double-filter prediction

Two model families are frozen:

```text
P  = source flora + plant geographic filter + area + climate
     + floristic status + lineage baseline propensity

DP = P + channel retention x lineage dependency
```

Primary holdout is leave-one-archipelago-out. Random island splits are not primary because nearby islands can share geography, source pools and observation structure.

Primary endpoint:

```text
held-out source-available genus entry
```

Primary metric:

```text
mean log loss
```

Secondary endpoints include Brier score, assemblage Bray-Curtis error and the downstream floral/reproductive response vector. The trait vector remains secondary: it cannot be used to construct channel retention or dependency.

N3 passes only if DP improves out-of-sample prediction and the gain is not driven by one archipelago or one channel.

## Routing

| Outcome | Route |
|---|---|
| N1 + N2 + N3 robust | Nature Ecology & Evolution first submission |
| N1 + N2 strong, N3 weak | Ecology Letters / Global Ecology and Biogeography |
| N1 only | frozen Chapter 1 route |
| N2 fails | retain H3: genus assembly explains the broad response; cause unidentified |

## Scope boundary with Chapter 2

This challenge stays at the assembly level. It does not import Izu or Campanula field data, global seed-set measurements or H5b within-lineage causal inference.

Chapter 1 challenge:

```text
partner geography -> lineage assembly -> emergent syndrome
```

Chapter 2 (`izu-core`):

```text
realized visitor community -> functional matching -> effective service
-> plant-specific reproductive / morphological response
```

The two chapters therefore remain complementary rather than competing for the same mechanism claim.
