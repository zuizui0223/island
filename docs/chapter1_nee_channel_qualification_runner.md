# Chapter 1 NEE channel-qualification execution boundary

## Purpose

This workflow is the last execution layer before N1. It does not fit N1 and does not
inspect plant outcomes. It combines only:

1. the completed primary `geo_k5` island × channel source-availability artifact; and
2. the five completed exact-island channel-observation artifacts.

It then calls the already merged `chapter1_nee_channel_inputs` validator to project
`retained / disrupted / unresolved` states and compute the frozen per-channel support
qualification receipt.

## Frozen upstream requirements

The trigger file must specify two completed workflow run IDs:

```text
source_proxy_run_id=<run id>
island_observation_run_id=<run id>
```

The workflow requires exactly one source-proxy artifact named
`chapter1-nee-source-proxies-from-full-scan-<run id>` and five observation artifacts
named `chapter1-nee-island-search-observation-<channel>-<run id>`.

A missing channel artifact is a hard stop. The workflow never substitutes another
channel, source mode, observation route, or support threshold.

## Projection rule

The existing frozen rule remains authoritative:

- source available + detected -> retained;
- source available + adequate non-detection -> disrupted;
- source available + insufficient effort/unresolved -> unresolved;
- unresolved source -> unresolved.

The Search-specific audit is checked again before projection: any
`adequate_non_detection` row must have `search_complete=true`. Thus a truncated Search
result cannot enter N1 as a disrupted channel.

## Qualification

The workflow writes the existing confirmatory/pilot/not-qualified support receipt.
Confirmatory support remains the only N1-gate-eligible tier. The workflow reports the
number of confirmatory channels but does not change the minimum-three-channel N1 rule
and does not fit a reduced rescue model.

## Output artifact

`chapter1-nee-channel-qualification-<workflow run id>` contains:

- the combined island-observation table;
- `nee_source_availability.csv`;
- `nee_island_channel_observation.csv`;
- `nee_projected_channel_state.csv`;
- `nee_channel_qualification_receipt.csv`;
- `nee_channel_qualification_receipt.json`.

These are the direct inputs to the separately frozen N1 execution wrapper.
