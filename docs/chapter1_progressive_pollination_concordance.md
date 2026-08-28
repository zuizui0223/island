# Progressive pollination-syndrome concordance contract

## Role

Pollination-syndrome concordance is a **fixed secondary interpretation family** in the progressive Chapter 1 analysis. It is rerun after the primary pollinator-name-free plant response on every new trait snapshot.

The sampled axes are:

- `large_bee_like`
- `butterfly_like`
- `bird_like`

They are non-exhaustive floral-architecture templates, not observations or classifications of realized pollinators. Moths, flies, bats, beetles, wasps and other channels are unmodelled, not absent.

## Fixed analysis rule

The underlying trait templates remain those in `config/chapter1_pr138_pollination_syndromes.yml`. The secondary family is declared in `config/chapter1_global_branching.yml` as `sampled_guild_concordance` and is analysed across both `analysis_regime` and `biogeographic_realm` layers.

BH-FDR is applied within the predeclared
`context_layer x axis_set x stratum x support_tier` family. Response support remains `<30` not promoted, `30-49` pilot, and `>=50` confirmatory. Direct-only evidence remains a sensitivity. The same sampled-guild scores are also carried into the source-adjusted/joint-source-selection layer in the progressive workflow.

A biological result is allowed to change when trait evidence changes. The templates, FDR family, context definitions and claim ceiling are not allowed to change after inspecting updated results.

## Baseline snapshot at PR #141 HEAD `660bc02bedca8e15548bb913bc69b6b904fb8d89`

The frozen raw syndrome/pathway artifact is run `32922758001`, artifact `pr138-attraction-selfing-32922758001`, digest `sha256:5c573fdd75ebab4b85dd97cd11a064d9764eb3772c8fa37e124236eb841e6770`.

The frozen realm artifact is run `32927060210`, artifact `pr138-realm-sensitivity-32927060210`, digest `sha256:6556a7938a52547ddf7b630a822d743b2ba5ca479589716919441b96dcfa0521`.

The baseline all-native confirmatory pattern is:

- northern mid-latitude: `large_bee_like = -0.04927`, secondary-family `q = 0.0320`; bird and butterfly axes do not survive the family correction;
- tropical: bird-, butterfly- and large-bee-like concordances are all positive and supported;
- Palearctic: bird-, butterfly- and large-bee-like concordances are all negative and supported;
- Neotropical: none of the three sampled guild axes is supported;
- Nearctic: confirmatory support is insufficient;
- southern extratropical: pilot only, with bird-like and large-bee-like concordance increasing together.

The simultaneous tropical response of all three templates and the multi-template Palearctic decline demonstrate why these scores must be interpreted as overlapping floral architectures rather than pollinator identities.

## Source-pool sensitivity

Using the frozen GIFT source-pool sensitivity (`32954909953`), the qualitative architecture contrast persists:

- tropical bird-, butterfly- and large-bee-like slopes remain positive across all four source definitions;
- Palearctic bird-, butterfly- and large-bee-like slopes remain negative across all four source definitions;
- Neotropical sampled-guild slopes remain unsupported;
- northern large-bee-like and bird-like negative slopes are stable across source definitions, while the northern butterfly-like axis is less stable.

This strengthens the statement that the observed regional floral architecture is not merely removed by the tested mainland-source expectation. It still does **not** establish a pollinator mechanism.

## Lineage boundary

PR #138/#139 show that much of floral-accessibility geography is represented by genus composition. Therefore the sampled-guild concordance must remain an assemblage-level compatibility statement unless a stricter lineage-residual analysis and independent channel measurement support a mechanism.

Allowed wording:

> Regional assemblage responses are concordant with different sampled floral architectures associated with candidate pollination-service channels.

Not allowed:

- `Bombus loss caused the Palearctic pattern`;
- `birds or butterflies replaced lost pollination service in tropical islands`;
- `a positive bird-like score means bird pollination`;
- `an unmodelled guild is absent`.

## Progressive update requirement

Every new trait snapshot must recompute the sampled-guild family under the same templates, support gates and multiplicity family, alongside the primary pollinator-name-free response. Changes are recorded as evidence evolution, not treated as CI failures.
