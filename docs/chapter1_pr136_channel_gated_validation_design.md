# Chapter 1 validation design anchored to PR #136

## Core hypothesis

This design treats PR #136 as the scientific mainline.

> **Island isolation does not create a universal floral syndrome; it creates a syndrome only when it removes a pollination channel that the regional flora could otherwise use.**

The Chapter 1 direct hypothesis remains **biogeographically contingent island floral filtering**. `Channel-gated island assembly` is the causal explanation to be evaluated with increasingly direct channel evidence.

The key distinction is:

```text
retained != disrupted != structurally absent
```

A pollination channel that was never available in the regional/source flora cannot be described as lost. Likewise, a lost channel that is effectively replaced need not produce the same plant filtering as an unreplaced loss.

## Causal decomposition

For island `i`, pollination channel `c`, and plant trait axis/category `t`:

```text
syndrome-forming pressure(i,c,t)
  ~ source_channel_exposure(i,c)
  x channel_disruption(i,c)
  x (1 - effective_replacement(i,c))
  x plant_functional_dependence(t,c)
```

This is a causal hypothesis, not an already identified regression equation.

Distance to continent remains a composite geography axis. It can encode dispersal limitation, connectivity, source supply and assembly history. It is never interpreted as a pure pollinator-loss treatment.

## Four-state natural experiment

The hypothesis is tested by separating four biological states rather than comparing `pollinator present` versus `pollinator absent`.

| State | Source/channel condition | Island condition | PR136 prediction |
|---|---|---|---|
| **Retained** | channel available to source flora | channel retained on island | little or no channel-specific syndrome pressure |
| **Disrupted** | channel available to source flora | channel reduced/lost with isolation | strongest predicted channel-specific filtering |
| **Structurally absent** | channel not available to source flora | channel absent by construction | absence alone must not create that channel-specific syndrome |
| **Replaced** | source channel disrupted | alternative channel provides comparable effective service | syndrome weakened, absent or redirected relative to unreplaced disruption |

`Structurally absent` is therefore a falsification regime, not a deficit category.

## Response definition: residual trait filtering, not raw trait composition

The primary biological response is not a raw island trait proportion.

For each island, floristic-status stratum and directly supported trait state:

```text
residual_trait(i,t)
  = observed_direct_trait_share(i,t)
    - same_genus_null_expectation(i,t)
```

The same-genus null preserves observed island species membership and genus composition while reassigning direct trait states within genera. It asks whether a geography/channel association remains beyond the trait composition expected from lineage turnover.

This is evaluated separately for:

- `all_native`
- `native_nonendemic`
- `endemic` where support permits
- `introduced` as a secondary negative-control stratum

Endemicity is also analysed separately as a response, not inserted as a nuisance covariate into a process that may generate it.

## Trait domains

### Required syndrome domains

A channel-specific **syndrome** is not declared from one trait. The minimum multivariate target is:

1. **access / mechanical architecture**
   - flower size
   - tube depth
   - symmetry
   - floral form
2. **reproductive assurance**
   - SI/SC
   - autonomous or delayed selfing
   - mating system
3. **display / colour**
   - flower colour is secondary and cannot establish a syndrome alone

Axes not yet passed through status, genus-fixed, support, FDR and coverage gates remain discovery/supporting axes.

### PR136 directional contrast for the Northern simplification/reproductive-assurance candidate

Expected positive direction under unreplaced loss of an important large/long-tongued channel:

```text
more small / very small
more open or shallow access
more actinomorphic / generalized-access architecture
more SC
more autonomous or delayed selfing
more predominantly/obligately selfing mating states
```

Expected negative direction:

```text
large / very large flowers
deep tubes
zygomorphic or mechanically restrictive forms
specialized deep/tubular/salverform/bilabiate/papilionaceous/spurred architectures
```

Colour is evaluated as an independent secondary domain and is not required to follow a single universal direction.

## Analysis ladder

### Stage 0 — observation boundary

Use the full frozen physical-island universe in the observation-selection audit.

Rules:

- no-record islands remain in the observation layer;
- no-record island != biologically empty island;
- missing trait != trait state 0;
- collapse GBIF multiplicity to unique `island x species` incidence for flora composition;
- fit outcome-blind selection from region, area, distance and climate;
- carry predeclared IPW / overlap-trimmed and region-specific selection sensitivities.

Observation models cannot manufacture ecological trait values for unobserved islands.

### Stage 1 — floristic-status pathway

Test geography against endemic share among resolved native flora:

```text
endemic_share
  ~ geography + log_area + climate_pc1..4
```

This quantifies the part of island restructuring expressed as floristic/endemic turnover before residual trait filtering is considered.

### Stage 2 — status-stratified genus-fixed residuals

Construct `residual_trait(i,t)` separately by floristic stratum.

Baseline residual model:

```text
residual_trait(i,t)
  ~ geography + log_area + climate_pc1..4
```

with frozen spatial-block cluster-robust uncertainty and direct-evidence support rules.

### Stage 3 — H1: universal-syndrome falsification

Question:

> Does one coherent residual floral/reproductive syndrome occur regardless of regional channel context?

Test pooled geography slopes and their heterogeneity after status and genus composition are represented.

A universal syndrome requires a coherent multivariate direction across biogeographic regimes. Failure is informative and motivates H2; it is not a failed study.

### Stage 4 — H2a: source-exposure necessary-condition test

This is the first direct test of the PR136 wording `the regional flora could otherwise use`.

For a focal channel, define **source-channel exposure** independently of island trait outcomes.

Primary interaction:

```text
residual_trait_vector
  ~ geography * source_channel_exposure
    + log_area + climate_pc1..4
```

Formal prediction:

- in `source_channel_exposure = 1`, isolation may produce the prespecified residual syndrome direction;
- in `source_channel_exposure = 0`, absence of that channel must not by itself generate the same coherent channel-specific syndrome;
- compare the vectors directly. `significant in one group / nonsignificant in another` is not sufficient.

This test does not require calling an island deficient yet. It tests the necessary source-exposure condition.

### Stage 5 — H2b: realized channel-disruption test within source-exposed systems

Restrict to `source_channel_exposure = 1`.

Estimate a channel state from direct or effort-aware channel evidence:

```text
expected_channel_availability - observed_channel_availability
```

The measurement hierarchy is:

1. direct functional interaction / effectiveness evidence where available;
2. occurrence or abundance with observation effort and adequate non-detection;
3. climatic compatibility only as a supporting opportunity diagnostic, never as deficit itself.

Test:

```text
channel_deficit
  ~ geography + climate + area/source-connectivity + observation_effort
```

then:

```text
residual_trait_vector
  ~ channel_deficit
    + geography
    + log_area
    + climate_pc1..4
```

Support for the PR136 channel mechanism requires the residual syndrome to track measured deficit beyond the total geography gradient.

### Stage 6 — H2c: retained versus disrupted versus structurally absent contrast

Compare the same predeclared residual syndrome contrasts among:

```text
retained
expected-but-deficient / disrupted
structurally absent
```

The killer prediction is:

> `disrupted` should differ from both `retained` and `structurally absent` in the predicted multivariate direction.

If `structurally absent` and `disrupted` show the same syndrome simply because the pollinator is absent, the channel-gated interpretation is weakened.

### Stage 7 — H3: effective replacement attenuation

This is the strongest version of the theory and requires direct alternative-channel data.

Define effective replacement from service rather than presence alone, ideally:

```text
replacement = visitation x per_visit_effectiveness
```

Test among source-exposed systems:

```text
residual_trait_vector
  ~ channel_deficit
    * effective_replacement
    + geography + area + climate
```

Prediction:

```text
high deficit + low replacement  -> strongest syndrome
high deficit + high replacement -> attenuated / redirected syndrome
```

If replacement data are not adequate, this stage remains explicitly untested and must not be inferred from bird/butterfly/moth-like floral traits.

## Multivariate syndrome gate

The term `syndrome` is promoted only after two conditions are met.

### Gate A — nonzero residual vector

A joint response-vector test rejects:

```text
H0: all supported residual trait slopes/contrasts = 0
```

### Gate B — prespecified directional coherence

The residual vector must also align with the frozen PR136 contrast across **at least architecture and reproductive-assurance domains**.

Recommended implementation:

- first test the full vector with a cluster-robust joint Wald / equivalent multivariate test;
- then test predeclared domain contrasts;
- require architecture and reproductive-assurance contrasts to point in the predicted direction;
- colour remains secondary;
- do not name a syndrome from a single surviving atomic category.

## Bombus as the first natural experiment

Bombus / large long-tongued bee function is the first channel case because PR136 already distinguishes source applicability from structural absence and has preserved climatic-compatibility diagnostics.

### Confirmatory hierarchy

1. **source-exposed versus structurally absent** is defined from biogeography/source evidence, outcome-blind;
2. within exposed systems, classify `retained` versus `expected-but-deficient` using occurrence/effort evidence where defensible;
3. climatic compatibility may condition whether establishment is environmentally plausible but cannot define deficit;
4. test the residual multivariate syndrome contrast after status and genus controls;
5. direct alternative-guild replacement is a stronger follow-up, not imputed from floral phenotype.

## Falsification table

| Result | Interpretation |
|---|---|
| same coherent syndrome in all regimes | supports a more universal island filter; weakens channel-gated necessity |
| syndrome only where source channel was available and disrupted | supports PR136 core prediction |
| structurally absent regime shows same syndrome as disrupted regime | absence alone or another common filter may suffice; weakens focal-channel explanation |
| source-exposed disruption occurs but residual traits remain null after genus control | channel loss may affect reproduction/demography without producing assemblage trait filtering |
| raw trait signal disappears after status/genus control | floristic/lineage turnover is sufficient for that apparent syndrome |
| deficit effect disappears after observation/IPW/coverage sensitivity | measurement/selection dependence; do not promote |
| deficit signal weakens under high measured replacement | supports the full channel-gated/replacement extension |
| SC-distance signal remains but channel-deficit signal does not | Baker/colonization assurance is sufficient; focal channel not required |

## Observation and support rules

- all 8,265 candidate islands stay in the observation universe;
- ecology models use only directly supported observed flora for the relevant response;
- no minimum mainland-distance cutoff;
- geography is continuous and can encode source-pool accessibility;
- direct trait evidence is preferred for confirmatory inference;
- `<30` islands per atomic response: not promoted;
- `30–49`: pilot;
- `>=50`: confirmatory count support;
- BH FDR is applied within declared stratum x region/channel-regime x trait-domain families;
- leave-one-spatial-block, information-weight and coverage sensitivities are mandatory for promoted results.

## Relationship of current post-PR136 analyses to this design

Recent `trait_probability`, raw `when_where` response-vector and all-island observation analyses are retained as **pattern/observation diagnostics**. They do not by themselves test the PR136 core mechanism because they do not define the primary response as a status- and genus-fixed residual channel contrast.

They become part of the primary chain only after their response is replaced by or explicitly linked to the PR136 residual trait construction.

## Claim ladder

### Chapter 1 can claim

- where geography is associated with floristic-status turnover;
- where status/genus-fixed residual trait filtering remains;
- whether a coherent residual syndrome is contingent on source-channel exposure;
- whether direct channel deficit, if adequately measured, tracks that residual syndrome.

### Chapter 1 cannot claim without stronger data

- historical pollinator loss from geography alone;
- functional replacement from alternative-guild presence alone;
- causal floral evolution;
- assembly versus in-situ evolution;
- pollinator identity from one colour/form category.

## Decision rule for the Chapter 1 headline

The PR136 hypothesis becomes the Chapter 1 headline only if the following chain survives:

```text
observation boundary
  -> status pathway quantified
  -> genus-fixed residual vector remains
  -> source-exposure contingency is supported
  -> disrupted differs from retained and structurally absent
  -> observation / coverage / spatial sensitivities do not reverse the result
```

If the chain breaks before source exposure or disruption, the correct Chapter 1 conclusion is a boundary result: geography restructures floristic status and lineage composition, but the proposed channel-specific floral syndrome is not recovered.

If the chain survives through direct deficit and replacement attenuation, the broader `channel-gated island assembly` mechanism is strongly supported and can be handed to Chapter 2 for plant-level functional identification.
