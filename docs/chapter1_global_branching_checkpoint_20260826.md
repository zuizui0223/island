# Chapter 1 global branching checkpoint

Date: 2026-08-26
Central hypothesis: **Dispersal-filtered pollination-service branching hypothesis**

## Result first

Chapter 1 now tests a pollinator-name-free plant response vector before consulting any
guild-labelled syndrome. The global result supports **regional branching**, not one
universal floral island syndrome and not a simple northern-selfing versus tropical-mobile-
pollinator dichotomy.

For all native species, the two-dimensional response vector differs directly between
northern mid-latitude and tropical islands (`q = 3.20e-6`):

- northern mid-latitude islands show increasing generalized accessibility with distance
  (`slope = +0.04755`, `q = 0.0282`), without supported strict reproductive assurance
  (`slope = +0.02238`, `q = 0.348`);
- tropical islands show increasing strict reproductive assurance
  (`slope = +0.10015`, `q = 0.0282`), without an axis-resolved accessibility change
  (`slope = -0.04314`, `q = 0.125`).

The native-nonendemic result is even less compatible with the simple dichotomy. Tropical
islands show both increasing reproductive assurance (`+0.13848`, `q = 1.44e-4`) and
increasing specialized-access architecture (`-0.05809` on generalized accessibility,
`q = 0.0376`). The north-versus-tropical vector difference remains supported
(`q = 3.63e-7`).

The supported conclusion is therefore:

> Mainland distance is associated with different combinations of reproductive assurance
> and floral accessibility in different biogeographic settings. Isolation does not map to
> one universal simplification or selfing response.

This establishes the **when/where plant-side pattern**. It does not identify effective
pollination service, pollinator loss, attraction relaxation, competition for scarce
pollinators, or functional replacement.

## Central causal concept

```text
pollination channels available in the source
x island distance and connectivity
x pollinator dispersal and establishment ability
    -> effective pollination service remaining on the island
       |-- reproductive-assurance branch
       |     -> SC, autonomous selfing, increased selfing
       `-- attraction/access branch
             |-- return to attraction collapses -> plainness/generalized access
             |-- competition for scarce service -> maintained or reinforced attraction
             `-- another guild supplies effective service -> maintenance or redirection
```

Chapter 1 observes only the right-hand plant responses and the geographic gradient. The
three upstream terms and the effective-service mediator remain latent. Consequently the
diagram is the theory being tested for compatibility, not an identified causal path.

## Predeclared evidence hierarchy

### Primary: universal plant response

The primary vector contains no pollinator name:

1. `accessibility_generalization`: the continuous `generalized_accessible` trait-
   concordance score;
2. `reproductive_assurance`: the strict `selfing_core` score, using SC, mating system and
   autonomous selfing but excluding floral-attraction traits.

The multivariate context vector is tested first. An axis is labelled only when both the
context-vector gate and its own FDR-adjusted test pass.

### Secondary: non-exhaustive guild concordance

Large-bee-like, butterfly-like and bird-like scores are placed in a separate secondary
test family. They are examples available from the current trait templates, not an
exhaustive channel set. Moths, flies, bats, beetles, wasps and other possible channels are
not represented; an unmodelled guild is never coded as an absent guild.

These scores measure multivariate floral-trait concordance. They do not observe visitor
identity, abundance, visitation, pollen delivery, seed production, mobility or
establishment.

## Where the branches are currently detected

### Broad analysis regimes, all native species

| Context | Testable islands / clusters | Multivariate vector | Axis-resolved branch |
|---|---:|---:|---|
| northern mid-latitude | 237 / 36 | supported, `q = 0.0380` | accessibility generalization |
| tropical | 134 / 45 | supported, `q = 0.00134` | reproductive assurance |
| northern high-latitude | below support | not testable | not testable |
| southern extratropical | below 50-island support | not testable | pilot vector unsupported |

The northern-midlatitude and tropical vectors differ directly (`q = 3.20e-6`). This is
the main global when/where result.

### Formal biogeographic realms

Palearctic all-native composition shows both increasing accessibility generalization
(`+0.10860`, `q = 5.57e-7`) and reproductive assurance (`+0.07843`, `q = 3.0e-6`).
Neotropical all-native composition also shows both (`+0.03339`, `q = 0.0337`; `+0.11077`,
`q = 0.0161`). Their all-native two-axis vectors do not differ at the predeclared threshold
(`q = 0.0574`).

In native non-endemics, Palearctic again shows both axes, whereas Neotropical retains only
reproductive assurance. The direct Palearctic-versus-Neotropical vector difference is then
supported (`q = 0.00170`). This status dependence is evidence that floristic composition
matters; it is not evidence that one pollinator guild is absent.

Nearctic remains below the 50-island confirmatory support gate. Its earlier 31-island
pilot non-replication is retained as a pilot result, not promoted to a global exception.

## What the sampled guild concordances say

The secondary scores are useful mainly as compatibility and specificity checks:

- northern mid-latitude large-bee-like concordance declines (`-0.04927`, `q = 0.0320`),
  which is compatible with the northern large-bee-decline side hypothesis;
- northern bird-like and butterfly-like slopes do not survive the secondary axis family;
- tropical bird-like, butterfly-like **and large-bee-like** concordances all increase;
- Palearctic bird-like, butterfly-like and large-bee-like concordances all decline;
- Neotropical has no supported sampled-guild axis;
- the southern-extratropical pilot increases bird-like and large-bee-like concordance
  together, rather than isolating a bird-specific signal.

Because several templates share colour, size, tube and form traits, the simultaneous
tropical increase in all three demonstrates that these are not pollinator-identity
classifiers. The current results do not support the claim that birds or butterflies
functionally replaced lost pollination service. The same logic applies to Bombus: the
northern/Palearctic pattern is compatible with the side hypothesis, but does not show that
Bombus loss caused it.

## Analysis universe and apparent attrition

The frozen geography universe still contains exactly **8,265 islands**. It was not cut in
half and unobserved islands were not coded as zero-trait islands.

- geography/isolation universe: 8,265 islands;
- any current island-level syndrome score: 424 islands;
- any derived branch score: 423 islands;
- syndrome artifact: 5,299 rows, which are repeated island x syndrome x floristic-stratum
  outcomes, not 5,299 different islands;
- derived branch artifact: 4,412 repeated island x axis x stratum rows.

The fitted sample is therefore an observation-support subset nested inside the retained
8,265-island universe. Coverage and selection sensitivity remain mandatory; missing trait
support is neither biological absence nor a zero score.

## Chapter boundary

### Chapter 1: global when/where

Identifies:

- whether the plant response is universal or regionally branched;
- where accessibility and reproductive-assurance responses occur along the distance
  gradient;
- whether region vectors differ directly;
- whether guild-labelled trait concordances are compatible with, or contradict, proposed
  regional side hypotheses.

Does not identify:

- source-channel exposure;
- a pure isolation effect separate from source supply and connectivity;
- pollinator mobility or establishment;
- effective service (`visitation x per-visit effectiveness`);
- competition, attraction elasticity, or functional replacement.

### Chapter 2: direct mechanism in `izu-core`

The direct test belongs in <https://github.com/zuizui0223/izu-core>:

```text
source-channel evidence x connectivity x mobility/establishment
    -> visitation and per-visit effectiveness
    -> reproductive assurance and quantitative attraction/access traits
```

The Izu system can test whether distance-related loss of functional matching and effective
service predicts the plant-side branches observed globally. It is the mechanism chapter;
Chapter 1 remains the global pattern and boundary chapter.

## Frozen inputs and reproducibility

- distance/regime run `29228212586`, artifact
  `purpose-shortest-distance-regime-29228212586`;
- syndrome/pathway run `32922758001`, artifact
  `pr138-attraction-selfing-32922758001`, digest
  `sha256:5c573fdd75ebab4b85dd97cd11a064d9764eb3772c8fa37e124236eb841e6770`;
- realm run `32927060210`, artifact `pr138-realm-sensitivity-32927060210`, digest
  `sha256:6556a7938a52547ddf7b630a822d743b2ba5ca479589716919441b96dcfa0521`.

The hypothesis contract is `config/chapter1_global_branching.yml`. The executable analysis
is `src/island_v2/chapter1_global_branching.py`; it writes branch scores, within-context
omnibus tests, direct between-context tests, classified primary branch states and a
fail-closed manifest.
