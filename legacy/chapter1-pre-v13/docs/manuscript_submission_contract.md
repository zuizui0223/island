# Chapter 1 manuscript submission contract

## Status

The Chapter 1 analysis is frozen as a **submission candidate** under
[`config/chapter1_submission_freeze_lock.json`](../config/chapter1_submission_freeze_lock.json).
The lock points to workflow run `33067239884`, artifact
`chapter1-area-capacity-moderation-33067239884`, and head
`6c2d7e267f872dc24983969cb6475ca2a8051975`.

The artifact contains the nested PR #138 pattern/robustness inputs, PR #139
source-matched lineage-representation outputs, and the continuous area-moderation
extension. The analysis is frozen, but the current GitHub Actions copy expires on
2026-11-25. A durable release/deposit and journal-formatted manuscript remain required
before submission.

PR #140 and later trait-acquisition products are post-freeze. They do not enter this
submission candidate unless a prospective update is declared before inspecting revised
Chapter 1 outcomes.

## Canonical question

> **When and where do floral and reproductive assemblage responses emerge along a
> mainland-distance/source-accessibility gradient, and how are those responses divided
> between source-pool lineage assembly and plant-side reproductive/accessibility branches?**

This is a global comparative analysis of assemblage composition. It does not identify a
historical pollinator-loss event, an effective pollination service, or in-situ floral
evolution from cross-sectional island data.

## Theory tested at the observable boundary

```text
mainland distance / source accessibility x island area / capacity
|
|-- plant source-pool filter
|     |-- source-available lineage entry
|     |-- within-lineage species loading
|     `-- Baker-compatible reproductive assurance
|
`-- pollination-channel filter
      `-- regional attraction/accessibility response
```

Chapter 1 observes geography, area, floristic composition, lineage representation, and
plant traits. Pollinator arrival, establishment, persistence, visitation, per-visit
effectiveness, and population maintenance remain latent. The second branch is therefore
a compatibility interpretation, not an identified causal pathway.

## Canonical analysis hierarchy

### 1. Global WHEN / WHERE pattern

Atomic floral/reproductive states are combined in grouped-binomial response-vector tests.
The primary within-context null is that every supported geographic slope in the vector is
zero. The direct between-context null is equality of the two geographic response vectors.
Responses require at least 50 islands for confirmatory inference; 30-49 islands are pilot
only.

The frozen result detects nonzero vectors in northern mid-latitude and tropical floras,
with a direct difference in both `all_native` and `native_nonendemic` strata. Northern
high-latitude and southern-extratropical floras remain confirmatorily unresolved.

### 2. Pollinator-name-free branch synthesis

Two prespecified plant responses summarize the higher-level branching test:

- `accessibility_generalization`;
- strict `reproductive_assurance`, excluding attraction traits.

The direct northern-midlatitude versus tropical difference remains supported across all
24 source-mode x outcome-blind selection-mode x primary-stratum scenarios. The broad
northern branch is not independently IPW-robust; the clearest positive branch is
Palearctic. This distinction prevents a latitude proxy from being relabelled as a
universal northern or Bombus effect.

### 3. Observation, support, and spatial safeguards

- all 8,265 islands remain in the observation-process universe;
- raw records are collapsed to unique island x species incidence;
- missing observation support is never coded as trait absence;
- information-weight, equal-island, trait-resolution, distance-form, and outcome-blind
  selection sensitivities are retained;
- spatial-block cluster-robust inference and leave-one-block-out influence tests are
  retained;
- status-stratified persistence is descriptive of assemblage strata, not a causal test of
  endemism.

### 4. Source pool and lineage decomposition

GIFT mainland membership and frozen outcome-blind source assignments define candidate
source availability. Source-availability matching holds the prevalence and source-richness
class structure fixed.

Three components are kept separate:

- `entry_enrichment`: which source-available genera are represented;
- `species_enrichment`: species-weighted lineage functional position;
- `loading_increment`: extra species weighting after genus entry.

Broad Palearctic enrichment is predominantly a genus-entry pattern: entry is positive and
FDR-supported in 8/8 primary source-mode x stratum scenarios, while loading is supported in
0/8. Exact-SI results shift toward species loading, but the final SI assemblage slope does
not exceed the genus-fixed source-pool null.

### 5. Continuous area / capacity moderation

Area is continuous; there is no small/large island threshold. Models include
distance, area, distance x area, climate PCs 1-4, and spatial-block cluster-robust
covariance.

The stable area result is source-pool-side: broad Palearctic genus-entry enrichment is
more strongly associated with distance on smaller islands in 8/8 scenarios. Species
loading has no equally stable area law. Neotropical reproductive assurance is
smaller-island amplified in both primary strata, but uses only 59 response-supported
islands and 14 spatial blocks and lacks a supported direct Palearctic-Neotropical
reproductive-assurance interaction difference.

## Manuscript-level results

The Results main line is:

1. no universal classical floral island syndrome is recovered;
2. northern-midlatitude and tropical floral/reproductive response vectors are both
   detectable and differ directly;
3. the clearest robust positive branch is Palearctic, while Nearctic remains pilot and
   tropical floras follow a different vector;
4. much of the floral-accessibility geography is lineage-associated;
5. source-matched Palearctic lineage representation is mainly a genus-entry response;
6. smaller island area amplifies this genus-entry response in all eight frozen scenarios.

The manuscript must not turn individually significant traits into the primary result.
Atomic states, named syndrome scores, exact-SI restrictions, genus deletion, and template
sensitivities decompose or stress-test the frozen vector result.

## Pollination-syndrome interpretation contract

Pollination-syndrome expectations enter only after the plant-side result is frozen.

Allowed:

- compare the regional plant response with literature-defined expectations for large
  bees, birds, Lepidoptera, hawkmoths, or other functional groups;
- describe concordance, mismatch, and alternative explanations;
- state that dispersal and establishment differences among pollinator guilds are candidate
  mechanisms for Chapter 2 or independent field tests.

Not allowed:

- infer actual pollinator identity from floral phenotype;
- claim Bombus loss caused the Palearctic branch;
- claim bird or Lepidoptera replacement caused a tropical or southern pattern;
- convert climatic compatibility or opportunistic nondetection into historical loss;
- describe a syndrome score as measured visitation, pollen delivery, seed set, mobility,
  establishment, or effective service.

Bombus niche/hypervolume, applicability, occurrence, and channel-score products remain
diagnostic or sensitivity products and are outside the primary Chapter 1 model.

## Claim ceiling

Supported:

> **Within opportunistically observed island floras, isolation/source accessibility is
> associated with biogeographically contingent floral and reproductive assemblage
> composition. The strongest source-matched area-conditioned result is greater broad
> Palearctic genus-entry enrichment with distance on smaller islands.**

Not established:

- pollinator loss, replacement, mobility, persistence, or effective service;
- genetic founder effects, founder number, drift, or a directly observed Baker process;
- historical mainland ancestry or colonisation probability;
- in-situ evolution, differential diversification, or within-lineage causal adaptation;
- a universal area-dependent floral island syndrome.

## Reproducibility and archival gate

The automated freeze validator must verify:

- GitHub artifact metadata and archive SHA-256;
- all nested PR #138/#139 inputs required by the area analysis;
- row counts and unique result keys;
- q-value bounds and hierarchical vector-before-axis classification;
- the 8/8 Palearctic genus-entry and bounded Neotropical results;
- the explicit causal claim ceiling.

Before journal submission, the expiring Actions artifact must be copied to a durable
release or repository with the same digest, and the manuscript must cite that permanent
identifier. Figure/table source files must be generated only from the locked artifact.

## Chapter 1 -> Chapter 2 boundary

Chapter 1 identifies where assemblage branches differ and how source-available lineages are
represented. Chapter 2 (`izu-core`) asks which directly measured interaction states produce
those conditional response signs. Its empirical gate requires measurable functional
position, visitor exposure/effectiveness, and reproductive outcome; Chapter 1 phenotype
scores alone cannot satisfy that gate.

## Noncanonical material

- Bombus-primary M0-M4 variants and coefficient-product mediation;
- climatic suitability or occurrence nondetection treated as realized service;
- PR #140 and later trait acquisitions unless prospectively promoted;
- interpretations based only on counts of individually significant atomic slopes;
- old latitude-only conclusions that do not retain the Palearctic/Nearctic distinction.
