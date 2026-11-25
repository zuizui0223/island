# GloBI five documented visitor channels: 88-family disaggregation and within-family null

**2026-11-11 · fd Draft PR #14 · exploratory original-data validation.** The response is the relative frequency of five GloBI *global interaction reference* categories attached to plant genera, NOT local island pollinator fauna, effective pollen transfer, or pollination syndromes.

## Original inputs and design

- Protected original `island` species Direct high/medium single-state flower trait archive, SHA256 `b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`.
- Frozen GloBI source-genus interaction-breadth artifact 10328696801, SHA256 `475dda8b0eef4ce22ad1aa736fa82cc049042d801217c65d52641489288ca254`; original source-of-truth is GloBI v0.9, published June 4, 2026, `10.5281/zenodo.20546682`. GloBI itself states its recorded interaction sources are geographically, temporally and taxonomically biased.
- Original frozen GBIF genus–family source-backed crosswalk, **484 genus IDs** uniquely assigned to **88 families** after removing 2 competing historical taxonomic labels, with **1,590** directly measured morphology species.
- Each botanical family is held out **in its entirety** when predicting documented genus×visitor reference proportions. Baseline is 3 training-only normalized documentation-effort covariates: `log1p(GloBI independent refs)`, `log1p(matched GloBI plant names)`, `log1p(direct morphological species in genus)`. Test adds **only fraction of open-radial flower form** in a genus. Fixed ridge alpha100 and Euclidean projection to five nonnegative shares summing to 1. This follows the validated 88-family leave-one-out result, without using a true local island pollinator response.

## Five guilds: genus- and family-equal gain are different estimands

Positive proper five-channel composition MSE gain means adding open floral form reduces prediction error compared with documentation effort alone.

| Documented global GloBI channel | Genus-equal MSE gain | Family-equal MSE gain | Family-equal conditional bootstrap 95% interval | Families positive |
| --- | ---: | ---: | --- | ---: |
| Bombus | +0.000272 | **−0.000265** | −0.001749 to +0.001172 | 42/88 |
| Other bees | +0.001348 | +0.001717 | −0.000947 to +0.004345 | 52/88 |
| Lepidoptera (butterflies AND moths) | +0.001054 | +0.000164 | −0.002718 to +0.002920 | 47/88 |
| **Flower-visiting birds** | +0.000726 | **+0.002091** | **+0.000927 to +0.003533** | 53/88 |
| Diptera | +0.002334 | +0.001741 | −0.000362 to +0.004009 | 55/88 |
| **Five-channel MSE average** | **+0.001147** | **+0.001090** | **+0.000122 to +0.002077** | 54/88 |

The conditional 4,999-draw *family-equal* bootstrap assigns each of the 88 family mean scores one vote; it does NOT imply 88 truly independent interaction-study providers or phylogenetic histories. Only the **documented bird-reference share** has an individual family-equal interval entirely positive. Relative visitor reference shares are compositionally dependent, so the pooled five-channel MSE is the more stable global metric.

The two largest families represent **107/484** genera. Scoring after removing them (without retraining) retains a positive average gain across each guild; the five-channel mean rises to +0.001571 over 377 genera. This is a support sensitivity rather than source-independent generalization.

## Within-family morphology-label negative control

To check whether the result only reflects each family's fixed average morphological composition, for every family with two or more genera and a genuinely varying open-form fraction, shuffle **only the genus-level open-form feature within family** while preserving its original family membership, GloBI response and three documentation covariates.

- **32/88 families / 308/484 genera** permit a nontrivial label exchange. The remaining families are single-genus or have no observable within-family open-form variation; the null therefore has limited discriminatory power.
- Repeat 999 shuffles, refitting all 88 family-heldout ridge models for each. This null preserves each family's recorded form-rate distribution and all source-effort features but is **NOT a fine-phylogeny or original data-provider null**.
- In the five-channel aggregate, observed genus-equal gain **+0.001147**, shuffled median **+0.000432**, shuffled 95% **[+0.000201,+0.000751]**, exploratory upper-tail fraction **0.001**.
- For documented **bird** interaction references, observed +0.000726, shuffled 95% [+0.000057,+0.000546], exploratory upper-tail **0.001**.
- For documented **Diptera**, observed +0.002334, shuffled 95% [+0.000395,+0.001561], exploratory upper-tail **0.001**.
- The other three channels do NOT all pass the individual family-equal uncertainty test. Source-dependent or compositional differences must not be exaggerated into a five-guild mechanism.

The 6 pooled/guild observed gains and all 999-shuffle null summaries were **recomputed with an independent second analysis function** against the same source inputs and matched numerically at machine precision (<1e−12). The public source/test code runs synthetic inputs only; this numerical replay is not independent ecological replication.

## Provider-independent interpretation remains blocked

The frozen 2,228-genus GloBI summary contains reference counts per pollinator channel but **NO original `provider_id` or `reference_id` on the plant genus×channel units**. The separately acquired ~6,928-row pollinator-taxon catalog has `reference_examples` but these are **pollinator-side species aggregates**, not a verified join to the scored 484 plant genera. They cannot reconstruct source-provider participation for the outcome being predicted. This is not a reason to invent provider labels.

GloBI's 2026-06-04 v0.9 Zenodo interpreted `interactions.tsv.gz` (SHA256 `095e57398420d1e0746e61233c796beedcfb8e03bd0c6a732aa17ae62043bd17`) includes original dataset citation and reference information and is about 2.7GB compressed. A genuine provider-disjoint test requires streaming that exact version, mapping accepted plant genus and five evidence types, removing duplicated references across providers and recomputing GloBI documented shares **after provider exclusion**. Existing public `src/fd/globi_provider_provenance_gate.py` explicitly rejects the current aggregated data as insufficient. No provider-heldout biological result is claimed.

## Claim ceiling

**Supported narrowly:** directly recorded coarse open-radial form structure is associated with a small out-of-family predictive improvement in the *global documentation shares* for some GloBI visitor groups; only bird references have a clearly positive family-equal interval among the five individual channels in this exploratory audit.

**Not supported:** the abundance, presence or absence of true bird/fly/bee/moth pollinators on an island, a functional pollination syndrome, effective flower visitation or pollen transfer, adaptive trait evolution or independent-source generality. Earlier exact-island five-channel GBIF record qualification lacks enough true nondetections for a balanced binary local guild turnover test, and Lepidoptera remains inseparable into moths and butterflies.

### Released artifacts

- `src/fd/globi_family_guild_audit.py`: leak-checked 88-family model, within-family negative control, original-record rights-safe group summaries.
- `tests/test_globi_family_guild_audit.py`: 3 synthetic data tests for correct family holdouts, missing/ambiguous family/morphology, unnormalized five-channel targets and absence of private identifiers.
- `results/20261111_floral_syndrome_coupling/per_guild_88_family_within_family_morphology_null.csv`: only 6 anonymous aggregate rows; no protected taxa, genus/family/source names, pollinator species or individual original citations.

**Next genuine scientific gate:** provider-heldout GloBI v0.9 target reconstruction and truly independent local pollinator reference/effectiveness validation, NOT more partial-absence claims from unresolved GBIF channels.
