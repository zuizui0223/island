# Source-archive deletion and exact-plant-species family holdout: 928 island floral assemblages

**12 November 2026 | New original-data analysis, completed; exploratory and source-provenance-limited.**
Source: [successful island Actions #43204982054](https://github.com/zuizui0223/island/actions/runs/43204982054), linked to [island PR #323](https://github.com/zuizui0223/island/pull/323). Original GloBI v0.9 `interactions.tsv.gz`, DOI `10.5281/zenodo.20546682`, SHA256 `095e57398420d1e0746e61233c796beedcfb8e03bd0c6a732aa17ae62043bd17`, was verified. Original GSHHG source SHA256 `8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41`. Source identifiers, exact georeferences, species, polygons, bibliography keys and provider identities stay in the restricted source job, not this public report.

## Scientific decision

The strongest independent **fd** ecological finding remains the directly observed, nonrandom **colour × flower-form × floral-symmetry coordination**, and the modest cross-family **structural morphology × documented global visitor-composition association** established in earlier genus-level models. **Neither this nor the new original GloBI locality overlay demonstrates pollinator shift, effective pollination, adaptation or phenotype evolution.**

This new analysis closes two prior evidence gaps:

1. In the original full 928-island geographic cohort, the same-island *exact accepted plant species* bibliographic records become much sparser after removing the single archive-first source label selected on the entire global plant sample, *before* inspecting locality outcomes.
2. An independent **exact plant-species** (not plant-genus aggregate) five-channel documentation-share prediction, holding out **entire botanical families**, gives no reliable incremental gain from flower colour/form/symmetry pairwise interactions, before or after removal of the dominant original archive label.

## Original full macro panel and source omission

The 928 eligible floral island lists are frozen before reading visitor outcomes. Of these, **927** have original WKB-identical GSHHG polygons and **one** is unassessable (never treated as a visitor absence). The Direct strict high/medium, single-state, triple-trait cohort has **2,135 accepted plant species** and **81,593 evaluable island × plant species list records**. Exact accepted species, original island polygon, original GloBI citation and original source archive label must agree for a documented local relation. Original v0.9 lacks event-time data.

| Observational unit / source-verified property | All cited sources | After removing largest global archive-first label |
| --- | ---: | ---: |
| Islands with ≥1 same-flora-plant geolocated citation | **109** | **28** |
| Distinct island × same-plant species pairs with citation | **552** | **166** |
| Strictly native-plant island × species pairs | **157** | **84** |
| Islands with strictly native-plant documented locality | **40** | **14** |
| Island × plant pairs with ≥2 **documented** visitor guilds | **213** | **96** |
| Strict native pairs with ≥2 documented guilds | **72** | **50** |

Original georeferenced relationship units matching exact island and species: **3,182**; distinct cited bibliography keys: **2,155**; distinct original archive-first labels in those local relation records: **12**. Only **47 of 552** island×plant source-cited pairs had ≥2 archive-first labels even before deletion (archive-label multiplicity is not study-level independence). There are **265** locality-cited pairs with **unresolved plant native status**; do not promote these to native. Bird + at least one other cited guild occurs in **10** island×plant pairs of any status, **5** strict-native pairs; these are bibliography/locality coincidences, not contemporaneous effective-visitor mixtures.

The dominant archive-first label contributes **84.7013%** of original globally cited species × guild × archive × bibliography units in the entire *pre-selected 2,135-plant* panel; it accounts for **65.7448%** of matched local cited units. This label may be an umbrella archive aggregating human studies. Omission is **source sensitivity**, not leave-one-independent-study-out and not a missing-visit experiment.

Documented original same-island guild pair counts (all origins → after deletion): Bombus **209→110**; other bees **251→95**; Lepidoptera (moths and butterflies combined) **222→32**; flower-visiting birds **25→6**; Diptera **190→70**. These overlap across guilds and must not be added to infer distinct plants.

## Separate exact-species whole-family-heldout floral prediction

Use ONLY original source-backed species-direct colour, floral form and symmetry, and original exact plant-species GloBI **global** documented five-guild bibliography-reference shares. Train a source-documentation-only model versus increasingly elaborate morphology features with **one botanical family at a time entirely excluded from training**. The same **243 species / 45 unambiguously named families** have ≥3 distinct cited reference keys *both before and after* deleting the globally dominant archive. The original global candidate comprises 2,135 floral species, of which 2,128 have unambiguous botanical-family assignment; small selected test support limits transportability. No local same-island visitation values enter the prediction target.

| Morphology added to source-documentation baseline | Before: species-equal MSE gain | After: species-equal MSE gain | Before: family-equal gain [conditional 95%] | After: family-equal gain [conditional 95%] |
| --- | ---: | ---: | ---: | ---: |
| Open flower form | **+0.000255** | **−0.000147** | +0.000479 [−0.000652,+0.001657] | +0.000206 [−0.001056,+0.001445] |
| Open form + symmetry | −0.000555 | −0.000865 | +0.000137 [−0.001055,+0.001348] | +0.000086 [−0.001238,+0.001368] |
| Open form + symmetry + inconspicuous colour | −0.000657 | −0.000736 | +0.000123 [−0.001047,+0.001269] | +0.000338 [−0.000968,+0.001628] |
| Three morphology terms + **all three pairwise interactions** | **−0.000958** | **−0.000928** | −0.000363 [−0.001880,+0.001087] | **+0.000005** [−0.001558,+0.001570] |

The baseline five-channel species-equal MSE is **0.027309 before** and **0.024267 after** (the *target documentation composition differs* after source omission; do not interpret their raw MSE difference as a longitudinal improvement). Positive gain means reduced heldout prediction MSE. The species-equal and family-equal estimands are distinct; **all family-equal 95% intervals above cross zero**. The 115-feature genus-morphology model and this exact-species three-term interaction model are different experiments and cannot be pooled as independent negative replications.

**Consequence:** simple directly recorded **floral architecture** can be a limited, source-sensitive *global documented visitor-composition proxy*, but colour×architecture coupling does **not** acquire a transferable, independently supported **pollination syndrome mechanism** through this new exact-species test. Earlier GloBI genus-level *documented bird share* results that survived the same broad source-label omission on 344 genera/77 families remain a separate narrow posthoc result; they cannot be promoted to island bird pollination, real bee→bird shifts or an all-channel species-level claim.

## Claim boundary and interpretation

- The original 928-island direct floral **within-island trait-coupling** result (909/928 positive under same-marginal colour shuffle) remains a **phenotypic coordination** finding. This test did not depend on GloBI source labels.
- Only **109/927** evaluable floral islands have any exact-flora-plant GloBI cited locality evidence and only **40/927** have strict-native plant evidence, falling to **28** and **14** islands after the largest archive is omitted. Never interpret the other islands as guild absence.
- Recorded `pollination` predicates do not quantify pollen deposition or successful seed set. GloBI v0.9 does not supply event timestamps; moth and butterfly relationships are pooled into Lepidoptera.
- A single accepted species has the same fixed floral trait label at every location in this dataset; an island–mainland floral difference is **species sorting by construction**, not local flower-colour/shape evolution.
- The family-heldout comparison is posthoc conditional on a selected documentation-rich 243-species set and can still carry source/provider, taxonomic and ecological ascertainment bias. A 45-family conditional interval is not a phylogenetic or causal evolutionary test.

## Reproducibility and status

The original runner `scripts/run_fd_globi_provider_v09_audit.py` emitted nested `macro_928_Direct_flora_exact_species_GloBI.original_exact_island_locality_coverage` and `macro_928_exact_species_flower_syndrome_documentation_transfer` inside the **successful source job #43204982054**, with completed original-input and privacy/synthetic safeguards. Its top-level legacy status `PROVIDER_COVERAGE_ONLY_NOT_FAMILY_PREDICTION` describes the original coverage payload and should **not** hide the nested model's explicit `COMPLETED_EXPLORATORY_EXACT_SPECIES_FAMILY_LOFO_DOCUMENTARY_SHARES` status. Treat this as a new descriptive falsification/evidence-gating result. Both `fd` PR #14 and `island` PR #323 remain **Draft**, not an independently validated evolutionary manuscript result.
