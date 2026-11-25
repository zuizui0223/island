# Frozen test: tropical radial × deep-tube flowers versus mapped pollinator guilds

**Status 2026-11-19:** post-hoc exploratory hypothesis frozen after inspecting the existing +10.99 pp descriptive result; module implemented with synthetic tests; **not executed on original plant×island or pollinator-map data**. Existing GIFT 4,505-island positive results are unchanged. This is an exploratory ecological association, not an observation of visits, pollen transfer, fitness or a historical pollinator shift.

## Primary estimand and single joint test

The existing original tropical actinomorphic × deep-tube, source-relative Direct motif excess is approximately +10.99 percentage points (666 eligible floras), driven largely by between-genus botanical composition. The new, independent test asks whether island-minus-eligible-mainland candidate SOURCE differences in this actual same-species pair motif co-vary with mapped, source-relative nectarivorous-bird and Sphingidae guild availability.

Freeze: all 4,505 original GIFT floras as parent; tropical analysis region; all-recorded botanical status; first five ranked GIFT native mainland candidate source regions ("geo_k5"); strict **same accepted species** Direct actinomorphic and deep-tube sector; require at least five measured same-species pairs both for island and source, plus 60 eligible islands, 15 disjoint 10-degree physical blocks, and 15 distinct first-source groups. No cross-products of separate species-level single-trait margins.

Response: Direct known-pair sector proportion on island minus the corresponding source pool proportion. Unmeasured species are **missing**, not sector-negative. The primary two independent covariates are log(1+mapped qualified nectarivorous-bird range species richness) and log(1+qualified Sphingidae modelled species richness), island-minus-source and standardized within the eligible tropical test cohort. Guild variables come only from independently acquired animal data, without focal flower traits. Neither zero GBIF detections nor generic Aves/Lepidoptera record totals count as qualified absence or moth/nectarivory measures.

Fit one joint 2-df null of no guild association, with both guild terms entered together, controlling log(1+distance), log(1+area), absolute latitude, temperature and precipitation PCs, and separately recorded original Direct island/source pair coverage. Original disjoint 10-degree physical blocks and first GIFT source group each receive equal cluster weight with separate CR1 F approximations and 499 group-resampled intervals. A directional positive *exploratory, post-hoc gate* requires both partial guild slopes nonnegative, joint p<0.05 by both group definitions, and the lower CI for at least one guild greater than zero under both groupings. This test was chosen after seeing the tropical +10.99 pp result and is therefore post-hoc exploratory, not prospectively preregistered. The joint p-values are exploratory, not family-wise confirmatory nor tests of causality.

## Before any original-data outcome: mandatory botanical genus baseline (correction 2026-11-19)

The first synthetic-only version omitted explicit adjustment for the botanical
genus composition, despite the source-verified tropical floral pairing
difference being largely **between-genus**. This is corrected **before any
original guild map is joined or guild outcome is scored**. The one primary
joint 2-df bird × Sphingidae test is unchanged, but it is now conditional on
`genus_reference_expected_delta`, an independently trained expected
island-minus-mainland-source radial × deep-tube fraction based on each flora's
observed genus composition and **different accepted congeneric plant names**.
It must be calculated without using the focal island's own plant-species trait
states, source-specific held-out plant names, or its flower-sector outcome.
The data builder must document its botanical evidence source, taxonomy,
spatial/source name holdouts, and first-source eligibility. Pass a separate
boolean `genus_reference_qualified`; reject missing or unmatched values,
rather than substituting genus richness or a within-island in-sample predictor.
Values must be real finite signed fractions in [-1,1].

This is deliberately a **conditional residual association**: a null after
genus adjustment does *not* falsify guild-associated taxonomic establishment,
which may be mediated by genus composition. A detected residual is also not
proof of actual visits or of independent within-lineage evolution. An
unadjusted contrast may be reported as descriptive only, never as the primary
guild-effect test. This is a pre-original-data correction to an exploratory
protocol, **not** prospective preregistration.

## Public animal-source acquisition status (19 Nov 2026)

The source-only BirdLife-2020 taxonomy file from AVONICHE
(`Avoniche1_Birdlife2020.csv`, Dryad DOI 10.5061/dryad.1zcrjdg56)
is independently discoverable, CC0-licensed, and contains published quantitative
nectar-use data for 10,999 birds. A strict public, eBird-free receipt/parser
is implemented in `fd.avoniche_nectar_preflight`, with five successful
synthetic tests. However, in original GitHub Actions run **44539952332**,
the two official original-file URLs returned **HTTP 403/401**, respectively.
Therefore **zero actual bird weights were acquired** by that run. Do not
describe species weights, spatial guild availability or the flower-guild
model as source-executed. A local copy of the original BirdLife2020 CSV can
be supplied to the same code with `--input-csv`, with SHA256 receipt.

BirdLife bird range polygons are a separate permission-controlled input;
species-level Sphingidae SDMs with validated GBIF sampling effort are also
not prepared. Until both spatial inputs and the name-disjoint genus baseline
are independently audited, output **NOT_IDENTIFIABLE**, not a positive or
negative guild result. No eBird occurrence source is used.

## Source-independent map requirements

**Birds:** use AVONICHE quantitative nectar-use fractions (https://doi.org/10.5061/dryad.1zcrjdg56) as fixed per-species weights, summed over native resident/breeding species with verified range polygons; do not choose a nectarivory threshold from floral outcomes on exact island geometry, not centroid radii. BirdLife polygons require noncommercial research data access/permission; do not publish them or species-specific locality/identity joins.

**Sphingidae:** require species-resolved adult nectar-feeding moth range/SDM surfaces with spatial and temporal GBIF detection-effort checks; sum species-level modeled occupancy probabilities as expected richness. Plain GBIF occurrence absence is not sufficient to assign zero richness; unknown must remain NA. Distinguish adult nectarivory from merely being in Sphingidae.

**Study population:** genuine original GIFT native source floras and original GSHHG polygon join, with independent 10-degree blocks. Map scores must be built blind to focal floral response. No eBird-based data stream is required.

## Recording coverage and limitations

The reported tropical symmetry coverage difference (island ~47%, source/mainland ~16%) is substantial, but is NOT automatically the same as Direct symmetry×depth **joint** coverage in the analysis. Use the genuine observed pair count and full accepted-species denominator separately for each island and source. The runner adjusts for measured coverage but explicitly returns worst-case bounds that permit every missing pair to be sector-positive or sector-negative. If zero lies in these bounds, a positive guild coefficient cannot establish robust direction of the total sector difference under arbitrary MNAR. A separate source-blind genus-conditioned sensitivity, phylogenetic and endemism checks are needed before interpreting guild-associated botanical sorting as lineage-level adaptation.

## Executable status and privacy

The complete 4,505-island row-count gate alone does not authenticate species/island/source identities: verify SHA-pinned original inputs separately upstream. The ready-to-run module accepts a **private** one-row-per-original-island CSV with exactly 4,505 rows and the columns enforced by the code. The private join must be prepared inside authorized source CI using original plant and separate independently reviewed guild maps. Unknown and non-qualified guild states must not be converted to zero. Execute:

    python -m fd.guild_sector_alignment --private-island-join /secure/original_4505_guild_pair_join.csv --bootstrap 499 --output /secure/anonymous_guild_pair.json

Only anonymous aggregate JSON, an input source manifest and code/test receipts may be published after real execution. A non-estimable support status is a valid result. No protected island, plant or animal identities should be committed to the public fd repository. Do not merge this mechanism into Main until original map and GIFT source hashes, model support, uncertainty and missingness bounds have all been checked.
