# Geological origin, establishment filtering and floral adaptation — scope and inference contract

Status: 2026-11-14 exploratory original-source macroecology. Three original datasets are not interchangeable: independently classified geological origin (Roeble et al. / Zizka), GIFT native-region potential source lists, and original species-direct floral trait evidence. User goal: distinguish establishment and adaptation of floral functional traits between **geologically oceanic and continental islands**.

## Geological axis and independence

- Never infer "oceanic" from LGM mainland disconnection. Geological origin and LGM connectivity are separate axes. Original conservative GSHHG crosswalk yields 176 geologically continental and 258 oceanic labels, plus 8 mixed and many unresolved island IDs; origin source is a conservative centroid/area crosswalk, not confirmed tectonic island population sampling.
- The stricter 928-island source-direct floral phenotype intersection gives **148 origin-classified, source-comparable complete floras**, with 73 northern-midlatitude (50 continental, 23 oceanic) and 64 tropical (21 continental, 43 oceanic), and smaller unsupported regions.
- Prior broad 37 origin-matched pair results share only 5 spatial/source components; later disjoint study selected just 10 origin pairs and no FDR-supported result. Current 148-flora *separate* matching uses area, mainland distance, species number, source genus count, four climate PCs, and forbids any geographic 10° block or verified source archipelago across separate selected pairs. Only **3 matched pairs each** in the two main regions under all-recorded origin cohorts; for WCVP-compatible 2 northern midlat and 0 tropical, for strict GIFT native 1 northern midlat and 0 tropical. **Do not report statistical significance, biological equivalence or an origin-specific evolutionary mechanism using these tiny disjoint samples**. Full original run [#43573857952](https://github.com/zuizui0223/island/actions/runs/43573857952) succeeded in two numerical seeds; these are repeated analyses of identical source data, not new islands.

## Model A — arrival, establishment and assemblage filtering

Potential source inventory from five geographically closest GIFT mainland native regions (not observed seed rain or proven ancestry). All accepted species occurring on the focal island are removed from donor pools before source-relative comparison. Exactly 30 high/medium species-direct flower colour, form, symmetry plants per flora; finite categorical Rao expectation is analytical, no 17D KDE bandwidth or pseudo continuous hypervolume.

Distinct outcomes:
1. Source-normalized marginal diversity (equal categorical trait Rao Q).
2. Joint actual colour × form × symmetry Rao Q and exact colour–architecture pairing departure.
3. Symmetric source vs island genus overlap (Sørensen, Jaccard) and 10 distinct-genera count-balanced expected overlap; compare with asymmetric overlap denominators separately.
4. Trait-associated genus/family turnover, with separately documented native provenance, source detection, area, isolation, climate and regional composition.

The primary biological interpretation is **present-day floristic sorting relative to an observed inventory**. Absence from a GIFT checklist does not establish colonization failure. Genus overlap differs with source vs island denominator; equalized ten-genus and Jaccard source contrasts currently do not show robust global oceanic-minus-continental separation. A source-normalized effect must survive geographic and taxonomic controls and an independently collected external island cohort before it can support geological source-selective establishment.

## Model B — adaptation and divergence after establishment

Original `chapter1_status_flora.csv.gz` independently contains `origin_status` (native/introduced/unresolved), `endemic_status` (endemic/nonendemic/unresolved), and combined `floristic_status` (endemic/native_nonendemic/native_endemism_unresolved/...). Earlier post-hoc output that classified endemism by looking for `origin_status=='endemic'` is **INVALID** and superseded. Strict potential endemic is `origin_status == 'native' && endemic_status == 'endemic' && floristic_status == 'endemic'`; strict native-nonendemic comparator likewise requires native + nonendemic + native_nonendemic. WCVP TDWG-L3 regional native ranges alone NEVER imply island endemism.

First-level evidence: count source-recorded endemic and nonendemic species whose flower colour/form/symmetry traits are genuinely species-direct and whose independently observed accepted congener occurs in the relevant GIFT native mainland source. A 0 count is meaningful only after reading the original separate endemic status. For supported strata, report island-specific endemic-minus-nonendemic categorical mismatch against nearest recorded source congener, with regional and source equal weighting, same-island paired controls and adequate distinct phylogenetic groups. This is **phenotypic divergence proxy, not historical adaptation**: nearest checklist congener is not necessarily a sister taxon/ancestor; plant traits are species-fixed once globally.

Progression toward an actual adaptation hypothesis requires:
- Plant island/mainland homologous population trait observations (colour, structure, sex system, floral size, nectar guides/UV, autonomous selfing), ideally multiple matched localities within species or sister lineages and verified native island origin.
- Independently reconstructed tree tips/branch distances (TPL/LCVP/WP published original true-tip backbones can filter candidate relatives, not prove local adaptation).
- Visitor guild identities, quantitative visitation/pollen deposition/seed set or pollen limitation, plus climate, botanical origin, lineage and sampling effort. Compare phenotype transitions and fitness within closely related lineages across geological classes; differentiate drift/founder and environmental vs pollinator filtering.
- WCVP regional-native candidates should be separately flagged, not conflated with indigenous island populations.

## Explicit identification table

| Outcome/causal language | Current evidence |
| --- | --- |
| Different observed floral assemblage composition by geographical region and source availability | Yes, descriptive observational evidence |
| Oceanic vs continental geological effect on floral trait richness after strong independent source/taxon controls | Unresolved; strict geographical matching has at most three pairs per region |
| Establishment success/failure of specific floral guild species from propagule rain | Not measured |
| Real species-level endemic floral divergence relative to source congeners | Requires correctly read `endemic_status`, source support and same-family/phylogenetic gating |
| Pollinator-mediated adaptive evolution, island phenotype change and real guild shifts | Not identified in the present fixed-species trait dataset |

Scientific PR remains Draft; no geology-treatment randomization, historical colonization, within-species trait evolutionary response or field-measured pollination effect may be claimed.
