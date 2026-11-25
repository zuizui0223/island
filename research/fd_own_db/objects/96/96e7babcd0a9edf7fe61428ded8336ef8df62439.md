# 14 November 2026 — Measured FCP flower photographs × geological island endemics

**ORIGINAL protected source job completed successfully.** [Island Actions #43578227883](https://github.com/zuizui0223/island/actions/runs/43578227883), anonymous artifact #14107327790. This is a cross-repository DATA OPPORTUNITY audit, not newly measured pictures, local floral adaptation or genetic polymorphism.

## Frozen original inputs and legal/provenance scope

- Independent FCP source: **4,475 measured plant species / 447,500 actual photographic ROI measurement attempts / 227,107 four-colour classifiable photographs**, [20261021 FCP original species ledger](https://github.com/zuizui0223/fcp/blob/c45a95420a96abfb3fdb8e9ad1a37c75288c63fb/results/fcp_np_all4475_unified_species_ledger_20261021.csv), verified SHA256 23d3daadc9799c2e7b5867d54ccc1ed98369ba79a81b72dc7f831e5197e63c49.
- GIFT original recorded island flora: original accepted plant names, separate independent native-origin and endemic-status fields, confirmed combined endemic or native_nonendemic status. WCVP regional range alone cannot prove native island endemism.
- Original Roeble/Zizka GSHHG geological island crosswalk: **434 oceanic/continental typed islands**, LGM mainland connection NOT treated as geological classification. Original GIFT five closest native mainland candidate floras are not the historical ancestors of island taxa.
- Implementation: [fd source audit code](../../src/fd/fcp_geological_endemic_opportunity.py), [synthetic tests](../../tests/test_fcp_geological_endemic_opportunity.py), [island original protected runner](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_fcp_geological_endemic_overlap.py).
- Only anonymous region × geological origin × source status count tables leave the original rights-safe workflow; no exact taxa/island/photo/GIFT source IDs are exported.

## Strict native endemic source overlap

**Every count below is the number of distinct original source accepted endemic plant names within that region/type, NOT evolutionary independent species or unique island populations.** FCP exact-match means the accepted source name exactly matches the FCP species spelling; unresolved synonyms are never silently mapped to accepted taxa.

| Geological type and region | Geological islands | Distinct source-verified endemic names | Species with an exact measured FCP match | FCP species with >=40 actually classified photos | FCP endemic species with a different source FCP-photographed same-genus mainland candidate |
| --- | ---: | ---: | ---: | ---: | ---: |
| Tropical oceanic | **160** | **1,232** | **0** | 0 | 0 |
| Tropical continental | **43** | **11,775** | **2** | 1 | 0 |
| Northern midlatitude oceanic | **70** | **442** | **2** | 2 | 2 |
| Northern midlatitude continental | **109** | **1,298** | **3** | 3 | 1 |

By contrast, strict native-nonendemic island botanical floras have considerably more original FCP measured species. Tropical oceanic **266** (of which **180** have >=40 classifiable photos), tropical continental **474** (322), northern-midlatitude oceanic **649** (457), northern-midlatitude continental **1,488** (1,086). The lack of an exact tropical-oceanic endemic match is NOT proof that no online flower image exists for these 1,232 accepted names; it is only this particular source's restricted exact species-name overlap and requires cautious synonym/measurement coverage audit.

## Evolutionary interpretation and next source gate

Stage 1 / establishment: island species-to-source botanical repertoire comparisons remain observable but source native checklists do NOT directly measure propagule arrival, establishment success or extinction. The earlier strict geology-paired independent dataset has only **3 separate oceanic/continental island pairs per major region**.

Stage 2 / adaptation: FCP's species-wide photo diversity D/W and 50-km photo organisation are **not island-specific population color measurements** and the current FCP cohort has 0 exact tropical-oceanic endemic names. These source data cannot establish a geographical origin by adaptive flower-colour coefficient. Source trait Direct flower-colour data independently have 13 tropical genera shared across geological origin with 9 colour-variable genera, but the geological within-genus effect was unsupported (−0.05267, 95% genus interval includes zero; exact sign p 0.67578).

A workable next acquisition design must (a) verify indigenous endemic accepted plant names and source synonym mappings without uncontrolled promotion, (b) deliberately target photo-measured tropical oceanic endemic taxa from distinct genera and independent archipelagos, (c) retain **island-local and mainland-local** photo coordinates and dates rather than only species-wide means, (d) acquire independently surveyed visitors, pollen deposition and reproductive outcomes, and (e) separate neutral population history, phenotypic plasticity, ecological sources, camera/observer and taxonomic sampling before calling it adaptation.

Keep all parent research PRs DRAFT. The current observation is a source-coverage blocker and research-acquisition target, not an ecological negative result for oceanic endemic colour evolution.
