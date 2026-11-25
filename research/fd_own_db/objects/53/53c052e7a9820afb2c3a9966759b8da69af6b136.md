# 14 November 2026 — Original iNaturalist public endemic photo opportunity on oceanic versus continental islands

**Genuine public-provider query executed; neither photos nor flower pixels were downloaded or measured in this pilot.** Original GitHub Actions [#43580304367](https://github.com/zuizui0223/island/actions/runs/43580304367) completed SUCCESS, anonymous source JSON artifact #14108395738. The prior fixed first-eight-genera pilot [#43579956771](https://github.com/zuizui0223/island/actions/runs/43579956771) also passed.

## Frozen pre-API biological sample

Select among independently Roeble/Zizka-typed **tropical oceanic and continental** native-endemic GIFT plants. Require original high/medium Direct species-specific flower primary colour, five original GIFT native candidate mainland source floras, and a separate source-listed species of the same genus with directly recorded colour. Do NOT inspect photo availability or actual floral colour state to choose species. Exclude taxa found on the focal island from all candidate source donor lists.

There are **13 botanical genera** with qualifying original endemic plants on BOTH geological island types. Deterministically select ONE different accepted endemic plant name per type per genus, 26 source taxa, using the same frozen source SHA order as the first eight-genera pilot. The selected cohort covers **9 oceanic and 5 continental 10-degree geographic blocks** and **11 independently named source archipelagos**. Distinct botanical genera are not independent historical colonizations or phylogenetically disjoint tree branches.

## Actual public API source returns

Rate-controlled iNaturalist v1 endpoint: first exact active species-level taxon lookup by accepted original species name, then research-grade photographs with public georeferencing and one of CC0, CC-BY, CC-BY-SA, CC-BY-NC, CC-BY-NC-SA licenses. First eligible observation page per exact taxon; reject obscured/privately georeferenced observations, positional accuracy greater than 5 km and misidentified/descendant taxa. Provider query failures MUST remain provider failures, not biological no-photograph states.

| Source outcome | Tropical oceanic endemic taxa | Tropical continental endemic taxa |
| --- | ---: | ---: |
| Preselected original source taxon queries | **13** | **13** |
| One unambiguous exact public iNat species taxon ID | **11** | **11** |
| At least one verified licensed, open georeferenced research-grade photographic observation in first page | **4** | **3** |
| Individual photographic observation records passing the same first-page checks | 27 | 20 |
| Provider API errors | 0 | 0 |

**One shared original botanical genus has at least one qualifying global public photographic record for the oceanic endemic species AND its continental endemic congener.** The first eight-genera pilot had 2 oceanic /1 continental photo-supported species and zero pairs supported on both; the full 13-congener extension reveals a single matched pair.

CAVEAT: qualified source photographic observations may show leaves, not a flower; their publicly visible coordinates need not fall inside the actual species' source island; the source photo may be cultivated despite original GIFT source status elsewhere; some photographs are on later API pages; the exact iNaturalist taxon ID can differ from a CoL-reviewed botanical concept. No photo pixels, pixel SHA, ROI colour metric, effective pollinator evidence or native island population genetic trait has been measured in this pilot. A photo-availability contrast is NOT a geology-specific adaptability result.

## Next essential geographic and image-quality filters

- Query the ORIGINAL target species with public geoprivacy open around the original independent GSHHG island centroid (source-area-scaled bounded radius), then require point-in-the-ORIGINAL GSHHG polygon before claiming island locality. Never reverse obscured-sensitive species locations; a centroid circle is NOT polygon proof.
- On genuinely source-island-local CC-licensed photographs, retrieve original image bytes, verify pixels, inspect flower presence and use the same calibrated FCP ROI-v4 / nine-composition measurement; retain source license/photographer attribution. Do NOT reclassify generic images as flowering or count ROI failures as negative pigment states.
- Increase genuinely independent botanical genera and geologically independent archipelagos before testing an oceanic-versus-continental adaptation coefficient. Botanical sister lineages, source population gene/phenotype and actual pollination/fitness measures are not yet present.

Code: [source-paired botanical/taxon and public photo pilot](../../src/fd/oceanic_endemic_inat_photo_pilot.py), [unit tests](../../tests/test_oceanic_endemic_inat_photo_pilot.py), [original private runner](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_oc_endemic_inat_photo_pilot.py), [workflow](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/.github/workflows/run-fd-oc-endemic-photo-pilot.yml). Keep original fd/island PRs DRAFT.
