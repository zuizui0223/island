# 14 November 2026 — Frozen FCP photo species and CoL XR synonym candidate audit

**Original-source completed and verified:** [island Actions #43579566802](https://github.com/zuizui0223/island/actions/runs/43579566802) SUCCESS; anonymous artifact 14107178629. This is a photo-coverage/taxonomic-candidate audit, NOT a measurement of oceanic/endemic adaptation or island-local flower colour.

## Frozen independent source frames

- Original Catalogue of Life Extended Release crosswalk: **115,328** unique original submitted botanical names, original 2026 CoL XR checklist UUID 7ddf754f-d193-4cc9-b351-99906754a03b, source file SHA256 079d9c27b0ca3aa51f0c7310f43ad4901b07a319d6063eab5a7d62b291bdd9e3. This is an AUTOMATIC CANDIDATE name normalization layer, not human-reviewed accepted taxonomy.
- Actual FCP measured-photo source: **4,475 species**, **447,500 measured photos**, **227,107 classifiable photo colour states**, original source-ledger SHA256 23d3daadc9799c2e7b5867d54ccc1ed98369ba79a81b72dc7f831e5197e63c49.
- Across the WHOLE 115,328 name source list (not the endemic set), **3,731** submitted botanical names already match an FCP photo taxon via exact spelling with whitespace-only normalization. A further **150** original names propose a CoL automatic EXACT synonym-to-FCP-accepted-name candidate. **129** candidates are excluded from unique name recovery for many-submitted-name-to-one-target FCP collision or existing exact name conflict. **21** remain one-to-one unpromoted *candidate* aliases, none automatically accepted as botanical identity.
- Strict candidate qualification: source exact synonym status, candidate species rank in Plantae, exact match type and confidence >=95, no flags/issues, and no target FCP label collisions. No name guesses, genus-prefix promotions, manual undocumented mapping or source native-status changes.
- Geological labels from independent Roeble/Zizka crosswalk of **434** original islands. Island biological endemic status must be independent original native, endemic and combined floristic status evidence; never inferred from WCVP TDWG regional range.

## Native endemic geological target — full strict synonym candidate result

| Geological region | Original distinct strict native endemic species | Exact photo-FCP match | One-to-one CoL synonym-to-FCP candidate |
| --- | ---: | ---: | ---: |
| **Tropical oceanic** | **1,232** | **0** | **0** |
| Tropical continental | 11,775 | 2 | 0 |
| Northern midlatitude oceanic | 442 | 2 | 0 |
| Northern midlatitude continental | 1,298 | 3 | 0 |
| Southern extratropical oceanic | 64 | 0 | 0 |
| Southern extratropical continental | 493 | 1 | 0 |
| Northern high-latitude oceanic | 2 | 0 | 0 |

Every geological-endemic region group has **zero UNAMBIGUOUS one-to-one CoL automatic candidate synonym rescue**. This does not prove zero images exist online, nor rule out unresolved source synonyms, the 129 many-to-one proposed aliases, nonexact botanical matches or names outside the CoL XR universe. In the current high-provenance automatic name gate, the original FCP species-wide flower colour D/W cannot be used to infer tropical-oceanic endemic island-local colour variation.

## Next acquisition, and interpretation limits

Proceed with separate provider-photo availability and exact species-ID verification among deliberately selected tropical endemic **oceanic/continental genus-paired** source candidates; select before observing photo response or flower-colour category. A licensed public georeferenced research-grade photograph is not itself a photograph of a flower, an island-local phenotype, observed pollen service, or evolutionary selection. Retain provider response status and source date and do NOT reinterpret missing photo data as biological absence.

[Source taxonomic source-audit module](../../src/fd/fcp_endemic_colxr_synonym_gate.py) · [unit tests](../../tests/test_fcp_endemic_colxr_synonym_gate.py) · [original protected runner](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_fcp_colxr_endemic_synonym.py) · [workflow](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/.github/workflows/run-fd-fcp-colxr-endemic-synonyms.yml). Scientific PR remains DRAFT.
