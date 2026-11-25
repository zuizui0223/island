# Source-provider robustness and support feasibility of island flower–visitor guild associations

Original data verified **18 November 2026** in [GIFT×GloBI CI #44301575888](https://github.com/zuizui0223/island/actions/runs/44301575888), with source archive and 4,505 accepted botanical flora SHA checks and successful historical 928/927/552 old cohort replay. The original GloBI source has **24,577,183** relation rows and the full source result includes **750 source-mode × region × provenance × genuine same-accepted-species paired-flower × visitor-guild cells** (2 × 5 × 3 × 5 × 5). **No cell is an independent experiment.**

## Explicit source-documentary support thresholds (exploratory, post hoc)

For each genuine same-species Direct flower predicate × explicitly documented GloBI guild, require the same region and botanical status panel to have:

- >= **10** GIFT source islands with at least one positive cited relation;
- >= **10** distinct accepted botanical species with the true Direct floral pair;
- >= **10** original source bibliography labels;
- >= **2** original archive-first provider labels.

Test once on full source and once after removing the **globally most frequent original provider archive**, selected without looking at flower or guild-specific outcomes. A cell only passes **provider-robust documentary candidate support** if both sets independently meet the same minimum, *on the selected cells' original botanical scope*. These thresholds were chosen as a practical feasibility screen **after seeing dataset support**, not preregistered biology; passing would still not verify any primary publication independence, effective pollination, animal nondetection, taxonomic lineage convergence or island floral evolutionary change. A failed gate is **NOT IDENTIFIABLE**, not ecological no-effect or absent pollinator.

## Actual original-data support counts

There are 25 original selected two-Direct-flower × GloBI animal-guild combinations per geographical/plant-status panel. The complete original anonymous 750-cell artifact from [Actions #44301575888](https://github.com/zuizui0223/island/actions/runs/44301575888) was evaluated against the above same four counts:

| Original geography & plant provenance | Original 25 candidate cells passing threshold on full cited source | Passing after globally largest archive excluded |
| --- | ---: | ---: |
| **Global all-recorded** | **14 /25** | **0 /25** |
| Global literal known-native | 1 /25 | **0 /25** |
| Northern midlatitude all-recorded | 11 /25 | **0 /25** |
| Northern midlatitude literal known-native | 0 /25 | **0 /25** |
| **Tropical all-recorded** | **0 /25** | **0 /25** |
| **Tropical literal known-native** | **0 /25** | **0 /25** |

**Important**: The full 750-cell tensor includes Diptera, which was NOT included in the concise initial source-action log printout; that subset had only 12 global and 10 northern-midlatitude passing combinations. The complete 25-cell totals are **14** and **11**, respectively, while the post-omission totals remain zero. A source's independent archive label is NOT itself an independently published field study. Thus even those 14 preliminary global cells cannot currently be promoted to valid tests of realized guild preference or floral adaptation. After omitting the global largest source, ZERO remain for this exploratory minimal gate; the result is source-availability fragility, not proof that flower shape and animal visitor guilds have no ecological association.

## Why tropical pollinator-syndrome displacement remains not identifiable

In the original true same-accepted-species source-linked tropical botanical floras: **actinomorphic × deep / Lepidoptera** has 7 positive citation relations on 7 actual GIFT islands/7 directly scored plant species, but all from **one archive-first provider**. All disappear after excluding the globally largest provider. Tropical actinomorphic × deep / birds has 2 original positive relations on 2 islands/2 plants with one provider; strict literal native is only one island/one plant. Tropical zygomorphic × deep / birds has 8 original cited records on 6 islands, but its record count drops to 1 island/plant after global-largest-provider omission.

GloBI's original five-source-guild ontology labels only combined **Lepidoptera**, not butterfly vs Sphingidae. A global flower-colour ontology may combine red and pink, or other distinct signal colours. Neither is sufficient to identify an independently measured bird/butterfly/hawkmoth flower syndrome. These positive citation records do NOT measure population-specific flower phenotype, actual pollen transfer, fruit/seed set or historical change.

## Independent data aggregators are not independent underlying ecology

Web of Life on exact full-parent geographic coverage supplies **78 network sites on 18 GIFT islands /2,529 positive botanical accepted plant×network×guild units**, after the original 75/16/2,469 baseline. GloBI provides **232 islands /1,799 accepted plant×island pairs with cited animal relationships**, but the two *sets of source island sites, plant identities, references and original studies* have NOT been de-duplicated across WOL, GloBI and Mangal. Do not add the island or interaction counts or calculate a pseudo-replicated pooled effect; their observed animal taxonomy resolution and study definitions also differ.

**Biological interpretation ceiling:** keep globally replicated botanical functional-repertoire composition, real same-species floral pair frequency and genus-composition source accounting as the main paper. Source-positive visitor links help identify where to obtain independent original site/time and animal effort studies, but cannot independently validate adaptive island pollinator-guild turnover. The new source-side gate software [visitor_documentation_support.py](../../src/fd/visitor_documentation_support.py) and [unit tests](../../tests/test_visitor_documentation_support.py) formalize these evidence qualifications; the original anonymous-artifact verification workflow is [island CI #44306963725](https://github.com/zuizui0223/island/actions/runs/44306963725). The second-stage source artifact replay [original CI #44308473532](https://github.com/zuizui0223/island/actions/runs/44308473532) **SUCCESS**: 3 strict synthetic source gate tests passed, both original GloBI and WOL GitHub artifact ZIP SHA256 digests reconciled, and full 750-cell original input evaluated (GLOBAL 14 before/0 after provider leaveout, TROPICS 0/0, all strict-native source gate 1/0). Anonymous result artifact **#14428704199**, signed ZIP SHA256 **701d0c5dfe7a8286d435423b56ca55dfcfdf9751aecab0667c2784b1178152b9**. This second-stage result reuses the SAME source data; it is NOT a new biological dataset or independent replication.

All original taxon/island/site/reference/provider identifiers are held within the protected original analysis, not part of this report. The main GIFT botanical population remains 4,505 and no absence is generated by documentation failure.
