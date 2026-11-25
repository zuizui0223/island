# Correction: EuPPollNet true SAME-plant source flower availability and within-source comparison support

**Original source verified 19 November 2026, NOT a pollinator preference estimate.** Source original protected [island Actions #44431962303](https://github.com/zuizui0223/island/actions/runs/44431962303) **SUCCESS**, 15 synthetic regressions passed, original 4,505 GIFT botanically accepted island-flora parent, 4,499 exact-hash GSHHG geographic polygons, previously fixed 13 EuPPollNet-visited GIFT island floras and source public positive visit/flower survey Git blobs replayed. New original anonymous source artifact #14487097731, GitHub signed ZIP archive SHA256 **4b0be59b2a42045346f0eb91cbf32f273f19bb921ccf2884b6bb93c4fc200638**. Private source taxa, study/site IDs, original flowering keys and coordinates remain in protected CI. This is a **quality-control correction** of the prior [flower count upper bound](FLOWER_ABUNDANCE_SOURCE_GATE.md), not a second independent animal study.

## Why the prior 11,041 source exposure rows were an upper bound

The prior source analysis required same published \`Flower_data_merger\` key, original study/network, same calendar day and a positive numerical flower abundance value with nonmissing \`Units\`. However, it did **not** independently require \`Flower_counts.Plant_species\` to have the SAME accepted botanical species binomial as \`Interaction_data.Plant_accepted_name\` on the actual GIFT original island. That is essential for a correct **species-specific** flower abundance denominator. A source key can link the same survey/study/day while referencing a different source plant; also unresolved taxonomic synonyms cannot be counted as certifiably the same species without an explicit backed taxonomic crosswalk.

The new source check applies strict SAME accepted genus+specific epithet matching, never imputes a synonym, morphologically unrelated plant or infraspecific taxon from another source. Additional original source waterfall:

| Original positive EuPPollNet source records on correctly accepted GIFT island plants | Source rows |
| --- | ---: |
| Previously verified positive visit and raw source flower-count unique key | 21,284 |
| Same published source study×network×DAY/month/year | 11,042 |
| Positive numerical flower count, known native survey unit | **11,041** |
| Positive count but source flower plant binomial **DIFFERS** (possibly true taxon mismatch, unresolved synonym or naming practice) | **632** |
| Original flower survey plant string **UNRESOLVED** at strict accepted binomial level | **10** |
| **CERTIFIED same source botanical plant in visitor event and flower survey** | **10,399** |

The 632 mismatches are **not proven different biological plant species**, because scientific synonym and original source taxonomy reconciliation has not been attempted. The strict result represents a lower, high-certainty directly matched subset; previously reported 11,041 is the broad upper source linkage count and must NOT be used as fully species-validated flower supply.

**The updated 10,399 positive visitor rows** are repeated original source observations on only **4 original GIFT island floras /4 independently labelled EuPPollNet source studies /174 distinct accepted botanical plant species /665 distinct original published flower survey keys**, not 10,399 independent pollinator choices or local plant populations. Original flower count measurement remains heterogeneous across **nine named source Units**, including cover percent, flowers/m², heads/m² and inflorescences/m². No cross-study raw abundance value is treated as commensurate without verified physical-unit calibration.

## Strongest possible within-source true phenotype contrast opportunities

After strict same-source accepted flower species matching and before any response model, group by the original **source study × network × date × GIFT island × unconverted flower-count Units**. This yields **166** original site/day/unit contexts overall, of which **122** contain >=2 accepted plant species, **124** have >=2 documented original visitor animal guilds, **121** span >=2 original botanical genera. These are separate marginal counts, **not guaranteed simultaneous**.

The actual joint comparison support requires at least one directly scored TRUE flower configuration and at least one explicitly scored FALSE configuration (both character axes recorded on the same accepted species), >=2 documented visitor guilds, and >=2 botanical genera in that identical original context:

| Genuine SAME-species Direct botanical flower profile (versus explicitly scored nonprofile) | Qualified original site/day/unit contexts | Original GIFT islands | Independent labelled source study IDs |
| --- | ---: | ---: | ---: |
| Actinomorphic × deep tube | **7** | **2** | **2** |
| Zygomorphic × deep tube | **2** | **1** | **1** |
| Tubular/bell × deep tube | **6** | **1** | **1** |
| White × tubular/bell | **0** | **0** | **0** |
| Red/pink × tubular/bell | **0** | **0** | **0** |

The seven radial×deep source contexts on two islands and two studies provide **at most** an exploratory opportunity to inspect within-study guild use after unit calibration, original source sampling completeness, and plant/visitor observation-effort harmonization. They are not enough to claim a transferable, source-independent macroecological pollination syndrome: the same plant may be visited by multiple guilds; positive literature interaction entries are not standardized true absence/nonvisitation, the original flora's species-Direct morphology is not observed at the original source field flower, and no pollen delivery / fitness is measured.

## Research decision and sources

This is a **new confirmed identification result:** genus-independent floral trait and animal guild observationally comparable source contexts become rare even when original 623,476 European visitor records, published flower counts, 4,505 original GIFT botanical island floras and each true same-accepted-plant phenotype are jointly included. Do NOT estimate a global or tropical butterfly/bird/hawkmoth plant adaptation effect from the surviving records. A species-matched source abundance denominator alone is not proof of a complete choice set of all flowers, animal non-detection, or a biologically comparable unit of pollination service.

The original global macroecological **functional floral repertoire recurrence + regionally asymmetric true floral architecture + genus-source assortment** remains the main positive finding. This source qualification is a supplementary data-support audit and should not shrink the global botanical inference population to four islands.

[Reusable original same-plant exposure matching audit](../../src/fd/euppollnet_species_identity_gate.py) · [synthetic tests](../../tests/test_euppollnet_species_identity_gate.py) · [source original protected runner](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_original_euppollnet_flower_supply_on_gift.py) · [verified original CI](https://github.com/zuizui0223/island/actions/runs/44431962303).

**Next falsification required:** audit published full site/day original flower censuses **including plant species with no recorded visitor edges**, at the same survey unit, to establish whether a genuine plant flower-availability *choice set* is sampled; then check whether study/network IDs locate one real geographic site and one botanical sampling scheme. Until that passes, no visitor preference coefficient or fitted positive biological effect.
