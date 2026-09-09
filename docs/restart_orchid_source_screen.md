# Orchid reproductive-biology source screen

## Executed follow-up: 2025 workbook and original Reunion study

The main acquisition pass retrieved the latest API-resolved record [18208949](https://zenodo.org/records/18208949), version 1.4.0. Its deposited publication date is 2025-01-10 despite the filename describing coverage through 2025; retain this inconsistency, do not infer chronology from it. `Pollination List Thru 2025.xlsx` MD5 `fac0a0e1d069bac0a1aec3761e8a7c7f` independently matches Zenodo. SHA256 `2d2baaf32eb9e2332117e67fbcc5f9d20a60f202b8bda7935b75df237b2d11a8`.

The species sheet contains 3,175 data rows. Only the genus heading is carried forward; species epithets, trait values and blanks are not filled. Exact binomial intersection with formally verified Run 34103212755 gives 432 rows with empty reproductive axes, of which 198 have an SC or SI flag. These are discovery candidates, **not 198 promoted cells**. No synonym/fuzzy expansion was applied.

The [original Reunion paper](https://www.biw.kuleuven.be/lbh/lbnl/ecology/pdf-files/pdf-art/hans%20j/JBI2005.pdf), DOI 10.1111/j.1365-2699.2005.01307.x, was downloaded (11 pages; SHA256 `4821aa20302ca818b83ff6e76e0c755f01b872b042e52f8d76a83e67458a77be`). Methods p.1753 and Appendix 1 pp.1759-1761 were read; all three appendix pages were visually inspected. Methods report bud-stage bagging and fruiting outcomes. They also make a broad all-species SC statement attributed to unpublished data, but do not give species-level hand-self/cross results. Do not claim those experimental details were verified.

Ten exact-master, previously empty reproduction cells have an explicit CP classification: Aeranthes arachnites, A. strangulata, A. tenella, Beclardia macrostachya, Cryptopus elatus, Graphorkis concolor, Habenaria sigillum, Jumellea recta, J. triquetra, Polystachya cultriformis. CP is defined in the table footnote as dependence on pollinators for fruit set, under the described bagging experiment. These ten are admitted as population-scoped autonomous_selfing_capacity=absent, High original-study evidence. This is not an inference from missing values. Per-species n and raw counts are unavailable and this limitation is retained in every receipt. No SC/SI conversion, pollinator-syndrome conversion, or genus training is performed.

The two empty-axis Cynorkis rows coded cleistogamous are held because the current ontology distinguishes facultative/obligate, but the table does not. AP alone is also held because fruiting after bagging does not distinguish sexual selfing from agamospermy. The workbook assigns SI to Oeoniella polystachys citing Agnew plus Jacquemyn, whereas the latter's general statement says all assessed taxa appeared SC (and prints Oeniella); retain an unresolved source/name audit, not a corrected claim. These checks demonstrate why the discovery workbook must not be bulk-promoted unchanged.

After these ten admissions, the saved pending queue has 422 rows representing 421 species, including 196 rows with SC/SI flags. Prasophyllum elatum appears twice with different citations and blank compatibility fields; it remains two source rows but only one species, never two resolved cells. Thus the original intersection of 432 rows represented 431 species. The saved queue retains literal SC/SI, selfing-evidence code, notes and references. Citation groups with the most remaining rows are Jacquemyn et al. 2005 (22), Wu et al. 2023 (21), Smith 1928 via Catling (13), Schlechter 1914 via Catling (12), and Ackerman 1995 (11). These are search priorities, not independent evidence or anticipated gains. Wu et al.'s bibliography link identifies DOI 10.1016/j.gecco.2023.e02778, a natural fruit-set study; do not infer SC from its abstract's natural fruiting. Workbook flags remain unpromoted pending the original-method check; the source screen is not a claim of 100% original-paper verification.

Screen date: 2026-09-07. Scope: source definitions and provenance only; no species evidence admitted, no trait assignments changed, no scientific result produced. Workbook/ledger intersection is a separate task.

## Source and field semantics

[Ackerman et al. (2023), Methods: Breeding systems](https://academic.oup.com/botlinnean/article/202/3/295/7076252), DOI [10.1093/botlinnean/boac082](https://doi.org/10.1093/botlinnean/boac082):

- SI: experimental hand-selfing gives 0% fruit set **or** <5% seed set, with cross-pollination comparisons. SC exceeds those thresholds; mechanism was not recorded. Threshold citation: Agnew (1986), *Plant Breeding* 97:183–186; terminology follows Neal & Anderson (2005).
- Autogamy combines autonomous self-pollination **and agamospermy**, seldom distinguished.
- Evidence `1`: pollinator-exclusion experiment. Evidence `2`: assumed from cleistogamy, exceptional fruiting without known pollinators, direct anther pollen-tube growth, absent rostellum, extra anthers, pollinarium bending onto stigma, or falling pollinia.
- Chasmogamy: open flowers without selfing evidence, or known pollinators alongside selfing. Mixed entries mark autogamy, chasmogamy and mixed columns.
- Literature cutoff for published analyses: 2020-12-31. Searches used compilations, databases, colleague networks and reference chasing; some photographic unpublished observations were added. Redundant references were not exhaustively retained.
- Geographic sampling and inferred autogamy introduce bias. Article license: CC BY-NC 4.0.

These definitions are article-level; literal workbook headers, blanks and row-specific notes still require inspection.

## Version and license pins

Live primary metadata: [Zenodo 7263689 API](https://zenodo.org/api/records/7263689) and [Zenodo 14601785 API](https://zenodo.org/api/records/14601785).

| Record | Version/date as deposited | Files relevant to lineage | License |
| --- | --- | --- | --- |
| [7263689](https://zenodo.org/records/7263689) | Version `29 Oct 2022`; metadata publication date `2022-03-13` (retain both, do not silently reconcile) | `Pollination_List_RLT_Data_For_Submission.xlsx`; `Pollination List.xlsx`; `Table S2 Pollinator List Literature Cited.docx`; `Pollination List Literature Cited Updated.docx` | CC BY 4.0 |
| [14601785](https://zenodo.org/records/14601785) | Version `1.3.0`; publication date `2025-01-05` | `Pollination List Thru 2024.xlsx`; `Pollination List Literature Cited thru 2024.docx` | CC BY 4.0 |

Both records belong to concept DOI `10.5281/zenodo.6350595`. Record 14601785 declares `isNewVersionOf` DOI `10.5281/zenodo.10471963`. Both queried records report `is_last=false`; neither should be described as the current latest release. This screen does not switch to a different version.

Publisher and deposited-data licenses differ: dataset reuse permission is not permission to redistribute the publisher article under CC BY alone.

Zenodo-declared checksums (MD5 metadata, not independently recomputed here):

- 7263689 submission workbook: `83cec2710db5a39049e37171c575517a`.
- 7263689 `Pollination List.xlsx`: `5c76d0b0aae942f05b6422f71279bdb1`.
- 7263689 Table S2 bibliography: `45af9e31d0794ccf61f52b3b050bd23f`.
- 7263689 updated bibliography: `4573382d1428a27006ef8884bbcbf37b`.
- 14601785 workbook: `a99bc5894b69c3a47f5663b3af118fca`.
- 14601785 bibliography: `0488ada69eccfe07e246e61edfe72b24`.

## Admission implications (screening judgments, not source claims)

1. Keep self-compatibility, autonomous sexual selfing, pollinator independence and agamospermy distinct. Do not derive one from another or treat family-level prevalence as species evidence.
2. A positive autogamy cell, even with experimental evidence `1`, does not resolve the source's selfing/agamospermy ambiguity. Require the species-linked original experiment before admitting specifically sexual autonomous selfing.
3. Chasmogamy or an unmarked autogamy cell is not a demonstrated autonomous-selfing negative. `2` remains inferred, not experimental.
4. Follow lineage as pinned workbook + sheet/row + literal field/note + cited author/year -> matching version's bibliography -> original publication and treatment/result locator. A row's citation list may cover several traits; attribution to the breeding-system result must be checked explicitly.
5. The companion bibliography is a discovery aid, not an independent replicate or substitute for species-level methods/results. Mixing the 2022 workbook with the 2024 bibliography without checking citation identity would obscure provenance.
6. Freeze version DOI and downloaded-file SHA-256 for any later extraction; preserve raw values, conflicting reports, uncertainty and taxonomic name matching. No automated synonym expansion or species-row promotion was performed in this screen.
7. Evidence codes `1`/`2` are defined for autonomous selfing, **not** as a general quality grade for SC/SI. An SC entry accompanied by autogamy evidence `2` is not thereby shown to be inferred SC. Conversely, the article's methods statement does not independently audit each SC cell, particularly entries added after its cutoff. No article statement allowing SC to be inferred directly from autogamy was found in the checked compatibility/methods passages; original treatment evidence remains necessary.

## Access and remaining checks

The publisher's normal HTML endpoint was intermittently blocked; its publisher-owned minimal article endpoint supplied readable methods and copyright text. Zenodo landing pages failed through the browsing tool, while the official record API returned metadata. These access differences are not biological evidence. Original species papers and workbook rows were not audited by this source-screen task.
