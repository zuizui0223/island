# Rights review: marginal cell clearance rather than raw source frequency

2026-11-26 JST. Actual computation on the canonical frozen rights-audit ZIP (artifact 10186293306, SHA256 6c5f2577e85e95e841b8659427e27c72ce631f59a3da5a335f54f505de5ea32a) and provenance audit (10154261753, SHA256 5d3e5ea80e1cf9757cc7b3bbb64db2ca67f0b8ecce7f7f52668472511f07ad36), independently downloaded again from authorized private Drive and matched to source hashes. No source rows exported. The original rights policy remains unchanged.

## Baseline

| Axis | Resolved | Already authorized |
|---|---:|---:|
| Floral structural complexity | 91,635 | 26,825 |
| Flower colour | 82,556 | 18,067 |
| Reproductive assurance | 48,497 | 1,382 |
| Total | 222,688 | 46,274 |

176,413 cells have one or more rights blockers; one has a separate unresolved rule-semantic discrepancy. 5,032 blocker families form 7,898 distinct blocker sets. 63,174 rights-review cells have exactly one blocker (12,755 reproductive). High/Medium/Low refer to frozen axis quality, not an exact Direct-only trait roster. The published derivative is strongly unbalanced; it is not an acceptable silent replacement for the scientific ledger.

## Counterfactual clearance

For each cell c let B_c be its exact original blocker set. For a *hypothetical* fully cleared family set S, additional releasable cells = sum of I(B_c ⊆ S) over rights-review-required cells only. No assumption that every family represents one licensor, that all members have the same licence, or that a clearance is legally attainable. Semantic/provenance gates and attribution/licence compatibility remain separate. Greedy ranking maximizes one-step completion of whole blocker sets; it is a heuristic, not global optimization.

| Scenario | Exact blocker labels | Additional cells | Potential authorized total | Additional reproductive cells |
|---|---:|---:|---:|---:|
| Three eFloras/FoC labels (separate permissions) | 3 | 14,001 | 60,275 | 38 |
| Four European database labels | 4 | 11,872 | 58,146 | 3,132 |
| Four item-license review candidates (WFO Malesiana, Baseflor, GBIF, GrassBase) | 4 | 11,136 | 57,410 | 3 |
| Top ten by raw incidence | 10 | 39,326 | 85,600 | 7,890 |
| Top ten by immediate singleton gain | 10 | 40,357 | 86,631 | 6,872 |
| Greedy ten, all axes | 10 | 42,856 | 89,130 | 10,082 |
| Greedy twenty, all axes | 20 | 68,436 | 114,710 | 16,299 |
| Greedy four, reproductive objective | 4 | 16,052 | 62,326 | 10,031 |
| Greedy eight, reproductive objective | 8 | 24,398 | 70,672 | 16,224 |

These are planning ceilings conditional on full legal clearance of the listed blockers. None have been granted. Clearing eight reproductive-priority labels could take the reproductive axis from 1,382 to 17,606 cells (36.3%), whereas morphology-heavy permissions alone yield almost no reproductive gain.

## Reproductive-first queue

| Step | Frozen blocker label | Newly eligible reproductive cells, conditional on previous steps |
|---:|---|---:|
| 1 | Go Botany | 3,799 |
| 2 | Ken Fern Useful Plants | 2,795 |
| 3 | FloraWeb/BiolFlor | 1,850 |
| 4 | Pladias | 1,587 |
| 5 | Unresolved citation hashes | 1,344 |
| 6 | Unresolved other DOIs | 1,527 |
| 7 | Europe PMC | 2,336 |
| 8 | Ecoflora | 986 |

The unresolved hash and DOI labels represent provenance-recovery work, not one rights holder or a promise that clearing every constituent record is feasible.

## Exact source scope and licensing evidence

- GBIF: historical provenance audit has 1,867 direct lineage records under domain:gbif.org, all species-page URLs, not occurrence downloads. GBIF's [dataset licensing guidance](https://www.gbif.org/publishing-data) and [taxon description guidance](https://ipt.gbif.org/manual/en/ipt/latest/best-practices-checklists) show why species-page provider/content identity must be recovered. No domain-wide CC0/CC BY clearance.
- Flora Malesiana: 3,393 mapped FLORML lineages trace to one specific [WFO content archive](https://files.worldfloraonline.org/files/Flora%20Malesiana/github/Flora%20Malesiana%20final%20files/Flora%20Malesiana%20final%20files.zip). The [WFO terms](https://www.worldfloraonline.org/termsOfUse) distinguish CC0 taxonomy from provider-specific content. Some individual [Naturalis publications](https://repository.naturalis.nl/pub/532515/) are explicitly CC BY, but do not establish that this imported treatment archive and all its records share that licence. Recover exact version and record-level licence metadata first.
- [Go Botany terms](https://gobotany.nativeplanttrust.org/terms-of-use/) restrict reuse of site content and identify a contact for uses beyond images. Request scope-specific written clarification for normalized trait codes, derived database distribution and editor/reviewer access, not blanket scraping or copying of site prose.
- [Pladias management rules](https://pladias.cz/en/homepage/rules) permit specified research uses but retain provider decisions and require board contact for substantial use. Distinguish individual source subsets, attribution/cooperation, sharing with reviewers and public release.
- [Useful Tropical Plants](https://tropical.theferns.info/donate.php) declares CC BY-NC-SA 3.0 for its database/code. This is not the same source as every PFAF version. Do not assume that current umbrella CC BY-SA covers NC material; verify actual snapshot and negotiate permission or a compliant separation.
- [BiolFlor provider](https://www.ufz.de/index.php?de=38567&m=0) and [FloraWeb original source descriptions](https://www.floraweb.de/ueberfloraweb/datenquellen.html) distinguish downstream portal from original trait data. Obtain exact version and rights for reproductive records.
- [Baseflor/eVeg licence notice](https://e-veg.net/eveg/web/en/about-eveg-project) identifies Baseflor as Open Database licensed. [ODbL terms](https://opendatacommons.org/licenses/odbl/1-0/index.html) require attention to adapted-database sharing; [CC compatible licences](https://creativecommons.org/compatible-licenses/) does not list ODbL as a BY-SA-compatible licence. Confirm the actual frozen acquisition and permissible distribution architecture; do not silently add it to the current umbrella.
- [Kew data terms](https://www.kew.org/science/collections-and-resources/data-and-digital/terms-of-use) require dataset-specific licence metadata. Recover GrassBase-specific terms before promoting its 1,427 sole-blocked cells.
- [Missouri Botanical Garden conditions](https://www.missouribotanicalgarden.org/conditions) restrict reuse of its hosted texts and acknowledge other publishers. The three eFloras-related families require exact treatment-provider scope and separate rights decisions; a common provider umbrella is not enough.

No legal clearance or source-policy changes are made from these public notices. The next practical goal is a source-level scope/attribution register and an editorial access arrangement, not reanalysis on a rights-selected sample. Author/institutional approval is needed before external messages or public releases.
