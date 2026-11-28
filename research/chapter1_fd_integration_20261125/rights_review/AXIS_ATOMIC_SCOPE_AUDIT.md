# Exact original axis provenance versus atomic Direct support

26 November 2026 JST. New local analysis of SHA256-pinned original archives 10058653212, 10042377593, and 10186293306. This is a correction to the scope of the earlier CORE_RIGHTS_SCOPE crosswalk, not a replacement scientific result. The new independent CI is pending.

## Exact axis identity closes the original rights crosswalk

The original analysis ZIP contains `snapshot/chapter1_trait_snapshot_species_axis.csv.gz`, `snapshot/chapter1_trait_ledger_direct_only.csv.gz`, and the corresponding state audit. The snapshot manifest pins v1.0 species-axis SHA256 a6eef8d731b2730a99e388ddd0683d7ac9b2af30d43f1f9569c857f89f664b2a. The 318,885 rows of the snapshot and original v1.0 axis ledger are identical after sorting by the unique species-axis key, in every column (composition, names, source groups, lineage sets, quality). The compressed bytes differ because of snapshot serialization, not records. The original source atomic ledger inside the analysis ZIP also has exactly the same SHA256 as the v1.0 atomic ledger.

All 12,401 informative single-state reproductive records of the frozen H1 Direct audit match the historical long snapshot exactly. They correspond to 9,186 species-axis cells (1,278 authorized and 7,908 review-required). At the trait-record level, 1,294 records correspond to authorized cells and 11,107 to rights-review-required cells. The previous 1,747 unmatched trait records were unmatched to the atomic ledger, not missing from the canonical axis rights audit.

## Material evidence-scope distinction

The exact implementation `src/island_v2/chapter1_trait_snapshot.py::materialize_long_ledger` copies the cell's quality, source_groups and source_lineages to every component trait. Therefore the historical Direct scope is an axis-level High/Medium restriction. It must not be described as guaranteed individually verified trait-specific evidence.

| Trait | H1 informative single state (axis Direct) | Same exact state in High/Medium atomic ledger | Atomic row absent |
|---|---:|---:|---:|
| SC/SI | 6671 | 6556 | 115 |
| Mating system | 2054 | 2047 | 7 |
| Autonomous capacity | 3676 | 2051 | 1625 |
| Colour | 26127 | 25524 | 603 |
| Form | 6934 | 6926 | 8 |
| Symmetry | 32039 | 32036 | 3 |
| Tube depth | 2566 | 2565 | 1 |

No same-key state/quality disagreements in these informative single-state rows. All 10,654 exact reproductive atomic lineage sets are subsets of the corresponding axis lineage sets; only 4,527 sets are equal. In particular, a Go Botany lineage attached to a reproductive axis does not establish that it is the direct source of a particular SC or autonomy value. The earlier individual-source trait counts remain appropriate for the *atomic* ledger; the axis blocker counts must now be interpreted on the full exact axis crosswalk.

For the H1 binary definitions, the three-informative-reproductive-component species roster is 564 at the axis scope, 504 with exact atomic support for all three. The informative binary symmetry/depth roster is 1,040 vs 1,039; this is **not** the full separate paired-trait analysis roster and should not be substituted for its 1,112 species without reproducing its definition. The missing atomic records are not proved false; reconstruct the integration rule and inspect evidence before drawing that conclusion.

## Rights and manuscript implications

The current frozen GEB candidate uses wording such as species-direct and directly observed components. This needs a scope clarification, and strict atomic-support sensitivity before strengthening the direct-evidence claims. Do not silently change the published scientific selector or discard original results. First reproduce the seven primary broad-flora H1 Direct fits using identical frozen settings and the stricter roster, and repeat the relevant H2 complete-three-component result and paired-trait support checks. The exact historical WCVP inputs remain missing; do not replace them. H3 remains an independent published experiment analysis; H4 reuse of Direct scores requires a separate scope audit if affected. Retain all original H1/H2 results and disclose the new sensitivity as post-hoc.

The Europe PMC statement audit run 45755161777 remains in progress at this checkpoint. The original 87 candidate cells are not authorized by this discovery; it changes provenance interpretation, not source permissions. No external correspondence, scientific fit, licence promotion or publication.
