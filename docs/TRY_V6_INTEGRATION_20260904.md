# TRY v6 private-request integration

This integration adds TRY as a **private local evidence provider** for the v2 trait ledger. It does not copy the raw TRY request into the repository and it does not treat TRY as one independent biological source. Original references are retained and collapsed to original-source lineages before direct-cell resolution.

## Accepted mapping contract

The adapter currently promotes only the following source-backed species-level mappings:

| TRY TraitID | TRY field | island trait | rule |
|---|---|---|---|
| 207 | Flower color | `flower_primary_color` | conservative multilingual colour mapping; multiple states within the same original source become `multicolored_variable` |
| 2935 | Flower symmetry type | `floral_symmetry` | `radial -> actinomorphic`; `bilateral -> zygomorphic`; boolean `yes` is positive evidence; boolean `no` is not converted to the opposite state |
| 211 | Flower sexual self-incompatibility | `self_incompatibility` | species-level SC/SI and explicit SI mechanisms only; genus/family records and inbreeding measures are rejected |
| 2936 | Flower corolla type | `corolla_fusion` sidecar | `free` or `partly_or_fully_fused`; **never mapped to `floral_form`** |

The first three traits enter the existing axes `flower_colour`, `floral_structural_complexity`, and `reproductive_assurance`. `corolla_fusion` stays outside the current three-axis strict coverage contract until it has its own ontology/analysis contract.

## Provenance and conflict rules

1. `AccSpeciesName` is matched to `accepted_species` first. Exact case-normalized names are treated as formatting-equivalent exact names. Optional exact synonym maps may be supplied.
2. TRY `SpeciesName`, `AccSpeciesName`, `DatasetID`, `Dataset`, `ObservationID`, `ObsDataID`, `DataID`, `DataName`, and original values remain in the private prepared outputs.
3. `Reference` is the preferred original-source lineage. DOI is used when available; otherwise a normalized citation fingerprint is used. Dataset identity is only the fallback when no citation exists.
4. Evidence is pre-aggregated at `accepted_species x trait x original-source lineage` before it reaches the shared resolver.
5. Multiple colour states inside one original source are retained as `multicolored_variable`. Conflicting symmetry or SI/SC states inside one original source are fail-closed and written to `try_exclusions.csv.gz`.
6. The shared `dedupe_direct_lineages()` and `resolve_direct_cells()` logic is still authoritative after TRY is added, so overlap with a paper already represented through another provider is not counted twice.
7. Genus/family SI/SC statements are never promoted to strict species-direct evidence.

## Local preparation

Keep the request ZIP outside the repository.

```bash
python -m island_v2.try_traits prepare \
  --try-zip /private/path/51897_04092026085543.zip \
  --master-csv data/v2/staging/gbif/collected/island_taxa.csv \
  --output-dir /private/path/try_51897_prepared
```

Optional exact synonym mapping:

```bash
  --synonym-map-csv /private/path/try_synonyms.csv
```

The synonym file must contain `accepted_species` plus one of `synonym`, `source_name`, `matched_source_name`, or `original_name`.

Preparation writes:

- `trait_candidates.csv.gz` — source-backed strict candidates;
- `try_common_direct_evidence.csv.gz` — rows already converted to the shared `EVIDENCE_COLUMNS` schema with `source_group=try`;
- `try_corolla_fusion_sidecar.csv.gz` — separate corolla-fusion evidence;
- `try_name_audit.csv.gz` — matched/unmatched TRY names;
- `try_exclusions.csv.gz` — rejected and conflicting source rows with reasons;
- `try_integration_summary.json` — counts and input SHA-256 hashes.

## Rebuild affected strict cells

Use the prepared common evidence with the current formal lineage/direct/Validated-Low/coverage artifacts:

```bash
python -m island_v2.try_trait_integration \
  --try-common /private/path/try_51897_prepared/try_common_direct_evidence.csv.gz \
  --formal-lineage /path/to/current/reviewed_source_lineages.csv.gz \
  --current-direct /path/to/current/resolved_direct_species_trait.csv.gz \
  --current-low /path/to/current/rebuilt_all_evidence_validated_low.csv.gz \
  --current-coverage /path/to/current/strict_species_axis_coverage.csv.gz \
  --master-genus-map data/v2/staging/gbif/collected/island_taxa.csv \
  --ontology config/trait_ontology.yml \
  --output /private/path/try_51897_integrated
```

If the current direct ledger contains post-formal common-evidence packages that are not in `--formal-lineage`, pass each package with repeated `--additional-common` arguments. This prevents an affected genus/trait recomputation from accidentally dropping a newer independent source.

The integration recomputes only genus x trait cells touched by TRY, but it reruns both direct conflict resolution and the corresponding Validated-Low rules. Unaffected cells are copied through unchanged.

Outputs include the post-TRY direct lineage ledger, duplicate-lineage audit, resolved direct ledger, rebuilt Validated-Low rows, strict species-axis coverage, affected genus rules, conflict audit, and a before/after integration summary.

## Data handling

The request ZIP and its row-level prepared outputs should remain private/local unless redistribution is separately confirmed. Commit the adapter, tests, mapping contract, and summary statistics—not the raw TRY request itself. This also keeps the public repository reproducible without turning a private TRY request into a redistributed dataset.

## 2026-09-04 adapter sanity check

Using the request's own 11,700 accepted names as a temporary all-TRY master (not the island master) only to validate parsing/mapping, the adapter produced:

- flower colour: 10,531 species;
- floral symmetry: 5,433 species;
- strict SI/SC: 1,444 species;
- corolla fusion sidecar: 4,765 species.

These are **not island coverage gains**. The definitive gain is the before/after result from `try_trait_integration` after matching against the current island species master and resolving overlap/conflicts with the existing evidence ledger.
