# All-evidence trait audit

`island-v2-all-evidence-trait-audit` is the strict, reproducible coverage and
genus-inference audit on 106,295 accepted angiosperms and three reporting axes.
It consumes the integrated artifact plus the latest pinned public-web artifact.
It never counts family inference or global fallback.

## Evidence and conflict contracts

Direct evidence is keyed by `accepted_species × trait_name`. Original-source
redistributions are collapsed at
`accepted_species × trait_name × source_lineage`; provider identity alone does
not create independent support. Direct evidence is resolved in High then Medium
order. Conflicting equal-quality evidence is classified as:

- independent-source agreement;
- true multistate or variable evidence;
- cultivar contamination;
- source ontology mismatch; or
- unresolved direct conflict.

Unresolved conflicts are excluded. Flower colour retains a state set and can be
`multicolored_variable`. Reproductive traits such as self-incompatibility,
autonomous selfing and mating system remain separate traits and are not used as
proxies for each other.

## Rebuilt Validated Low

Genus rules are rebuilt from every resolved direct High/Medium cell. Both rule
construction and application use `genus × trait_name`; the broad axis is
reporting metadata only. The primary Low setting requires at least three direct
species, axis-specific dominance, leave-one-species-out accuracy and
leave-one-source-lineage-out accuracy. The relaxed `min_species=3` setting is a
secondary robustness candidate. Both `min_species=2` settings are diagnostic
only.

## Analysis tracks

- Confirmatory uses source-lineage-deduplicated species/synonym-direct evidence,
  excludes unresolved direct conflicts, and retains multistate composition.
- Secondary robustness adds validated probabilistic genus inference. A species'
  imputed value is shared across all islands within one imputation, and
  uncertainty must be propagated by multiple imputation or measurement error.
- Sensitivity only may inspect family inference and global fallback. These
  values cannot supply headline effects, significance or posterior
  probabilities.

Missing evidence is never converted to zero or biological absence. The
reproductive-assurance pathway uses one common island set across its equations,
and analysis must address spatial dependence and source-pool or lineage
composition.

## Artifact contents

The workflow uploads the base integrated files and:

- `all_evidence_trait_coverage_summary.json`;
- `before_after_coverage_summary.json`;
- `strict_species_axis_coverage.csv.gz`;
- `rebuilt_all_evidence_validated_low.csv.gz`;
- `trait_specific_genus_rule_audit.csv.gz`;
- `unresolved_reason_by_species_axis.csv.gz`;
- `threshold_sensitivity_report.json` and its summary CSV;
- `source_lineage_conflicts.csv.gz` and trait-specific duplicate table;
- `prioritized_acquisition_queue.csv.gz`;
- confirmatory and secondary inputs plus the three-tier analysis manifest;
- per-island endpoint coverage and the common reproductive-pathway island set;
  and
- `all_evidence_source_run_manifest.json`, including hashes of every local
  input.

The acquisition queue ranks genus-trait cells with two agreeing direct species,
many unresolved congeners, high island occurrence and record counts, and gives
a recommended source/query. It does not initiate a new web scan.
