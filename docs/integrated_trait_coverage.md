# Integrated three-axis trait coverage

`island-v2-integrated-trait-coverage` produces a coverage ledger on the fixed
`106,295 accepted species × 3 axes` denominator. It is intentionally separate
from the general fill cascade: family inference and global fallback are never
eligible here.

## Evidence precedence

Each `accepted_species × axis` key keeps one best record in this order:

1. species/synonym-direct High;
2. species/synonym-direct Medium;
3. Validated Low from the validated genus-consensus artifact.

The three-pass artifacts, every artifact-confirmed expanded wave, targeted
reacquisition, and Validated Low form the `before` ledger. The audited bulk and
promoted-public-web records form the integration increment.

## Bulk acceptance gate

Bulk rows may retain `review_status=pending`; that status is not treated as an
approval. A row is counted only when all deterministic source-contract checks
pass:

- the provider has an explicit contract in
  `integrated_trait_coverage.py`;
- `candidate_kind=source_backed`;
- evidence is species-direct or an exact synonym-to-accepted match;
- name matching is exact or author-stripped exact;
- the trait and controlled value map to one of the three axes;
- URL, source record ID, citation, and excerpt are present;
- no inference rule is present; and
- confidence is High or Medium.

eFloras uses its corresponding species-direct, source-backed extraction
contract. Promoted public-web records must already be source-backed High or
Medium rows in the all-master evidence ledger. Every accepted and rejected row
is retained in `bulk_candidate_acceptance_audit.csv.gz`.

## Source-lineage deduplication

Redistributing providers do not create independent evidence. Lineage keys
prefer an original DOI, then known origin aliases (for example, EOL records
whose citation identifies USDA PLANTS), provider component datasets, normalized
original citations, and finally canonical source URLs. Only one record per
`accepted_species × axis × source_lineage` is retained. Duplicate rows and
same-lineage value disagreements remain visible in
`source_lineage_duplicates.csv.gz`.

## Reproduction

The pinned run and artifact inventory is
`config/integrated_trait_coverage_sources.json`. The workflow
`build-integrated-trait-coverage.yml` resolves every artifact to its live
artifact ID, digest, expiry, run status, and run URL before downloading it. It
also verifies the expanded-wave summaries are non-overlapping and cover the
artifact-confirmed interval `[0,107000)`.

Run the workflow manually after every pinned source run is successful, or call
the CLI with the equivalent artifact directories:

```bash
island-v2-integrated-trait-coverage \
  --three-pass-high-dir /path/to/pass1 \
  --three-pass-medium-dir /path/to/pass2 \
  --three-pass-unresolved-dir /path/to/pass3 \
  --expanded-artifact-dir /path/to/expanded-run-1 \
  --expanded-artifact-dir /path/to/expanded-run-2 \
  --targeted-dir /path/to/targeted \
  --validated-low-dir /path/to/validated-low \
  --all-master-dir /path/to/all-master \
  --source-manifest-json config/integrated_trait_coverage_sources.json \
  --output-dir /path/to/output
```

The summary records before/after quality and axis counts, fill rates, exact
0/1/2/3-axis species counts, unresolved species-axis keys, Validated Low
upgrades, source-specific marginal gains, and per-run direct contributions.
The artifact also contains the selected coverage table, full accepted-lineage
ledger, unresolved tables, source-lineage duplicates, and the resolved source
manifest.
