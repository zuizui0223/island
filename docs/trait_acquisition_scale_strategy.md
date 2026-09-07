# Trait acquisition scale strategy

## Purpose

The Chapter 1 universe is 106,295 analysis-applicable species. Acquisition must therefore optimize **net unresolved species-axis cells recovered per source and per unit effort**, not the number of pages inspected or commits produced.

The current bottleneck is `reproductive_assurance`. Easy flower-colour gains must not displace reproductive acquisition simply because they are faster to inspect.

## 1. Source-scale acquisition is the default unit

Acquire a source as a whole whenever it contains a multi-species table, monograph, flora treatment, review appendix, supplementary dataset, thesis table, or other structured multi-taxon evidence.

Do not stop after the first few matching species. For each selected source:

1. inventory the full taxon scope;
2. extract all strict trait statements that can be mapped without cross-trait inference;
3. reconcile source names to the fixed analysis universe;
4. retain source lineage, page/table/row location and exact provenance;
5. only then compare the complete extracted packet with unresolved cells.

Species-by-species page hunting is reserved for conflict repair or exceptionally high-leverage missing evidence, not routine acquisition.

## 2. Rank sources before extraction

Every candidate source should be screened before detailed extraction. Record at least:

- source identity and type;
- estimated number of taxa represented;
- overlap with currently unresolved species, separately for reproductive assurance, structure and colour;
- expected number of strict extractable species-trait rows;
- taxonomic reconciliation burden;
- access/extraction format (CSV/table/text/OCR/manual);
- estimated extraction effort;
- expected net unresolved cells and expected cells per unit effort.

Prioritization is lexicographic rather than threshold gaming:

1. reproductive-assurance net-cell opportunity;
2. total expected strict net cells;
3. expected cells per unit effort;
4. lower reconciliation/provenance risk.

A source that is easy but mostly fills already-strong flower colour should rank below a comparably tractable source that materially advances reproductive assurance.

## 3. Evaluate existing Validated Low capacity in parallel

Additional acquisition and Low evaluation are separate lanes that run in parallel.

Use the complete direct-evidence ledger to re-evaluate all `genus × individual_trait` rules under the **existing current thresholds** only. Do not lower minimum species, dominance, species-LOO, source-lineage-LOO, or masked-accuracy requirements to create coverage.

For every eligible rule, report:

- trait-specific supporting species and independent source lineages;
- dominant state and counterexamples;
- dominance;
- species leave-one-out accuracy;
- source-lineage leave-one-out accuracy;
- number of unresolved species-trait cells that would be filled;
- conflicts or taxonomic-identity blockers.

Family inference and global fallback do not count as strict Chapter 1 acquisition.

## 4. Batch formal integration

Reviewed source packets may accumulate without re-running formal coverage after each edit.

Formal integration should occur only after a meaningful source-scale packet has been completed and reviewed. The integration run then validates the complete batch, rebuilds relevant direct cells and current-threshold Low rules, and reports net species-axis gain.

The `Integrate reviewed restart evidence` workflow is manual-only for this reason. Its dispatch `batch_label` must identify the completed source packet being integrated.

## 5. Reporting contract

For each formal batch report, in order:

1. source(s) fully processed and source-scale extraction completeness;
2. reviewed direct rows and unique species;
3. net newly resolved cells by axis;
4. reproductive-assurance gain first;
5. current-threshold Validated Low gains or invalidations;
6. conflicts/rejections;
7. remaining unresolved cells;
8. source ROI observations for selecting the next batch.

Do not report candidate rows, page hits, or machine discoveries as coverage gains.

## PR151 pivot checkpoint

Run `34107096317`, artifact `10012805875`, verifies the latest Agalmyla increment relative to Run `34106558478` as exactly **18 species-axis cells: 14 flower colour + 4 floral structure + 0 reproductive assurance**.

Across the recoverable public restart baseline (222,375 cells), the current integrated public artifact has 222,482 filled cells: **+107 total = +54 colour, +22 structure, +31 reproduction**. This public restart state is not a reconstruction of the historical private TRY + Wave55 checkpoint and must remain labelled separately.

The 18-cell Agalmyla increment is valid evidence, but its axis composition demonstrates why future acquisition must be selected at source scale with reproductive assurance as the primary bottleneck.
