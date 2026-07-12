# Fill-first trait cascade

## Purpose

The operational priority is **trait yield**, not measurement. Direct
source-backed evidence covers well under 1% of the master for most fields, so
`island-v2-trait-fill-cascade` fills every eligible species-trait by descending a
source-blind fallback ladder until the 9-column `unknown` is driven toward zero:

```text
species_direct -> synonym_direct -> genus_inference -> family_inference -> global_fallback
```

## Angiosperm eligibility gate (before the cascade)

The target is 100% fill **within the eligible flowering-plant universe**, not
the raw Tracheophyta acquisition universe. Non-seed vascular plants (ferns,
lycophytes), gymnosperms, and fossil lineages have no flowers or comparable
pollen-vector states, so a family-based angiosperm gate runs *before* the
cascade (`config/angiosperm_scope.yml`, `island-v2-angiosperm-scope`). Out-of-scope
species are marked in `angiosperm_scope_by_species.csv` and are **never** filled
with floral or pollination priors. The fill denominator is the angiosperm-
eligible set. On the current master that is 106,295 of 115,328 species; the
9,033 out-of-scope split into non_seed_vascular, gymnosperm, fossil, and
family_unresolved (missing family). Two fern-sounding families that are actually
angiosperms — Cardiopteridaceae, Pteridophyllaceae — stay eligible.

Low-resolution fills are never discarded and never hidden. Every fill carries
`fill_tier`, `evidence_scope`, and `confidence`, so a genus/family/global fill is
always separable for sensitivity analysis and never masquerades as accepted
direct evidence. Fills are written to staging, never to curated. A cascade fill
is candidate coverage, never a biological absence.

## Guardrails

- **Traits are filled independently.** One trait is never derived from another;
  self-compatibility never fills autonomous selfing.
- **Inference keeps its distribution.** Genus/family/global fills record the
  supporting `value_distribution`, so a modal fill can be checked against a
  distribution-aware draw rather than trusted as a hard code.
- **Support thresholds.** `min_genus_support` / `min_family_support` gate the
  inference tiers; below family support the global fallback fills.
- **Direct evidence wins.** A species with its own value always keeps it. Curated
  evidence outweighs machine candidates in the modal vote.

## Run it

```bash
island-v2-trait-fill-cascade run --output-dir data/v2/staging/traits/fill_cascade
```

Outputs:

- `fill_coverage_summary.json` — eligible denominator, out-of-scope counts by
  group, and per-trait fill rate, species-direct count, unknown remaining, and
  tier composition.
- `benchmark_sample.csv` — the deterministic 100-species benchmark (eligible
  species), one column per trait, each cell `value [fill_tier]`.
- `angiosperm_scope_by_species.csv` — every master species with its
  `taxonomic_group` and `angiosperm_analysis_eligible` flag (regenerable,
  git-ignored).
- `trait_fills.csv.gz` — the full long fill table over eligible species
  (regenerable, git-ignored).

## Sources

All evidence channels are unioned before the cascade, source-blind: the machine
wave index, search-enabled LLM batches, bulk trait joins (`candidate_long`
layout), and curated accepted evidence. A newly landed bulk source raises
species-direct coverage and sharpens every inference tier with no code change.

## Reading the current baseline

With direct evidence alone at 0.02–1.45% per field, the cascade takes every
field to 100% fill (zero unknown) across the 106,295 angiosperm-eligible
species — dominated by genus/family/global inference, with only ~3,700
species-direct cells. That is the intended yield-first state: maximal coverage
within the flowering-plant universe now, resolution made explicit so quality is
a downstream sensitivity axis, not an acquisition gate. The lever that upgrades
resolution is landing more direct evidence (bulk sources, un-gated colour
extraction); rerun the cascade after each landing and watch species-direct climb
and the global-fallback share fall.
