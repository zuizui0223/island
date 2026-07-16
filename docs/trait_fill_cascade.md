# Evidence-grounded trait cascade

## Purpose

The operational priority is a **complete rectangular form without inventing
biological values**, while direct evidence keeps increasing. Every eligible
species-trait has a row, but `island-v2-trait-fill-cascade` assigns a value only
through this evidence ladder:

```text
species_direct -> synonym_direct -> qualified genus_inference
               -> qualified family_inference -> unresolved_no_evidence
```

## Angiosperm eligibility gate (before the cascade)

The target is 100% cell presence **within the eligible flowering-plant universe**, not
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

Low-resolution fills are never hidden. Every row carries `fill_tier`,
`evidence_scope`, and `confidence`. Genus/family values also retain support and
consensus diagnostics. `unresolved_no_evidence` carries `filled_value=unresolved`,
blank reported/inferred channels, and `analysis_eligible=false`. It is cell
accounting, not a biological state or absence.

## Guardrails

- **Traits are filled independently.** One trait is never derived from another;
  self-compatibility never fills autonomous selfing.
- **Inference keeps its basis.** Genus/family fills record `support_n` for the
  winning value, `total_support_n`, `supporting_taxa`, `supporting_genera_n`,
  `supporting_genera`, `value_distribution`, `winner_share`, and
  `inference_basis`. A tied mode is never randomly resolved. Family inference
  additionally requires support across the configured minimum number of genera.
- **Direct conflicts are not silently decided.** Equal-weight species-level
  states are removed from `reported_value`, retained in
  `direct_conflict_distribution`, and rescued only through the explicitly
  inferred channel. Multiple reported flower colours are retained as
  `multicolored_variable`, because colour is intentionally a plural trait.
- **One species, one taxonomic vote.** After resolving each species' direct value,
  genus and family distributions count each species once. Duplicate
  databases or repeated source rows cannot dominate inference merely by volume.
- **Ontology is fail-closed.** Source-specific direct-evidence adapters require
  species-direct scope and source provenance. Values are normalized to
  `config/trait_ontology.yml`; an out-of-domain value fails materialization
  instead of spreading through a fallback tier.
- **Support and consensus thresholds.** `min_genus_support` /
  `min_family_support` and `min_genus_consensus` / `min_family_consensus` gate
  inference. Below either gate the next taxonomic tier is tried, then the cell
  remains unresolved.
- **Direct evidence wins.** A species with its own value always keeps it. Curated
  evidence outweighs machine candidates in the modal vote.
- **Reported and inferred values remain distinguishable.** `species_direct` and
  `synonym_direct` are reported-evidence tiers; genus and family tiers are
  explicitly inferred sensitivity values. Every shard also has mutually
  exclusive `reported_value` and `inferred_value` columns alongside the selected
  `filled_value`; consumers do not have to infer origin from the value itself. A
  model-only answer is not loaded as direct evidence and can never train the
  reported genus/family summaries.
- **No guild completion requirement.** Pollinator functional guild is not a
  target of the exhaustive cascade. The target is the flower itself: colour,
  form, symmetry, tube depth, size, display and reward, plus pollen-vector and
  reproductive architecture/assurance traits.

## Sharded, checkpointed, resumable

The cascade does not recompute all ~115k species every run. Eligible species are
partitioned into deterministic shards by a stable hash of the accepted name, and
each shard is checkpointed independently. The shared model (species/genus/family
distributions) is built once from all evidence and applied per shard, so a
sharded run is bit-for-bit identical to a single whole-universe pass.

A shard is recomputed only when its inputs change:

- **model fingerprint** — any evidence row or cascade parameter change flips it
  and correctly invalidates every shard;
- **species fingerprint** — the shard's exact species membership, so master
  growth only touches the shards whose membership changed.

Repeat runs with unchanged inputs process **zero** shards and just re-aggregate
the cached per-shard summaries. `--max-shards` bounds work per invocation for a
scheduled 30-minute lane; rerun to resume and finalize.

```bash
# full build (128 shards), resuming any already-current shards
island-v2-trait-fill-cascade run --output-dir data/v2/staging/traits/fill_cascade

# bounded lane: advance at most 32 stale shards, then rerun to continue
island-v2-trait-fill-cascade run --output-dir DIR --max-shards 32

# how many shards are current for the present inputs, without filling
island-v2-trait-fill-cascade status --output-dir DIR
```

Outputs:

- `fill_coverage_summary.json` — eligible denominator, out-of-scope counts by
  group, per-trait fill rate/species-direct/unknown/tier composition, and shard
  completion (`shards_complete`, `n_shards_current`). Aggregated from the compact
  per-shard summaries, which are disjoint so the sums are exact.
- `cascade_manifest.json` — contract version, shard count, model fingerprint, and
  completion state.
- `benchmark_sample.csv` — the deterministic 100-species benchmark (eligible
  species), one column per trait, each cell `value [fill_tier]`. Reads only the
  shards holding the sample species.
- `shards/shard_XXXXX/` — per-shard `fills.csv.gz` + `shard_checkpoint.json`
  (regenerable, git-ignored).
- `angiosperm_scope_by_species.csv` — every master species with its
  `taxonomic_group` and `angiosperm_analysis_eligible` flag (regenerable,
  git-ignored).

## Sources

All reported-evidence channels are unioned before the cascade: the source-backed
machine-wave index, evidence-locked LLM extraction with exact source quotes,
bulk trait joins (`candidate_long` layout), and curated accepted evidence.
Model-only LLM output is deliberately excluded from this direct-evidence union.
A newly landed reported source raises
species-direct coverage and sharpens every inference tier with no code change.

The materialization workflow first restores the newest durable all-master
artifact, including its reported-evidence ledger and generated eFloras, EOL,
GIFT, USDA and AusTraits candidate files. It reruns after a successful all-master
build. The all-master workflow explicitly dispatches materialization with its
exact run id after artifact upload, avoiding workflow-chain depth limits. Bulk
and eFloras refresh completions also trigger a new all-master build, so their
direct evidence does not wait for a later public-web run. This prevents evidence acquired in Actions from disappearing merely
because the generated source files are not committed to Git. Workflow-triggered
runs use the triggering default-branch artifact; other runs resolve only a
successful default-branch artifact, and record that run id in the cascade
manifest. Generated sources are read from a clean temporary overlay, while the
tracked machine index and quote-locked LLM bundles remain separate inputs.
Pull-request validation may use the latest successful artifact from that exact
head branch; merged materialization is always rebuilt from the triggering main
all-master run.

## Reading the current baseline

The cascade retains 100% of the 106,295 x target-trait cells. A cell is resolved
only by species-direct evidence or a qualified genus/family distribution; every
other cell is explicit unresolved and excluded from analysis. Each new bulk
source or quote-locked extraction upgrades cells to `species_direct`, can rescue
related unresolved cells through auditable taxonomic support, and never
overwrites direct evidence.
