# Source-scale trait acquisition and probabilistic genus Low (2026-08-09)

This wave starts from formal artifact Run `31290142917`. Its fixed universe is
106,295 species and 318,885 species-axis cells; 119,571 cells are resolved and
72,570 species-trait cells have direct High or Medium evidence. Rows already in
that artifact are not reacquired merely to inflate the package.

The immediate direct-evidence addition is the published seed-plant mating-system
database archived at Dryad (`10.5061/dryad.577q1`; Moeller et al. 2017). The
source has 624 study rows for 492 raw taxa. Exact binomial matching to the fixed
master and the formal novelty baseline yields 206 retained source rows for 161
accepted species and 173 previously absent species-trait pairs:

- `mating_system`: the continuous multilocus outcrossing rate (`mean.tm`) is
  retained in the quote and mapped to predominantly selfing (`<0.2`), mixed
  (`0.2-0.8`), or predominantly outcrossing (`>0.8`);
- `self_incompatibility`: the source's explicit `si` field is retained as SI or
  SC and is never substituted for mating system or autonomous selfing.

Each study row keeps its underlying publication citation as `source_lineage`.
Multiple studies or values for one species remain separate evidence rows so the
common conflict resolver can distinguish agreement, biological variation, and
unresolved conflict. Non-exact names, non-species ranks, cultivar/hybrid text,
and taxa outside the fixed master fail closed.

The complete incremental package contains 208 evidence rows and 175
species-trait pairs, all absent from the exact novelty baseline reconstructed
from the integrated baseline evidence plus the cumulative public-Web ledger.
Two additional exact records from PlantNET and Micheneau et al. remain in the
audited package but below the ten-row per-trait production gate. The deterministic
record audit reviews all 208 rows and the common promotion gate selects the 206
Moeller rows in the mating-system and self-incompatibility strata. The package
summary records every input and output hash.

## Separate probabilistic Low artifact

The strict Validated Low contract is unchanged: the repository's shared
trait-specific implementation joins `genus × trait_name`, requires at least
three direct species, and applies the current trait-specific dominance and
masked-validation thresholds. It never joins on `genus × axis` alone.

`analysis/build_secondary_probabilistic_genus_low_20260809.py` materializes two
additional, non-strict analysis inputs from the same audited rules:

- relaxed-min3 secondary predictions: 2,498 eligible rules, 13,493 additional
  species-trait predictions, and 5,781 currently unresolved species-axis cells
  potentially covered. Species leave-one-out accuracy is 0.9538 and
  source-lineage leave-one-out accuracy is 0.9494;
- min-species-2 diagnostic predictions: 20,161 potentially covered cells, but
  these are sensitivity-only and are not recommended for the main analysis.

Both files use `quality=probabilistic_low`, set
`strict_coverage_eligible=false`, preserve the empirical value distribution,
and require one species-level draw to be shared across all islands within an
imputation. Family inference, global fallback, and axis-only joins are absent.
The strict 119,571-cell coverage is not changed by these artifacts.

## Trait and provenance guardrails

Flower colours remain state sets. `pollen_vector_mode` and `reward_type` remain
independent traits and are not mapped to the strict floral-structure axis.
Self-incompatibility, autonomous selfing, mating system, and cleistogamy remain
separate traits. Search snippets are not evidence; every direct row stores the
source URL, exact record or quote, retrieval metadata, match method, and source
lineage.

The workflow chains from `reviewed-open-web-evidence-31290142917`, promotes the
reviewed direct package, rebuilds the strict all-evidence ledger and Validated
Low using the shared implementation, then emits the secondary and sensitivity
genus-prediction artifacts. GitHub's final artifact is authoritative; local
figures are preflight checks only.
