# PR138 exact-SI genus-fixed + source-pool checkpoint

Status: **lineage decomposition / claim-boundary checkpoint; not a new primary estimand**.

Frozen-input deterministic run: 1,000 within-genus permutations, seed `20260827`. The dedicated Actions workflow is committed separately.

## Question

The exact-SI analysis already showed that the Palearctic attraction/accessibility shift does not require measured selfing capacity. This checkpoint asks the next stricter question: **does the SI-restricted shift exceed what is expected from the observed genus composition alone?**

For each of `large_bee_like` and `generalized_accessible`, syndrome scores are permuted only among scored exact-SI species within the same genus. Island species membership, floristic status, missing-score locations, GIFT mainland membership, source-region eligibility, outcome-blind source assignments, area, climate, distance and spatial blocks remain fixed. Island and mainland source scores are recomputed for every draw before source subtraction.

## Estimand reproduction

The all-analysis observed exact-SI Palearctic slopes reproduce the frozen PR138 restriction artifact exactly: maximum absolute coefficient difference `2.9e-16` (SE `2.1e-16`; p-value `6.3e-19`). The lineage null therefore attacks the existing SI estimand rather than replacing it with a different response.

## Main result: all-analysis-eligible evidence

### `all_native`

- observed source-adjusted SI slope: `+0.0714` to `+0.0788` across four frozen source definitions;
- genus-fixed null mean slope: `+0.0752` to `+0.0841`;
- observed-minus-null residual slope: `-0.0065` to `-0.0038`;
- empirical one-sided P(observed > genus-null): `0.739` to `0.838`;
- no positive residual is recovered.

### `native_nonendemic`

- observed source-adjusted SI slope: `+0.0643` to `+0.0687` across four frozen source definitions;
- genus-fixed null mean slope: `+0.0708` to `+0.0744`;
- observed-minus-null residual slope: `-0.0065` to `-0.0042`;
- empirical one-sided P(observed > genus-null): `0.775` to `0.845`;
- no positive residual is recovered.

The null is not merely wide enough to include the observed effect: its mean is slightly **more positive** than the observed SI slope. The source-adjusted Palearctic SI response therefore does not exceed the between-genus expectation.

## Strict direct-only sensitivity

Using both direct SI status and direct-only floral syndrome evidence, the observed slopes remain positive (`+0.0543` to `+0.0646`), but the genus-null means are essentially identical (`+0.0556` to `+0.0650`). Residual slopes are only `-0.0015` to `+0.0001`, all residual FDR q >= `0.964`, and empirical one-sided q >= `0.616`.

Thus the lineage result is not created by genus-consensus SI gap filling or by genus-consensus floral evidence.

## Evidence audit

- `all_analysis_eligible` / `large_bee_like`: 3,879 exact-SI species in the ledger; 827 axis-scored SI species across 310 genera; 101 genera / 618 species contribute within-genus permutations; 521 GIFT source regions remain eligible.
- `all_analysis_eligible` / `generalized_accessible`: 3,879 exact-SI species in the ledger; 1,765 axis-scored SI species across 276 genera; 97 genera / 1,586 species contribute within-genus permutations; 571 GIFT source regions remain eligible.
- `direct_only` / `large_bee_like`: 1,881 exact-SI species in the ledger; 468 axis-scored SI species across 280 genera; 73 genera / 261 species contribute within-genus permutations; 467 GIFT source regions remain eligible.
- `direct_only` / `generalized_accessible`: 1,881 exact-SI species in the ledger; 380 axis-scored SI species across 230 genera; 61 genera / 211 species contribute within-genus permutations; 458 GIFT source regions remain eligible.

## Updated scientific boundary

Supported:

1. Palearctic isolation is associated with an assemblage-level attraction/accessibility shift.
2. The shift is not a simple measured-selfing by-product: it remains positive in exact-SI floras.
3. **However, the SI-restricted source-adjusted slope is fully reproducible under a genus-composition-preserving null. No positive residual attraction/accessibility slope remains after genus structure is fixed.**

Therefore the current data do **not** support an interaction-mediated floral filter acting beyond measured genus composition. The strongest supported ordering is now `isolation -> biogeographic / lineage assembly -> floral syndrome composition`. A pollination filter could still act by filtering whole lineages, so disappearance under the genus null does not demonstrate that pollination is irrelevant; it only removes the claim of residual within-genus / beyond-genus filtering.

This result strengthens the Chapter 1 pattern-first boundary: PR138 identifies geographically contingent assemblage reorganization, while causal pollination-channel attribution remains unresolved and belongs to independent channel evidence / Chapter 2-style mechanism identification.
