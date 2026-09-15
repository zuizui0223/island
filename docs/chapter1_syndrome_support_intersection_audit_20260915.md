# Chapter 1 syndrome support-intersection audit — 2026-09-15

Status: **completed diagnostic; minimum informative-trait requirement is not the explanation**.

## Provenance

- workflow run: `34964326934`
- artifact: `chapter1-syndrome-support-audit-34964326934`
- artifact ID: `10394368256`
- digest: `sha256:131529ce68a8edeba4287087a21702d94a6d1b181d745a4fad4ff1f50cbc0901`

## Question

The historical two-axis route first builds species-level syndrome-concordance scores and, in the frozen definition, requires at least two informative traits per species. Does that support intersection explain why the all-observed two-axis North--Tropical contrast is weak while the atomic probability route is supported?

The audit reran the same historical syndrome definitions twice, changing only:

- frozen: `minimum_informative_traits = 2`;
- diagnostic: `minimum_informative_traits = 1`.

Everything else was held fixed: trait weights, evidence scope, island covariates, contexts, clustering, branching model and multiple-testing structure.

## Result

| evidence scope | min informative traits | islands | blocks | p | q | supported |
|---|---:|---:|---:|---:|---:|---|
| all-analysis | 2 | 3,253 | 184 | 0.73998 | 0.88797 | no |
| all-analysis | 1 | 3,620 | 187 | 0.50197 | 0.50197 | no |
| direct-only | 2 | 3,173 | 183 | 0.06142 | 0.13094 | no |
| direct-only | 1 | 3,563 | 187 | 0.85119 | 0.85119 | no |

Lowering the requirement to one informative trait recovers nearly the full atomic island support, but **does not restore the two-axis North--Tropical contrast**. In direct evidence it becomes weaker rather than stronger.

## Interpretation

The weak historical two-axis route is therefore not caused simply by losing species or islands through the `minimum_informative_traits >= 2` rule.

Combined with the separate response-compression audit, the evidence now localizes the discrepancy further:

1. six atomic beta-binomial vector: supported;
2. matched grouped-binomial vector: supported;
3. linear projection of those same atomic slopes to two family means: supported;
4. species-level syndrome-concordance route: unsupported;
5. relaxing the species-level minimum informative-trait intersection: still unsupported.

Thus the remaining difference lies in the **species-first concordance/averaging construction itself** (including how heterogeneous trait availability is normalized within species and how species then receive equal weight in island means), not merely in dimensionality, model family, or the two-trait support gate.

## Next gate

Use a trait-first family summary or an equivalent constrained atomic model that preserves trait-specific denominators, then compare it directly against the species-first concordance score on the same island frame. Do not change the biological trait definitions or inspect new signs when choosing the contrast.
