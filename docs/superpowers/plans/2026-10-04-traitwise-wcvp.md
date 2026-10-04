# Traitwise H1 and WCVP-only origin sensitivity — 2026-10-04

User-approved post-result redesign. Remove the H1 directional composite; do not construct three domain scores. Retain seven previously defined binary outcomes, model each in four regions, with corrected isolation, area and climate controls. Broad flora remains primary; WCVP regional-native-compatible flora is the only active origin sensitivity. Direct-only evidence remains an evidence-quality sensitivity. Preserve H2–H4 and historical analyses.

1. Test and implement separate beta-binomial isolation slopes with cluster-robust t inference (G-1 df), pointwise 95% CIs and two-sided Holm adjustment over the fixed 28 region-trait tests per flora/evidence scope. Holm is valid without assumptions about correlation among traits; failed tests retain family size with p=1 and are marked untestable.
2. Rerun four scopes using frozen input artifacts, corrected geography and WCVP mapping; record input hashes, support, convergence, and all estimates. No selection by significance.
3. Publish a traitwise results table and clear scope/provenance report. Change active H1 selector and current documentation only after successful replay. Archive old aggregate results unchanged.
4. Run focused regression tests and a final independent review; report results and any remaining constraints.

Ruling: use Holm rather than BH for 28 correlated tests to control family-wise error without a dependency assumption. Domain labels organize interpretation only. This redesign is retrospective, not a new preregistered test.

## Completion ledger

- 112/112 separate fits converged from pinned artifacts; WCVP coverage matched 513,320 records / 2,372 islands.
- New tests first failed on missing module; implementation then passed; final focused suite: 21 passed.
- Independent reviewer confirmed numerical inference/receipts and found two documentation defects, both fixed.
- Active H1 selector, README, manuscript, SI and figure captions switched; old main/SI/captions saved verbatim under submission/chapter1_current/history/pre_traitwise_20261004.
- Administrative package files carry explicit supersession notices; no claim of submission readiness.
- H2–H4 original result files remain unchanged. Poster rendering is outside this bounded reanalysis.
