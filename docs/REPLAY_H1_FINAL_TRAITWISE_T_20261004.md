# Replay traitwise H1 — 2026-10-04

Run from the repository root with Python >=3.11 and project dependencies installed (`pip install -e ".[dev]"`). Authentication must permit reading the pinned GitHub Actions artifacts. These artifacts have retention limits; archive the exact inputs for long-term replay. SHA256 values are recorded in `results/h1_final_traitwise_t_20261004/manifest.json`.

```powershell
New-Item -ItemType Directory -Force work/traitwise/source, work/traitwise/status, work/traitwise/wcvp | Out-Null
gh run download 34232450884 --repo zuizui0223/island -n chapter1-progressive-analysis-34232450884 -D work/traitwise/source
gh run download 32559322028 --repo zuizui0223/island -n official-wcvp-corroborated-status-32559322028 -D work/traitwise/status
gh run download 37093097826 --repo zuizui0223/island -n chapter1-wcvp-native-corrected-37093097826 -D work/traitwise/wcvp
$env:PYTHONPATH='src'
$env:OPENBLAS_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
python -m island_v2.chapter1_h1_traitwise --flora work/traitwise/source/fixed/canonical/input/chapter1_status_flora.csv.gz --all-audit work/traitwise/source/snapshot/chapter1_trait_state_audit_all.csv.gz --direct-audit work/traitwise/source/snapshot/chapter1_trait_state_audit_direct.csv.gz --covariates results/geography_20260924/corrected_geography_covariates.csv --wcvp work/traitwise/wcvp/ch1-wcvp/ranges/wcvp_native_range_summary.csv.gz --mapping work/traitwise/status/resolved/island_tdwg_l3_mapping.csv --config config/chapter1_h1_final_traitwise_t_20261004.yml --output work/traitwise/replay
python -m pytest tests/test_chapter1_h1_traitwise.py tests/test_chapter1_wcvp_native_compatibility.py tests/test_chapter1_submission_surface.py -q
```

The CLI replays broad + WCVP flora for all-analysis + Direct-only evidence. Expect 112 rows, 28 per scope. Compare input hashes before comparing estimates; do not substitute newer WCVP ranges silently. The committed replay used Python 3.12.14, NumPy 2.3.5, SciPy 1.18.1. Numerical comparisons across platforms should allow optimizer tolerance rather than require byte identity.

`python scripts/render_h1_final_traitwise_t.py` renders the committed results to PDF/SVG/PNG. Every symbol is a fitted coefficient; none are illustrative data. The figure uses pointwise intervals and independently marks nominal unadjusted p < .05 coefficients.

H2–H4 files are unchanged. H2 functional scores remain conditioning variables and H4 matching targets; this H1 request does not remove those distinct definitions. Historical directional-score modules remain for replay but are not the active H1 entry point.
