# Corrected submission replay

The current baseline is `config/chapter1_submission_current.json`. Do not use old v14 locks as the active result selection. This directory reuses the original analysis modules: the correction changes geographic exposure and removes the verified continental fragment, not trait values or model specifications.

Use Python >=3.11 and install the repository (`pip install -e .`). From the repository root:

```powershell
$env:PYTHONPATH='src'
python scripts/geography_correction/replay.py --sources C:/path/to/extracted-sources --output C:/path/to/new-replay
```

`--verify-only` checks hashes without running models. The output must be new or empty. The runner checks every external file in `config/chapter1_corrected_input_files_20260924.json` and all committed corrected-result hashes, reproduces old H2/H3/exact-H4/raw-pattern results, then runs corrected H1 (all three strata), H2, H3 sensitivities, exact/atomic H4 and raw patterns in both evidence scopes. Numerical platform differences are expected at floating-point tolerance. A completed stage receipt is not a claim that every hypothesis is supported.

## Input sources and layout

Download immutable artifacts from `zuizui0223/island` and extract under these directory names. Artifact SHA256 provenance is also retained in the parent v14 lock; **per-file checks in the corrected input manifest are mandatory**, even if the archive was downloaded earlier.

| Directory | Artifact ID | Workflow run |
|---|---:|---:|
| `v14-artifact` | 10535020072 | 35314955780 |
| `progressive-input` | 10058653212 | 34232450884 |
| `all-data-primary` | 10447959173 | 35100991898 |
| `glopl-global` | 10444156159 | 35090599662 |
| `h4-ra` | 10444749163 | 35093274622 |
| `h4-arch` | 10445189257 | 35094588521 |
| `raw-coupling` | 10495522252 | 35216672430 |

The manifest's `island--*.yml` files are corresponding `config/*.yml` files from parent commit `8f36b62fae99fe5bbb227768ca330b0bfd267e78`, copied with the prefix. `chapter1_v13_functional_bridge_v1.yml` retains its unprefixed name. Retrieve exact bytes from that commit; do not edit numerical settings to pass checks. If GitHub artifact retention removes an input, replay must stop until the exact matching archive is available; no newer dataset is substituted.

## Geometry reproduction

The published distances were calculated from the original GSHHG archive with SHA256 `8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41`. Obtain it from the [official GMT GSHHG 2.3.7 release](https://github.com/GenericMappingTools/gshhg-gmt/releases/download/2.3.7/gshhg-shp-2.3.7.zip). The locked island GPKG is in artifact `8066083419`, run `28659417688`, path `islands/gshhg/prepared/islands_v2.gpkg` (archive SHA256 `bee33e14672ec7ff4ed1f7acaea36b32cc6170a45c3d14b2caeb5f43bfaa623b`).

```powershell
python scripts/geography_correction/recompute_geometry.py --archive C:/sources/gshhg-shp-2.3.7.zip --islands C:/sources/islands_v2.gpkg --sites results/geography_20260924/glopl_corrected_site_distances.csv --output C:/new-geometry
```

This can be slow. It recomputes minimum great-circle arc distance on a mean-radius sphere, not an ellipsoid, and verifies the excluded component's exact geometry. The default model replay uses the hash-locked committed distance tables. Recomputed geometry should be compared by island/site identity with numerical tolerance before changing the submission contract. Do not overwrite the lock automatically.

## Interpretation and scope

This is an authorized post-hoc measurement correction selected as primary submission baseline. It is not new prospective confirmation. Corrected all-analysis H1 joint tests remain supported in four regions, but southern shallow/open tube effects are negative. Tropical direct-only H2 accessibility is not FDR-supported; H3 supplemental-only and H4 supplemental-only selfing remain weak. Read the full tables and `docs/chapter1_corrected_submission_20260924.md`.
