from __future__ import annotations

import tomllib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LEGACY = ROOT / "legacy" / "chapter1-pre-v13"

CANONICAL_V13 = (
    "README.md",
    "docs/PAPER_PIPELINE.md",
    "docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md",
    "config/chapter1_v13_unified_island_syndrome_result_lock.json",
    "config/chapter1_v13_functional_bridge_result_lock.json",
    "docs/chapter1_v13_submission_figure_sync_20260917.md",
    "docs/chapter1_unified_hypothesis_20260917.md",
    "docs/chapter1_v13_raw_colour_audit_20260917.md",
    "docs/chapter1_v13_raw_colour_coupling_audit_20260917.md",
    "docs/chapter1_submission_freeze_v13_20260917.md",
)

# Sentinels for retired publication/mechanism branches.  The archive retains the
# exact old relative path so provenance is discoverable without polluting the
# active v13 surface.
RETIRED_SENTINELS = (
    "src/island_v2/bombus_regime.py",
    "src/island_v2/chapter1_h5_bombus_upstream_bridge.py",
    "src/island_v2/chapter1_pr138_palearctic_restricted_ipw.py",
    "docs/chapter1_submission_freeze_20260915.md",
    "docs/manuscript_submission_contract.md",
    "config/chapter1_h5_bombus_upstream_bridge_v1.yml",
    "tests/test_chapter1_h5_bombus_upstream_bridge.py",
    ".github/workflows/run-chapter1-h5-bombus-upstream-bridge.yml",
)

RETIRED_CLI_TARGET_TOKENS = (
    "island_v2.bombus_",
    "island_v2.trait_bombus",
    "island_v2.m1_m3_bombus",
    "island_v2.chapter1_h5_bombus",
    "island_v2.chapter1_pr138_palearctic",
    "island_v2.chapter1_submission_freeze",
)


def test_canonical_v13_surface_stays_active() -> None:
    missing = [path for path in CANONICAL_V13 if not (ROOT / path).exists()]
    assert not missing, f"canonical v13 files missing from active surface: {missing}"


def test_retired_storyline_files_live_only_in_legacy() -> None:
    still_active = [path for path in RETIRED_SENTINELS if (ROOT / path).exists()]
    missing_archive = [path for path in RETIRED_SENTINELS if not (LEGACY / path).exists()]
    assert not still_active, f"retired pre-v13 files still active: {still_active}"
    assert not missing_archive, f"retired files missing from archive: {missing_archive}"


def test_active_cli_surface_does_not_expose_retired_storylines() -> None:
    data = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    scripts = data["project"]["scripts"]
    offenders = {
        name: target
        for name, target in scripts.items()
        if any(token in str(target) for token in RETIRED_CLI_TARGET_TOKENS)
    }
    assert not offenders, f"retired CLI commands remain active: {offenders}"


def test_active_handoff_points_to_archive_boundary() -> None:
    for path in (ROOT / "README.md", ROOT / "docs" / "PAPER_PIPELINE.md"):
        text = path.read_text(encoding="utf-8")
        assert "legacy/chapter1-pre-v13/" in text
