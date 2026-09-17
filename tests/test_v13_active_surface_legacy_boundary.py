from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "legacy" / "chapter1-pre-v13"

# These branches are historical scientific provenance after the v13 (#229-#234)
# canonicalization. They may exist under legacy/, but must not re-enter the active
# package, test, workflow, configuration, or publication-facing document surface.
RETIRED_ACTIVE_PATTERNS = (
    "src/island_v2/bombus_*.py",
    "src/island_v2/trait_bombus_analysis.py",
    "src/island_v2/m1_m3_bombus_channel_input.py",
    "src/island_v2/chapter1_h5_bombus_*.py",
    "src/island_v2/chapter1_pr138_*.py",
    "tests/test_bombus*.py",
    "tests/test_trait_bombus_analysis.py",
    "tests/test_m1_m3_bombus_channel_input.py",
    "tests/test_chapter1_h5_bombus_*.py",
    "tests/test_chapter1_pr138_*.py",
    ".github/workflows/*bombus*.yml",
    ".github/workflows/*pr138*.yml",
    "config/chapter1_h5_bombus_*",
    "docs/chapter1_h5_bombus_*.md",
    "docs/chapter1_pr138_*.md",
    "docs/chapter1_submission_freeze_20260914.md",
    "docs/chapter1_submission_freeze_20260915.md",
    "docs/chapter1_p0_*.md",
    "docs/chapter1_p1_*.md",
    "docs/manuscript_submission_contract.md",
    "scripts/build_chapter1_v10_*.py",
    "scripts/promote_chapter1_v11_*.py",
)

CANONICAL_V13 = (
    "docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md",
    "config/chapter1_v13_unified_island_syndrome_result_lock.json",
    "config/chapter1_v13_functional_bridge_result_lock.json",
    "docs/chapter1_unified_hypothesis_20260917.md",
    "docs/chapter1_v13_raw_colour_audit_20260917.md",
    "docs/chapter1_v13_raw_colour_coupling_audit_20260917.md",
    "docs/chapter1_submission_freeze_v13_20260917.md",
)

# Frozen parent result locks are deliberately retained in active config because the
# v13 result lock validates their contracts directly. Historical implementation code
# does not need to remain importable just because the frozen parent result does.
PROTECTED_PARENT_LOCKS = (
    "config/chapter1_all_data_route_result_lock.json",
    "config/chapter1_v12_two_panel_result_lock.json",
    "config/chapter1_v12_h5_glopl_extension_result_lock.json",
)


def test_retired_pre_v13_branches_are_not_on_active_surface() -> None:
    offenders: list[str] = []
    for pattern in RETIRED_ACTIVE_PATTERNS:
        offenders.extend(str(path.relative_to(ROOT)) for path in ROOT.glob(pattern))
    assert not sorted(set(offenders)), "retired pre-v13 files remain active: " + ", ".join(
        sorted(set(offenders))
    )


def test_archive_manifest_exists() -> None:
    assert (ARCHIVE / "ARCHIVE_MANIFEST.md").is_file()


def test_canonical_v13_and_required_parent_locks_remain_active() -> None:
    for relative in (*CANONICAL_V13, *PROTECTED_PARENT_LOCKS):
        assert (ROOT / relative).is_file(), relative
