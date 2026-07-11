"""One-shot migration that stops recovery workflows from blocking primary screening."""

from __future__ import annotations

import argparse
from pathlib import Path

PRIMARY = Path(".github/workflows/drive-global-trait-campaign.yml")
LATEST_RECOVERY = Path(
    ".github/workflows/recover-verified-multilingual-lead-full.yml"
)
OBSOLETE_RECOVERIES = (
    Path(".github/workflows/recover-multilingual-wikimedia.yml"),
    Path(".github/workflows/recover-verified-multilingual-wikimedia.yml"),
)

PRIMARY_OLD = "  group: global-trait-campaign\n  cancel-in-progress: false"
PRIMARY_NEW = "  group: global-trait-campaign-primary\n  cancel-in-progress: false"
RECOVERY_OLD = "  group: global-trait-campaign\n  cancel-in-progress: false"
RECOVERY_NEW = "  group: multilingual-recovery-v2\n  cancel-in-progress: true"


def transform_primary(text: str) -> str:
    if text.count(PRIMARY_OLD) != 1:
        raise ValueError("primary concurrency marker was not unique")
    return text.replace(PRIMARY_OLD, PRIMARY_NEW)


def transform_latest_recovery(text: str) -> str:
    push_start = text.find("  push:\n")
    permissions_start = text.find("\npermissions:\n", push_start)
    if push_start < 0 or permissions_start < 0:
        raise ValueError("latest recovery push block boundaries were not found")
    without_push = text[:push_start] + text[permissions_start + 1 :]
    if without_push.count(RECOVERY_OLD) != 1:
        raise ValueError("latest recovery concurrency marker was not unique")
    transformed = without_push.replace(RECOVERY_OLD, RECOVERY_NEW)
    if "  push:\n" in transformed:
        raise ValueError("latest recovery still has an automatic push trigger")
    return transformed


def validate_inputs() -> tuple[str, str]:
    primary = PRIMARY.read_text(encoding="utf-8")
    latest = LATEST_RECOVERY.read_text(encoding="utf-8")
    for path in OBSOLETE_RECOVERIES:
        if not path.exists():
            raise ValueError(f"obsolete recovery workflow already missing: {path}")
        if "workflow_run:" not in path.read_text(encoding="utf-8"):
            raise ValueError(f"obsolete recovery lacks expected workflow_run: {path}")
    return transform_primary(primary), transform_latest_recovery(latest)


def apply() -> None:
    primary, latest = validate_inputs()
    PRIMARY.write_text(primary, encoding="utf-8")
    LATEST_RECOVERY.write_text(latest, encoding="utf-8")
    for path in OBSOLETE_RECOVERIES:
        path.unlink()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    primary, latest = validate_inputs()
    if args.apply:
        PRIMARY.write_text(primary, encoding="utf-8")
        LATEST_RECOVERY.write_text(latest, encoding="utf-8")
        for path in OBSOLETE_RECOVERIES:
            path.unlink()
    else:
        assert "group: global-trait-campaign-primary" in primary
        assert "group: multilingual-recovery-v2" in latest
        assert "  push:\n" not in latest


if __name__ == "__main__":
    main()
