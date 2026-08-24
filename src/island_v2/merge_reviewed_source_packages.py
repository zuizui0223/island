"""Merge independently audited trait source packages deterministically."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

from island_v2.open_web_common import reviewed_source_package_evidence


def _sha256(path: Path) -> str:
    payload = path.read_bytes()
    if path.suffix.casefold() in {".csv", ".json", ".md", ".txt"}:
        payload = payload.replace(b"\r\n", b"\n")
    return hashlib.sha256(payload).hexdigest()


def _write_gzip_csv(frame: pd.DataFrame, path: Path) -> None:
    payload = frame.to_csv(index=False, lineterminator="\n").encode("utf-8")
    with path.open("wb") as raw, gzip.GzipFile(
        filename="", mode="wb", fileobj=raw, mtime=0
    ) as compressed:
        compressed.write(payload)


def merge_packages(
    evidence_paths: list[Path],
    audit_paths: list[Path],
    output_dir: Path,
    output_stem: str,
    generated_at_utc: str,
) -> dict[str, object]:
    """Merge packages while re-running the shared trait-specific review gate."""

    if not evidence_paths or len(evidence_paths) != len(audit_paths):
        raise ValueError("supply one audit path for every evidence path")

    evidence = pd.concat(
        [pd.read_csv(path, dtype=str).fillna("") for path in evidence_paths],
        ignore_index=True,
        sort=False,
    ).fillna("")
    audit = pd.concat(
        [pd.read_csv(path, dtype=str).fillna("") for path in audit_paths],
        ignore_index=True,
        sort=False,
    ).fillna("")
    selected, trait_audit, gate = reviewed_source_package_evidence(evidence, audit)
    if not gate["package_gate_passed"]:
        raise ValueError("merged package did not pass the shared source-package gate")

    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_output = output_dir / f"{output_stem}_evidence.csv.gz"
    audit_output = output_dir / f"{output_stem}_audit.csv.gz"
    trait_output = output_dir / f"{output_stem}_trait_audit.csv"
    manifest_output = output_dir / f"{output_stem}_manifest.json"

    _write_gzip_csv(evidence, evidence_output)
    _write_gzip_csv(audit, audit_output)
    trait_audit.to_csv(trait_output, index=False, lineterminator="\n")

    manifest: dict[str, object] = {
        "contract": "merged_reviewed_source_package_v1",
        "generated_at_utc": generated_at_utc,
        "input_packages": [
            {
                "evidence_path": path.as_posix(),
                "evidence_sha256": _sha256(path),
                "audit_path": audit_path.as_posix(),
                "audit_sha256": _sha256(audit_path),
            }
            for path, audit_path in zip(evidence_paths, audit_paths, strict=True)
        ],
        "candidate_evidence_rows": len(evidence),
        "audit_rows": len(audit),
        "selected_evidence_rows": len(selected),
        "selected_species_trait": int(
            selected[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
        ),
        "source_package_gate": gate,
        "output_hashes": {
            evidence_output.name: _sha256(evidence_output),
            audit_output.name: _sha256(audit_output),
            trait_output.name: _sha256(trait_output),
        },
    }
    manifest_output.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path, action="append", required=True)
    parser.add_argument("--audit", type=Path, action="append", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--output-stem", default="merged_reviewed_source_package"
    )
    parser.add_argument(
        "--generated-at-utc",
        default="",
        help="Pinned ISO timestamp; defaults to the current UTC time.",
    )
    args = parser.parse_args()
    timestamp = args.generated_at_utc or datetime.now(UTC).isoformat().replace(
        "+00:00", "Z"
    )
    manifest = merge_packages(
        args.evidence,
        args.audit,
        args.output_dir,
        args.output_stem,
        timestamp,
    )
    print(json.dumps(manifest, ensure_ascii=False))


if __name__ == "__main__":
    main()
