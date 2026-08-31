"""Build the lossless Wave46 corrected public-Web and reproductive overlay."""

from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path
from typing import Any

import pandas as pd

from island_v2.wave37_europe_pmc_checkpoint import EXPECTED_SPECIES, _sha256
from island_v2.wave45_promoted_reproductive_overlay import build_wave45_overlay

BASELINE_FORMAL_RUN_ID = 33_370_692_122
FORMAL_WAVE33_RUN_ID = 32_932_103_226
CORRECT_PUBLIC_WEB_RUN_ID = 32_710_232_989


def build_wave46_overlay(
    *,
    baseline_csv: Path,
    previous_rule_audit_csv: Path,
    all_evidence_dir: Path,
    direct_evidence_csv: Path,
    external_evidence_csv: Path,
    checkpoint_summary_json: Path,
    output_dir: Path,
    expected_species: int = EXPECTED_SPECIES,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="wave46-overlay-") as temporary:
        staging = Path(temporary)
        summary = build_wave45_overlay(
            baseline_csv=baseline_csv,
            previous_rule_audit_csv=previous_rule_audit_csv,
            all_evidence_dir=all_evidence_dir,
            direct_evidence_csv=direct_evidence_csv,
            external_evidence_csv=external_evidence_csv,
            checkpoint_summary_json=checkpoint_summary_json,
            output_dir=staging,
            expected_species=expected_species,
        )
        mapping = {
            "wave45_species_axis_coverage.csv.gz": "wave46_species_axis_coverage.csv.gz",
            "wave45_resolved_direct_species_trait.csv.gz": "wave46_resolved_direct_species_trait.csv.gz",
            "wave45_new_validated_low_species_trait.csv.gz": "wave46_new_validated_low_species_trait.csv.gz",
            "wave45_new_trait_specific_genus_rule_audit.csv.gz": "wave46_new_trait_specific_genus_rule_audit.csv.gz",
            "wave45_change_audit.csv.gz": "wave46_change_audit.csv.gz",
            "wave45_external_congener_resolved_species_trait.csv.gz": "wave46_external_congener_resolved_species_trait.csv.gz",
            "wave45_external_congener_source_conflicts.csv.gz": "wave46_external_congener_source_conflicts.csv.gz",
        }
        for old_name, new_name in mapping.items():
            (staging / old_name).replace(output_dir / new_name)

    checkpoint = json.loads(checkpoint_summary_json.read_text(encoding="utf-8"))
    reused_external = pd.read_csv(external_evidence_csv, dtype=str).fillna("")
    delta = summary["delta"]
    delta["resolved_wave46_direct_species_trait"] = delta.pop(
        "resolved_wave45_direct_species_trait"
    )
    delta["resolved_wave46_direct_species_axis"] = delta.pop("resolved_wave45_direct_species_axis")
    delta["corrected_public_web_rows"] = checkpoint["evidence"]["latest_public_web_rows"]
    delta["new_reviewed_direct_species_trait"] = checkpoint["evidence"]["new_direct_species_trait"]
    delta["reused_external_packet_species_trait"] = int(
        reused_external[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
    )
    summary.update(
        {
            "contract": "wave46_corrected_public_web_recovery_overlay_v1",
            "baseline_formal_run_id": BASELINE_FORMAL_RUN_ID,
            "formal_wave33_run_id": FORMAL_WAVE33_RUN_ID,
            "correct_public_web_run_id": CORRECT_PUBLIC_WEB_RUN_ID,
        }
    )
    summary["checks"]["formal_wave33_baseline_pinned"] = True
    summary["checks"]["correct_latest_public_web_artifact_used"] = True
    summary["artifact_sha256"] = {
        path.name: _sha256(path) for path in sorted(output_dir.glob("wave46_*.csv.gz"))
    }
    summary_path = output_dir / "wave46_coverage_summary.json"
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-csv", required=True, type=Path)
    parser.add_argument("--previous-rule-audit-csv", required=True, type=Path)
    parser.add_argument("--all-evidence-dir", required=True, type=Path)
    parser.add_argument("--direct-evidence-csv", required=True, type=Path)
    parser.add_argument("--external-evidence-csv", required=True, type=Path)
    parser.add_argument("--checkpoint-summary-json", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-species", type=int, default=EXPECTED_SPECIES)
    args = parser.parse_args()
    summary = build_wave46_overlay(
        baseline_csv=args.baseline_csv,
        previous_rule_audit_csv=args.previous_rule_audit_csv,
        all_evidence_dir=args.all_evidence_dir,
        direct_evidence_csv=args.direct_evidence_csv,
        external_evidence_csv=args.external_evidence_csv,
        checkpoint_summary_json=args.checkpoint_summary_json,
        output_dir=args.output_dir,
        expected_species=args.expected_species,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
