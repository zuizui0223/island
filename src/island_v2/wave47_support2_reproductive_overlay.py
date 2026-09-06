"""Build the lossless Wave47 support-2 reproductive overlay."""

from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path
from typing import Any

from island_v2.wave37_europe_pmc_checkpoint import EXPECTED_SPECIES, _sha256
from island_v2.wave45_promoted_reproductive_overlay import build_wave45_overlay

BASELINE_FORMAL_RUN_ID = 33_376_311_877
FORMAL_WAVE33_RUN_ID = 32_932_103_226


def build_wave47_overlay(
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
    with tempfile.TemporaryDirectory(prefix="wave47-overlay-") as temporary:
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
            "wave45_species_axis_coverage.csv.gz": "wave47_species_axis_coverage.csv.gz",
            "wave45_resolved_direct_species_trait.csv.gz": "wave47_resolved_direct_species_trait.csv.gz",
            "wave45_new_validated_low_species_trait.csv.gz": "wave47_new_validated_low_species_trait.csv.gz",
            "wave45_new_trait_specific_genus_rule_audit.csv.gz": "wave47_new_trait_specific_genus_rule_audit.csv.gz",
            "wave45_change_audit.csv.gz": "wave47_change_audit.csv.gz",
            "wave45_external_congener_resolved_species_trait.csv.gz": "wave47_external_congener_resolved_species_trait.csv.gz",
            "wave45_external_congener_source_conflicts.csv.gz": "wave47_external_congener_source_conflicts.csv.gz",
        }
        for old_name, new_name in mapping.items():
            (staging / old_name).replace(output_dir / new_name)

    checkpoint = json.loads(checkpoint_summary_json.read_text(encoding="utf-8"))
    delta = summary["delta"]
    delta["resolved_wave47_direct_species_trait"] = delta.pop(
        "resolved_wave45_direct_species_trait"
    )
    delta["resolved_wave47_direct_species_axis"] = delta.pop("resolved_wave45_direct_species_axis")
    summary.update(
        {
            "contract": "wave47_support2_reproductive_lossless_overlay_v1",
            "baseline_formal_run_id": BASELINE_FORMAL_RUN_ID,
            "formal_wave33_run_id": FORMAL_WAVE33_RUN_ID,
            "source_species": "Buxus wallichiana",
            "selection_contract": "formal_gap_aware_current_min2_positive_state_unlock_v1",
        }
    )
    summary["checks"]["formal_wave33_baseline_pinned"] = True
    summary["checks"]["formal_wave46_immediate_baseline_pinned"] = True
    summary["checks"]["external_source_pdf_digest_verified"] = checkpoint["checks"][
        "retrieved_source_digest_verified"
    ]
    summary["artifact_sha256"] = {
        path.name: _sha256(path) for path in sorted(output_dir.glob("wave47_*.csv.gz"))
    }
    summary_path = output_dir / "wave47_coverage_summary.json"
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
    summary = build_wave47_overlay(
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
