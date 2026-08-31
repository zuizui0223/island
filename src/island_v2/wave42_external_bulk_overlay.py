"""Build the lossless Wave42 overlay from external bulk reproductive support."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from island_v2.wave37_europe_pmc_checkpoint import EXPECTED_SPECIES
from island_v2.wave41_external_congener_overlay import build_overlay
from island_v2.wave42_external_bulk_checkpoint import SHARED_LINEAGE

BASELINE_FORMAL_RUN_ID = 33355297355
BASELINE_FORMAL_ARTIFACT = (
    "wave41-external-congener-reproduction-33355297355"
)


def build_wave42_overlay(
    *,
    baseline_csv: Path,
    previous_rule_audit_csv: Path,
    all_evidence_dir: Path,
    external_evidence_csv: Path,
    checkpoint_summary_json: Path,
    output_dir: Path,
    expected_species: int = EXPECTED_SPECIES,
) -> dict[str, Any]:
    summary = build_overlay(
        baseline_csv=baseline_csv,
        previous_rule_audit_csv=previous_rule_audit_csv,
        all_evidence_dir=all_evidence_dir,
        external_evidence_csv=external_evidence_csv,
        checkpoint_summary_json=checkpoint_summary_json,
        output_dir=output_dir,
        expected_species=expected_species,
        wave_label="wave42",
        baseline_formal_run_id=BASELINE_FORMAL_RUN_ID,
        contract="wave42_external_bulk_reproduction_lossless_overlay_v1",
    )
    checkpoint = json.loads(checkpoint_summary_json.read_text(encoding="utf-8"))
    all_evidence = json.loads(
        (all_evidence_dir / "all_evidence_trait_coverage_summary.json").read_text(
            encoding="utf-8"
        )
    )
    combined = all_evidence["source_lineage_audit"]["external_congener_support"]
    delta = summary["delta"]
    delta["new_external_input_species_trait"] = checkpoint["evidence"][
        "species_trait"
    ]
    delta["combined_external_evidence_rows"] = combined["rows"]
    delta["combined_external_resolved_species_trait"] = combined[
        "resolved_species_trait"
    ]
    delta["combined_external_direct_conflicts"] = sum(
        int(value)
        for key, value in combined["cell_resolution_classification_counts"].items()
        if key == "unresolved_direct_conflict"
    )
    delta.pop("external_resolved_species_trait", None)
    delta.pop("external_direct_conflicts", None)
    summary["baseline_formal_artifact"] = BASELINE_FORMAL_ARTIFACT
    summary["external_source_lineage_policy"] = SHARED_LINEAGE
    summary["checks"]["shared_meta_database_lineage_guard"] = (
        checkpoint["evidence"]["shared_redistribution_guard_lineages"] == 1
    )
    output_path = output_dir / "wave42_coverage_summary.json"
    output_path.write_text(
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
    parser.add_argument("--external-evidence-csv", required=True, type=Path)
    parser.add_argument("--checkpoint-summary-json", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-species", type=int, default=EXPECTED_SPECIES)
    args = parser.parse_args()
    summary = build_wave42_overlay(
        baseline_csv=args.baseline_csv,
        previous_rule_audit_csv=args.previous_rule_audit_csv,
        all_evidence_dir=args.all_evidence_dir,
        external_evidence_csv=args.external_evidence_csv,
        checkpoint_summary_json=args.checkpoint_summary_json,
        output_dir=args.output_dir,
        expected_species=args.expected_species,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
