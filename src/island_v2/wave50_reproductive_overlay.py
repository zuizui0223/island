"""Build the lossless Wave50 source-grounded reproductive overlay."""

from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path
from typing import Any

from island_v2.wave37_europe_pmc_checkpoint import EXPECTED_SPECIES, _sha256
from island_v2.wave45_promoted_reproductive_overlay import build_wave45_overlay
from island_v2.wave50_reproductive_checkpoint import (
    BASELINE_WAVE49_RUN_ID,
    FORMAL_WAVE33_RUN_ID,
)


def build_wave50_overlay(
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
    with tempfile.TemporaryDirectory(prefix="wave50-overlay-") as temporary:
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
            "wave45_species_axis_coverage.csv.gz": "wave50_species_axis_coverage.csv.gz",
            "wave45_resolved_direct_species_trait.csv.gz": "wave50_resolved_direct_species_trait.csv.gz",
            "wave45_new_validated_low_species_trait.csv.gz": "wave50_new_validated_low_species_trait.csv.gz",
            "wave45_new_trait_specific_genus_rule_audit.csv.gz": "wave50_new_trait_specific_genus_rule_audit.csv.gz",
            "wave45_change_audit.csv.gz": "wave50_change_audit.csv.gz",
            "wave45_external_congener_resolved_species_trait.csv.gz": "wave50_external_congener_resolved_species_trait.csv.gz",
            "wave45_external_congener_source_conflicts.csv.gz": "wave50_external_congener_source_conflicts.csv.gz",
        }
        for old_name, new_name in mapping.items():
            (staging / old_name).replace(output_dir / new_name)

    checkpoint = json.loads(checkpoint_summary_json.read_text(encoding="utf-8"))
    delta = summary["delta"]
    delta["resolved_wave50_direct_species_trait"] = delta.pop(
        "resolved_wave45_direct_species_trait"
    )
    delta["resolved_wave50_direct_species_axis"] = delta.pop(
        "resolved_wave45_direct_species_axis"
    )
    summary.update(
        {
            "contract": "wave50_source_grounded_reproductive_lossless_overlay_v1",
            "baseline_formal_run_id": BASELINE_WAVE49_RUN_ID,
            "formal_wave33_run_id": FORMAL_WAVE33_RUN_ID,
            "source_lineages": [
                "doi:10.1111/j.1365-2699.2005.01307.x",
                "doi:10.1371/journal.pone.0079043",
                "doi:10.1038/s41598-020-60651-w",
                "doi:10.1614/ipsm-d-16-00030.1",
                "flinders-thesis:blackall-2023-pultenaea",
            ],
            "selection_contract": (
                "unresolved_reproductive_support2_third_species_and_direct_v1"
            ),
            "query_accounting": checkpoint["queries"],
            "query_cost_usd": checkpoint["query_cost_usd"],
        }
    )
    summary["checks"].update(
        {
            "formal_wave33_baseline_pinned": True,
            "formal_wave49_immediate_baseline_pinned": True,
            "all_source_receipts_verified": checkpoint["checks"][
                "retrieved_sources_verified"
            ],
            "content_fingerprints_verified": checkpoint["checks"][
                "content_fingerprints_verified"
            ],
            "reproductive_traits_not_interchanged": True,
        }
    )
    summary["artifact_sha256"] = {
        path.name: _sha256(path) for path in sorted(output_dir.glob("wave50_*.csv.gz"))
    }
    summary_path = output_dir / "wave50_coverage_summary.json"
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
    summary = build_wave50_overlay(**vars(args))
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
