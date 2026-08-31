"""Validate the frozen Wave50 reproductive evidence and source receipts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from island_v2.wave48_reproductive_checkpoint import (
    EXPECTED_SPECIES,
    validate_packet as validate_reproductive_packet,
)

FORMAL_WAVE33_RUN_ID = 32_932_103_226
BASELINE_WAVE49_RUN_ID = 33_397_671_718
DIRECT_ROWS = 5
EXTERNAL_ROWS = 0
REVIEW_ROWS = 5
IDENTITY_ROWS = 5
REJECTED_ROWS = 9


def validate_packet(
    *,
    packet_dir: Path,
    target_coverage_csv: Path,
    retrieved_source_dir: Path,
    output_dir: Path,
    output_json: Path,
    expected_species: int = EXPECTED_SPECIES,
) -> dict[str, Any]:
    return validate_reproductive_packet(
        packet_dir=packet_dir,
        target_coverage_csv=target_coverage_csv,
        retrieved_source_dir=retrieved_source_dir,
        output_dir=output_dir,
        output_json=output_json,
        expected_species=expected_species,
        packet_label="wave50",
        baseline_formal_run_id=BASELINE_WAVE49_RUN_ID,
        expected_direct_rows=DIRECT_ROWS,
        expected_external_rows=EXTERNAL_ROWS,
        expected_review_rows=REVIEW_ROWS,
        expected_identity_rows=IDENTITY_ROWS,
        expected_rejected_rows=REJECTED_ROWS,
        contract="wave50_source_grounded_reproductive_checkpoint_v1",
        baseline_check_label="immediate_wave49_baseline_pinned",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--packet-dir", required=True, type=Path)
    parser.add_argument("--target-coverage-csv", required=True, type=Path)
    parser.add_argument("--retrieved-source-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    parser.add_argument("--expected-species", type=int, default=EXPECTED_SPECIES)
    args = parser.parse_args()
    summary = validate_packet(
        packet_dir=args.packet_dir,
        target_coverage_csv=args.target_coverage_csv,
        retrieved_source_dir=args.retrieved_source_dir,
        output_dir=args.output_dir,
        output_json=args.output_json,
        expected_species=args.expected_species,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
