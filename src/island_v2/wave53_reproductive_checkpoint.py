"""Validate the frozen Wave53 support-two reproductive packet."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from island_v2.wave48_reproductive_checkpoint import EXPECTED_SPECIES
from island_v2.wave48_reproductive_checkpoint import validate_packet as validate_reproductive_packet

FORMAL_WAVE33_RUN_ID = 32_932_103_226
BASELINE_WAVE52_RUN_ID = 33_470_056_127
DIRECT_ROWS = 1
EXTERNAL_ROWS = 1
REVIEW_ROWS = 2
IDENTITY_ROWS = 2
REJECTED_ROWS = 2


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
        packet_label="wave53",
        baseline_formal_run_id=BASELINE_WAVE52_RUN_ID,
        expected_direct_rows=DIRECT_ROWS,
        expected_external_rows=EXTERNAL_ROWS,
        expected_review_rows=REVIEW_ROWS,
        expected_identity_rows=IDENTITY_ROWS,
        expected_rejected_rows=REJECTED_ROWS,
        contract="wave53_support_two_reproductive_checkpoint_v1",
        baseline_check_label="immediate_wave52_baseline_pinned",
        allow_completed_direct_axis_enrichment=False,
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
    summary = validate_packet(**vars(args))
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
