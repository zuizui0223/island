"""Rebuild only Wave53-touched trait-specific reproductive rules."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from island_v2.wave37_europe_pmc_checkpoint import EXPECTED_SPECIES
from island_v2.wave48_incremental_all_evidence import build_incremental_audit

EXPECTED_NEW_RULES = frozenset(
    {
        ("Gomphrena", "reproductive_assurance", "self_incompatibility"),
        ("Marcgravia", "reproductive_assurance", "self_incompatibility"),
    }
)


def build_wave53_incremental_audit(
    *,
    master_csv: Path,
    ontology_yaml: Path,
    baseline_coverage_csv: Path,
    previous_rule_audit_csv: Path,
    previous_resolved_direct_csv: Path,
    previous_external_resolved_csv: Path,
    previous_external_conflicts_csv: Path,
    previous_rebuilt_low_csv: Path,
    new_direct_evidence_csv: Path,
    new_external_evidence_csv: Path,
    output_dir: Path,
    expected_species: int = EXPECTED_SPECIES,
) -> dict[str, Any]:
    return build_incremental_audit(
        master_csv=master_csv,
        ontology_yaml=ontology_yaml,
        baseline_coverage_csv=baseline_coverage_csv,
        previous_rule_audit_csv=previous_rule_audit_csv,
        previous_resolved_direct_csv=previous_resolved_direct_csv,
        previous_external_resolved_csv=previous_external_resolved_csv,
        previous_external_conflicts_csv=previous_external_conflicts_csv,
        previous_rebuilt_low_csv=previous_rebuilt_low_csv,
        new_direct_evidence_csv=new_direct_evidence_csv,
        new_external_evidence_csv=new_external_evidence_csv,
        output_dir=output_dir,
        expected_species=expected_species,
        expected_direct_rows=1,
        expected_external_rows=1,
        expected_new_rules=EXPECTED_NEW_RULES,
        expected_blocked_rules=frozenset(),
        expected_counterexample_rules=frozenset(),
        output_label="wave53",
        contract="wave53_incremental_all_evidence_touched_rule_rebuild_v1",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--master-csv", required=True, type=Path)
    parser.add_argument("--ontology-yaml", required=True, type=Path)
    parser.add_argument("--baseline-coverage-csv", required=True, type=Path)
    parser.add_argument("--previous-rule-audit-csv", required=True, type=Path)
    parser.add_argument("--previous-resolved-direct-csv", required=True, type=Path)
    parser.add_argument("--previous-external-resolved-csv", required=True, type=Path)
    parser.add_argument("--previous-external-conflicts-csv", required=True, type=Path)
    parser.add_argument("--previous-rebuilt-low-csv", required=True, type=Path)
    parser.add_argument("--new-direct-evidence-csv", required=True, type=Path)
    parser.add_argument("--new-external-evidence-csv", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-species", type=int, default=EXPECTED_SPECIES)
    args = parser.parse_args()
    summary = build_wave53_incremental_audit(**vars(args))
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
