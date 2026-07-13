"""Materialize the preregistered analysis contract as JSON metadata."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml


def build_analysis_contract(contract_yaml: Path, output_json: Path) -> None:
    with contract_yaml.open("r", encoding="utf-8") as handle:
        contract = yaml.safe_load(handle)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract-yaml", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    args = parser.parse_args()
    build_analysis_contract(args.contract_yaml, args.output_json)


if __name__ == "__main__":
    main()
