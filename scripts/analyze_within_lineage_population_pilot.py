"""CLI wrapper for the within-lineage population pilot."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from island_v2.within_lineage_population_pilot import write_outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--arabidopsis-csv", type=Path, required=True)
    parser.add_argument("--phormium-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    summary = write_outputs(
        args.arabidopsis_csv,
        args.phormium_csv,
        args.output_dir,
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
