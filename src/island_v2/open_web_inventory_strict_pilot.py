"""Strict reporting wrapper for discovered-domain inventory backfills."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import typer

from island_v2.open_web_inventory_pilot import run_inventory
from island_v2.open_web_network_pilot import audit_sample

app = typer.Typer(add_completion=False, no_args_is_help=True)


def run_strict(
    *,
    baseline_dir: Path,
    master_csv: Path,
    config_yaml: Path,
    registry_yaml: Path,
    output_dir: Path,
    domain: str,
    max_species_pages: int,
    inventory_pause_seconds: float,
    page_pause_seconds: float,
    completed_species_pages_csv: Path | None = None,
) -> dict[str, object]:
    report = run_inventory(
        baseline_dir=baseline_dir,
        master_csv=master_csv,
        config_yaml=config_yaml,
        registry_yaml=registry_yaml,
        output_dir=output_dir,
        domain=domain,
        max_species_pages=max_species_pages,
        inventory_pause_seconds=inventory_pause_seconds,
        page_pause_seconds=page_pause_seconds,
        completed_species_pages_csv=completed_species_pages_csv,
    )
    all_candidates = pd.read_csv(
        output_dir / "all_inventory_candidates.csv", dtype=str
    ).fillna("")
    if all_candidates.empty:
        contract = all_candidates.copy()
        discovery_only = all_candidates.copy()
    else:
        contract = all_candidates.loc[
            all_candidates["domain_review_status"].eq("approved_registry_trait")
            & all_candidates["already_direct_in_baseline"].astype(str).str.casefold().ne("true")
        ].copy()
        discovery_only = all_candidates.loc[
            ~all_candidates.index.isin(contract.index)
        ].copy()

    contract.to_csv(output_dir / "contract_eligible_novel_candidates.csv", index=False)
    discovery_only.to_csv(output_dir / "discovery_only_or_existing_candidates.csv", index=False)
    contract.to_csv(output_dir / "novel_inventory_candidates.csv", index=False)
    audit_sample(contract, limit=100).to_csv(
        output_dir / "manual_audit_sample.csv", index=False
    )

    report.update(
        {
            "contract_eligible_novel_rows": len(contract),
            "contract_eligible_novel_species": (
                int(contract["accepted_species"].nunique()) if not contract.empty else 0
            ),
            "contract_eligible_novel_species_trait": (
                len(contract.drop_duplicates(["accepted_species", "trait_name"]))
                if not contract.empty
                else 0
            ),
            "contract_eligible_novel_species_axis": (
                len(contract.drop_duplicates(["accepted_species", "axis"]))
                if not contract.empty
                else 0
            ),
            "contract_eligible_by_trait": (
                {
                    str(trait): int(count)
                    for trait, count in contract.drop_duplicates(
                        ["accepted_species", "trait_name"]
                    )["trait_name"].value_counts().items()
                }
                if not contract.empty
                else {}
            ),
            "discovery_only_or_existing_rows": len(discovery_only),
            "promotion_input": "contract_eligible_novel_candidates.csv",
        }
    )
    (output_dir / "inventory_pilot_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(report, ensure_ascii=False, sort_keys=True))
    return report


@app.command()
def main(
    baseline_dir: Path = typer.Option(..., exists=True, file_okay=False),
    master_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    output_dir: Path = typer.Option(..., file_okay=False),
    domain: str = typer.Option(...),
    config_yaml: Path = typer.Option(
        Path("config/open_web_trait_acquisition.yml"), exists=True, dir_okay=False
    ),
    registry_yaml: Path = typer.Option(
        Path("config/open_web_domain_registry.yml"), exists=True, dir_okay=False
    ),
    max_species_pages: int = typer.Option(250, min=1, max=5000),
    inventory_pause_seconds: float = typer.Option(0.5, min=0.1, max=10),
    page_pause_seconds: float = typer.Option(0.35, min=0.1, max=10),
    completed_species_pages_csv: Path | None = typer.Option(
        None,
        exists=True,
        dir_okay=False,
        help="Prior discovered_species_pages.csv; listed URLs are never fetched again.",
    ),
) -> None:
    run_strict(
        baseline_dir=baseline_dir,
        master_csv=master_csv,
        config_yaml=config_yaml,
        registry_yaml=registry_yaml,
        output_dir=output_dir,
        domain=domain,
        max_species_pages=max_species_pages,
        inventory_pause_seconds=inventory_pause_seconds,
        page_pause_seconds=page_pause_seconds,
        completed_species_pages_csv=completed_species_pages_csv,
    )


if __name__ == "__main__":
    app()
