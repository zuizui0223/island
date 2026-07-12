"""Family-based angiosperm eligibility gate, applied before floral trait fills.

The GBIF acquisition universe is Tracheophyta, but non-seed vascular plants and
gymnosperms have no flowers or comparable pollen-vector states. Assigning them
global floral priors would be nonsense, so this gate classifies every master
species by family into an angiosperm-eligible flowering-plant universe versus an
out-of-scope remainder (gymnosperm, non-seed vascular, fossil non-angiosperm, or
family_unresolved when the family is missing). It is deliberately upstream of the
trait cascade: out-of-scope species are never filled, only marked.

This is a coarse, source-blind family gate for coverage bookkeeping, not the
reviewed per-species taxon-scope adjudication in ``taxon_scope``. It never
promotes a species into a confirmatory denominator; it only decides who is
eligible to receive a floral fill.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)


@app.callback()
def main() -> None:
    """Classify master species into the angiosperm-eligible flowering-plant universe."""


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def load_config(path: Path) -> dict[str, Any]:
    """Load and minimally validate the versioned angiosperm-scope configuration."""
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    required = {"family_column", "eligible_group", "unresolved_group", "non_angiosperm_families"}
    if not isinstance(config, dict) or not required.issubset(config):
        raise typer.BadParameter(
            "angiosperm-scope config must contain family_column, eligible_group, "
            "unresolved_group, and non_angiosperm_families"
        )
    return config


def build_family_group_map(config: dict[str, Any]) -> dict[str, str]:
    """Flatten the per-group non-angiosperm family lists into family -> group."""
    mapping: dict[str, str] = {}
    for group, families in config["non_angiosperm_families"].items():
        for family in families or []:
            mapping[_text(family)] = str(group)
    return mapping


def classify_scope(master: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    """Return per-species taxonomic_group and angiosperm eligibility.

    A species is angiosperm-eligible when its family is present and not on any
    non-angiosperm list. A missing family is unresolved and out of scope.
    """
    family_column = config["family_column"]
    if "accepted_species" not in master.columns or family_column not in master.columns:
        raise typer.BadParameter(f"master must contain accepted_species and {family_column!r}")

    group_map = build_family_group_map(config)
    eligible_group = str(config["eligible_group"])
    unresolved_group = str(config["unresolved_group"])

    families = master[family_column].map(_text)

    def classify(family: str) -> str:
        if not family:
            return unresolved_group
        return group_map.get(family, eligible_group)

    groups = families.map(classify)
    out = pd.DataFrame(
        {
            "accepted_species": master["accepted_species"].map(_text),
            "family": families,
            "taxonomic_group": groups,
            "angiosperm_analysis_eligible": groups.eq(eligible_group),
        }
    )
    return out.loc[out["accepted_species"].ne("")].drop_duplicates("accepted_species").reset_index(drop=True)


def scope_summary(scope: pd.DataFrame) -> dict[str, Any]:
    """Summarise the eligible flowering-plant universe versus out-of-scope groups."""
    n_total = int(len(scope))
    n_eligible = int(scope["angiosperm_analysis_eligible"].sum())
    by_group = {str(k): int(v) for k, v in scope["taxonomic_group"].value_counts().items()}
    return {
        "version": "1.0",
        "n_master_species": n_total,
        "n_angiosperm_eligible": n_eligible,
        "angiosperm_eligible_rate": n_eligible / n_total if n_total else 0.0,
        "n_out_of_scope": n_total - n_eligible,
        "by_taxonomic_group": by_group,
        "interpretation": (
            "Family-based angiosperm eligibility for floral/pollination denominators. "
            "Out-of-scope species are never assigned floral priors; they are marked, "
            "not filled. Coarse family gate, not reviewed per-species taxon scope."
        ),
    }


@app.command("run")
def run(
    master_taxa_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/angiosperm_scope.yml"), exists=True),
) -> None:
    """Write per-species angiosperm scope classification and a summary."""
    config = load_config(config_path)
    master = pd.read_csv(master_taxa_csv, dtype=str).fillna("")
    scope = classify_scope(master, config)
    summary = scope_summary(scope)

    output_dir.mkdir(parents=True, exist_ok=True)
    scope.to_csv(output_dir / "angiosperm_scope_by_species.csv", index=False)
    (output_dir / "angiosperm_scope_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(
        f"Classified {summary['n_master_species']} species: "
        f"{summary['n_angiosperm_eligible']} angiosperm-eligible "
        f"({summary['angiosperm_eligible_rate']:.1%}), "
        f"{summary['n_out_of_scope']} out-of-scope {summary['by_taxonomic_group']}."
    )


if __name__ == "__main__":
    app()
