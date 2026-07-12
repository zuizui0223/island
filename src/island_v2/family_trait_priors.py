"""Generate low-confidence likely floral-trait candidates from declared family priors.

These are deliberately separated from reported evidence. The mappings below are
broad architectural regularities useful for sensitivity analyses and acquisition
triage, not species-level observations. A prior is emitted only when the species
still lacks a reported candidate for that trait.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import typer

from island_v2.search_result_trait_candidates import OUTPUT_COLUMNS

app = typer.Typer(add_completion=False, no_args_is_help=True)

# Only relatively diagnostic family-level floral architectures are included.
# Ambiguous families are intentionally omitted rather than forced into a state.
FAMILY_TRAIT_PRIORS: dict[str, tuple[tuple[str, str, str], ...]] = {
    "Asteraceae": (("floral_form", "composite_head", "family_typical_composite_head"),),
    "Poaceae": (("floral_form", "reduced_wind", "family_typical_reduced_grass_flower"),),
    "Cyperaceae": (("floral_form", "reduced_wind", "family_typical_reduced_sedge_flower"),),
    "Orchidaceae": (("floral_symmetry", "zygomorphic", "family_typical_orchid_zygomorphy"),),
    "Lamiaceae": (
        ("floral_symmetry", "zygomorphic", "family_typical_lamiaceae_zygomorphy"),
        ("floral_form", "tubular", "family_typical_fused_tubular_corolla"),
    ),
    "Acanthaceae": (
        ("floral_symmetry", "zygomorphic", "family_typical_acanthaceae_zygomorphy"),
        ("floral_form", "tubular", "family_typical_fused_tubular_corolla"),
    ),
    "Gesneriaceae": (
        ("floral_symmetry", "zygomorphic", "family_typical_gesneriaceae_zygomorphy"),
        ("floral_form", "tubular", "family_typical_fused_tubular_corolla"),
    ),
    "Campanulaceae": (("floral_form", "bell_campanulate", "family_typical_campanulate_corolla"),),
    "Brassicaceae": (
        ("floral_symmetry", "actinomorphic", "family_typical_brassicaceae_symmetry"),
        ("floral_form", "open_radial", "family_typical_open_cruciform_flower"),
    ),
    "Rosaceae": (("floral_symmetry", "actinomorphic", "family_typical_rosaceae_symmetry"),),
    "Apocynaceae": (
        ("floral_symmetry", "actinomorphic", "family_typical_apocynaceae_symmetry"),
        ("floral_form", "tubular", "family_typical_fused_corolla_tube"),
    ),
    "Rubiaceae": (
        ("floral_symmetry", "actinomorphic", "family_typical_rubiaceae_symmetry"),
        ("floral_form", "tubular", "family_typical_fused_corolla_tube"),
    ),
    "Melastomataceae": (("floral_symmetry", "actinomorphic", "family_typical_melastomataceae_symmetry"),),
    "Lauraceae": (("floral_symmetry", "actinomorphic", "family_typical_lauraceae_symmetry"),),
    "Primulaceae": (("floral_symmetry", "actinomorphic", "family_typical_primulaceae_symmetry"),),
    "Annonaceae": (("floral_symmetry", "actinomorphic", "family_typical_annonaceae_symmetry"),),
    "Araceae": (("floral_form", "other_described", "family_typical_spadix_spathe_architecture"),),
}


def _text(value: object) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return " ".join(str(value).strip().split())


def generate_family_trait_priors(
    species_table: pd.DataFrame,
    existing_candidates: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, dict[str, object]]:
    if "accepted_species" not in species_table.columns:
        raise typer.BadParameter("species table missing accepted_species")
    if "family" not in species_table.columns:
        raise typer.BadParameter("species table missing family")

    reported_pairs: set[tuple[str, str]] = set()
    if existing_candidates is not None and not existing_candidates.empty:
        table = existing_candidates.fillna("")
        if {"accepted_species", "trait_name"}.issubset(table.columns):
            if "candidate_class" in table.columns:
                table = table.loc[table["candidate_class"].eq("reported")]
            reported_pairs = {
                (_text(row.get("accepted_species")), _text(row.get("trait_name")))
                for row in table.to_dict("records")
                if _text(row.get("accepted_species")) and _text(row.get("trait_name"))
            }

    rows: list[dict[str, str]] = []
    for record in species_table.fillna("").to_dict("records"):
        species = _text(record.get("accepted_species"))
        family = _text(record.get("family"))
        if not species:
            continue
        for trait_name, value, rule in FAMILY_TRAIT_PRIORS.get(family, ()):
            if (species, trait_name) in reported_pairs:
                continue
            rows.append(
                {
                    "accepted_species": species,
                    "trait_name": trait_name,
                    "provisional_candidate_value": value,
                    "candidate_class": "proxy",
                    "inference_status": "likely",
                    "confidence": "low",
                    "evidence_quote": f"family={family}",
                    "inference_rule": f"likely_{rule}",
                    "searched_name": species,
                    "search_query": "",
                    "search_engine": "",
                    "search_rank": "",
                    "source_url": "",
                    "source_title": "",
                    "source_type": "declared_family_trait_prior",
                    "evidence_scope": "proxy_not_reported",
                    "review_status": "unreviewed_proxy",
                }
            )

    frame = pd.DataFrame(rows, columns=OUTPUT_COLUMNS).drop_duplicates()
    report: dict[str, object] = {
        "n_species_input": int(species_table["accepted_species"].astype(str).str.strip().ne("").sum()),
        "n_prior_rows": int(len(frame)),
        "n_species_with_prior": int(frame["accepted_species"].nunique()) if len(frame) else 0,
        "coverage_by_trait": {
            str(k): int(v)
            for k, v in frame.groupby("trait_name")["accepted_species"].nunique().sort_index().items()
        }
        if len(frame)
        else {},
        "policy": (
            "All rows are low-confidence likely proxies from declared family-level floral architecture. "
            "They are never reported species-level evidence and are emitted only for traits lacking a "
            "reported candidate in the supplied candidate table."
        ),
    }
    return frame, report


@app.command()
def generate(
    species_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    existing_candidates_csv: Path | None = typer.Option(None, exists=True),
) -> None:
    species = pd.read_csv(species_csv, dtype=str).fillna("")
    existing = (
        pd.read_csv(existing_candidates_csv, dtype=str).fillna("")
        if existing_candidates_csv is not None
        else None
    )
    frame, report = generate_family_trait_priors(species, existing)
    output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_dir / "family_trait_likely_candidates.csv", index=False)
    (output_dir / "family_trait_likely_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
