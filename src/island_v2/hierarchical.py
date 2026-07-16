"""Build pending-review genus/family trait distributions from accepted evidence."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)


def load_csv(path: Path, required: set[str]) -> pd.DataFrame:
    data = pd.read_csv(path, dtype=str).fillna("")
    absent = required - set(data.columns)
    if absent:
        raise typer.BadParameter(f"Missing columns: {', '.join(sorted(absent))}")
    return data


def get_rule(policy: dict[str, Any], trait: str, level: str) -> dict[str, Any]:
    return policy.get("trait_rules", {}).get(trait, {}).get(level, policy["default"][level])


def make_candidates(
    taxa: pd.DataFrame,
    evidence: pd.DataFrame,
    ontology: dict[str, Any],
    policy: dict[str, Any],
) -> pd.DataFrame:
    direct_scopes = set(policy["accepted_direct_scopes"])
    source_grades = set(policy["accepted_source_reliability"])
    records = evidence.loc[
        evidence["review_status"].eq("accepted")
        & evidence["evidence_scope"].isin(direct_scopes)
        & evidence["source_reliability"].isin(source_grades)
        & evidence["standardized_value"].ne("unresolved")
    ].copy()
    records = records.merge(taxa[["accepted_species", "genus", "family"]], on="accepted_species", how="left")
    already_direct = set(zip(records["accepted_species"], records["trait_name"], strict=False))
    rows: list[dict[str, Any]] = []

    for focal in taxa.drop_duplicates("accepted_species").to_dict("records"):
        for trait_name in ontology["traits"]:
            if (focal["accepted_species"], trait_name) in already_direct:
                continue
            for level, scope in (("genus", "genus_inference"), ("family", "family_inference")):
                group = records.loc[
                    (records["trait_name"].eq(trait_name))
                    & (records[level].eq(focal[level]))
                    & (records["accepted_species"].ne(focal["accepted_species"]))
                ].drop_duplicates(subset=["accepted_species", "trait_name"])
                if group.empty:
                    continue
                n_species = group["accepted_species"].nunique()
                probs = group["standardized_value"].value_counts(normalize=True).sort_index().to_dict()
                modal_value = max(probs, key=probs.get)
                modal_probability = float(probs[modal_value])
                rule = get_rule(policy, trait_name, level)
                if not bool(rule.get("enabled", True)):
                    continue
                if n_species < int(rule["minimum_supporting_species"]):
                    continue
                if modal_probability < float(rule["minimum_modal_probability"]):
                    continue
                rows.append(
                    {
                        "accepted_species": focal["accepted_species"],
                        "genus": focal["genus"],
                        "family": focal["family"],
                        "trait_name": trait_name,
                        "candidate_kind": "hierarchical_inference",
                        "evidence_scope": scope,
                        "standardized_value": modal_value,
                        "category_probabilities_json": json.dumps(probs, ensure_ascii=False),
                        "modal_probability": modal_probability,
                        "supporting_taxa_json": json.dumps(sorted(group["accepted_species"].unique()), ensure_ascii=False),
                        "supporting_evidence_ids_json": json.dumps(sorted(group["evidence_id"].unique()), ensure_ascii=False),
                        "inference_rule": f"{level} distribution from {n_species} accepted supporting species",
                        "confidence": "medium" if level == "genus" else "low",
                        "review_status": "pending",
                        "analysis_track": f"expanded_{level}_distribution",
                    }
                )
                break
    return pd.DataFrame(rows)


@app.command()
def run(
    taxa_csv: Path = typer.Option(..., exists=True),
    accepted_evidence_csv: Path = typer.Option(..., exists=True),
    output_csv: Path = typer.Option(...),
    ontology_path: Path = typer.Option(Path("config/trait_ontology.yml")),
    policy_path: Path = typer.Option(Path("config/hierarchical_inference.yml")),
) -> None:
    """Write probability-preserving inferred candidates; never overwrite evidence."""
    taxa = load_csv(taxa_csv, {"accepted_species", "genus", "family"})
    evidence = load_csv(
        accepted_evidence_csv,
        {"evidence_id", "accepted_species", "trait_name", "standardized_value", "evidence_scope", "source_reliability", "review_status"},
    )
    ontology = yaml.safe_load(ontology_path.read_text(encoding="utf-8"))
    policy = yaml.safe_load(policy_path.read_text(encoding="utf-8"))
    output = make_candidates(taxa, evidence, ontology, policy)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(output_csv, index=False)
    typer.echo(f"Wrote {len(output)} inference candidates to {output_csv}")


if __name__ == "__main__":
    app()
