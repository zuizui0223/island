from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RESULTS = ROOT / "results" / "geography_20260924"
FINAL_H1 = ROOT / "results" / "h1_final_traitwise_t_20261004" / "traitwise_results.csv"
DEFAULT_OUTPUT = ROOT / "submission" / "chapter1_current" / "supplement"

SOURCE_FILES = [
    "all/h2_decomposition_models.csv",
    "direct/h2_decomposition_models.csv",
    "h3_original_corrected_comparison.json",
    "h3_corrected_effect_rows.csv.gz",
    "h3_corrected_measurement_cells.csv.gz",
    "h4_exact_corrected.csv",
    "h4_atomic_corrected.csv",
]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _clean_float(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    out = frame.copy()
    for column in columns:
        if column in out.columns:
            out[column] = pd.to_numeric(out[column], errors="coerce")
    return out


def build_s1_summary(results: Path) -> pd.DataFrame:
    repo_root = results.parents[1]
    db = yaml.safe_load(
        (repo_root / "config/chapter1_database_versions/v1.0.0.yml").read_text(
            encoding="utf-8"
        )
    )["database"]
    selector = json.loads(
        (repo_root / "config/chapter1_submission_current.json").read_text(
            encoding="utf-8"
        )
    )
    h3 = json.loads(
        (results / "h3_original_corrected_comparison.json").read_text(encoding="utf-8")
    )["corrected"]["global_gradient"]
    effect_rows = pd.read_csv(results / "h3_corrected_effect_rows.csv.gz")
    measurement_cells = pd.read_csv(results / "h3_corrected_measurement_cells.csv.gz")

    denominator_species = int(db["denominator_species"])
    axes = list(db["axes"])
    rows = [
        ("geography", "analysis_universe", selector["population"]["universe"], "island_units", "config/chapter1_submission_current.json"),
        ("geography", "broad_H1_union", selector["population"]["broad_H1_union"], "island_units", "config/chapter1_submission_current.json"),
        ("traits", "accepted_angiosperm_species", denominator_species, "species", "config/chapter1_database_versions/v1.0.0.yml"),
        ("traits", "possible_species_axis_cells", denominator_species * len(axes), "species_axis_cells", "config/chapter1_database_versions/v1.0.0.yml"),
        ("traits", "resolved_species_axis_cells", db["resolved_cells"], "species_axis_cells", "config/chapter1_database_versions/v1.0.0.yml"),
    ]
    for axis in axes:
        rows.append(
            (
                "traits",
                f"resolved_{axis}_cells",
                db["axis_resolved_cells"][axis],
                "species_axis_cells",
                "config/chapter1_database_versions/v1.0.0.yml",
            )
        )
    rows.extend(
        [
            ("GloPL", "effect_rows", len(effect_rows), "experimental_rows", "results/geography_20260924/h3_corrected_effect_rows.csv.gz"),
            ("GloPL", "measurement_cells", len(measurement_cells), "measurement_cells", "results/geography_20260924/h3_corrected_measurement_cells.csv.gz"),
            ("GloPL", "sites", h3["n_sites"], "sites", "results/geography_20260924/h3_original_corrected_comparison.json"),
            ("GloPL", "publications", h3["n_publications"], "publications", "results/geography_20260924/h3_original_corrected_comparison.json"),
            ("GloPL", "true_continental_zero_sites", selector["population"]["true_continental_site_zeros"], "sites", "config/chapter1_submission_current.json"),
        ]
    )
    return pd.DataFrame(rows, columns=["section", "metric", "value", "unit", "source"])


def build_h2(results: Path) -> pd.DataFrame:
    out = pd.concat(
        [
            pd.read_csv(results / "all/h2_decomposition_models.csv"),
            pd.read_csv(results / "direct/h2_decomposition_models.csv"),
        ],
        ignore_index=True,
    )
    return out[
        [
            "evidence_scope",
            "context",
            "response",
            "analysis_role",
            "model_family",
            "status",
            "n_islands",
            "n_clusters",
            "distance_estimate",
            "distance_se",
            "distance_p",
            "primary_H2b_q",
            "selfing_core_estimate",
            "selfing_core_se",
            "selfing_core_p",
            "kappa",
            "optimizer_success",
        ]
    ]


def build_h3(results: Path) -> pd.DataFrame:
    payload = json.loads(
        (results / "h3_original_corrected_comparison.json").read_text(encoding="utf-8")
    )
    corrected = payload["corrected"]
    rows = []
    for analysis, result in [
        ("primary", corrected["global_gradient"]),
        ("supplemental_only", corrected["sensitivities"]["supplemental_only"]),
        ("no_zero_constant", corrected["sensitivities"]["no_zero_constant"]),
    ]:
        rows.append(
            {
                "analysis": analysis,
                "estimate": result["distance_slope"],
                "se": result["distance_slope_se"],
                "two_sided_p": result["two_sided_p"],
                "one_sided_positive_p": result["one_sided_positive_p"],
                "n_measurement_cells": result["n_cells"],
                "n_publications": result["n_publications"],
                "n_sites": result["n_sites"],
            }
        )
    return pd.DataFrame(rows)


def build_h4_scores(results: Path) -> pd.DataFrame:
    return pd.read_csv(results / "h4_exact_corrected.csv")[
        [
            "family",
            "H2_syndrome",
            "score_name",
            "analysis",
            "inferential_role",
            "evaluable",
            "estimate",
            "se",
            "z",
            "two_sided_p",
            "one_sided_negative_p",
            "n_cells",
            "n_publications",
            "n_sites",
            "n_species",
        ]
    ]


def build_h4_atomic(results: Path) -> pd.DataFrame:
    frame = pd.read_csv(results / "h4_atomic_corrected.csv")
    return frame[
        [
            "family",
            "trait",
            "analysis",
            "inferential_role",
            "evaluable",
            "trait_state_estimate",
            "trait_state_se",
            "trait_state_z",
            "trait_state_two_sided_p",
            "trait_state_one_sided_negative_p",
            "n_cells",
            "n_publications",
            "n_sites",
            "n_species",
            "reason",
        ]
    ]


def write_tables(results: Path, output: Path) -> dict[str, object]:
    output.mkdir(parents=True, exist_ok=True)
    tables = {
        "Table_S1_data_summary.csv": build_s1_summary(results),
        "Table_S3_H2_decomposition.csv": build_h2(results),
        "Table_S5_H3_pollen_limitation.csv": build_h3(results),
        "Table_S6a_H4_scores.csv": build_h4_scores(results),
        "Table_S6b_H4_atomic.csv": build_h4_atomic(results),
    }
    output_meta = {}
    for name, frame in tables.items():
        path = output / name
        frame.to_csv(path, index=False, float_format="%.12g")
        output_meta[name] = {
            "rows": len(frame),
            "columns": list(frame.columns),
            "sha256": _sha256(path),
        }

    h1_name = "Table_S2_H1_traitwise.csv"
    h1_target = output / h1_name
    h1_target.write_bytes(FINAL_H1.read_bytes())
    h1_frame = pd.read_csv(FINAL_H1)
    output_meta[h1_name] = {
        "rows": len(h1_frame),
        "columns": list(h1_frame.columns),
        "sha256": _sha256(h1_target),
    }

    source_meta = {rel: _sha256(results / rel) for rel in SOURCE_FILES}
    source_meta["results/h1_final_traitwise_t_20261004/traitwise_results.csv"] = _sha256(FINAL_H1)
    repo_root = results.parents[1]
    rel = "config/chapter1_database_versions/v1.0.0.yml"
    source_meta[rel] = _sha256(repo_root / rel)
    manifest = {
        "contract": "chapter1_submission_supplement_tables_v1",
        "scientific_surface": "corrected_geography_20260924_plus_final_traitwise_H1_20261004",
        "sources": source_meta,
        "outputs": output_meta,
        "notes": {
            "Table_S2": (
                "Final active H1: seven separate traits in four regions, with broad All "
                "primary plus WCVP and Direct-only sensitivities. No pooled/domain score."
            ),
            "Table_S4": (
                "Full raw colour and colour-by-architecture results remain in "
                "results/geography_20260924/{all,direct}/raw_patterns/ rather than "
                "being duplicated into one very large submission CSV."
            ),
        },

    }
    manifest_path = output / "SUPPLEMENT_TABLES_MANIFEST.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    manifest = write_tables(args.results, args.output)
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
