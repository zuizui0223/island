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
FINAL_AUDIT_DIR = ROOT / "results" / "h1_final_directional_20261003"
FINAL_H2_ALL = FINAL_AUDIT_DIR / "h2_all_finite_cluster_audit.csv"
FINAL_H2_DIRECT = FINAL_AUDIT_DIR / "h2_direct_finite_cluster_audit.csv"
FINAL_H3 = FINAL_AUDIT_DIR / "h3_finite_publication_audit.csv"
FINAL_H4 = FINAL_AUDIT_DIR / "h4_finite_publication_audit.csv"
DEFAULT_OUTPUT = ROOT / "submission" / "chapter1_current" / "supplement"

SOURCE_FILES = [
    "h3_original_corrected_comparison.json",
    "h3_corrected_effect_rows.csv.gz",
    "h3_corrected_measurement_cells.csv.gz",
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


def write_tables(results: Path, output: Path) -> dict[str, object]:
    output.mkdir(parents=True, exist_ok=True)
    tables = {
        "Table_S1_data_summary.csv": build_s1_summary(results),
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

    copied = {
        "Table_S2_H1_traitwise.csv": FINAL_H1,
        "Table_S3a_H2_all_finite_cluster.csv": FINAL_H2_ALL,
        "Table_S3b_H2_direct_finite_cluster.csv": FINAL_H2_DIRECT,
        "Table_S5_H3_pollen_limitation.csv": FINAL_H3,
        "Table_S6a_H4_scores.csv": FINAL_H4,
    }
    for name, source in copied.items():
        target = output / name
        target.write_bytes(source.read_bytes())
        frame = pd.read_csv(source)
        output_meta[name] = {
            "rows": len(frame),
            "columns": list(frame.columns),
            "sha256": _sha256(target),
        }

    source_meta = {rel: _sha256(results / rel) for rel in SOURCE_FILES}
    source_meta["results/h1_final_traitwise_t_20261004/traitwise_results.csv"] = _sha256(FINAL_H1)
    source_meta["results/h1_final_directional_20261003/h2_all_finite_cluster_audit.csv"] = _sha256(FINAL_H2_ALL)
    source_meta["results/h1_final_directional_20261003/h2_direct_finite_cluster_audit.csv"] = _sha256(FINAL_H2_DIRECT)
    source_meta["results/h1_final_directional_20261003/h3_finite_publication_audit.csv"] = _sha256(FINAL_H3)
    source_meta["results/h1_final_directional_20261003/h4_finite_publication_audit.csv"] = _sha256(FINAL_H4)
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
            "Table_S3": (
                "Active H2 tables are copied from the final finite-cluster t audits, "
                "separately for All and Direct-only evidence."
            ),
            "Table_S5_S6": (
                "Active H3 and exact-score H4 tables are copied from the final "
                "finite-publication t audits."
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
