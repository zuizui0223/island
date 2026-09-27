from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RESULTS = ROOT / "results" / "geography_20260924"
DEFAULT_OUTPUT = ROOT / "submission" / "chapter1_current" / "supplement"

SOURCE_FILES = [
    "all/beta_binomial_within_slopes.csv",
    "direct/beta_binomial_within_slopes.csv",
    "all/beta_binomial_within_omnibus.csv",
    "direct/beta_binomial_within_omnibus.csv",
    "all/h2_decomposition_models.csv",
    "direct/h2_decomposition_models.csv",
    "h3_original_corrected_comparison.json",
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


def build_h1_atomic(results: Path) -> pd.DataFrame:
    pieces = []
    for scope, rel in (
        ("all_analysis_eligible", "all/beta_binomial_within_slopes.csv"),
        ("direct_only", "direct/beta_binomial_within_slopes.csv"),
    ):
        frame = pd.read_csv(results / rel)
        frame.insert(0, "evidence_scope", scope)
        pieces.append(frame)
    out = pd.concat(pieces, ignore_index=True)
    out = out[
        [
            "evidence_scope",
            "stratum",
            "context",
            "outcome",
            "n_islands",
            "n_species_trials",
            "geography_slope_log_odds",
            "cluster_robust_se",
            "p_value",
            "kappa",
            "optimizer_success",
        ]
    ]
    out["submission_note"] = ""
    mask = (
        out["evidence_scope"].eq("direct_only")
        & out["stratum"].eq("all_observed")
        & out["context"].eq("northern_high_latitude")
        & out["outcome"].eq("shallow_open_tube")
    )
    out.loc[mask, "submission_note"] = (
        "Frozen optimizer flag audited separately; converged retry changed slope "
        "by 2.93e-06 and the six-response sensitivity remained supported."
    )
    return out


def build_h1_joint(results: Path) -> pd.DataFrame:
    pieces = []
    for scope, rel in (
        ("all_analysis_eligible", "all/beta_binomial_within_omnibus.csv"),
        ("direct_only", "direct/beta_binomial_within_omnibus.csv"),
    ):
        frame = pd.read_csv(results / rel)
        frame.insert(0, "evidence_scope", scope)
        pieces.append(frame)
    out = pd.concat(pieces, ignore_index=True)
    return out[
        [
            "evidence_scope",
            "stratum",
            "context",
            "status",
            "n_retained_outcomes",
            "retained_outcomes",
            "n_unique_islands",
            "n_clusters",
            "joint_wald_chisq",
            "joint_df",
            "p_value",
            "q_value",
            "all_optimizers_converged",
            "vector_supported",
        ]
    ]


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
        "Table_S2a_H1_atomic.csv": build_h1_atomic(results),
        "Table_S2b_H1_joint.csv": build_h1_joint(results),
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

    source_meta = {
        rel: _sha256(results / rel)
        for rel in SOURCE_FILES
    }
    manifest = {
        "contract": "chapter1_submission_supplement_tables_v1",
        "scientific_surface": "corrected_geography_20260924",
        "sources": source_meta,
        "outputs": output_meta,
        "notes": {
            "Table_S2a": (
                "The frozen Direct-only northern-high-latitude shallow/open-tube "
                "optimizer flag is preserved and annotated; the separate convergence "
                "audit closes that numerical warning without rewriting frozen results."
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
