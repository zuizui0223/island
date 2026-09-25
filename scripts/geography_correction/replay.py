"""Replay all corrected models from hash-verified, frozen parent inputs."""

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

from island_v2.corrected_submission import verify_files, verify_replay_tables


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[2]
    manifest = json.loads(
        (repo / "config/chapter1_corrected_input_files_20260924.json").read_text(encoding="utf8")
    )
    contract = json.loads(
        (repo / "config/chapter1_submission_current.json").read_text(encoding="utf8")
    )
    verify_files(args.sources, manifest["files"])
    verify_files(repo, contract["files"])
    if args.verify_only:
        print("Locked inputs and submission results verified. No models run.")
        return
    out = args.output.resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError(
            "Replay output must be new or empty; committed results cannot be overwritten"
        )
    out.mkdir(parents=True, exist_ok=True)
    geometry = repo / contract["results_directory"]
    shutil.copy2(geometry / "glopl_corrected_site_distances.csv", out)
    env = os.environ.copy()
    env.update(
        ISLAND_CORRECTION_SOURCES=str(args.sources.resolve()),
        ISLAND_CORRECTION_OUTPUT=str(out),
        ISLAND_CORRECTION_GEOMETRY=str(geometry),
        PYTHONPATH=str(repo / "src"),
        PYTHONIOENCODING="utf-8",
    )
    for stage in [
        "reproduce_h2",
        "refit_corrected_geography",
        "refit_glopl",
        "refit_raw_flower_patterns",
        "refit_atomic",
    ]:
        subprocess.run(
            [sys.executable, str(Path(__file__).parent / (stage + ".py"))], env=env, check=True
        )
    replay_tables = [
        "corrected_geography_covariates.csv",
        "glopl_corrected_site_distances.csv",
        "h2_original_reproduced_all.csv",
        "h2_original_reproduced_direct.csv",
        "h4_atomic_corrected.csv",
        "h4_atomic_original_reproduced.csv",
        "h4_exact_corrected.csv",
        "h4_exact_original_reproduced.csv",
        "all/beta_binomial_between_omnibus.csv",
        "all/beta_binomial_between_slopes.csv",
        "all/beta_binomial_within_omnibus.csv",
        "all/beta_binomial_within_slopes.csv",
        "all/h2_decomposition_models.csv",
        "all/raw_patterns/raw_colour_architecture_model_results.csv",
        "all/raw_patterns/raw_colour_conditioned_architecture_model_results.csv",
        "all/raw_patterns/raw_colour_joint_omnibus.csv",
        "all/raw_patterns/raw_colour_model_results.csv",
        "direct/beta_binomial_between_omnibus.csv",
        "direct/beta_binomial_between_slopes.csv",
        "direct/beta_binomial_within_omnibus.csv",
        "direct/beta_binomial_within_slopes.csv",
        "direct/h2_decomposition_models.csv",
        "direct/raw_patterns/raw_colour_architecture_model_results.csv",
        "direct/raw_patterns/raw_colour_conditioned_architecture_model_results.csv",
        "direct/raw_patterns/raw_colour_joint_omnibus.csv",
        "direct/raw_patterns/raw_colour_model_results.csv",
    ]
    verification = verify_replay_tables(out, geometry, replay_tables)
    (out / "repository_replay_verification.json").write_text(
        json.dumps(verification, indent=2) + "\n",
        encoding="utf8",
    )
    (out / "REPLAY_COMPLETE.json").write_text(
        json.dumps({"status": "all_stages_completed_and_verified", "contract": contract["contract"]}, indent=2)
        + "\n",
        encoding="utf8",
    )


if __name__ == "__main__":
    main()
