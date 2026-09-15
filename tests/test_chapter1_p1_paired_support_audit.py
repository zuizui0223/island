import json
from pathlib import Path

import pandas as pd
import yaml

from island_v2.chapter1_p1_paired_support_audit import run_audit


def _write_scope(root: Path, scope: str) -> None:
    out = root / "taxonomic-depth" / scope
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for stratum in ("all_native", "native_nonendemic"):
        for source_mode in ("geo_k5", "geo_k10", "geo_k20", "geo50_climate10"):
            for syndrome in ("generalized_accessible", "selfing_core"):
                rows.append(
                    {
                        "island_id": "1",
                        "syndrome": syndrome,
                        "stratum": stratum,
                        "source_mode": source_mode,
                        "observed_score": 0.4,
                        "family_expected": 0.1,
                        "genus_expected": 0.3,
                        "after_family_residual": 0.3,
                        "after_genus_residual": 0.1,
                        "n_species": 10,
                    }
                )
    decomp = pd.DataFrame(rows)
    decomp.to_csv(out / "decomposition.csv", index=False)
    long = decomp.melt(
        id_vars=["island_id", "syndrome", "stratum", "source_mode", "n_species"],
        value_vars=["observed_score", "after_family_residual", "after_genus_residual"],
        var_name="taxonomic_stage",
        value_name="syndrome_score",
    )
    long["syndrome"] = long["taxonomic_stage"] + "__" + long["syndrome"]
    long.to_csv(out / "long_scores.csv", index=False)


def test_paired_support_audit_passes_common_support(tmp_path: Path) -> None:
    artifact = tmp_path / "artifact"
    _write_scope(artifact, "direct_only")
    _write_scope(artifact, "all_analysis_eligible")
    config = {
        "contract": "chapter1_p1_assembly_depth_defense_v1",
        "pinned_primary_artifact": {
            "workflow_run_id": 1,
            "artifact_id": 2,
            "artifact_digest": "sha256:test",
        },
        "primary_context": {
            "response_axes": ["generalized_accessible", "selfing_core"],
            "strata": ["all_native", "native_nonendemic"],
            "source_modes": ["geo_k5", "geo_k10", "geo_k20", "geo50_climate10"],
        },
    }
    config_path = tmp_path / "config.yml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    out = tmp_path / "out"
    manifest = run_audit(artifact, config_path, out)
    assert manifest["primary_direct_only_pass"] is True
    saved = json.loads((out / "chapter1_p1a_paired_support_manifest.json").read_text())
    assert saved["new_biological_model_fitted"] is False
    cells = pd.read_csv(out / "p1a_paired_support_cells.csv")
    assert len(cells) == 32
    assert cells["same_support_all_three_stages"].all()
