from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from island_v2.chapter1_v9_figure3_p1_defense_submission import Figure3P1Error, _load_p1_inputs


def _write_inputs(tmp_path: Path, *, valid_permutations: int = 2000, any_ci_positive: bool = False):
    p1c = tmp_path / "p1c"
    p1d = tmp_path / "p1d"
    p1c.mkdir()
    p1d.mkdir()
    result = {
        "status": "true_genus_exceeds_matched_complexity_null",
        "valid_permutations": valid_permutations,
        "permutations_regenerated": False,
        "observed_primary_median_conditional_attenuation": 0.72,
        "one_sided_randomization_p_value": 0.029,
    }
    (p1c / "chapter1_p1c_matched_genus_null_result.json").write_text(json.dumps(result))
    pd.DataFrame({
        "permutation_id": range(2000),
        "valid": [True] * 2000,
        "primary_median_conditional_attenuation": [0.2] * 2000,
    }).to_csv(p1c / "p1c_matched_genus_null_permutations.csv.gz", index=False)
    rows = []
    for stratum in ["all_native", "native_nonendemic"]:
        for mode in ["geo50_climate10", "geo_k5", "geo_k10", "geo_k20"]:
            rows.append({
                "evidence_scope": "direct_only",
                "stratum": stratum,
                "source_mode": mode,
                "observed_family_to_genus_extra_attenuation": 0.5,
                "family_to_genus_extra_attenuation_ci_low": -0.1,
                "family_to_genus_extra_attenuation_ci_high": 0.9,
                "extra_attenuation_ci_above_zero": any_ci_positive,
            })
    pd.DataFrame(rows).to_csv(p1d / "p1b_paired_bootstrap_summary.csv", index=False)
    return p1c, p1d


def test_load_p1_inputs_accepts_frozen_contract(tmp_path: Path):
    p1c, p1d = _write_inputs(tmp_path)
    out = _load_p1_inputs(p1c, p1d)
    assert len(out["permutations"]) == 2000
    assert len(out["p1d"]) == 8


def test_load_p1_inputs_rejects_incomplete_permutations(tmp_path: Path):
    p1c, p1d = _write_inputs(tmp_path, valid_permutations=1999)
    with pytest.raises(Figure3P1Error):
        _load_p1_inputs(p1c, p1d)


def test_load_p1_inputs_rejects_changed_p1d_claim(tmp_path: Path):
    p1c, p1d = _write_inputs(tmp_path, any_ci_positive=True)
    with pytest.raises(Figure3P1Error):
        _load_p1_inputs(p1c, p1d)
