import hashlib
from pathlib import Path

import pytest

from island_v2.corrected_submission import compare_csv_table, verify_files


def test_hash_gate_rejects_changed_or_missing_inputs(tmp_path):
    p = tmp_path / "a.csv"
    p.write_bytes(b"locked\n")
    manifest = {"a.csv": hashlib.sha256(p.read_bytes()).hexdigest()}
    verify_files(tmp_path, manifest)
    p.write_bytes(b"changed\n")
    with pytest.raises(ValueError, match="SHA256"):
        verify_files(tmp_path, manifest)
    p.unlink()
    with pytest.raises(FileNotFoundError):
        verify_files(tmp_path, manifest)


def test_hash_gate_rejects_path_escape(tmp_path):
    with pytest.raises(ValueError, match="outside"):
        verify_files(tmp_path, {"../file": "x"})



def test_replay_table_comparison_accepts_numeric_tolerance_and_rejects_content_changes(tmp_path):
    expected = tmp_path / "expected.csv"
    observed = tmp_path / "observed.csv"
    expected.write_text("id,value,label\na,1.0,x\nb,2.0,y\n", encoding="utf-8")
    observed.write_text("id,value,label\na,1.0000001,x\nb,2.0,y\n", encoding="utf-8")
    report = compare_csv_table(observed, expected)
    assert report["rows"] == 2
    assert report["numeric_columns_checked"] == ["value"]

    observed.write_text("id,value,label\na,1.0,x\nb,2.0,z\n", encoding="utf-8")
    with pytest.raises(ValueError, match="categorical mismatch"):
        compare_csv_table(observed, expected)


def test_replay_table_comparison_rejects_numeric_drift(tmp_path):
    expected = tmp_path / "expected.csv"
    observed = tmp_path / "observed.csv"
    expected.write_text("value\n1.0\n", encoding="utf-8")
    observed.write_text("value\n1.1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="numeric mismatch"):
        compare_csv_table(observed, expected)


def test_current_submission_selects_complete_corrected_results():
    import json

    import pandas as pd

    root = Path(__file__).resolve().parents[1]
    lock = json.loads((root / "config/chapter1_submission_current.json").read_text(encoding="utf8"))
    assert lock["status"] == "primary_submission_baseline" and lock["primary"]
    verify_files(root, lock["files"])
    result = root / lock["results_directory"]
    dist = pd.read_csv(result / "gshhg_spherical_distances_all.csv")
    assert len(dist) == 8264 and dist.island_id.is_unique
    assert dist.spherical_coast_distance_km.gt(0).all()
    zero = pd.read_csv(result / "formerly_zero_islands_recalculated.csv")
    assert len(zero) == 1113
    sites = pd.read_csv(result / "glopl_corrected_site_distances.csv")
    assert len(sites) == 1248 and sites.site_key.is_unique
    assert sites.spherical_distance_km.eq(0).sum() == 996
    for scope in ["all", "direct"]:
        h1 = pd.read_csv(result / scope / "beta_binomial_within_omnibus.csv")
        broad = h1[h1.stratum.eq("all_observed")]
        assert len(broad) == 4 and broad.q_value.lt(0.05).all()
        assert {"all_observed", "all_native", "native_nonendemic"} == set(h1.stratum)
        assert len(pd.read_csv(result / scope / "h2_decomposition_models.csv")) == 12
        assert len(list((result / scope / "raw_patterns").glob("raw_*.csv"))) == 4
    h3 = json.loads((result / "h3_original_corrected_comparison.json").read_text())["corrected"]
    assert h3["global_gradient"]["distance_slope"] == pytest.approx(0.09191043596699681)
    assert h3["sensitivities"]["supplemental_only"]["two_sided_p"] > 0.05
    h4 = pd.read_csv(result / "h4_exact_corrected.csv")
    assert len(h4) == 6
