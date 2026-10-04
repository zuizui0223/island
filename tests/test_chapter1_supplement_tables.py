from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "submission" / "build_chapter1_supplement_tables.py"
RESULTS = ROOT / "results" / "geography_20260924"
COMMITTED = ROOT / "submission" / "chapter1_current" / "supplement"

OUTPUTS = [
    "Table_S1_data_summary.csv",
    "Table_S2_H1_traitwise.csv",
    "Table_S3a_H2_all_finite_cluster.csv",
    "Table_S3b_H2_direct_finite_cluster.csv",
    "Table_S5_H3_pollen_limitation.csv",
    "Table_S6a_H4_scores.csv",
    "SUPPLEMENT_TABLES_MANIFEST.json",
]


def test_supplement_tables_rebuild_byte_for_byte(tmp_path: Path) -> None:
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--results",
            str(RESULTS),
            "--output",
            str(tmp_path),
        ],
        check=True,
    )
    for name in OUTPUTS:
        assert (tmp_path / name).read_bytes() == (COMMITTED / name).read_bytes(), name


def test_supplement_table_shapes_and_primary_values() -> None:
    manifest = json.loads(
        (COMMITTED / "SUPPLEMENT_TABLES_MANIFEST.json").read_text(encoding="utf-8")
    )
    assert manifest["contract"] == "chapter1_submission_supplement_tables_v1"
    assert (
        manifest["scientific_surface"]
        == "corrected_geography_20260924_plus_final_traitwise_H1_20261004"
    )

    expected_rows = {
        "Table_S1_data_summary.csv": 13,
        "Table_S2_H1_traitwise.csv": 112,
        "Table_S3a_H2_all_finite_cluster.csv": 12,
        "Table_S3b_H2_direct_finite_cluster.csv": 12,
        "Table_S5_H3_pollen_limitation.csv": 3,
        "Table_S6a_H4_scores.csv": 6,
    }
    for name, rows in expected_rows.items():
        assert manifest["outputs"][name]["rows"] == rows

    s1 = pd.read_csv(COMMITTED / "Table_S1_data_summary.csv")
    values = dict(zip(s1["metric"], s1["value"], strict=True))
    assert int(values["analysis_universe"]) == 8264
    assert int(values["broad_H1_union"]) == 4379
    assert int(values["accepted_angiosperm_species"]) == 106295
    assert int(values["resolved_species_axis_cells"]) == 222688
    assert int(values["effect_rows"]) == 2969
    assert int(values["sites"]) == 1248
    assert int(values["publications"]) == 919

    h1 = pd.read_csv(COMMITTED / "Table_S2_H1_traitwise.csv")
    assert len(h1) == 112
    assert h1["optimizer_success"].all()
    primary_sc = h1.loc[
        h1["evidence_scope"].eq("all")
        & h1["flora_scope"].eq("broad")
        & h1["outcome"].eq("self_compatibility")
    ]
    assert len(primary_sc) == 4
    assert (primary_sc["estimate"] > 0).all()
    assert (primary_sc["p_two_sided"] < 0.05).all()

    h2 = pd.read_csv(COMMITTED / "Table_S3a_H2_all_finite_cluster.csv")
    tropical_access = h2.loc[
        h2["context"].eq("tropical")
        & h2["response"].eq("generalized_accessible")
    ].iloc[0]
    assert float(tropical_access["primary_H2b_q_finite_cluster_t"]) == pytest.approx(
        0.019131173936551844
    )

    h3 = pd.read_csv(COMMITTED / "Table_S5_H3_pollen_limitation.csv")
    primary_h3 = h3.loc[h3["analysis"].eq("primary")].iloc[0]
    assert float(primary_h3["estimate"]) == pytest.approx(0.091910435967)
    assert float(primary_h3["finite_publication_t_two_sided_p"]) == pytest.approx(
        0.01594350716910782
    )

    h4 = pd.read_csv(COMMITTED / "Table_S6a_H4_scores.csv")
    selfing = h4.loc[
        h4["family"].eq("reproductive_assurance")
        & h4["analysis"].eq("primary")
    ].iloc[0]
    access = h4.loc[
        h4["family"].eq("accessibility_generalization")
        & h4["analysis"].eq("primary")
    ].iloc[0]
    assert float(selfing["estimate"]) == pytest.approx(-0.298300568155)
    assert float(selfing["finite_publication_t_two_sided_p"]) == pytest.approx(
        0.004166318114023793
    )
    assert float(access["estimate"]) == pytest.approx(-0.295659590051)
    assert float(access["finite_publication_t_two_sided_p"]) == pytest.approx(
        0.0233429161772711
    )
