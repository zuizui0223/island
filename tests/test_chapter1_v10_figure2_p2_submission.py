from __future__ import annotations

from pathlib import Path

import pandas as pd

from island_v2.chapter1_v10_figure2_p2_submission import _profile_label, _read


def test_profile_label() -> None:
    assert _profile_label("direct_only", "native_nonendemic") == "Direct · NNE"
    assert _profile_label("all_analysis_eligible", "all_native") == "All · native"


def test_read_requires_file(tmp_path: Path) -> None:
    path = tmp_path / "x.csv"
    pd.DataFrame({"a": [1]}).to_csv(path, index=False)
    out = _read(tmp_path, "x.csv")
    assert out.loc[0, "a"] == 1
