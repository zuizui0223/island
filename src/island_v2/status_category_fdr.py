"""Multiple-testing correction for status-aware category isolation tests."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _bh(p_values: pd.Series) -> pd.Series:
    p = pd.to_numeric(p_values, errors="coerce")
    valid = p.notna() & np.isfinite(p)
    out = pd.Series(np.nan, index=p.index, dtype=float)
    if not valid.any():
        return out
    values = p.loc[valid]
    order = values.sort_values().index
    ranked = values.loc[order].to_numpy(dtype=float)
    m = len(ranked)
    adjusted = ranked * m / np.arange(1, m + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    out.loc[order] = adjusted
    return out


def add_category_fdr(isolation: pd.DataFrame) -> pd.DataFrame:
    required = {"stratum", "regime", "domain", "category", "p_value"}
    missing = required - set(isolation.columns)
    if missing:
        raise typer.BadParameter(f"category isolation table missing columns: {sorted(missing)}")
    out = isolation.copy()
    out["q_value_bh"] = np.nan
    for _, idx in out.groupby(["stratum", "regime", "domain"], sort=True).groups.items():
        out.loc[idx, "q_value_bh"] = _bh(out.loc[idx, "p_value"])
    out["fdr_supported"] = out["q_value_bh"].lt(0.05)
    return out


@app.command("run")
def run(
    isolation_csv: Path = typer.Option(..., exists=True),
    output_csv: Path = typer.Option(...),
) -> None:
    result = add_category_fdr(pd.read_csv(isolation_csv))
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_csv, index=False)


if __name__ == "__main__":
    app()
