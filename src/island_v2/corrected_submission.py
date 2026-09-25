"""Fail-closed provenance checks for the corrected Chapter 1 submission."""

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd


def verify_files(root: Path, files: dict[str, str]) -> None:
    root = root.resolve()
    if not files:
        raise ValueError("Empty input manifest")
    for name, expected in files.items():
        path = (root / name).resolve()
        if not path.is_relative_to(root):
            raise ValueError(f"Input outside declared root: {name}")
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1048576), b""):
                digest.update(block)
        if digest.hexdigest() != expected:
            raise ValueError(f"SHA256 mismatch: {name}")


def compare_csv_table(
    observed_path: Path,
    expected_path: Path,
    *,
    absolute_tolerance: float = 1e-6,
    relative_tolerance: float = 1e-5,
) -> dict[str, object]:
    """Compare one replayed table with its frozen counterpart."""

    observed = pd.read_csv(observed_path)
    expected = pd.read_csv(expected_path)
    if list(observed.columns) != list(expected.columns):
        raise ValueError(f"Replay columns differ: {observed_path.name}")
    if len(observed) != len(expected):
        raise ValueError(f"Replay row count differs: {observed_path.name}")

    numeric_columns: list[str] = []
    for column in expected.columns:
        left = observed[column]
        right = expected[column]
        if pd.api.types.is_numeric_dtype(left) and pd.api.types.is_numeric_dtype(right):
            x = pd.to_numeric(left, errors="coerce").to_numpy(float)
            y = pd.to_numeric(right, errors="coerce").to_numpy(float)
            if not np.allclose(
                x,
                y,
                atol=absolute_tolerance,
                rtol=relative_tolerance,
                equal_nan=True,
            ):
                delta = np.abs(x - y)
                finite = np.isfinite(delta)
                maximum = float(delta[finite].max()) if finite.any() else float("nan")
                raise ValueError(
                    f"Replay numeric mismatch: {observed_path.name}:{column}; max_abs={maximum}"
                )
            numeric_columns.append(column)
        else:
            x = left.astype("string").fillna("<NA>")
            y = right.astype("string").fillna("<NA>")
            if not x.equals(y):
                raise ValueError(f"Replay categorical mismatch: {observed_path.name}:{column}")

    return {
        "file": observed_path.as_posix(),
        "rows": int(len(observed)),
        "numeric_columns_checked": numeric_columns,
    }


def verify_replay_tables(
    observed_root: Path,
    expected_root: Path,
    relative_paths: list[str],
    *,
    absolute_tolerance: float = 1e-6,
    relative_tolerance: float = 1e-5,
) -> dict[str, object]:
    """Fail if any primary replay table differs from the frozen corrected result."""

    observed_root = observed_root.resolve()
    expected_root = expected_root.resolve()
    reports: list[dict[str, object]] = []
    for relative in relative_paths:
        observed = (observed_root / relative).resolve()
        expected = (expected_root / relative).resolve()
        if not observed.is_relative_to(observed_root) or not expected.is_relative_to(expected_root):
            raise ValueError(f"Replay table outside declared root: {relative}")
        report = compare_csv_table(
            observed,
            expected,
            absolute_tolerance=absolute_tolerance,
            relative_tolerance=relative_tolerance,
        )
        report["file"] = relative
        reports.append(report)
    return {
        "status": "pass",
        "absolute_tolerance": absolute_tolerance,
        "relative_tolerance": relative_tolerance,
        "tables": reports,
    }
