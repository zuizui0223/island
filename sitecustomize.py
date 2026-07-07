"""Temporary compatibility for the compressed GBIF island-species snapshot.

Some manually triggered legacy workflows still request the historical
``island_species_occurrences.csv`` path. During the repository migration, only
that exact missing path is redirected to its gzip-compressed successor. No other
file name or pandas operation is changed.

Remove this shim once every legacy workflow explicitly requests the ``.csv.gz``
path.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

try:
    import pandas as _pandas
except ImportError:
    _pandas = None

if _pandas is not None:
    _original_read_csv = _pandas.read_csv

    def _read_csv_with_gbif_snapshot_fallback(filepath_or_buffer: Any, *args: Any, **kwargs: Any):
        if isinstance(filepath_or_buffer, (str, Path)):
            path = Path(filepath_or_buffer)
            if path.name == "island_species_occurrences.csv" and not path.exists():
                compressed = path.with_name(f"{path.name}.gz")
                if compressed.exists():
                    filepath_or_buffer = compressed
        return _original_read_csv(filepath_or_buffer, *args, **kwargs)

    _pandas.read_csv = _read_csv_with_gbif_snapshot_fallback
