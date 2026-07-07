"""Resolve the repository's cumulative GBIF island-species snapshot path.

The snapshot is gzip-compressed to keep GitHub commits below large-file limits.
A legacy uncompressed file is accepted only during one-time migration.
"""

from __future__ import annotations

from pathlib import Path

COMPRESSED_NAME = "island_species_occurrences.csv.gz"
LEGACY_NAME = "island_species_occurrences.csv"


def resolve_snapshot(directory: str | Path) -> Path:
    directory = Path(directory)
    compressed = directory / COMPRESSED_NAME
    if compressed.exists():
        return compressed
    legacy = directory / LEGACY_NAME
    if legacy.exists():
        return legacy
    raise FileNotFoundError(
        f"No GBIF island-species snapshot found in {directory}; expected {COMPRESSED_NAME}."
    )
