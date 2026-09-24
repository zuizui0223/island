"""Fail-closed provenance checks for the corrected Chapter 1 submission."""

import hashlib
from pathlib import Path


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
