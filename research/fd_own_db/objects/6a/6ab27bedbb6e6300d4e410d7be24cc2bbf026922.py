#!/usr/bin/env python3
"""Fetch only the PR321 public aggregate surfaces at an immutable source commit.

No raw taxon-to-trait ledger, individual species profile or proprietary
source-citation table is fetched or redistributed.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import urllib.request

COMMIT = "d4a5dafab5f0737ea6c1c1c0263fc20decf6035c"
BASE = f"https://raw.githubusercontent.com/zuizui0223/island/{COMMIT}/results"
SCOPES = ("all", "direct", "all_single_state", "direct_single_state")
RELATIVE_PATHS = [
    "floral_space_20261106/all/profile_coordinates.csv",
    "floral_space_pairs_20261106/pairwise_observed_and_null.csv.gz",
]
RELATIVE_PATHS += [
    f"floral_space_20261106/{scope}/{name}"
    for scope in SCOPES
    for name in ("island_profile_carriers.csv.gz", "island_metrics.csv")
]


def fetch(root: Path, *, max_bytes: int = 80 * 1024 * 1024) -> dict:
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    manifest = {"upstream_repository": "zuizui0223/island", "upstream_commit": COMMIT,
                "rights_boundary": "Public aggregate only; no name-level trait records",
                "files": {}}
    for relative in RELATIVE_PATHS:
        url = f"{BASE}/{relative}"
        destination = root / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        temp = destination.with_name(destination.name + ".download")
        digest = hashlib.sha256()
        size = 0
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "fd-floral-analysis/1.0"})
            with urllib.request.urlopen(request, timeout=90) as response, temp.open("wb") as handle:
                if response.status != 200:
                    raise RuntimeError(f"HTTP {response.status} for {relative}")
                while chunk := response.read(1 << 20):
                    size += len(chunk)
                    if size > max_bytes:
                        raise ValueError(f"Upstream aggregate exceeds size limit: {relative}")
                    digest.update(chunk)
                    handle.write(chunk)
            if not size:
                raise ValueError(f"Empty aggregate: {relative}")
            # Validate headers before atomic replacement; a 200 HTML error is not CSV.
            opener = gzip.open if destination.suffix == ".gz" else open
            with opener(temp, "rt", encoding="utf-8") as reader:
                header = reader.readline().rstrip("\r\n")
            if "," not in header or "island_id" not in header and "profile_id" not in header and "encoding" not in header:
                raise ValueError(f"Unexpected aggregate header: {relative}")
            temp.replace(destination)
        finally:
            if temp.exists():
                temp.unlink()
        manifest["files"][relative] = {"url": url, "sha256": digest.hexdigest(), "bytes": size}
    (root / "source_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("work/island_public"))
    args = parser.parse_args()
    manifest = fetch(args.output)
    print(f"Fetched {len(manifest['files'])} source-commit-pinned public aggregates")


if __name__ == "__main__":
    main()
