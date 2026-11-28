#!/usr/bin/env python3
"""Restore the pinned Chapter 1 original ZIP from its private Drive parts.

Standard library only; no network, upload, extraction or changes to input files.
Use an existing, authorized private output directory.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile
import zipfile

ORIGINAL = "10058653212-chapter1-progressive-analysis-34232450884.zip"
SIZE = 105740636
SHA = "b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999"
PARTS = [
    (ORIGINAL + ".part01", 67108864, "eeac1b25389cd99d46915e7df918acce65a4a172c2e4a529a051cd82b0acc28b"),
    (ORIGINAL + ".part02", 38631772, "f45cc430f265f5cb79b5797384562ab460bb1e1fdde4668f1c5fc1fc43ba0797"),
]

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1048576), b""):
            h.update(b)
    return h.hexdigest()


def verify(path, size, digest):
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"Missing regular file: {path.name}")
    if path.stat().st_size != size or sha256(path) != digest:
        raise ValueError(f"Size or SHA256 mismatch: {path.name}")


def restore(source, destination):
    m = json.loads((source / (ORIGINAL + ".parts.json")).read_text(encoding="utf-8"))
    if (m["original"], m["original_bytes"], m["original_sha256"]) != (ORIGINAL, SIZE, SHA):
        raise ValueError("Original archive pin mismatch")
    if [(p["name"], p["bytes"], p["sha256"]) for p in m["parts"]] != PARTS:
        raise ValueError("Part order or pin mismatch")
    paths = []
    for name, size, digest in PARTS:
        p = source / name
        # Some Drive downloads append .bin; the content pin is unchanged.
        if not p.exists():
            p = source / (name + ".bin")
        verify(p, size, digest)
        paths.append(p)
    output = destination / ORIGINAL
    if not output.exists():
        fd, tmp = tempfile.mkstemp(prefix=".island-restore-", dir=destination)
        tmp = Path(tmp)
        try:
            with os.fdopen(fd, "wb") as out:
                for p in paths:
                    with p.open("rb") as f:
                        for b in iter(lambda: f.read(1048576), b""):
                            out.write(b)
            verify(tmp, SIZE, SHA)
            with zipfile.ZipFile(tmp) as z:
                if z.testzip() is not None:
                    raise ValueError("ZIP CRC mismatch")
            # Never overwrite an existing archive, even if created concurrently.
            try:
                os.link(tmp, output)
            except FileExistsError:
                verify(output, SIZE, SHA)
        finally:
            tmp.unlink(missing_ok=True)
    verify(output, SIZE, SHA)
    with zipfile.ZipFile(output) as z:
        if z.testzip() is not None:
            raise ValueError("ZIP CRC mismatch")
    return {"original": ORIGINAL, "bytes": SIZE, "sha256": SHA,
            "parts_verified": len(paths), "crc_verified": True,
            "network_used": False, "input_files_changed": False}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    os.umask(0o077)
    src = args.source_dir.resolve(strict=True)
    dst = args.output_dir.resolve(strict=True)
    if not src.is_dir() or not dst.is_dir():
        raise ValueError("Both directories must already exist")
    print(json.dumps(restore(src, dst), indent=2))
