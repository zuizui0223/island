from __future__ import annotations

import hashlib
import json
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
LOCK_PATH = ROOT / "config/chapter1_route_config_separation_lock.json"


def _git_blob_sha(path: Path) -> str:
    content = path.read_bytes()
    payload = f"blob {len(content)}\0".encode() + content
    return hashlib.sha1(payload, usedforsecurity=False).hexdigest()


def test_canonical_config_blobs_are_frozen_and_route_configs_are_separate() -> None:
    lock = json.loads(LOCK_PATH.read_text(encoding="utf-8"))

    for item in lock["canonical_configs"]:
        path = ROOT / item["path"]
        assert _git_blob_sha(path) == item["git_blob_sha"]
        config = yaml.safe_load(path.read_text(encoding="utf-8"))
        assert config["contract"] == item["contract"]

    for item in lock["all_data_route_configs"]:
        path = ROOT / item["path"]
        config = yaml.safe_load(path.read_text(encoding="utf-8"))
        assert config["contract"] == item["contract"]


def test_all_data_workflows_do_not_rebind_canonical_config_paths() -> None:
    lock = json.loads(LOCK_PATH.read_text(encoding="utf-8"))
    forbidden = {item["path"] for item in lock["canonical_configs"]}
    for relpath in lock["all_data_workflows"]:
        text = (ROOT / relpath).read_text(encoding="utf-8")
        for path in forbidden:
            assert path not in text, f"{relpath} still references frozen canonical path {path}"
