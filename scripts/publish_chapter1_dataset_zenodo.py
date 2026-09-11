from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from urllib.parse import quote

import requests


def require_ok(response: requests.Response, expected: set[int], context: str) -> None:
    if response.status_code not in expected:
        raise RuntimeError(
            f"{context} failed: HTTP {response.status_code}: {response.text[:2000]}"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package-dir", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--publish", action="store_true")
    parser.add_argument("--base-url", default="https://zenodo.org")
    args = parser.parse_args()

    token = os.environ.get("ZENODO_TOKEN", "").strip()
    if not token:
        raise RuntimeError("ZENODO_TOKEN is required")

    package_dir = args.package_dir.resolve()
    metadata_path = args.metadata.resolve()
    if not package_dir.is_dir():
        raise FileNotFoundError(package_dir)
    if not metadata_path.is_file():
        raise FileNotFoundError(metadata_path)

    files = sorted(p for p in package_dir.iterdir() if p.is_file())
    required = {
        "species_axis_coverage.public.csv.gz",
        "PUBLIC_SUBSET_MANIFEST.json",
        "README.md",
    }
    missing = sorted(required - {p.name for p in files})
    if missing:
        raise RuntimeError(f"Zenodo package missing required files: {missing}")

    payload = json.loads(metadata_path.read_text(encoding="utf-8"))
    if set(payload) != {"metadata"}:
        raise RuntimeError("metadata JSON must contain exactly one top-level 'metadata' object")
    metadata = payload["metadata"]
    if metadata.get("upload_type") != "dataset":
        raise RuntimeError("Zenodo upload_type must be dataset")
    if metadata.get("license") != "cc-by-sa-4.0":
        raise RuntimeError("Zenodo dataset licence must be cc-by-sa-4.0")

    headers = {"Authorization": f"Bearer {token}"}
    json_headers = {**headers, "Content-Type": "application/json"}
    api = args.base_url.rstrip("/") + "/api"

    create = requests.post(
        f"{api}/deposit/depositions",
        headers=json_headers,
        data=json.dumps({}),
        timeout=60,
    )
    require_ok(create, {201}, "create deposition")
    deposit = create.json()
    deposition_id = int(deposit["id"])
    bucket_url = deposit["links"]["bucket"]

    try:
        for path in files:
            with path.open("rb") as handle:
                upload = requests.put(
                    f"{bucket_url}/{quote(path.name)}",
                    headers=headers,
                    data=handle,
                    timeout=600,
                )
            require_ok(upload, {200, 201}, f"upload {path.name}")

        update = requests.put(
            f"{api}/deposit/depositions/{deposition_id}",
            headers=json_headers,
            data=json.dumps(payload),
            timeout=60,
        )
        require_ok(update, {200}, "update metadata")
        deposit = update.json()

        result = {
            "deposition_id": deposition_id,
            "draft_url": deposit.get("links", {}).get("html"),
            "reserved_doi": deposit.get("metadata", {}).get("prereserve_doi", {}).get("doi"),
            "published": False,
        }

        if args.publish:
            publish = requests.post(
                f"{api}/deposit/depositions/{deposition_id}/actions/publish",
                headers=headers,
                timeout=60,
            )
            require_ok(publish, {202}, "publish deposition")
            published = publish.json()
            result.update(
                {
                    "published": True,
                    "doi": published.get("doi")
                    or published.get("metadata", {}).get("doi"),
                    "doi_url": published.get("doi_url"),
                    "record_url": published.get("record_url")
                    or published.get("links", {}).get("record_html"),
                }
            )

        print(json.dumps(result, indent=2, sort_keys=True))
    except Exception:
        # Keep the draft for forensic inspection instead of deleting it silently.
        print(
            json.dumps(
                {
                    "deposition_id": deposition_id,
                    "draft_url": deposit.get("links", {}).get("html"),
                    "published": False,
                    "status": "failed_after_draft_creation",
                },
                indent=2,
                sort_keys=True,
            )
        )
        raise


if __name__ == "__main__":
    main()
