"""Acquire the Goodwillie et al. 2010 Dryad files through ECU ScholarShip.

ECU's public record for the Dryad DOI is a DataRecord wrapper.  The actual
bitstreams live on the item referenced by ``relation.hasHubDataset``.  This
adapter follows that explicit relation, requires exact filenames, downloads
only public ORIGINAL bitstreams, and records the full HAL provenance chain.
It never guesses a bitstream UUID and never changes trait coverage itself.
"""
from __future__ import annotations

import argparse
import json
import urllib.request
from pathlib import Path

DATARECORD_ID = "bc729ea8-46f7-48c3-b013-af5cb93e9888"
SOURCE_DOI = "10.5061/dryad.bn5v3cm5"
TARGET_XLS = "floral display for outcrossing rate species.xls"
TARGET_README = "README_for_floral display for outcrossing rate species.txt"
XLS_MAGIC = bytes.fromhex("d0cf11e0a1b11ae1")


def get_json(url: str) -> dict:
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "island-source-scale-acquisition/1.0",
            "Accept": "application/hal+json, application/json",
        },
    )
    with urllib.request.urlopen(req, timeout=120) as response:
        return json.load(response)


def get_bytes(url: str) -> bytes:
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "island-source-scale-acquisition/1.0"},
    )
    with urllib.request.urlopen(req, timeout=120) as response:
        return response.read()


def embedded(payload: dict, key: str) -> list[dict]:
    value = payload.get("_embedded", {}).get(key, [])
    return value if isinstance(value, list) else []


def relation_values(item: dict, key: str) -> list[str]:
    rows = item.get("metadata", {}).get(key, [])
    return [str(row.get("value", "")).strip() for row in rows if str(row.get("value", "")).strip()]


def add_size(url: str, size: int) -> str:
    return url + ("&" if "?" in url else "?") + f"size={size}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--api-base", default="https://thescholarship.ecu.edu/server/api")
    parser.add_argument("--datarecord-id", default=DATARECORD_ID)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    base = args.api_base.rstrip("/")
    wrapper = get_json(f"{base}/core/items/{args.datarecord_id}")
    hubs = relation_values(wrapper, "relation.hasHubDataset")
    if len(hubs) != 1:
        raise ValueError(f"expected exactly one relation.hasHubDataset, got {hubs}")
    hub_id = hubs[0]
    hub = get_json(f"{base}/core/items/{hub_id}")

    bundles_url = (
        (hub.get("_links", {}).get("bundles", {}) or {}).get("href")
        or f"{base}/core/items/{hub_id}/bundles"
    )
    bundles_payload = get_json(add_size(bundles_url, 100))
    bundles = embedded(bundles_payload, "bundles")
    originals = [b for b in bundles if str(b.get("name", "")).casefold() == "original"]
    if len(originals) != 1:
        raise ValueError(f"expected exactly one ORIGINAL bundle on hub item, got {[b.get('name') for b in bundles]}")
    original = originals[0]
    bitstreams_url = (
        (original.get("_links", {}).get("bitstreams", {}) or {}).get("href")
        or f"{base}/core/bundles/{original.get('uuid') or original.get('id')}/bitstreams"
    )
    bits_payload = get_json(add_size(bitstreams_url, 200))
    bits = embedded(bits_payload, "bitstreams")

    targets = {TARGET_XLS.casefold(): None, TARGET_README.casefold(): None}
    for bit in bits:
        name = str(bit.get("name", "")).strip()
        key = name.casefold()
        if key not in targets:
            continue
        content = (bit.get("_links", {}).get("content", {}) or {}).get("href")
        if not content:
            uuid = bit.get("uuid") or bit.get("id")
            if uuid:
                content = f"{base}/core/bitstreams/{uuid}/content"
        if not content:
            raise ValueError(f"target bitstream lacks content URL: {name}")
        targets[key] = {"metadata": bit, "content_url": content}

    missing = [key for key, value in targets.items() if value is None]
    if missing:
        raise ValueError(
            f"target filenames absent from ORIGINAL bundle: {missing}; "
            f"available={[b.get('name') for b in bits]}"
        )

    xls_info = targets[TARGET_XLS.casefold()]
    readme_info = targets[TARGET_README.casefold()]
    assert xls_info is not None and readme_info is not None
    xls = get_bytes(str(xls_info["content_url"]))
    readme = get_bytes(str(readme_info["content_url"]))
    if not xls.startswith(XLS_MAGIC):
        raise ValueError(f"downloaded workbook is not OLE2 XLS; first8={xls[:8].hex()}")
    try:
        readme.decode("utf-8")
    except UnicodeDecodeError:
        readme.decode("latin-1")

    xls_path = args.output / "goodwillie2010_source.xls"
    readme_path = args.output / "goodwillie2010_README.txt"
    xls_path.write_bytes(xls)
    readme_path.write_bytes(readme)

    provenance = {
        "contract": "goodwillie2010_ecu_public_mirror_acquisition_v1",
        "source_doi": SOURCE_DOI,
        "ecu_datarecord_id": args.datarecord_id,
        "ecu_hub_dataset_id": hub_id,
        "wrapper_handle": wrapper.get("handle"),
        "hub_handle": hub.get("handle"),
        "relation_key": "relation.hasHubDataset",
        "original_bundle_id": original.get("uuid") or original.get("id"),
        "workbook": {
            "name": TARGET_XLS,
            "bitstream_id": xls_info["metadata"].get("uuid") or xls_info["metadata"].get("id"),
            "content_url": xls_info["content_url"],
            "size_bytes": len(xls),
        },
        "readme": {
            "name": TARGET_README,
            "bitstream_id": readme_info["metadata"].get("uuid") or readme_info["metadata"].get("id"),
            "content_url": readme_info["content_url"],
            "size_bytes": len(readme),
        },
        "wrapper_item": wrapper,
        "hub_item": hub,
        "original_bundle": original,
        "original_bitstreams": bits,
    }
    (args.output / "ecu_mirror_provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "ecu_hub_dataset_id": hub_id,
        "workbook_bytes": len(xls),
        "readme_bytes": len(readme),
        "workbook_bitstream_id": provenance["workbook"]["bitstream_id"],
        "readme_bitstream_id": provenance["readme"]["bitstream_id"],
    }, indent=2))


if __name__ == "__main__":
    main()
