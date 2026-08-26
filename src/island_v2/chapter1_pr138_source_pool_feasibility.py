from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path

TARGETS = ("bio1", "bio4", "bio5", "bio6", "bio12", "bio15", "bio18", "bio19", "elevation")


def _text(value: object) -> str:
    return "" if value is None else str(value).strip()


def _num(value: object) -> int | None:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def _load(path: Path) -> list[dict]:
    payload = json.loads(path.read_text())
    if not isinstance(payload, list) or not payload:
        raise ValueError(f"unexpected or empty GIFT payload: {path}")
    return payload


def _write_csv(path: Path, rows: list[dict], fallback: str) -> None:
    fields = sorted({key for row in rows for key in row}) if rows else [fallback]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def run(lists_path: Path, taxonomy_path: Path, regions_path: Path, raster_path: Path, misc_path: Path, output_dir: Path) -> dict:
    lists = _load(lists_path)
    taxonomy = _load(taxonomy_path)
    regions = _load(regions_path)
    raster = _load(raster_path)
    misc = _load(misc_path)

    candidates = []
    for row in lists:
        if "mainland" not in _text(row.get("entity_class")).lower():
            continue
        if _num(row.get("native_indicated")) != 1:
            continue
        if _num(row.get("restricted")) not in (0, None):
            continue
        candidates.append({key: row.get(key) for key in (
            "list_ID", "ref_ID", "entity_ID", "geo_entity", "entity_class",
            "suit_geo", "native_indicated", "restricted", "taxon_ID"
        )})

    entity_ids = {_num(row.get("entity_ID")) for row in candidates}
    entity_ids.discard(None)
    region_by_id = {_num(row.get("entity_ID")): row for row in regions}
    matched_regions = [region_by_id[entity_id] for entity_id in entity_ids if entity_id in region_by_id]

    env_rows = []
    for row in raster:
        haystack = " ".join(_text(row.get(key)) for key in ("dataset", "layer_name", "layer", "description")).lower()
        layer_name = _text(row.get("layer_name")).lower()
        matched: set[str] = set()
        for target in TARGETS:
            if re.search(rf"(?<![a-z0-9]){re.escape(target)}(?![a-z0-9])", haystack):
                matched.add(target)
            match = re.fullmatch(r"bio(\d+)", target)
            if match and re.search(rf"(?:_|bio)0?{int(match.group(1))}(?:$|_)", layer_name):
                matched.add(target)
        if matched:
            enriched = dict(row)
            enriched["target_inputs"] = "|".join(sorted(matched))
            env_rows.append(enriched)

    output_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(output_dir / "gift_mainland_native_candidate_lists.csv", candidates, "list_ID")
    _write_csv(output_dir / "gift_candidate_mainland_regions.csv", matched_regions, "entity_ID")
    _write_csv(output_dir / "gift_candidate_pca_environment_layers.csv", env_rows, "target_inputs")

    audit = {
        "contract": "pr138_source_pool_feasibility_v2",
        "n_list_rows": len(lists),
        "n_candidate_mainland_native_lists": len(candidates),
        "n_candidate_suit_geo_1": sum(_num(row.get("suit_geo")) == 1 for row in candidates),
        "n_unique_candidate_mainland_entities": len(entity_ids),
        "n_candidate_entities_with_region_metadata": len(matched_regions),
        "n_region_rows": len(regions),
        "n_env_raster_layers": len(raster),
        "n_env_misc_variables": len(misc),
        "pca_target_inputs": list(TARGETS),
        "pca_target_layers_found": sorted({target for row in env_rows for target in row["target_inputs"].split("|")}),
        "n_candidate_pca_environment_layer_rows": len(env_rows),
        "entity_class_counts": dict(Counter(_text(row.get("entity_class")) for row in lists)),
        "taxonomy_payload_rows": len(taxonomy),
        "outcome_blind": True,
        "current_island_trait_outcomes_used": False,
        "claim_boundary": "Metadata feasibility only; no source region or source weight is assigned to a current island.",
    }
    (output_dir / "source_pool_feasibility.json").write_text(json.dumps(audit, indent=2) + "\n")
    (output_dir / "README.md").write_text(
        "# PR138 source-pool feasibility audit\n\n"
        f"Mainland/native/unrestricted lists: **{len(candidates)}**\n\n"
        f"Unique candidate mainland entities: **{len(entity_ids)}**\n\n"
        f"Entities with region metadata: **{len(matched_regions)}**\n\n"
        "PCA target variables found: **" + ", ".join(audit["pca_target_layers_found"]) + "**\n\n"
        "No current island trait values or fitted ecological results are used.\n"
    )
    return audit


def main() -> None:
    parser = argparse.ArgumentParser()
    for name in ("lists", "taxonomy", "regions", "raster", "misc"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.lists, args.taxonomy, args.regions, args.raster, args.misc, args.output_dir), indent=2))


if __name__ == "__main__":
    main()
