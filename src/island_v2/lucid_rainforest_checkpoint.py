"""Acquire strict floral evidence from CSIRO Australian Rainforest Plants.

The source exposes one public current-name index and static species fact sheets.
The index is intersected locally with the fixed strict universe, and only species
with at least one unresolved strict axis are fetched.  Completed pages are cached,
so reruns never request them again.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import hashlib
import json
import re
from pathlib import Path
from typing import Any

import pandas as pd
import requests
from bs4 import BeautifulSoup

from analysis.acquire_flora_of_australia_traits_20260809 import extract_description
from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS
from island_v2.trait_measurements import load_rules

BASE_URL = "https://apps.lucidcentral.org/rainforest/text/entities"
INDEX_URL = f"{BASE_URL}/index.htm?v=20260811"
SOURCE_RUN_ID = "csiro-rainforest-plants-source-scale-20260811"
SOURCE_ARTIFACT = "csiro-rainforest-plants-reviewed-source-package"
SOURCE_FILE = "csiro_rainforest_plants_evidence.csv.gz"
SOURCE_PROVIDER = "csiro_australian_tropical_rainforest_plants"
SOURCE_CITATION = "CSIRO Australian Tropical Rainforest Plants, online edition 8"
ACCEPTANCE_CONTRACT = "official_flora_exact_species_family_flower_section_v1"
USER_AGENT = "island-trait-research/1.0 (github.com/zuizui0223/island)"
AUDIT_SIZE = 200
AUDIT_SEED = "csiro-rainforest-audit-20260811"

AXIS = {
    "flower_primary_color": "flower_colour",
    "floral_symmetry": "floral_structural_complexity",
    "floral_form": "floral_structural_complexity",
    "flower_size_class": "floral_structural_complexity",
    "tube_depth_class": "floral_structural_complexity",
    "inflorescence_display": "floral_structural_complexity",
}


def _text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def _sha(value: str | bytes) -> str:
    payload = value if isinstance(value, bytes) else value.encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _write_gzip(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.GzipFile(filename=str(path), mode="wb", mtime=0) as handle:
        handle.write(payload)


def _read_gzip(path: Path) -> bytes:
    with gzip.open(path, "rb") as handle:
        return handle.read()


def _fetch(url: str, path: Path) -> tuple[bytes, bool]:
    if path.exists():
        return _read_gzip(path), False
    response = requests.get(url, headers={"User-Agent": USER_AGENT}, timeout=45)
    response.raise_for_status()
    _write_gzip(path, response.content)
    return response.content, True


def acquire_index(cache_dir: Path) -> tuple[pd.DataFrame, bool, str]:
    payload, fetched = _fetch(INDEX_URL, cache_dir / "index.htm.gz")
    soup = BeautifulSoup(payload, "html.parser")
    rows: list[dict[str, str]] = []
    for anchor in soup.find_all("a", href=True):
        href = _text(anchor.get("href"))
        name = _text(anchor.get_text(" ", strip=True))
        if not href.endswith(".htm") or "/" in href or len(name.split()) < 2:
            continue
        rows.append({"accepted_species": name, "href": href})
    frame = pd.DataFrame(rows).drop_duplicates(["accepted_species", "href"])
    return frame.reset_index(drop=True), fetched, _sha(payload)


def _page_path(cache_dir: Path, href: str) -> Path:
    return cache_dir / "pages" / f"{Path(href).stem}.html.gz"


def acquire_pages(
    selected: pd.DataFrame, cache_dir: Path, workers: int
) -> dict[str, int]:
    def one(row: Any) -> tuple[bool, bool]:
        try:
            _, fetched = _fetch(
                f"{BASE_URL}/{row.href}", _page_path(cache_dir, row.href)
            )
            return fetched, False
        except requests.RequestException:
            return False, True

    fetched = failed = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        for was_fetched, was_failed in executor.map(
            one, selected.itertuples(index=False)
        ):
            fetched += int(was_fetched)
            failed += int(was_failed)
    return {"selected": len(selected), "fetched": fetched, "failed": failed}


def _page_fields(payload: bytes) -> dict[str, str]:
    soup = BeautifulSoup(payload, "html.parser")
    title = _text(soup.select_one("h1.bannertitle").get_text(" ", strip=True)) if soup.select_one("h1.bannertitle") else ""
    family_node = soup.select_one("div#family div.content")
    flower_node = soup.select_one("div#flowers div.content")
    return {
        "page_title": title,
        "family": _text(family_node.get_text(" ", strip=True)) if family_node else "",
        "flower_text": _text(flower_node.get_text(" ", strip=True)) if flower_node else "",
    }


def build_evidence(
    selected: pd.DataFrame,
    master: pd.DataFrame,
    cache_dir: Path,
    baseline_pairs: set[tuple[str, str]],
    rules: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    master_family = dict(zip(master.accepted_species, master.family, strict=True))
    evidence_rows: list[dict[str, str]] = []
    inventory_rows: list[dict[str, Any]] = []
    for row in selected.sort_values("accepted_species").itertuples(index=False):
        path = _page_path(cache_dir, row.href)
        if not path.exists():
            continue
        payload = _read_gzip(path)
        fields = _page_fields(payload)
        species = _text(row.accepted_species)
        expected_family = _text(master_family.get(species))
        title_ok = bool(re.match(rf"^{re.escape(species)}(?:\s|$)", fields["page_title"]))
        family_ok = bool(expected_family and fields["family"] == expected_family)
        cultivar_excluded = bool(
            re.search(
                r"\b(?:cultivar|hybrid|grex)\b|[×✕]",
                fields["flower_text"],
                re.IGNORECASE,
            )
        )
        identity_ok = title_ok and family_ok
        extracted = (
            extract_description(fields["flower_text"], rules)
            if identity_ok and fields["flower_text"] and not cultivar_excluded
            else {}
        )
        novel = {
            trait: item
            for trait, item in extracted.items()
            if (species, trait) not in baseline_pairs
        }
        url = f"{BASE_URL}/{row.href}"
        lineage = f"provider_treatment:csiro_rainforest_key:{species}"
        for trait, item in novel.items():
            value = _text(item["value"])
            quote = _text(item["quote"])
            record_id = f"csiro-rainforest-{_sha(f'{species}|{trait}|{value}|{quote}')[:24]}"
            evidence_rows.append(
                {
                    "accepted_species": species,
                    "axis": AXIS[trait],
                    "trait_name": trait,
                    "normalized_value": value,
                    "quality": "high",
                    "source_group": "latest_public_web",
                    "source_provider": SOURCE_PROVIDER,
                    "source_url": url,
                    "source_record_id": record_id,
                    "source_citation": SOURCE_CITATION,
                    "source_excerpt": quote,
                    "evidence_scope": "species_direct",
                    "name_match_method": "accepted_name_exact_plus_family_exact",
                    "source_lineage": lineage,
                    "lineage_method": "official_provider_species_treatment",
                    "source_run_id": SOURCE_RUN_ID,
                    "source_artifact": SOURCE_ARTIFACT,
                    "source_file": SOURCE_FILE,
                    "acceptance_contract": ACCEPTANCE_CONTRACT,
                }
            )
        inventory_rows.append(
            {
                "accepted_species": species,
                "href": row.href,
                "source_url": url,
                "page_title": fields["page_title"],
                "page_family": fields["family"],
                "master_family": expected_family,
                "identity_gate_passed": identity_ok,
                "cultivar_excluded": cultivar_excluded,
                "flower_text_chars": len(fields["flower_text"]),
                "all_extracted_traits": "|".join(sorted(extracted)),
                "novel_extracted_traits": "|".join(sorted(novel)),
                "content_sha256": _sha(payload),
                "source_lineage": lineage,
            }
        )
    evidence = pd.DataFrame(evidence_rows, columns=EVIDENCE_COLUMNS)
    if not evidence.empty:
        evidence = evidence.drop_duplicates("source_record_id").sort_values(
            ["trait_name", "accepted_species", "source_record_id"], kind="stable"
        )
    return evidence.reset_index(drop=True), pd.DataFrame(inventory_rows)


def deterministic_audit_queue(evidence: pd.DataFrame) -> pd.DataFrame:
    if len(evidence) < AUDIT_SIZE:
        raise ValueError(f"only {len(evidence)} candidates available for 200-row audit")
    pool = evidence.copy()
    pool["_audit_hash"] = pool.source_record_id.map(
        lambda value: _sha(f"{AUDIT_SEED}|{value}")
    )
    selected = []
    for _, group in pool.groupby("trait_name", sort=True):
        selected.append(group.sort_values("_audit_hash").head(min(25, len(group))))
    audit = pd.concat(selected).drop_duplicates("source_record_id")
    remaining = pool.loc[~pool.index.isin(audit.index)].sort_values("_audit_hash")
    if len(audit) > AUDIT_SIZE:
        audit = audit.sort_values("_audit_hash").head(AUDIT_SIZE)
    elif len(audit) < AUDIT_SIZE:
        audit = pd.concat([audit, remaining.head(AUDIT_SIZE - len(audit))])
    audit = audit.sort_values(["trait_name", "_audit_hash"]).rename(
        columns={
            "source_record_id": "candidate_id",
            "source_provider": "provider",
            "source_excerpt": "exact_supporting_quote",
        }
    )
    audit["accepted_correct"] = ""
    audit["cultivar_status"] = "wild_or_naturalised_official_flora_treatment"
    audit["reviewer"] = ""
    audit["reviewed_at_utc"] = ""
    audit["audit_reason"] = ""
    return audit[
        [
            "candidate_id",
            "accepted_species",
            "trait_name",
            "normalized_value",
            "provider",
            "source_url",
            "source_lineage",
            "source_citation",
            "exact_supporting_quote",
            "name_match_method",
            "cultivar_status",
            "accepted_correct",
            "reviewer",
            "reviewed_at_utc",
            "audit_reason",
        ]
    ].reset_index(drop=True)


def _hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run(
    *,
    master_csv: Path,
    strict_coverage_csv: Path,
    baseline_evidence_csv: Path,
    cache_dir: Path,
    output_dir: Path,
    measurement_config: Path,
    workers: int,
) -> dict[str, Any]:
    strict = pd.read_csv(strict_coverage_csv, dtype=str).fillna("")
    universe_names = set(strict.accepted_species)
    if len(universe_names) != 106_295:
        raise ValueError("coverage universe must contain exactly 106,295 species")
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    master = master.loc[master.accepted_species.isin(universe_names)].copy()
    if len(master) != 106_295 or master.accepted_species.nunique() != 106_295:
        raise ValueError("master does not uniquely cover the strict universe")
    baseline = pd.read_csv(baseline_evidence_csv, dtype=str).fillna("")
    baseline_pairs = set(zip(baseline.accepted_species, baseline.trait_name, strict=True))
    unresolved_names = set(strict.loc[strict.quality.eq(""), "accepted_species"])
    index, index_fetched, index_hash = acquire_index(cache_dir)
    selected = index.loc[
        index.accepted_species.isin(universe_names & unresolved_names)
    ].drop_duplicates("accepted_species")
    fetch = acquire_pages(selected, cache_dir, workers)
    evidence, inventory = build_evidence(
        selected,
        master,
        cache_dir,
        baseline_pairs,
        load_rules(measurement_config),
    )
    audit = deterministic_audit_queue(evidence)
    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / SOURCE_FILE
    inventory_path = output_dir / "csiro_rainforest_plants_treatment_inventory.csv.gz"
    audit_path = output_dir / "csiro_rainforest_plants_audit_queue_200.csv"
    index_path = output_dir / "csiro_rainforest_plants_index.csv.gz"
    deterministic_gzip = {"method": "gzip", "mtime": 0}
    evidence.to_csv(evidence_path, index=False, compression=deterministic_gzip)
    inventory.to_csv(inventory_path, index=False, compression=deterministic_gzip)
    audit.to_csv(audit_path, index=False)
    index.to_csv(index_path, index=False, compression=deterministic_gzip)
    summary = {
        "contract": "csiro_rainforest_plants_source_scale_acquisition_v1",
        "source_run_id": SOURCE_RUN_ID,
        "index_url": INDEX_URL,
        "index_sha256": index_hash,
        "index_fetched_this_run": index_fetched,
        "index_species": int(index.accepted_species.nunique()),
        "strict_exact_species": int(index.accepted_species.isin(universe_names).sum()),
        "selected_unresolved_species": len(selected),
        "fetch": fetch,
        "identity_gate_passed": int(inventory.identity_gate_passed.sum()),
        "pages_with_flower_text": int(inventory.flower_text_chars.gt(0).sum()),
        "baseline_species_trait_pairs": len(baseline_pairs),
        "novel_candidate_rows": len(evidence),
        "novel_species_trait": int(
            evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
        ),
        "novel_species_axis": int(
            evidence[["accepted_species", "axis"]].drop_duplicates().shape[0]
        ),
        "unique_species": int(evidence.accepted_species.nunique()),
        "by_trait": evidence.groupby("trait_name").size().sort_values(ascending=False).to_dict(),
        "search_queries": 0,
        "api_cost_usd": 0.0,
        "audit_rows": len(audit),
        "audit_seed": AUDIT_SEED,
        "guardrails": {
            "accepted_name_exact": True,
            "family_exact": True,
            "flower_section_only": True,
            "snippet_as_evidence": False,
            "cultivar_excluded": True,
            "family_inference": False,
            "global_fallback": False,
        },
    }
    (output_dir / "csiro_rainforest_plants_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
    )
    manifest = []
    for path in sorted(output_dir.iterdir()):
        if path.is_file() and path.name != "file_manifest.sha256":
            manifest.append(f"{_hash_file(path)}  {path.name}")
    (output_dir / "file_manifest.sha256").write_text(
        "\n".join(manifest) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--master-csv",
        type=Path,
        default=Path("data/v2/staging/gbif/collected/island_taxa.csv"),
    )
    parser.add_argument("--strict-coverage-csv", type=Path, required=True)
    parser.add_argument("--baseline-evidence-csv", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--measurement-config",
        type=Path,
        default=Path("config/measurement_classification.yml"),
    )
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    if not 1 <= args.workers <= 16:
        raise ValueError("workers must be between 1 and 16")
    print(json.dumps(run(**vars(args)), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
