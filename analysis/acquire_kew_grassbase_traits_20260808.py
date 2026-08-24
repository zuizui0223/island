from __future__ import annotations

import argparse
import hashlib
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from threading import Lock

import pandas as pd
import requests
from bs4 import BeautifulSoup

INDEX_URL = "https://static1.kew.org/data/grasses-db/sppindex.htm"
USER_AGENT = (
    "island-trait-research/1.0 "
    "(academic citation-preserving evidence acquisition; "
    "https://github.com/zuizui0223/island)"
)


HEADINGS = {
    "HABIT",
    "INFLORESCENCE",
    "FERTILE SPIKELETS",
    "GLUMES",
    "FLORETS",
    "FLOWER",
    "FRUIT",
    "DISTRIBUTION",
    "NOTES",
}

INFLORESCENCE_PATTERNS = (
    (
        r"\b(?:racemes?|racemose|spikes?|spicate|panicles?|paniculate|"
        + r"thyrses?|thyrsoid|catkins?)\b",
        "raceme_spike_panicle",
    ),
    (r"\b(?:umbels?|umbellate|corymbs?|corymbose|cymes?|cymose)\b", "umbel_corymb"),
    (r"\b(?:capitulum|capitula|capitate|flower heads?)\b", "composite_display"),
    (r"\b(?:brush|puff)\b", "brush_puff_display"),
)


def text(value: object) -> str:
    return " ".join(str(value or "").replace("\xa0", " ").split())


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def parse_sections(page: bytes) -> tuple[str, dict[str, str], str]:
    soup = BeautifulSoup(page, "html.parser")
    heading = soup.find("h1", class_="scientificName")
    matched_name = text(heading.get_text(" ", strip=True) if heading else "")
    description = None
    for bold in soup.find_all(["b", "strong"]):
        if text(bold.get_text(" ", strip=True)).upper() == "HABIT":
            description = bold.find_parent("p")
            break
    if description is None:
        return matched_name, {}, ""

    description_text = text(description.get_text(" ", strip=True))
    heading_pattern = re.compile(
        r"\b("
        + "|".join(sorted((re.escape(item) for item in HEADINGS), key=len, reverse=True))
        + r")\b"
    )
    matches = list(heading_pattern.finditer(description_text))
    sections: dict[str, str] = {}
    for index, match in enumerate(matches):
        end = matches[index + 1].start() if index + 1 < len(matches) else len(description_text)
        sections[match.group(1)] = text(description_text[match.end() : end])
    return matched_name, sections, description_text


def inflorescence_value(excerpt: str) -> str:
    states = {
        state
        for pattern, state in INFLORESCENCE_PATTERNS
        if re.search(pattern, excerpt, re.IGNORECASE)
    }
    return "|".join(sorted(states))


def cleistogamy_value(description: str) -> str:
    if not re.search(r"\bcleistogam", description, re.IGNORECASE):
        return ""
    if re.search(r"\b(?:only|entirely|obligate(?:ly)?)\s+cleistogam", description, re.IGNORECASE):
        return "obligate"
    return "facultative"


def fetch_index(cache_path: Path) -> pd.DataFrame:
    if cache_path.exists():
        page = cache_path.read_bytes()
    else:
        response = requests.get(INDEX_URL, timeout=45, headers={"User-Agent": USER_AGENT})
        response.raise_for_status()
        page = response.content
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_bytes(page)
    soup = BeautifulSoup(page, "html.parser")
    rows = []
    for anchor in soup.select('a[href*="www/"][href$=".htm"]'):
        species = text(anchor.get_text(" ", strip=True))
        if len(species.split()) != 2:
            continue
        href = str(anchor.get("href", "")).lstrip("./")
        rows.append(
            {
                "grassbase_species": species,
                "source_url": "https://static1.kew.org/data/grasses-db/" + href,
            }
        )
    result = pd.DataFrame(rows).drop_duplicates("grassbase_species")
    if len(result) < 10_000:
        raise ValueError(f"GrassBase index returned only {len(result)} species")
    return result.sort_values("grassbase_species").reset_index(drop=True)


def fetch_one(row: dict[str, str], cache_dir: Path) -> dict[str, str]:
    url = row["source_url"]
    page_id = Path(url).stem
    cache = cache_dir / f"{page_id}.html"
    data = b""
    status = ""
    error = ""
    if cache.exists():
        data = cache.read_bytes()
        status = "cached"
    else:
        for attempt in range(3):
            try:
                response = requests.get(
                    url,
                    timeout=30,
                    headers={"User-Agent": USER_AGENT},
                )
                response.raise_for_status()
                data = response.content
                cache.write_bytes(data)
                status = f"http_{response.status_code}"
                break
            except Exception as exc:  # noqa: BLE001
                error = f"{type(exc).__name__}: {exc}"
                time.sleep(0.5 * (attempt + 1))
    matched_name, sections, description = parse_sections(data) if data else ("", {}, "")
    inflorescence = sections.get("INFLORESCENCE", "")
    return {
        **row,
        "matched_page_name": matched_name,
        "identity_exact": str(matched_name == row["grassbase_species"]),
        "fetch_status": status or "failed",
        "fetch_error": error,
        "content_sha256": sha256(data) if data else "",
        "inflorescence_excerpt": inflorescence,
        "inflorescence_value": inflorescence_value(inflorescence),
        "cleistogamy_excerpt": description if cleistogamy_value(description) else "",
        "cleistogamy_value": cleistogamy_value(description),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", type=Path)
    parser.add_argument("--strict", type=Path, required=True)
    parser.add_argument("--master", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--include-low", action="store_true")
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    cache_dir = args.output / "raw_html"
    cache_dir.mkdir(exist_ok=True)
    index = (
        pd.read_csv(args.index, dtype=str).fillna("")
        if args.index
        else fetch_index(args.output / "grassbase_species_index.html")
    )
    index.to_csv(args.output / "grassbase_index.csv", index=False)
    strict = pd.read_csv(args.strict, dtype=str).fillna("")
    master = pd.read_csv(args.master, dtype=str).fillna("")
    target_quality = {""} | ({"low"} if args.include_low else set())
    unresolved = strict.loc[
        strict["axis"].eq("floral_structural_complexity") & strict["quality"].isin(target_quality),
        ["accepted_species", "quality"],
    ].rename(columns={"quality": "baseline_quality"})
    targets = index.merge(
        unresolved,
        left_on="grassbase_species",
        right_on="accepted_species",
        how="inner",
        validate="one_to_one",
    ).merge(
        master[["accepted_species", "genus", "family", "n_islands", "n_records"]],
        on="accepted_species",
        how="left",
        validate="one_to_one",
    )
    targets["n_islands_numeric"] = pd.to_numeric(targets["n_islands"], errors="coerce").fillna(0)
    targets = targets.sort_values(
        ["n_islands_numeric", "n_records", "accepted_species"],
        ascending=[False, False, True],
        kind="stable",
    ).drop(columns="n_islands_numeric")
    if args.limit > 0:
        targets = targets.head(args.limit).copy()
    targets.to_csv(args.output / "acquisition_targets.csv", index=False)

    rows: list[dict[str, str]] = []
    lock = Lock()
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [
            pool.submit(fetch_one, row, cache_dir) for row in targets.astype(str).to_dict("records")
        ]
        for i, future in enumerate(as_completed(futures), start=1):
            result = future.result()
            with lock:
                rows.append(result)
                if i % 50 == 0 or i == len(futures):
                    pd.DataFrame(rows).to_csv(args.output / "fetch_inventory.csv.gz", index=False)
                    print(f"completed={i}/{len(futures)}", flush=True)
    result = pd.DataFrame(rows).sort_values("accepted_species", kind="stable")
    result.to_csv(args.output / "fetch_inventory.csv.gz", index=False)
    print(result[["fetch_status", "identity_exact"]].value_counts().to_string())
    print("inflorescence species", int(result["inflorescence_value"].ne("").sum()))
    print("cleistogamy species", int(result["cleistogamy_value"].ne("").sum()))


if __name__ == "__main__":
    main()
