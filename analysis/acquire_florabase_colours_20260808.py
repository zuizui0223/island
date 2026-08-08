"""Acquire explicit flower-colour evidence from the official WA Florabase.

The script is intentionally resumable: every search/profile response is cached,
and reruns only request missing pages.  Only current, species-rank records whose
binomial exactly matches the strict universe are eligible.  A result is emitted
only when the taxon profile itself contains a floral description and an explicit
colour token; the search-result snippet is never used as evidence.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import math
import re
import time
from pathlib import Path
from urllib.parse import urlencode

import pandas as pd
import requests
from bs4 import BeautifulSoup

BASE = "https://florabase.dbca.wa.gov.au"
COLORS = (
    "black",
    "blue",
    "brown",
    "cream",
    "green",
    "grey",
    "orange",
    "pink",
    "purple",
    "red",
    "violet",
    "white",
    "yellow",
)
USER_AGENT = (
    "island-trait-research/1.0 "
    "(academic evidence acquisition; https://github.com/zuizui0223/island)"
)
COLOR_PATTERNS = {
    "black": r"\b(?:black|blackish)\b",
    "blue": r"\b(?:blue|bluish|azure)\b",
    "brown": r"\b(?:brown|brownish|maroon)\b",
    "cream": r"\b(?:cream|creamy)\b",
    "green": r"\b(?:green|greenish)\b",
    "grey": r"\b(?:grey|gray|greyish|grayish)\b",
    "orange": r"\b(?:orange|orange-red|orange-yellow)\b",
    "pink": r"\b(?:pink|pinkish|rose(?:-pink)?)\b",
    "purple": r"\b(?:purple|purplish|mauve|magenta)\b",
    "red": r"\b(?:red|reddish|scarlet|crimson)\b",
    "violet": r"\b(?:violet|violet-blue)\b",
    "white": r"\b(?:white|whitish)\b",
    "yellow": r"\b(?:yellow|yellowish|golden)\b",
}


def request(url: str, path: Path, retries: int = 4) -> tuple[str, str]:
    if path.exists() and path.stat().st_size > 1000:
        return "cached", path.read_text(encoding="utf-8", errors="replace")
    error = ""
    for attempt in range(retries):
        try:
            response = requests.get(url, headers={"User-Agent": USER_AGENT}, timeout=45)
            if response.status_code == 200 and len(response.content) > 1000:
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(response.content)
                return "http_200", response.text
            error = f"http_{response.status_code}"
        except requests.RequestException as exc:
            error = type(exc).__name__
        time.sleep(1.5 * (attempt + 1))
    return error or "request_failed", ""


def result_count(html: str) -> int:
    text = " ".join(BeautifulSoup(html, "html.parser").get_text(" ", strip=True).split())
    match = re.search(r"Displaying records \S+ of ([0-9,]+)", text)
    return int(match.group(1).replace(",", "")) if match else 0


def parse_search_page(html: str, color: str, page_url: str) -> list[dict]:
    soup = BeautifulSoup(html, "html.parser")
    rows: list[dict] = []
    for tr in soup.find_all("tr"):
        if "monomial" in (tr.get("class") or []):
            continue
        cells = tr.find_all("td")
        profile = tr.select_one('a[href^="/browse/profile/"]')
        if len(cells) < 2 or profile is None:
            continue
        italic = [" ".join(i.get_text(" ", strip=True).split()) for i in cells[1].find_all("i")]
        # Species rank only; infraspecific records have three italic components.
        if len(italic) != 2:
            continue
        species = " ".join(italic)
        rows.append(
            {
                "florabase_species": species,
                "colour_query": color,
                "profile_id": profile["href"].rstrip("/").split("/")[-1],
                "profile_url": BASE + profile["href"],
                "search_url": page_url,
            }
        )
    return rows


def profile_identity(soup: BeautifulSoup) -> str:
    heading = soup.find("h1")
    if heading is None:
        return ""
    italic = [" ".join(i.get_text(" ", strip=True).split()) for i in heading.find_all("i")]
    return " ".join(italic[:2]) if len(italic) >= 2 else ""


def extract_description(soup: BeautifulSoup) -> tuple[str, str]:
    block = soup.select_one("main#main-content .text-block.brief")
    if block is None:
        return "", ""
    paragraph = block.find("p")
    footer = block.find("footer")
    description = " ".join(paragraph.get_text(" ", strip=True).split()) if paragraph else ""
    attribution = " ".join(footer.get_text(" ", strip=True).split()) if footer else ""
    return description, attribution


def floral_context(description: str) -> bool:
    # Fail closed: colour words elsewhere may describe leaves, fruits, or soils.
    return bool(
        re.search(
            r"(?:\bfl\.|\bflowers?\b|\bcorolla\b|\bperianth\b|\bpetals?\b|\btepals?\b)",
            description,
            flags=re.IGNORECASE,
        )
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--strict", required=True)
    parser.add_argument("--master", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--workers", type=int, default=5)
    args = parser.parse_args()

    out = Path(args.output)
    search_cache = out / "raw_search"
    profile_cache = out / "raw_profile"
    out.mkdir(parents=True, exist_ok=True)

    strict = pd.read_csv(args.strict)
    strict_color = strict[(strict["axis"] == "flower_colour") & (strict["quality"].isna())].copy()
    universe = set(strict_color["accepted_species"].astype(str))
    master = pd.read_csv(
        args.master, usecols=["accepted_species", "family", "n_islands", "n_records"]
    )
    master = master.drop_duplicates("accepted_species")

    first_pages: dict[str, str] = {}
    jobs: list[tuple[str, int, str, Path]] = []
    for color in COLORS:
        url = f"{BASE}/search/advanced?{urlencode({'current': 'yes', 'colour': color})}"
        _, html = request(url, search_cache / f"{color}-0001.html")
        first_pages[color] = html
        count = result_count(html)
        for page in range(2, math.ceil(count / 20) + 1):
            page_url = url + f"&page={page}"
            jobs.append((color, page, page_url, search_cache / f"{color}-{page:04d}.html"))

    def fetch_search(job: tuple[str, int, str, Path]) -> tuple[str, int, str, str, str]:
        color, page, url, path = job
        status, html = request(url, path)
        return color, page, url, status, html

    pages: list[tuple[str, str, str]] = []
    for color, html in first_pages.items():
        pages.append((color, f"{BASE}/search/advanced?current=yes&colour={color}", html))
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        for i, (color, page, url, status, html) in enumerate(pool.map(fetch_search, jobs), 1):
            if status in {"cached", "http_200"}:
                pages.append((color, url, html))
            if i % 100 == 0 or i == len(jobs):
                print(f"search_pages={i}/{len(jobs)}", flush=True)

    search_rows = [row for color, url, html in pages for row in parse_search_page(html, color, url)]
    search_df = pd.DataFrame(search_rows).drop_duplicates(
        ["florabase_species", "colour_query", "profile_id"]
    )
    search_df.to_csv(out / "florabase_colour_index.csv.gz", index=False)

    targets = search_df[search_df["florabase_species"].isin(universe)].copy()
    grouped = (
        targets.groupby(["florabase_species", "profile_id", "profile_url"], as_index=False)
        .agg(
            indexed_colours=("colour_query", lambda x: "|".join(sorted(set(x)))),
            search_urls=("search_url", lambda x: "|".join(sorted(set(x)))),
        )
        .rename(columns={"florabase_species": "accepted_species"})
        .merge(master, on="accepted_species", how="left", validate="one_to_one")
        .sort_values(["n_islands", "accepted_species"], ascending=[False, True])
    )
    targets = grouped
    targets.to_csv(out / "acquisition_targets.csv", index=False)

    def fetch_profile(row: dict) -> dict:
        status, html = request(row["profile_url"], profile_cache / f"{row['profile_id']}.html")
        result = dict(row)
        result.update(
            {
                "fetch_status": status,
                "matched_page_name": "",
                "identity_exact": False,
                "description": "",
                "description_attribution": "",
                "explicit_colours": "",
                "supporting_excerpt": "",
                "content_sha256": hashlib.sha256(html.encode("utf-8")).hexdigest() if html else "",
            }
        )
        if not html:
            return result
        soup = BeautifulSoup(html, "html.parser")
        identity = profile_identity(soup)
        description, attribution = extract_description(soup)
        result["matched_page_name"] = identity
        result["identity_exact"] = identity == row["accepted_species"]
        result["description"] = description
        result["description_attribution"] = attribution
        if not result["identity_exact"] or not description or not floral_context(description):
            return result
        flower_clause_match = re.search(r"\bFl\.\s*([^.]*)", description, re.IGNORECASE)
        flower_clause = flower_clause_match.group(0) if flower_clause_match else ""
        indexed = row["indexed_colours"].split("|")
        explicit = [
            colour
            for colour in indexed
            if flower_clause and re.search(COLOR_PATTERNS[colour], flower_clause, re.IGNORECASE)
        ]
        result["explicit_colours"] = "|".join(explicit)
        if explicit:
            result["supporting_excerpt"] = flower_clause
        return result

    records = targets.to_dict("records")
    fetched: list[dict] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        for i, result in enumerate(pool.map(fetch_profile, records), 1):
            fetched.append(result)
            if i % 100 == 0 or i == len(records):
                print(f"profiles={i}/{len(records)}", flush=True)
    inventory = pd.DataFrame(fetched)
    inventory.to_csv(out / "fetch_inventory.csv.gz", index=False)
    accepted = inventory[inventory["identity_exact"] & inventory["explicit_colours"].astype(bool)]
    print(f"search_rows={len(search_df)}")
    print(f"unresolved_exact_targets={len(targets)}")
    print(f"accepted_species={accepted['accepted_species'].nunique()}")
    print(
        f"high_occurrence={accepted.loc[accepted['n_islands'] >= 11, 'accepted_species'].nunique()}"
    )


if __name__ == "__main__":
    main()
