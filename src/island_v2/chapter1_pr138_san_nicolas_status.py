from __future__ import annotations

import argparse
import json
import re
from difflib import SequenceMatcher
from html.parser import HTMLParser
from pathlib import Path

import pandas as pd

SCI_RE = re.compile(r"^\s*([A-Z][A-Za-z-]+)\s+([a-z][A-Za-z-]+)\b")
CEDROS_SCI_RE = re.compile(
    r"(?<![A-Za-z-])(\*?)\s*([A-Z][A-Za-z-]+)\s+([a-z][A-Za-z-]+)\b"
)
STATUS_TOKEN_RE = re.compile(r"^(N|E)(?:\d+)?$")
CEDROS_START_MARKER = "Appendix 2."
CEDROS_PROSE_STARTS = {
    "collected",
    "distribution",
    "endemic",
    "found",
    "growing",
    "grows",
    "inhabits",
    "isla",
    "known",
    "listed",
    "location",
    "observed",
    "occurs",
    "preliminary",
    "present",
    "reported",
    "species",
    "this",
    "widespread",
}
CEDROS_FUZZY_MIN_PART = 0.84
CEDROS_FUZZY_MIN_MEAN = 0.88
CEDROS_FUZZY_MIN_MARGIN = 0.06
SAN_NICOLAS_STATUS_OVERRIDES = {
    # Table E-1 marks this N1, but footnote 1 explicitly states it is introduced
    # to San Nicolas Island despite being native to California.
    "encelia californica": "introduced",
}


class _CCH2ChecklistHTMLParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.records: list[tuple[str, str]] = []
        self._container_depth = 0
        self._note_depth = 0
        self._capture_name = False
        self._capture_note = False
        self._name_parts: list[str] = []
        self._note_parts: list[str] = []

    @staticmethod
    def _classes(attrs: list[tuple[str, str | None]]) -> set[str]:
        values = dict(attrs).get("class") or ""
        return set(values.split())

    def handle_starttag(
        self, tag: str, attrs: list[tuple[str, str | None]]
    ) -> None:
        classes = self._classes(attrs)
        if tag == "div" and not self._container_depth and "taxon-container" in classes:
            self._container_depth = 1
            self._name_parts = []
            self._note_parts = []
            return
        if not self._container_depth:
            return
        if tag == "div":
            self._container_depth += 1
            if "note-div" in classes:
                self._capture_note = True
                self._note_depth = self._container_depth
        elif tag == "span" and "taxon-span" in classes:
            self._capture_name = True

    def handle_endtag(self, tag: str) -> None:
        if not self._container_depth:
            return
        if tag == "span" and self._capture_name:
            self._capture_name = False
        if tag != "div":
            return
        if self._capture_note and self._container_depth == self._note_depth:
            self._capture_note = False
            self._note_depth = 0
        self._container_depth -= 1
        if self._container_depth == 0:
            name = " ".join(" ".join(self._name_parts).split())
            note = " ".join(" ".join(self._note_parts).split())
            if name:
                self.records.append((name, note))

    def handle_data(self, data: str) -> None:
        if self._capture_name:
            self._name_parts.append(data)
        if self._capture_note:
            self._note_parts.append(data)


def species_key(name: str) -> str:
    parts = str(name).strip().lstrip("*\"").split()
    if len(parts) < 2:
        return ""
    return f"{parts[0]} {parts[1]}".lower()


def _collapse_rows(rows: list[dict[str, str]]) -> pd.DataFrame:
    if not rows:
        return pd.DataFrame(
            columns=["source_species", "species_key", "origin_status", "status_conflict"]
        )
    frame = pd.DataFrame(rows).drop_duplicates()
    out: list[dict[str, object]] = []
    for key, group in frame.groupby("species_key", sort=True):
        statuses = sorted(set(group["origin_status"]))
        conflict = len(statuses) != 1
        out.append(
            {
                "source_species": "|".join(sorted(set(group["source_species"]))),
                "species_key": key,
                "origin_status": statuses[0] if len(statuses) == 1 else "unresolved",
                "status_conflict": conflict,
            }
        )
    return pd.DataFrame(out)


def parse_nps_san_nicolas_text(text: str) -> pd.DataFrame:
    """Parse NPS San Nicolas Table E-1 where N=native and E=exotic.

    Source footnotes override table tokens where they explicitly change island origin.
    This prevents N1 for Encelia californica from being misread as island-native.
    """
    rows: list[dict[str, str]] = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        match = SCI_RE.match(line)
        if not match:
            continue
        statuses = [token for token in line.split()[2:] if STATUS_TOKEN_RE.match(token)]
        if not statuses:
            continue
        status_token = statuses[-1][0]
        source_name = f"{match.group(1)} {match.group(2)}"
        key = species_key(source_name)
        origin = "native" if status_token == "N" else "introduced"
        origin = SAN_NICOLAS_STATUS_OVERRIDES.get(key, origin)
        rows.append(
            {
                "source_species": source_name,
                "species_key": key,
                "origin_status": origin,
            }
        )
    return _collapse_rows(rows)


def parse_cedros_oberbauer_text(text: str) -> pd.DataFrame:
    """Parse Oberbauer (1993) Cedros Appendix 2.

    The source explicitly states that taxa denoted with an asterisk are presumably
    introduced. Unstarred species entries in the annotated island list are therefore
    treated as native for this status sensitivity. The paper mentions ``Appendix 2``
    earlier in the prose, so parsing deliberately starts at the *last* marker, which is
    the actual annotated checklist heading.
    """
    marker_index = text.rfind(CEDROS_START_MARKER)
    if marker_index < 0:
        marker_index = text.rfind("Appendix 2")
    if marker_index < 0:
        raise ValueError("Cedros Appendix 2 marker not found")
    appendix = text[marker_index:]
    rows: list[dict[str, str]] = []
    for raw_line in appendix.splitlines():
        # ``pdftotext -layout`` emits this four-column appendix as one physical
        # line. Searching the whole line is therefore necessary; matching only
        # its start silently retained the leftmost column and discarded the
        # other three.
        for match in CEDROS_SCI_RE.finditer(raw_line):
            genus = match.group(2)
            if genus.lower() in CEDROS_PROSE_STARTS:
                continue
            source_name = f"{genus} {match.group(3)}"
            rows.append(
                {
                    "source_species": source_name,
                    "species_key": species_key(source_name),
                    "origin_status": "introduced" if match.group(1) else "native",
                }
            )
    return _collapse_rows(rows)


def parse_cch2_san_nicolas_html(html: str) -> pd.DataFrame:
    """Parse the author-maintained CCH2 San Nicolas island checklist.

    The checklist displays structural species parents above infraspecific terminal
    taxa. Only terminal entries are status observations. Explicit ``Non-native``
    or island-introduction notes are introduced; the page states that otherwise
    unmarked checklist entries are presumed native. An explicitly unclear island
    status remains unresolved.
    """
    parser = _CCH2ChecklistHTMLParser()
    parser.feed(html)
    if not parser.records:
        raise ValueError("CCH2 San Nicolas checklist contains no taxon records")

    grouped: dict[str, list[tuple[str, str]]] = {}
    for name, note in parser.records:
        grouped.setdefault(species_key(name), []).append((name, note))

    rows: list[dict[str, str]] = []
    for key, records in grouped.items():
        terminal = [record for record in records if len(record[0].split()) > 2]
        for source_name, note in terminal or records:
            lower_note = note.lower()
            if "status on island unclear" in lower_note:
                origin = "unresolved"
            elif "non-native" in lower_note or "introduc" in lower_note:
                origin = "introduced"
            else:
                origin = "native"
            rows.append(
                {
                    "source_species": source_name,
                    "species_key": key,
                    "origin_status": origin,
                }
            )
    return _collapse_rows(rows)


def _name_similarity(source_key: str, target_key: str) -> tuple[float, float]:
    source_parts = source_key.split()
    target_parts = target_key.split()
    if len(source_parts) != 2 or len(target_parts) != 2:
        return 0.0, 0.0
    genus = SequenceMatcher(None, source_parts[0], target_parts[0]).ratio()
    epithet = SequenceMatcher(None, source_parts[1], target_parts[1]).ratio()
    return min(genus, epithet), (genus + epithet) / 2


def _unique_fuzzy_target(source_key: str, target_keys: list[str]) -> str:
    scored: list[tuple[float, float, str]] = []
    for target_key in target_keys:
        min_part, mean = _name_similarity(source_key, target_key)
        scored.append((mean, min_part, target_key))
    scored.sort(reverse=True)
    if not scored:
        return ""
    top_mean, top_min, top_key = scored[0]
    second_mean = scored[1][0] if len(scored) > 1 else 0.0
    if (
        top_min < CEDROS_FUZZY_MIN_PART
        or top_mean < CEDROS_FUZZY_MIN_MEAN
        or top_mean - second_mean < CEDROS_FUZZY_MIN_MARGIN
    ):
        return ""
    return top_key


def _map_status_to_flora(
    parsed_status: pd.DataFrame,
    target_keys: list[str],
    *,
    allow_fuzzy_name_match: bool,
) -> pd.DataFrame:
    target_set = set(target_keys)
    mapped: list[dict[str, object]] = []
    for record in parsed_status.to_dict(orient="records"):
        source_key = str(record["species_key"])
        target_key = source_key if source_key in target_set else ""
        method = "exact_species_key" if target_key else ""
        if not target_key and allow_fuzzy_name_match:
            target_key = _unique_fuzzy_target(source_key, target_keys)
            method = "ocr_fuzzy_unique" if target_key else ""
        if not target_key:
            continue
        mapped.append(
            {
                "species_key": target_key,
                "origin_status": record["origin_status"],
                "status_conflict": bool(record["status_conflict"]),
                "source_species": record["source_species"],
                "name_match_method": method,
            }
        )

    if not mapped:
        return pd.DataFrame(
            columns=[
                "species_key",
                "origin_status",
                "status_conflict",
                "source_species",
                "name_match_method",
            ]
        )

    out: list[dict[str, object]] = []
    for key, group in pd.DataFrame(mapped).groupby("species_key", sort=True):
        statuses = set(group["origin_status"])
        source_conflict = bool(group["status_conflict"].any())
        conflict = source_conflict or len(statuses) != 1 or "unresolved" in statuses
        out.append(
            {
                "species_key": key,
                "origin_status": next(iter(statuses)) if not conflict else "unresolved",
                "status_conflict": conflict,
                "source_species": "|".join(sorted(set(group["source_species"]))),
                "name_match_method": "|".join(
                    sorted(set(group["name_match_method"]))
                ),
            }
        )
    return pd.DataFrame(out)


def build_status_ledger(
    island_species: pd.DataFrame,
    parsed_status: pd.DataFrame,
    *,
    island_id: str,
    status_source: str,
    status_reference: str,
    evidence_prefix: str,
    allow_fuzzy_name_match: bool = False,
) -> pd.DataFrame:
    required = {"island_id", "accepted_species"}
    missing = required.difference(island_species.columns)
    if missing:
        raise ValueError(f"island species missing columns: {sorted(missing)}")

    flora = island_species.loc[
        island_species["island_id"].astype(str).eq(str(island_id)),
        ["island_id", "accepted_species"],
    ].drop_duplicates().copy()
    flora["species_key"] = flora["accepted_species"].map(species_key)
    source = _map_status_to_flora(
        parsed_status,
        sorted(set(flora["species_key"]) - {""}),
        allow_fuzzy_name_match=allow_fuzzy_name_match,
    )
    out = flora.merge(source, on="species_key", how="left")
    out["origin_status"] = out["origin_status"].fillna("unresolved")
    out["endemic_status"] = "unresolved"
    out["status_source"] = status_source
    out["status_reference"] = status_reference
    out["status_evidence_id"] = out.apply(
        lambda r: f"{evidence_prefix}:{r['species_key']}"
        if r["origin_status"] != "unresolved"
        else "",
        axis=1,
    )
    columns = [
        "island_id",
        "accepted_species",
        "origin_status",
        "endemic_status",
        "status_source",
        "status_reference",
        "status_evidence_id",
        "source_species",
        "name_match_method",
    ]
    return out[columns].sort_values("accepted_species").reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdftotext", required=True)
    parser.add_argument("--island-species", required=True)
    parser.add_argument("--island-id", required=True)
    parser.add_argument(
        "--source-kind",
        choices=["san_nicolas_nps", "san_nicolas_cch2", "cedros_oberbauer"],
        required=True,
    )
    parser.add_argument("--status-reference", required=True)
    parser.add_argument("--output-csv", required=True)
    parser.add_argument("--manifest-json", required=True)
    args = parser.parse_args()

    text = Path(args.pdftotext).read_text(encoding="utf-8", errors="replace")
    if args.source_kind == "san_nicolas_nps":
        parsed = parse_nps_san_nicolas_text(text)
        status_source = (
            "NPS San Nicolas Island Integrated Natural Resources Management Plan 2010, Table E-1"
        )
        evidence_prefix = "san-nicolas-2010"
    elif args.source_kind == "san_nicolas_cch2":
        parsed = parse_cch2_san_nicolas_html(text)
        status_source = (
            "CCH2 author-maintained Checklist of the Vascular Flora of San Nicolas Island"
        )
        evidence_prefix = "san-nicolas-cch2-2025"
    else:
        parsed = parse_cedros_oberbauer_text(text)
        status_source = (
            "Oberbauer 1993 Preliminary Annotated List of Vascular Plants of Isla de Cedros"
        )
        evidence_prefix = "cedros-oberbauer-1993"

    flora = pd.read_csv(args.island_species)
    ledger = build_status_ledger(
        flora,
        parsed,
        island_id=args.island_id,
        status_source=status_source,
        status_reference=args.status_reference,
        evidence_prefix=evidence_prefix,
        allow_fuzzy_name_match=args.source_kind == "cedros_oberbauer",
    )
    output = Path(args.output_csv)
    output.parent.mkdir(parents=True, exist_ok=True)
    ledger.to_csv(output, index=False)

    resolved = ledger["origin_status"].ne("unresolved")
    manifest = {
        "contract": "chapter1_pr138_source_backed_island_status_v2",
        "source_kind": args.source_kind,
        "island_id": args.island_id,
        "n_source_species_keys": len(parsed),
        "n_island_species": len(ledger),
        "n_resolved": int(resolved.sum()),
        "n_exact_name_matches": int(
            ledger["name_match_method"].fillna("").str.contains("exact_species_key").sum()
        ),
        "n_fuzzy_name_matches": int(
            ledger["name_match_method"].fillna("").str.contains("ocr_fuzzy_unique").sum()
        ),
        "n_native": int(ledger["origin_status"].eq("native").sum()),
        "n_introduced": int(ledger["origin_status"].eq("introduced").sum()),
        "n_unresolved": int(ledger["origin_status"].eq("unresolved").sum()),
        "status_policy": (
            "source-backed island checklist only; unmatched or conflicting names remain unresolved"
        ),
        "status_reference": args.status_reference,
    }
    manifest_path = Path(args.manifest_json)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
