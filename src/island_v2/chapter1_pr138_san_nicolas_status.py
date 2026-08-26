from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd

SCI_RE = re.compile(r"^\s*([A-Z][A-Za-z-]+)\s+([a-z][A-Za-z-]+)\b")
STATUS_TOKEN_RE = re.compile(r"^(N|E)(?:\d+)?$")
CEDROS_START_MARKER = "Appendix 2."
SAN_NICOLAS_STATUS_OVERRIDES = {
    # Table E-1 marks this N1, but footnote 1 explicitly states it is introduced
    # to San Nicolas Island despite being native to California.
    "encelia californica": "introduced",
}


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
    treated as native for this status sensitivity. Parsing begins only after Appendix 2
    so vegetation-community prose cannot create status records.
    """
    marker_index = text.find(CEDROS_START_MARKER)
    if marker_index < 0:
        marker_index = text.find("Appendix 2")
    if marker_index < 0:
        raise ValueError("Cedros Appendix 2 marker not found")
    appendix = text[marker_index:]
    rows: list[dict[str, str]] = []
    for raw_line in appendix.splitlines():
        stripped = raw_line.strip()
        introduced = stripped.startswith("*")
        line = stripped.lstrip("*\"").strip()
        match = SCI_RE.match(line)
        if not match:
            continue
        source_name = f"{match.group(1)} {match.group(2)}"
        if match.group(1) in {
            "Endemic",
            "Found",
            "Occurs",
            "Present",
            "Reported",
            "Known",
            "Listed",
            "Growing",
            "Grows",
            "Inhabits",
            "Widespread",
            "Preliminary",
            "Species",
            "This",
            "Isla",
            "Distribution",
            "Location",
        }:
            continue
        rows.append(
            {
                "source_species": source_name,
                "species_key": species_key(source_name),
                "origin_status": "introduced" if introduced else "native",
            }
        )
    return _collapse_rows(rows)


def build_status_ledger(
    island_species: pd.DataFrame,
    parsed_status: pd.DataFrame,
    *,
    island_id: str,
    status_source: str,
    status_reference: str,
    evidence_prefix: str,
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
    source = parsed_status[
        ["species_key", "origin_status", "status_conflict", "source_species"]
    ].copy()
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
    ]
    return out[columns].sort_values("accepted_species").reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdftotext", required=True)
    parser.add_argument("--island-species", required=True)
    parser.add_argument("--island-id", required=True)
    parser.add_argument(
        "--source-kind",
        choices=["san_nicolas_nps", "cedros_oberbauer"],
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
    )
    output = Path(args.output_csv)
    output.parent.mkdir(parents=True, exist_ok=True)
    ledger.to_csv(output, index=False)

    resolved = ledger["origin_status"].ne("unresolved")
    manifest = {
        "contract": "chapter1_pr138_source_backed_island_status_v2",
        "source_kind": args.source_kind,
        "island_id": args.island_id,
        "n_source_species_keys": int(len(parsed)),
        "n_island_species": int(len(ledger)),
        "n_resolved": int(resolved.sum()),
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
