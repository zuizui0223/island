"""Audit the Chamaenerion/Epilobium fleischeri reproductive source chain.

Diagnostic only: no promotion. The audit resolves the exact Meyer workbook row,
then inspects the Gamba & Muchhala 2020 Dryad dataset and its source/reference
fields before deciding whether any strict reproductive trait is actually stated.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from urllib.parse import quote

import pandas as pd
import requests

TARGET_LABELS = {
    "Chamaenerion fleischeri",
    "Chamerion fleischeri",
    "Epilobium fleischeri",
}
EXPECTED_MEYER_SHA256 = "d6f30a8ad728a3b1b3c8ac5a5143fb935ffc64a38ce06708f6f2df4331737647"
GAMBA_DATASET_DOI = "10.5061/dryad.d2547d819"
GAMBA_ARTICLE_DOI = "10.1111/mec.15575"


def text(x: object) -> str:
    if x is None or pd.isna(x):
        return ""
    return " ".join(str(x).strip().split())


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def norm_species(x: object) -> str:
    return text(x).replace("_", " ")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--meyer-xlsx", type=Path, required=True)
    ap.add_argument("--coverage", type=Path, required=True)
    ap.add_argument("--direct-ledger", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    observed = sha256(args.meyer_xlsx)
    if observed != EXPECTED_MEYER_SHA256:
        raise ValueError(f"Meyer hash mismatch: {observed}")

    coverage = pd.read_csv(args.coverage, dtype=str).fillna("")
    rep = coverage.loc[coverage.axis.eq("reproductive_assurance")].copy()
    target_cov = rep.loc[rep.accepted_species.isin(TARGET_LABELS)].copy()
    target_cov.to_csv(args.output / "target_current_coverage.csv", index=False)

    direct = pd.read_csv(args.direct_ledger, dtype=str).fillna("")
    target_direct = direct.loc[direct.accepted_species.isin(TARGET_LABELS)].copy()
    target_direct.to_csv(args.output / "target_current_direct_rows.csv", index=False)

    book = pd.ExcelFile(args.meyer_xlsx, engine="openpyxl")
    target_rows = []
    source_blocks = []
    for sheet in book.sheet_names:
        frame = pd.read_excel(args.meyer_xlsx, sheet_name=sheet, dtype=object, engine="openpyxl").fillna("")
        lookup = {str(c).strip().casefold(): c for c in frame.columns}
        if "genus_species" not in lookup:
            continue
        spcol = lookup["genus_species"]
        mask = frame[spcol].map(norm_species).isin(TARGET_LABELS)
        if not mask.any():
            continue
        for idx, row in frame.loc[mask].iterrows():
            rec = {"sheet": sheet, "source_row": int(idx + 2)}
            rec.update({str(k): text(v) for k, v in row.to_dict().items()})
            target_rows.append(rec)
            refcol = lookup.get("ref")
            if refcol is not None:
                ref = text(row[refcol])
                block = frame.loc[frame[refcol].map(text).eq(ref)].copy()
                for j, brow in block.iterrows():
                    brec = {"sheet": sheet, "source_row": int(j + 2), "target_ref": ref}
                    brec.update({str(k): text(v) for k, v in brow.to_dict().items()})
                    source_blocks.append(brec)

    pd.DataFrame(target_rows).to_csv(args.output / "meyer_target_rows.csv", index=False)
    pd.DataFrame(source_blocks).to_csv(args.output / "meyer_target_source_blocks.csv", index=False)

    headers = {"User-Agent": "island-source-audit/1.0 (reproductive trait verification)"}
    dryad_meta = {}
    dryad_files = []
    dryad_hits = []
    dryad_error = ""
    try:
        encoded = quote(f"doi:{GAMBA_DATASET_DOI}", safe="")
        meta = requests.get(
            f"https://datadryad.org/api/v2/datasets/{encoded}",
            headers=headers, timeout=30,
        )
        meta.raise_for_status()
        dryad_meta = meta.json()
        version = (dryad_meta.get("_links", {}).get("stash:version", {}).get("href") or "").rstrip("/").split("/")[-1]
        files_url = f"https://datadryad.org/api/v2/versions/{version}/files" if version else ""
        if files_url:
            fr = requests.get(files_url, headers=headers, timeout=30)
            fr.raise_for_status()
            payload = fr.json()
            dryad_files = payload.get("_embedded", {}).get("stash:files", [])
            for f in dryad_files:
                name = text(f.get("path"))
                href = (f.get("_links", {}).get("stash:download", {}) or {}).get("href") or ""
                if href.startswith("/"):
                    href = "https://datadryad.org" + href
                if not href or not name.casefold().endswith(".csv"):
                    continue
                rr = requests.get(href, headers=headers, timeout=30)
                rr.raise_for_status()
                local = args.output / ("dryad_" + Path(name).name)
                local.write_bytes(rr.content)
                try:
                    df = pd.read_csv(local, dtype=str).fillna("")
                except Exception:
                    continue
                hitmask = pd.Series(False, index=df.index)
                for c in df.columns:
                    vals = df[c].astype(str).str.replace("_", " ", regex=False)
                    hitmask |= vals.isin(TARGET_LABELS)
                    hitmask |= vals.str.contains("fleischeri", case=False, na=False)
                if hitmask.any():
                    for idx, row in df.loc[hitmask].iterrows():
                        dryad_hits.append({"file": name, "row": int(idx + 2), **{str(k): text(v) for k, v in row.to_dict().items()}})
    except Exception as exc:
        dryad_error = repr(exc)

    (args.output / "gamba_dryad_metadata.json").write_text(json.dumps(dryad_meta, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (args.output / "gamba_dryad_files.json").write_text(json.dumps(dryad_files, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    pd.DataFrame(dryad_hits).to_csv(args.output / "gamba_dryad_target_hits.csv", index=False)

    summary = {
        "contract": "chamaenerion_fleischeri_source_chain_diagnostic_v1",
        "meyer_sha256": observed,
        "gamba_article_doi": GAMBA_ARTICLE_DOI,
        "gamba_dataset_doi": GAMBA_DATASET_DOI,
        "current_target_coverage_rows": int(len(target_cov)),
        "current_target_direct_rows": int(len(target_direct)),
        "meyer_target_rows": int(len(target_rows)),
        "meyer_target_refs": sorted({text(r.get("ref")) for r in target_rows if text(r.get("ref"))}),
        "meyer_target_mating_system_values": sorted({text(r.get("mating_system")) for r in target_rows if text(r.get("mating_system"))}),
        "meyer_target_quant_values": sorted({text(r.get("quant")) for r in target_rows if text(r.get("quant"))}),
        "dryad_files_seen": int(len(dryad_files)),
        "dryad_target_hits": int(len(dryad_hits)),
        "dryad_target_files": sorted({text(r.get("file")) for r in dryad_hits}),
        "dryad_error": dryad_error,
        "formal_gain": 0,
        "promotion_allowed": False,
        "next_gate": "inspect exact Gamba/Dryad mating-system field and its cited original source; do not recode mixed mating or selfing ability into SC unless compatibility is explicitly supported",
    }
    (args.output / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
