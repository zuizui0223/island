"""Build a strict self-compatibility source packet for Justicia betonica.

Evidence contract:
- The pinned Meyer workbook must contain exactly one Layek et al. (2022) row for
  the exact fixed-universe species label, with raw mating_system=SC.
- The species must still be an unresolved reproductive-assurance cell and must
  not already contain direct self-incompatibility evidence.
- The original Plant Species Biology article explicitly reports the species as
  self-compatible from a study of floral biology, breeding system and pollination
  ecology (DOI 10.1111/1442-1984.12380).

No genus/family inference or autonomous-selfing proxy is allowed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

AXIS = "reproductive_assurance"
TRAIT = "self_incompatibility"
TARGET = "Justicia betonica"
SOURCE_REF = "Layek et al (2022)"
SOURCE_DOI = "10.1111/1442-1984.12380"
MEYER_MIRROR_DOI = "10.5061/dryad.cc2fqz6hr"
EXPECTED_MEYER_SHA256 = "d6f30a8ad728a3b1b3c8ac5a5143fb935ffc64a38ce06708f6f2df4331737647"
REVIEW_STATUS = "source_methodology_reviewed_reference_backed_strict_direct"


def text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--meyer-xlsx", type=Path, required=True)
    p.add_argument("--coverage", type=Path, required=True)
    p.add_argument("--direct-ledger", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    observed = sha256(args.meyer_xlsx)
    if observed != EXPECTED_MEYER_SHA256:
        raise ValueError(f"pinned Meyer workbook hash mismatch: {observed}")

    coverage = pd.read_csv(args.coverage, dtype=str).fillna("")
    rep = coverage.loc[coverage["axis"].eq(AXIS), ["accepted_species", "quality"]].drop_duplicates()
    if len(rep) != 106295 or rep["accepted_species"].duplicated().any():
        raise ValueError(f"expected fixed 106295-species reproductive coverage, got {len(rep)}")
    quality = dict(zip(rep["accepted_species"], rep["quality"], strict=True))
    if TARGET not in quality:
        raise ValueError(f"target missing from fixed universe: {TARGET}")
    if quality[TARGET]:
        raise ValueError(f"target is no longer a reproductive gap: {TARGET} quality={quality[TARGET]}")

    direct = pd.read_csv(args.direct_ledger, dtype=str).fillna("")
    existing_si = direct.loc[
        direct["accepted_species"].eq(TARGET)
        & direct["trait_name"].eq(TRAIT)
    ]
    if len(existing_si):
        raise ValueError(f"target already has direct self-incompatibility evidence:\n{existing_si.to_string(index=False)}")

    book = pd.ExcelFile(args.meyer_xlsx, engine="openpyxl")
    hits = []
    for sheet in book.sheet_names:
        frame = pd.read_excel(args.meyer_xlsx, sheet_name=sheet, dtype=object, engine="openpyxl").fillna("")
        lookup = {str(c).strip().casefold(): c for c in frame.columns}
        if not {"ref", "genus_species", "mating_system"}.issubset(lookup):
            continue
        subset = frame.loc[
            frame[lookup["ref"]].map(text).eq(SOURCE_REF)
            & frame[lookup["genus_species"]].map(text).str.replace("_", " ", regex=False).eq(TARGET)
        ].copy()
        if len(subset):
            subset.insert(0, "source_sheet", sheet)
            subset.insert(1, "source_row", subset.index + 2)
            subset["mating_system_raw"] = subset[lookup["mating_system"]].map(text)
            hits.append(subset)
    if len(hits) != 1 or len(hits[0]) != 1:
        raise ValueError(f"expected exactly one pinned Layek row for {TARGET}, got {[len(x) for x in hits]}")
    source_row = hits[0].iloc[0]
    if text(source_row["mating_system_raw"]) != "SC":
        raise ValueError(f"pinned Layek state changed for {TARGET}: {source_row['mating_system_raw']!r}")
    hits[0].to_csv(args.output / "justicia_betonica_meyer_source_row.csv", index=False)

    packet = pd.DataFrame([{
        "accepted_species": TARGET,
        "axis": AXIS,
        "source_reference_raw": f"Layek et al. 2022 DOI {SOURCE_DOI}; Meyer mirror DOI {MEYER_MIRROR_DOI}",
        "source_url": f"https://doi.org/{SOURCE_DOI}",
        "source_article_url": f"https://doi.org/{SOURCE_DOI}",
        "source_lineage": f"primary-article:{SOURCE_DOI};pinned-compilation:{MEYER_MIRROR_DOI}",
        "name_match_method": "exact_fixed_universe_binomial_and_exact_primary_article_species",
        "unresolved_reproductive_before_source": True,
        "promotion_allowed": False,
        "genus_rule_training_allowed": False,
        "trait_name": TRAIT,
        "normalized_value": "SC",
        "quality": "high",
        "source_column": "Layek_2022_breeding_system_explicit_self_compatible",
        "source_raw_value": "Original article explicitly states that Justicia betonica is self-compatible while pollinator-dependent for fruit and seed set",
        "review_status": REVIEW_STATUS,
        "acceptance_basis": (
            "The original Plant Species Biology study explicitly investigated the breeding system and reports "
            "Justicia betonica as self-compatible. The pinned Meyer row independently carries SC for the exact "
            "same species. Pollinator dependence is kept separate from compatibility; no autonomous-selfing "
            "proxy or higher-taxon inference is used."
        ),
    }])
    packet.to_csv(
        args.output / "justicia_betonica_unresolved_sc_batch.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )

    summary = {
        "contract": "justicia_betonica_primary_explicit_sc_v1",
        "accepted_species": TARGET,
        "source_doi": SOURCE_DOI,
        "meyer_mirror_doi": MEYER_MIRROR_DOI,
        "meyer_workbook_sha256": observed,
        "pinned_compilation_state": "SC",
        "primary_evidence": "original article explicitly reports Justicia betonica as self-compatible in a breeding-system study",
        "reviewed_strict_packet_rows": 1,
        "reviewed_strict_packet_species": [TARGET],
        "reviewed_strict_packet_state_counts": {"SC": 1},
        "formal_gain": 0,
        "formal_gain_pending_batch_integration": 1,
        "promotion_allowed": False,
        "genus_rule_training_allowed": False,
        "claim_limit": "species-level self-compatibility only; pollinator dependence remains a separate biological property",
    }
    (args.output / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
