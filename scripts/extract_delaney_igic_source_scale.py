"""Build a reviewed current-gap SC/SI packet from Delaney & Igic (2022).

Delaney & Igic compiled a peer-reviewed species-level Fabaceae breeding-system
database.  Meyer et al.'s pinned public Input_Data_MS.xlsx mirror retains that
large-review block and its source `mating_system` labels unchanged.  This
adapter accepts only literal SC/SI states, exact fixed-universe binomials, and
currently unresolved reproductive cells.  DIO, DIC, BSO, D, HET, blanks and any
within-species SC/SI conflict fail closed.

The Delaney compilation is one source lineage here; it cannot train genus rules
or masquerade as independent underlying studies.  Quality is medium because the
public Meyer mirror does not retain per-row underlying reference keys for this
block even though the source paper is a peer-reviewed species-level synthesis.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd

AXIS = "reproductive_assurance"
TRAIT = "self_incompatibility"
REF = "Delaney and Igic"
SOURCE_REPO = "elenacicada/interesting_flowers"
SOURCE_COMMIT = "3d57315fe77097fe582025b884917550139801e8"
SOURCE_FILE = "Input_Data_MS.xlsx"
SOURCE_SHA256 = "d6f30a8ad728a3b1b3c8ac5a5143fb935ffc64a38ce06708f6f2df4331737647"
MEYER_DOI = "10.5061/dryad.cc2fqz6hr"
ARTICLE_DOI = "10.1086/717329"
SOURCE_URL = f"https://doi.org/{ARTICLE_DOI}"
SOURCE_LINEAGE = f"doi:{ARTICLE_DOI}"
REVIEW_STATUS = "source_methodology_reviewed_reference_backed_strict_direct"
BINOMIAL = re.compile(r"^[A-Z][A-Za-z.-]+ [a-z][A-Za-z.-]+$")
EXPECTED_PACKET = {
    "Hedysarum spinosissimum": "SC",
    "Pueraria phaseoloides": "SC",
    "Scorpiurus subvillosus": "SC",
}
EXPECTED_CONFLICT = {"Bauhinia pauletia"}


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


def locate_block(path: Path) -> tuple[str, pd.DataFrame]:
    book = pd.ExcelFile(path, engine="openpyxl")
    hits: list[tuple[str, pd.DataFrame]] = []
    for sheet in book.sheet_names:
        frame = pd.read_excel(path, sheet_name=sheet, dtype=object, engine="openpyxl").fillna("")
        lookup = {str(c).strip().casefold(): c for c in frame.columns}
        if "ref" not in lookup or "genus_species" not in lookup or "mating_system" not in lookup:
            continue
        subset = frame.loc[frame[lookup["ref"]].map(text).eq(REF)].copy()
        if len(subset):
            hits.append((sheet, subset))
    if len(hits) != 1:
        raise ValueError(f"expected one {REF!r} block, found {[x[0] for x in hits]}")
    return hits[0]


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input-xlsx", type=Path, required=True)
    p.add_argument("--coverage", type=Path, required=True)
    p.add_argument("--direct-ledger", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    observed = sha256(args.input_xlsx)
    if observed != SOURCE_SHA256:
        raise ValueError(f"pinned Meyer workbook hash mismatch: {observed}")

    sheet, source = locate_block(args.input_xlsx)
    lookup = {str(c).strip().casefold(): c for c in source.columns}
    species_col = lookup["genus_species"]
    state_col = lookup["mating_system"]

    source = source.copy()
    source["source_sheet"] = sheet
    source["source_row_in_sheet"] = source.index + 2
    source["accepted_species_candidate"] = source[species_col].map(text).str.replace("_", " ", regex=False)
    source["raw_mating_system"] = source[state_col].map(text)
    source["literal_strict_state"] = source["raw_mating_system"].where(
        source["raw_mating_system"].isin({"SC", "SI"}), ""
    )
    source["exact_binomial"] = source["accepted_species_candidate"].map(lambda x: bool(BINOMIAL.fullmatch(x)))

    coverage = pd.read_csv(args.coverage, dtype=str).fillna("")
    rep = coverage.loc[coverage["axis"].eq(AXIS), ["accepted_species", "quality"]].drop_duplicates()
    if len(rep) != 106295 or rep["accepted_species"].duplicated().any():
        raise ValueError(f"expected fixed 106295-species reproductive coverage, got {len(rep)}")
    quality = dict(zip(rep["accepted_species"], rep["quality"], strict=True))
    universe = set(quality)
    source["in_fixed_universe_exact"] = source["accepted_species_candidate"].isin(universe)
    source["current_reproductive_quality"] = source["accepted_species_candidate"].map(quality).fillna("")
    source["currently_unresolved_reproductive"] = (
        source["in_fixed_universe_exact"] & source["current_reproductive_quality"].eq("")
    )

    direct = pd.read_csv(args.direct_ledger, dtype=str).fillna("")
    current_direct_pairs = set(zip(direct["accepted_species"], direct["trait_name"], strict=False))

    strict_rows = source.loc[
        source["literal_strict_state"].ne("")
        & source["exact_binomial"]
        & source["in_fixed_universe_exact"]
        & source["currently_unresolved_reproductive"]
    ].copy()

    value_counts = strict_rows.groupby("accepted_species_candidate")["literal_strict_state"].nunique()
    conflict_species = set(value_counts.loc[value_counts.gt(1)].index)
    clean = strict_rows.loc[~strict_rows["accepted_species_candidate"].isin(conflict_species)].copy()
    clean = clean.sort_values("source_row_in_sheet").drop_duplicates(
        ["accepted_species_candidate", "literal_strict_state"], keep="first"
    )

    source_reference = (
        f"Delaney & Igic (2022), International Journal of Plant Sciences, DOI {ARTICLE_DOI}; "
        f"Meyer et al. public mirror DOI {MEYER_DOI}, GitHub {SOURCE_REPO}@{SOURCE_COMMIT}"
    )
    packet_rows: list[dict[str, object]] = []
    for row in clean.to_dict("records"):
        species = text(row["accepted_species_candidate"])
        state = text(row["literal_strict_state"])
        if (species, TRAIT) in current_direct_pairs:
            raise ValueError(f"current direct ledger already contains {species} {TRAIT}")
        packet_rows.append({
            "accepted_species": species,
            "source_species_name": species,
            "axis": AXIS,
            "source_reference_raw": source_reference,
            "source_url": SOURCE_URL,
            "source_article_url": SOURCE_URL,
            "source_lineage": SOURCE_LINEAGE,
            "name_match_method": "exact_accepted_species_in_fixed_universe",
            "unresolved_reproductive_before_source": True,
            "promotion_allowed": False,
            "genus_rule_training_allowed": False,
            "trait_name": TRAIT,
            "normalized_value": state,
            "quality": "medium",
            "source_column": "Delaney_and_Igic_species_level_mating_system",
            "source_raw_value": (
                f"sheet={sheet}; row={int(row['source_row_in_sheet'])}; "
                f"genus_species={text(row[species_col])}; mating_system={text(row[state_col])}"
            ),
            "review_status": REVIEW_STATUS,
            "acceptance_basis": (
                "Delaney & Igic is a peer-reviewed species-level Fabaceae breeding-system synthesis. "
                "The pinned Meyer input retains the Delaney block's literal mating_system labels; "
                "only unambiguous SC/SI rows are accepted. DIO/DIC/BSO/D/HET/blank states and "
                "within-species SI/SC conflicts are excluded. The compilation is one lineage and "
                "cannot train genus rules."
            ),
        })

    packet = pd.DataFrame(packet_rows)
    observed_packet = dict(zip(packet["accepted_species"], packet["normalized_value"], strict=True)) if len(packet) else {}
    if observed_packet != EXPECTED_PACKET:
        raise ValueError(f"Delaney strict current-gap packet changed: {observed_packet}")
    if conflict_species != EXPECTED_CONFLICT:
        raise ValueError(f"Delaney current-gap conflict set changed: {sorted(conflict_species)}")

    source["selection_reason"] = ""
    source.loc[~source["exact_binomial"], "selection_reason"] = "not_exact_binomial"
    source.loc[source["exact_binomial"] & ~source["in_fixed_universe_exact"], "selection_reason"] = "not_fixed_universe_exact"
    source.loc[source["in_fixed_universe_exact"] & ~source["currently_unresolved_reproductive"], "selection_reason"] = "already_resolved_reproductive_axis"
    source.loc[
        source["currently_unresolved_reproductive"] & source["literal_strict_state"].eq(""),
        "selection_reason",
    ] = "non_strict_or_missing_source_state"
    source.loc[source["accepted_species_candidate"].isin(conflict_species), "selection_reason"] = "within_source_sc_si_conflict"
    source.loc[
        source["accepted_species_candidate"].isin(EXPECTED_PACKET), "selection_reason"
    ] = "selected_reviewed_strict_sc_si"

    packet_columns = [
        "accepted_species", "source_species_name", "axis", "source_reference_raw",
        "source_url", "source_article_url", "source_lineage", "name_match_method",
        "unresolved_reproductive_before_source", "promotion_allowed",
        "genus_rule_training_allowed", "trait_name", "normalized_value", "quality",
        "source_column", "source_raw_value", "review_status", "acceptance_basis",
    ]
    packet = packet[packet_columns].sort_values("accepted_species").reset_index(drop=True)
    source.to_csv(
        args.output / "delaney_igic_full_source_audit.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    packet.to_csv(
        args.output / "delaney_igic_unresolved_si_sc_batch.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    strict_rows.to_csv(args.output / "delaney_igic_current_gap_strict_rows.csv", index=False)

    summary = {
        "contract": "delaney_igic_2022_source_scale_strict_sc_si_v1",
        "source_article_doi": ARTICLE_DOI,
        "meyer_mirror_doi": MEYER_DOI,
        "source_repo": SOURCE_REPO,
        "source_commit": SOURCE_COMMIT,
        "source_file": SOURCE_FILE,
        "source_sha256": observed,
        "source_sheet": sheet,
        "source_rows": int(len(source)),
        "raw_unique_taxon_labels": int(source["accepted_species_candidate"].nunique()),
        "exact_fixed_universe_species": int(source.loc[source["in_fixed_universe_exact"], "accepted_species_candidate"].nunique()),
        "exact_current_unresolved_reproductive_species": int(source.loc[source["currently_unresolved_reproductive"], "accepted_species_candidate"].nunique()),
        "current_gap_literal_sc_si_rows": int(len(strict_rows)),
        "current_gap_literal_sc_si_species": int(strict_rows["accepted_species_candidate"].nunique()),
        "within_source_current_gap_conflicts": sorted(conflict_species),
        "reviewed_strict_packet_rows": int(len(packet)),
        "reviewed_strict_packet_species": packet["accepted_species"].tolist(),
        "reviewed_strict_packet_state_counts": packet["normalized_value"].value_counts().sort_index().astype(int).to_dict(),
        "formal_gain": 0,
        "formal_gain_pending_batch_integration": int(len(packet)),
        "promotion_allowed": False,
        "genus_rule_training_allowed": False,
        "claim_limit": "Peer-reviewed species-level compilation; literal SC/SI only; no DIO/BSO/other proxy mapping and no underlying-study independence claim.",
    }
    (args.output / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
