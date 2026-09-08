"""Build a strict SI packet for the remaining Grossenbacher-2017 gap.

The current fixed universe spells the taxon ``Hypochoeris radicata``.  The
Grossenbacher et al. (2017) source block retained in the pinned Meyer workbook
contains that exact fixed-universe label with mating_system=SI.  Primary
experimental literature uses the standard orthographic form ``Hypochaeris
radicata`` and directly tests self- versus cross-pollination.

Evidence hierarchy used here:
  * Ortiz et al. 2006 (AJB, DOI 10.3732/ajb.93.2.234): direct compatibility
    experiments; 205 plants across 21 populations, with the species classified
    as self-incompatible despite occasional SC individuals.
  * Picó et al. 2004 (Bot. J. Linn. Soc., DOI
    10.1111/j.1095-8339.2004.00330.x): experimental self and outcross crosses;
    selfing dramatically reduced seed set, supporting self-incompatibility.
  * Grossenbacher et al. 2017 (New Phytologist, DOI 10.1111/nph.14534):
    species-level SI/SC compilation; conflicting species and species with
    among-population breeding-system variation were excluded, while occasional
    SC individuals within otherwise SI species were explicitly permitted.

Only the single exact current reproductive gap is emitted.  No genus/family
inference or genus-rule training is allowed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

AXIS = "reproductive_assurance"
TRAIT = "self_incompatibility"
TARGET = "Hypochoeris radicata"
SOURCE_TAXON = "Hypochaeris radicata"
GROSSENBACHER_DOI = "10.1111/nph.14534"
ORTIZ_DOI = "10.3732/ajb.93.2.234"
PICO_DOI = "10.1111/j.1095-8339.2004.00330.x"
MEYER_DOI = "10.5061/dryad.cc2fqz6hr"
MEYER_REF = "Grossenbacher et al (2017)"
EXPECTED_MEYER_SHA256 = "d6f30a8ad728a3b1b3c8ac5a5143fb935ffc64a38ce06708f6f2df4331737647"
REVIEW_STATUS = "primary_experiment_plus_species_compilation_reviewed_strict_direct"


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
    existing = direct.loc[
        direct["accepted_species"].eq(TARGET)
        & direct["trait_name"].isin(["self_incompatibility", "mating_system"])
    ]
    if len(existing):
        raise ValueError(f"target already has direct reproductive evidence:\n{existing.to_string(index=False)}")

    book = pd.ExcelFile(args.meyer_xlsx, engine="openpyxl")
    hits: list[pd.DataFrame] = []
    for sheet in book.sheet_names:
        frame = pd.read_excel(args.meyer_xlsx, sheet_name=sheet, dtype=object, engine="openpyxl").fillna("")
        lookup = {str(c).strip().casefold(): c for c in frame.columns}
        if not {"ref", "genus_species", "mating_system"}.issubset(lookup):
            continue
        subset = frame.loc[
            frame[lookup["ref"]].map(text).eq(MEYER_REF)
            & frame[lookup["genus_species"]].map(text).str.replace("_", " ", regex=False).eq(TARGET)
        ].copy()
        if len(subset):
            subset.insert(0, "source_sheet", sheet)
            subset.insert(1, "source_row", subset.index + 2)
            subset["mating_system_raw"] = subset[lookup["mating_system"]].map(text)
            hits.append(subset)
    if len(hits) != 1 or len(hits[0]) != 1:
        raise ValueError(f"expected exactly one pinned Grossenbacher row for {TARGET}, got {[len(x) for x in hits]}")
    source_row = hits[0].iloc[0]
    if text(source_row["mating_system_raw"]) != "SI":
        raise ValueError(f"pinned Grossenbacher state changed for {TARGET}: {source_row['mating_system_raw']!r}")

    hits[0].to_csv(args.output / "hypochoeris_grossenbacher_source_row.csv", index=False)

    source_reference = (
        f"Grossenbacher et al. 2017 DOI {GROSSENBACHER_DOI}; "
        f"Ortiz et al. 2006 DOI {ORTIZ_DOI}; Picó et al. 2004 DOI {PICO_DOI}; "
        f"Meyer mirror DOI {MEYER_DOI}"
    )
    packet = pd.DataFrame([
        {
            "accepted_species": TARGET,
            "source_species_name": SOURCE_TAXON,
            "axis": AXIS,
            "source_reference_raw": source_reference,
            "source_url": f"https://doi.org/{GROSSENBACHER_DOI}",
            "source_article_url": f"https://doi.org/{ORTIZ_DOI}",
            "source_lineage": f"primary-doi-set:{ORTIZ_DOI}|{PICO_DOI};compilation:{GROSSENBACHER_DOI}",
            "name_match_method": "fixed_universe_exact_Grossenbacher_label_with_documented_Hypochaeris_orthographic_source_form",
            "unresolved_reproductive_before_source": True,
            "promotion_allowed": False,
            "genus_rule_training_allowed": False,
            "trait_name": TRAIT,
            "normalized_value": "SI",
            "quality": "high",
            "source_column": "Grossenbacher_2017_species_SI_plus_primary_hand_pollination_experiments",
            "source_raw_value": (
                "Grossenbacher species state=SI; Ortiz 2006: H. radicata classified SI after direct "
                "compatibility tests across 21 populations/205 plants, with occasional SC individuals; "
                "Pico 2004: selfing dramatically reduced seed set"
            ),
            "review_status": REVIEW_STATUS,
            "acceptance_basis": (
                "The pinned Grossenbacher 2017 species-level compilation gives the exact fixed-universe "
                "label as SI and explicitly excludes conflicting species/among-population breeding-system "
                "variation. Independent primary controlled-pollination studies directly support species-level "
                "self-incompatibility. Occasional SC individuals do not overturn the source-defined SI species "
                "state under Grossenbacher's stated policy. No proxy, autonomous-selfing shortcut, or taxonomic "
                "higher-level inference is used."
            ),
        }
    ])
    packet.to_csv(
        args.output / "hypochoeris_radicata_unresolved_si_batch.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )

    summary = {
        "contract": "hypochoeris_radicata_primary_experiment_backed_si_v1",
        "accepted_species": TARGET,
        "primary_source_species_name": SOURCE_TAXON,
        "grossenbacher_doi": GROSSENBACHER_DOI,
        "ortiz_doi": ORTIZ_DOI,
        "pico_doi": PICO_DOI,
        "meyer_mirror_doi": MEYER_DOI,
        "meyer_workbook_sha256": observed,
        "pinned_compilation_state": "SI",
        "primary_evidence": {
            "Ortiz_2006": "direct self/cross compatibility experiments; species classified SI; 205 plants across 21 populations; occasional SC individuals retained as exceptions",
            "Pico_2004": "experimental self/outcross crosses; selfing dramatically reduced seed set, indicating SI",
        },
        "reviewed_strict_packet_rows": 1,
        "reviewed_strict_packet_species": [TARGET],
        "reviewed_strict_packet_state_counts": {"SI": 1},
        "formal_gain": 0,
        "formal_gain_pending_batch_integration": 1,
        "promotion_allowed": False,
        "genus_rule_training_allowed": False,
        "claim_limit": "species-level SI only; spelling bridge is confined to the Grossenbacher exact fixed-universe label and the same taxon's primary-literature orthographic form",
    }
    (args.output / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
