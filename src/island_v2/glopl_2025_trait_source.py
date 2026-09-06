"""Build direct trait evidence from the 2025 peer-reviewed GloPL table.

The source workbook reports species-level autofertility, self compatibility,
floral symmetry, and several related but distinct traits.  Only autofertility,
self compatibility, and symmetry enter the frozen strict three-axis ledger.
Floral specialisation, reward, sexual system, heterostyly, dichogamy, and
herkogamy are preserved in a separate ledger and cannot silently resolve a
strict axis.

The workbook is a redistribution of the GloPL/Burns trait compilation.  Its
lineage therefore identifies that compilation, not the Nature provider URL,
so the same underlying observation cannot be counted as independent evidence
when it appears through another provider.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
import typer

from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS

app = typer.Typer(add_completion=False, no_args_is_help=True)

ARTICLE_DOI = "10.1038/s41467-025-61032-5"
ARTICLE_URL = f"https://doi.org/{ARTICLE_DOI}"
SUPPLEMENT_URL = (
    "https://www.ebi.ac.uk/europepmc/webservices/rest/"
    "PMC12216600/supplementaryFiles"
)
SOURCE_ARTIFACT = "41467_2025_61032_MOESM7_ESM.xlsx"
SOURCE_SHA256 = "e6fe5c2d798269c8b41a678d029ccc25dc600e1f0626ed834cf45c494b8a9775"
RETRIEVAL_DATE = "2026-08-27"
SOURCE_GROUP = "glopl_2025_trait_compilation"
SOURCE_PROVIDER = "Lin et al. 2025 GloPL Supplementary Data"
SOURCE_CITATION = (
    "Lin et al. (2025), Global meta-analysis shows that threatened flowering "
    "plants have higher pollination deficits, Nature Communications 16:5882, "
    f"DOI {ARTICLE_DOI}, Source Data Fig. 2."
)
SOURCE_RUN_ID = "europe-pmc:PMC12216600:supplementary-files"

STRICT_MAP = {
    "AF": {
        "trait_name": "autonomous_selfing_capacity",
        "axis": "reproductive_assurance",
        "values": {"yes": "autonomous", "no": "absent"},
    },
    "BS1": {
        "trait_name": "self_incompatibility",
        "axis": "reproductive_assurance",
        "values": {"sc": "SC", "si": "SI", "p": "mixed_or_variable"},
    },
    "flower_symmetry_1z_0a": {
        "trait_name": "floral_symmetry",
        "axis": "floral_structural_complexity",
        "values": {"0": "actinomorphic", "1": "zygomorphic"},
    },
}

INDEPENDENT_MAP = {
    "floral_structure": {
        "trait_name": "floral_specialisation_class",
        "values": {"g": "generalized", "s": "specialized"},
    },
    "type_of_reward": {
        "trait_name": "reward_type",
        "values": {"p": "pollen", "n": "nectar_and_pollen", "s": "specialized"},
    },
    "accessibility_1hard_0easy": {
        "trait_name": "reward_accessibility",
        "values": {"0": "easy", "1": "hard"},
    },
    "BS2": {
        "trait_name": "heterostyly",
        "values": {"yes": "present", "no": "absent"},
    },
    "BS3": {
        "trait_name": "sex_system",
        "values": {"hm": "bisexual_hermaphroditic", "none": "not_bisexual"},
    },
    "BS4": {
        "trait_name": "dichogamy",
        "values": {"yes": "present", "none": "absent"},
    },
    "BS5": {
        "trait_name": "herkogamy",
        "values": {"yes": "present", "no": "absent"},
    },
}

FAMILY_ALIASES = {
    "Compositae": "Asteraceae",
    "Leguminosae": "Fabaceae",
    "Papilionaceae": "Fabaceae",
    "Gramineae": "Poaceae",
    "Labiatae": "Lamiaceae",
    "Cruciferae": "Brassicaceae",
    "Umbelliferae": "Apiaceae",
}

AUDIT_COLUMNS = [
    "source_row",
    "submitted_species",
    "source_family",
    "master_family",
    "trait_name",
    "raw_value",
    "normalized_value",
    "decision",
    "reason",
]


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _record_id(*parts: object) -> str:
    payload = "\n".join(_text(part) for part in parts)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:24]


def _family(value: object) -> str:
    value = _text(value)
    return FAMILY_ALIASES.get(value, value)


def _source_excerpt(rows: pd.DataFrame, column: str) -> str:
    excerpts: list[str] = []
    for row in rows.sort_values("source_row").itertuples(index=False):
        record = row._asdict()
        excerpts.append(
            f"Source Data Fig.2 row {record['source_row']}: "
            f"Species_accepted_names={record['Species_accepted_names']}; "
            f"Family={record['Family']}; {column}={record[column]}; "
            f"study={record['Author']} ({record['Year']}); DOI={record['DOI']}"
        )
    return " || ".join(excerpts)


def _write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        frame.to_csv(
            path,
            index=False,
            compression={"method": "gzip", "mtime": 0},
        )
    else:
        frame.to_csv(path, index=False)


def _prepare_source(source: pd.DataFrame) -> pd.DataFrame:
    required = {
        "Author",
        "Year",
        "DOI",
        "Species_accepted_names",
        "Family",
        *STRICT_MAP,
        *INDEPENDENT_MAP,
    }
    missing = required.difference(source.columns)
    if missing:
        raise ValueError(f"GloPL Source Data missing columns: {sorted(missing)}")
    result = source.copy().fillna("")
    result.insert(0, "source_row", range(2, len(result) + 2))
    result["submitted_species"] = result["Species_accepted_names"].map(
        lambda value: _text(value).replace("_", " ")
    )
    return result


def build_package(
    source: pd.DataFrame,
    master: pd.DataFrame,
    universe_species: set[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return strict evidence, independent-trait evidence, and full audit."""

    source = _prepare_source(source)
    if missing := {"accepted_species", "family"}.difference(master.columns):
        raise ValueError(f"master list missing columns: {sorted(missing)}")
    master = master[["accepted_species", "family"]].fillna("").drop_duplicates(
        "accepted_species"
    )
    if universe_species is not None:
        master = master.loc[master["accepted_species"].isin(universe_species)].copy()
    master_family = {
        _text(row.accepted_species): _family(row.family)
        for row in master.itertuples(index=False)
    }

    audit_rows: list[dict[str, str | int]] = []
    accepted_rows: list[dict[str, Any]] = []
    mappings = [(True, column, spec) for column, spec in STRICT_MAP.items()]
    mappings += [(False, column, spec) for column, spec in INDEPENDENT_MAP.items()]

    for record in source.to_dict("records"):
        submitted = _text(record["submitted_species"])
        source_family = _family(record["Family"])
        species_rank = bool(
            re.fullmatch(r"[A-Z][A-Za-z.-]+ [a-z][A-Za-z.-]+", submitted)
        )
        for strict, column, spec in mappings:
            raw = _text(record[column])
            if not raw:
                continue
            normalized = spec["values"].get(raw, "")
            reason = "selected"
            if not species_rank:
                reason = "not_exact_species_rank"
            elif submitted not in master_family:
                reason = "not_exact_island_accepted_species"
            elif source_family != master_family[submitted]:
                reason = "family_conflict"
            elif not normalized:
                reason = "unmapped_source_value"
            audit_rows.append(
                {
                    "source_row": record["source_row"],
                    "submitted_species": submitted,
                    "source_family": source_family,
                    "master_family": master_family.get(submitted, ""),
                    "trait_name": spec["trait_name"],
                    "raw_value": raw,
                    "normalized_value": normalized,
                    "decision": "accept" if reason == "selected" else "exclude",
                    "reason": reason,
                }
            )
            if reason == "selected":
                accepted_rows.append(
                    {
                        **record,
                        "accepted_species": submitted,
                        "strict_three_axis_included": strict,
                        "source_column": column,
                        "trait_name": spec["trait_name"],
                        "axis": spec.get("axis", ""),
                        "normalized_value": normalized,
                    }
                )

    accepted = pd.DataFrame(accepted_rows).fillna("")
    audit = pd.DataFrame(audit_rows, columns=AUDIT_COLUMNS).fillna("")
    evidence_rows: list[dict[str, Any]] = []
    independent_rows: list[dict[str, Any]] = []
    keys = [
        "accepted_species",
        "strict_three_axis_included",
        "source_column",
        "trait_name",
        "axis",
    ]
    for key, rows in accepted.groupby(keys, sort=True, dropna=False):
        species, strict, column, trait, axis = key
        values = sorted(set(rows["normalized_value"]))
        if len(values) != 1:
            mask = (
                audit["submitted_species"].eq(species)
                & audit["trait_name"].eq(trait)
                & audit["decision"].eq("accept")
            )
            audit.loc[mask, "decision"] = "exclude"
            audit.loc[mask, "reason"] = "within_source_direct_conflict"
            continue
        value = values[0]
        lineage = f"glopl-burns-trait-compilation:{species.casefold()}:{trait}"
        common = {
            "accepted_species": species,
            "axis": axis,
            "trait_name": trait,
            "normalized_value": value,
            "quality": "high",
            "source_group": SOURCE_GROUP,
            "source_provider": SOURCE_PROVIDER,
            "source_url": ARTICLE_URL,
            "source_record_id": "glopl2025:" + _record_id(species, trait, value),
            "source_citation": SOURCE_CITATION,
            "source_excerpt": _source_excerpt(rows, column),
            "evidence_scope": "species_direct",
            "name_match_method": "accepted_name_exact",
            "source_lineage": lineage,
            "lineage_method": "underlying_glopl_burns_trait_compilation",
            "source_run_id": SOURCE_RUN_ID,
            "source_artifact": SOURCE_ARTIFACT,
            "source_file": SOURCE_ARTIFACT,
            "acceptance_contract": (
                "peer_reviewed_species_trait_table_exact_identity_v1"
            ),
        }
        if strict:
            evidence_rows.append(common)
        else:
            independent_rows.append(
                {
                    **common,
                    "strict_three_axis_included": False,
                    "strict_exclusion_reason": "independent_trait_not_in_frozen_axis",
                }
            )

    strict_evidence = pd.DataFrame(evidence_rows, columns=EVIDENCE_COLUMNS)
    independent_columns = [
        *EVIDENCE_COLUMNS,
        "strict_three_axis_included",
        "strict_exclusion_reason",
    ]
    independent = pd.DataFrame(independent_rows, columns=independent_columns)
    return strict_evidence, independent, audit


@app.command("build")
def build(
    source_xlsx: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    master_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    universe_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    output_dir: Annotated[Path, typer.Option(file_okay=False)],
) -> None:
    """Build immutable strict and independent trait ledgers."""

    source_hash = _sha256(source_xlsx)
    if source_hash != SOURCE_SHA256:
        raise ValueError(
            f"source hash mismatch: expected {SOURCE_SHA256}, observed {source_hash}"
        )
    source = pd.read_excel(source_xlsx, sheet_name="Fig.2", dtype=str).fillna("")
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    universe = pd.read_csv(
        universe_csv, usecols=["accepted_species"], dtype=str
    ).fillna("")
    universe_species = set(universe["accepted_species"])
    if len(universe_species) != 106_295:
        raise ValueError(
            "strict universe must contain exactly 106295 accepted species; "
            f"observed {len(universe_species)}"
        )
    strict, independent, audit = build_package(
        source, master, universe_species=universe_species
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "glopl_2025_reviewed_direct_evidence.csv.gz": strict,
        "glopl_2025_independent_trait_evidence.csv.gz": independent,
        "glopl_2025_extraction_audit.csv.gz": audit,
    }
    hashes: dict[str, str] = {}
    for filename, frame in outputs.items():
        path = output_dir / filename
        _write_csv(frame, path)
        hashes[filename] = _sha256(path)

    summary = {
        "contract": "glopl_2025_trait_source_package_v1",
        "created_at": RETRIEVAL_DATE,
        "source": {
            "article_doi": ARTICLE_DOI,
            "supplement_url": SUPPLEMENT_URL,
            "artifact": SOURCE_ARTIFACT,
            "sha256": source_hash,
        },
        "input": {
            "source_rows": len(source),
            "source_species": int(
                source["Species_accepted_names"].replace("", pd.NA).nunique()
            ),
            "master_species": len(universe_species),
        },
        "strict_direct": {
            "rows": len(strict),
            "species": int(strict["accepted_species"].nunique()),
            "species_trait": int(
                strict[["accepted_species", "trait_name"]]
                .drop_duplicates()
                .shape[0]
            ),
            "by_trait": strict["trait_name"].value_counts().sort_index().to_dict(),
        },
        "independent_traits": {
            "rows": len(independent),
            "species_trait": int(
                independent[["accepted_species", "trait_name"]]
                .drop_duplicates()
                .shape[0]
            ),
            "by_trait": independent["trait_name"]
            .value_counts()
            .sort_index()
            .to_dict(),
            "strict_three_axis_included": False,
        },
        "exclusions": audit.loc[audit["decision"].eq("exclude"), "reason"]
        .value_counts()
        .sort_index()
        .to_dict(),
        "checks": {
            "accepted_names_exact_only": True,
            "family_conflicts_excluded": True,
            "within_source_conflicts_excluded": True,
            "reward_not_in_structure_axis": True,
            "floral_specialisation_not_substituted_for_floral_form": True,
            "family_inference_absent": True,
            "global_fallback_absent": True,
            "source_lineage_underlying_not_provider_url": True,
        },
        "artifact_hashes": hashes,
    }
    summary_path = output_dir / "glopl_2025_trait_source_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    app()
