"""Build a strict all-species inventory from official WFO-hosted flora archives.

This is an acquisition-stage program.  It resolves provider WFO identifiers to
the fixed 106,295-species universe, records complete source descriptions, and
does not promote extracted traits into the accepted evidence ledger.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import html
import io
import json
import re
import tarfile
import zipfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator

import pandas as pd


EXPECTED_UNIVERSE_SPECIES = 106_295
BACKBONE_MEMBER = "classification.csv"
BACKBONE_BYTES = 121_131_289
BACKBONE_SHA256 = "d0472c3814f3b5bed0af84a4fb7716d36e9200cd585dd1fc9c7b0f244e723ee6"


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def text(value: Any) -> str:
    if value is None:
        return ""
    return str(value).strip()


def normalized_description(value: str) -> str:
    unescaped = html.unescape(text(value))
    without_markup = re.sub(r"<[^>]+>", " ", unescaped)
    return " ".join(without_markup.split())


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def fixed_universe(coverage_path: Path, master_path: Path) -> pd.DataFrame:
    coverage = pd.read_csv(
        coverage_path,
        usecols=["accepted_species"],
        dtype=str,
    ).drop_duplicates()
    if len(coverage) != EXPECTED_UNIVERSE_SPECIES:
        raise ValueError(
            f"fixed universe mismatch: {len(coverage):,} != "
            f"{EXPECTED_UNIVERSE_SPECIES:,}"
        )
    master = pd.read_csv(
        master_path,
        usecols=["accepted_species", "genus", "family", "n_islands", "n_records"],
        dtype=str,
    ).fillna("")
    master = master.loc[
        master["accepted_species"].isin(set(coverage["accepted_species"]))
    ].drop_duplicates("accepted_species")
    result = coverage.merge(master, on="accepted_species", how="left", validate="1:1")
    if result["family"].eq("").any():
        raise ValueError("fixed universe contains species without family metadata")
    return result


def iter_backbone(backbone_zip: Path) -> Iterator[dict[str, str]]:
    if backbone_zip.stat().st_size != BACKBONE_BYTES:
        raise ValueError(
            f"backbone byte mismatch: {backbone_zip.stat().st_size} != {BACKBONE_BYTES}"
        )
    digest = sha256_file(backbone_zip)
    if digest != BACKBONE_SHA256:
        raise ValueError(f"backbone SHA-256 mismatch: {digest}")
    with zipfile.ZipFile(backbone_zip) as archive:
        bad = archive.testzip()
        if bad:
            raise ValueError(f"corrupt backbone ZIP member: {bad}")
        with archive.open(BACKBONE_MEMBER) as raw:
            handle = io.TextIOWrapper(
                raw,
                encoding="utf-8-sig",
                errors="replace",
                newline="",
            )
            yield from csv.DictReader(handle, delimiter="\t")


def build_wfo_identity_map(
    backbone_zip: Path,
    universe: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, dict[str, str]], pd.DataFrame, dict[str, int]]:
    universe_by_name = universe.set_index("accepted_species").to_dict("index")
    matches: list[dict[str, str]] = []
    status_counts: Counter[str] = Counter()
    backbone_rows = 0
    for row in iter_backbone(backbone_zip):
        backbone_rows += 1
        scientific_name = text(row.get("scientificName"))
        if scientific_name not in universe_by_name:
            continue
        status = text(row.get("taxonomicStatus"))
        status_counts[status] += 1
        master = universe_by_name[scientific_name]
        family_match = text(row.get("family")).casefold() == text(master["family"]).casefold()
        rank_match = text(row.get("taxonRank")).casefold() == "species"
        matches.append(
            {
                "accepted_species": scientific_name,
                "master_family": text(master["family"]),
                "master_genus": text(master["genus"]),
                "wfo_id": text(row.get("taxonID")).casefold(),
                "wfo_scientific_name": scientific_name,
                "wfo_family": text(row.get("family")),
                "wfo_rank": text(row.get("taxonRank")),
                "wfo_taxonomic_status": status,
                "wfo_accepted_name_usage_id": text(
                    row.get("acceptedNameUsageID")
                ).casefold(),
                "family_match": str(family_match).lower(),
                "rank_match": str(rank_match).lower(),
            }
        )

    match_table = pd.DataFrame(matches)
    eligible = match_table.loc[
        match_table["family_match"].eq("true")
        & match_table["rank_match"].eq("true")
        & match_table["wfo_taxonomic_status"].isin(["Accepted", "Synonym"])
    ].copy()

    candidates: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    for row in eligible.to_dict("records"):
        status = row["wfo_taxonomic_status"]
        direct_method = (
            "wfo_id_species_rank_accepted_name_exact"
            if status == "Accepted"
            else "wfo_id_species_rank_synonym_name_exact"
        )
        candidates[row["wfo_id"]].append(
            {
                "accepted_species": row["accepted_species"],
                "name_match_method": direct_method,
                "matched_wfo_name": row["wfo_scientific_name"],
                "wfo_name_status": status,
                "master_family": row["master_family"],
            }
        )
        if status == "Synonym" and row["wfo_accepted_name_usage_id"]:
            candidates[row["wfo_accepted_name_usage_id"]].append(
                {
                    "accepted_species": row["accepted_species"],
                    "name_match_method": "wfo_accepted_usage_from_exact_synonym",
                    "matched_wfo_name": row["wfo_scientific_name"],
                    "wfo_name_status": status,
                    "master_family": row["master_family"],
                }
            )

    identity_map: dict[str, dict[str, str]] = {}
    ambiguous_rows: list[dict[str, str]] = []
    for wfo_id, rows in candidates.items():
        species = sorted({row["accepted_species"] for row in rows})
        if len(species) != 1:
            for row in rows:
                ambiguous_rows.append(
                    {
                        "wfo_id": wfo_id,
                        "accepted_species": row["accepted_species"],
                        "name_match_method": row["name_match_method"],
                        "ambiguity": "|".join(species),
                    }
                )
            continue
        methods = sorted(
            rows,
            key=lambda row: (
                0
                if row["name_match_method"]
                == "wfo_id_species_rank_accepted_name_exact"
                else 1
                if row["name_match_method"]
                == "wfo_id_species_rank_synonym_name_exact"
                else 2,
                row["matched_wfo_name"],
            ),
        )
        identity_map[wfo_id] = methods[0]

    ambiguous = pd.DataFrame(ambiguous_rows)
    stats = {
        "backbone_rows": backbone_rows,
        "master_exact_name_rows": len(match_table),
        "master_exact_name_species": int(match_table["accepted_species"].nunique()),
        "eligible_exact_name_rows": len(eligible),
        "eligible_exact_name_species": int(eligible["accepted_species"].nunique()),
        "resolved_wfo_ids": len(identity_map),
        "ambiguous_wfo_ids": int(
            ambiguous["wfo_id"].nunique() if not ambiguous.empty else 0
        ),
        **{f"exact_status_{key}": value for key, value in sorted(status_counts.items())},
    }
    return match_table, identity_map, ambiguous, stats


def dict_rows_from_zip(
    path: Path,
    member: str,
    *,
    header: bool = True,
    fieldnames: list[str] | None = None,
    delimiter: str = "\t",
) -> Iterator[dict[str, str]]:
    with zipfile.ZipFile(path) as archive, archive.open(member) as raw:
        handle = io.TextIOWrapper(
            raw,
            encoding="utf-8-sig",
            errors="replace",
            newline="",
        )
        if header:
            yield from csv.DictReader(handle, delimiter=delimiter)
        else:
            reader = csv.reader(handle, delimiter=delimiter)
            names = fieldnames or []
            for row in reader:
                yield dict(zip(names, row))


def remap_rows(
    rows: Iterator[dict[str, str]],
    *,
    taxon_id: str = "taxonID",
    description: str = "description",
    description_type: str = "type",
    language: str = "language",
    source: str = "source",
    license_name: str = "license",
    rights: str = "rights",
    defaults: dict[str, str] | None = None,
) -> Iterator[dict[str, str]]:
    defaults = defaults or {}
    for row in rows:
        yield {
            "taxonID": text(row.get(taxon_id)),
            "description": text(row.get(description)),
            "type": text(row.get(description_type) or defaults.get("type")),
            "language": text(row.get(language) or defaults.get("language")),
            "source": text(row.get(source)),
            "license": text(row.get(license_name) or defaults.get("license")),
            "rights": text(row.get(rights)),
        }


def fna_rows(path: Path) -> Iterator[dict[str, str]]:
    with path.open(encoding="utf-8-sig", errors="replace", newline="") as handle:
        reader = csv.reader(handle)
        for row in reader:
            if len(row) < 3:
                continue
            yield {
                "taxonID": row[0],
                "type": row[1],
                "language": "en",
                "description": row[2],
                "source": row[3] if len(row) > 3 else "",
                "license": "provider_archive_terms",
            }


def provider_rows(base: Path) -> Iterable[tuple[str, str, Iterator[dict[str, str]]]]:
    yield (
        "brazilian_flora_2020",
        "https://files.worldfloraonline.org/Files/Brazilian%20Flora%202020/"
        "dwca-lista_especies_flora_brasil.zip",
        dict_rows_from_zip(base / "brazilian_flora_2020.zip", "description.txt"),
    )
    yield (
        "flora_north_america",
        "https://files.worldfloraonline.org/Files/FNA/fna_description.txt",
        fna_rows(base / "fna_description.txt"),
    )
    yield (
        "eflora_thailand",
        "https://files.worldfloraonline.org/Files/eFloraOfThailand/FloraOfThailand.zip",
        dict_rows_from_zip(base / "eflora_thailand.zip", "description.txt"),
    )
    yield (
        "flore_du_gabon",
        "https://files.worldfloraonline.org/Files/Flore%20du%20Gabon/FloreduGabon.zip",
        dict_rows_from_zip(base / "flore_gabon.zip", "description.txt"),
    )
    yield (
        "plants_of_nepal",
        "https://files.worldfloraonline.org/Files/PlantsOfNepal/PlantsOfNepal.zip",
        dict_rows_from_zip(base / "plants_nepal.zip", "description.txt"),
    )
    yield (
        "turkey_flora",
        "https://files.worldfloraonline.org/Files/Turkey/Turkey.zip",
        dict_rows_from_zip(
            base / "turkey_flora.zip",
            "description.txt",
            header=False,
            fieldnames=["taxonID", "description"],
        ),
    )


def provider_rows_round2(
    base: Path,
) -> Iterable[tuple[str, str, Iterator[dict[str, str]]]]:
    dwca_fields = [
        "taxonID",
        "description",
        "type",
        "source",
        "language",
        "created",
        "contributor",
        "audience",
        "rights",
        "license",
        "rightsHolder",
        "creator",
    ]
    yield (
        "flora_of_australia",
        "https://files.worldfloraonline.org/Files/Australia/FoA/"
        "WFO_FoA_2020_12_03.zip",
        dict_rows_from_zip(
            base / "australia_foa.zip",
            "WFO_FoA_description_2020_12_03.dsv",
            delimiter="|",
        ),
    )
    yield (
        "flora_of_colombia",
        "https://files.worldfloraonline.org/Files/Colombia/Flora_Of_Colombia.zip",
        remap_rows(
            dict_rows_from_zip(
                base / "flora_colombia.zip",
                "description.txt",
                header=False,
                fieldnames=[
                    "taxonID",
                    "language",
                    "description",
                    "contributor",
                    "extra1",
                    "extra2",
                    "extra3",
                    "extra4",
                    "extra5",
                    "extra6",
                    "extra7",
                    "extra8",
                    "extra9",
                    "extra10",
                    "extra11",
                    "extra12",
                ],
            ),
            defaults={"type": "general"},
        ),
    )
    yield (
        "flora_d_afrique_centrale",
        "https://files.worldfloraonline.org/Files/Flora%20d'Afrique/WFO-fdac.zip",
        dict_rows_from_zip(
            base / "flora_afrique.zip",
            "description_table.csv",
            delimiter=",",
        ),
    )
    yield (
        "flora_helvetica",
        "https://files.worldfloraonline.org/Files/Geneva/"
        "FloraHelvetica-DE-withWFOID.zip",
        dict_rows_from_zip(
            base / "flora_helvetica.zip",
            "description.csv",
            delimiter=",",
        ),
    )
    for provider, filename, archive_url in (
        (
            "flora_mesoamericana",
            "flora_mesoamericana.zip",
            "https://files.worldfloraonline.org/Files/MBG/Flora_Mesoamericana/"
            "Flora_Mesoamericana.zip",
        ),
        (
            "flora_of_nicaragua",
            "flora_nicaragua.zip",
            "https://files.worldfloraonline.org/Files/MBG/Flora_Of_Nicaragua/"
            "Flora_Of_Nicaragua.zip",
        ),
        (
            "flora_of_panama",
            "flora_panama.zip",
            "https://files.worldfloraonline.org/Files/MBG/Flora_Of_Panama/"
            "Flora_Of_Panama.zip",
        ),
        (
            "manual_plants_costa_rica",
            "flora_costa_rica.zip",
            "https://files.worldfloraonline.org/Files/MBG/"
            "Manual_de_Plantas_de_Costa_Rica/"
            "Manual_de_Plantas_de_Costa_Rica.zip",
        ),
    ):
        yield (
            provider,
            archive_url,
            dict_rows_from_zip(
                base / filename,
                "description.txt",
                header=False,
                fieldnames=dwca_fields,
            ),
        )
    for provider, filename, archive_url in (
        (
            "nybg_brittonia",
            "nybg_brittonia.zip",
            "https://files.worldfloraonline.org/Files/NYBG/"
            "nybg_wfo_dwca_brittonia/nybg_wfo_dwca_brittonia.zip",
        ),
        (
            "nybg_floraneotropica",
            "nybg_floraneotropica.zip",
            "https://files.worldfloraonline.org/Files/NYBG/"
            "nybg_wfo_dwca_floraneotropica/nybg_wfo_dwca_floraneotropica.zip",
        ),
        (
            "nybg_memoirs",
            "nybg_memoirs.zip",
            "https://files.worldfloraonline.org/Files/NYBG/"
            "nybg_wfo_dwca_memoirsnybg/nybg_wfo_dwca_memoirsnybg.zip",
        ),
        (
            "nybg_neflora",
            "nybg_neflora.zip",
            "https://files.worldfloraonline.org/Files/NYBG/"
            "nybg_wfo_dwca_neflora/nybg_wfo_dwca_neflora.zip",
        ),
    ):
        yield (
            provider,
            archive_url,
            dict_rows_from_zip(
                base / filename,
                "description.txt",
                delimiter=",",
            ),
        )
    yield (
        "rhododendron_monograph",
        "https://files.worldfloraonline.org/Files/Edinburgh/Rhododendron/"
        "Rhododendron.zip",
        dict_rows_from_zip(
            base / "rhododendron.zip",
            "RhododendronDescriptions.txt",
        ),
    )
    yield (
        "solanaceae_source",
        "https://files.worldfloraonline.org/Files/Solanaceae/WFO/Solanaceae.zip",
        remap_rows(
            dict_rows_from_zip(
                base / "solanaceae.zip",
                "description.txt",
            ),
            taxon_id="wfoID",
            description="descritpion",
        ),
    )


def provider_rows_round3(
    base: Path,
) -> Iterable[tuple[str, str, Iterator[dict[str, str]]]]:
    yield (
        "kew_flora_tropical_east_africa",
        "https://files.worldfloraonline.org/Files/Kew/FTEA/ftea_dwca.zip",
        remap_rows(
            dict_rows_from_zip(
                base / "kew_ftea.zip",
                "description.txt",
            ),
            taxon_id="taxonId",
            defaults={"language": "en"},
        ),
    )
    yield (
        "kew_flora_zambesiaca_current",
        "https://files.worldfloraonline.org/Files/Kew/FZ/fz_dwca.zip",
        remap_rows(
            dict_rows_from_zip(
                base / "kew_fz.zip",
                "description.txt",
                header=False,
                fieldnames=[
                    "taxonId",
                    "type",
                    "source",
                    "description",
                    "identifier",
                    "reference",
                    "coreId",
                ],
            ),
            taxon_id="taxonId",
            defaults={"language": "en"},
        ),
    )
    yield (
        "kew_flora_west_tropical_africa",
        "https://files.worldfloraonline.org/Files/Kew/FWTA/kew_wfo_dwca.zip",
        remap_rows(
            dict_rows_from_zip(
                base / "kew_fwta_wfo.zip",
                "description.txt",
                header=False,
                fieldnames=[
                    "taxonID",
                    "type",
                    "source",
                    "description",
                    "identifier",
                    "reference",
                    "coreId",
                    "scientificNameID",
                    "extra1",
                    "extra2",
                ],
            ),
            defaults={"language": "en"},
        ),
    )


ELIGIBLE_DESCRIPTION_TYPES = {
    "brazilian_flora_2020": {"morphology"},
    "flora_north_america": {"morphology"},
    "eflora_thailand": {"general"},
    "flore_du_gabon": {"morphology"},
    "plants_of_nepal": set(),
    "turkey_flora": {""},
    "flora_of_australia": {"morphology", "general", "diagnostic"},
    "flora_of_colombia": {"general"},
    "flora_d_afrique_centrale": {"morphology"},
    "flora_helvetica": {"morphology"},
    "flora_mesoamericana": {"general"},
    "flora_of_nicaragua": {"general"},
    "flora_of_panama": {"general"},
    "manual_plants_costa_rica": {"general"},
    "nybg_brittonia": {"general"},
    "nybg_floraneotropica": {"general"},
    "nybg_memoirs": {"general"},
    "nybg_neflora": {"general"},
    "rhododendron_monograph": {"morphology"},
    "solanaceae_source": {"morphology", "general"},
    "kew_flora_tropical_east_africa": {"morphology"},
    "kew_flora_west_tropical_africa": {"morphology"},
    "kew_flora_zambesiaca_current": {"morphology"},
}


def acquire_descriptions(
    base: Path,
    identity_map: dict[str, dict[str, str]],
    round2_base: Path | None = None,
    round3_base: Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    records: list[dict[str, str]] = []
    provider_stats: list[dict[str, Any]] = []
    rejected: list[dict[str, str]] = []
    providers = list(provider_rows(base))
    if round2_base is not None:
        providers.extend(provider_rows_round2(round2_base))
    if round3_base is not None:
        providers.extend(provider_rows_round3(round3_base))
    for provider, archive_url, rows in providers:
        types: Counter[str] = Counter()
        languages: Counter[str] = Counter()
        total_rows = 0
        eligible_rows = 0
        matched_rows = 0
        matched_species: set[str] = set()
        matched_ids: set[str] = set()
        seen: set[tuple[str, str, str]] = set()
        for row in rows:
            total_rows += 1
            description_type = text(row.get("type"))
            description_type_key = description_type.casefold()
            language = text(row.get("language")) or (
                "tr" if provider == "turkey_flora" else ""
            )
            types[description_type or "(blank)"] += 1
            languages[language or "(blank)"] += 1
            if description_type_key not in ELIGIBLE_DESCRIPTION_TYPES[provider]:
                continue
            eligible_rows += 1
            wfo_id = text(row.get("taxonID") or row.get("WFOID")).casefold()
            identity = identity_map.get(wfo_id)
            if identity is None:
                continue
            description = normalized_description(text(row.get("description")))
            if not description:
                rejected.append(
                    {
                        "provider": provider,
                        "wfo_id": wfo_id,
                        "reason": "empty_description",
                    }
                )
                continue
            fingerprint = hashlib.sha256(description.encode("utf-8")).hexdigest()
            duplicate_key = (provider, wfo_id, fingerprint)
            if duplicate_key in seen:
                continue
            seen.add(duplicate_key)
            matched_rows += 1
            matched_ids.add(wfo_id)
            matched_species.add(identity["accepted_species"])
            records.append(
                {
                    "accepted_species": identity["accepted_species"],
                    "wfo_id": wfo_id,
                    "matched_wfo_name": identity["matched_wfo_name"],
                    "wfo_name_status": identity["wfo_name_status"],
                    "name_match_method": identity["name_match_method"],
                    "provider": provider,
                    "source_record": text(row.get("source")),
                    "description_type": description_type,
                    "language": language,
                    "description": description,
                    "description_sha256": fingerprint,
                    "source_url": f"https://list.worldfloraonline.org/{wfo_id}",
                    "archive_url": archive_url,
                    "license": text(row.get("license")),
                    "rights": text(row.get("rights")),
                    "retrieved_at_utc": now_utc(),
                    "source_lineage": f"provider_treatment:{provider}:{wfo_id}",
                    "review_status": "source_description_acquired_unreviewed",
                }
            )
        provider_stats.append(
            {
                "provider": provider,
                "archive_url": archive_url,
                "description_rows_total": total_rows,
                "description_rows_trait_eligible": eligible_rows,
                "matched_description_rows": matched_rows,
                "matched_wfo_ids": len(matched_ids),
                "matched_master_species": len(matched_species),
                "description_type_counts": json.dumps(
                    dict(sorted(types.items())),
                    ensure_ascii=False,
                    sort_keys=True,
                ),
                "language_counts": json.dumps(
                    dict(sorted(languages.items())),
                    ensure_ascii=False,
                    sort_keys=True,
                ),
                "api_queries": 0,
                "cost_usd": 0.0,
            }
        )
    return (
        pd.DataFrame(records),
        pd.DataFrame(provider_stats),
        pd.DataFrame(rejected),
    )


def write_file_manifest(output_dir: Path) -> None:
    rows = []
    for path in sorted(output_dir.iterdir()):
        if not path.is_file() or path.name == "file_manifest.sha256":
            continue
        rows.append(f"{sha256_file(path)}  {path.name}")
    (output_dir / "file_manifest.sha256").write_text(
        "\n".join(rows) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--backbone-zip", type=Path, required=True)
    parser.add_argument("--round2-base", type=Path)
    parser.add_argument("--round3-base", type=Path)
    parser.add_argument("--coverage", type=Path, required=True)
    parser.add_argument("--master", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    universe = fixed_universe(args.coverage, args.master)
    matches, identity_map, ambiguous, backbone_stats = build_wfo_identity_map(
        args.backbone_zip,
        universe,
    )
    descriptions, provider_stats, rejected = acquire_descriptions(
        args.base,
        identity_map,
        args.round2_base,
        args.round3_base,
    )

    matches.to_csv(
        args.output / "wfo_master_name_matches.csv.gz",
        index=False,
        compression="gzip",
    )
    ambiguous.to_csv(
        args.output / "ambiguous_wfo_identity_map.csv",
        index=False,
    )
    descriptions.to_csv(
        args.output / "matched_master_species_descriptions.csv.gz",
        index=False,
        compression="gzip",
    )
    provider_stats.to_csv(args.output / "provider_coverage.csv", index=False)
    rejected.to_csv(args.output / "description_rejections.csv", index=False)

    summary = {
        "contract": "wfo_official_bulk_description_inventory_v1",
        "created_at_utc": now_utc(),
        "formal_promotion_status": "not_promoted_acquisition_stage_only",
        "fixed_universe_species": EXPECTED_UNIVERSE_SPECIES,
        "fixed_species_axis_denominator": EXPECTED_UNIVERSE_SPECIES * 3,
        "backbone": {
            "file": args.backbone_zip.name,
            "bytes": args.backbone_zip.stat().st_size,
            "sha256": sha256_file(args.backbone_zip),
            **backbone_stats,
        },
        "descriptions": {
            "rows": len(descriptions),
            "species": int(
                descriptions["accepted_species"].nunique()
                if not descriptions.empty
                else 0
            ),
            "wfo_ids": int(
                descriptions["wfo_id"].nunique() if not descriptions.empty else 0
            ),
            "providers": int(
                descriptions["provider"].nunique() if not descriptions.empty else 0
            ),
            "languages": sorted(
                descriptions["language"].dropna().unique().tolist()
                if not descriptions.empty
                else []
            ),
        },
        "provider_coverage": provider_stats.to_dict("records"),
        "evidence_contract": {
            "accepted_rows": "none_at_this_stage",
            "identity": "species-rank exact accepted or exact synonym, family-consistent",
            "source_lineage": "provider treatment plus WFO taxon identifier",
            "next_gate": "trait extraction, exact quote validation, conflict and cultivar review",
        },
    }
    write_json(args.output / "run_manifest.json", summary)
    write_file_manifest(args.output)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
