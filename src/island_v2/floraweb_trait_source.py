"""Build a provenance-complete FloraWeb/BiolFlor trait source package.

The input is the frozen 2026-06-24 FloraWeb scrape published by taxifydb.
Controlled BiolFlor reproductive and flower-type fields are mapped to the
repository ontology.  Free-text Rothmaler descriptions contribute flower
colour only when a colour is explicitly attached to a flower, corolla, petal,
or perianth phrase.  Pollination vector, reward, dichogamy, and sex system are
preserved as independent traits and are never projected onto the strict three
axes.
"""

from __future__ import annotations

import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

RAW_SHA256 = "a32b5779b2d7a92284fdfdf0f77ec5a7364283bb9a8f277009d6e29116e12382"
PROVENANCE_SHA256 = "b6a074c8f9c5e62b053e44c2afdcfb0aff417625820bfb34195c82456943afee"
SNAPSHOT_RELEASE = "scrape-snapshots-2026.06"
SNAPSHOT_URL = (
    "https://github.com/gcol33/taxifydb/releases/download/"
    f"{SNAPSHOT_RELEASE}/floraweb_raw_2026-06-24.csv"
)
RETRIEVAL_DATE = "2026-06-24"
SOURCE_CITATION = (
    "FloraWeb. Daten und Informationen zu Wildpflanzen und zur Vegetation "
    "Deutschlands. Bundesamt fuer Naturschutz, Bonn. https://www.floraweb.de/ "
    "(accessed 2026-06-24); reproductive fields derive from Klotz, Kuehn & "
    "Durka (eds) 2002, BIOLFLOR, Schriftenreihe fuer Vegetationskunde 38."
)

STRICT_EVIDENCE_COLUMNS = [
    "accepted_species",
    "trait_name",
    "normalized_value",
    "quality",
    "source_group",
    "source_provider",
    "source_url",
    "source_record_id",
    "source_citation",
    "source_excerpt",
    "evidence_scope",
    "name_match_method",
    "source_lineage",
    "lineage_method",
    "source_run_id",
    "source_artifact",
    "acceptance_contract",
    "inference_rule",
]

INDEPENDENT_EVIDENCE_COLUMNS = [
    *STRICT_EVIDENCE_COLUMNS,
    "strict_three_axis_included",
]

MATING_VALUES = {
    "xenogam": "predominantly_outcrossing",
    "fakultativ xenogam": "predominantly_outcrossing",
    "gemischte Befruchtung": "mixed_mating",
    "fakultativ autogam": "predominantly_selfing",
    "obligat Autogam": "obligate_selfing",
}

SI_VALUES = {
    "selbstkompatibel": "SC",
    "selbstinkompatibel": "SI",
    "± selbstinkompatibel": "mixed_or_variable",
    "± selbstkompatibel": "mixed_or_variable",
}

AUTONOMY_VALUES = {
    "immer Selbstbestäubung": "autonomous",
    "in der Regel Selbstbestäubung": "autonomous",
    "häufig Selbstbestäubung": "autonomous",
    "bei ausbleib. Fremdbest. Selbstbestäubung": "facilitated",
    "selten Selbstbestäubung": "mixed_or_variable",
    "möglich Selbstbestäubung": "mixed_or_variable",
    "nie Selbstbestäubung": "absent",
}

CLEISTOGAMY_VALUES = {
    "in der Regel Kleistogamie": "facultative",
    "häufig Kleistogamie": "facultative",
    "selten Kleistogamie": "facultative",
    "nie Kleistogamie": "absent",
}

DICHOGAMY_VALUES = {
    "homogam": "absent",
    "proterandrisch": "protandry",
    "ausgeprägt protandrisch": "protandry",
    "leicht proterandrisch": "protandry",
    "protogyn": "protogyny",
    "ausgeprägt protogyn": "protogyny",
    "leicht protogyn": "protogyny",
}

FORM_PATTERNS: tuple[tuple[str, str, str], ...] = (
    ("Scheibenblumen", "floral_form", "open_radial"),
    ("Köpfchenblumen", "floral_form", "composite_head"),
    ("Stieltellerblumen", "floral_form", "salverform"),
    ("Glockenblumen", "floral_form", "bell_campanulate"),
    ("Trichterblumen", "floral_form", "funnel_trumpet"),
    ("Schmetterlingsblumen", "floral_form", "papilionaceous"),
    ("Pinsel- und Bürstenblumen", "floral_form", "brush_puff"),
    ("Eigentliche Lippenblumen", "floral_form", "bilabiate"),
    ("Lippenblumen, Rachenblumen", "floral_form", "bilabiate"),
    ("Lippenblumen, Maskenblumen", "floral_form", "bilabiate"),
    ("Kolbenblumen", "inflorescence_display", "composite_display"),
)

COLOUR_PATTERNS: tuple[tuple[re.Pattern[str], str], ...] = (
    (re.compile(r"\b(?:rein)?wei(?:ß|ss)(?:lich|l\.)?\b|\bcreme(?:farben)?\b", re.IGNORECASE), "white"),
    (
        re.compile(
            r"\b(?:hell|dunkel|dkl)?grün(?:lich|l\.)?\b|"
            r"\b(?:schwarz|rot)?braun(?:violett|lich)?\b|"
            r"\bbräunlich\b|\bunscheinbar\b",
            re.IGNORECASE,
        ),
        "green_brown_inconspicuous",
    ),
    (
        re.compile(
            r"\b(?:hell|blass|gold|zitronen)?gelb(?:lich|l\.|e[snmr]?)?\b|"
            r"\b(?:rot)?orange(?:farben)?\b",
            re.IGNORECASE,
        ),
        "yellow_orange",
    ),
    (
        re.compile(
            r"\b(?:scharlach|blut|purpur|karmin|rosen|dunkel|dkl|hell)?rot(?:lich|l\.)?\b|"
            r"\brötlich\b|\b(?:zart|blass)?rosa\b|\bpink\b|\bpurpur(?:n|farben|rosa)?\b",
            re.IGNORECASE,
        ),
        "red_pink",
    ),
    (
        re.compile(
            r"\b(?:hell|himmel|azur|dunkel|dkl)?blau(?:violett|lich|l\.)?\b|"
            r"\b(?:hell|dunkel|dkl)?violett\b|\blila(?:blau)?\b",
            re.IGNORECASE,
        ),
        "blue_purple",
    ),
)

FLOWER_SUBJECT_RE = re.compile(
    r"\b(?:Blüten|Blüte)(?!stand|köpf|kopf|hüll|zeit)\b|"
    r"\b(?:Blütenkrone|Krone|Kronblätter|Perigon(?:blätter)?|Blütenhülle|"
    r"Zungenblüten|Röhrenblüten|Randblüten)\b",
    re.IGNORECASE,
)
COLOUR_CLAUSE_STOP_RE = re.compile(
    r"\b(?:Staub(?:blatt|blätter|beutel|fäden?)|Antheren?|Narbe|Griffel|"
    r"Frucht(?:knoten|stand)?|Samen|Laubblätter|Stgblätter|Grundblätter|"
    r"Hochblätter)\b",
    re.IGNORECASE,
)
FADED_RE = re.compile(r"\b(?:welke|verblühte|abgeblüht|nach der Blüte)\b", re.IGNORECASE)


def _text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _stable_id(*parts: object) -> str:
    payload = "\x1f".join(_text(part) for part in parts)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:24]


def _write_gzip(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, compression={"method": "gzip", "mtime": 0})


def _floraweb_url(page: str, name_usage_id: str) -> str:
    return f"https://www.floraweb.de/php/{page}.php?name-use-id={name_usage_id}"


def _lineage(page: str, label: str, name_usage_id: str) -> tuple[str, str]:
    if page == "morphologie" and label == "Beschreibung":
        return (
            f"rothmaler:floraweb:name-use-id:{name_usage_id}:description",
            "underlying_source_trait",
        )
    slug = re.sub(r"[^a-z0-9]+", "_", label.casefold()).strip("_")
    return (
        f"biolflor:floraweb:name-use-id:{name_usage_id}:trait:{slug}",
        "underlying_source_trait",
    )


def _colour_states(description: str) -> tuple[list[str], str]:
    """Extract flower-bound colours and return states plus the exact quote."""
    states: set[str] = set()
    quotes: list[str] = []
    for subject in FLOWER_SUBJECT_RE.finditer(description):
        start = subject.start()
        prefix = description[max(0, start - 18) : subject.end()]
        if FADED_RE.search(prefix):
            continue
        end = min(len(description), subject.end() + 115)
        segment = description[start:end]
        period = segment.find(".", 1)
        if period >= 0:
            segment = segment[: period + 1]
        stop = COLOUR_CLAUSE_STOP_RE.search(segment, pos=max(0, subject.end() - start))
        if stop:
            segment = segment[: stop.start()]
        if FADED_RE.search(segment):
            continue
        found = {
            state
            for pattern, state in COLOUR_PATTERNS
            if pattern.search(segment)
        }
        if found:
            states.update(found)
            quotes.append(_text(segment))

    # Handle adjective + floral organ constructions such as "gelbe Kronblätter".
    for pattern, state in COLOUR_PATTERNS:
        for match in pattern.finditer(description):
            tail = description[match.start() : min(len(description), match.end() + 28)]
            context = description[max(0, match.start() - 18) : match.end() + 28]
            if FADED_RE.search(context):
                continue
            if re.search(
                r"\b(?:Kronblätter|Perigon(?:blätter)?|Blüten)(?!stand|köpf)\b",
                tail,
                re.IGNORECASE,
            ):
                states.add(state)
                quotes.append(_text(tail))

    return sorted(states), " | ".join(dict.fromkeys(quotes))


def _independent_states(label: str, values: list[str]) -> dict[str, set[str]]:
    result: dict[str, set[str]] = defaultdict(set)
    if label == "Bestäubung (Pollenvektoren)":
        for value in values:
            folded = value.casefold()
            if "windbestäubung" in folded:
                result["pollen_vector_mode"].add("abiotic_wind")
            if "insektenbestäubung" in folded:
                result["pollen_vector_mode"].add("biotic")
    elif label == "Belohnung für Bestäuber":
        for value in values:
            folded = value.casefold()
            if folded.startswith(("vorhanden nektar", "reichlich nektar", "wenig nektar")):
                result["reward_type"].add("nectar")
            elif folded.startswith(("vorhanden pollen", "reichlich pollen", "wenig pollen")):
                result["reward_type"].add("pollen")
            elif folded.startswith("vorhanden öl"):
                result["reward_type"].add("oil")
            elif folded.startswith("vorhanden täuschung"):
                result["reward_type"].add("rewardless")
    elif label == "Dichogamie (zeitliche Geschlechtertrennung)":
        for value in values:
            if value in DICHOGAMY_VALUES:
                result["dichogamy"].add(DICHOGAMY_VALUES[value])
    elif label == "Diklinie (räumliche Geschlechtertrennung)":
        for value in values:
            folded = value.casefold()
            if "synözisch" in folded:
                result["sex_system"].add("hermaphroditic")
            elif "gynodiözisch" in folded:
                result["sex_system"].add("gynodioecious")
            elif "androdiözisch" in folded:
                result["sex_system"].add("androdioecious")
            elif "diözisch" in folded:
                result["sex_system"].add("dioecious")
            elif "monözisch" in folded:
                if any(prefix in folded for prefix in ("gyno", "andro", "tri")):
                    result["sex_system"].add("polygamous_or_other")
                else:
                    result["sex_system"].add("monoecious")
            elif "triözisch" in folded:
                result["sex_system"].add("polygamous_or_other")
    return result


def _coalesce_states(trait: str, states: set[str]) -> str:
    if trait == "self_incompatibility" and (
        "mixed_or_variable" in states or {"SI", "SC"}.issubset(states)
    ):
        return "mixed_or_variable"
    if trait == "mating_system" and len(states) > 1:
        return "mixed_mating"
    if trait == "autonomous_selfing_capacity" and len(states) > 1:
        return "mixed_or_variable"
    if trait == "cleistogamy" and len(states) > 1:
        return "facultative"
    if trait == "pollen_vector_mode" and len(states) > 1:
        return "mixed"
    if trait == "reward_type" and len(states) > 1:
        return "mixed_or_multiple"
    if trait == "dichogamy" and len(states) > 1:
        return "variable"
    if trait == "sex_system" and len(states) > 1:
        return "polygamous_or_other"
    return "|".join(sorted(states))


def _evidence_row(
    *,
    species: str,
    trait: str,
    states: set[str],
    quality: str,
    page: str,
    label: str,
    name_usage_id: str,
    excerpt: str,
    acceptance_contract: str,
) -> dict[str, str]:
    source_lineage, lineage_method = _lineage(page, label, name_usage_id)
    return {
        "accepted_species": species,
        "trait_name": trait,
        "normalized_value": _coalesce_states(trait, states),
        "quality": quality,
        "source_group": "floraweb_biolflor_bulk",
        "source_provider": "FloraWeb BfN via taxifydb frozen scrape",
        "source_url": _floraweb_url(page, name_usage_id),
        "source_record_id": f"floraweb:name-use-id:{name_usage_id}:{page}:{label}",
        "source_citation": SOURCE_CITATION,
        "source_excerpt": excerpt,
        "evidence_scope": "species_direct",
        "name_match_method": "accepted_name_exact",
        "source_lineage": source_lineage,
        "lineage_method": lineage_method,
        "source_run_id": f"external_release:{SNAPSHOT_RELEASE}",
        "source_artifact": f"{SNAPSHOT_RELEASE}:floraweb_raw_2026-06-24.csv",
        "acceptance_contract": acceptance_contract,
        "inference_rule": "",
    }


def build_source_package(
    raw: pd.DataFrame,
    master: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    required_raw = {"name_usage_id", "canonical_name", "page", "label", "value"}
    required_master = {"accepted_species", "family"}
    if missing := required_raw.difference(raw.columns):
        raise ValueError(f"raw snapshot missing columns: {sorted(missing)}")
    if missing := required_master.difference(master.columns):
        raise ValueError(f"master missing columns: {sorted(missing)}")

    raw = raw.fillna("").copy()
    master = master.fillna("").copy()
    master_names = set(master["accepted_species"])
    source_names = sorted(set(raw["canonical_name"].map(_text)))
    name_audit = pd.DataFrame(
        {
            "source_name": source_names,
            "accepted_species": [name if name in master_names else "" for name in source_names],
            "name_match_method": [
                "accepted_name_exact" if name in master_names else "unmatched_no_synonym_claim"
                for name in source_names
            ],
            "accepted": [name in master_names for name in source_names],
        }
    )
    matched = raw.loc[raw["canonical_name"].isin(master_names)].copy()

    strict_rows: list[dict[str, str]] = []
    independent_rows: list[dict[str, str]] = []
    extraction_audit: list[dict[str, Any]] = []

    for (species, name_usage_id, page, label), group in matched.groupby(
        ["canonical_name", "name_usage_id", "page", "label"],
        sort=True,
    ):
        values = sorted({_text(value) for value in group["value"] if _text(value)})
        excerpt = " | ".join(f"{label}: {value}" for value in values)
        mapped: dict[str, set[str]] = defaultdict(set)
        quality = "high"
        contract = "floraweb_exact_species_controlled_trait_v1"

        if label == "Befruchtungstyp":
            for value in values:
                if value in MATING_VALUES:
                    mapped["mating_system"].add(MATING_VALUES[value])
        elif label == "SI-Reaktion":
            for value in values:
                if value in SI_VALUES:
                    mapped["self_incompatibility"].add(SI_VALUES[value])
        elif label == "Bestäubung (Pollenvektoren)":
            for value in values:
                if value in AUTONOMY_VALUES:
                    mapped["autonomous_selfing_capacity"].add(AUTONOMY_VALUES[value])
                if value in CLEISTOGAMY_VALUES:
                    mapped["cleistogamy"].add(CLEISTOGAMY_VALUES[value])
        elif label == "Blumentyp (nach Kugler 1970)":
            for value in values:
                for marker, trait, state in FORM_PATTERNS:
                    if marker in value:
                        mapped[trait].add(state)
        elif label == "Beschreibung":
            states: set[str] = set()
            quotes: list[str] = []
            for value in values:
                value_states, quote = _colour_states(value)
                states.update(value_states)
                if quote:
                    quotes.append(quote)
            if states:
                mapped["flower_primary_color"].update(states)
                excerpt = " | ".join(quotes)
                quality = "medium"
                contract = "floraweb_exact_species_flower_bound_colour_v1"

        for trait, states in sorted(mapped.items()):
            if not states:
                continue
            strict_rows.append(
                _evidence_row(
                    species=species,
                    trait=trait,
                    states=states,
                    quality=quality,
                    page=page,
                    label=label,
                    name_usage_id=_text(name_usage_id),
                    excerpt=excerpt,
                    acceptance_contract=contract,
                )
            )

        independent = _independent_states(label, values)
        for trait, states in sorted(independent.items()):
            if not states:
                continue
            row = _evidence_row(
                species=species,
                trait=trait,
                states=states,
                quality="high",
                page=page,
                label=label,
                name_usage_id=_text(name_usage_id),
                excerpt=excerpt,
                acceptance_contract="floraweb_exact_species_independent_trait_v1",
            )
            row["strict_three_axis_included"] = "false"
            independent_rows.append(row)

        if mapped or independent:
            extraction_audit.append(
                {
                    "accepted_species": species,
                    "name_usage_id": _text(name_usage_id),
                    "page": page,
                    "label": label,
                    "raw_values": json.dumps(values, ensure_ascii=False),
                    "strict_traits": json.dumps(
                        {key: sorted(value) for key, value in mapped.items()},
                        ensure_ascii=False,
                        sort_keys=True,
                    ),
                    "independent_traits": json.dumps(
                        {key: sorted(value) for key, value in independent.items()},
                        ensure_ascii=False,
                        sort_keys=True,
                    ),
                    "source_url": _floraweb_url(page, _text(name_usage_id)),
                    "accepted": bool(mapped or independent),
                }
            )

    strict = pd.DataFrame(strict_rows, columns=STRICT_EVIDENCE_COLUMNS)
    independent = pd.DataFrame(
        independent_rows,
        columns=INDEPENDENT_EVIDENCE_COLUMNS,
    )
    audit = pd.DataFrame(extraction_audit)
    for frame in (strict, independent, audit, name_audit):
        if len(frame):
            frame.sort_values(list(frame.columns[: min(4, len(frame.columns))]), inplace=True)
            frame.reset_index(drop=True, inplace=True)
    return strict, independent, audit, name_audit


@app.command("build")
def build(
    raw_csv: Annotated[Path, typer.Option("--raw-csv", exists=True, dir_okay=False)],
    provenance_json: Annotated[
        Path, typer.Option("--provenance-json", exists=True, dir_okay=False)
    ],
    master_csv: Annotated[Path, typer.Option("--master-csv", exists=True, dir_okay=False)],
    output_dir: Annotated[Path, typer.Option("--output-dir", file_okay=False)],
    universe_csv: Annotated[
        Path | None, typer.Option("--universe-csv", exists=True, dir_okay=False)
    ] = None,
    expected_species: Annotated[int, typer.Option("--expected-species", min=1)] = 106_295,
) -> None:
    if _sha256(raw_csv) != RAW_SHA256:
        raise typer.BadParameter("FloraWeb raw snapshot SHA-256 mismatch")
    if _sha256(provenance_json) != PROVENANCE_SHA256:
        raise typer.BadParameter("FloraWeb provenance SHA-256 mismatch")
    provenance = json.loads(provenance_json.read_text(encoding="utf-8"))
    if provenance.get("accessed") != RETRIEVAL_DATE:
        raise typer.BadParameter("unexpected FloraWeb snapshot access date")

    raw = pd.read_csv(raw_csv, dtype=str).fillna("")
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    if universe_csv is not None:
        universe = pd.read_csv(
            universe_csv,
            usecols=["accepted_species"],
            dtype=str,
        ).fillna("")
        universe_names = set(universe["accepted_species"])
        if len(universe_names) != expected_species:
            raise typer.BadParameter(
                f"fixed universe has {len(universe_names)} species, expected {expected_species}"
            )
        master = master.loc[master["accepted_species"].isin(universe_names)].copy()
        if master["accepted_species"].nunique() != expected_species:
            raise typer.BadParameter("master does not contain the complete fixed universe")
    strict, independent, extraction_audit, name_audit = build_source_package(raw, master)
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "floraweb_reviewed_direct_evidence.csv.gz": strict,
        "floraweb_independent_trait_evidence.csv.gz": independent,
        "floraweb_extraction_audit.csv.gz": extraction_audit,
        "floraweb_name_match_audit.csv.gz": name_audit,
    }
    for name, frame in outputs.items():
        _write_gzip(frame, output_dir / name)

    summary = {
        "contract": "floraweb_biolflor_trait_source_package_v1",
        "source": {
            "snapshot_url": SNAPSHOT_URL,
            "snapshot_release": SNAPSHOT_RELEASE,
            "snapshot_access_date": RETRIEVAL_DATE,
            "raw_sha256": RAW_SHA256,
            "provenance_sha256": PROVENANCE_SHA256,
            "citation": SOURCE_CITATION,
        },
        "input": {
            "raw_rows": len(raw),
            "source_species": int(raw["canonical_name"].nunique()),
            "master_species": int(master["accepted_species"].nunique()),
        },
        "name_matching": {
            "exact_species": int(name_audit["accepted"].sum()),
            "unmatched_without_synonym_claim": int((~name_audit["accepted"]).sum()),
            "fuzzy_matches_accepted": 0,
        },
        "strict_direct": {
            "rows": len(strict),
            "species": int(strict["accepted_species"].nunique()),
            "species_trait": int(strict[["accepted_species", "trait_name"]].drop_duplicates().shape[0]),
            "by_trait": {
                str(key): int(value)
                for key, value in strict["trait_name"].value_counts().sort_index().items()
            },
            "by_quality": {
                str(key): int(value)
                for key, value in strict["quality"].value_counts().sort_index().items()
            },
        },
        "independent_traits": {
            "rows": len(independent),
            "species_trait": int(
                independent[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
            ),
            "by_trait": {
                str(key): int(value)
                for key, value in independent["trait_name"].value_counts().sort_index().items()
            },
            "strict_three_axis_included": False,
        },
        "checks": {
            "source_hashes_fixed": True,
            "accepted_names_exact_only": True,
            "source_lineage_underlying_not_provider_url": True,
            "pollen_vector_not_in_structure_axis": True,
            "reward_not_in_structure_axis": True,
            "family_inference_absent": True,
            "global_fallback_absent": True,
        },
    }
    summary["artifact_hashes"] = {
        name: _sha256(output_dir / name) for name in outputs
    }
    (output_dir / "floraweb_trait_source_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    app()
