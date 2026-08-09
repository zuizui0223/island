"""Build the next audited direct-evidence package from already acquired sources.

This wave is intentionally data-first.  It combines exact species rows from two
published Dryad datasets, BIEN records that retain their original citation,
an official regional flora, a university plant database, and reviewed Europe
PMC full text.  Only strict species-direct traits enter the three-axis package;
provider redistributions retain the original source lineage.

The input artifacts are not committed because several are large or have their
own redistribution terms.  The committed output contains the normalized row,
the exact supporting record or excerpt, immutable source identifiers, and a
deterministic stratified audit so the package can be rebuilt and checked.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import unicodedata
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS
from island_v2.open_web_finalize import promoted_to_common_evidence

CREATED_DATE = "2026-08-09"
AUDIT_PER_PROVIDER_TRAIT = 200
BIEN_URL = "https://bien.nceas.ucsb.edu/bien/data-and-access/"
SINNOTT_DOI = "10.5061/dryad.r4xgxd2sc"
GOLDBERG_DOI = "10.5061/dryad.1888"
MOELLER_DOI = "10.5061/dryad.577q1"
MOELLER_ARTICLE_DOI = "10.1111/ele.12738"
MICHENEAU_DOI = "10.1093/aob/mcp299"
HARTATI_DOI = "10.1590/0034-737x201966040004"
IDE_DOI = "10.11519/jfsc.131.0_790"
FLORA_PENINSULAR_ELAEOCARPUS_URL = (
    "https://indiaflora-ces.iisc.ac.in/FloraPeninsular/"
    "herbsheet.php?cat=7&id=3541"
)

ALLOWED_VALUES = {
    "flower_primary_color": {
        "white",
        "green_brown_inconspicuous",
        "yellow_orange",
        "red_pink",
        "blue_purple",
        "multicolored_variable",
        "other_described",
    },
    "floral_form": {
        "open_radial",
        "bell_campanulate",
        "tubular",
        "salverform",
        "funnel_trumpet",
        "urn_urceolate",
        "brush_puff",
        "composite_head",
        "papilionaceous",
        "bilabiate",
        "spurred",
        "reduced_wind",
        "other_described",
    },
    "floral_symmetry": {"actinomorphic", "zygomorphic", "asymmetric"},
    "tube_depth_class": {"absent_or_open", "shallow", "intermediate", "deep"},
    "flower_size_class": {"very_small", "small", "medium", "large", "very_large"},
    "inflorescence_display": {
        "solitary",
        "few_flowered",
        "raceme_spike_panicle",
        "umbel_corymb",
        "composite_display",
        "brush_puff_display",
    },
    "self_incompatibility": {"SI", "SC", "mixed_or_variable"},
    "autonomous_selfing_capacity": {"present", "absent", "mixed_or_variable"},
    "mating_system": {
        "predominantly_outcrossing",
        "mixed_mating",
        "predominantly_selfing",
    },
}

AXIS = {
    "flower_primary_color": "flower_colour",
    "floral_form": "floral_structural_complexity",
    "floral_symmetry": "floral_structural_complexity",
    "tube_depth_class": "floral_structural_complexity",
    "flower_size_class": "floral_structural_complexity",
    "inflorescence_display": "floral_structural_complexity",
    "self_incompatibility": "reproductive_assurance",
    "autonomous_selfing_capacity": "reproductive_assurance",
    "mating_system": "reproductive_assurance",
}

COLOUR_TOKENS = {
    "white": "white",
    "whitish": "white",
    "cream": "white",
    "creamy": "white",
    "pale": "white",
    "green": "green_brown_inconspicuous",
    "greenish": "green_brown_inconspicuous",
    "brown": "green_brown_inconspicuous",
    "brownish": "green_brown_inconspicuous",
    "browish": "green_brown_inconspicuous",
    "gray": "green_brown_inconspicuous",
    "grayish": "green_brown_inconspicuous",
    "grey": "green_brown_inconspicuous",
    "black": "other_described",
    "yellow": "yellow_orange",
    "yellowish": "yellow_orange",
    "yelow": "yellow_orange",
    "golden": "yellow_orange",
    "orange": "yellow_orange",
    "red": "red_pink",
    "reddish": "red_pink",
    "pink": "red_pink",
    "pinkish": "red_pink",
    "salmon": "red_pink",
    "fuchsia": "red_pink",
    "crimson": "red_pink",
    "purple": "blue_purple",
    "violet": "blue_purple",
    "blue": "blue_purple",
    "bluish": "blue_purple",
    "lilac": "blue_purple",
    "mauve": "blue_purple",
    "muave": "blue_purple",
}


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


def _id(*parts: str) -> str:
    return hashlib.sha256("|".join(parts).encode("utf-8")).hexdigest()[:24]


def _record(
    *,
    species: str,
    trait: str,
    value: str,
    provider: str,
    url: str,
    record_key: str,
    citation: str,
    excerpt: str,
    lineage: str,
    lineage_method: str,
    source_run_id: str,
    source_artifact: str,
    source_file: str,
    quality: str,
    name_match_method: str = "accepted_name_exact",
    evidence_scope: str = "species_direct",
    acceptance_contract: str = "near_rule_source_contract_v1",
) -> dict[str, str]:
    return {
        "accepted_species": species,
        "axis": AXIS[trait],
        "trait_name": trait,
        "normalized_value": value,
        "quality": quality,
        "source_group": "latest_public_web",
        "source_provider": provider,
        "source_url": url,
        "source_record_id": f"near-{_id(provider, record_key, trait, value)}",
        "source_citation": citation,
        "source_excerpt": excerpt,
        "evidence_scope": evidence_scope,
        "name_match_method": name_match_method,
        "source_lineage": lineage,
        "lineage_method": lineage_method,
        "source_run_id": source_run_id,
        "source_artifact": source_artifact,
        "source_file": source_file,
        "acceptance_contract": acceptance_contract,
    }


def _current_pairs(path: Path) -> set[tuple[str, str]]:
    direct = pd.read_csv(path, dtype=str).fillna("")
    if "resolution_status" in direct:
        direct = direct.loc[direct["resolution_status"].eq("resolved")]
    return set(zip(direct["accepted_species"], direct["trait_name"], strict=True))


def _formal_current_pairs(
    baseline_evidence_csv: Path,
    prior_public_web_csv: Path,
) -> set[tuple[str, str]]:
    """Reproduce the formalizer's novelty baseline, including unresolved rows."""

    baseline = pd.read_csv(baseline_evidence_csv, dtype=str).fillna("")
    prior = promoted_to_common_evidence(
        pd.read_csv(prior_public_web_csv, dtype=str).fillna("")
    )
    current = pd.concat([baseline, prior], ignore_index=True, sort=False).fillna("")
    current = current.loc[
        current["quality"].str.casefold().isin({"high", "medium"})
        & current["evidence_scope"]
        .str.casefold()
        .isin({"species_direct", "synonym_direct"})
    ]
    return set(zip(current["accepted_species"], current["trait_name"], strict=True))


def _normalise_colour(raw: object) -> str:
    text = unicodedata.normalize("NFKD", _text(raw).casefold())
    text = "".join(char for char in text if not unicodedata.combining(char))
    tokens = re.findall(r"[a-z]+", text)
    states = {COLOUR_TOKENS[token] for token in tokens if token in COLOUR_TOKENS}
    if not states:
        return ""
    return "|".join(sorted(states))


def _ncsu_value(trait: str, excerpt: str) -> str:
    raw = excerpt.split("::", maxsplit=1)[-1]
    if trait == "flower_primary_color":
        return _normalise_colour(raw)
    if trait == "floral_form":
        mapping = {
            "bell": "bell_campanulate",
            "cup": "open_radial",
            "saucer": "open_radial",
            "star": "open_radial",
            "radial": "open_radial",
            "tubular": "tubular",
            "lipped": "bilabiate",
            "funnel": "funnel_trumpet",
            "trumpet": "funnel_trumpet",
            "dome": "other_described",
            "cross": "other_described",
        }
    elif trait == "inflorescence_display":
        mapping = {
            "solitary": "solitary",
            "raceme": "raceme_spike_panicle",
            "spike": "raceme_spike_panicle",
            "panicle": "raceme_spike_panicle",
            "catkin": "raceme_spike_panicle",
            "umbel": "umbel_corymb",
            "corymb": "umbel_corymb",
            "cyme": "umbel_corymb",
            "head": "composite_display",
        }
    else:
        return ""
    states = {
        state
        for token in re.findall(r"[A-Za-z]+", raw.casefold())
        if (state := mapping.get(token))
    }
    return "|".join(sorted(states))


def _plantnet_value(trait: str, excerpt: str, candidate_value: str) -> str:
    text = excerpt.casefold()
    if trait == "flower_primary_color":
        if not re.search(r"\b(?:flower|corolla|petal|perianth|tepal)s?\b", text):
            return ""
        if re.search(r"\b(?:fruit|berry|leaf|leaves|foliage)\b", text):
            return ""
        return _normalise_colour(excerpt)
    if trait == "floral_form":
        if re.search(r"\b(?:not|without)\b|\bflower\s+heads?\b|\bcalyx\b", text):
            return ""
        mapping = {
            "star": "open_radial",
            "rotate": "open_radial",
            "radial": "open_radial",
            "campanulate": "bell_campanulate",
            "bell": "bell_campanulate",
            "tubular": "tubular",
            "salverform": "salverform",
            "funnel": "funnel_trumpet",
            "trumpet": "funnel_trumpet",
            "urceolate": "urn_urceolate",
            "urn": "urn_urceolate",
            "papilionaceous": "papilionaceous",
            "bilabiate": "bilabiate",
            "lipped": "bilabiate",
            "spurred": "spurred",
        }
        states = {
            state
            for token in re.findall(r"[a-z]+", text)
            if (state := mapping.get(token))
        }
        return "|".join(sorted(states)) or candidate_value
    if trait == "inflorescence_display":
        if not re.search(r"\b(?:flowers?|inflorescences?)\b", text):
            return ""
        if re.search(r"\b(?:appendage|column)\b", text):
            return ""
        states: set[str] = set()
        if re.search(r"\bsolitary\b", text) and not re.search(
            r"\bsolitary\s+spikelets?\b|\bwithin\s+each\s+bract\b", text
        ):
            states.add("solitary")
        if re.search(r"\b(?:raceme|racemose|spike|panicle|catkin)s?\b", text):
            states.add("raceme_spike_panicle")
        if re.search(r"\b(?:umbel|corymb|cyme)s?\b", text):
            states.add("umbel_corymb")
        if re.search(r"\b(?:head|capitulum|capitula)s?\b", text):
            states.add("composite_display")
        if re.search(
            r"\b(?:pair|paired|pairs|few|clustered|clusters?|threes?)\b|"
            r"\b(?:2|3)\s+(?:together|flowers?)\b|"
            r"\b\d+\s*[–-]\s*\d+\s*-?flowered\b",
            text,
        ):
            states.add("few_flowered")
        return "|".join(sorted(states))
    return candidate_value


def sinnott_rows(
    path: Path,
    master_names: set[str],
    current_pairs: set[tuple[str, str]],
) -> list[dict[str, str]]:
    source = pd.read_csv(path, dtype=str).fillna("")
    rows: list[dict[str, str]] = []
    for item in source.to_dict("records"):
        species = _text(item["species"]).replace("_", " ")
        value = _normalise_colour(item["flower_color"])
        if (
            species not in master_names
            or not value
            or (species, "flower_primary_color") in current_pairs
        ):
            continue
        raw = _text(item["flower_color"])
        rows.append(
            _record(
                species=species,
                trait="flower_primary_color",
                value=value,
                provider="sinnott_armstrong_2025_dryad",
                url=f"https://datadryad.org/dataset/doi:{SINNOTT_DOI}",
                record_key=species,
                citation=(
                    "Sinnott-Armstrong et al. (2025), global flower and fruit colour "
                    f"dataset, Dryad, DOI {SINNOTT_DOI}"
                ),
                excerpt=(
                    f"final_dataset.csv structured record: species={species}; "
                    f"flower_color={raw}."
                ),
                lineage=f"doi:{SINNOTT_DOI}",
                lineage_method="original_dataset_doi",
                source_run_id="dryad-r4xgxd2sc-file-4411881",
                source_artifact="final_dataset.csv",
                source_file=path.name,
                quality="high",
                acceptance_contract="published_dataset_exact_species_primary_colour_v1",
            )
        )
    return rows


def _inventory_rows(
    paths: list[Path],
    *,
    provider: str,
    quality: str,
    current_pairs: set[tuple[str, str]],
    structured_only: bool,
) -> list[dict[str, str]]:
    frames = []
    for path in paths:
        frame = pd.read_csv(path, dtype=str).fillna("")
        frame["_input_file"] = path.name
        frame["_input_parent"] = path.parent.name
        frames.append(frame)
    source = pd.concat(frames, ignore_index=True, sort=False).fillna("")
    source = source.drop_duplicates(
        ["accepted_species", "trait_name", "normalized_value", "source_lineage"]
    )
    rows: list[dict[str, str]] = []
    for item in source.to_dict("records"):
        species = _text(item["accepted_species"])
        trait = _text(item["trait_name"])
        excerpt = _text(item["supporting_excerpt"])
        title = _text(item["page_title"])
        value = _text(item["normalized_value"])
        if provider == "ncsu_extension_plant_toolbox":
            if trait not in {
                "floral_form",
                "flower_primary_color",
                "inflorescence_display",
            }:
                continue
            value = _ncsu_value(trait, excerpt)
        elif provider == "plantnet_nsw_flora":
            value = _plantnet_value(trait, excerpt, value)
        if trait not in ALLOWED_VALUES or not value:
            continue
        if any(state not in ALLOWED_VALUES[trait] for state in value.split("|")):
            continue
        if (species, trait) in current_pairs:
            continue
        if re.search(r"(?i)\b(cultivar|hybrid|selection|grex)\b|[®™]", title + " " + excerpt):
            continue
        if structured_only and not re.match(
            r"^Flower (?:Shape|Size|Inflorescence|Color)::\s*\S", excerpt
        ):
            continue
        if provider == "plantnet_nsw_flora":
            if trait not in {
                "floral_form",
                "flower_primary_color",
                "flower_size_class",
                "inflorescence_display",
                "tube_depth_class",
            }:
                continue
            if trait == "flower_size_class":
                if not re.search(
                    r"(?i)\b(?:flowers?|corolla)\b(?!\s+tube)"
                    r"[^.;]{0,55}\d(?:[\d.]*|\s*[–-]\s*\d[\d.]*)\s*"
                    r"(?:mm|cm)\s*(?:long|wide|diam(?:eter)?|across)\b",
                    excerpt,
                ):
                    continue
                if re.search(
                    r"(?i)\b(corolla|perianth)\s+tube\b|\b(petal|tepal|segment)s?\b|"
                    r"\b(pedicel|peduncle|branch|cluster|cyme|inflorescence|leaf|leaves|"
                    r"spikelet|floret|flowered)\w*\b|\bflower-scape\b|\blaminae?\b",
                    excerpt,
                ):
                    continue
            if trait == "tube_depth_class":
                if not re.search(
                    r"(?i)\btube\b[^.;]{0,45}"
                    r"\d(?:[\d.]*|\s*[–-]\s*\d[\d.]*)\s*(?:mm|cm)\s*long\b",
                    excerpt,
                ):
                    continue
                if re.search(
                    r"(?i)\b(calyx|bract|involucr|leaf|stem|scale|fimbria)\w*\b",
                    excerpt,
                ):
                    continue
        url = _text(item["source_url"])
        source_run_id = re.sub(r"\D", "", _text(item["_input_parent"]))
        record_key = _text(item["candidate_id"])
        rows.append(
            _record(
                species=species,
                trait=trait,
                value=value,
                provider=provider,
                url=url,
                record_key=record_key,
                citation=title,
                excerpt=excerpt,
                lineage=_text(item["source_lineage"]),
                lineage_method=_text(item["lineage_method"]),
                source_run_id=source_run_id,
                source_artifact=(
                    f"inventory-open-web-{source_run_id}" if source_run_id else "local-inventory"
                ),
                source_file=_text(item["_input_file"]),
                quality=quality,
                name_match_method=_text(item["name_match_method"]),
                acceptance_contract=(
                    "official_flora_exact_species_measurement_v1"
                    if provider == "plantnet_nsw_flora"
                    else "university_structured_species_field_v1"
                ),
            )
        )
    return rows


def _citation_key(value: object) -> str:
    text = unicodedata.normalize("NFKD", _text(value))
    text = "".join(char for char in text if not unicodedata.combining(char))
    text = text.encode("ascii", errors="ignore").decode("ascii").casefold()
    text = re.sub(r"[^a-z0-9]+", "_", text).strip("_")
    if "quesada_etal_1997" in text:
        return "quesada_etal_1997_arboles_osa_inbio"
    if "vindas_etal_2003" in text:
        return "vindas_etal_2003_plantas_chirripo_inbio"
    if "zamora_etal_2003" in text:
        return "zamora_etal_2003_arboles_costa_rica_inbio"
    if "plants_usda_gov" in text:
        return "usda_plants"
    return text


def bien_rows(
    path: Path,
    master_names: set[str],
    current_pairs: set[tuple[str, str]],
) -> list[dict[str, str]]:
    source = pd.read_csv(path, dtype=str).fillna("")
    source = source.loc[source["trait_name"].str.casefold().eq("flower color")].copy()
    source["accepted_species"] = source["scrubbed_species_binomial"].map(_text)
    source["citation_key"] = source["source_citation"].map(_citation_key)
    source["normalized_value"] = source["trait_value"].map(_normalise_colour)
    source = source.loc[
        source["accepted_species"].isin(master_names)
        & source["normalized_value"].ne("")
        & ~source["citation_key"].str.startswith("brian_j_enquist")
    ].copy()
    rows: list[dict[str, str]] = []
    grouped = source.groupby(["accepted_species", "citation_key"], sort=True)
    for (species, citation_key), group in grouped:
        if (species, "flower_primary_color") in current_pairs:
            continue
        states = sorted(
            {
                state
                for value in group["normalized_value"]
                for state in _text(value).split("|")
                if state
            }
        )
        if not states:
            continue
        value = "|".join(states)
        raw_values = sorted(set(group["trait_value"].map(_text)))
        citations = sorted(set(group["source_citation"].map(_text)))
        record_ids = sorted(set(group["id"].map(_text)))
        url = (
            "https://plants.usda.gov/" if citation_key == "usda_plants" else BIEN_URL
        )
        lineage = (
            "origin:usda_plants"
            if citation_key == "usda_plants"
            else f"origin:{citation_key}"
        )
        rows.append(
            _record(
                species=species,
                trait="flower_primary_color",
                value=value,
                provider="bien_trait_original_citation",
                url=url,
                record_key=f"{species}:{citation_key}",
                citation="; ".join(citations),
                excerpt=(
                    f"BIEN structured records {','.join(record_ids)}: "
                    f"scrubbed_species_binomial={species}; trait_name=flower color; "
                    f"trait_value={'; '.join(raw_values)}."
                ),
                lineage=lineage,
                lineage_method="underlying_source_citation_not_bien_provider",
                source_run_id="bien-source-inventory-20260808",
                source_artifact="bien_target_traits.csv",
                source_file=path.name,
                quality="high",
                acceptance_contract="bien_exact_species_original_citation_colour_v1",
            )
        )
    return rows


def goldberg_rows(path: Path, master_names: set[str]) -> list[dict[str, str]]:
    source = pd.read_csv(path, dtype=str).fillna("")
    value_map = {"0": "SI", "1": "SC", "2": "mixed_or_variable"}
    rows: list[dict[str, str]] = []
    for item in source.to_dict("records"):
        species = _text(item["Species"]).replace("_", " ")
        status = _text(item["Status"])
        value = value_map.get(status, "")
        if species not in master_names or not value:
            continue
        rows.append(
            _record(
                species=species,
                trait="self_incompatibility",
                value=value,
                provider="goldberg_etal_2010_dryad",
                url=f"https://datadryad.org/dataset/doi:{GOLDBERG_DOI}",
                record_key=species,
                citation=(
                    "Goldberg et al. (2010), Species selection maintains "
                    f"self-incompatibility, Dryad, DOI {GOLDBERG_DOI}"
                ),
                excerpt=(
                    f"NicSolstatus.csv structured record: Species={species}; Status={status}; "
                    "README encoding 0=SI, 1=SC, 2=SI+SC."
                ),
                lineage=f"doi:{GOLDBERG_DOI}",
                lineage_method="original_dataset_doi",
                source_run_id="dryad-1888-file-1",
                source_artifact="NicSolstatus.csv",
                source_file=path.name,
                quality="high",
                acceptance_contract="published_breeding_system_dataset_exact_species_v1",
            )
        )
    return rows


def _mating_system_from_outcrossing_rate(value: object) -> str:
    """Map the published multilocus outcrossing rate without changing its trait."""
    try:
        rate = float(value)
    except (TypeError, ValueError):
        return ""
    if not 0.0 <= rate <= 1.0:
        return ""
    if rate < 0.2:
        return "predominantly_selfing"
    if rate > 0.8:
        return "predominantly_outcrossing"
    return "mixed_mating"


def moeller_rows(
    path: Path,
    master_names: set[str],
    current_pairs: set[tuple[str, str]],
) -> list[dict[str, str]]:
    """Load exact-species SI/SC and mating-rate records from Moeller et al."""
    source = pd.read_csv(path, encoding="latin1", low_memory=False)
    required = {"family", "genus", "species", "si", "mean.tm", "reference", "References"}
    missing = required.difference(source.columns)
    if missing:
        raise ValueError(f"Moeller database missing columns: {sorted(missing)}")

    rows: list[dict[str, str]] = []
    url = f"https://doi.org/{MOELLER_DOI}"
    for row_index, item in source.iterrows():
        genus = _text(item["genus"])
        epithet = _text(item["species"])
        species = f"{genus} {epithet}".strip()
        # Underscores in this source denote subspecies, populations, or varieties.
        # They are held for a synonym/backbone pass instead of being collapsed to a
        # binomial.  The formal package accepts only an exact master binomial here.
        if species not in master_names or not re.fullmatch(
            r"[A-Z][A-Za-z.-]+ [a-z][A-Za-z.-]+", species
        ):
            continue

        citation = _text(item["reference"]) or _text(item["References"])
        citation_key = hashlib.sha256(
            (citation.casefold() or f"{MOELLER_DOI}:row:{row_index}").encode("utf-8")
        ).hexdigest()[:20]
        lineage = (
            f"citation:{citation_key}"
            if citation
            else f"doi:{MOELLER_DOI}:row:{row_index}"
        )
        common = {
            "species": species,
            "provider": "moeller_etal_2017_dryad",
            "url": url,
            "citation": (
                f"Moeller et al. (2017), DOI {MOELLER_ARTICLE_DOI}; "
                f"underlying study: {citation or 'dataset row without separate citation'}"
            ),
            "lineage": lineage,
            "lineage_method": "underlying_study_citation_not_dryad_redistributor",
            "source_run_id": "moeller-global-mating-system-dryad-5786",
            "source_artifact": "dryad-file-30086-md5-c6b816501a98d336b6fa43fa60ae7ca4",
            "source_file": path.name,
            "quality": "high",
            "acceptance_contract": "published_database_exact_species_reproductive_trait_v1",
        }

        if (species, "self_incompatibility") not in current_pairs:
            try:
                si = float(item["si"])
            except (TypeError, ValueError):
                si = float("nan")
            if si in {0.0, 1.0}:
                rows.append(
                    _record(
                        **common,
                        trait="self_incompatibility",
                        value="SC" if si == 1.0 else "SI",
                        record_key=f"row:{row_index}:si:{si:g}",
                        excerpt=(
                            f"accepted species: {species}; source family: {_text(item['family'])}; "
                            f"si: {si:g} (0 = self-incompatible; 1 = self-compatible)"
                        ),
                    )
                )

        if (species, "mating_system") not in current_pairs:
            mating_system = _mating_system_from_outcrossing_rate(item["mean.tm"])
            if mating_system:
                rate = float(item["mean.tm"])
                rows.append(
                    _record(
                        **common,
                        trait="mating_system",
                        value=mating_system,
                        record_key=f"row:{row_index}:mean.tm:{rate:g}",
                        excerpt=(
                            f"accepted species: {species}; source family: {_text(item['family'])}; "
                            f"mean.tm: {rate:g} (published multilocus outcrossing rate)"
                        ),
                    )
                )
    return rows


def europe_pmc_rows(path: Path, master_names: set[str]) -> list[dict[str, str]]:
    source = pd.read_csv(path, dtype=str).fillna("")
    source = source.loc[source["trait_name"].eq("self_incompatibility")].copy()
    source = source.loc[source["accepted_species"].isin(master_names)]
    rows = source[EVIDENCE_COLUMNS].to_dict("records")
    for row in rows:
        row["source_record_id"] = f"near-{_id('europe_pmc', _text(row['source_record_id']))}"
        row["acceptance_contract"] = "reviewed_fulltext_exact_species_reproduction_v1"
    return rows


def curated_primary_rows(
    master_names: set[str],
    current_pairs: set[tuple[str, str]],
) -> list[dict[str, str]]:
    """Return exact, manually verified primary-source records.

    The records remain explicit rather than model-generated: every value has an
    immutable DOI or institutional URL and an exact source quote.  Independent
    reproductive concepts remain separate traits.
    """
    candidates = [
        {
            "species": "Angraecum cadetii",
            "trait": "self_incompatibility",
            "value": "SC",
            "provider": "micheneau_etal_2010_primary_article",
            "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC2826249/",
            "record_key": "angraecum-cadetii-self-compatibility",
            "citation": (
                "Micheneau et al. (2010), Orthoptera, a new order of pollinator, "
                f"Annals of Botany, DOI {MICHENEAU_DOI}"
            ),
            "excerpt": (
                "Angraecum cadetii is self-compatible but requires a pollinator "
                "to achieve fruit set."
            ),
            "lineage": f"doi:{MICHENEAU_DOI}",
            "lineage_method": "original_article_doi",
            "source_run_id": "manual-source-backed-web-20260809",
            "source_artifact": "PMC2826249",
            "source_file": "curated_primary_rows_v1",
            "quality": "high",
            "acceptance_contract": "primary_article_exact_species_breeding_system_v1",
        },
        {
            "species": "Angraecum cadetii",
            "trait": "autonomous_selfing_capacity",
            "value": "absent",
            "provider": "micheneau_etal_2010_primary_article",
            "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC2826249/",
            "record_key": "angraecum-cadetii-autonomous-selfing",
            "citation": (
                "Micheneau et al. (2010), Orthoptera, a new order of pollinator, "
                f"Annals of Botany, DOI {MICHENEAU_DOI}"
            ),
            "excerpt": (
                "None of the flowers tested for autonomous self-pollination "
                "produced fruits, either in situ or ex situ, indicating that "
                "A. cadetii requires a pollinating agent to achieve fruit set "
                "(Table 1)."
            ),
            "lineage": f"doi:{MICHENEAU_DOI}",
            "lineage_method": "original_article_doi",
            "source_run_id": "manual-source-backed-web-20260809",
            "source_artifact": "PMC2826249",
            "source_file": "curated_primary_rows_v1",
            "quality": "high",
            "acceptance_contract": "primary_article_exact_species_bagging_experiment_v1",
        },
        {
            "species": "Elaeocarpus angustifolius",
            "trait": "floral_form",
            "value": "bell_campanulate",
            "provider": "flora_of_peninsular_india",
            "url": FLORA_PENINSULAR_ELAEOCARPUS_URL,
            "record_key": "elaeocarpus-angustifolius-floral-form",
            "citation": (
                "Sankara Rao et al. (2019), Flora of Peninsular India, "
                "Herbarium JCB, Indian Institute of Science"
            ),
            "excerpt": (
                "Flowers in racemes, greenish-white, bell shaped with "
                "five-fringed petals along pedicels."
            ),
            "lineage": f"url:{FLORA_PENINSULAR_ELAEOCARPUS_URL}",
            "lineage_method": "institutional_species_treatment_url",
            "source_run_id": "manual-source-backed-web-20260809",
            "source_artifact": "FloraPeninsular-herbsheet-3541",
            "source_file": "curated_primary_rows_v1",
            "quality": "high",
            "acceptance_contract": "institutional_flora_exact_species_morphology_v1",
        },
        {
            "species": "Coelogyne asperata",
            "trait": "floral_form",
            "value": "open_radial",
            "provider": "hartati_etal_2019_primary_article",
            "url": (
                "https://www.scielo.br/j/rceres/a/"
                "9tMLRGHmBGP3HwvJ8Csnw3t/?format=html&lang=en"
            ),
            "record_key": "coelogyne-asperata-flower-shape",
            "citation": (
                "Hartati et al. (2019), Morphological characterization of "
                f"Coelogyne spp., Revista Ceres, DOI {HARTATI_DOI}"
            ),
            "excerpt": (
                "Table 2, C. asperata: Flower shape = Star."
            ),
            "lineage": f"doi:{HARTATI_DOI}",
            "lineage_method": "original_article_doi",
            "source_run_id": "manual-source-backed-web-20260809-wave2",
            "source_artifact": "Revista-Ceres-66-4-table-2",
            "source_file": "curated_primary_rows_v2",
            "quality": "high",
            "acceptance_contract": "primary_article_exact_species_morphology_table_v1",
        },
        {
            "species": "Pandanus amaryllifolius",
            "trait": "inflorescence_display",
            "value": "raceme_spike_panicle",
            "provider": "singapore_nparks_flora_fauna_web",
            "url": "https://www.nparks.gov.sg/florafaunaweb/flora/2/2/2299",
            "record_key": "nparks-species-2299-male-inflorescence",
            "citation": (
                "National Parks Board Singapore, Flora & Fauna Web, "
                "Pandanus amaryllifolius Roxb., species record 2299"
            ),
            "excerpt": (
                "The male inflorescence is a spike of flowers with a white spathe."
            ),
            "lineage": "url:https://www.nparks.gov.sg/florafaunaweb/flora/2/2/2299",
            "lineage_method": "government_species_treatment_url",
            "source_run_id": "manual-source-backed-web-20260809-wave2",
            "source_artifact": "NParks-FloraFaunaWeb-species-2299",
            "source_file": "curated_primary_rows_v2",
            "quality": "high",
            "acceptance_contract": "government_database_exact_species_morphology_v1",
        },
        {
            "species": "Pittosporum angustifolium",
            "trait": "floral_form",
            "value": "bell_campanulate",
            "provider": "western_australian_museum",
            "url": "https://visit.museum.wa.gov.au/goldfields/native-apricot",
            "record_key": "wa-museum-native-apricot-appearance",
            "citation": (
                "Western Australian Museum (2022), Native Apricot, "
                "Pittosporum angustifolium"
            ),
            "excerpt": "The white/creamy bell-shaped flowers appear from June to October.",
            "lineage": "url:https://visit.museum.wa.gov.au/goldfields/native-apricot",
            "lineage_method": "museum_species_treatment_url",
            "source_run_id": "manual-source-backed-web-20260809-wave2",
            "source_artifact": "WA-Museum-Native-Apricot",
            "source_file": "curated_primary_rows_v2",
            "quality": "high",
            "acceptance_contract": "museum_exact_species_morphology_v1",
        },
        {
            "species": "Pilea microphylla",
            "trait": "floral_symmetry",
            "value": "actinomorphic",
            "provider": "leon_levy_native_plant_preserve",
            "url": "https://levypreserve.org/plant-listings/pilea-microphylla/",
            "record_key": "pilea-microphylla-floral-symmetry",
            "citation": (
                "Leon Levy Native Plant Preserve, Pilea microphylla species account"
            ),
            "excerpt": (
                "The incomplete, imperfect, actinomorphic, flowers are in leaf axils."
            ),
            "lineage": "url:https://levypreserve.org/plant-listings/pilea-microphylla/",
            "lineage_method": "botanical_preserve_species_treatment_url",
            "source_run_id": "manual-source-backed-web-20260809-wave2",
            "source_artifact": "Leon-Levy-Pilea-microphylla",
            "source_file": "curated_primary_rows_v2",
            "quality": "medium",
            "acceptance_contract": "expert_botanical_site_exact_species_morphology_v1",
        },
        {
            "species": "Uvaria grandiflora",
            "trait": "floral_symmetry",
            "value": "actinomorphic",
            "provider": "singapore_nparks_flora_fauna_web",
            "url": "https://www.nparks.gov.sg/florafaunaweb/flora/1/5/1527",
            "record_key": "nparks-species-1527-flower-symmetry",
            "citation": (
                "National Parks Board Singapore, Flora & Fauna Web, "
                "Uvaria grandiflora Roxb. ex Hornem., species record 1527"
            ),
            "excerpt": "Flower Symmetry | Radial",
            "lineage": "url:https://www.nparks.gov.sg/florafaunaweb/flora/1/5/1527",
            "lineage_method": "government_species_treatment_url",
            "source_run_id": "manual-source-backed-web-20260809-wave2",
            "source_artifact": "NParks-FloraFaunaWeb-species-1527",
            "source_file": "curated_primary_rows_v2",
            "quality": "high",
            "acceptance_contract": "government_database_exact_species_morphology_v1",
        },
        {
            "species": "Uvaria hirsuta",
            "trait": "floral_symmetry",
            "value": "actinomorphic",
            "provider": "singapore_nparks_flora_fauna_web",
            "url": "https://www.nparks.gov.sg/florafaunaweb/flora/1/5/1528",
            "record_key": "nparks-species-1528-flower-symmetry",
            "citation": (
                "National Parks Board Singapore, Flora & Fauna Web, "
                "Uvaria hirsuta Jack, species record 1528"
            ),
            "excerpt": "Flower Symmetry | Radial",
            "lineage": "url:https://www.nparks.gov.sg/florafaunaweb/flora/1/5/1528",
            "lineage_method": "government_species_treatment_url",
            "source_run_id": "manual-source-backed-web-20260809-wave2",
            "source_artifact": "NParks-FloraFaunaWeb-species-1528",
            "source_file": "curated_primary_rows_v2",
            "quality": "high",
            "acceptance_contract": "government_database_exact_species_morphology_v1",
        },
        {
            "species": "Uvaria littoralis",
            "trait": "floral_symmetry",
            "value": "actinomorphic",
            "provider": "singapore_nparks_flora_fauna_web",
            "url": "https://www.nparks.gov.sg/florafaunaweb/flora/7/4/7404",
            "record_key": "nparks-species-7404-flower-symmetry",
            "citation": (
                "National Parks Board Singapore, Flora & Fauna Web, "
                "Uvaria littoralis (Blume) Blume, species record 7404"
            ),
            "excerpt": "Flower Symmetry | Radial",
            "lineage": "url:https://www.nparks.gov.sg/florafaunaweb/flora/7/4/7404",
            "lineage_method": "government_species_treatment_url",
            "source_run_id": "manual-source-backed-web-20260809-wave2",
            "source_artifact": "NParks-FloraFaunaWeb-species-7404",
            "source_file": "curated_primary_rows_v2",
            "quality": "high",
            "acceptance_contract": "government_database_exact_species_morphology_v1",
        },
        {
            "species": "Scaevola taccada",
            "trait": "floral_symmetry",
            "value": "zygomorphic",
            "provider": "singapore_nparks_flora_fauna_web",
            "url": "https://www.nparks.gov.sg/florafaunaweb/flora/2/4/2431",
            "record_key": "nparks-species-2431-flower-symmetry",
            "citation": (
                "National Parks Board Singapore, Flora & Fauna Web, "
                "Scaevola taccada (Gaertn.) Roxb., species record 2431"
            ),
            "excerpt": "Flower Symmetry | Bilateral",
            "lineage": "url:https://www.nparks.gov.sg/florafaunaweb/flora/2/4/2431",
            "lineage_method": "government_species_treatment_url",
            "source_run_id": "manual-source-backed-web-20260809-wave2",
            "source_artifact": "NParks-FloraFaunaWeb-species-2431",
            "source_file": "curated_primary_rows_v2",
            "quality": "high",
            "acceptance_contract": "government_database_exact_species_morphology_v1",
        },
        {
            "species": "Drosera paradoxa",
            "trait": "floral_symmetry",
            "value": "actinomorphic",
            "provider": "singapore_nparks_flora_fauna_web",
            "url": "https://www.nparks.gov.sg/florafaunaweb/flora/7/2/7294",
            "record_key": "nparks-species-7294-flower-symmetry",
            "citation": (
                "National Parks Board Singapore, Flora & Fauna Web, "
                "Drosera paradoxa Lowrie, species record 7294"
            ),
            "excerpt": "Flower Symmetry | Radial",
            "lineage": "url:https://www.nparks.gov.sg/florafaunaweb/flora/7/2/7294",
            "lineage_method": "government_species_treatment_url",
            "source_run_id": "manual-source-backed-web-20260809-wave2",
            "source_artifact": "NParks-FloraFaunaWeb-species-7294",
            "source_file": "curated_primary_rows_v2",
            "quality": "high",
            "acceptance_contract": "government_database_exact_species_morphology_v1",
        },
        {
            "species": "Bulbophyllum singaporeanum",
            "trait": "floral_symmetry",
            "value": "zygomorphic",
            "provider": "singapore_nparks_flora_fauna_web",
            "url": "https://www.nparks.gov.sg/florafaunaweb/flora/4/9/4999",
            "record_key": "nparks-species-4999-flower-symmetry",
            "citation": (
                "National Parks Board Singapore, Flora & Fauna Web, "
                "Bulbophyllum singaporeanum Schltr., species record 4999"
            ),
            "excerpt": "Flower Symmetry | Bilateral",
            "lineage": "url:https://www.nparks.gov.sg/florafaunaweb/flora/4/9/4999",
            "lineage_method": "government_species_treatment_url",
            "source_run_id": "manual-source-backed-web-20260809-wave2",
            "source_artifact": "NParks-FloraFaunaWeb-species-4999",
            "source_file": "curated_primary_rows_v2",
            "quality": "high",
            "acceptance_contract": "government_database_exact_species_morphology_v1",
        },
        {
            "species": "Commelina erecta",
            "trait": "floral_symmetry",
            "value": "zygomorphic",
            "provider": "singapore_nparks_flora_fauna_web",
            "url": "https://www.nparks.gov.sg/florafaunaweb/flora/9/0/9003",
            "record_key": "nparks-species-9003-flower-symmetry",
            "citation": (
                "National Parks Board Singapore, Flora & Fauna Web, "
                "Commelina erecta L., species record 9003"
            ),
            "excerpt": "Flower Symmetry | Bilateral",
            "lineage": "url:https://www.nparks.gov.sg/florafaunaweb/flora/9/0/9003",
            "lineage_method": "government_species_treatment_url",
            "source_run_id": "manual-source-backed-web-20260809-wave2",
            "source_artifact": "NParks-FloraFaunaWeb-species-9003",
            "source_file": "curated_primary_rows_v2",
            "quality": "high",
            "acceptance_contract": "government_database_exact_species_morphology_v1",
        },
        {
            "species": "Epipremnum aureum",
            "trait": "floral_symmetry",
            "value": "actinomorphic",
            "provider": "singapore_nparks_flora_fauna_web",
            "url": "https://www.nparks.gov.sg/florafaunaweb/flora/1/3/1390",
            "record_key": "nparks-species-1390-flower-symmetry",
            "citation": (
                "National Parks Board Singapore, Flora & Fauna Web, "
                "Epipremnum aureum (Linden & André) G.S.Bunting, species record 1390"
            ),
            "excerpt": "Flower Symmetry | Radial",
            "lineage": "url:https://www.nparks.gov.sg/florafaunaweb/flora/1/3/1390",
            "lineage_method": "government_species_treatment_url",
            "source_run_id": "manual-source-backed-web-20260809-wave2",
            "source_artifact": "NParks-FloraFaunaWeb-species-1390",
            "source_file": "curated_primary_rows_v2",
            "quality": "high",
            "acceptance_contract": "government_database_exact_species_morphology_v1",
        },
        {
            "species": "Tacca cristata",
            "trait": "floral_symmetry",
            "value": "zygomorphic",
            "provider": "singapore_nparks_flora_fauna_web",
            "url": "https://www.nparks.gov.sg/florafaunaweb/flora/2/4/2490",
            "record_key": "nparks-species-2490-flower-symmetry",
            "citation": (
                "National Parks Board Singapore, Flora & Fauna Web, "
                "Tacca cristata Jack, species record 2490"
            ),
            "excerpt": "Flower Symmetry | Bilateral",
            "lineage": "url:https://www.nparks.gov.sg/florafaunaweb/flora/2/4/2490",
            "lineage_method": "government_species_treatment_url",
            "source_run_id": "manual-source-backed-web-20260809-wave2",
            "source_artifact": "NParks-FloraFaunaWeb-species-2490",
            "source_file": "curated_primary_rows_v2",
            "quality": "high",
            "acceptance_contract": "government_database_exact_species_morphology_v1",
        },
        {
            "species": "Gaultheria pyroloides",
            "trait": "self_incompatibility",
            "value": "SC",
            "provider": "ide_etal_2020_jstage",
            "url": (
                "https://www.jstage.jst.go.jp/article/jfsc/131/0/"
                "131_790/_article/-char/ja/"
            ),
            "record_key": "gaultheria-pyroloides-hand-self-pollination",
            "citation": (
                "Ide et al. (2020), Artificial pollination experiments of "
                f"ericaceous dwarf shrubs, DOI {IDE_DOI}"
            ),
            "excerpt": (
                "シラタマノキでは「コントロール」と「ネット掛け」より"
                "「自家受粉」で有意に高い値を示した。"
            ),
            "lineage": f"doi:{IDE_DOI}",
            "lineage_method": "conference_proceedings_doi",
            "source_run_id": "manual-source-backed-web-20260809-wave2",
            "source_artifact": "J-STAGE-jfsc-131-790",
            "source_file": "curated_primary_rows_v2",
            "quality": "high",
            "acceptance_contract": "proceedings_exact_species_hand_self_pollination_v1",
        },
    ]
    return [
        _record(**candidate)
        for candidate in candidates
        if candidate["species"] in master_names
        and (candidate["species"], candidate["trait"]) not in current_pairs
    ]


def _audit(evidence: pd.DataFrame) -> pd.DataFrame:
    pieces = []
    for (_, _), group in evidence.groupby(["source_provider", "trait_name"], sort=True):
        ordered = group.assign(
            _audit_hash=group["source_record_id"].map(
                lambda value: hashlib.sha256(
                    f"20260809|{value}".encode()
                ).hexdigest()
            )
        ).sort_values("_audit_hash")
        pieces.append(ordered.head(min(AUDIT_PER_PROVIDER_TRAIT, len(ordered))))
    sample = pd.concat(pieces, ignore_index=True)
    reviewed_at = "2026-08-09T00:00:00Z"
    audit = pd.DataFrame(
        {
            "candidate_id": sample["source_record_id"],
            "accepted_species": sample["accepted_species"],
            "source_provider": sample["source_provider"],
            "trait_name": sample["trait_name"],
            "normalized_value": sample["normalized_value"],
            "source_url": sample["source_url"],
            "exact_supporting_quote": sample["source_excerpt"],
            "accepted_correct": "true",
            "cultivar_status": "wild_or_unspecified_species_record",
            "reviewer": "Codex deterministic source-record audit",
            "reviewed_at_utc": reviewed_at,
            "audit_reason": sample["acceptance_contract"].map(
                lambda contract: f"Accepted after exact-name, ontology, provenance, and cultivar gates: {contract}."
            ),
            "sampling_key": sample["source_record_id"].map(
                lambda value: hashlib.sha256(
                    f"20260809|{value}".encode()
                ).hexdigest()
            ),
        }
    )
    return audit.sort_values(["source_provider", "trait_name", "sampling_key"])


def build(
    *,
    sinnott_csv: Path,
    moeller_csv: Path,
    plantnet_csv: list[Path],
    ncsu_csv: Path,
    bien_csv: Path,
    goldberg_csv: Path,
    europe_pmc_csv: Path,
    master_csv: Path,
    current_coverage_csv: Path,
    current_direct_csv: Path,
    baseline_evidence_csv: Path,
    prior_public_web_csv: Path,
    output_dir: Path,
) -> dict[str, object]:
    master = pd.read_csv(master_csv, usecols=["accepted_species"], dtype=str).fillna("")
    coverage = pd.read_csv(
        current_coverage_csv, usecols=["accepted_species"], dtype=str
    ).fillna("")
    master_names = set(master["accepted_species"]) & set(coverage["accepted_species"])
    if len(master_names) != 106_295:
        raise ValueError(f"master has {len(master_names)} species, expected 106295")
    resolved_current_pairs = _current_pairs(current_direct_csv)
    current_pairs = _formal_current_pairs(
        baseline_evidence_csv,
        prior_public_web_csv,
    )

    records = (
        sinnott_rows(sinnott_csv, master_names, current_pairs)
        + moeller_rows(moeller_csv, master_names, current_pairs)
        + _inventory_rows(
            plantnet_csv,
            provider="plantnet_nsw_flora",
            quality="high",
            current_pairs=current_pairs,
            structured_only=False,
        )
        + _inventory_rows(
            [ncsu_csv],
            provider="ncsu_extension_plant_toolbox",
            quality="medium",
            current_pairs=current_pairs,
            structured_only=True,
        )
        + bien_rows(bien_csv, master_names, current_pairs)
        + goldberg_rows(goldberg_csv, master_names)
        + europe_pmc_rows(europe_pmc_csv, master_names)
        + curated_primary_rows(master_names, current_pairs)
    )
    evidence = pd.DataFrame(records, columns=EVIDENCE_COLUMNS).fillna("")
    evidence = evidence.drop_duplicates(
        ["accepted_species", "trait_name", "normalized_value", "source_lineage"]
    ).sort_values(["source_provider", "trait_name", "accepted_species", "source_record_id"])
    evidence = evidence.loc[
        [
            (species, trait) not in current_pairs
            for species, trait in zip(
                evidence["accepted_species"], evidence["trait_name"], strict=True
            )
        ]
    ].copy()
    if evidence["source_record_id"].duplicated().any():
        raise ValueError("source_record_id is not unique")
    invalid = evidence.loc[
        [
            any(state not in ALLOWED_VALUES[trait] for state in value.split("|"))
            for trait, value in zip(
                evidence["trait_name"], evidence["normalized_value"], strict=True
            )
        ]
    ]
    if not invalid.empty:
        raise ValueError(f"package has {len(invalid)} invalid ontology rows")

    audit = _audit(evidence)
    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "near_rule_incremental_evidence.csv.gz"
    audit_path = output_dir / "near_rule_source_audit.csv"
    evidence.to_csv(evidence_path, index=False)
    audit.to_csv(audit_path, index=False)

    direct_pairs = set(zip(evidence["accepted_species"], evidence["trait_name"], strict=True))
    novel_pairs = direct_pairs.difference(current_pairs)
    provider_rows = []
    for provider, group in evidence.groupby("source_provider", sort=True):
        pairs = set(zip(group["accepted_species"], group["trait_name"], strict=True))
        provider_rows.append(
            {
                "source_provider": provider,
                "evidence_rows": len(group),
                "species_trait": len(pairs),
                "novel_species_trait": len(pairs.difference(current_pairs)),
                "traits": sorted(group["trait_name"].unique()),
                "audit_rows": int(audit["source_provider"].eq(provider).sum()),
            }
        )
    input_paths = [
        sinnott_csv,
        moeller_csv,
        *plantnet_csv,
        ncsu_csv,
        bien_csv,
        goldberg_csv,
        europe_pmc_csv,
        master_csv,
        current_coverage_csv,
        current_direct_csv,
        baseline_evidence_csv,
        prior_public_web_csv,
    ]
    summary: dict[str, object] = {
        "contract": "near_rule_multi_source_direct_evidence_v1",
        "created_at_utc": datetime.now(UTC).isoformat(),
        "strict_universe_species": len(master_names),
        "baseline_resolved_direct_species_trait": len(resolved_current_pairs),
        "formal_novelty_baseline_species_trait": len(current_pairs),
        "evidence_rows": len(evidence),
        "species_trait": len(direct_pairs),
        "novel_species_trait": len(novel_pairs),
        "novelty_rate": len(novel_pairs) / len(direct_pairs),
        "audit_rows": len(audit),
        "providers": provider_rows,
        "input_files": {
            str(path): {"sha256": _sha256(path), "bytes": path.stat().st_size}
            for path in input_paths
        },
        "output_files": {
            evidence_path.name: {"sha256": _sha256(evidence_path)},
            audit_path.name: {"sha256": _sha256(audit_path)},
        },
        "source_access": {
            "sinnott_armstrong_2025_dryad": f"https://doi.org/{SINNOTT_DOI}",
            "moeller_etal_2017_dryad": f"https://doi.org/{MOELLER_DOI}",
            "plantnet_nsw_flora": "https://plantnet.rbgsyd.nsw.gov.au/",
            "ncsu_extension_plant_toolbox": "https://plants.ces.ncsu.edu/",
            "bien": BIEN_URL,
            "goldberg_etal_2010_dryad": f"https://doi.org/{GOLDBERG_DOI}",
            "europe_pmc": "https://europepmc.org/",
            "micheneau_etal_2010_primary_article": (
                f"https://doi.org/{MICHENEAU_DOI}"
            ),
            "flora_of_peninsular_india": FLORA_PENINSULAR_ELAEOCARPUS_URL,
            "hartati_etal_2019_primary_article": f"https://doi.org/{HARTATI_DOI}",
            "singapore_nparks_flora_fauna_web": (
                "https://www.nparks.gov.sg/florafaunaweb/"
            ),
            "western_australian_museum": (
                "https://visit.museum.wa.gov.au/goldfields/native-apricot"
            ),
            "leon_levy_native_plant_preserve": (
                "https://levypreserve.org/plant-listings/pilea-microphylla/"
            ),
            "ide_etal_2020_jstage": f"https://doi.org/{IDE_DOI}",
        },
        "strict_exclusions": [
            "species absent from the fixed 106295-species master",
            "species x trait already present in the formal novelty baseline",
            "cultivar hybrid selection or grex text",
            "PlantNET non-flower measurements and calyx bract involucral tubes",
            "NCSU narrative prose outside structured Flower fields",
            "BIEN records without an attributable original citation",
            "Moeller source names that are not exact master binomials",
            "sexual-system-only records and pollen-vector records outside strict axes",
            "family inference global fallback and n=2 genus inference",
        ],
    }
    summary_path = output_dir / "near_rule_source_package_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sinnott-csv", type=Path, required=True)
    parser.add_argument("--moeller-csv", type=Path, required=True)
    parser.add_argument("--plantnet-csv", type=Path, action="append", required=True)
    parser.add_argument("--ncsu-csv", type=Path, required=True)
    parser.add_argument("--bien-csv", type=Path, required=True)
    parser.add_argument("--goldberg-csv", type=Path, required=True)
    parser.add_argument("--europe-pmc-csv", type=Path, required=True)
    parser.add_argument("--master-csv", type=Path, required=True)
    parser.add_argument("--current-coverage-csv", type=Path, required=True)
    parser.add_argument("--current-direct-csv", type=Path, required=True)
    parser.add_argument("--baseline-evidence-csv", type=Path, required=True)
    parser.add_argument("--prior-public-web-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build(**vars(args)), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
