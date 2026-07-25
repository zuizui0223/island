"""Free, source-backed search and extraction for the legacy v1 trait table.

The collector emulates the useful part of the former manual ChatGPT workflow:
search several public sources with trait-specific intent, retain the supporting
text, and then collapse only explicit matches into the nine-column v1 schema.
It does not call an LLM or a paid API.
"""

from __future__ import annotations

import hashlib
import json
import re
import time
from collections.abc import Callable, Iterable
from concurrent.futures import ThreadPoolExecutor
from html.parser import HTMLParser
from pathlib import Path
from typing import Any
from urllib.parse import quote, quote_plus, urljoin
from xml.etree import ElementTree

import pandas as pd
import typer

from island_v2 import reported_ecology_machine, trait_descriptive_scout, trait_web_reported_scout
from island_v2.europe_pmc_discovery import (
    discover_europe_pmc_full_text_sources,
    discover_europe_pmc_sources,
)
from island_v2.v1_category_traits import (
    OUTPUT_COLUMNS,
    derive_v1_categories,
    validate_result_table,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

SOURCE_COLUMNS = [
    "accepted_species",
    "source_text",
    "source_url",
    "source_citation",
    "source_type",
    "evidence_scope",
]
TRAIT_TEXT_MARKERS = re.compile(
    r"flower|floral|corolla|petal|blossom|pollinat|bee|bumble|fly|moth|"
    r"butterfl|bird|wind|selfing|outcross|mating|breeding system|"
    r"self[- ](?:in)?compatib|autogam|cleistogam",
    flags=re.IGNORECASE,
)

WEB_DESCRIPTION_SOURCES = [
    {
        "name": "ncsu_plant_toolbox",
        "source_type": "horticulture_site",
        "url_style": "slug",
        "url_template": "https://plants.ces.ncsu.edu/plants/{name}/",
        "sitemap_url": "https://plants.ces.ncsu.edu/sitemap.xml",
        "required_marker": "Plant Detail",
    },
    {
        "name": "nzpcn_flora",
        "source_type": "flora_or_monograph",
        "url_style": "slug",
        "url_template": "https://www.nzpcn.org.nz/flora/species/{name}/",
        "sitemap_url": "https://www.nzpcn.org.nz/sitemap.xml",
        "required_marker": "Detailed description",
    },
    {
        "name": "pfaf",
        "source_type": "horticulture_site",
        "url_style": "query",
        "url_template": "https://pfaf.org/user/Plant.aspx?LatinName={name}",
        "sitemap_url": "https://pfaf.org/sitemap.xml",
        "required_marker": "Physical Characteristics",
        "visible_text_mode": "aspnet_form",
    },
    {
        "name": "useful_tropical_plants",
        "source_type": "horticulture_site",
        "url_style": "query",
        "url_template": "https://tropical.theferns.info/viewtropical.php?id={name}",
        "required_marker": "General Information",
        "reject_marker": "We have no record for",
    },
]
WEB_DESCRIPTION_INDEX_VERSION = "public_plant_site_sitemaps_v1"

EVIDENCE_COLUMNS = [
    "species",
    "field",
    "value",
    "matched_term",
    "source_type",
    "source_url",
    "source_citation",
    "source_excerpt",
    "evidence_type",
    "confidence",
    "source_tier",
    "rule_id",
]

LIKELY_EVIDENCE_COLUMNS = [
    "species",
    "field",
    "value",
    "basis_fields",
    "basis_values",
    "basis_source_urls",
    "rule_id",
    "confidence",
]

LIKELY_WIDE_COLUMNS = [
    "species",
    "likely_pollination_guild",
    "likely_selfing_status",
    "likely_self_incompatibility",
    "likely_basis",
    "likely_rule_ids",
    "likely_confidence",
]

PRIORITY_COLUMNS = [
    "species",
    "flower_color",
    "flower_shape",
    "pollination_guild",
    "pollination_guild_mode",
    "mating_system",
    "selfing_status",
    "selfing_status_mode",
    "self_incompatibility",
    "self_incompatibility_mode",
    "basis_rule_ids",
]

V2_REPORTED_COLUMNS = [
    "accepted_species",
    "trait_layer",
    "trait_name",
    "provisional_candidate_value",
    "matched_term",
    "candidate_class",
    "source",
    "source_type",
    "source_url",
    "raw_description",
    "evidence_scope",
    "wild_or_cultivated",
    "review_status",
]

V2_FIELD_VALUE_MAP: dict[tuple[str, str], tuple[str, str, str]] = {
    ("flower_color", "white"): ("M0_floral_phenotype", "flower_primary_color", "white"),
    ("flower_color", "cream"): ("M0_floral_phenotype", "flower_primary_color", "white"),
    ("flower_color", "yellow"): ("M0_floral_phenotype", "flower_primary_color", "yellow_orange"),
    ("flower_color", "orange"): ("M0_floral_phenotype", "flower_primary_color", "yellow_orange"),
    ("flower_color", "red"): ("M0_floral_phenotype", "flower_primary_color", "red_pink"),
    ("flower_color", "pink"): ("M0_floral_phenotype", "flower_primary_color", "red_pink"),
    ("flower_color", "purple"): ("M0_floral_phenotype", "flower_primary_color", "blue_purple"),
    ("flower_color", "blue"): ("M0_floral_phenotype", "flower_primary_color", "blue_purple"),
    ("flower_color", "green"): (
        "M0_floral_phenotype",
        "flower_primary_color",
        "green_brown_inconspicuous",
    ),
    ("flower_color", "brown"): (
        "M0_floral_phenotype",
        "flower_primary_color",
        "green_brown_inconspicuous",
    ),
    ("flower_shape", "actinomorphic / radially symmetric"): (
        "M0_floral_phenotype",
        "floral_symmetry",
        "actinomorphic",
    ),
    ("flower_shape", "zygomorphic / bilaterally symmetric"): (
        "M0_floral_phenotype",
        "floral_symmetry",
        "zygomorphic",
    ),
    ("flower_shape", "tubular"): ("M0_floral_phenotype", "floral_form", "tubular"),
    ("flower_shape", "bell-shaped"): ("M0_floral_phenotype", "floral_form", "bell_campanulate"),
    ("flower_shape", "funnel / trumpet-shaped"): (
        "M0_floral_phenotype",
        "floral_form",
        "funnel_trumpet",
    ),
    ("flower_shape", "open / bowl / cup-shaped"): (
        "M0_floral_phenotype",
        "floral_form",
        "open_radial",
    ),
    ("flower_shape", "composite head / capitulum"): (
        "M0_floral_phenotype",
        "floral_form",
        "composite_head",
    ),
    ("flower_shape", "papilionaceous"): ("M0_floral_phenotype", "floral_form", "papilionaceous"),
    ("flower_shape", "spurred"): ("M0_floral_phenotype", "floral_form", "spurred"),
    ("flower_shape", "urn-shaped"): ("M0_floral_phenotype", "floral_form", "urn_urceolate"),
    ("flower_shape", "globose / spherical flower head"): (
        "M0_floral_phenotype",
        "floral_form",
        "composite_head",
    ),
    ("flower_shape", "spike / elongated inflorescence"): (
        "M0_floral_phenotype",
        "inflorescence_display",
        "raceme_spike_panicle",
    ),
    ("pollination_guild", "bumblebees"): (
        "M2_pollinator_functional_evidence",
        "pollination_functional_guild",
        "bumblebees",
    ),
    ("pollination_guild", "bees"): (
        "M2_pollinator_functional_evidence",
        "pollination_functional_guild",
        "other_bees",
    ),
    ("pollination_guild", "flies"): (
        "M2_pollinator_functional_evidence",
        "pollination_functional_guild",
        "flies",
    ),
    ("pollination_guild", "butterflies"): (
        "M2_pollinator_functional_evidence",
        "pollination_functional_guild",
        "butterflies",
    ),
    ("pollination_guild", "moths"): (
        "M2_pollinator_functional_evidence",
        "pollination_functional_guild",
        "moths",
    ),
    ("pollination_guild", "birds"): (
        "M2_pollinator_functional_evidence",
        "pollination_functional_guild",
        "birds",
    ),
    ("pollination_guild", "wind"): (
        "M1_reproductive_assurance_pollen_vector",
        "pollen_vector_mode",
        "abiotic_wind",
    ),
    ("pollination_guild", "self"): (
        "M1_reproductive_assurance_pollen_vector",
        "autonomous_selfing_capacity",
        "autonomous",
    ),
    ("mating_system", "obligate_outcrossing"): (
        "M1_reproductive_assurance_pollen_vector",
        "mating_system",
        "predominantly_outcrossing",
    ),
    ("mating_system", "mixed_mating"): (
        "M1_reproductive_assurance_pollen_vector",
        "mating_system",
        "mixed_mating",
    ),
    ("mating_system", "mainly_selfing"): (
        "M1_reproductive_assurance_pollen_vector",
        "mating_system",
        "predominantly_selfing",
    ),
    ("mating_system", "obligate_selfing"): (
        "M1_reproductive_assurance_pollen_vector",
        "mating_system",
        "obligate_selfing",
    ),
    ("self_incompatibility", "SI"): (
        "M1_reproductive_assurance_pollen_vector",
        "self_incompatibility",
        "SI",
    ),
    ("self_incompatibility", "SC"): (
        "M1_reproductive_assurance_pollen_vector",
        "self_incompatibility",
        "SC",
    ),
    ("selfing_signal", "autonomous_selfing"): (
        "M1_reproductive_assurance_pollen_vector",
        "autonomous_selfing_capacity",
        "autonomous",
    ),
}

SOURCE_TIER_ORDER = {
    "A_flora": 0,
    "A_literature": 0,
    "B_database": 1,
    "B_encyclopedia": 2,
    "C_horticulture": 3,
    "D_inference": 4,
}

Rule = tuple[str, str, str, str]

COLOR_RULES: list[Rule] = [
    ("flower_color", "white", r"\b(?:white|whitish|ivory)\b", "color_white"),
    ("flower_color", "cream", r"\b(?:cream|creamy)\b", "color_cream"),
    ("flower_color", "yellow", r"\b(?:yellow|yellowish|golden)\b", "color_yellow"),
    ("flower_color", "orange", r"\borange\b", "color_orange"),
    ("flower_color", "red", r"\b(?:red|reddish|crimson|scarlet)\b", "color_red"),
    ("flower_color", "pink", r"\b(?:pink|pinkish|rose[- ]colou?red)\b", "color_pink"),
    (
        "flower_color",
        "purple",
        r"\b(?:purple|purplish|violet|mauve|lilac|lavender)\b",
        "color_purple",
    ),
    ("flower_color", "blue", r"\b(?:blue|bluish)\b", "color_blue"),
    ("flower_color", "green", r"\b(?:green|greenish)\b", "color_green"),
    ("flower_color", "brown", r"\b(?:brown|brownish)\b", "color_brown"),
]

SHAPE_RULES: list[Rule] = [
    (
        "flower_shape",
        "actinomorphic / radially symmetric",
        r"\b(?:actinomorphic|radially symmetric|radial symmetry)\b",
        "shape_actinomorphic",
    ),
    (
        "flower_shape",
        "zygomorphic / bilaterally symmetric",
        r"\b(?:zygomorphic|bilaterally symmetric|bilateral symmetry|bilabiate|two[- ]lipped)\b",
        "shape_zygomorphic",
    ),
    (
        "flower_shape",
        "tubular",
        r"\b(?:tubular|tube[- ]shaped|salverform|salver[- ]shaped|corolla[- ]tube)\b",
        "shape_tubular",
    ),
    (
        "flower_shape",
        "bell-shaped",
        r"\b(?:campanulate|bell[- ]shaped)\b",
        "shape_bell",
    ),
    (
        "flower_shape",
        "funnel / trumpet-shaped",
        r"\b(?:funnelform|funnel[- ]shaped|trumpet[- ]shaped|infundibuliform)\b|"
        r"\bflower shape\s*:?\s*(?:funnel|trumpet)\b",
        "shape_funnel",
    ),
    (
        "flower_shape",
        "open / bowl / cup-shaped",
        r"\b(?:bowl[- ]shaped|cup[- ]shaped|rotate flower|open[- ]faced flower)\b",
        "shape_open",
    ),
    (
        "flower_shape",
        "papilionaceous",
        r"\b(?:papilionaceous|pea[- ]shaped)\b",
        "shape_papilionaceous",
    ),
    (
        "flower_shape",
        "composite head / capitulum",
        r"\b(?:capitulum|capitula|composite flower head|flower head)\b",
        "shape_head",
    ),
    (
        "flower_shape",
        "spurred",
        r"\b(?:spurred flower|floral spur|nectar spur)\b",
        "shape_spurred",
    ),
    (
        "flower_shape",
        "globose / spherical flower head",
        r"\b(?:(?:spherical|globose|globbose|globe[- ]shaped) "
        r"(?:flower )?(?:heads?|balls?)|(?:flower )?(?:heads?|balls?) "
        r"(?:are )?(?:spherical|globose|globbose|globe[- ]shaped))\b",
        "shape_globose_head",
    ),
    (
        "flower_shape",
        "spike / elongated inflorescence",
        r"\b(?:inflorescences? (?:are )?spiciform|spiciform inflorescences?|"
        r"spicate inflorescences?|floral spikes?|flower spikes?)\b",
        "shape_spike_inflorescence",
    ),
]

GUILD_RULES: list[Rule] = [
    (
        "pollination_guild",
        "bumblebees",
        r"\b(?:pollinated by bumblebees?|bumblebee[- ]pollinated|Bombus[- ]pollinated)\b",
        "guild_bumblebees",
    ),
    (
        "pollination_guild",
        "bees",
        r"\b(?:pollinated by bees?|bee[- ]pollinated|melittophil(?:y|ous))\b",
        "guild_bees",
    ),
    (
        "pollination_guild",
        "bees",
        r"\b(?:attracts?|visited by) (?:honeybees?|bees?)\b",
        "guild_attracts_bees",
    ),
    (
        "pollination_guild",
        "flies",
        r"\b(?:pollinated by fl(?:y|ies)|fl(?:y|ies)[- ]pollinated|myi?ophil(?:y|ous))\b",
        "guild_flies",
    ),
    (
        "pollination_guild",
        "butterflies",
        r"\b(?:pollinated by (?:bees? (?:and|or) )?butterflies|"
        r"butterfl(?:y|ies)[- ]pollinated|psychophil(?:y|ous))\b",
        "guild_butterflies",
    ),
    (
        "pollination_guild",
        "butterflies",
        r"\b(?:attracts?|visited by) (?:adult )?butterflies\b",
        "guild_attracts_butterflies",
    ),
    (
        "pollination_guild",
        "moths",
        r"\b(?:pollinated by moths?|moth[- ]pollinated|phalaenophil(?:y|ous))\b",
        "guild_moths",
    ),
    (
        "pollination_guild",
        "moths",
        r"\b(?:attracts?|visited by) moths?\b",
        "guild_attracts_moths",
    ),
    (
        "pollination_guild",
        "birds",
        r"\b(?:pollinated by birds?|bird[- ]pollinated|hummingbird[- ]pollinated|ornithophil(?:y|ous))\b",
        "guild_birds",
    ),
    (
        "pollination_guild",
        "birds",
        r"\b(?:attracts?|visited by) (?:birds?|hummingbirds?)\b",
        "guild_attracts_birds",
    ),
    (
        "pollination_guild",
        "wind",
        r"\b(?:pollinated by (?:the )?wind|wind[- ]pollinated|wind pollination|anemophil(?:y|ous))\b",
        "guild_wind",
    ),
    (
        "pollination_guild",
        "self",
        r"\b(?:self[- ]pollinated|self[- ]pollinating|automatic self[- ]pollination)\b",
        "guild_self",
    ),
]

MATING_RULES: list[Rule] = [
    (
        "mating_system",
        "obligate_outcrossing",
        r"\b(?:obligate(?:ly)? outcross(?:ing|er)|obligate allogam(?:y|ous))\b",
        "mating_obligate_outcrossing",
    ),
    (
        "mating_system",
        "mixed_mating",
        r"\b(?:mixed mating|facultative selfing|both selfing and outcrossing)\b",
        "mating_mixed",
    ),
    (
        "mating_system",
        "mainly_selfing",
        r"\b(?:predominantly|primarily|mainly|largely) self(?:ing|[- ]fertili[sz]ing)\b",
        "mating_mainly_selfing",
    ),
    (
        "mating_system",
        "obligate_selfing",
        r"\b(?:obligate(?:ly)? self(?:ing|er)|obligate autogam(?:y|ous))\b",
        "mating_obligate_selfing",
    ),
]

SI_RULES: list[Rule] = [
    (
        "self_incompatibility",
        "SI",
        r"\bself[- ]incompatib(?:le|ility)\b",
        "si_explicit",
    ),
    (
        "self_incompatibility",
        "SC",
        r"\bself[- ]compatib(?:le|ility)\b",
        "sc_explicit",
    ),
    (
        "self_incompatibility",
        "likely_SC",
        r"\b(?:autonomous|spontaneous) selfing\b|\bautogam(?:y|ous)\b",
        "sc_from_explicit_autogamy",
    ),
]

VISITATION_RULES: list[Rule] = [
    (
        "pollinator_visitation",
        guild,
        pattern,
        f"visitation_{guild}",
    )
    for guild, pattern in (
        ("bees", r"\b(?:attracts?|visited by) (?:honeybees?|bees?)\b"),
        ("butterflies", r"\b(?:attracts?|visited by) (?:adult )?butterflies\b"),
        ("moths", r"\b(?:attracts?|visited by) moths?\b"),
        ("birds", r"\b(?:attracts?|visited by) (?:birds?|hummingbirds?)\b"),
    )
]

SELFING_SIGNAL_RULES: list[Rule] = [
    (
        "selfing_signal",
        "autonomous_selfing",
        r"\b(?:autonomous|spontaneous) selfing\b|\bautogam(?:y|ous)\b|"
        r"\bautofert(?:ile|ility)\b",
        "selfing_signal_autonomous",
    ),
    (
        "selfing_signal",
        "cleistogamy",
        r"\bcleistogam(?:y|ous)\b|\bcleistogamous flowers?\b",
        "selfing_signal_cleistogamy",
    ),
]

ALL_RULES = [
    *COLOR_RULES,
    *SHAPE_RULES,
    *GUILD_RULES,
    *MATING_RULES,
    *SI_RULES,
    *VISITATION_RULES,
    *SELFING_SIGNAL_RULES,
]

BULK_VALUE_MAP: dict[str, dict[str, tuple[str, str]]] = {
    "flower_primary_color": {
        "white": ("flower_color", "white"),
        "green_brown_inconspicuous": ("flower_color", "green/brown/inconspicuous"),
        "yellow_orange": ("flower_color", "yellow/orange"),
        "red_pink": ("flower_color", "red/pink"),
        "blue_purple": ("flower_color", "blue/purple"),
        "multicolored_variable": ("flower_color", "multicolored/variable"),
    },
    "floral_symmetry": {
        "actinomorphic": ("flower_shape", "actinomorphic / radially symmetric"),
        "zygomorphic": ("flower_shape", "zygomorphic / bilaterally symmetric"),
    },
    "floral_form": {
        "tubular": ("flower_shape", "tubular"),
        "bell_campanulate": ("flower_shape", "bell-shaped"),
        "funnel_trumpet": ("flower_shape", "funnel / trumpet-shaped"),
        "open_radial": ("flower_shape", "open / bowl / cup-shaped"),
        "composite_head": ("flower_shape", "composite head / capitulum"),
        "urn_urceolate": ("flower_shape", "urn-shaped"),
    },
    "pollination_functional_guild": {
        "other_bees": ("pollination_guild", "bees"),
        "bumblebees": ("pollination_guild", "bumblebees"),
        "flies": ("pollination_guild", "flies"),
        "butterflies": ("pollination_guild", "butterflies"),
        "moths": ("pollination_guild", "moths"),
        "birds": ("pollination_guild", "birds"),
        "mixed_or_generalist": ("pollination_guild", "mixed"),
    },
    "pollination_vector_reported": {
        "bee": ("pollination_guild", "bees"),
        "bird": ("pollination_guild", "birds"),
        "butterfly_moth": ("pollination_guild", "mixed"),
        "fly": ("pollination_guild", "flies"),
        "wind": ("pollination_guild", "wind"),
    },
    "pollen_vector_mode": {
        "abiotic_wind": ("pollination_guild", "wind"),
    },
    "mating_system": {
        "mixed_mating": ("mating_system", "mixed_mating"),
        "predominantly_selfing": ("mating_system", "mainly_selfing"),
        "obligate_selfing": ("mating_system", "obligate_selfing"),
    },
    "mating_system_reported": {
        "mixed_mating": ("mating_system", "mixed_mating"),
    },
    "self_incompatibility": {
        "SI": ("self_incompatibility", "SI"),
        "SC": ("self_incompatibility", "SC"),
    },
    "self_incompatibility_reported": {
        "self_incompatible": ("self_incompatibility", "SI"),
    },
    "self_compatibility_reported": {
        "self_compatible": ("self_incompatibility", "SC"),
    },
}

FLORAL_CONTEXT = re.compile(
    r"\b(?:flower|flowers|floral|floret|florets|corolla|corollas|petal|petals|sepal|"
    r"sepals|tepal|tepals|perianth|inflorescence|inflorescences|capitulum|capitula)\b",
    re.IGNORECASE,
)
NON_FLORAL_CONTEXT = re.compile(
    r"\b(?:leaf|leaves|foliage|phyllode|phyllodes|fruit|fruits|berry|berries|capsule|"
    r"capsules|seed|seeds|bark|stem|wood|spine|spines|anther|anthers|pollen|style|"
    r"styles|stigma|stigmas|ovary|ovaries|ripe|ripening|powdery bloom)\b",
    re.IGNORECASE,
)
CULTIVAR_CONTEXT = re.compile(
    r"\b(?:cultivar|cultivars|variety|varieties|cv\.)\b",
    re.IGNORECASE,
)
NAME_CONTEXT = re.compile(
    r"\b(?:called|known as|commonly known|named|name|meaning)\b",
    re.IGNORECASE,
)
NEGATION = re.compile(r"\b(?:not|no|without|lack(?:s|ed|ing)?)\b", re.IGNORECASE)


@app.callback()
def main() -> None:
    """Build source-backed v1 trait rows without an LLM or paid API."""


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def _bounded_trait_text(text: str, species: str, max_chars: int = 4000) -> str:
    normalized = _text(text)
    if len(normalized) <= max_chars:
        return normalized
    anchors = [match.start() for match in TRAIT_TEXT_MARKERS.finditer(normalized)]
    species_key = _text(species)
    if species_key:
        species_match = re.search(re.escape(species_key), normalized, re.IGNORECASE)
        if species_match:
            anchors.append(species_match.start())
    windows: list[tuple[int, int]] = []
    for anchor in sorted(set(anchors)):
        start = max(0, anchor - 350)
        end = min(len(normalized), anchor + 650)
        if windows and start <= windows[-1][1]:
            windows[-1] = (windows[-1][0], max(windows[-1][1], end))
        else:
            windows.append((start, end))
    pieces: list[str] = []
    used = 0
    for start, end in windows:
        remaining = max_chars - used
        if remaining <= 0:
            break
        piece = normalized[start:end].strip()[:remaining]
        if piece:
            pieces.append(piece)
            used += len(piece)
    compact = " … ".join(pieces)
    return compact[:max_chars] if compact else normalized[:max_chars]


def compact_source_texts(sources: pd.DataFrame) -> pd.DataFrame:
    """Keep bounded evidence-oriented text plus the original retrieval hash."""
    compact = pd.DataFrame(sources, columns=SOURCE_COLUMNS).fillna("").copy()
    for index, row in compact.iterrows():
        citation = _text(row["source_citation"])
        if "retrieved_content_sha256=" in citation:
            continue
        full_text = _text(row["source_text"])
        stored = _bounded_trait_text(full_text, _text(row["accepted_species"]))
        content_sha = hashlib.sha256(full_text.encode("utf-8")).hexdigest()
        compact.loc[index, "source_text"] = stored
        compact.loc[index, "source_citation"] = (
            f"{citation} | retrieved_content_sha256={content_sha} | "
            f"stored_excerpt_chars={len(stored)}"
        ).strip(" |")
    return compact


def _excerpt(text: str, start: int, end: int, pad: int = 110) -> str:
    lo = max(0, start - pad)
    hi = min(len(text), end + pad)
    return _text(text[lo:hi])


def _local_context(text: str, start: int, end: int, pad: int = 90) -> str:
    return text[max(0, start - pad) : min(len(text), end + pad)]


def _has_floral_context(text: str, start: int, end: int) -> bool:
    sentence_start = max(
        text.rfind(".", 0, start),
        text.rfind(";", 0, start),
        text.rfind("\n", 0, start),
    )
    sentence_end_candidates = [
        boundary
        for boundary in (
            text.find(".", end),
            text.find(";", end),
            text.find("\n", end),
        )
        if boundary >= 0
    ]
    sentence_end = min(sentence_end_candidates) if sentence_end_candidates else len(text)
    local = text[sentence_start + 1 : sentence_end]
    floral = list(FLORAL_CONTEXT.finditer(local))
    if not floral:
        return False
    center = start - sentence_start - 1
    floral_distance = min(abs(center - match.start()) for match in floral)
    names = list(NAME_CONTEXT.finditer(local))
    if names and min(abs(center - match.start()) for match in names) <= floral_distance:
        return False
    non_floral = list(NON_FLORAL_CONTEXT.finditer(local))
    if not non_floral:
        return True
    return floral_distance <= min(abs(center - match.start()) for match in non_floral)


def _has_local_negation(text: str, start: int) -> bool:
    return bool(NEGATION.search(text[max(0, start - 45) : start]))


def _has_cultivar_context(text: str, start: int) -> bool:
    return bool(CULTIVAR_CONTEXT.search(text[max(0, start - 100) : start]))


def _source_evidence(source_type: str, citation: str) -> tuple[str, str, str]:
    source = source_type.casefold()
    citation_l = citation.casefold()
    if source == "gbif_species_description":
        return "flora", "medium", "A_flora"
    if source == "flora_or_monograph":
        return "flora", "medium", "A_flora"
    if source == "horticulture_site":
        return "horticulture", "medium", "C_horticulture"
    if source == "openalex_title_abstract":
        if re.search(r"\b(?:review|meta-analysis|systematic review)\b", citation_l):
            return "review", "medium", "A_literature"
        if re.search(
            r"\b(?:field study|field experiment|experimental study|hand pollination|"
            r"controlled pollination)\b",
            citation_l,
        ):
            return "field_study", "medium", "A_literature"
        return "inference", "low", "D_inference"
    if source in {
        "europe_pmc_full_text_xml",
        "europe_pmc_title_abstract_direct",
    }:
        return "field_study", "medium", "A_literature"
    if source == "gbif_checklist_description":
        return "flora", "medium", "A_flora"
    if source == "general_web_scientific":
        return "field_study", "medium", "A_literature"
    if source in {"general_web_institutional", "general_web_curated"}:
        return "flora", "medium", "A_flora"
    if source in {"wikipedia_extract", "wikidata_description"}:
        return "inference", "low", "B_encyclopedia"
    return "inference", "low", "D_inference"


def extract_evidence_from_sources(sources: pd.DataFrame) -> pd.DataFrame:
    """Extract conservative v1 field candidates from source-backed text rows."""
    missing = set(SOURCE_COLUMNS).difference(sources.columns)
    if missing:
        raise typer.BadParameter(f"source table is missing columns: {sorted(missing)}")
    rows: list[dict[str, str]] = []
    seen: set[tuple[str, str, str, str, str]] = set()
    for record in sources.fillna("").to_dict("records"):
        species = _text(record.get("accepted_species"))
        text = _text(record.get("source_text"))
        scope = _text(record.get("evidence_scope"))
        if not species or not text or scope != "species_direct":
            continue
        source_type = _text(record.get("source_type")) or "source_text"
        citation = _text(record.get("source_citation"))
        evidence_type, confidence, source_tier = _source_evidence(source_type, citation)
        for field, value, pattern, rule_id in ALL_RULES:
            for match in re.finditer(pattern, text, flags=re.IGNORECASE):
                if field in {"flower_color", "flower_shape"} and not _has_floral_context(
                    text, match.start(), match.end()
                ):
                    continue
                if field == "flower_color" and _has_cultivar_context(text, match.start()):
                    continue
                if _has_local_negation(text, match.start()):
                    continue
                key = (
                    species,
                    field,
                    value,
                    _text(record.get("source_url")),
                    rule_id,
                )
                if key in seen:
                    break
                seen.add(key)
                candidate_evidence_type = evidence_type
                candidate_confidence = confidence
                if rule_id.startswith("guild_attracts_"):
                    candidate_evidence_type = "inference"
                    candidate_confidence = "low"
                rows.append(
                    {
                        "species": species,
                        "field": field,
                        "value": value,
                        "matched_term": match.group(0),
                        "source_type": source_type,
                        "source_url": _text(record.get("source_url")),
                        "source_citation": citation,
                        "source_excerpt": _excerpt(text, match.start(), match.end()),
                        "evidence_type": candidate_evidence_type,
                        "confidence": candidate_confidence,
                        "source_tier": source_tier,
                        "rule_id": rule_id,
                    }
                )
                break
    return pd.DataFrame(rows, columns=EVIDENCE_COLUMNS)


def evidence_from_bulk_candidates(candidates: pd.DataFrame) -> pd.DataFrame:
    """Translate existing EOL/campaign candidates into conservative v1 fields."""
    required = {"accepted_species", "trait_name", "candidate_value"}
    missing = required.difference(candidates.columns)
    if missing:
        raise typer.BadParameter(f"candidate table is missing columns: {sorted(missing)}")
    rows: list[dict[str, str]] = []
    for record in candidates.fillna("").to_dict("records"):
        species = _text(record.get("accepted_species"))
        trait = _text(record.get("trait_name"))
        value = _text(record.get("candidate_value"))
        mapped = BULK_VALUE_MAP.get(trait, {}).get(value)
        if not species or mapped is None:
            continue
        field, v1_value = mapped
        excerpt = _text(
            record.get("source_excerpt")
            or record.get("evidence_excerpt")
            or record.get("raw_description")
        )
        rows.append(
            {
                "species": species,
                "field": field,
                "value": v1_value,
                "matched_term": _text(record.get("matched_term")) or value,
                "source_type": _text(record.get("source_type")) or "bulk_trait_candidate",
                "source_url": _text(record.get("source_url")),
                "source_citation": _text(record.get("source_citation")),
                "source_excerpt": excerpt or f"{trait} = {value}",
                "evidence_type": "inference",
                "confidence": "low",
                "source_tier": "B_database",
                "rule_id": f"bulk:{trait}:{value}",
            }
        )
    return pd.DataFrame(rows, columns=EVIDENCE_COLUMNS).drop_duplicates()


def v2_reported_candidates(evidence: pd.DataFrame) -> pd.DataFrame:
    """Translate explicit extracted evidence to the canonical v2 long format."""
    rows: list[dict[str, str]] = []
    for record in evidence.fillna("").to_dict("records"):
        rule_id = _text(record.get("rule_id"))
        if rule_id.startswith("guild_attracts_") or rule_id.startswith("visitation_"):
            continue
        mapped = V2_FIELD_VALUE_MAP.get((_text(record.get("field")), _text(record.get("value"))))
        if mapped is None:
            continue
        trait_layer, trait_name, value = mapped
        rows.append(
            {
                "accepted_species": _text(record.get("species")),
                "trait_layer": trait_layer,
                "trait_name": trait_name,
                "provisional_candidate_value": value,
                "matched_term": _text(record.get("matched_term")),
                "candidate_class": "reported",
                "source": _text(record.get("source_type")),
                "source_type": _text(record.get("source_type")),
                "source_url": _text(record.get("source_url")),
                "raw_description": _text(record.get("source_excerpt")),
                "evidence_scope": "species_direct",
                "wild_or_cultivated": "unknown",
                "review_status": "unreviewed",
            }
        )
    return pd.DataFrame(rows, columns=V2_REPORTED_COLUMNS).drop_duplicates()


def _unique(values: Iterable[object]) -> list[str]:
    return list(dict.fromkeys(_text(value) for value in values if _text(value)))


def _collapse_field(rows: pd.DataFrame, field: str) -> str:
    values = _unique(rows.loc[rows["field"].eq(field), "value"])
    if not values:
        return "unknown"
    if field == "pollination_guild":
        reduced = set(values)
        if reduced == {"bees", "bumblebees"}:
            return "bumblebees"
        return values[0] if len(values) == 1 else "mixed"
    if field == "mating_system":
        return values[0] if len(values) == 1 else "mixed_mating"
    if field == "self_incompatibility":
        explicit = [value for value in values if value in {"SI", "SC"}]
        if len(set(explicit)) > 1:
            return "unknown"
        if explicit:
            return explicit[0]
        return values[0]
    return ", ".join(values)


def _best_tier_rows(rows: pd.DataFrame, field: str) -> pd.DataFrame:
    selected = rows.loc[rows["field"].eq(field)]
    if selected.empty:
        return selected
    ranks = selected["source_tier"].map(SOURCE_TIER_ORDER).fillna(99)
    return selected.loc[ranks.eq(ranks.min())]


def collapse_v1_rows(species: list[str], evidence: pd.DataFrame) -> pd.DataFrame:
    """Collapse long evidence into the exact legacy nine-column contract."""
    output: list[dict[str, str]] = []
    empty = evidence.iloc[0:0]
    by_species = (
        {_text(name): rows for name, rows in evidence.groupby("species", sort=False)}
        if len(evidence)
        else {}
    )
    for name in species:
        rows = by_species.get(name, empty)
        if rows.empty:
            output.append(
                {
                    "species": name,
                    "flower_color": "unknown",
                    "flower_shape": "unknown",
                    "pollination_guild": "unknown",
                    "pollination_notes": "unknown",
                    "mating_system": "unknown",
                    "self_incompatibility": "unknown",
                    "evidence_type": "inference",
                    "confidence": "low",
                }
            )
            continue
        fields = (
            "flower_color",
            "flower_shape",
            "pollination_guild",
            "mating_system",
            "self_incompatibility",
        )
        selected_by_field = {field: _best_tier_rows(rows, field) for field in fields}
        values = {field: _collapse_field(selected_by_field[field], field) for field in fields}
        selected_rows = (
            pd.concat(
                [selected for selected in selected_by_field.values() if not selected.empty],
                ignore_index=True,
            )
            if any(not selected.empty for selected in selected_by_field.values())
            else empty
        )
        known = [field for field, value in values.items() if value != "unknown"]
        evidence_types = _unique(selected_rows["evidence_type"])
        notes = selected_rows.loc[
            selected_rows["field"].isin(
                {"pollination_guild", "mating_system", "self_incompatibility"}
            ),
            "source_excerpt",
        ]
        note_values = _unique(notes)[:3]
        confidences = set(selected_rows["confidence"])
        confidence = "low"
        if known and "medium" in confidences:
            confidence = "high" if len(known) >= 3 and confidences == {"medium"} else "medium"
        output.append(
            {
                "species": name,
                "flower_color": values["flower_color"],
                "flower_shape": values["flower_shape"],
                "pollination_guild": values["pollination_guild"],
                "pollination_notes": " | ".join(note_values) if note_values else "unknown",
                "mating_system": values["mating_system"],
                "self_incompatibility": values["self_incompatibility"],
                "evidence_type": "|".join(evidence_types) if evidence_types else "inference",
                "confidence": confidence,
            }
        )
    return validate_result_table(
        pd.DataFrame(output, columns=OUTPUT_COLUMNS), expected_species=species
    )


def infer_likely_traits(raw: pd.DataFrame, evidence: pd.DataFrame) -> pd.DataFrame:
    """Build auditable proxy traits without changing the direct nine-column table."""
    rows: list[dict[str, str]] = []
    seen: set[tuple[str, str, str, str]] = set()
    empty = evidence.iloc[0:0]
    by_species = (
        {_text(name): group for name, group in evidence.groupby("species", sort=False)}
        if len(evidence)
        else {}
    )

    def add(
        species: str,
        field: str,
        value: str,
        basis_fields: str,
        basis_values: str,
        source_urls: Iterable[object],
        rule_id: str,
    ) -> None:
        key = (species, field, value, rule_id)
        if key in seen:
            return
        seen.add(key)
        rows.append(
            {
                "species": species,
                "field": field,
                "value": value,
                "basis_fields": basis_fields,
                "basis_values": basis_values,
                "basis_source_urls": "|".join(_unique(source_urls)),
                "rule_id": rule_id,
                "confidence": "low",
            }
        )

    for record in raw.to_dict("records"):
        species = _text(record["species"])
        species_evidence = by_species.get(species, empty)
        direct_guild = _text(record["pollination_guild"])
        direct_mating = _text(record["mating_system"])
        direct_si = _text(record["self_incompatibility"])
        color = _text(record["flower_color"]).casefold()
        shape = _text(record["flower_shape"]).casefold()

        if direct_guild == "unknown":
            visitation = _best_tier_rows(species_evidence, "pollinator_visitation")
            for visitor in _unique(visitation["value"]):
                add(
                    species,
                    "pollination_guild",
                    f"likely_{visitor}",
                    "pollinator_visitation",
                    visitor,
                    visitation.loc[visitation["value"].eq(visitor), "source_url"],
                    f"likely_from_visitation_{visitor}",
                )

            floral_urls = species_evidence.loc[
                species_evidence["field"].isin({"flower_color", "flower_shape"}),
                "source_url",
            ]
            composite = any(term in shape for term in ("composite head", "capitulum"))
            tube = not composite and any(term in shape for term in ("tubular", "funnel", "trumpet"))
            warm = any(term in color for term in ("red", "pink", "orange"))
            bee_color = any(term in color for term in ("blue", "purple", "violet", "yellow"))
            bee_form = any(
                term in shape for term in ("zygomorphic", "bilateral", "papilionaceous", "spurred")
            )
            bee_form = bee_form or (not composite and "tubular" in shape)
            if warm and tube:
                add(
                    species,
                    "pollination_guild",
                    "likely_birds_or_butterflies",
                    "flower_color|flower_shape",
                    f"{record['flower_color']}|{record['flower_shape']}",
                    floral_urls,
                    "likely_guild_warm_tubular",
                )
            if bee_color and bee_form:
                add(
                    species,
                    "pollination_guild",
                    "likely_bees_or_bumblebees",
                    "flower_color|flower_shape",
                    f"{record['flower_color']}|{record['flower_shape']}",
                    floral_urls,
                    "likely_guild_bee_syndrome",
                )

        direct_selfing = (
            direct_mating
            in {
                "mixed_mating",
                "mainly_selfing",
                "obligate_selfing",
            }
            or direct_guild == "self"
        )
        selfing_signals = _best_tier_rows(species_evidence, "selfing_signal")
        if not direct_selfing and not selfing_signals.empty:
            signal_values = _unique(selfing_signals["value"])
            add(
                species,
                "selfing_status",
                "likely_selfing_capable",
                "selfing_signal",
                "|".join(signal_values),
                selfing_signals["source_url"],
                "likely_selfing_from_autogamy_or_cleistogamy",
            )
        if direct_si in {"unknown", "likely_SC"} and not selfing_signals.empty:
            add(
                species,
                "self_incompatibility",
                "likely_SC",
                "selfing_signal",
                "|".join(_unique(selfing_signals["value"])),
                selfing_signals["source_url"],
                "likely_sc_from_autogamy_or_cleistogamy",
            )

    return pd.DataFrame(rows, columns=LIKELY_EVIDENCE_COLUMNS)


def collapse_likely_traits(species: list[str], likely: pd.DataFrame) -> pd.DataFrame:
    """Collapse proxy evidence while retaining the explicit likely namespace."""
    output: list[dict[str, str]] = []
    empty = likely.iloc[0:0]
    by_species = (
        {_text(name): group for name, group in likely.groupby("species", sort=False)}
        if len(likely)
        else {}
    )
    field_columns = {
        "pollination_guild": "likely_pollination_guild",
        "selfing_status": "likely_selfing_status",
        "self_incompatibility": "likely_self_incompatibility",
    }
    for name in species:
        rows = by_species.get(name, empty)
        values = {
            column: "|".join(_unique(rows.loc[rows["field"].eq(field), "value"])) or "unknown"
            for field, column in field_columns.items()
        }
        output.append(
            {
                "species": name,
                **values,
                "likely_basis": " || ".join(_unique(rows["basis_values"])[:4]) or "unknown",
                "likely_rule_ids": "|".join(_unique(rows["rule_id"])) or "unknown",
                "likely_confidence": "low",
            }
        )
    return pd.DataFrame(output, columns=LIKELY_WIDE_COLUMNS)


def build_priority_traits(raw: pd.DataFrame, likely_wide: pd.DataFrame) -> pd.DataFrame:
    """Resolve direct-first core traits for coverage and downstream sensitivity runs."""
    likely_by_species = likely_wide.set_index("species").to_dict("index")
    output: list[dict[str, str]] = []
    for record in raw.to_dict("records"):
        species = _text(record["species"])
        inferred = likely_by_species[species]
        direct_guild = _text(record["pollination_guild"])
        likely_guild = _text(inferred["likely_pollination_guild"])
        guild = direct_guild if direct_guild != "unknown" else likely_guild
        guild_mode = (
            "direct"
            if direct_guild != "unknown"
            else ("likely" if likely_guild != "unknown" else "unknown")
        )

        mating = _text(record["mating_system"])
        direct_selfing = {
            "obligate_outcrossing": "outcrossing",
            "mixed_mating": "mixed_selfing",
            "mainly_selfing": "selfing",
            "obligate_selfing": "selfing",
        }.get(mating, "selfing" if direct_guild == "self" else "unknown")
        likely_selfing = _text(inferred["likely_selfing_status"])
        selfing = direct_selfing if direct_selfing != "unknown" else likely_selfing
        selfing_mode = (
            "direct"
            if direct_selfing != "unknown"
            else ("likely" if likely_selfing != "unknown" else "unknown")
        )

        direct_si = _text(record["self_incompatibility"])
        likely_si = _text(inferred["likely_self_incompatibility"])
        si = direct_si if direct_si != "unknown" else likely_si
        si_mode = (
            "likely"
            if direct_si.startswith("likely_")
            else "direct"
            if direct_si != "unknown"
            else ("likely" if likely_si != "unknown" else "unknown")
        )
        output.append(
            {
                "species": species,
                "flower_color": _text(record["flower_color"]),
                "flower_shape": _text(record["flower_shape"]),
                "pollination_guild": guild,
                "pollination_guild_mode": guild_mode,
                "mating_system": mating,
                "selfing_status": selfing,
                "selfing_status_mode": selfing_mode,
                "self_incompatibility": si,
                "self_incompatibility_mode": si_mode,
                "basis_rule_ids": _text(inferred["likely_rule_ids"]),
            }
        )
    return pd.DataFrame(output, columns=PRIORITY_COLUMNS)


def _json_getter(
    client: Any,
    *,
    max_attempts: int = 3,
) -> Callable[[str, dict[str, Any]], dict[str, Any]]:
    def getter(url: str, params: dict[str, Any]) -> dict[str, Any]:
        for attempt in range(max_attempts):
            response = client.get(url, params=params)
            status = int(getattr(response, "status_code", 0))
            retryable = status in {408, 425, 429, 500, 502, 503, 504}
            if not retryable or attempt + 1 >= max_attempts:
                response.raise_for_status()
                return response.json()
            retry_after = _text(getattr(response, "headers", {}).get("Retry-After"))
            delay = float(retry_after) if retry_after.replace(".", "", 1).isdigit() else 2**attempt
            time.sleep(min(delay, 10.0))
        raise RuntimeError("unreachable JSON retry state")

    return getter


def _wikipedia_exact_sources(
    species: list[str],
    getter: Callable[[str, dict[str, Any]], dict[str, Any]],
    pause_seconds: float,
    wikipedia_languages: list[str] | None = None,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Fetch exact scientific-name pages without the extra Wikidata round trips."""
    rows: list[dict[str, str]] = []
    errors: list[dict[str, str]] = []
    languages = trait_web_reported_scout.normalize_wikipedia_languages(wikipedia_languages)
    for name in species:
        for language in languages:
            try:
                extract, url, resolved = trait_web_reported_scout.fetch_wikipedia_extract(
                    getter,
                    name,
                    language,
                )
            except Exception as exc:  # noqa: BLE001
                errors.append(
                    {
                        "species": name,
                        "source": f"wikipedia_{language}",
                        "error": str(exc),
                    }
                )
                continue
            if not extract:
                continue
            scope = (
                "species_direct"
                if name.casefold() in extract.casefold()
                or trait_web_reported_scout._scope_for(name, resolved) == "species_direct"  # noqa: SLF001
                else "species_indirect"
            )
            rows.append(
                {
                    "accepted_species": name,
                    "source_text": _text(extract),
                    "source_url": url,
                    "source_citation": f"{resolved or name} Wikipedia extract ({language})",
                    "source_type": "wikipedia_extract",
                    "evidence_scope": scope,
                }
            )
            break
        if pause_seconds > 0:
            time.sleep(pause_seconds)
    return rows, errors


class _VisibleTextParser(HTMLParser):
    """Extract visible HTML text while preserving block boundaries."""

    _SKIP_TAGS = {
        "aside",
        "footer",
        "form",
        "header",
        "nav",
        "noscript",
        "script",
        "style",
        "svg",
        "template",
    }
    _BLOCK_TAGS = {
        "article",
        "br",
        "div",
        "h1",
        "h2",
        "h3",
        "h4",
        "h5",
        "li",
        "p",
        "section",
        "td",
        "th",
        "tr",
    }

    def __init__(self, *, include_form_text: bool = False) -> None:
        super().__init__(convert_charrefs=True)
        self._skip_tags = self._SKIP_TAGS - ({"form"} if include_form_text else set())
        self._skip_depth = 0
        self._parts: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        del attrs
        if tag in self._skip_tags:
            self._skip_depth += 1
        elif tag in self._BLOCK_TAGS and not self._skip_depth:
            self._parts.append(". ")

    def handle_endtag(self, tag: str) -> None:
        if tag in self._skip_tags and self._skip_depth:
            self._skip_depth -= 1
        elif tag in self._BLOCK_TAGS and not self._skip_depth:
            self._parts.append(". ")

    def handle_data(self, data: str) -> None:
        if not self._skip_depth and data.strip():
            self._parts.append(data)
            self._parts.append(" ")

    def text(self) -> str:
        text = "".join(self._parts)
        text = re.sub(r"(?:\s*\.\s*){2,}", ". ", text)
        return _text(text)


def html_visible_text(html: str, *, include_form_text: bool = False) -> str:
    parser = _VisibleTextParser(include_form_text=include_form_text)
    parser.feed(html)
    parser.close()
    return parser.text()


class _WorldFloraTaxonLinkParser(HTMLParser):
    def __init__(self, species: str) -> None:
        super().__init__(convert_charrefs=True)
        self.species = species.casefold()
        self.href = ""

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag != "a" or self.href:
            return
        attributes = {key: _text(value) for key, value in attrs}
        title = attributes.get("title", "").casefold()
        href = attributes.get("href", "")
        if title == self.species and href.startswith("/taxon/"):
            self.href = href.split(";", 1)[0]


def _world_flora_sources(
    species: list[str],
    client: Any,
    pause_seconds: float,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    rows: list[dict[str, str]] = []
    errors: list[dict[str, str]] = []
    base_url = "https://www.worldfloraonline.org"
    for name in species:
        try:
            search_url = f"{base_url}/search?query={quote_plus(name)}"
            response = client.get(search_url)
            response.raise_for_status()
            parser = _WorldFloraTaxonLinkParser(name)
            parser.feed(response.text)
            parser.close()
            if not parser.href:
                continue
            taxon_url = urljoin(base_url, parser.href)
            response = client.get(taxon_url)
            response.raise_for_status()
            if len(response.content) > 3_000_000:
                raise ValueError("World Flora Online response exceeded 3 MB")
            text = html_visible_text(response.text)
            if name.casefold() not in text.casefold() or len(text) < 300:
                continue
            rows.append(
                {
                    "accepted_species": name,
                    "source_text": text,
                    "source_url": str(response.url),
                    "source_citation": f"{name} - World Flora Online",
                    "source_type": "flora_or_monograph",
                    "evidence_scope": "species_direct",
                }
            )
        except Exception as exc:  # noqa: BLE001
            errors.append({"species": name, "source": "world_flora_online", "error": str(exc)})
        if pause_seconds > 0:
            time.sleep(pause_seconds)
    return rows, errors


def _web_source_url(source: dict[str, str], species: str) -> str:
    if source["url_style"] == "slug":
        name = re.sub(r"[^a-z0-9]+", "-", species.casefold()).strip("-")
        encoded = quote(name, safe="-")
    else:
        encoded = quote_plus(species)
    return source["url_template"].format(name=encoded)


def _sitemap_locations(content: bytes) -> tuple[str, list[str]]:
    if len(content) > 8_000_000:
        raise ValueError("sitemap response exceeded 8 MB")
    root = ElementTree.fromstring(content)
    root_type = root.tag.rsplit("}", 1)[-1]
    locations = [
        _text(element.text)
        for element in root.iter()
        if element.tag.rsplit("}", 1)[-1] == "loc" and _text(element.text)
    ]
    return root_type, locations


def _fetch_sitemap_urls(client: Any, sitemap_url: str) -> set[str]:
    """Fetch one sitemap or one bounded sitemap-index level."""
    response = client.get(sitemap_url)
    response.raise_for_status()
    root_type, locations = _sitemap_locations(response.content)
    if root_type != "sitemapindex":
        return {location.casefold() for location in locations}
    if len(locations) > 50:
        raise ValueError("sitemap index exceeded 50 child sitemaps")
    indexed_urls: set[str] = set()
    for child_url in locations:
        child = client.get(child_url)
        child.raise_for_status()
        child_type, child_locations = _sitemap_locations(child.content)
        if child_type == "sitemapindex":
            raise ValueError("nested sitemap index depth exceeded")
        indexed_urls.update(location.casefold() for location in child_locations)
        if len(indexed_urls) > 250_000:
            raise ValueError("sitemap index exceeded 250,000 URLs")
    return indexed_urls


def build_web_description_index(
    client: Any,
    *,
    max_attempts: int = 3,
    backoff_seconds: float = 1.0,
) -> dict[str, Any]:
    """Fetch configured plant-site indexes once for reuse by all shards."""
    if max_attempts < 1:
        raise ValueError("max_attempts must be positive")
    sources: dict[str, dict[str, Any]] = {}
    for source in WEB_DESCRIPTION_SOURCES:
        sitemap_url = source.get("sitemap_url", "")
        if not sitemap_url:
            continue
        error = ""
        urls: list[str] = []
        attempts = 0
        for attempts in range(1, max_attempts + 1):
            try:
                urls = sorted(_fetch_sitemap_urls(client, sitemap_url))
                error = ""
                break
            except Exception as exc:  # noqa: BLE001
                error = f"{type(exc).__name__}:{exc}"
                if attempts < max_attempts and backoff_seconds > 0:
                    time.sleep(backoff_seconds * (2 ** (attempts - 1)))
        if error:
            sources[source["name"]] = {
                "status": "retry",
                "sitemap_url": sitemap_url,
                "n_urls": 0,
                "urls": [],
                "error": error,
                "attempts": attempts,
            }
        else:
            sources[source["name"]] = {
                "status": "completed",
                "sitemap_url": sitemap_url,
                "n_urls": len(urls),
                "urls": urls,
                "error": "",
                "attempts": attempts,
            }
    return {
        "contract_version": WEB_DESCRIPTION_INDEX_VERSION,
        "sources": sources,
    }


def load_web_description_index(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("contract_version") != WEB_DESCRIPTION_INDEX_VERSION:
        raise typer.BadParameter(f"unsupported web-description index contract in {path}")
    if not isinstance(payload.get("sources"), dict):
        raise typer.BadParameter(f"web-description index has no sources mapping: {path}")
    return payload


def _web_description_sources(
    species: list[str],
    client: Any,
    pause_seconds: float,
    source_types: set[str] | None = None,
    source_index: dict[str, Any] | None = None,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    configured_sources = [
        source
        for source in WEB_DESCRIPTION_SOURCES
        if source_types is None or source["source_type"] in source_types
    ]
    if not configured_sources:
        return [], []

    def fetch_source(
        source: dict[str, str],
    ) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
        source_rows: list[dict[str, str]] = []
        source_errors: list[dict[str, str]] = []
        indexed_urls: set[str] | None = None
        sitemap_url = source.get("sitemap_url", "")
        if sitemap_url:
            try:
                if source_index is None:
                    indexed_urls = _fetch_sitemap_urls(client, sitemap_url)
                else:
                    entry = (source_index.get("sources") or {}).get(source["name"])
                    if not isinstance(entry, dict):
                        raise ValueError("source missing from frozen sitemap index")
                    if entry.get("status") != "completed":
                        raise ValueError(
                            f"frozen sitemap index unavailable:{entry.get('error', '')}"
                        )
                    indexed_urls = {
                        _text(url).casefold() for url in entry.get("urls", []) if _text(url)
                    }
            except Exception as exc:  # noqa: BLE001
                source_errors.append(
                    {
                        "species": "",
                        "source": source["name"],
                        "error": f"sitemap:{exc}",
                    }
                )
                return source_rows, source_errors
        for name in species:
            url = _web_source_url(source, name)
            if indexed_urls is not None and url.casefold() not in indexed_urls:
                continue
            try:
                response = client.get(url)
                if response.status_code == 404:
                    continue
                response.raise_for_status()
                if len(response.content) > 3_000_000:
                    raise ValueError("HTML response exceeded 3 MB")
                text = html_visible_text(
                    response.text,
                    include_form_text=source.get("visible_text_mode") == "aspnet_form",
                )
                if name.casefold() not in text.casefold():
                    continue
                required = source.get("required_marker", "")
                rejected = source.get("reject_marker", "")
                if required and required.casefold() not in text.casefold():
                    continue
                if rejected and rejected.casefold() in text.casefold():
                    continue
                source_rows.append(
                    {
                        "accepted_species": name,
                        "source_text": text,
                        "source_url": str(response.url),
                        "source_citation": f"{name} - {source['name']}",
                        "source_type": source["source_type"],
                        "evidence_scope": "species_direct",
                    }
                )
            except Exception as exc:  # noqa: BLE001
                source_errors.append(
                    {
                        "species": name,
                        "source": source["name"],
                        "error": str(exc),
                    }
                )
            if pause_seconds > 0:
                time.sleep(pause_seconds)
        return source_rows, source_errors

    rows: list[dict[str, str]] = []
    errors: list[dict[str, str]] = []
    # One sequential worker per site keeps per-host traffic bounded while allowing
    # independent sites to run in parallel.
    with ThreadPoolExecutor(max_workers=len(configured_sources)) as executor:
        results = executor.map(fetch_source, configured_sources)
        for source_rows, source_errors in results:
            rows.extend(source_rows)
            errors.extend(source_errors)
    return rows, errors


def _strict_binomial_identity(value: object) -> str:
    """Return a comparison key only for an unadorned two-part scientific name."""
    parts = _text(value).split()
    if len(parts) != 2:
        return ""
    for part in parts:
        letters = part.replace("-", "")
        if not letters or not letters.isalpha() or part.startswith("×"):
            return ""
    return " ".join(parts).casefold()


def _gbif_taxon_name(payload: dict[str, Any]) -> str:
    canonical_name = _text(payload.get("canonicalName"))
    if _strict_binomial_identity(canonical_name):
        return canonical_name
    scientific_name = _text(payload.get("scientificName"))
    candidate = " ".join(scientific_name.split()[:2])
    return candidate if _strict_binomial_identity(candidate) else ""


def _validated_gbif_species_match(
    getter: Callable[[str, dict[str, Any]], dict[str, Any]],
    name: str,
) -> dict[str, str] | None:
    """Accept only an exact species match or a verified synonym-to-species mapping."""
    requested_identity = _strict_binomial_identity(name)
    if not requested_identity:
        return None

    payload = getter(
        f"{trait_descriptive_scout.GBIF_API}/species/match",
        {"name": name},
    )
    match_type = _text(payload.get("matchType")).upper()
    rank = _text(payload.get("rank")).upper()
    status = _text(payload.get("status")).upper()
    matched_name = _gbif_taxon_name(payload)
    matched_identity = _strict_binomial_identity(matched_name)
    matched_key = _text(payload.get("usageKey") or payload.get("nubKey"))
    if (
        match_type != "EXACT"
        or rank != "SPECIES"
        or matched_identity != requested_identity
        or not matched_key
    ):
        return None

    resolved_key = matched_key
    resolved_name = matched_name
    if status == "SYNONYM":
        resolved_key = _text(payload.get("acceptedUsageKey"))
        if not resolved_key:
            return None
        accepted = getter(
            f"{trait_descriptive_scout.GBIF_API}/species/{resolved_key}",
            {},
        )
        accepted_key = _text(accepted.get("key") or accepted.get("nubKey"))
        accepted_rank = _text(accepted.get("rank")).upper()
        accepted_status = _text(accepted.get("taxonomicStatus") or accepted.get("status")).upper()
        resolved_name = _gbif_taxon_name(accepted)
        if (
            accepted_key != resolved_key
            or accepted_rank != "SPECIES"
            or accepted_status != "ACCEPTED"
            or not resolved_name
        ):
            return None
    elif status != "ACCEPTED":
        return None

    return {
        "query_name": _text(name),
        "matched_name": matched_name,
        "resolved_name": resolved_name,
        "matched_key": matched_key,
        "resolved_key": resolved_key,
        "match_type": match_type,
        "rank": rank,
        "status": status,
        "confidence": _text(payload.get("confidence")) or "unknown",
    }


def _gbif_match_provenance(match: dict[str, str]) -> str:
    return "; ".join(
        [
            "GBIF taxon match",
            f"query={match['query_name']}",
            f"matched={match['matched_name']}",
            f"resolved={match['resolved_name']}",
            f"matchType={match['match_type']}",
            f"rank={match['rank']}",
            f"status={match['status']}",
            f"confidence={match['confidence']}",
            f"matchedKey={match['matched_key']}",
            f"resolvedKey={match['resolved_key']}",
        ]
    )


def _gbif_sources(
    species: list[str],
    getter: Callable[[str, dict[str, Any]], dict[str, Any]],
    max_descriptions: int,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    rows: list[dict[str, str]] = []
    errors: list[dict[str, str]] = []
    wanted_types = {"morphology", "description", "general", "diagnostic", "reproduction"}
    for name in species:
        try:
            match = _validated_gbif_species_match(getter, name)
            key = match["resolved_key"] if match else ""
            descriptions = (
                trait_descriptive_scout.fetch_descriptions(getter, key, max_descriptions)
                if key
                else []
            )
            for document in descriptions:
                description_type = _text(document.get("type")).casefold()
                text = _text(document.get("description"))
                if not text or (description_type and description_type not in wanted_types):
                    continue
                citation = _text(document.get("source"))
                provenance = _gbif_match_provenance(match) if match else ""
                rows.append(
                    {
                        "accepted_species": name,
                        "source_text": text,
                        "source_url": f"https://www.gbif.org/species/{key}",
                        "source_citation": " | ".join(
                            part for part in (citation, provenance) if part
                        ),
                        "source_type": "gbif_species_description",
                        "evidence_scope": "species_direct",
                    }
                )
        except Exception as exc:  # noqa: BLE001
            errors.append({"species": name, "source": "gbif", "error": str(exc)})
    return rows, errors


def _missing_species_for_fields(
    species: list[str],
    source_rows: list[dict[str, str]],
    seed_evidence: pd.DataFrame,
    fields: set[str],
) -> list[str]:
    source_table = pd.DataFrame(source_rows, columns=SOURCE_COLUMNS)
    text_evidence = extract_evidence_from_sources(source_table)
    evidence = pd.concat([seed_evidence, text_evidence], ignore_index=True)
    covered = {
        (row.species, row.field) for row in evidence[["species", "field"]].itertuples(index=False)
    }
    return [name for name in species if any((name, field) not in covered for field in fields)]


def _europe_pmc_sources(
    species: list[str],
    pause_seconds: float,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Adapt exact Europe PMC title/abstract discovery to the shared source schema."""
    if not species:
        return [], []
    batch = discover_europe_pmc_sources(
        species,
        min_interval_seconds=max(1.0, pause_seconds),
    )
    full_text = discover_europe_pmc_full_text_sources(
        batch.sources,
        min_interval_seconds=max(0.25, pause_seconds),
    )
    rows = [
        {
            "accepted_species": record["accepted_species"],
            "source_text": record["source_text"],
            "source_url": record["source_url"],
            "source_citation": record["source_citation"],
            "source_type": "europe_pmc_title_abstract",
            "evidence_scope": "species_direct",
        }
        for record in batch.sources
    ]
    rows.extend(
        {
            "accepted_species": record["accepted_species"],
            "source_text": record["source_text"],
            "source_url": record["source_url"],
            "source_citation": record["source_citation"],
            "source_type": "europe_pmc_full_text_xml",
            "evidence_scope": "species_direct",
        }
        for record in full_text.sources
    )
    errors = [
        {
            "species": record["accepted_species"],
            "source": "europe_pmc",
            "error": ":".join(
                value
                for value in (
                    record["error_class"],
                    record["error_code"],
                    record["message"],
                )
                if value
            ),
        }
        for record in batch.errors
    ]
    errors.extend(
        {
            "species": record["accepted_species"],
            "source": "europe_pmc_full_text",
            "error": ":".join(
                value
                for value in (
                    record["error_class"],
                    record["error_code"],
                    record["message"],
                )
                if value
            ),
        }
        for record in full_text.errors
    )
    return rows, errors


def collect_free_sources(
    species: list[str],
    *,
    include_gbif: bool = True,
    include_wikimedia: bool = True,
    include_openalex: bool = False,
    include_europe_pmc: bool = False,
    include_web_descriptions: bool = True,
    include_world_flora: bool = False,
    max_descriptions: int = 20,
    openalex_per_page: int = 5,
    pause_seconds: float = 0.05,
    wikipedia_languages: list[str] | None = None,
    seed_evidence: pd.DataFrame | None = None,
    web_description_index_path: Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Retrieve public text in reliability order, querying only remaining gaps."""
    import ssl

    import httpx
    import truststore

    source_rows: list[dict[str, str]] = []
    error_rows: list[dict[str, str]] = []
    seed = seed_evidence if seed_evidence is not None else pd.DataFrame(columns=EVIDENCE_COLUMNS)
    all_fields = {
        "flower_color",
        "flower_shape",
        "mating_system",
        "self_incompatibility",
    }
    floral_fields = {"flower_color", "flower_shape"}
    reproductive_fields = {"mating_system", "self_incompatibility"}
    web_description_index = (
        load_web_description_index(web_description_index_path)
        if include_web_descriptions and web_description_index_path is not None
        else None
    )
    with httpx.Client(
        timeout=45.0,
        follow_redirects=True,
        verify=truststore.SSLContext(ssl.PROTOCOL_TLS_CLIENT),
        headers={
            "User-Agent": "island-floral-v2/0.1 (https://github.com/zuizui0223/island)",
        },
    ) as client:
        getter = _json_getter(client)
        if include_gbif:
            missing = _missing_species_for_fields(species, source_rows, seed, floral_fields)
            rows, errors = _gbif_sources(missing, getter, max_descriptions)
            source_rows.extend(rows)
            error_rows.extend(errors)
        if include_world_flora:
            missing = _missing_species_for_fields(species, source_rows, seed, all_fields)
            rows, errors = _world_flora_sources(missing, client, pause_seconds)
            source_rows.extend(rows)
            error_rows.extend(errors)
        if include_web_descriptions:
            missing = _missing_species_for_fields(species, source_rows, seed, all_fields)
            rows, errors = _web_description_sources(
                missing,
                client,
                pause_seconds,
                source_types={"flora_or_monograph"},
                source_index=web_description_index,
            )
            source_rows.extend(rows)
            error_rows.extend(errors)
        if include_wikimedia:
            missing = _missing_species_for_fields(species, source_rows, seed, all_fields)
            rows, errors = _wikipedia_exact_sources(
                missing,
                getter,
                pause_seconds,
                wikipedia_languages=wikipedia_languages,
            )
            source_rows.extend(rows)
            error_rows.extend(errors)
        if include_openalex:
            missing = _missing_species_for_fields(
                species,
                source_rows,
                seed,
                reproductive_fields,
            )
            species_frame = pd.DataFrame({"accepted_species": missing})
            try:
                rows = reported_ecology_machine.openalex_text_sources(
                    species_frame,
                    len(missing),
                    openalex_per_page,
                    pause_seconds,
                )
                source_rows.extend(rows.to_dict("records"))
            except Exception as exc:  # noqa: BLE001
                error_rows.append({"species": "", "source": "openalex", "error": str(exc)})
        if include_europe_pmc:
            missing = _missing_species_for_fields(
                species,
                source_rows,
                seed,
                reproductive_fields,
            )
            rows, errors = _europe_pmc_sources(missing, pause_seconds)
            source_rows.extend(rows)
            error_rows.extend(errors)
        if include_web_descriptions:
            missing = _missing_species_for_fields(
                species,
                source_rows,
                seed,
                floral_fields,
            )
            rows, errors = _web_description_sources(
                missing,
                client,
                pause_seconds,
                source_types={"horticulture_site"},
                source_index=web_description_index,
            )
            source_rows.extend(rows)
            error_rows.extend(errors)
    sources = pd.DataFrame(source_rows, columns=SOURCE_COLUMNS).drop_duplicates()
    errors = pd.DataFrame(error_rows, columns=["species", "source", "error"])
    return sources.reset_index(drop=True), errors.reset_index(drop=True)


def _species_list(table: pd.DataFrame, column: str, max_taxa: int) -> list[str]:
    if column not in table.columns:
        raise typer.BadParameter(f"species column {column!r} not found")
    return _unique(table[column])[:max_taxa]


def write_search_outputs(
    output_dir: Path,
    species: list[str],
    sources: pd.DataFrame,
    errors: pd.DataFrame,
    bulk_candidates: pd.DataFrame | None = None,
) -> dict[str, Any]:
    sources = compact_source_texts(sources)
    text_evidence = extract_evidence_from_sources(sources)
    bulk_evidence = (
        evidence_from_bulk_candidates(bulk_candidates)
        if bulk_candidates is not None
        else pd.DataFrame(columns=EVIDENCE_COLUMNS)
    )
    evidence = pd.concat([text_evidence, bulk_evidence], ignore_index=True).drop_duplicates()
    raw = collapse_v1_rows(species, evidence)
    derived = derive_v1_categories(raw)
    likely_evidence = infer_likely_traits(raw, evidence)
    likely_wide = collapse_likely_traits(species, likely_evidence)
    priority = build_priority_traits(raw, likely_wide)
    v2_candidates = v2_reported_candidates(evidence)
    output_dir.mkdir(parents=True, exist_ok=True)
    sources.to_csv(output_dir / "source_texts.csv", index=False)
    evidence.to_csv(output_dir / "v1_category_evidence.csv", index=False)
    raw.to_csv(output_dir / "v1_category_traits.csv", index=False)
    derived.to_csv(output_dir / "v1_category_traits_derived.csv", index=False)
    likely_evidence.to_csv(output_dir / "v1_likely_traits.csv", index=False)
    likely_wide.to_csv(output_dir / "v1_likely_traits_wide.csv", index=False)
    priority.to_csv(output_dir / "v1_priority_traits.csv", index=False)
    v2_candidates.to_csv(output_dir / "v2_reported_candidates.csv", index=False)
    errors.to_csv(output_dir / "lookup_errors.csv", index=False)

    fields = [
        "flower_color",
        "flower_shape",
        "pollination_guild",
        "mating_system",
        "self_incompatibility",
    ]
    coverage = {
        field: {
            "n_species": int(raw[field].ne("unknown").sum()),
            "rate": round(float(raw[field].ne("unknown").mean()), 4),
        }
        for field in fields
    }
    complete = raw[fields].ne("unknown").all(axis=1)
    any_known = raw[fields].ne("unknown").any(axis=1)
    core_masks = {
        "flower_color": priority["flower_color"].ne("unknown"),
        "flower_shape": priority["flower_shape"].ne("unknown"),
        "color_and_shape": priority[["flower_color", "flower_shape"]].ne("unknown").all(axis=1),
        "selfing_direct": priority["selfing_status_mode"].eq("direct"),
        "selfing_direct_or_likely": priority["selfing_status"].ne("unknown"),
        "pollination_direct_or_likely": priority["pollination_guild"].ne("unknown"),
        "self_incompatibility_direct_or_likely": priority["self_incompatibility"].ne("unknown"),
    }
    core_masks["color_shape_and_selfing"] = (
        core_masks["color_and_shape"] & core_masks["selfing_direct_or_likely"]
    )
    core_coverage = {
        name: {"n_species": int(mask.sum()), "rate": round(float(mask.mean()), 4)}
        for name, mask in core_masks.items()
    }
    v2_species_by_trait = {
        trait: set(rows["accepted_species"])
        for trait, rows in v2_candidates.groupby("trait_name", sort=False)
    }
    v2_color = v2_species_by_trait.get("flower_primary_color", set())
    v2_form = v2_species_by_trait.get("floral_form", set())
    v2_symmetry = v2_species_by_trait.get("floral_symmetry", set())
    v2_shape = v2_form | v2_symmetry
    v2_reproductive = set().union(
        *(
            v2_species_by_trait.get(trait, set())
            for trait in (
                "autonomous_selfing_capacity",
                "mating_system",
                "self_incompatibility",
            )
        )
    )
    v2_core_coverage = {
        "flower_color": len(v2_color),
        "flower_form_or_symmetry": len(v2_shape),
        "color_and_shape": len(v2_color & v2_shape),
        "any_reproductive_assurance": len(v2_reproductive),
    }
    report = {
        "version": "1.0",
        "method": "free source search plus deterministic explicit-text extraction",
        "n_species": len(species),
        "n_source_rows": int(len(sources)),
        "n_evidence_rows": int(len(evidence)),
        "n_text_evidence_rows": int(len(text_evidence)),
        "n_bulk_evidence_rows": int(len(bulk_evidence)),
        "n_likely_evidence_rows": int(len(likely_evidence)),
        "n_v2_reported_candidate_rows": int(len(v2_candidates)),
        "evidence_rows_by_source_tier": {
            str(tier): int(count)
            for tier, count in evidence["source_tier"].value_counts().sort_index().items()
        },
        "n_lookup_errors": int(len(errors)),
        "coverage": coverage,
        "core_coverage": core_coverage,
        "v2_core_coverage_n_species": v2_core_coverage,
        "n_species_any_trait": int(any_known.sum()),
        "any_trait_rate": round(float(any_known.mean()), 4),
        "n_species_all_five_traits": int(complete.sum()),
        "all_five_traits_rate": round(float(complete.mean()), 4),
        "policy": (
            "Only species-direct explicit text is collapsed into the legacy table. "
            "Proxy results use a separate likely_* namespace and never use island "
            "Bombus presence or absence. Unknown stays unknown; all evidence excerpts "
            "and URLs remain auditable."
        ),
    }
    (output_dir / "v1_category_search_manifest.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return report


@app.command("build-web-index")
def build_web_index_command(
    output_json: Path = typer.Option(...),
) -> None:
    """Freeze plant-site sitemap indexes once for reuse across all shards."""
    import ssl

    import httpx
    import truststore

    with httpx.Client(
        timeout=45.0,
        follow_redirects=True,
        verify=truststore.SSLContext(ssl.PROTOCOL_TLS_CLIENT),
        headers={
            "User-Agent": "island-floral-v2/0.1 (https://github.com/zuizui0223/island)",
        },
    ) as client:
        payload = build_web_description_index(client)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    failed = {
        name: entry.get("error", "")
        for name, entry in payload["sources"].items()
        if entry.get("status") != "completed"
    }
    if failed:
        raise typer.BadParameter(f"plant-site sitemap index failures: {failed}")
    typer.echo(
        json.dumps(
            {
                "contract_version": WEB_DESCRIPTION_INDEX_VERSION,
                "output_json": str(output_json),
                "n_sources": len(payload["sources"]),
                "n_urls": sum(entry["n_urls"] for entry in payload["sources"].values()),
            },
            ensure_ascii=False,
        )
    )


@app.command("search")
def search_command(
    species_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    species_column: str = typer.Option("accepted_species"),
    max_taxa: int = typer.Option(25, min=1, max=500),
    max_descriptions: int = typer.Option(20, min=1, max=100),
    openalex_per_page: int = typer.Option(5, min=1, max=25),
    pause_seconds: float = typer.Option(0.05, min=0.0),
    wikipedia_language: list[str] | None = typer.Option(None, "--wikipedia-language"),
    gbif: bool = typer.Option(True, "--gbif/--no-gbif"),
    wikimedia: bool = typer.Option(True, "--wikimedia/--no-wikimedia"),
    openalex: bool = typer.Option(False, "--openalex/--no-openalex"),
    europe_pmc: bool = typer.Option(False, "--europe-pmc/--no-europe-pmc"),
    web_descriptions: bool = typer.Option(
        True,
        "--web-descriptions/--no-web-descriptions",
    ),
    world_flora: bool = typer.Option(False, "--world-flora/--no-world-flora"),
    candidate_csv: Path | None = typer.Option(None, exists=True),
    web_description_index: Path | None = typer.Option(None, exists=True),
) -> None:
    """Search public sources and write v1 rows plus auditable evidence."""
    table = pd.read_csv(species_csv, dtype=str).fillna("")
    species = _species_list(table, species_column, max_taxa)
    if not species:
        raise typer.BadParameter("no species selected")
    bulk_candidates = (
        pd.read_csv(candidate_csv, dtype=str).fillna("") if candidate_csv is not None else None
    )
    seed_evidence = (
        evidence_from_bulk_candidates(bulk_candidates)
        if bulk_candidates is not None
        else pd.DataFrame(columns=EVIDENCE_COLUMNS)
    )
    sources, errors = collect_free_sources(
        species,
        include_gbif=gbif,
        include_wikimedia=wikimedia,
        include_openalex=openalex,
        include_europe_pmc=europe_pmc,
        include_web_descriptions=web_descriptions,
        include_world_flora=world_flora,
        max_descriptions=max_descriptions,
        openalex_per_page=openalex_per_page,
        pause_seconds=pause_seconds,
        wikipedia_languages=wikipedia_language,
        seed_evidence=seed_evidence,
        web_description_index_path=web_description_index,
    )
    report = write_search_outputs(
        output_dir,
        species,
        sources,
        errors,
        bulk_candidates=bulk_candidates,
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


@app.command("extract-text")
def extract_text_command(
    source_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    species_csv: Path | None = typer.Option(None, exists=True),
    species_column: str = typer.Option("accepted_species"),
    candidate_csv: Path | None = typer.Option(None, exists=True),
) -> None:
    """Build v1 rows from cached source text without making network calls."""
    sources = pd.read_csv(source_csv, dtype=str).fillna("")
    if species_csv is None:
        species = _unique(sources["accepted_species"])
    else:
        table = pd.read_csv(species_csv, dtype=str).fillna("")
        species = _species_list(table, species_column, len(table))
    errors = pd.DataFrame(columns=["species", "source", "error"])
    bulk_candidates = (
        pd.read_csv(candidate_csv, dtype=str).fillna("") if candidate_csv is not None else None
    )
    report = write_search_outputs(
        output_dir,
        species,
        sources,
        errors,
        bulk_candidates=bulk_candidates,
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
