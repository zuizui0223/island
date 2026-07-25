from __future__ import annotations

import httpx
import pandas as pd

from island_v2.broad_web_medium_pilot import (
    FETCHED_SOURCE_COLUMNS,
    _canonical_url,
    _direct_evidence,
    _gobotany_sitemap_index,
    extract_direct_passages,
    parse_duckduckgo_urls,
    parse_gobotany_evidence,
    parse_pladias_evidence,
    select_all_unresolved_targets,
    select_pilot_targets,
)


def test_duckduckgo_parser_returns_destinations_without_snippets() -> None:
    html = """
    <div class="result">
      <a class="result__a"
         href="//duckduckgo.com/l/?uddg=https%3A%2F%2Fexample.edu%2Fflora%2Fone">
         Result one
      </a>
      <a class="result__snippet">Do not retain this snippet.</a>
    </div>
    <div class="result">
      <a class="result__a" href="https://example.gov/taxon/two?tracking=1">Two</a>
    </div>
    """
    assert parse_duckduckgo_urls(html, max_results=10) == [
        "https://example.edu/flora/one",
        "https://example.gov/taxon/two",
    ]


def test_extract_direct_passages_requires_exact_taxon_page_or_same_sentence() -> None:
    accepted = extract_direct_passages(
        accepted_species="Plantus alpha",
        search_name="Plantus alpha",
        target_axis="flower_colour",
        title="Plantus alpha - University Flora",
        text="Leaves are green. Flowers are bright yellow with five petals.",
    )
    assert accepted == ["Flowers are bright yellow with five petals."]

    rejected = extract_direct_passages(
        accepted_species="Plantus alpha",
        search_name="Plantus alpha",
        target_axis="flower_colour",
        title="Regional flora",
        text="Plantus beta has yellow flowers. Plantus alpha occurs nearby.",
    )
    assert rejected == []


def test_select_pilot_is_axis_balanced_and_genus_diverse() -> None:
    axes = [
        "floral_structural_complexity",
        "flower_colour",
        "reproductive_assurance",
    ]
    species = [
        "Alpha one",
        "Alpha two",
        "Beta one",
        "Gamma one",
        "Delta one",
        "Epsilon one",
        "Zeta one",
    ]
    missing_axes = {
        "Alpha one": axes[0],
        "Alpha two": axes[0],
        "Beta one": axes[0],
        "Gamma one": axes[1],
        "Delta one": axes[1],
        "Epsilon one": axes[2],
        "Zeta one": axes[2],
    }
    coverage_rows = []
    for name in species:
        for axis in axes:
            coverage_rows.append(
                {
                    "accepted_species": name,
                    "axis": axis,
                    "after_quality": "" if axis == missing_axes[name] else "medium",
                }
            )
    master = pd.DataFrame(
        {
            "accepted_species": species,
            "genus": [name.split()[0] for name in species],
            "family": ["Exampleaceae"] * len(species),
            "n_records": ["100", "90", "80", "70", "60", "50", "40"],
        }
    )
    direct = pd.DataFrame(
        {
            "accepted_species": ["Alpha donor", "Alpha donor2"],
            "axis": [axes[0], axes[0]],
            "quality": ["medium", "medium"],
            "evidence_scope": ["species_direct", "species_direct"],
        }
    )
    selected = select_pilot_targets(
        pd.DataFrame(coverage_rows),
        master,
        direct,
        limit=6,
    )
    assert len(selected) == 6
    assert selected["target_axis"].value_counts().to_dict() == {
        axes[0]: 2,
        axes[1]: 2,
        axes[2]: 2,
    }
    assert selected["accepted_species"].nunique() == len(selected)
    assert selected.loc[selected["target_axis"].eq(axes[0]), "genus"].nunique() == 2


def test_direct_evidence_accepts_institutional_page_and_rejects_generic_page() -> None:
    common = {
        "accepted_species": "Plantus alpha",
        "source_text": "Plantus alpha flowers are bright yellow.",
        "source_url": "",
        "source_citation": "University flora",
        "evidence_scope": "species_direct",
        "search_name": "Plantus alpha",
        "name_role": "accepted",
        "name_resolution_lineage": "",
        "target_axis": "flower_colour",
        "source_host": "",
        "retrieved_at_epoch": "1",
        "raw_content_sha256": "abc",
        "rights_status": "unknown",
        "robots_status": "allowed",
    }
    trusted = {
        **common,
        "source_url": "https://example.edu/flora/plantus-alpha",
        "source_type": "general_web_institutional",
        "source_host": "example.edu",
    }
    untrusted = {
        **common,
        "source_url": "https://example.com/plantus-alpha",
        "source_type": "general_web_unreviewed",
        "source_host": "example.com",
    }
    sources = pd.DataFrame([trusted, untrusted], columns=FETCHED_SOURCE_COLUMNS)
    _, direct = _direct_evidence(sources)
    assert len(direct) == 1
    assert direct.iloc[0]["axis"] == "flower_colour"
    assert direct.iloc[0]["normalized_value"] == "yellow_orange"
    assert direct.iloc[0]["evidence_quality"] == "medium"
    assert direct.iloc[0]["source_lineage"] == ("url:https://example.edu/flora/plantus-alpha")


def test_canonical_url_drops_query_and_fragment() -> None:
    assert _canonical_url("HTTPS://Example.ORG/a//b/?utm_source=x#part") == (
        "https://example.org/a/b"
    )


def test_pladias_parser_uses_only_structured_feature_values() -> None:
    html = """
    <html>
      <head><title>Plantus alpha – Pladias</title></head>
      <body>
        <li class="featureDetail">
          <div class="d-flex"><div>
            <span class="font-italic">Flower colour</span>
            <span class="font-italic">:</span>
            <b>green, white</b>
          </div></div>
        </li>
        <li class="featureDetail">
          <div class="d-flex"><div>
            <span class="font-italic">Flower symmetry</span>
            <span class="font-italic">:</span>
            <b>actinomorphic</b>
          </div></div>
        </li>
        <li class="featureDetail">
          <div class="d-flex"><div>
            <span class="font-italic">Generative reproduction type</span>
            <span class="font-italic">:</span>
            <b>allogamy</b>
          </div></div>
        </li>
        <aside>Help: flower colours include red, blue, yellow, and white.</aside>
      </body>
    </html>
    """
    evidence = parse_pladias_evidence(
        html,
        accepted_species="Plantus alpha",
        source_url="https://pladias.cz/en/taxon/data/Plantus+alpha",
    )
    assert set(zip(evidence["trait_name"], evidence["normalized_value"])) == {
        ("flower_primary_color", "multicolored_variable"),
        ("floral_symmetry", "actinomorphic"),
        ("mating_system", "predominantly_outcrossing"),
    }
    assert evidence["evidence_quality"].eq("medium").all()
    assert evidence["evidence_scope"].eq("species_direct").all()


def test_gobotany_parser_extracts_structured_colour_symmetry_and_size() -> None:
    html = """
    <html>
      <head><title>Plantus alpha (example plant): Go Botany</title></head>
      <body>
        <dl><dt>Flower petal color</dt><dd>white</dd></dl>
        <dl>
          <dt>Flower symmetry</dt>
          <dd>there are two or more ways to evenly divide the flower
              (the flower is radially symmetrical)</dd>
        </dl>
        <dl><dt>Flower diameter</dt><dd>12–20 mm</dd></dl>
      </body>
    </html>
    """
    evidence = parse_gobotany_evidence(
        html,
        accepted_species="Plantus alpha",
        source_url="https://gobotany.nativeplanttrust.org/species/plantus/alpha/",
    )
    assert set(zip(evidence["trait_name"], evidence["normalized_value"])) == {
        ("flower_primary_color", "white"),
        ("floral_symmetry", "actinomorphic"),
        ("flower_size_class", "medium"),
    }


def test_gobotany_parser_rejects_ambiguous_symmetry_value() -> None:
    html = """
    <html>
      <head><title>Plantus alpha (example plant): Go Botany</title></head>
      <body>
        <dl>
          <dt>Flower symmetry</dt>
          <dd>there are two or more ways to evenly divide some flowers and
              only one way to divide other flowers</dd>
        </dl>
      </body>
    </html>
    """
    evidence = parse_gobotany_evidence(
        html,
        accepted_species="Plantus alpha",
        source_url="https://gobotany.nativeplanttrust.org/species/plantus/alpha/",
    )
    assert evidence.empty


def test_structured_parsers_reject_nonmatching_species_page() -> None:
    html = """
    <html><head><title>Plantus beta – Pladias</title></head>
    <body>
      <li class="featureDetail"><div class="d-flex"><div>
        <span class="font-italic">Flower colour</span><b>yellow</b>
      </div></div></li>
    </body></html>
    """
    evidence = parse_pladias_evidence(
        html,
        accepted_species="Plantus alpha",
        source_url="https://pladias.cz/en/taxon/data/Plantus+alpha",
    )
    assert evidence.empty


def test_gobotany_sitemap_index_keeps_exact_species_pages_only() -> None:
    response = httpx.Response(
        200,
        request=httpx.Request(
            "GET",
            "https://gobotany.nativeplanttrust.org/sitemap.txt",
        ),
        text=(
            "https://gobotany.nativeplanttrust.org/species/plantus/alpha/\n"
            "https://gobotany.nativeplanttrust.org/species/plantus/beta/subsp-one/\n"
            "https://gobotany.nativeplanttrust.org/simple/plantus-alpha/\n"
        ),
    )

    class Client:
        def get(self, _url: str) -> httpx.Response:
            return response

    assert _gobotany_sitemap_index(Client()) == {
        "plantus alpha": ("https://gobotany.nativeplanttrust.org/species/plantus/alpha/")
    }


def test_select_all_unresolved_targets_keeps_each_species_once() -> None:
    coverage = pd.DataFrame(
        [
            {"accepted_species": "Alpha one", "axis": "flower_colour", "after_quality": ""},
            {
                "accepted_species": "Alpha one",
                "axis": "reproductive_assurance",
                "after_quality": "",
            },
            {"accepted_species": "Beta one", "axis": "flower_colour", "after_quality": "medium"},
            {
                "accepted_species": "Beta one",
                "axis": "reproductive_assurance",
                "after_quality": "",
            },
        ]
    )
    master = pd.DataFrame(
        {
            "accepted_species": ["Alpha one", "Beta one"],
            "genus": ["Alpha", "Beta"],
            "family": ["A", "B"],
            "n_records": ["10", "20"],
        }
    )
    selected = select_all_unresolved_targets(coverage, master)
    assert selected["accepted_species"].tolist() == ["Alpha one", "Beta one"]
    assert selected.set_index("accepted_species").loc["Alpha one", "target_axis"] == (
        "flower_colour|reproductive_assurance"
    )
