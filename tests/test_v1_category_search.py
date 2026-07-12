from __future__ import annotations

import pandas as pd

from island_v2.v1_category_search import (
    SOURCE_COLUMNS,
    _json_getter,
    _web_description_sources,
    _wikipedia_exact_sources,
    _world_flora_sources,
    build_priority_traits,
    collapse_likely_traits,
    collapse_v1_rows,
    evidence_from_bulk_candidates,
    extract_evidence_from_sources,
    infer_likely_traits,
    v2_reported_candidates,
)


class _Response:
    def __init__(self, html: str, status_code: int, url: str) -> None:
        self.text = html
        self.content = html.encode()
        self.status_code = status_code
        self.url = url

    def raise_for_status(self) -> None:
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")


class _Client:
    def get(self, url: str, params: object = None) -> _Response:
        del params
        if url.endswith("sitemap.xml"):
            return _Response(
                """
                <urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">
                  <url><loc>https://plants.ces.ncsu.edu/plants/plantus-hortensis/</loc></url>
                </urlset>
                """,
                200,
                url,
            )
        if "ncsu" in url:
            return _Response(
                """
                <html><body><h1>Plantus hortensis</h1>
                <h2>Plant Detail</h2>
                <dl><dt>Flower Color:</dt><dd>Pink</dd>
                <dt>Flower Shape:</dt><dd>Trumpet</dd></dl>
                <p>The flowers are pollinated by bees.</p>
                <p>The nectar attracts hummingbirds.</p>
                <script>Flowers are green and pollinated by wind.</script>
                </body></html>
                """,
                200,
                url,
            )
        return _Response("<html><body>Not found</body></html>", 404, url)


class _WorldFloraClient:
    def get(self, url: str, params: object = None) -> _Response:
        del params
        if "/search?" in url:
            return _Response(
                '<a title="Plantus florensis" href="/taxon/wfo-1;jsessionid=x">taxon</a>',
                200,
                url,
            )
        return _Response(
            "<html><body><h1>Plantus florensis</h1>"
            "<p>" + ("Botanical description. " * 20) + "</p>"
            "<p>Flowers have a white corolla tube and are pollinated by bees.</p>"
            "</body></html>",
            200,
            url,
        )


class _JsonResponse:
    def __init__(
        self,
        status_code: int,
        payload: dict[str, object],
        headers: dict[str, str] | None = None,
    ) -> None:
        self.status_code = status_code
        self.payload = payload
        self.headers = headers or {}

    def raise_for_status(self) -> None:
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")

    def json(self) -> dict[str, object]:
        return self.payload


class _JsonClient:
    def __init__(self, responses: list[_JsonResponse]) -> None:
        self.responses = responses
        self.calls = 0

    def get(self, url: str, params: object = None) -> _JsonResponse:
        del url, params
        response = self.responses[self.calls]
        self.calls += 1
        return response


def _sources(*rows: dict[str, str]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=SOURCE_COLUMNS).fillna("")


def test_json_getter_retries_rate_limit_with_retry_after(monkeypatch) -> None:
    client = _JsonClient(
        [
            _JsonResponse(429, {}, {"Retry-After": "0"}),
            _JsonResponse(200, {"ok": True}),
        ]
    )
    delays: list[float] = []
    monkeypatch.setattr("island_v2.v1_category_search.time.sleep", delays.append)

    payload = _json_getter(client)("https://example.org/api", {"q": "Plantus"})

    assert payload == {"ok": True}
    assert client.calls == 2
    assert delays == [0.0]


def test_wikipedia_exact_source_uses_scientific_name_without_wikidata() -> None:
    calls: list[tuple[str, dict[str, object]]] = []

    def getter(url: str, params: dict[str, object]) -> dict[str, object]:
        calls.append((url, params))
        return {
            "query": {
                "pages": {
                    "1": {
                        "title": "Plantus exactus",
                        "extract": "Plantus exactus has red tubular flowers.",
                    }
                }
            }
        }

    rows, errors = _wikipedia_exact_sources(
        ["Plantus exactus"],
        getter,
        pause_seconds=0,
    )

    assert not errors
    assert len(calls) == 1
    assert calls[0][1]["titles"] == "Plantus exactus"
    assert "wikidata.org" not in calls[0][0]
    assert rows[0]["evidence_scope"] == "species_direct"
    assert rows[0]["source_type"] == "wikipedia_extract"


def test_extracts_explicit_traits_with_floral_context() -> None:
    sources = _sources(
        {
            "accepted_species": "Plantus exampleii",
            "source_text": (
                "The leaves are green. Flowers are white and tubular. "
                "The species is pollinated by bumblebees and is self-incompatible."
            ),
            "source_url": "https://example.org/plant",
            "source_citation": "Example flora",
            "source_type": "gbif_species_description",
            "evidence_scope": "species_direct",
        }
    )

    evidence = extract_evidence_from_sources(sources)
    observed = set(zip(evidence["field"], evidence["value"], strict=True))

    assert ("flower_color", "white") in observed
    assert ("flower_color", "green") not in observed
    assert ("flower_shape", "tubular") in observed
    assert ("pollination_guild", "bumblebees") in observed
    assert ("self_incompatibility", "SI") in observed
    assert set(evidence["evidence_type"]) == {"flora"}


def test_extracts_globose_heads_and_spiciform_inflorescences() -> None:
    sources = _sources(
        {
            "accepted_species": "Plantus capitatus",
            "source_text": "The flowers are borne in spherical heads on short peduncles.",
            "source_url": "https://example.org/capitatus",
            "source_citation": "Example flora",
            "source_type": "flora_or_monograph",
            "evidence_scope": "species_direct",
        },
        {
            "accepted_species": "Plantus spicatus",
            "source_text": "The female inflorescences are spiciform and densely flowered.",
            "source_url": "https://example.org/spicatus",
            "source_citation": "Example flora",
            "source_type": "flora_or_monograph",
            "evidence_scope": "species_direct",
        },
    )

    evidence = extract_evidence_from_sources(sources)
    observed = set(zip(evidence["species"], evidence["value"], strict=True))

    assert ("Plantus capitatus", "globose / spherical flower head") in observed
    assert ("Plantus spicatus", "spike / elongated inflorescence") in observed


def test_rejects_nonfloral_powdery_bloom_color_near_flowering_time() -> None:
    sources = _sources(
        {
            "accepted_species": "Plantus powderii",
            "source_text": (
                "The phyllodes are distinctly covered with a white powdery bloom, "
                "and flowers from August to September."
            ),
            "source_url": "https://example.org/plantus-powderii",
            "source_citation": "Example flora",
            "source_type": "flora_or_monograph",
            "evidence_scope": "species_direct",
        },
        {
            "accepted_species": "Plantus ripeningii",
            "source_text": (
                "The flower balls are white, then form hooked spikes which when ripe "
                "brown off."
            ),
            "source_url": "https://example.org/plantus-ripeningii",
            "source_citation": "Example flora",
            "source_type": "flora_or_monograph",
            "evidence_scope": "species_direct",
        },
        {
            "accepted_species": "Plantus commonnameii",
            "source_text": (
                "Plantus commonnameii, the Japanese red maple, has a local name meaning "
                "flower maple and is a rare shrub."
            ),
            "source_url": "https://example.org/plantus-commonnameii",
            "source_citation": "Example encyclopedia",
            "source_type": "wikipedia_extract",
            "evidence_scope": "species_direct",
        },
        {
            "accepted_species": "Plantus medicinalis",
            "source_text": "An infusion of the open flower is drunk as a headache remedy.",
            "source_url": "https://example.org/plantus-medicinalis",
            "source_citation": "Example horticulture account",
            "source_type": "horticulture_site",
            "evidence_scope": "species_direct",
        },
        {
            "accepted_species": "Plantus spinosa",
            "source_text": "The flowers are white. The mature capitula have yellow spines.",
            "source_url": "https://example.org/plantus-spinosa",
            "source_citation": "Example flora",
            "source_type": "flora_or_monograph",
            "evidence_scope": "species_direct",
        },
        {
            "accepted_species": "Plantus capsularis",
            "source_text": "The plant has white flowers and broadly cup-shaped capsules.",
            "source_url": "https://example.org/plantus-capsularis",
            "source_citation": "Example flora",
            "source_type": "flora_or_monograph",
            "evidence_scope": "species_direct",
        },
    )

    evidence = extract_evidence_from_sources(sources)

    powder = evidence.loc[evidence["species"].eq("Plantus powderii")]
    ripening = evidence.loc[evidence["species"].eq("Plantus ripeningii")]
    common_name = evidence.loc[evidence["species"].eq("Plantus commonnameii")]
    medicinal = evidence.loc[evidence["species"].eq("Plantus medicinalis")]
    spiny = evidence.loc[evidence["species"].eq("Plantus spinosa")]
    capsular = evidence.loc[evidence["species"].eq("Plantus capsularis")]
    assert "flower_color" not in set(powder["field"])
    assert set(ripening.loc[ripening["field"].eq("flower_color"), "value"]) == {"white"}
    assert "flower_color" not in set(common_name["field"])
    assert "flower_shape" not in set(medicinal["field"])
    assert set(spiny.loc[spiny["field"].eq("flower_color"), "value"]) == {"white"}
    assert set(capsular.loc[capsular["field"].eq("flower_color"), "value"]) == {"white"}
    assert "flower_shape" not in set(capsular["field"])


def test_rejects_internal_organ_and_cultivar_colors() -> None:
    sources = _sources(
        {
            "accepted_species": "Plantus antheralis",
            "source_text": (
                "The flowers are red with purple anthers and a black stigma."
            ),
            "source_url": "https://example.org/plantus-antheralis",
            "source_citation": "Example flora",
            "source_type": "flora_or_monograph",
            "evidence_scope": "species_direct",
        },
        {
            "accepted_species": "Plantus cultivarensis",
            "source_text": (
                "The flowers are pink to red. Cultivars / Varieties: 'Alba'. "
                "White flowers."
            ),
            "source_url": "https://example.org/plantus-cultivarensis",
            "source_citation": "Example horticulture account",
            "source_type": "horticulture_site",
            "evidence_scope": "species_direct",
        },
    )

    evidence = extract_evidence_from_sources(sources)

    antheralis = evidence.loc[
        evidence["species"].eq("Plantus antheralis")
        & evidence["field"].eq("flower_color")
    ]
    cultivarensis = evidence.loc[
        evidence["species"].eq("Plantus cultivarensis")
        & evidence["field"].eq("flower_color")
    ]
    assert set(antheralis["value"]) == {"red"}
    assert set(cultivarensis["value"]) == {"pink", "red"}


def test_ignores_indirect_and_negated_claims() -> None:
    sources = _sources(
        {
            "accepted_species": "Plantus exampleii",
            "source_text": "Flowers are not pollinated by bees and have red petals.",
            "source_url": "https://example.org/indirect",
            "source_citation": "Indirect page",
            "source_type": "wikipedia_extract",
            "evidence_scope": "species_indirect",
        },
        {
            "accepted_species": "Plantus exampleii",
            "source_text": "The flowers are not pollinated by bees.",
            "source_url": "https://example.org/direct",
            "source_citation": "Direct page",
            "source_type": "wikipedia_extract",
            "evidence_scope": "species_direct",
        },
    )

    evidence = extract_evidence_from_sources(sources)

    assert evidence.empty


def test_collapses_to_legacy_schema_and_keeps_unknowns() -> None:
    sources = _sources(
        {
            "accepted_species": "Plantus exampleii",
            "source_text": (
                "Flowers have red petals and are bilaterally symmetric. "
                "They are pollinated by bees. Mixed mating and self-compatible."
            ),
            "source_url": "https://example.org/article",
            "source_citation": "Field experiment 2024",
            "source_type": "openalex_title_abstract",
            "evidence_scope": "species_direct",
        }
    )
    evidence = extract_evidence_from_sources(sources)

    result = collapse_v1_rows(["Plantus exampleii", "Absentia nulla"], evidence)

    first = result.iloc[0]
    assert first["flower_color"] == "red"
    assert first["flower_shape"] == "zygomorphic / bilaterally symmetric"
    assert first["pollination_guild"] == "bees"
    assert first["mating_system"] == "mixed_mating"
    assert first["self_incompatibility"] == "SC"
    assert first["evidence_type"] == "field_study"
    assert first["confidence"] == "high"
    assert "pollinated by bees" in first["pollination_notes"]

    second = result.iloc[1]
    assert second["flower_color"] == "unknown"
    assert second["pollination_guild"] == "unknown"
    assert second["evidence_type"] == "inference"
    assert second["confidence"] == "low"


def test_conflicting_guilds_collapse_to_mixed() -> None:
    sources = _sources(
        {
            "accepted_species": "Plantus generalis",
            "source_text": "Flowers are pollinated by bees and pollinated by flies.",
            "source_url": "https://example.org/generalist",
            "source_citation": "Pollination study",
            "source_type": "openalex_title_abstract",
            "evidence_scope": "species_direct",
        }
    )

    evidence = extract_evidence_from_sources(sources)
    result = collapse_v1_rows(["Plantus generalis"], evidence)

    assert result.loc[0, "pollination_guild"] == "mixed"


def test_translates_eol_style_bulk_candidates_without_overclaiming() -> None:
    candidates = pd.DataFrame(
        [
            {
                "accepted_species": "Plantus bulkus",
                "trait_name": "flower_primary_color",
                "candidate_value": "red_pink",
                "source_url": "https://eol.org/pages/1/data",
                "source_excerpt": "EOL TraitBank: flower color = red",
            },
            {
                "accepted_species": "Plantus bulkus",
                "trait_name": "self_incompatibility",
                "candidate_value": "SI",
                "source_url": "https://eol.org/pages/1/data",
                "source_excerpt": "EOL TraitBank: self incompatibility = SI",
            },
            {
                "accepted_species": "Plantus bulkus",
                "trait_name": "sex_system",
                "candidate_value": "dioecious",
                "source_url": "https://eol.org/pages/1/data",
            },
        ]
    )

    evidence = evidence_from_bulk_candidates(candidates)
    result = collapse_v1_rows(["Plantus bulkus"], evidence)

    assert len(evidence) == 2
    assert result.loc[0, "flower_color"] == "red/pink"
    assert result.loc[0, "self_incompatibility"] == "SI"
    assert result.loc[0, "mating_system"] == "unknown"
    assert result.loc[0, "evidence_type"] == "inference"
    assert result.loc[0, "confidence"] == "low"


def test_horticulture_page_is_verified_and_extracted_with_provenance() -> None:
    sources, errors = _web_description_sources(
        ["Plantus hortensis"],
        _Client(),
        pause_seconds=0,
    )
    source_table = pd.DataFrame(sources, columns=SOURCE_COLUMNS)
    evidence = extract_evidence_from_sources(source_table)
    result = collapse_v1_rows(["Plantus hortensis"], evidence)

    assert not errors
    assert len(sources) == 1
    assert result.loc[0, "flower_color"] == "pink"
    assert result.loc[0, "flower_shape"] == "funnel / trumpet-shaped"
    assert result.loc[0, "pollination_guild"] == "mixed"
    assert result.loc[0, "evidence_type"] == "horticulture|inference"
    assert "green" not in set(evidence["value"])
    inferred = evidence.loc[evidence["rule_id"].eq("guild_attracts_birds")].iloc[0]
    assert inferred["confidence"] == "low"
    assert inferred["source_tier"] == "C_horticulture"


def test_warm_tubular_flower_is_only_emitted_as_likely_guild() -> None:
    sources = _sources(
        {
            "accepted_species": "Plantus rubritubus",
            "source_text": "The flowers are red with a long tubular corolla.",
            "source_url": "https://example.org/flora/rubritubus",
            "source_citation": "Example flora",
            "source_type": "flora_or_monograph",
            "evidence_scope": "species_direct",
        }
    )
    evidence = extract_evidence_from_sources(sources)
    raw = collapse_v1_rows(["Plantus rubritubus"], evidence)
    likely = infer_likely_traits(raw, evidence)
    wide = collapse_likely_traits(["Plantus rubritubus"], likely)

    assert raw.loc[0, "pollination_guild"] == "unknown"
    assert wide.loc[0, "likely_pollination_guild"] == "likely_birds_or_butterflies"
    row = likely.loc[likely["rule_id"].eq("likely_guild_warm_tubular")].iloc[0]
    assert row["basis_fields"] == "flower_color|flower_shape"
    assert row["confidence"] == "low"
    assert row["basis_source_urls"] == "https://example.org/flora/rubritubus"


def test_composite_tubular_florets_do_not_emit_warm_tube_likely_guild() -> None:
    sources = _sources(
        {
            "accepted_species": "Plantus compositus",
            "source_text": "Capitula have pink tubular florets and white ray florets.",
            "source_url": "https://example.org/flora/compositus",
            "source_citation": "Example flora",
            "source_type": "flora_or_monograph",
            "evidence_scope": "species_direct",
        }
    )

    evidence = extract_evidence_from_sources(sources)
    raw = collapse_v1_rows(["Plantus compositus"], evidence)
    likely = infer_likely_traits(raw, evidence)

    assert "tubular" in raw.loc[0, "flower_shape"]
    assert "composite head" in raw.loc[0, "flower_shape"]
    assert "likely_guild_warm_tubular" not in set(likely["rule_id"])


def test_autogamy_yields_likely_selfing_and_likely_sc_only() -> None:
    sources = _sources(
        {
            "accepted_species": "Plantus autogamus",
            "source_text": "Flowers show autonomous selfing under cultivation.",
            "source_url": "https://example.org/autogamy",
            "source_citation": "Cultivation account",
            "source_type": "horticulture_site",
            "evidence_scope": "species_direct",
        }
    )
    evidence = extract_evidence_from_sources(sources)
    raw = collapse_v1_rows(["Plantus autogamus"], evidence)
    likely = infer_likely_traits(raw, evidence)
    wide = collapse_likely_traits(["Plantus autogamus"], likely)
    priority = build_priority_traits(raw, wide)

    assert raw.loc[0, "mating_system"] == "unknown"
    assert raw.loc[0, "self_incompatibility"] == "likely_SC"
    assert wide.loc[0, "likely_selfing_status"] == "likely_selfing_capable"
    assert wide.loc[0, "likely_self_incompatibility"] == "likely_SC"
    assert priority.loc[0, "selfing_status_mode"] == "likely"
    assert priority.loc[0, "self_incompatibility_mode"] == "likely"


def test_v2_export_does_not_promote_attraction_to_pollination_guild() -> None:
    sources = _sources(
        {
            "accepted_species": "Plantus visitatus",
            "source_text": ("The flowers are pollinated by bees and also attract hummingbirds."),
            "source_url": "https://example.org/visitatus",
            "source_citation": "Example flora",
            "source_type": "flora_or_monograph",
            "evidence_scope": "species_direct",
        }
    )
    evidence = extract_evidence_from_sources(sources)
    raw = collapse_v1_rows(["Plantus visitatus"], evidence)
    v2 = v2_reported_candidates(evidence)

    assert raw.loc[0, "pollination_guild"] == "mixed"
    guild = v2.loc[v2["trait_name"].eq("pollination_functional_guild")]
    assert list(guild["provisional_candidate_value"]) == ["other_bees"]


def test_v2_export_maps_color_shape_and_autogamy_to_canonical_traits() -> None:
    sources = _sources(
        {
            "accepted_species": "Plantus v2us",
            "source_text": ("Flowers are red and tubular and show autonomous selfing."),
            "source_url": "https://example.org/v2us",
            "source_citation": "Example flora",
            "source_type": "flora_or_monograph",
            "evidence_scope": "species_direct",
        }
    )

    v2 = v2_reported_candidates(extract_evidence_from_sources(sources))
    observed = set(zip(v2["trait_name"], v2["provisional_candidate_value"], strict=True))

    assert ("flower_primary_color", "red_pink") in observed
    assert ("floral_form", "tubular") in observed
    assert ("autonomous_selfing_capacity", "autonomous") in observed
    assert ("self_incompatibility", "SC") not in observed
    assert set(v2["candidate_class"]) == {"reported"}


def test_unrelated_trumpet_flower_name_is_not_a_shape_claim() -> None:
    sources = _sources(
        {
            "accepted_species": "Plantus vinea",
            "source_text": (
                "Consider other vines such as Trumpet Flower. "
                "The flowers of Plantus vinea are pink."
            ),
            "source_url": "https://example.org/plantus-vinea",
            "source_citation": "Plantus vinea - horticulture guide",
            "source_type": "horticulture_site",
            "evidence_scope": "species_direct",
        }
    )

    evidence = extract_evidence_from_sources(sources)

    assert "flower_shape" not in set(evidence["field"])


def test_world_flora_precedes_lower_tier_horticulture_values() -> None:
    flora_rows, errors = _world_flora_sources(
        ["Plantus florensis"],
        _WorldFloraClient(),
        pause_seconds=0,
    )
    horticulture = {
        "accepted_species": "Plantus florensis",
        "source_text": "The flowers of Plantus florensis are red.",
        "source_url": "https://example.org/horticulture",
        "source_citation": "Plantus florensis horticulture page",
        "source_type": "horticulture_site",
        "evidence_scope": "species_direct",
    }
    sources = pd.DataFrame([*flora_rows, horticulture], columns=SOURCE_COLUMNS)
    evidence = extract_evidence_from_sources(sources)
    result = collapse_v1_rows(["Plantus florensis"], evidence)

    assert not errors
    assert set(evidence["source_tier"]) == {"A_flora", "C_horticulture"}
    assert result.loc[0, "flower_color"] == "white"
    assert result.loc[0, "flower_shape"] == "tubular"
    assert result.loc[0, "pollination_guild"] == "bees"
    assert result.loc[0, "evidence_type"] == "flora"
