import pandas as pd

from island_v2.regional_flora_network_acquisition import (
    SITES,
    _terminal_task_rows,
    parse_family_index,
    parse_species_page,
)


def test_family_index_parser_is_provider_generic() -> None:
    payload = b"""
    <html><h1>Checklist: Acanthaceae</h1>
    <a href="../species.php?species_id=152940">Thunbergia alata</a>
    <a href="../species.php?species_id=152940">Thunbergia alata</a>
    <a href="../species.php?species_id=152950">*Thunbergia annua</a>
    <a href="../species.php?species_id=1">Thunbergia alata var. alba</a>
    </html>
    """
    family, rows = parse_family_index(
        payload,
        site=SITES["zambia"],
        family_id=1,
    )
    assert family == "Acanthaceae"
    assert [(row["source_scientific_name"], row["species_id"]) for row in rows] == [
        ("Thunbergia alata", "152940"),
        ("Thunbergia annua", "152950"),
    ]
    assert all(row["domain"] == "zambiaflora.com" for row in rows)


def test_expanded_african_flora_network_is_registered() -> None:
    site = SITES["zimbabwe"]
    assert site.provider == "Flora of Zimbabwe"
    assert site.domain == "zimbabweflora.co.zw"
    assert site.family_url(1).startswith("https://www.zimbabweflora.co.zw/")
    assert {
        "caprivi",
        "burundi",
        "drcongo",
        "rwanda",
        "angola",
        "kenya",
        "tanzania",
        "uganda",
    }.issubset(SITES)
    assert len({site.domain for site in SITES.values()}) == len(SITES)


def _species_html(*, name: str = "Thunbergia alata", family: str = "Acanthaceae") -> bytes:
    return f"""
    <html><head><title>Flora: Species information: {name}</title></head><body>
    <a href="family.php?family_id=1">{family}</a>
    <h1>{name}<span class="author"> Bojer ex Sims</span></h1>
    <table><tr><td><a href="about.php#descr">Description:</a></td><td>
    Flowers solitary or in pairs. Corolla campanulate, 25-35 mm long,
    pale yellow to deep orange with a deep purple throat.
    </td></tr></table>
    </body></html>
    """.encode()


def test_species_page_requires_exact_identity_and_preserves_excerpts() -> None:
    audit, rows = parse_species_page(
        _species_html(),
        site=SITES["zambia"],
        expected_species="Thunbergia alata",
        expected_family="Acanthaceae",
        source_url="https://www.zambiaflora.com/speciesdata/species.php?species_id=152940",
        missing_traits=(
            "flower_primary_color",
            "floral_form",
            "flower_size_class",
            "inflorescence_display",
        ),
        retrieved_at_utc="2026-08-14T00:00:00Z",
    )
    assert audit["status"] == "candidate_extracted"
    evidence = {row["trait_name"]: row for row in rows}
    assert evidence["flower_primary_color"]["normalized_value"] == "multicolored_variable"
    assert evidence["floral_form"]["normalized_value"] == "bell_campanulate"
    assert evidence["flower_size_class"]["normalized_value"] == "medium"
    assert evidence["inflorescence_display"]["normalized_value"] == "solitary"
    assert all(row["review_status"] == "pending_individual_review" for row in rows)
    assert all(row["source_excerpt"] for row in rows)
    assert all(len(row["content_sha256"]) == 64 for row in rows)


def test_species_page_rejects_name_or_family_conflict() -> None:
    name_audit, name_rows = parse_species_page(
        _species_html(name="Thunbergia annua"),
        site=SITES["malawi"],
        expected_species="Thunbergia alata",
        expected_family="Acanthaceae",
        source_url="https://www.malawiflora.com/speciesdata/species.php?species_id=1",
        missing_traits=("flower_primary_color",),
        retrieved_at_utc="2026-08-14T00:00:00Z",
    )
    assert name_audit["status"] == "rejected_name_mismatch"
    assert name_rows == []

    family_audit, family_rows = parse_species_page(
        _species_html(family="Lamiaceae"),
        site=SITES["botswana"],
        expected_species="Thunbergia alata",
        expected_family="Acanthaceae",
        source_url="https://www.botswanaflora.com/speciesdata/species.php?species_id=1",
        missing_traits=("flower_primary_color",),
        retrieved_at_utc="2026-08-14T00:00:00Z",
    )
    assert family_audit["status"] == "rejected_family_mismatch"
    assert family_rows == []


def test_resume_keeps_permanent_no_hits_but_retries_transient_failures() -> None:
    tasks = pd.DataFrame(
        {
            "status": ["success", "fetch_failed", "fetch_failed", "fetch_failed", ""],
            "error": ["", "http_status_404", "http_status_503", "ReadTimeout:x", ""],
        }
    )
    assert _terminal_task_rows(tasks).tolist() == [True, True, False, False, False]
