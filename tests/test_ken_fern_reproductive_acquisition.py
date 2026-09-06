from analysis.acquire_ken_fern_reproductive_traits_20260808 import parse_page


def _page(*, pollinators: str, cultivation: str = "Wild") -> bytes:
    return f"""
    <div class="latin_name"><h1>Alpha one</h1></div>
    <div class="family"><h4>Exampleaceae</h4></div>
    <table class="tableProperties">
      <tr><td>Pollinators</td><td>{pollinators}</td></tr>
      <tr><td>Self-fertile</td><td>Yes</td></tr>
      <tr><td>Cultivation Status</td><td>{cultivation}</td></tr>
    </table>
    """.encode()


def test_wild_sexual_self_fertile_maps_only_to_sc() -> None:
    row = parse_page(
        "Alpha one", "Exampleaceae", "https://example.test/a", _page(pollinators="Bees, Self")
    )
    assert row["status"] == "accepted_wild_self_fertile"
    assert row["trait_name"] == "self_incompatibility"
    assert row["normalized_value"] == "SC"


def test_apomictic_self_fertile_is_not_mapped_to_sc() -> None:
    row = parse_page(
        "Alpha one", "Exampleaceae", "https://example.test/a", _page(pollinators="Apomictic")
    )
    assert row["status"] == "nonsexual_reproduction_category_error"
    assert row["trait_name"] == ""
    assert row["normalized_value"] == ""
