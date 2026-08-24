from pathlib import Path

import pandas as pd

from island_v2.lucid_rainforest_checkpoint import (
    _page_fields,
    _write_gzip,
    build_evidence,
    deterministic_audit_queue,
)
from island_v2.trait_measurements import load_rules

PAGE = b"""
<html><body>
<h1 class="bannertitle">Glochidion testum Author</h1>
<div id="family"><div class="content">Phyllanthaceae</div></div>
<div id="flowers"><div class="content">
Flowers white, tubular, solitary, about 4 mm diam.
</div></div>
</body></html>
"""


def test_page_fields_use_only_identity_family_and_flower_section() -> None:
    fields = _page_fields(PAGE)
    assert fields == {
        "page_title": "Glochidion testum Author",
        "family": "Phyllanthaceae",
        "flower_text": "Flowers white, tubular, solitary, about 4 mm diam.",
    }


def test_build_evidence_is_exact_family_gated_and_incremental(tmp_path: Path) -> None:
    cache = tmp_path / "cache"
    _write_gzip(cache / "pages" / "glochidion_testum.html.gz", PAGE)
    selected = pd.DataFrame(
        [{"accepted_species": "Glochidion testum", "href": "glochidion_testum.htm"}]
    )
    master = pd.DataFrame(
        [{"accepted_species": "Glochidion testum", "family": "Phyllanthaceae"}]
    )
    evidence, inventory = build_evidence(
        selected,
        master,
        cache,
        {("Glochidion testum", "flower_primary_color")},
        load_rules(Path("config/measurement_classification.yml")),
    )
    assert inventory.loc[0, "identity_gate_passed"]
    assert "flower_primary_color" not in set(evidence["trait_name"])
    assert {"floral_form", "flower_size_class", "inflorescence_display"}.issubset(
        set(evidence["trait_name"])
    )
    assert evidence["source_excerpt"].str.contains("Flowers white").all()


def test_audit_queue_is_deterministic_and_keeps_review_fields_blank() -> None:
    evidence = pd.DataFrame(
        [
            {
                "source_record_id": f"candidate-{index:03d}",
                "accepted_species": f"Example species{index}",
                "trait_name": "floral_form",
                "normalized_value": "tubular",
                "source_provider": "provider",
                "source_url": f"https://example.org/{index}",
                "source_lineage": f"lineage:{index}",
                "source_citation": "citation",
                "source_excerpt": "Flowers tubular.",
                "name_match_method": "accepted_name_exact_plus_family_exact",
            }
            for index in range(220)
        ]
    )
    first = deterministic_audit_queue(evidence)
    second = deterministic_audit_queue(evidence)
    pd.testing.assert_frame_equal(first, second)
    assert len(first) == 200
    assert first["accepted_correct"].eq("").all()
    assert first["reviewer"].eq("").all()
