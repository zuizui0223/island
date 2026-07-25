from __future__ import annotations

import pytest

from island_v2.plazi_trait_reacquisition import (
    build_plazi_fetch_url,
    extract_plazi_treatment_evidence,
)


def _xml() -> bytes:
    return b"""<document
      ID-DOI="10.1000/plant.1"
      masterDocId="MASTER-1"
      masterDocTitle="A revision of Plantus"
      zenodo-license-treatments="CC-BY">
      <treatment>
        <taxonomicName>Plantus alpha Author</taxonomicName>
        <subSubSection type="description">
          <paragraph>
            Flowers sessile; corolla 4.1-7.8 mm long; corolla tube 8-12 mm long.
          </paragraph>
        </subSubSection>
      </treatment>
    </document>"""


def test_plazi_exact_treatment_maps_controlled_measurement_classes() -> None:
    evidence = extract_plazi_treatment_evidence(
        _xml(),
        accepted_species="Plantus alpha",
        treatment_uuid="03A5194AFF97FFF5B0E22FDF5E59F94E",
    )

    values = set(zip(evidence["trait_name"], evidence["normalized_value"], strict=True))
    assert values == {
        ("flower_size_class", "small"),
        ("tube_depth_class", "intermediate"),
    }
    assert set(evidence["evidence_quality"]) == {"medium"}
    assert set(evidence["source_lineage"]) == {"doi:10.1000/plant.1"}
    assert set(evidence["source_license"]) == {"CC-BY"}
    assert set(evidence["inference_rule"]) == {""}


def test_plazi_treatment_rejects_wrong_species_and_invalid_uuid() -> None:
    with pytest.raises(ValueError, match="exact accepted species"):
        extract_plazi_treatment_evidence(
            _xml(),
            accepted_species="Othera beta",
            treatment_uuid="03A5194AFF97FFF5B0E22FDF5E59F94E",
        )
    with pytest.raises(ValueError, match="32 hexadecimal"):
        build_plazi_fetch_url("not-a-uuid")
