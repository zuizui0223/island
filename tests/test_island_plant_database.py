from __future__ import annotations

import json
from pathlib import Path

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import Polygon

from island_v2.island_plant_database import build_islands_core, validate_bundle


CONTRACT = Path("config/island_plant_database_v2.yml")


def _write_bundle(tmp_path: Path, *, with_traits: bool = False) -> Path:
    bundle = tmp_path / "bundle"
    bundle.mkdir()

    pd.DataFrame(
        [
            {
                "island_id": "isl-1",
                "source_island_id": "gshhg-1",
                "island_name": "Example Island",
                "archipelago": "",
                "country_or_territory": "",
                "area_km2": "100",
                "centroid_lat": "1",
                "centroid_lon": "2",
                "geometry_source": "GSHHG",
                "geometry_version": "2.3.7",
                "geometry_sha256": "0" * 64,
                "release_status": "redistributable",
                "source_license": "public-source",
            }
        ]
    ).to_csv(bundle / "islands.csv", index=False)

    pd.DataFrame(
        [
            {
                "taxon_id": "tax-1",
                "accepted_name": "Alpha one",
                "authorship": "",
                "taxonomic_rank": "species",
                "family": "Alphaaceae",
                "genus": "Alpha",
                "backbone_name": "GBIF",
                "backbone_key": "1",
                "backbone_version": "test",
                "taxonomic_status": "accepted",
                "release_status": "redistributable",
                "source_license": "CC-BY-4.0",
            }
        ]
    ).to_csv(bundle / "taxa.csv", index=False)

    pd.DataFrame(
        [
            {
                "island_id": "isl-1",
                "taxon_id": "tax-1",
                "membership_status": "candidate",
                "establishment_status": "unknown",
                "endemic_status": "unknown",
                "first_record_year": "2000",
                "last_record_year": "2025",
                "occurrence_record_count": "3",
                "specimen_record_count": "1",
                "evidence_count": "1",
                "review_status": "unreviewed",
                "release_status": "redistributable",
                "provenance_id": "prov-1",
            }
        ]
    ).to_csv(bundle / "island_taxa.csv", index=False)

    pd.DataFrame(
        [
            {
                "evidence_id": "ev-1",
                "island_id": "isl-1",
                "taxon_id": "tax-1",
                "source_type": "gbif_download",
                "source_record_id": "occ-1",
                "source_url": "https://example.org/occ-1",
                "dataset_key_or_doi": "doi:test",
                "basis_of_record": "PRESERVED_SPECIMEN",
                "event_date": "2025-01-01",
                "coordinate_uncertainty_m": "10",
                "evidence_role": "supports_presence",
                "source_license": "CC-BY-4.0",
                "rights_status": "redistributable",
                "review_status": "unreviewed",
            }
        ]
    ).to_csv(bundle / "evidence.csv", index=False)

    pd.DataFrame(
        [
            {
                "object_type": "evidence",
                "object_id": "ev-1",
                "source_id": "occ-1",
                "source_url": "https://example.org/occ-1",
                "source_license": "CC-BY-4.0",
                "rights_status": "redistributable",
                "rights_evidence": "source metadata",
                "redistributable_value": "true",
                "redistributable_provenance": "true",
                "review_status": "accepted",
            }
        ]
    ).to_csv(bundle / "RIGHTS_LEDGER.csv", index=False)

    if with_traits:
        pd.DataFrame(
            [
                {
                    "trait_record_id": "tr-1",
                    "taxon_id": "tax-1",
                    "island_id": "",
                    "trait_name": "flower_colour",
                    "trait_value": "red",
                    "trait_unit": "",
                    "trait_ontology_version": "1",
                    "evidence_id": "ev-1",
                    "quality": "high",
                    "source_license": "CC-BY-4.0",
                    "rights_status": "redistributable",
                    "release_status": "redistributable",
                }
            ]
        ).to_csv(bundle / "traits.csv", index=False)

    return bundle


def _write_island_source(tmp_path: Path, backend: str) -> tuple[Path, Path]:
    gpkg = tmp_path / "islands.gpkg"
    geometry = Polygon([(140.0, 35.0), (140.2, 35.0), (140.2, 35.2), (140.0, 35.2)])
    frame = gpd.GeoDataFrame(
        [
            {
                "island_id": "gshhg_2.3.7_h_abc",
                "source_label": "gshhg_2.3.7_h",
                "parent_feature_id": "42",
                "part_index": 1,
                "island_name": "Example Island",
                "area_km2": 400.0,
                "geometry_sha256": "legacy-short-hash",
                "landmass_rule": "5 <= area_km2 <= 7000000",
            }
        ],
        geometry=[geometry],
        crs=4326,
    )
    frame.to_file(gpkg, layer="islands", driver="GPKG")
    policy = tmp_path / "source_policy.json"
    policy.write_text(
        json.dumps(
            {
                "source_backend": backend,
                "source_version": "2.3.7",
            }
        ),
        encoding="utf-8",
    )
    return gpkg, policy


def test_valid_alpha1_bundle(tmp_path: Path) -> None:
    bundle = _write_bundle(tmp_path, with_traits=True)
    report = validate_bundle(bundle, CONTRACT)
    assert report == {
        "database_id": "global_island_plant_database",
        "version": "2.0.0-alpha1",
        "islands": 1,
        "taxa": 1,
        "island_taxa": 1,
        "evidence": 1,
        "traits": 1,
        "rights_ledger_present": True,
        "valid": True,
    }


def test_unknown_island_taxon_reference_fails(tmp_path: Path) -> None:
    bundle = _write_bundle(tmp_path)
    frame = pd.read_csv(bundle / "island_taxa.csv", dtype=str)
    frame.loc[0, "taxon_id"] = "tax-missing"
    frame.to_csv(bundle / "island_taxa.csv", index=False)
    with pytest.raises(ValueError, match="unknown taxa"):
        validate_bundle(bundle, CONTRACT)


def test_trait_requires_existing_evidence(tmp_path: Path) -> None:
    bundle = _write_bundle(tmp_path, with_traits=True)
    frame = pd.read_csv(bundle / "traits.csv", dtype=str)
    frame.loc[0, "evidence_id"] = "ev-missing"
    frame.to_csv(bundle / "traits.csv", index=False)
    with pytest.raises(ValueError, match="unknown evidence"):
        validate_bundle(bundle, CONTRACT)


def test_gshhg_island_export_is_public_and_hash_locked(tmp_path: Path) -> None:
    gpkg, policy = _write_island_source(tmp_path, "gshhg")
    frame = build_islands_core(gpkg, policy, CONTRACT)
    assert len(frame) == 1
    row = frame.iloc[0]
    assert row["source_island_id"] == "gshhg_2.3.7_h:42"
    assert row["geometry_source"] == "GSHHG"
    assert row["source_license"] == "LGPL-3.0-or-later"
    assert row["release_status"] == "redistributable"
    assert len(row["geometry_sha256"]) == 64
    assert float(row["centroid_lat"]) > 35.0
    assert float(row["centroid_lon"]) > 140.0


def test_natural_earth_fallback_is_marked_public_domain(tmp_path: Path) -> None:
    gpkg, policy = _write_island_source(tmp_path, "natural_earth_10m_fallback")
    frame = build_islands_core(gpkg, policy, CONTRACT)
    assert frame.loc[0, "geometry_source"] == "Natural Earth 10m"
    assert frame.loc[0, "source_license"] == "PUBLIC-DOMAIN"
