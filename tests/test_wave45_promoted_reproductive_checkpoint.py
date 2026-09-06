from __future__ import annotations

import pandas as pd
import pytest

from island_v2.wave45_promoted_reproductive_checkpoint import _validate_evidence


def _row() -> dict[str, str]:
    return {
        "accepted_species": "Example plant",
        "axis": "reproductive_assurance",
        "trait_name": "self_incompatibility",
        "normalized_value": "SC",
        "quality": "high",
        "source_group": "wave45_primary_reproductive_articles",
        "source_provider": "Example primary article",
        "source_url": "https://doi.org/10.1234/example",
        "source_record_id": "wave45:example:self-incompatibility",
        "source_citation": "Example et al. 2026",
        "source_excerpt": "the species is self-compatible",
        "evidence_scope": "species_direct",
        "name_match_method": "accepted_name_exact",
        "source_lineage": "doi:10.1234/example",
        "lineage_method": "original_article_doi",
        "source_run_id": "web:retrieved-20260831",
        "source_artifact": "example",
        "source_file": "html_sha256:example",
        "acceptance_contract": "primary_article_exact_species_v1",
    }


def test_wave45_checkpoint_keeps_reproductive_traits_distinct() -> None:
    _validate_evidence(pd.DataFrame([_row()]), label="test")

    wrong = _row()
    wrong["trait_name"] = "floral_form"
    with pytest.raises(ValueError, match="interchanges distinct reproductive traits"):
        _validate_evidence(pd.DataFrame([wrong]), label="test")


def test_wave45_checkpoint_keeps_auxiliary_traits_outside_strict_axes() -> None:
    wrong = _row()
    wrong["axis"] = "floral_structural_complexity"
    wrong["trait_name"] = "reward_type"
    with pytest.raises(ValueError, match="independent auxiliary trait"):
        _validate_evidence(pd.DataFrame([wrong]), label="test")
