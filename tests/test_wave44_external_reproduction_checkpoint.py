from __future__ import annotations

import pandas as pd
import pytest

from island_v2.wave44_external_reproduction_checkpoint import (
    DIRECT_CONTRACT,
    EXTERNAL_CONTRACT,
    _validate_evidence,
)


def _evidence_row(*, external: bool) -> dict[str, str]:
    species = "Other plant" if external else "Target plant"
    provider = "Lin et al. 2025 GloPL Supplementary Data"
    return {
        "accepted_species": species,
        "axis": "reproductive_assurance",
        "trait_name": "self_incompatibility",
        "normalized_value": "SC",
        "quality": "high",
        "source_group": "wave44_glopl_2025_reproduction",
        "source_provider": provider,
        "source_url": "https://doi.org/10.1038/s41559-025-02784-9",
        "source_record_id": f"wave44:{species}",
        "source_citation": "Lin et al. 2025",
        "source_excerpt": "Species row; self-compatible=1",
        "evidence_scope": (
            "external_congener_species_direct" if external else "synonym_direct"
        ),
        "name_match_method": "strict_wfo_gbif_two_backbone",
        "source_lineage": (
            "glopl-burns-trait-compilation:"
            f"{species.casefold()}:self_incompatibility"
        ),
        "lineage_method": "underlying_glopl_burns_trait_compilation",
        "source_run_id": "europe-pmc:PMC12216600:supplementary-files",
        "source_artifact": "41467_2025_61032_MOESM7_ESM.xlsx",
        "source_file": "41467_2025_61032_MOESM7_ESM.xlsx",
        "acceptance_contract": EXTERNAL_CONTRACT if external else DIRECT_CONTRACT,
    }


def test_wave44_evidence_keeps_direct_and_external_scopes_separate() -> None:
    target = {"Target plant"}
    _validate_evidence(
        pd.DataFrame([_evidence_row(external=False)]),
        target_species=target,
        external=False,
    )
    _validate_evidence(
        pd.DataFrame([_evidence_row(external=True)]),
        target_species=target,
        external=True,
    )


def test_wave44_external_evidence_cannot_enter_direct_scope() -> None:
    row = _evidence_row(external=True)
    row["evidence_scope"] = "synonym_direct"
    with pytest.raises(ValueError, match="external evidence violates"):
        _validate_evidence(
            pd.DataFrame([row]),
            target_species={"Target plant"},
            external=True,
        )
