"""Shared controlled vocabularies for v2 trait evidence tables.

These enums label evidence provenance and taxonomic scope. They are used by
bulk imports and reviewable evidence workflows.
"""

from __future__ import annotations

from enum import StrEnum


class EvidenceScope(StrEnum):
    """Taxonomic resolution of the biological statement."""

    SPECIES_DIRECT = "species_direct"
    SPECIES_INDIRECT = "species_indirect"
    GENUS_INFERENCE = "genus_inference"
    FAMILY_INFERENCE = "family_inference"
    UNRESOLVED = "unresolved"


class SourceType(StrEnum):
    """Where the candidate came from; source quality is stored separately."""

    PRIMARY_STUDY = "primary_study"
    FLORA_OR_MONOGRAPH = "flora_or_monograph"
    TAXONOMIC_DATABASE = "taxonomic_database"
    REVIEW = "review"
    CURATED_TRAIT_DATABASE = "curated_trait_database"
    INSTITUTIONAL_WEB = "institutional_web"
    CURATED_SPECIALIST_WEB = "curated_specialist_web"
    COMMUNITY_OR_UNVETTED_WEB = "community_or_unvetted_web"
    NONE_FOUND = "none_found"


class SourceReliability(StrEnum):
    """Reliability of the particular source, not confidence in the trait itself."""

    A_PRIMARY_OR_MONOGRAPH = "A_primary_or_monograph"
    B_CURATED_DATABASE_OR_INSTITUTION = "B_curated_database_or_institution"
    C_CURATED_SPECIALIST_WEB = "C_curated_specialist_web"
    D_UNVETTED_WEB = "D_unvetted_web"
    NONE = "none"


class CandidateKind(StrEnum):
    """Whether the record reports a source statement or a declared inference."""

    SOURCE_BACKED = "source_backed"
    HIERARCHICAL_INFERENCE = "hierarchical_inference"
    UNRESOLVED = "unresolved"


class Confidence(StrEnum):
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"
    UNRESOLVED = "unresolved"
