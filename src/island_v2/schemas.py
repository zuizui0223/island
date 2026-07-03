"""Schemas for LLM-generated v2 trait candidates.

The pipeline is designed for broad web retrieval *and* controlled taxonomic
inference. Each record is a reviewable candidate, never an automatically
adjudicated biological fact.
"""

from __future__ import annotations

from enum import StrEnum
from typing import Literal

from pydantic import BaseModel, Field, HttpUrl, field_validator


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
    """Whether the record reports a source statement or a declared imputation."""

    SOURCE_BACKED = "source_backed"
    HIERARCHICAL_INFERENCE = "hierarchical_inference"
    UNRESOLVED = "unresolved"


class Confidence(StrEnum):
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"
    UNRESOLVED = "unresolved"


class TraitEvidenceCandidate(BaseModel):
    """One trait value candidate with an auditable evidence or inference trail."""

    trait_name: str = Field(
        description="Canonical v2 trait name from config/trait_ontology.yml."
    )
    raw_value: str = Field(
        description="Source wording or a faithful short paraphrase; never a fabricated quotation."
    )
    standardized_value: str = Field(
        description="One allowed ontology value, or unresolved."
    )
    candidate_kind: CandidateKind
    evidence_scope: EvidenceScope
    source_type: SourceType
    source_reliability: SourceReliability
    source_citation: str = Field(
        description="Searchable bibliographic title, database record name, or none_found."
    )
    source_url: str | None = Field(
        default=None,
        description="Stable URL retrieved during this run when available; otherwise null.",
    )
    evidence_excerpt: str = Field(
        description="Short faithful source excerpt/paraphrase, or explicit statement of missing evidence."
    )
    supporting_taxa: list[str] = Field(
        default_factory=list,
        description="Congeners/family taxa supporting an inference; empty for direct evidence."
    )
    inference_rule: str | None = Field(
        default=None,
        description="Declared reason for genus/family inference; null for source-backed evidence."
    )
    confidence: Confidence
    inference_note: str = Field(
        description="Why the coding is direct, web-derived, inferred, conflicting, or unresolved."
    )
    needs_human_review: bool = Field(
        description="Must be true for every candidate produced by the LLM pipeline."
    )

    @field_validator("needs_human_review")
    @classmethod
    def require_review(cls, value: bool) -> bool:
        if value is not True:
            raise ValueError("LLM candidates must always require human review.")
        return value

    @field_validator("source_url")
    @classmethod
    def validate_url(cls, value: str | None) -> str | None:
        if value is None or value == "":
            return None
        HttpUrl(value)
        return value


class SpeciesTraitBatch(BaseModel):
    """The atomic output object from one controlled web-search batch."""

    accepted_species: str
    taxon_id: str | None = None
    trait_candidates: list[TraitEvidenceCandidate]
    batch_notes: str = Field(
        description="Ambiguities, conflicts, missing traits, and recommended review targets."
    )
    no_evidence_found: bool = Field(
        description="True only if no source-backed or declared hierarchical candidate could be recovered."
    )


class TraitExtractionResponse(BaseModel):
    """Schema enforced through OpenAI Structured Outputs."""

    prompt_version: Literal["v2.1"]
    model_interpretation_note: str
    species: list[SpeciesTraitBatch]
