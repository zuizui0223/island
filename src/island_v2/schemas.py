"""Strict schemas for LLM-generated trait evidence candidates.

These models deliberately represent *candidate evidence*, not adjudicated trait values.
No caller may promote an LLM record directly to the curated database without review.
"""

from __future__ import annotations

from enum import StrEnum
from typing import Literal

from pydantic import BaseModel, Field, HttpUrl, field_validator


class EvidenceScope(StrEnum):
    """Taxonomic resolution of the cited evidence."""

    SPECIES_DIRECT = "species_direct"
    SPECIES_INDIRECT = "species_indirect"
    GENUS_INFERENCE = "genus_inference"
    FAMILY_INFERENCE = "family_inference"
    UNRESOLVED = "unresolved"


class SourceType(StrEnum):
    """Source type, retained separately from inference status."""

    PRIMARY_STUDY = "primary_study"
    FLORA_OR_MONOGRAPH = "flora_or_monograph"
    TAXONOMIC_DATABASE = "taxonomic_database"
    REVIEW = "review"
    CURATED_TRAIT_DATABASE = "curated_trait_database"
    OTHER = "other"
    NONE_FOUND = "none_found"


class Confidence(StrEnum):
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"
    UNRESOLVED = "unresolved"


class TraitEvidenceCandidate(BaseModel):
    """One trait claim with a directly auditable evidence trail."""

    trait_name: str = Field(
        description="Canonical v2 trait name from config/trait_ontology.yml."
    )
    raw_value: str = Field(
        description="The source wording or a faithful short transcription; never a fabricated quotation."
    )
    standardized_value: str = Field(
        description="One allowed ontology value, or unresolved."
    )
    evidence_scope: EvidenceScope
    source_type: SourceType
    source_citation: str = Field(
        description="Short bibliographic identifier, database record name, or none_found."
    )
    source_url: str | None = Field(
        default=None,
        description="Stable source URL only when available; otherwise null.",
    )
    evidence_excerpt: str = Field(
        description="Short supporting excerpt or an explicit statement that no usable evidence was found."
    )
    confidence: Confidence
    inference_note: str = Field(
        description="Why the coding is direct, indirect, inferred, or unresolved."
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
        # Avoid coercing opaque identifiers into URLs while still rejecting malformed URLs.
        HttpUrl(value)
        return value


class SpeciesTraitBatch(BaseModel):
    """The atomic output object from one controlled LLM batch."""

    accepted_species: str
    taxon_id: str | None = None
    trait_candidates: list[TraitEvidenceCandidate]
    batch_notes: str = Field(
        description="Concise notes on ambiguity, conflicts, or missing evidence."
    )
    no_evidence_found: bool = Field(
        description="True only if no trait-level evidence could be recovered for this species."
    )


class TraitExtractionResponse(BaseModel):
    """Schema enforced through OpenAI Structured Outputs."""

    prompt_version: Literal["v2.0"]
    model_interpretation_note: str
    species: list[SpeciesTraitBatch]
