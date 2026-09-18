import pytest

from island_v2.chapter1_v14_submission_promotion import promote


def test_promote_marks_verified_result_canonical():
    preflight = {
        "contract": "chapter1_v14_preflight_result_lock_v1",
        "status": "local_reproduction_pending_ci_promotion",
        "canonical": False,
    }
    result = {"contract": "chapter1_v14_reordered_hypotheses_result_v1"}
    verification = {"verified": True}
    promoted = promote(
        result,
        preflight,
        verification,
        ci_run_id=123,
        ci_artifact_id=456,
        ci_artifact_digest="sha256:" + "a" * 64,
        head_sha="b" * 40,
    )
    assert promoted["canonical"] is True
    assert promoted["status"] == "canonical_v14_reproduced"
    assert promoted["promotion_provenance"]["ci_run_id"] == 123
    assert promoted["promotion_provenance"]["ci_artifact_id"] == 456


def test_promote_rejects_unverified_result():
    with pytest.raises(ValueError, match="unverified"):
        promote(
            {"contract": "chapter1_v14_reordered_hypotheses_result_v1"},
            {"canonical": False},
            {"verified": False},
            ci_run_id=1,
            ci_artifact_id=2,
            ci_artifact_digest="sha256:" + "a" * 64,
            head_sha="b" * 40,
        )


def test_promote_requires_full_sha_and_digest():
    with pytest.raises(ValueError, match="sha256"):
        promote(
            {"contract": "chapter1_v14_reordered_hypotheses_result_v1"},
            {"canonical": False},
            {"verified": True},
            ci_run_id=1,
            ci_artifact_id=2,
            ci_artifact_digest="not-a-digest",
            head_sha="b" * 40,
        )
    with pytest.raises(ValueError, match="40-character"):
        promote(
            {"contract": "chapter1_v14_reordered_hypotheses_result_v1"},
            {"canonical": False},
            {"verified": True},
            ci_run_id=1,
            ci_artifact_id=2,
            ci_artifact_digest="sha256:" + "a" * 64,
            head_sha="short",
        )
