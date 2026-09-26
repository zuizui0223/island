from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

HISTORICAL_WORKFLOWS = [
    ".github/workflows/audit-chapter1-effect-fingerprint.yml",
    ".github/workflows/audit-chapter1-p1-paired-support.yml",
    ".github/workflows/audit-chapter1-v13-global-only.yml",
    ".github/workflows/audit-chapter1-v13-submission.yml",
    ".github/workflows/promote-chapter1-v11-p3.yml",
    ".github/workflows/promote-chapter1-v11-submission-surface.yml",
    ".github/workflows/promote-chapter1-v13-submission-surface.yml",
    ".github/workflows/render-chapter1-figure3-response-fingerprint.yml",
    ".github/workflows/render-chapter1-figure3-submission.yml",
    ".github/workflows/render-chapter1-v10-figure2-p2.yml",
    ".github/workflows/render-chapter1-v11-figure4.yml",
    ".github/workflows/render-chapter1-v8-figure1.yml",
    ".github/workflows/render-chapter1-v8-figure2.yml",
    ".github/workflows/render-chapter1-v8-figure3-submission.yml",
    ".github/workflows/render-chapter1-v8-figure3.yml",
    ".github/workflows/render-chapter1-v8-figure4.yml",
    ".github/workflows/render-chapter1-v9-figure3-p1-defense.yml",
    ".github/workflows/run-chapter1-all-data-probability.yml",
    ".github/workflows/run-chapter1-h3-conditional-area.yml",
    ".github/workflows/run-chapter1-h3-equal-island.yml",
    ".github/workflows/run-chapter1-h3-grouped-sensitivity.yml",
    ".github/workflows/run-chapter1-h3-support-guardrail.yml",
    ".github/workflows/run-chapter1-h5-channel-state-overlap-audit.yml",
    ".github/workflows/run-chapter1-h5-glopl-floral-architecture-moderation.yml",
    ".github/workflows/run-chapter1-h5-glopl-global-distance.yml",
    ".github/workflows/run-chapter1-h5-glopl-global-shape-audit.yml",
    ".github/workflows/run-chapter1-h5-glopl-influence-audit.yml",
    ".github/workflows/run-chapter1-h5-glopl-island-overlap-preflight.yml",
    ".github/workflows/run-chapter1-h5-glopl-pollen-limitation.yml",
    ".github/workflows/run-chapter1-h5-glopl-reproductive-assurance-moderation.yml",
    ".github/workflows/run-chapter1-h5-identity-matched-two-route.yml",
    ".github/workflows/run-chapter1-h5-multichannel-two-route.yml",
    ".github/workflows/run-chapter1-h5-pollination-triangulation-v2.yml",
    ".github/workflows/run-chapter1-h5-within-island-identity-filter.yml",
    ".github/workflows/run-chapter1-hypothesis-ladder-v2.yml",
    ".github/workflows/run-chapter1-latest-robustness.yml",
    ".github/workflows/run-chapter1-latest-trait-checkpoint.yml",
    ".github/workflows/run-chapter1-response-compression-audit.yml",
    ".github/workflows/run-chapter1-status-resolution-audit.yml",
    ".github/workflows/run-chapter1-trait-first-family-audit.yml",
    ".github/workflows/run-chapter1-v13-functional-bridge.yml",
    ".github/workflows/run-chapter1-v13-raw-colour-audit.yml",
    ".github/workflows/run-chapter1-v14-reordered-hypotheses.yml",
    ".github/workflows/run-chapter1-wcvp-native-compatibility.yml",
    ".github/workflows/validate-chapter1-explanation-gap-contract.yml",
]


def _on_block(text: str) -> str:
    before_permissions = text.split("\npermissions:", 1)[0]
    return before_permissions.split("\non:\n", 1)[1].strip()


def test_superseded_publication_workflows_are_manual_only() -> None:
    for relative in HISTORICAL_WORKFLOWS:
        text = (ROOT / relative).read_text(encoding="utf-8")
        assert text.startswith("name: Historical replay —"), relative
        assert _on_block(text) == "workflow_dispatch:", relative
        assert "\n  push:" not in text.split("\npermissions:", 1)[0], relative
        assert "\n  pull_request:" not in text.split("\npermissions:", 1)[0], relative
        assert "\n  schedule:" not in text.split("\npermissions:", 1)[0], relative
