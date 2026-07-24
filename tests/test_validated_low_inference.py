import pandas as pd

from island_v2.validated_low_inference import Thresholds, infer_low_evidence


def evidence(rows):
    return pd.DataFrame(rows, columns=[
        "accepted_species", "trait_name", "normalized_value",
        "evidence_quality", "evidence_scope",
    ])


def test_accepts_consistent_structural_genus_rule():
    direct = evidence([
        ("Alpha one", "floral_symmetry", "actinomorphic", "high", "species_direct"),
        ("Alpha two", "floral_symmetry", "actinomorphic", "medium", "species_direct"),
        ("Alpha three", "floral_symmetry", "actinomorphic", "high", "synonym_direct"),
        ("Alpha four", "floral_symmetry", "actinomorphic", "high", "species_direct"),
    ])
    tasks = pd.DataFrame([{"accepted_species": "Alpha target", "axis": "floral_structural_complexity"}])
    accepted, remaining, rules = infer_low_evidence(
        tasks, direct, Thresholds(min_species=3, min_masked_accuracy=0.75)
    )
    assert len(accepted) == 1
    assert accepted.iloc[0]["normalized_value"] == "actinomorphic"
    assert accepted.iloc[0]["evidence_scope"] == "genus_consensus"
    assert remaining.empty
    assert bool(rules.iloc[0]["eligible"])


def test_rejects_variable_flower_colour():
    direct = evidence([
        ("Beta one", "flower_color", "red", "high", "species_direct"),
        ("Beta two", "flower_color", "white", "high", "species_direct"),
        ("Beta three", "flower_color", "red", "high", "species_direct"),
        ("Beta four", "flower_color", "white", "high", "species_direct"),
    ])
    tasks = pd.DataFrame([{"accepted_species": "Beta target", "axis": "flower_colour"}])
    accepted, remaining, rules = infer_low_evidence(tasks, direct)
    assert accepted.empty
    assert len(remaining) == 1
    assert not bool(rules.iloc[0]["eligible"])


def test_never_uses_family_or_global_fallback():
    direct = evidence([
        ("Gamma one", "mating_system", "self_compatible", "high", "family_inference"),
        ("Gamma two", "mating_system", "self_compatible", "medium", "global_fallback"),
        ("Gamma three", "mating_system", "self_compatible", "high", "species_direct"),
    ])
    tasks = pd.DataFrame([{"accepted_species": "Gamma target", "axis": "reproductive_assurance"}])
    accepted, remaining, rules = infer_low_evidence(tasks, direct)
    assert accepted.empty
    assert len(remaining) == 1
    assert rules.iloc[0]["n_direct_species"] == 1


def test_reads_real_acquisition_raw_candidate_value():
    direct = pd.DataFrame([
        ("Delta one", "floral_form", "tubular", "high", "species_direct"),
        ("Delta two", "floral_form", "tubular", "medium", "species_direct"),
        ("Delta three", "floral_form", "tubular", "high", "synonym_direct"),
        ("Delta four", "floral_form", "tubular", "high", "species_direct"),
    ], columns=[
        "accepted_species", "trait_name", "raw_candidate_value",
        "evidence_quality", "evidence_scope",
    ])
    tasks = pd.DataFrame([{"accepted_species": "Delta target", "axis": "floral_structural_complexity"}])
    accepted, remaining, rules = infer_low_evidence(tasks, direct)
    assert len(accepted) == 1
    assert accepted.iloc[0]["normalized_value"] == "tubular"
    assert remaining.empty
    assert bool(rules.iloc[0]["eligible"])


def test_empty_usable_evidence_returns_schema_stable_outputs():
    direct = pd.DataFrame([{
        "accepted_species": "Epsilon one",
        "trait_name": "floral_form",
        "raw_candidate_value": "tubular",
        "evidence_quality": "low",
        "evidence_scope": "family_inference",
    }])
    tasks = pd.DataFrame([{"accepted_species": "Epsilon target", "axis": "floral_structural_complexity"}])
    accepted, remaining, rules = infer_low_evidence(tasks, direct)
    assert accepted.empty
    assert list(accepted.columns)
    assert len(remaining) == 1
    assert remaining.iloc[0]["reason"] == "no_usable_direct_evidence"
    assert list(rules.columns)
