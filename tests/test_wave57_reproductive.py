"""Regression tests for the Wave56 genus generalization and unreviewed grades."""
import importlib.util
from pathlib import Path
import unittest

P = Path(__file__).parents[1] / "scripts" / "acquire_wave57_reproductive.py"
spec = importlib.util.spec_from_file_location("wave57", P)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class Gates(unittest.TestCase):
    def test_genus_hybrid_and_cultivar_rejected(self):
        for name in ["Cymbalaria", "Ixora sp", "Dombeya acutangula x delislei", "Lindernia Sudley"]:
            self.assertFalse(m.is_species(name))
        self.assertTrue(m.is_species("Cymbalaria muralis"))

    def test_actual_genus_generalization_fails_closed(self):
        rows = [{"accepted_species": "Cymbalaria", "trait_name": "autonomous_selfing_capacity", "evidence_quality": "high", "excerpt": "Urbanization selects autogamous species."}]
        audit = m.audit_prior(rows, {"Cymbalaria"})
        self.assertEqual(audit[0]["decision"], "invalid_species_name")
        self.assertEqual(audit[0]["promotion_allowed"], "false")

    def test_off_axis_separated(self):
        row = {"accepted_species": "Pilea pumila", "trait_name": "flower_primary_color"}
        self.assertEqual(m.audit_prior([row], {"Pilea pumila"})[0]["decision"], "off_reproductive_axis")

    def test_document_cooccurrence_never_high_or_promoted(self):
        paper = {"text": "Cymbalaria muralis was studied. Other plants are self-compatible.", "lineage": "doi:example", "url": "https://example.org", "provider": "test"}
        row = m.paper_leads(paper, {"Cymbalaria muralis"}, {})[0]
        self.assertEqual(row["evidence_quality"], "unreviewed")
        self.assertEqual(row["normalized_value"], "")
        self.assertEqual(row["species_claim_verified"], "false")
        self.assertEqual(row["promotion_allowed"], "false")

    def test_synonym_preserves_frozen_universe_without_rule_vote(self):
        paper = {"text": "Lindernia micrantha is self-compatible.", "lineage": "doi:example", "url": "https://example.org", "provider": "test"}
        rows = m.paper_leads(paper, {"Vandellia micrantha"}, {})
        self.assertEqual(rows[0]["accepted_species"], "Vandellia micrantha")
        self.assertEqual(rows[0]["genus_rule_training_allowed"], "false")

    def test_absent_breeding_context_or_scope_yields_nothing(self):
        p = {"text": "Pilea pumila has green leaves.", "lineage": "x", "url": "https://example.org", "provider": "test"}
        self.assertEqual(m.paper_leads(p, {"Pilea pumila"}, {}), [])
        p["text"] = "Pilea pumila is self-compatible."
        self.assertEqual(m.paper_leads(p, {"Pilea fontana"}, {}), [])

    def test_si_negation_is_not_normalized(self):
        p = {"text": "Durio graveolens is not self-compatible.", "lineage": "x", "url": "https://example.org", "provider": "test"}
        row = m.paper_leads(p, {"Durio graveolens"}, {})[0]
        self.assertEqual(row["normalized_value"], "")
        self.assertEqual(row["conflicted_genus_rule_blocked"], "true")


if __name__ == "__main__":
    unittest.main()
