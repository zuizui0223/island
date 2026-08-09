import hashlib
from pathlib import Path

import pandas as pd

from analysis.acquire_flora_of_australia_traits_20260809 import (
    AXIS,
    _source_lineage,
    deterministic_audit_queue,
    extract_description,
)
from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS
from island_v2.open_web_common import reviewed_source_package_evidence
from island_v2.trait_measurements import load_rules

COMMITTED_PACKAGE = (
    Path("data/v2/staging/traits/direct_llm_pilot")
    / "20260809_flora_of_australia_source_acquisition"
)


def _rules():
    return load_rules(Path("config/measurement_classification.yml"))


def test_six_traits_are_extracted_from_one_species_description() -> None:
    description = (
        "Flowers solitary, zygomorphic, tubular, 12-18 mm long, red and yellow; "
        "corolla tube 8-10 mm long. Leaves green. Inflorescences in racemes."
    )
    result = extract_description(description, _rules())

    assert result["flower_primary_color"]["value"] == "red_pink|yellow_orange"
    assert result["floral_symmetry"]["value"] == "zygomorphic"
    assert result["floral_form"]["value"] == "tubular"
    assert result["flower_size_class"]["value"] == "medium"
    assert result["tube_depth_class"]["value"] == "intermediate"
    assert result["inflorescence_display"]["value"] == (
        "raceme_spike_panicle|solitary"
    )
    assert "Leaves green" not in result["flower_primary_color"]["quote"]


def test_non_target_calyx_and_leaf_descriptions_fail_closed() -> None:
    description = (
        "Calyx tubular and green, 8 mm long. Leaves red beneath. "
        "Petals white; flowers 4 mm across."
    )
    result = extract_description(description, _rules())

    assert "floral_form" not in result
    assert result["flower_primary_color"]["value"] == "white"
    assert result["flower_size_class"]["value"] == "small"
    assert "tube_depth_class" not in result


def test_botanical_circa_abbreviation_does_not_detach_subject() -> None:
    result = extract_description(
        "Petals obovate, c. 15 mm long, rounded apically, yellow.",
        _rules(),
    )

    assert result["flower_primary_color"]["value"] == "yellow_orange"


def test_numeric_new_sentence_does_not_borrow_previous_flower_subject() -> None:
    result = extract_description(
        "Flowers white. 1.5-2.6 mm long, consisting of three fruitlets.",
        _rules(),
    )

    assert result["flower_primary_color"]["value"] == "white"
    assert "flower_size_class" not in result


def test_decimal_measurement_is_not_split_as_a_new_sentence() -> None:
    result = extract_description("Corolla 1.5-2.5 mm long, white.", _rules())

    assert result["flower_size_class"]["value"] == "very_small"
    assert result["flower_primary_color"]["value"] == "white"


def test_comparison_object_does_not_steal_postcomma_flower_value() -> None:
    result = extract_description(
        "Petals shorter than calyx, yellow. "
        "Corolla tube longer than sepals, light green, urceolate, 4-8 mm long.",
        _rules(),
    )

    assert result["flower_primary_color"]["value"] == (
        "green_brown_inconspicuous|yellow_orange"
    )
    assert result["floral_form"]["value"] == "urn_urceolate"
    assert result["tube_depth_class"]["value"] == "intermediate"


def test_explicit_spike_is_display_and_its_colour_is_floral() -> None:
    result = extract_description(
        "Leaves green. Spikes 2.5-4.5 cm long, golden. Flowers 5-merous.",
        _rules(),
    )

    assert result["inflorescence_display"]["value"] == "raceme_spike_panicle"
    assert "flower_primary_color" not in result


def test_leaf_asymmetry_and_petal_base_form_are_not_whole_flower_traits() -> None:
    result = extract_description(
        "Leaves asymmetric at base. Petals linear, with a tubular base c. 6 mm long.",
        _rules(),
    )

    assert "floral_symmetry" not in result
    assert "floral_form" not in result


def test_negated_form_flower_head_and_lobe_symmetry_fail_closed() -> None:
    result = extract_description(
        "Flowers not campanulate. Flower heads narrowly campanulate. "
        "Corolla lobes asymmetric. Plants usually drying orange.",
        _rules(),
    )

    assert "floral_form" not in result
    assert "floral_symmetry" not in result
    assert "flower_primary_color" not in result


def test_infructescence_and_comparison_measurements_are_excluded() -> None:
    result = extract_description(
        "Infructescence a raceme 20 cm long. "
        "Pedicels shorter than corolla tube, 1 mm long. "
        "Corolla tube 7-12 mm long.",
        _rules(),
    )

    assert "inflorescence_display" not in result
    assert result["tube_depth_class"]["value"] == "intermediate"


def test_comparison_does_not_transfer_scale_measurement_to_corolla_tube() -> None:
    result = extract_description(
        "Scales equal to or exceeding the corolla tube, abundantly fimbriate, "
        "fimbriae to 0.5 mm long.",
        _rules(),
    )

    assert "tube_depth_class" not in result


def test_solitary_within_bract_is_not_a_solitary_flower_display() -> None:
    result = extract_description(
        "Inflorescences many-flowered spikes. "
        "Flowers solitary within each bract, white.",
        _rules(),
    )

    assert result["inflorescence_display"]["value"] == "raceme_spike_panicle"


def test_multi_flowered_inflorescence_does_not_become_solitary() -> None:
    result = extract_description(
        "Inflorescence diffuse, 2-5(-10)-flowered; flowers solitary or loosely clustered.",
        _rules(),
    )

    assert "inflorescence_display" not in result


def test_object_flower_and_nonfloral_colour_nouns_fail_closed() -> None:
    result = extract_description(
        "Inflorescence bracts overtopping flowers, pinkish green. "
        "Corolla pink to orange. Male flowers with an apical fringe of white hairs.",
        _rules(),
    )

    assert result["flower_primary_color"]["value"] == "red_pink|yellow_orange"


def test_umbel_ray_indumentum_is_not_flower_colour() -> None:
    result = extract_description(
        "Central flowers sessile; rays of umbel densely white-tomentose.",
        _rules(),
    )

    assert "flower_primary_color" not in result


def test_measurement_owner_is_not_transferred_to_flower() -> None:
    result = extract_description(
        "Inflorescence peduncles erect in flower stage, 1.3-3.3 cm long. "
        "Bracteoles shorter than the flower, 1-1.3 mm long. "
        "Male flowers: androecium 1.8-2 mm long. "
        "Flower heads depressed-globular, 2.8-4 mm long. "
        "Male flowers with spathe 10-20 mm long. "
        "Corolla white; operculum 8-10 mm long. "
        "Stamens inserted near corolla mouth; filaments 1-3 mm long. "
        "Female flowers present. Capsule 18-30 mm long. "
        "Capsule without calyx and corolla 3-5 mm long.",
        _rules(),
    )

    assert "flower_size_class" not in result


def test_flower_size_rejects_numbered_next_sentence_and_nested_display_parts() -> None:
    result = extract_description(
        "Flowers up to c. 8 mm across, tepals of male flowers 4, of female "
        "flowers 5. Fruits c. 2 cm across. "
        "Inflorescence has widely spaced flowers, fertile portion 25 mm long. "
        "Flowers borne on pedicels of open flowers 4-7 mm long. "
        "Supplementary inflorescences have 1-2 flowers on an axis 2.5 cm long. "
        "Pedicels in flower 2.6 cm long. "
        "Racemes dense in flower, 3 cm long in fruit. "
        "Female flowers on 8.5 cm long peduncles. "
        "Green-flower tree 5 m high. "
        "Flowers on a wiry scape 400 mm high. "
        "Spikes with 3-5 flowers, 9 mm long in fruit. "
        "Flowers with lower glumes and lemmas 12 mm long. "
        "Flowers with heads 12 mm long. "
        "Upper lip of corolla 14 mm long. "
        "Corolla blue; limb 15 mm across. "
        "Calyx tube in female flowers 4 mm long. "
        "Capsule enclosed by corolla remains 4 mm long. "
        "Irregular flowers borne on 8 cm long unbranched peduncles. "
        "Three bracts to a flower, 5 mm long. "
        "Perianth of six bristles 3 mm long. "
        "Upper portion of perianth 15 mm long. "
        "Panicles with clusters of flowers; female 180 mm long. "
        "Fruit clearly exceeding perianth, 3 mm long. "
        "Flowers in umbels; buds 12 mm long. "
        "Perianth 5 cm long; perianth claws 8 mm long. "
        "Cymes with bracts looking as if flowers solitary, bracts 3 mm long. "
        "Corolla yellow; standard 14 mm long. "
        "Flowers on short, 26 mm long, straight peduncles. "
        "Female flowers with staminodes on a gynophore 0.3 mm long. "
        "Perianth 4-5 cm long; perianth limbs 5-6 mm long. "
        "Flowers in terminal flower heads, up to 200 mm long.",
        _rules(),
    )

    assert result["flower_size_class"]["value"] == "large|small"
    assert "medium" not in result["flower_size_class"]["value"]
    assert "very_large" not in result["flower_size_class"]["value"]


def test_cross_semicolon_corolla_tube_is_depth_not_whole_flower_size() -> None:
    result = extract_description(
        "Corolla yellow; tube c. 2 mm long, bearded at throat.",
        _rules(),
    )

    assert result["tube_depth_class"]["value"] == "shallow"
    assert "flower_size_class" not in result


def test_parenthetical_tube_range_uses_core_range_and_rejects_lip_length() -> None:
    ranged = extract_description(
        "Corolla 20-25(-30) mm long; tube 12-16(-18) mm long, with basal "
        "inflation 4 mm diam.",
        _rules(),
    )
    lip_only = extract_description(
        "Corolla 7 mm long; tube straight; upper lip 1 mm long.",
        _rules(),
    )

    assert ranged["tube_depth_class"]["value"] == "intermediate"
    assert ranged["tube_depth_class"]["raw_measurements"] == ("12-16 mm",)
    assert "tube_depth_class" not in lip_only


def test_source_lineage_uses_original_source_attribute_not_profile_url() -> None:
    first = _source_lineage("Adapted from Smith (1984).", "profile-1")
    second = _source_lineage("  Adapted FROM Smith (1984) ", "profile-2")
    fallback = _source_lineage("", "profile-3")

    assert first == second
    assert first[1] == "normalized_source_attribute_fingerprint"
    assert fallback == (
        "flora-of-australia-profile:profile-3",
        "official_profile_uuid_fallback",
    )


def test_audit_is_deterministic_and_represents_every_trait() -> None:
    rows = []
    for trait, axis in AXIS.items():
        for index in range(40):
            row = {column: "x" for column in EVIDENCE_COLUMNS}
            row.update(
                {
                    "accepted_species": f"Genus species{index}",
                    "axis": axis,
                    "trait_name": trait,
                    "normalized_value": "white",
                    "source_record_id": f"{trait}-{index}",
                    "source_provider": "flora_of_australia_official",
                    "source_excerpt": "Flowers white.",
                }
            )
            rows.append(row)
    evidence = pd.DataFrame(rows)

    first = deterministic_audit_queue(evidence)
    second = deterministic_audit_queue(evidence.sample(frac=1, random_state=7))

    assert len(first) == 200
    assert first["candidate_id"].tolist() == second["candidate_id"].tolist()
    assert first.groupby("trait_name").size().min() >= 20
    assert first["accepted_correct"].eq("").all()


def test_committed_source_package_passes_manual_audit_gate() -> None:
    evidence = pd.read_csv(
        COMMITTED_PACKAGE / "flora_of_australia_evidence.csv.gz",
        dtype=str,
    ).fillna("")
    audit = pd.read_csv(
        COMMITTED_PACKAGE / "flora_of_australia_manual_audit_200.csv",
        dtype=str,
    ).fillna("")

    selected, scopes, summary = reviewed_source_package_evidence(evidence, audit)

    assert summary["reviewed"] == 200
    assert summary["accepted_correct"] == 200
    assert summary["precision"] == 1.0
    assert summary["cultivar_contamination_rate"] == 0.0
    assert summary["package_gate_passed"] is True
    assert summary["selected_evidence_rows"] == 4_060
    assert set(scopes.loc[scopes["production_approved"], "trait_name"]) == set(AXIS)
    assert len(selected) == 4_060


def test_committed_source_package_manifest_is_exact() -> None:
    entries = {}
    for line in (COMMITTED_PACKAGE / "file_manifest.sha256").read_text(
        encoding="ascii"
    ).splitlines():
        digest, filename = line.split("  ", maxsplit=1)
        entries[filename] = digest

    files = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in COMMITTED_PACKAGE.iterdir()
        if path.is_file() and path.name != "file_manifest.sha256"
    }
    assert entries == files
