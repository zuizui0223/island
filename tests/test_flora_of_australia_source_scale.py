import hashlib
from pathlib import Path

import pandas as pd

from analysis.acquire_flora_of_australia_traits_20260809 import (
    AXIS,
    _cultivar_or_hybrid_treatment,
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


def test_cultivar_gate_distinguishes_dimension_cross_from_hybrid_name() -> None:
    description = "Flowers white, 5.5-14 × 2-6 cm; corolla 12 mm long."

    assert not _cultivar_or_hybrid_treatment("Plantus alba", description)
    assert _cultivar_or_hybrid_treatment("Plantus × hybrida", description)
    assert _cultivar_or_hybrid_treatment(
        "Plantus alba", "The flowers described here are from a red cultivar."
    )


def test_bud_form_and_receptacle_tube_do_not_define_open_flower_form() -> None:
    bud = extract_description(
        "Flowers white, slightly urceolate in bud; petals obovate.",
        _rules(),
    )
    receptacle = extract_description(
        "Flowers unisexual; receptacle-tube campanulate, enclosing the ovary.",
        _rules(),
    )

    assert "floral_form" not in bud
    assert "floral_form" not in receptacle


def test_flower_form_rejects_component_and_comparison_object_shapes() -> None:
    descriptions = (
        "Corolla-tube 1.4 cm long; labellum funnel-shaped, white.",
        "Calyx enclosing the corolla tube, broadly campanulate.",
        "Copious panicles of purple campanulate flower-heads.",
        "Corolla with a narrower cylindrical-tubular portion at fruiting-time.",
    )

    for description in descriptions:
        assert "floral_form" not in extract_description(description, _rules())


def test_flower_form_preserves_corolla_transition_across_lobe_phrase() -> None:
    result = extract_description(
        "Corolla campanulate with entire lobes to funnel-shaped with ciliate lobes.",
        _rules(),
    )

    assert result["floral_form"]["value"] == "bell_campanulate|funnel_trumpet"


def test_colour_excludes_nonfloral_drying_and_comparison_contexts() -> None:
    axis = extract_description(
        "Racemes with 2-3 flowers above the straw-coloured glabrescent axis.",
        _rules(),
    )
    comparison = extract_description(
        "Handsome white sweet-scented flowers somewhat resembling the wild rose.",
        _rules(),
    )
    drying = extract_description(
        "White flowers (turning red-brown on drying) in abbreviated spikes. "
        "Petals white, drying red-brown.",
        _rules(),
    )

    assert "flower_primary_color" not in axis
    assert comparison["flower_primary_color"]["value"] == "white"
    assert drying["flower_primary_color"]["value"] == "white"


def test_colour_rejects_supporting_parts_indumentum_and_fruiting_perianth() -> None:
    ancillary = extract_description(
        "Flowers green with crimson stamens. Flowers cream with pink pedicels.",
        _rules(),
    )
    indumentum = extract_description(
        "Female flowers with tips white-ciliate and petals with white hairs. "
        "Corolla white, greyish pubescent outside.",
        _rules(),
    )
    fruiting = extract_description(
        "Achene enclosed in the persistent, red, fleshy, accrescent perianth. "
        "Achene surrounded by the orange, accrescent perianth.",
        _rules(),
    )
    nonfloral_parts = extract_description(
        "Aerial parts deep coral pink, perianth paler, subterranean parts whitish.",
        _rules(),
    )

    assert ancillary["flower_primary_color"]["value"] == (
        "green_brown_inconspicuous|white"
    )
    assert indumentum["flower_primary_color"]["value"] == "white"
    assert "flower_primary_color" not in fruiting
    assert "flower_primary_color" not in nonfloral_parts


def test_colour_excludes_pressed_dry_state_colours() -> None:
    result = extract_description(
        "Perianth segments pale greenish to white, in the dry state pale "
        "brownish to blackish.",
        _rules(),
    )

    assert result["flower_primary_color"]["value"] == (
        "green_brown_inconspicuous|white"
    )


def test_colour_rejects_bract_hair_pollen_and_oil_dot_colours() -> None:
    descriptions = (
        "Male spikes are subtended by three large white leafy bracts.",
        "Stellate hairs on the outer surface of the perianth cream or golden.",
        "Petals orbicular, with inconspicuous oil dots.",
        "Flowers fragrant; pollen white; staminal filaments fleshy.",
    )

    for description in descriptions:
        assert "flower_primary_color" not in extract_description(description, _rules())


def test_numeric_two_lipped_form_and_regular_symmetry_are_preserved() -> None:
    form = extract_description(
        "Corolla tubular below and 2-lipped at the apex.",
        _rules(),
    )
    symmetry = extract_description(
        "Flowers regular or slightly zygomorphic.",
        _rules(),
    )

    assert form["floral_form"]["value"] == "bilabiate|tubular"
    assert symmetry["floral_symmetry"]["value"] == "actinomorphic|zygomorphic"


def test_colour_keeps_coordinated_and_distant_flower_states() -> None:
    coordinated = extract_description(
        "Flowers reddish and green with yellow anthers.",
        _rules(),
    )
    long_perianth = extract_description(
        "Perianth-tube greenish, 9-15 cm long, curved above the erect involucre; "
        "segments green in bud, but during anthesis pure white with a broad "
        "sharply bordered purple dorsal streak.",
        _rules(),
    )

    assert coordinated["flower_primary_color"]["value"] == (
        "green_brown_inconspicuous|red_pink"
    )
    assert long_perianth["flower_primary_color"]["value"] == (
        "blue_purple|green_brown_inconspicuous|white"
    )


def test_flower_size_rejects_inflorescence_and_reproductive_part_measurements() -> None:
    descriptions = (
        "Flowers 2-6, pseudo-umbellate or with rhachis up to 4 mm long.",
        "Flowers forming apparent verticillasters 1-15 cm long.",
        "Style 5 mm long in female flowers with capitate stigma 1 mm long.",
        (
            "Flowers sweetly scented; sepals rounded; tube 17-23 mm long; "
            "lobes 19-30 mm long; corolla mouth open; disk 2-4.5 mm high."
        ),
        "Catkins of flowers 3-3.5 cm long.",
        "Long-styled flowers, stigma arms 2-2.5 mm long.",
        "Flowers in dense terminal cymes 1.5-4.5 cm wide.",
        (
            "Anthers present in both flower forms, 0.75 mm long in functionally "
            "male flowers, 0.5 mm long in functionally female flowers."
        ),
        (
            "Filaments 2 mm long in long-styled flowers, 1.5 mm long in "
            "short-styled flowers."
        ),
        "Flowers fragrant; iInflorescence about 40 cm long; stamens 18-22.",
    )

    for description in descriptions:
        assert "flower_size_class" not in extract_description(description, _rules())


def test_tube_depth_does_not_use_overall_corolla_length() -> None:
    result = extract_description(
        "Petals fused to form a short corolla tube, corolla about 3 mm long overall.",
        _rules(),
    )

    assert "tube_depth_class" not in result


def test_display_preserves_solitary_and_grouped_multistates() -> None:
    cyme = extract_description(
        "Flowers terminal in few-flowered cymes or solitary.",
        _rules(),
    )
    fascicle = extract_description(
        "Flowers axillary, solitary or 2-3 in a fascicle or short raceme.",
        _rules(),
    )

    assert cyme["inflorescence_display"]["value"] == (
        "few_flowered|solitary|umbel_corymb"
    )
    assert fascicle["inflorescence_display"]["value"] == (
        "few_flowered|raceme_spike_panicle|solitary"
    )


def test_display_handles_inflected_terms_and_rejects_incomplete_or_fruit_context() -> None:
    inflected = extract_description(
        "Inflorescence racemose, spicate and few-flowered.",
        _rules(),
    )
    unsupported_fascicle = extract_description(
        "Flowers in dense axillary fascicles, sometimes flowers solitary.",
        _rules(),
    )
    fruit = extract_description(
        "Achenes 45-100 per capitulum, 1.2-3 mm long.",
        _rules(),
    )
    grouped_solitary = extract_description(
        "Flowers in dense cymes forming corymbs either terminal and solitary "
        "or grouped in panicles.",
        _rules(),
    )

    assert inflected["inflorescence_display"]["value"] == (
        "few_flowered|raceme_spike_panicle"
    )
    assert "inflorescence_display" not in unsupported_fascicle
    assert "inflorescence_display" not in fruit
    assert grouped_solitary["inflorescence_display"]["value"] == (
        "raceme_spike_panicle|umbel_corymb"
    )


def test_display_preserves_heads_clusters_long_solo_and_corymbose_states() -> None:
    heads = extract_description(
        "White florets in small heads in terminal panicles.",
        _rules(),
    )
    clusters = extract_description(
        "Flowers fragrant, solitary or in clusters of three together.",
        _rules(),
    )
    long_solo = extract_description(
        "Flowers greenish-white, pedicels up to 2 mm long, in short "
        "few-flowered lateral racemes up to about 1 cm long or solitary.",
        _rules(),
    )
    corymbose = extract_description(
        "Racemes many-flowered, pedunculate and corymbose.",
        _rules(),
    )
    paniculate = extract_description(
        "The yellow florets are held in capitula in paniculate inflorescences.",
        _rules(),
    )

    assert heads["inflorescence_display"]["value"] == (
        "composite_display|raceme_spike_panicle"
    )
    assert clusters["inflorescence_display"]["value"] == "few_flowered|solitary"
    assert long_solo["inflorescence_display"]["value"] == (
        "few_flowered|raceme_spike_panicle|solitary"
    )
    assert corymbose["inflorescence_display"]["value"] == (
        "raceme_spike_panicle|umbel_corymb"
    )
    assert paniculate["inflorescence_display"]["value"] == (
        "composite_display|raceme_spike_panicle"
    )


def test_display_rejects_fruit_leaf_and_stamen_arrangement_contexts() -> None:
    descriptions = (
        "Figs borne on long panicles on the main stem.",
        "Leaves reaching at least to the base of the spike.",
        "Stem 4 mm in diameter at the base of the spike.",
        (
            "Male inflorescences with spadices; stamens racemose, the individual "
            "racemes 1 cm long."
        ),
    )

    for description in descriptions:
        assert "inflorescence_display" not in extract_description(description, _rules())


def test_display_keeps_capitate_bracteate_head_and_small_cluster_states() -> None:
    capitate = extract_description(
        "Inflorescences capitate, few-flowered or many-flowered.",
        _rules(),
    )
    heads = extract_description(
        "Flowers in few-flowered bracteate heads.",
        _rules(),
    )
    clusters = extract_description(
        "Male flowers solitary or in small axillary clusters.",
        _rules(),
    )
    numbered = extract_description(
        "Inflorescence 2-flowered; flowers sometimes solitary.",
        _rules(),
    )
    variable = extract_description(
        "Flowers in few- to many-flowered corymbose cymes.",
        _rules(),
    )

    assert capitate["inflorescence_display"]["value"] == (
        "composite_display|few_flowered"
    )
    assert heads["inflorescence_display"]["value"] == (
        "composite_display|few_flowered"
    )
    assert clusters["inflorescence_display"]["value"] == "few_flowered|solitary"
    assert numbered["inflorescence_display"]["value"] == "few_flowered|solitary"
    assert variable["inflorescence_display"]["value"] == (
        "few_flowered|umbel_corymb"
    )


def test_display_rejects_large_cluster_range_and_comparison_only_panicle() -> None:
    large_cluster = extract_description(
        "Inflorescence in clusters of 3-12 flowers.",
        _rules(),
    )
    comparison = extract_description(
        "The inflorescence simulates a panicle and is commonly mistaken for Panicum.",
        _rules(),
    )
    capitate_hairs = extract_description(
        "A perennial 150 cm high including the inflorescence; hairs capitate-glandular.",
        _rules(),
    )

    assert "inflorescence_display" not in large_cluster
    assert "inflorescence_display" not in comparison
    assert "inflorescence_display" not in capitate_hairs


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
