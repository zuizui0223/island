import numpy as np
import pandas as pd
import pytest

shapely = pytest.importorskip("shapely")
from shapely.geometry import box  # noqa: E402

from fd.range_map_nectarivore import (  # noqa: E402
    OTHER,
    PARTIAL,
    SPECIALIST,
    build_species_table,
    crosswalk_names,
    diet_class,
    diet_counts,
    filter_range_records,
    island_incidence,
    normalise_binomial,
)


def islands():
    # three 10x10 islands at x = 0, 100, 200
    return ["a", "b", "c"], [box(0, 0, 10, 10), box(100, 0, 110, 10), box(200, 0, 210, 10)]


def test_normalise_binomial_drops_authority_and_subspecies():
    assert normalise_binomial("Myzomela  rubratra saffordi (Wetmore, 1917)") == "myzomela rubratra"
    assert normalise_binomial("Vestiaria_coccinea") == "vestiaria coccinea"
    with pytest.raises(ValueError):
        normalise_binomial("Myzomela")


def test_filter_range_records_codes():
    r = pd.DataFrame(
        {
            "sci_name": list("abcdef"),
            "presence": [1, 2, 3, 1, 1, 1],
            "origin": [1, 1, 1, 3, 2, 1],
            "seasonal": [1, 2, 1, 1, 1, 3],
        }
    )
    assert filter_range_records(r)["sci_name"].tolist() == ["a", "b", "e"]
    wide = filter_range_records(r, seasonal=frozenset({1, 2, 3, 4}))
    assert wide["sci_name"].tolist() == ["a", "b", "e", "f"]


def test_cover_threshold_and_union_of_polygons():
    ids, geoms = islands()
    species = ["Aa aa", "Aa aa", "Bb bb", "Cc cc"]
    ranges = [
        box(0, 0, 3, 10),  # 30 % of a
        box(3, 0, 6, 10),  # +30 % of a after union -> 60 %
        box(100, 0, 102, 10),  # 20 % of b
        box(-50, -50, 500, 50),  # covers all
    ]
    inc = island_incidence(ids, geoms, species, ranges, min_cover=0.5)
    sp = list(inc.species)
    m = inc.matrix.toarray()
    assert sp == ["aa aa", "bb bb", "cc cc"]
    assert m[:, sp.index("aa aa")].tolist() == [True, False, False]
    assert m[:, sp.index("bb bb")].tolist() == [False, False, False]
    assert m[:, sp.index("cc cc")].tolist() == [True, True, True]
    np.testing.assert_allclose(inc.cover.toarray()[0, 0], 0.6)
    assert inc.at_threshold(0.1).matrix.toarray()[1, sp.index("bb bb")]
    assert inc.at_threshold(0).matrix.toarray()[1, sp.index("bb bb")]


def test_streaming_matches_in_memory_when_species_split_across_chunks():
    from fd.range_map_nectarivore import island_incidence_streaming

    ids, geoms = islands()
    species = ["Aa aa", "Aa aa", "Bb bb", "Cc cc"]
    ranges = [box(0, 0, 6, 10), box(3, 0, 8, 10), box(100, 0, 102, 10), box(-50, -50, 500, 50)]
    full = island_incidence(ids, geoms, species, ranges, min_cover=0.5)
    chunks = [(species[:1], ranges[:1]), (species[1:3], ranges[1:3]), (species[3:], ranges[3:])]
    stream = island_incidence_streaming(ids, geoms, chunks, min_cover=0.5)
    assert stream.species == full.species
    np.testing.assert_allclose(stream.cover.toarray(), full.cover.toarray())
    # overlapping polygons of one species are unioned, not summed: a is 80 % covered, not 110 %
    np.testing.assert_allclose(stream.cover.toarray()[0, 0], 0.8)


def test_nectarivore_clades():
    from fd.range_map_nectarivore import CLADES, nectarivore_clade

    assert CLADES[nectarivore_clade("Trochilidae", "Eulampis")] == "trochilidae"
    assert CLADES[nectarivore_clade("Fringillidae", "Vestiaria")] == "drepanidini"
    assert CLADES[nectarivore_clade("Meliphagidae", "Myzomela")] == "meliphagidae"
    assert CLADES[nectarivore_clade("Psittaculidae", "Vini")] == "other"


def test_buffer_keeps_unbuffered_denominator():
    ids, geoms = islands()
    # range sits 1 unit off the coast of island a
    inc = island_incidence(ids, geoms, ["Aa aa"], [box(11, 0, 30, 10)], min_cover=0.0)
    assert inc.matrix.nnz == 0
    buf = island_incidence(ids, geoms, ["Aa aa"], [box(11, 0, 30, 10)], min_cover=0.0, island_buffer=2)
    assert buf.matrix.toarray()[0, 0]
    assert buf.cover.toarray()[0, 0] <= 1.0


def test_crosswalk_exact_then_synonym_no_fuzzy():
    cw = crosswalk_names(
        ["Aa aa", "Bb bb", "Cc cc", "Cc cx"],
        ["Aa aa", "Bb bz", "Cc cc"],
        pd.DataFrame({"from": ["Bb bb"], "to": ["Bb bz"]}),
    )
    assert cw.mapping == {"aa aa": "aa aa", "bb bb": "bb bz", "cc cc": "cc cc"}
    assert cw.unmatched == ("cc cx",)
    assert cw.via_synonym == 1
    with pytest.raises(ValueError):
        crosswalk_names(["Aa aa"], ["Aa aa"], pd.DataFrame({"from": ["Bb bb", "Bb bb"], "to": ["X x", "Y y"]}))


def test_diet_classes():
    assert diet_class([0, 10, 49.9, 50, 100]).tolist() == [OTHER, PARTIAL, PARTIAL, SPECIALIST, SPECIALIST]
    assert diet_class([20], specialist_cut=20).tolist() == [SPECIALIST]
    with pytest.raises(ValueError):
        diet_class([101])


def _inc_with_species(names):
    ids, geoms = islands()
    return island_incidence(ids, geoms, names, [box(-1, -1, 500, 20)] * len(names), min_cover=0.5)


def test_species_table_fails_closed_on_unmatched_occurring_species():
    inc = _inc_with_species(["Aa aa", "Bb bb", "Cc cc"])
    traits = pd.DataFrame({"scientific": ["Aa aa", "Bb bb"], "diet_nect": [60, 0], "family": ["F1", "F2"]})
    cw = crosswalk_names(inc.species, traits["scientific"])
    with pytest.raises(ValueError, match="unmatched"):
        build_species_table(inc, traits, cw)
    sub, table, report = build_species_table(inc, traits, cw, max_unmatched_share=0.5)
    assert report["unmatched_occurring"] == 1
    assert table.species == ("aa aa", "bb bb")
    assert table.diet.tolist() == [SPECIALIST, OTHER]
    assert sub.matrix.shape == (3, 2)


def test_diet_counts_rows_sum_to_richness():
    inc = _inc_with_species(["Aa aa", "Bb bb", "Cc cc"])
    counts = diet_counts(inc.matrix, np.array([SPECIALIST, PARTIAL, OTHER]))
    assert counts.tolist() == [[1, 1, 1]] * 3
