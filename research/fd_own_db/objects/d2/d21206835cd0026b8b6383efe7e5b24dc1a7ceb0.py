from fd.dopi_real_motif_multiresolution import (
    make_real_multiresolution_motifs,botanical_axis_coverage,
)


def inputs():
    return {
        "flower_primary_color":{
            "Alpha one":"white","Alpha two":"red_pink","Beta three":"yellow_orange",
        },
        "floral_symmetry":{
            "Alpha one":"actinomorphic","Alpha two":"zygomorphic",
            "Beta three":"actinomorphic","Gamma no_color":"actinomorphic",
        },
        "floral_form":{
            "Alpha one":"tubular","Alpha two":"open_radial",
            "Beta three":"bell_campanulate","Gamma no_color":"open_radial",
        },
        "tube_depth_class":{
            "Alpha one":"deep","Alpha two":"shallow",
            "Beta three":"deep","Gamma no_color":"shallow",
        },
    }


def test_all_27_genuine_motifs_without_global_complete_carrier_filter():
    z=make_real_multiresolution_motifs(inputs())
    assert len(z)==27
    assert botanical_axis_coverage(z)=={1:4,2:11,3:9,4:3}
    assert "Gamma no_color" in z["deep_tube"]["known"]
    assert "Gamma no_color" in z["actinomorphic_deep"]["known"]
    assert "Gamma no_color" not in z["white_radial_tubular_deep"]["known"]
    assert "Alpha one" in z["white_radial_tubular_deep"]["joint"]
    assert "Alpha one" in z["tubular_bell_deep"]["joint"]
    assert "Alpha two" not in z["tubular_bell_deep"]["joint"]


def test_partial_species_must_not_be_imputed_negative_for_missing_axes():
    z=make_real_multiresolution_motifs(inputs())
    known=z["white_radial_deep"]["known"]
    assert known==frozenset({"Alpha one","Alpha two","Beta three"})
    assert "Gamma no_color" not in known
    assert "Gamma no_color" not in z["white_radial_deep"]["joint"]
    assert z["open_radial"]["axis_count"]==1
