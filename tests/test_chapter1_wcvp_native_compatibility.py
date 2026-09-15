import pandas as pd

from island_v2.chapter1_wcvp_native_compatibility import build_regional_native_compatible_flora


def test_wcvp_only_upgrades_unresolved_compatible_rows():
    flora = pd.DataFrame([
        {"island_id":"i1","accepted_species":"sp1","origin_status":"unresolved","floristic_status":"unresolved"},
        {"island_id":"i1","accepted_species":"sp2","origin_status":"introduced","floristic_status":"introduced"},
        {"island_id":"i1","accepted_species":"sp3","origin_status":"native","floristic_status":"native_nonendemic"},
        {"island_id":"i2","accepted_species":"sp1","origin_status":"unresolved","floristic_status":"unresolved"},
    ])
    ranges = pd.DataFrame([
        {"accepted_species":"sp1","native_l3_codes":"AAA|BBB"},
        {"accepted_species":"sp2","native_l3_codes":"AAA"},
        {"accepted_species":"sp3","native_l3_codes":"AAA"},
    ])
    mapping = pd.DataFrame([
        {"island_id":"i1","tdwg_l3_code":"AAA","tdwg_match_status":"accepted"},
        {"island_id":"i2","tdwg_l3_code":"CCC","tdwg_match_status":"accepted"},
    ])
    out, audit = build_regional_native_compatible_flora(flora, ranges, mapping)
    x=out.set_index(["island_id","accepted_species"])
    assert x.loc[("i1","sp1"),"origin_status"] == "native"
    assert x.loc[("i1","sp1"),"floristic_status"] == "native_endemism_unresolved"
    assert x.loc[("i1","sp2"),"origin_status"] == "introduced"
    assert x.loc[("i1","sp3"),"origin_status"] == "native"
    assert x.loc[("i2","sp1"),"origin_status"] == "unresolved"
    assert audit["n_wcvp_compatible_upgrades"] == 1
