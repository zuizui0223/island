"""Build the Wave 20 reviewed support-two evidence checkpoint.

Wave 20 records exact-species statements from primary literature, official
floras, institutional species profiles, and one morphology-only nursery page.
It does not create genus inference; the shared all-evidence rebuild remains the
sole producer of Validated Low rows.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from island_v2.targeted_support2_wave15_checkpoint import (
    AUDIT_COLUMNS,
    EVIDENCE_COLUMNS,
    _sha,
)
from island_v2.targeted_support2_wave15_checkpoint import (
    _row as _wave15_row,
)
from island_v2.targeted_support2_wave15_checkpoint import (
    build_audit as _wave15_build_audit,
)

SOURCE_GROUP = "targeted_support2_wave20_checkpoint_20260821"
CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/targeted_support2_wave20_checkpoint_20260821"
)
PRIOR = Path(
    "data/v2/staging/traits/open_web_pilot/targeted_support2_wave19_checkpoint_20260820"
)
RETRIEVED_AT = "2026-08-21T05:00:00Z"


def _row(*args: str, **kwargs: str) -> dict[str, str]:
    row = _wave15_row(*args, **kwargs)
    row["source_group"] = SOURCE_GROUP
    row["retrieved_at_utc"] = RETRIEVED_AT
    return row


def primary_rows() -> list[dict[str, str]]:
    """Return 24 individually reviewed, exact-species evidence rows."""
    helicteres_excerpt = (
        "After being dusted with self-pollen, the flowers of both species were "
        "observed over several days after which the ovaries were dissected and "
        "found to have set seed. This confirms that the flowers of M. arboreus "
        "are self-compatible, as shown in previous studies (Webb, 1983), as are "
        "those of H. guazumifolia. Visitors to target flowers of Malvaviscus "
        "arboreus and Helicteres guazumifolia over the day were recorded."
    )
    rows = [
        _row(
            "Helicteres guazumifolia",
            "reproductive_assurance",
            "self_incompatibility",
            "self-pollen treatment set seed; self-compatible",
            "SC",
            "high",
            "University of St Andrews Research Repository",
            (
                "https://research-repository.st-andrews.ac.uk/bitstream/handle/"
                "10023/3219/CarolineKingPhDThesis.pdf?isAllowed=y&sequence=3"
            ),
            "Putting Pollination Quality into Analyses of Floral Ecology",
            (
                "King 2012. PhD thesis, University of St Andrews, p. 96; "
                "original controlled self-pollen observation."
            ),
            helicteres_excerpt,
            "hdl:10023/3219:king-2012-original-self-pollen-experiment",
            "A",
            "original_university_thesis_controlled_self_pollination",
            "en",
            '"Helicteres guazumifolia" self compatible self pollen seed set',
            content_sha256=(
                "3d46d524bdf798d099f218e6a70a8a309d60bd2ea46be05d11b7c8ae3c132f22"
            ),
            content_sha256_basis="retrieved_original_thesis_pdf_bytes",
        ),
        _row(
            "Patersonia sericea",
            "floral_structural_complexity",
            "tube_depth_class",
            "perianth tube 1.5-3 cm long",
            "deep",
            "high",
            "Royal Botanic Gardens Victoria / VicFlora",
            (
                "https://vicflora.rbg.vic.gov.au/flora/taxon/"
                "74f4ad73-9c9c-4ad1-9bb1-505dec4e6382"
            ),
            "Patersonia sericea - VicFlora",
            (
                "Conn 1994. Iridaceae. Flora of Victoria 2:686-716; "
                "official VicFlora accepted-species treatment."
            ),
            (
                "Perianth tube 1.5-3 cm long, basal half hairy; outer perianth "
                "lobes broadly ovate, obtuse, 2-3 cm long, 1.5-2.5 cm wide; "
                "inner perianth lobes ovate."
            ),
            "flora-treatment:flora-of-victoria:Patersonia_sericea",
            "A",
            "official_flora_exact_species_treatment",
            "en",
            '"Patersonia sericea" "Perianth tube"',
            content_sha256=(
                "54e0ff80fe1aa49b58d3c1e536e3965697afe46439d0cd2f41d3bad2ad6a9473"
            ),
            content_sha256_basis="retrieved_official_species_page_bytes",
        ),
        _row(
            "Bocquillonia castaneifolia",
            "flower_colour",
            "flower_primary_color",
            "les fleurs sont rouges",
            "red_pink",
            "medium",
            "Association pour la sauvegarde du patrimoine minier et historique",
            "https://endemia.nc/files/Tiebaghi-V-2022.pdf",
            "Tiébaghi: de la flore exceptionnelle de la montagne",
            (
                "ASPMHNC 2022. Tiébaghi flora guide, p. 44; exact species "
                "profile redistributed by Endemia."
            ),
            (
                "BOCQUILLONIA CASTANEIFOLIA. FAMILLE : Euphorbiaceae. "
                "Les feuilles sont densément groupées au sommet des tiges. "
                "Les jeunes feuilles sont rouges ainsi que les fleurs."
            ),
            "regional-guide:aspmhnc:tiebaghi-2022:Bocquillonia_castaneifolia",
            "B",
            "regional_conservation_guide_exact_species_profile",
            "fr",
            '"Bocquillonia castaneifolia" fleurs rouges',
        ),
        _row(
            "Scabiosa triandra",
            "floral_structural_complexity",
            "flower_size_class",
            "corolla ca. 1 cm long",
            "small",
            "medium",
            "Provincia di Prato",
            "https://www.provincia.prato.it/pagina2509_dipsacaceae.html",
            "Dipsacaceae - Provincia di Prato",
            "Provincia di Prato. Scheda botanica di Scabiosa triandra L.",
            (
                "Scabiosa triandra L. I fiori - corolla, da color lillà a rosa "
                "e lunga ca. 1 cm - sono riuniti in capolini numerosi, del "
                "diametro di ca. 2-3 cm. Habitat: radure boschive, prati."
            ),
            "official-regional-flora:provincia-prato:Scabiosa_triandra",
            "B",
            "official_regional_flora_exact_species_page",
            "it",
            '"Scabiosa triandra" corolla 1 cm',
        ),
        _row(
            "Schizomeria ovata",
            "flower_colour",
            "flower_primary_color",
            "white flowers",
            "white",
            "medium",
            "Burringbar Rainforest Nursery",
            (
                "https://burringbarrainforestnursery.com.au/plant-type/"
                "rainforest-height-10-15m/?orderby=price-desc"
            ),
            "Rainforest height 10-15m - Burringbar Rainforest Nursery",
            "Burringbar Rainforest Nursery exact native-species listing.",
            (
                "Schizomeria ovata - CRAB APPLE. Attractive tree to 12m with "
                "bushy crown, showy bronze pink new shoots, white flowers & "
                "white fruit. NE.Qld. to SE.NSW."
            ),
            "provider-treatment:burringbar-rainforest-nursery:Schizomeria_ovata",
            "C",
            "native_nursery_exact_species_morphology_only",
            "en",
            '"Schizomeria ovata" "white flowers"',
            content_sha256=(
                "eab8dc8c148e26ce75215eab39f678261b4681a83a2221653c5e544df980cab8"
            ),
            content_sha256_basis="retrieved_source_page_bytes",
        ),
    ]

    institutional = [
        (
            "Cassipourea barteri",
            "floral_structural_complexity",
            "floral_symmetry",
            "radiär-symmetrisch",
            "actinomorphic",
            "6023",
            "Cassipourea barteri (Hook. f. ex Oliv.) N.E. Br.",
            (
                "Cassipourea barteri (Hook. f. ex Oliv.) N.E. Br. Blüte: "
                "Farbe weiß; Symmetrie radiär-symmetrisch; Blütenstand büschelig."
            ),
            "c70466a8883f9c6c8882d51aa258a23b09a098635edec82bec0f06dbcce7a76c",
        ),
        (
            "Cassipourea barteri",
            "flower_colour",
            "flower_primary_color",
            "Blüte Farbe weiß",
            "white",
            "6023",
            "Cassipourea barteri (Hook. f. ex Oliv.) N.E. Br.",
            (
                "Cassipourea barteri (Hook. f. ex Oliv.) N.E. Br. Blüte: "
                "Farbe weiß; Symmetrie radiär-symmetrisch; Blütenstand büschelig."
            ),
            "c70466a8883f9c6c8882d51aa258a23b09a098635edec82bec0f06dbcce7a76c",
        ),
        (
            "Maerua angolensis",
            "floral_structural_complexity",
            "floral_symmetry",
            "radiär-symmetrisch",
            "actinomorphic",
            "1069",
            "Maerua angolensis DC.",
            (
                "Maerua angolensis DC. Blüte: Farbe weiß gelb; Symmetrie "
                "radiär-symmetrisch; Blütenstand einzel- oder wenigblütig."
            ),
            "44acde6ddca63da3d4a1caed7f3d65c7a312bba97970c44482334763e31e6930",
        ),
        (
            "Acridocarpus longifolius",
            "floral_structural_complexity",
            "inflorescence_display",
            "Blütenstand traubig",
            "raceme_spike_panicle",
            "5952",
            "Acridocarpus longifolius (G. Don) Hook. f.",
            (
                "Acridocarpus longifolius (G. Don) Hook. f. Blüte: Farbe gelb; "
                "Symmetrie zygomorph; Blütenstand traubig."
            ),
            "68667f6cfc280b6cde5024fa221aeb336f93e90947265708a7b7a3b37023d57c",
        ),
        (
            "Acridocarpus longifolius",
            "floral_structural_complexity",
            "floral_symmetry",
            "Symmetrie zygomorph",
            "zygomorphic",
            "5952",
            "Acridocarpus longifolius (G. Don) Hook. f.",
            (
                "Acridocarpus longifolius (G. Don) Hook. f. Blüte: Farbe gelb; "
                "Symmetrie zygomorph; Blütenstand traubig."
            ),
            "68667f6cfc280b6cde5024fa221aeb336f93e90947265708a7b7a3b37023d57c",
        ),
        (
            "Amischotolype tenuis",
            "floral_structural_complexity",
            "floral_symmetry",
            "radiär-symmetrisch",
            "actinomorphic",
            "5830",
            "Amischotolype tenuis (C.B. Clarke) R.S. Rao",
            (
                "Amischotolype tenuis (C.B. Clarke) R.S. Rao. Blüte: Farbe "
                "rosa lila rot; Symmetrie radiär-symmetrisch; Blütenstand ährig."
            ),
            "6a42677e6f1a89b8b0c0d02d191b2c55bdd832a68d65b1aa6ddb7a2828144822",
        ),
        (
            "Amischotolype tenuis",
            "flower_colour",
            "flower_primary_color",
            "Blüte Farbe rosa lila rot",
            "blue_purple|red_pink",
            "5830",
            "Amischotolype tenuis (C.B. Clarke) R.S. Rao",
            (
                "Amischotolype tenuis (C.B. Clarke) R.S. Rao. Blüte: Farbe "
                "rosa lila rot; Symmetrie radiär-symmetrisch; Blütenstand ährig."
            ),
            "6a42677e6f1a89b8b0c0d02d191b2c55bdd832a68d65b1aa6ddb7a2828144822",
        ),
    ]
    for species, axis, trait, raw, value, record, name, excerpt, digest in institutional:
        rows.append(
            _row(
                species,
                axis,
                trait,
                raw,
                value,
                "medium",
                "Senckenberg African Plants Photo Guide",
                f"https://africanplants.senckenberg.de/index.php/item/{record}",
                f"{name} | African Plants Photo Guide",
                "Senckenberg institutional exact-species profile.",
                excerpt,
                f"institutional-species-profile:senckenberg-african-plants:{record}",
                "B",
                "institutional_exact_species_profile",
                "de",
                f'"{species}" {raw}',
                content_sha256=digest,
                content_sha256_basis="retrieved_institutional_species_page_bytes",
            )
        )

    catesbaea_excerpt = (
        "SMALL-FLOWERED LILY THORN. Catesbya parviflora Swartz. Synonyms: "
        "Echinodendrum parviflorum (Sw.) A. Rich.; Catesbaea parviflora Sw. "
        "Flowers less than 13 mm long, white, tubular with 4 spreading lobes, "
        "borne on spur shoots. Habitat: Pine rocklands, hammocks, coastal "
        "strands, back dunes."
    )
    for trait, raw, value in (
        ("flower_primary_color", "flowers white", "white"),
        ("floral_form", "flowers tubular with 4 spreading lobes", "tubular"),
        ("flower_size_class", "flowers less than 13 mm long", "small"),
    ):
        axis = "flower_colour" if trait == "flower_primary_color" else "floral_structural_complexity"
        rows.append(
            _row(
                "Catesbaea parviflora",
                axis,
                trait,
                raw,
                value,
                "high",
                "Florida Natural Areas Inventory",
                "https://www.fnai.org/PDFs/FieldGuides/Catesbya_parviflora.pdf",
                "Small-flowered lily thorn - Catesbaea parviflora",
                "Florida Natural Areas Inventory species field guide (2000).",
                catesbaea_excerpt,
                "official-species-profile:fnai:Catesbaea_parviflora",
                "A",
                "official_conservation_species_profile",
                "en",
                '"Catesbaea parviflora" flowers white tubular',
                matched_name="Catesbya parviflora",
                scope="synonym_direct",
                name_match_method="exact_synonym",
                name_resolution_lineage=(
                    "source_prints_Catesbya_parviflora_and_exact_accepted_synonym_"
                    "Catesbaea_parviflora"
                ),
                content_sha256=(
                    "f8375a9aa5a0fceef4bc93ab60a9de46b5a89d145773d1955c13f3ed70ca6a8e"
                ),
                content_sha256_basis="retrieved_official_field_guide_pdf_bytes",
            )
        )

    chrysanthemum_excerpt = (
        "3. Chrysanthemum indicum Linnaeus, Sp. Pl. 2: 889. 1753. "
        "Synflorescence a lax terminal flat-topped cyme. Capitula many or few. "
        "Phyllaries in 5 rows. Ray floret lamina yellow, 1-1.3 cm."
    )
    rows.append(
        _row(
            "Chrysanthemum indicum",
            "floral_structural_complexity",
            "floral_form",
            "capitula with ray florets",
            "composite_head",
            "high",
            "Flora of China / eFloras",
            "https://www.efloras.org/florataxon.aspx?flora_id=2&taxon_id=220002857",
            "Chrysanthemum indicum in Flora of China",
            "Flora of China 20-21:669-671; exact accepted-species treatment.",
            chrysanthemum_excerpt,
            "flora-treatment:flora-of-china:220002857",
            "A",
            "official_flora_exact_species_treatment",
            "en",
            '"Chrysanthemum indicum" capitula ray floret',
            content_sha256=(
                "c17665807cab9a60c2956eda2ae7fa7e3c374985e5400722f32104e0efdfc590"
            ),
            content_sha256_basis="retrieved_official_flora_page_bytes",
        )
    )

    critonia_excerpt = (
        "Critonia macropoda. Glabrous shrubs to 2 m; leaves opposite on "
        "1.5-2.0 cm petioles; inflorescence corymbose, heads with 5 tubular "
        "florets in 3-6 fascicles; involucral bracts in 4-5-series."
    )
    for trait, raw, value in (
        ("floral_form", "heads with 5 tubular florets", "tubular"),
        (
            "inflorescence_display",
            "corymbose inflorescence; heads in 3-6 fascicles",
            "composite_display|umbel_corymb",
        ),
    ):
        rows.append(
            _row(
                "Critonia macropoda",
                "floral_structural_complexity",
                trait,
                raw,
                value,
                "high",
                "Smithsonian Institution",
                "https://www.govinfo.gov/content/pkg/GOVPUB-SI-PURL-gpo111574/pdf/GOVPUB-SI-PURL-gpo111574.pdf",
                "Flora of Dominica, Part 2: Dicotyledoneae",
                (
                    "Nicolson et al. 1991. Smithsonian Contributions to Botany "
                    "77:38; exact species treatment."
                ),
                critonia_excerpt,
                "flora-treatment:smithsonian-flora-dominica-1991:Critonia_macropoda",
                "A",
                "official_flora_exact_species_treatment",
                "en",
                f'"Critonia macropoda" {raw}',
                content_sha256=(
                    "60317020f4720e5945b8f6b28e05d079fcf3b982b9cbfbb9a87503e4e4332ec2"
                ),
                content_sha256_basis="retrieved_smithsonian_flora_pdf_bytes",
            )
        )

    aeonium_excerpt = (
        "Aeonium appendiculatum Bañares, sp. nova. Inflorescence dome-shaped, "
        "25-40 x 20-30 cm, glabrous; branches 5-15 cm with lanceolate, "
        "acuminate bracts and 25-90 flowers; pedicels 1.5-3 mm long, glabrous. "
        "Flowers 8-merous, 1.4 cm in diameter; calyx 4 mm wide, base rounded, "
        "glabrous; petals lanceolate, entire, 6-6.3 x 2-2.2 mm, white, median "
        "portion pink variegated, glabrous."
    )
    rows.append(
        _row(
            "Aeonium appendiculatum",
            "floral_structural_complexity",
            "flower_size_class",
            "flowers 1.4 cm in diameter",
            "medium",
            "high",
            "Botanic Garden and Botanical Museum Berlin / Willdenowia",
            "https://www.bgbm.org/sites/default/files/documents/w29Banares.pdf",
            "Notes on the taxonomy of Aeonium urbicum and A. appendiculatum",
            (
                "Bañares Baudet 1999. Willdenowia 29:95-103, pp. 98-101; "
                "original species description based on wild La Gomera material."
            ),
            aeonium_excerpt,
            "doi:10.3372/wi.29.2908:Aeonium_appendiculatum",
            "A",
            "primary_taxonomic_species_description",
            "en",
            '"Aeonium appendiculatum" flower diameter',
            content_sha256=(
                "98ce31080062bba5477c882b46d609840f39b7a70aeec48a4036363ac5ce771e"
            ),
            content_sha256_basis="retrieved_primary_taxonomic_pdf_bytes",
        )
    )

    rows.append(
        _row(
            "Arytera microphylla",
            "flower_colour",
            "flower_primary_color",
            "flowers white",
            "white",
            "medium",
            "Native Plants Queensland",
            (
                "https://npq.org.au/kingaroy-districts/dryrainforestplants/"
                "small-leavedcoogera/"
            ),
            "Small-Leaved Coogera - Arytera microphylla",
            (
                "Native Plants Queensland exact-species dry-rainforest profile; "
                "page updated 22 January 2025."
            ),
            (
                "Arytera microphylla. Inflorescences axillary, thyrses, "
                "2.5-8 cm long; peduncles sparsely hairy. Flowers functionally "
                "unisexual and plants monoecious, white, 2-3 mm diam. Habitat "
                "and Distribution: In DRf; north from Yarraman, Qld, to the "
                "Bundaberg district, Qld."
            ),
            "specialist-regional-profile:native-plants-queensland:Arytera_microphylla",
            "B",
            "specialist_regional_native_plant_profile",
            "en",
            '"Arytera microphylla" flowers white',
            content_sha256=(
                "6837079a69150bbde3f58a501054a18fb21c0fd055526d3f32c91ce73d551204"
            ),
            content_sha256_basis="retrieved_source_page_bytes",
        )
    )

    paraboea_excerpt = (
        "Paraboea barnettiae C. Puglisi. "
        "ช่อดอกแบบช่อกระจุกด้านเดียวเชิงประกอบ มี ๒-๖ ช่อ ออกตามปลายกิ่ง "
        "ก้านช่อดอกยาว ๓-๑๕ ซม. ดอกสีม่วงอมฟ้า"
    )
    for trait, raw, value, axis in (
        (
            "inflorescence_display",
            "ช่อดอกแบบช่อกระจุกด้านเดียวเชิงประกอบ (compound unilateral cyme)",
            "umbel_corymb",
            "floral_structural_complexity",
        ),
        (
            "flower_primary_color",
            "ดอกสีม่วงอมฟ้า (flowers purple-blue)",
            "blue_purple",
            "flower_colour",
        ),
    ):
        rows.append(
            _row(
                "Paraboea barnettiae",
                axis,
                trait,
                raw,
                value,
                "medium",
                "NECTEC Thai Plant Taxonomy Encyclopedia",
                (
                    "https://lst.nectec.or.th/encyclopedia/wikipedia/plant_taxonomy/"
                    "view/index.php?cid=189&namelink=%E0%B8%8A%E0%B9%89%E0%B8%B2%"
                    "E0%B8%AB%E0%B8%99%E0%B8%B2%E0%B8%94%E0%B9%80%E0%B8%82%E0%B8%B2"
                ),
                "Paraboea barnettiae C. Puglisi - ช้าหนาดเขา",
                (
                    "Thai Plant Taxonomy Encyclopedia exact-species profile by "
                    "Dr Pramote Triboun; wild limestone habitat in southern Thailand."
                ),
                paraboea_excerpt,
                "institutional-species-profile:nectec-thai-taxonomy:cid-189",
                "A",
                "national_institutional_taxonomic_species_profile",
                "th",
                f'"Paraboea barnettiae" {raw}',
                content_sha256=(
                    "d5e44ad61877e22bc1ff3ad1ab1aca693bb2740babfcf0583cc7b8c51c084912"
                ),
                content_sha256_basis="retrieved_institutional_species_page_bytes",
            )
        )

    xanthosoma_excerpt = (
        "Xanthosoma brasiliense (Desf.) Engl. Inflorescences 1-2 per axil, "
        "peduncle 20-25 x .5 cm; spathe 18-19 cm long, tube 5 x 2.5 cm, green "
        "outside, white inside, lamina 13-14 x 3 cm long, white in both sides; "
        "spadix 14-16 cm long, fertile male portion white 10 x 1-1.5 cm, acute "
        "at apex, sterile male portion 3.5 x 1.1 cm, white, weakly dimorphic, "
        "female portion 2-3 x 1 cm, pale yellow."
    )
    for trait, raw, value, axis in (
        (
            "flower_primary_color",
            "spathe and male portions white; female portion pale yellow",
            "white|yellow_orange",
            "flower_colour",
        ),
        (
            "inflorescence_display",
            "spadix 14-16 cm long",
            "raceme_spike_panicle",
            "floral_structural_complexity",
        ),
    ):
        rows.append(
            _row(
                "Xanthosoma brasiliense",
                axis,
                trait,
                raw,
                value,
                "medium",
                "Royal Botanic Gardens, Kew / CATE Araceae",
                (
                    "https://powo.science.kew.org/taxon/urn:lsid:ipni.org:"
                    "names:269304-2/general-information"
                ),
                "Xanthosoma brasiliense - Plants of the World Online",
                (
                    "CATE Araceae 2011 exact-species treatment distributed by "
                    "Royal Botanic Gardens, Kew; species-level cultivated material."
                ),
                xanthosoma_excerpt,
                "provider-treatment:cate-araceae:Xanthosoma_brasiliense:2011",
                "A",
                "institutional_monograph_species_treatment",
                "en",
                f'"Xanthosoma brasiliense" {raw}',
                status="species_level_cultivated_material_no_named_cultivar_or_hybrid",
            )
        )

    if len(rows) != 24:
        raise AssertionError("Wave20 must define exactly 24 direct rows")
    return rows


def build_audit(evidence: pd.DataFrame) -> pd.DataFrame:
    audit = _wave15_build_audit(evidence)
    audit["reviewer"] = "Codex Wave20 individual source audit"
    audit["reviewed_at_utc"] = RETRIEVED_AT
    nursery = audit["source_tier"].eq("C")
    audit.loc[nursery, "decision_reason"] = (
        "Accepted as limited Medium for exact-species flower colour only; the "
        "native nursery listing names no cultivar or hybrid and is not used for "
        "reproductive evidence or High-tier inference."
    )
    catesbaea = audit["accepted_species"].eq("Catesbaea parviflora")
    audit.loc[catesbaea, "decision_reason"] = (
        "Accepted after the official FNAI page image was rendered and visually "
        "checked: it prints both the historic spelling and exact accepted name, "
        "then states flower size, colour, and tubular form in the field description."
    )
    xanthosoma = audit["accepted_species"].eq("Xanthosoma brasiliense")
    audit.loc[xanthosoma, "decision_reason"] = (
        "Accepted as Medium species-direct morphology only. The institutional "
        "treatment explicitly identifies Xanthosoma brasiliense and records the "
        "spathe/spadix states, but the observed material is cultivated; no named "
        "cultivar or hybrid is transferred and no reproductive inference is made."
    )
    return audit.loc[:, AUDIT_COLUMNS]


def build_checkpoint(output_dir: Path = CHECKPOINT) -> dict[str, Any]:
    evidence = pd.DataFrame(primary_rows(), columns=EVIDENCE_COLUMNS).sort_values(
        ["accepted_species", "trait_name", "source_lineage"], kind="stable"
    )
    evidence = evidence.reset_index(drop=True)
    audit = build_audit(evidence)
    if len(evidence) != 24 or len(audit) != 24:
        raise ValueError("Wave20 must contain exactly 24 individually reviewed rows")
    if evidence["candidate_id"].duplicated().any():
        raise ValueError("Wave20 candidate IDs must be unique")
    if evidence.duplicated(["accepted_species", "trait_name"]).any():
        raise ValueError("Wave20 species x trait pairs must be unique")

    prior_evidence = pd.read_csv(
        PRIOR / "combined_curated_evidence_20260820.csv", dtype=str
    ).fillna("")
    prior_audit = pd.read_csv(
        PRIOR / "combined_curated_manual_audit_20260820.csv", dtype=str
    ).fillna("")
    combined_evidence = pd.concat([prior_evidence, evidence], ignore_index=True)
    combined_audit = pd.concat([prior_audit, audit], ignore_index=True)
    if combined_evidence["candidate_id"].duplicated().any():
        raise ValueError("combined evidence candidate IDs must be unique")
    if combined_audit["candidate_id"].duplicated().any():
        raise ValueError("combined audit candidate IDs must be unique")

    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "evidence": output_dir / "targeted_support2_wave20_evidence_20260821.csv",
        "audit": output_dir / "targeted_support2_wave20_manual_audit_20260821.csv",
        "combined_evidence": output_dir / "combined_curated_evidence_20260821.csv",
        "combined_audit": output_dir / "combined_curated_manual_audit_20260821.csv",
        "manifest": output_dir / "source_acquisition_manifest_wave20.json",
    }
    evidence.to_csv(paths["evidence"], index=False, lineterminator="\n")
    audit.to_csv(paths["audit"], index=False, lineterminator="\n")
    combined_evidence.to_csv(paths["combined_evidence"], index=False, lineterminator="\n")
    combined_audit.to_csv(paths["combined_audit"], index=False, lineterminator="\n")

    manifest = {
        "checkpoint": SOURCE_GROUP,
        "built_at_utc": RETRIEVED_AT,
        "baseline_formal_run_id": 32347770192,
        "accepted_evidence_rows": len(evidence),
        "accepted_species_trait": int(
            evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
        ),
        "accepted_species": int(evidence["accepted_species"].nunique()),
        "recorded_queries": int(evidence["query"].nunique()),
        "formal_search_api_queries": 0,
        "search_cost_usd": 0.0,
        "targeted_support2_rules": [
            "Acridocarpus|inflorescence_display",
            "Amischotolype|floral_symmetry",
            "Arytera|flower_primary_color",
            "Bocquillonia|flower_primary_color",
            "Cassipourea|floral_symmetry",
            "Catesbaea|flower_primary_color",
            "Chrysanthemum|floral_form",
            "Critonia|floral_form",
            "Helicteres|self_incompatibility",
            "Maerua|floral_symmetry",
            "Patersonia|tube_depth_class",
            "Paraboea|inflorescence_display",
            "Scabiosa|flower_size_class",
            "Schizomeria|flower_primary_color",
        ],
        "theoretical_rule_cells_touched": 179,
        "guardrails": {
            "search_snippet_as_evidence": False,
            "family_inference": False,
            "global_fallback": False,
            "min_species_two_production": False,
            "cross_trait_substitution": False,
            "genus_axis_only_join": False,
            "pollen_vector_or_reward_mixed_into_structure": False,
            "self_fertility_silently_mapped_to_autonomous_selfing": False,
            "cultivar_or_hybrid_transferred_to_wild_species": False,
        },
        "output_sha256": {
            key: _sha(path.read_text(encoding="utf-8"))
            for key, path in paths.items()
            if key != "manifest"
        },
        "notes": (
            "The 179 cells are a support-two queue ceiling, not observed coverage "
            "gain. Formal all-evidence rebuilding and species/lineage leave-one-out "
            "validation determine realized direct and Validated Low changes. The "
            "remaining rows preserve additional species-direct traits from the same "
            "reviewed treatments without expanding their allowed evidence tier."
        ),
    }
    paths["manifest"].write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=CHECKPOINT)
    args = parser.parse_args()
    print(json.dumps(build_checkpoint(args.output_dir), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
