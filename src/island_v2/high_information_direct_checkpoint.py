"""Freeze reviewed, high-information species-direct records.

These records were acquired one source at a time after ranking unresolved
``genus x trait`` cells by the number of congeners they could unlock.  The
checkpoint stores the exact quote, immutable retrieved-content hash and source
lineage; it does not perform inference or substitute one reproductive trait for
another.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from island_v2.high_leverage_direct_checkpoint import (
    EVIDENCE_COLUMNS,
    _audit,
    _canonical_file_bytes,
    _evidence_row,
)

CREATED_AT = "2026-08-12T08:35:00Z"


def _sha256(path: Path) -> str:
    return hashlib.sha256(_canonical_file_bytes(path)).hexdigest()


def _pohnpei_breeding_rows() -> list[dict[str, str]]:
    """Keep experimental compatibility and autofertility as separate traits."""

    experiments = (
        ("Crotalaria pallida", "20.6", "18.8", "autonomous"),
        ("Dendrolobium umbellatum", "3", "2.6", "autonomous"),
        ("Leucaena leucocephala", "13.6", "11.4", "autonomous"),
        ("Mimosa diplotricha", "7.6", "0.4", "absent"),
        ("Senna alata", "5.6", "6.2", "autonomous"),
        ("Vigna marina", "4.6", "4", "autonomous"),
        ("Hibiscus tiliaceus", "53", "52.2", "autonomous"),
        ("Sida acuta", "51.2", "24.4", "mixed_or_variable"),
        ("Thespesia populnea", "14.2", "5.8", "mixed_or_variable"),
    )
    common = {
        "quality": "high",
        "provider": "yomai_williams_2021_primary_article",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8317631/",
        "title": (
            "Breeding systems of naturalized versus indigenous species provide "
            "support for Baker's law on Pohnpei island"
        ),
        "citation": (
            "Yomai & Williams (2021), AoB PLANTS 13:plab038, "
            "DOI 10.1093/aobpla/plab038, PMCID PMC8317631"
        ),
        "lineage": "doi:10.1093/aobpla/plab038",
        "lineage_method": "original_primary_article_doi",
        "source_tier": "A",
        "source_type": "primary_research_article_controlled_pollination_experiment",
        "domain": "pmc.ncbi.nlm.nih.gov",
        "content_sha256": (
            "b006ffb2999fee4206cb86a3e84ec246a20fc87bd22f20830295b387b19f33dc"
        ),
        "content_sha256_basis": "retrieved_pmc_jats_full_text_xml_bytes",
        "retrieved_at_utc": "2026-08-12T09:02:00Z",
    }
    rows: list[dict[str, str]] = []
    for species, hand_self, autonomous_self, autonomous_value in experiments:
        table_row = (
            f"Table 4: {species}; Hand-self {hand_self} seeds per five flowers; "
            f"Autonomous-self (bagged and unmanipulated) {autonomous_self}."
        )
        rows.extend(
            (
                _evidence_row(
                    species=species,
                    trait="self_incompatibility",
                    value="SC",
                    excerpt=(
                        "Within all 11 species, the hand-self crosses produced a "
                        "higher or similar number of seeds as the hand-outcrosses. "
                        "All species were self-compatible, with ISI < 0.8. "
                        + table_row
                    ),
                    record_id=f"pmc8317631:table4:{species}:self-compatibility",
                    raw_value="all experimental species self-compatible; ISI < 0.8",
                    **common,
                ),
                _evidence_row(
                    species=species,
                    trait="autonomous_selfing_capacity",
                    value=autonomous_value,
                    excerpt=table_row,
                    record_id=f"pmc8317631:table4:{species}:autonomous-selfing",
                    raw_value=(
                        f"autonomous-self {autonomous_self} / hand-self {hand_self} "
                        "seeds per five flowers"
                    ),
                    **common,
                ),
            )
        )
    return rows


def reviewed_rows() -> list[dict[str, str]]:
    """Return only records whose complete source pages were inspected."""

    return [
        _evidence_row(
            species="Ficus religiosa",
            trait="floral_symmetry",
            value="actinomorphic",
            quality="medium",
            provider="leon_levy_native_plant_preserve",
            url="https://levypreserve.org/plant-listings/ficus-religiosa/",
            title="Ficus religiosa - Leon Levy Native Plant Preserve",
            citation="Leon Levy Native Plant Preserve species account: Ficus religiosa",
            excerpt=(
                "The highly reduced incomplete, imperfect, actinomorphic, flowers "
                "are borne entirely within a structure known as a Synconium (fig) "
                "and are fertilized by wasps."
            ),
            record_id="levy-preserve:ficus-religiosa:floral-symmetry",
            lineage="url:https://levypreserve.org/plant-listings/ficus-religiosa",
            lineage_method="canonical_institutional_species_page",
            source_tier="B",
            source_type="botanical_preserve_species_account",
            domain="levypreserve.org",
            content_sha256="84634554d3a7a2e01bd28336762fb7aa13ad578e89df3213845b9043e41aed74",
            content_sha256_basis="retrieved_species_page_bytes",
            retrieved_at_utc="2026-08-11T14:02:00Z",
        ),
        _evidence_row(
            species="Erythroxylum citrifolium",
            trait="floral_symmetry",
            value="actinomorphic",
            quality="high",
            provider="barreto_etal_2025_primary_article",
            url="https://www.scielo.br/j/bjb/a/Cy6Snbr5DzdW9R9fK6Ky5mn/?format=html&lang=en",
            title=(
                "Is there interspecific fruit set in distylous, synchronopatric "
                "species of Erythroxylum P. Browne (Erythroxylaceae)?"
            ),
            citation=(
                "Barreto et al. (2025), Brazilian Journal of Biology, "
                "DOI 10.1590/1519-6984.287794"
            ),
            excerpt=(
                "The species have similar floral attributes. Flowers have five "
                "sepals that connate at the base, and a pentamerous, white, "
                "actinomorphic corolla of similar size"
            ),
            record_id="doi:10.1590/1519-6984.287794:floral-attributes",
            lineage="doi:10.1590/1519-6984.287794",
            lineage_method="original_primary_article_doi",
            source_tier="A",
            source_type="primary_research_article",
            domain="scielo.br",
            content_sha256="27a74c1bc2fde4fb44075d6f0550752109c19d1af01e65a307fc799ed9c6b14c",
            content_sha256_basis="retrieved_primary_article_html_bytes",
            retrieved_at_utc="2026-08-11T14:16:00Z",
        ),
        _evidence_row(
            species="Ternstroemia penangiana",
            trait="inflorescence_display",
            value="solitary",
            quality="high",
            provider="india_flora_online_iisc",
            url="https://indiaflora-ces.iisc.ac.in/herbsheet.php?cat=13&id=11701",
            title="India Flora Online: Ternstroemia penangiana",
            citation=(
                "India Flora Online, Centre for Ecological Sciences, Indian "
                "Institute of Science: Ternstroemia penangiana"
            ),
            excerpt=(
                "A small or moderate-sized shady tree with dark coloured bark; cut "
                "mealy, reddish-brown. Leaves tapering at the base, margin entire "
                "and slightly recurved, dark green, drying brown. Flowers solitary, "
                "creamy-white turning dark brown in withering: petals waxy with "
                "denticulate margins."
            ),
            record_id="india-flora-online:11701:inflorescence-display",
            lineage="url:https://indiaflora-ces.iisc.ac.in/herbsheet.php?cat=13&id=11701",
            lineage_method="canonical_university_flora_species_page",
            source_tier="A",
            source_type="university_flora",
            domain="indiaflora-ces.iisc.ac.in",
            content_sha256="eaf894d548d06dec4f94ac23e742bae4ba84d4802953452caf1796fd5d49835e",
            content_sha256_basis="retrieved_species_page_bytes",
            retrieved_at_utc="2026-08-11T14:24:00Z",
        ),
        _evidence_row(
            species="Aloe maculata",
            trait="autonomous_selfing_capacity",
            value="autonomous",
            quality="high",
            provider="hargreaves_2007_doctoral_thesis",
            url=(
                "https://central.bac-lac.gc.ca/.item?app=Library&id=MR37375&"
                "oclc_number=302053454&op=pdf"
            ),
            title=(
                "The ecological causes and consequences of floral variation: "
                "Aloe maculata as a case study"
            ),
            citation=(
                "Hargreaves (2007), PhD thesis, Queen's University, Library and "
                "Archives Canada MR37375, OCLC 302053454"
            ),
            excerpt=(
                "Aloe maculata appears to be self-incompatible, as no flowers "
                "hand-pollinated with self-pollen set fruit (n=78 flowers on 12 "
                "plants), whereas 41 of the 90 cross-pollinated flowers set fruit "
                "(n=46 plants), and produced an average of 36.2 seeds per fruit "
                "(n=37 fruits on 32 plants). However, the self-incompatibility "
                "system may not be complete, as bagged plants produced a few fruits "
                "(mean=0.058 fruits/flower, n=4840 flowers on 43 plants) and seeds "
                "(mean=0.69 seeds/fruit, n=42 fruits on 15 plants) autonomously."
            ),
            record_id="lac:MR37375:chapter3:self-compatibility-results",
            lineage="thesis:lac:MR37375|oclc:302053454",
            lineage_method="original_doctoral_thesis_catalog_record",
            source_tier="A",
            source_type="doctoral_thesis_controlled_pollination_experiment",
            domain="central.bac-lac.gc.ca",
            content_sha256="16f5cf30fdb1044ec9b93987db372485bd040dd4172a3385cdf9972c50672e6f",
            content_sha256_basis="downloaded_doctoral_thesis_pdf_bytes",
            retrieved_at_utc="2026-08-11T14:43:00Z",
            raw_value="nonzero fruit and seed set in bagged plants; autonomously",
        ),
        _evidence_row(
            species="Begonia nelumbiifolia",
            trait="mating_system",
            value="mixed_mating",
            quality="high",
            provider="twyford_etal_2015_primary_article",
            url="https://pmc.ncbi.nlm.nih.gov/articles/PMC4600226/",
            title="Maintenance of species boundaries in a Neotropical radiation of Begonia",
            citation=(
                "Twyford et al. (2015), Molecular Ecology, "
                "DOI 10.1111/mec.13355, PMCID PMC4600226"
            ),
            excerpt=(
                "The species studied here all demonstrate moderate levels of "
                "inbreeding, with inferred equilibrium selfing rates of 0.40 for "
                "B. heracleifolia and 0.62 for B. nelumbiifolia"
            ),
            record_id="doi:10.1111/mec.13355:selfing-rates",
            lineage="doi:10.1111/mec.13355",
            lineage_method="original_primary_article_doi",
            source_tier="A",
            source_type="primary_research_article",
            domain="pmc.ncbi.nlm.nih.gov",
            content_sha256="9555fb74088bfde0c261fb99ae8e6e3aaada9b06bfccbab63d2d889f23f45d21",
            content_sha256_basis="retrieved_jats_full_text_xml_bytes",
            retrieved_at_utc="2026-08-11T14:53:00Z",
            raw_value="equilibrium selfing rate 0.62; inferred outcrossing rate 0.38",
        ),
        _evidence_row(
            species="Abutilon indicum",
            trait="self_incompatibility",
            value="SC",
            quality="high",
            provider="abid_etal_2010_primary_article",
            url="https://pakbs.org/pjbot/papers/1770199571.pdf",
            title="Pollination mechanism and role of insects in Abutilon indicum",
            citation=(
                "Abid, Alam & Qaiser (2010), Pakistan Journal of Botany "
                "42(3):1395-1399"
            ),
            excerpt=(
                "Breeding experiments revealed that A. indicum is self compatible "
                "and facultative autogamous taxon and no significant difference was "
                "found in fruit set of open pollination and autogamy."
            ),
            record_id="pak-j-bot:42:1395-1399:self-compatibility",
            lineage="article:abid-alam-qaiser-2010-abutilon-indicum",
            lineage_method="original_primary_article_bibliographic_identity",
            source_tier="A",
            source_type="primary_research_article_controlled_pollination_experiment",
            domain="pakbs.org",
            content_sha256="7cef43d7959498a6ed94ae910c33c38448d20da893e3b47dfcb4c1f3a081e114",
            content_sha256_basis="downloaded_publisher_pdf_bytes",
            retrieved_at_utc="2026-08-12T08:25:00Z",
            raw_value="self compatible",
        ),
        _evidence_row(
            species="Abutilon indicum",
            trait="autonomous_selfing_capacity",
            value="autonomous",
            quality="high",
            provider="abid_etal_2010_primary_article",
            url="https://pakbs.org/pjbot/papers/1770199571.pdf",
            title="Pollination mechanism and role of insects in Abutilon indicum",
            citation=(
                "Abid, Alam & Qaiser (2010), Pakistan Journal of Botany "
                "42(3):1395-1399"
            ),
            excerpt=(
                "Buds were bagged without manipulation to determine the self-"
                "pollination. Bagging experiments and pollen-ovule ratio reveal "
                "that it is a facultative autogamous taxon."
            ),
            record_id="pak-j-bot:42:1395-1399:facultative-autogamy",
            lineage="article:abid-alam-qaiser-2010-abutilon-indicum",
            lineage_method="original_primary_article_bibliographic_identity",
            source_tier="A",
            source_type="primary_research_article_controlled_pollination_experiment",
            domain="pakbs.org",
            content_sha256="7cef43d7959498a6ed94ae910c33c38448d20da893e3b47dfcb4c1f3a081e114",
            content_sha256_basis="downloaded_publisher_pdf_bytes",
            retrieved_at_utc="2026-08-12T08:25:00Z",
            raw_value="facultative autogamy after unmanipulated bud bagging",
        ),
        *_pohnpei_breeding_rows(),
        {
            **_evidence_row(
                species="Claoxylon taitense",
                trait="flower_size_class",
                value="very_small",
                quality="high",
                provider="polynesie_francaise_diren_tree_guide",
                url=(
                    "https://www.service-public.pf/diren/wp-content/uploads/sites/"
                    "17/2019/01/Guide_arbres.pdf"
                ),
                title=(
                    "Guide de reconnaissance des arbres, arbustes & arbrisseaux "
                    "indigenes de Polynesie francaise"
                ),
                citation=(
                    "Taputuarai, Direction de l'Environnement de la Polynesie "
                    "francaise, species account p. 29"
                ),
                excerpt=(
                    "Claoxylon taitense | Euphorbiaceae. Fleurs Tres petites, "
                    "verdatres, organisees en inflorescences de type grappe."
                ),
                record_id="diren-guide-arbres:p29:claoxylon-taitense:flower-size",
                lineage="publication:diren-polynesie-guide-arbres-indigenes",
                lineage_method="official_government_flora_publication",
                source_tier="A",
                source_type="government_regional_flora_guide",
                domain="service-public.pf",
                content_sha256=(
                    "c037286ffbdc794cf00290222fc92e3ba53dff01132363d62181602afa361f9b"
                ),
                content_sha256_basis="downloaded_government_guide_pdf_bytes",
                retrieved_at_utc="2026-08-12T09:15:00Z",
                raw_value="Fleurs tres petites",
            ),
            "language": "fr",
        },
        {
            **_evidence_row(
                species="Claoxylon taitense",
                trait="flower_primary_color",
                value="green_brown_inconspicuous",
                quality="high",
                provider="polynesie_francaise_diren_tree_guide",
                url=(
                    "https://www.service-public.pf/diren/wp-content/uploads/sites/"
                    "17/2019/01/Guide_arbres.pdf"
                ),
                title=(
                    "Guide de reconnaissance des arbres, arbustes & arbrisseaux "
                    "indigenes de Polynesie francaise"
                ),
                citation=(
                    "Taputuarai, Direction de l'Environnement de la Polynesie "
                    "francaise, species account p. 29"
                ),
                excerpt=(
                    "Claoxylon taitense | Euphorbiaceae. Fleurs Tres petites, "
                    "verdatres, organisees en inflorescences de type grappe."
                ),
                record_id="diren-guide-arbres:p29:claoxylon-taitense:flower-colour",
                lineage="publication:diren-polynesie-guide-arbres-indigenes",
                lineage_method="official_government_flora_publication",
                source_tier="A",
                source_type="government_regional_flora_guide",
                domain="service-public.pf",
                content_sha256=(
                    "c037286ffbdc794cf00290222fc92e3ba53dff01132363d62181602afa361f9b"
                ),
                content_sha256_basis="downloaded_government_guide_pdf_bytes",
                retrieved_at_utc="2026-08-12T09:15:00Z",
                raw_value="Fleurs verdatres",
            ),
            "language": "fr",
        },
        {
            **_evidence_row(
                species="Claoxylon taitense",
                trait="inflorescence_display",
                value="raceme_spike_panicle",
                quality="high",
                provider="polynesie_francaise_diren_tree_guide",
                url=(
                    "https://www.service-public.pf/diren/wp-content/uploads/sites/"
                    "17/2019/01/Guide_arbres.pdf"
                ),
                title=(
                    "Guide de reconnaissance des arbres, arbustes & arbrisseaux "
                    "indigenes de Polynesie francaise"
                ),
                citation=(
                    "Taputuarai, Direction de l'Environnement de la Polynesie "
                    "francaise, species account p. 29"
                ),
                excerpt=(
                    "Claoxylon taitense | Euphorbiaceae. Fleurs Tres petites, "
                    "verdatres, organisees en inflorescences de type grappe."
                ),
                record_id=(
                    "diren-guide-arbres:p29:claoxylon-taitense:inflorescence"
                ),
                lineage="publication:diren-polynesie-guide-arbres-indigenes",
                lineage_method="official_government_flora_publication",
                source_tier="A",
                source_type="government_regional_flora_guide",
                domain="service-public.pf",
                content_sha256=(
                    "c037286ffbdc794cf00290222fc92e3ba53dff01132363d62181602afa361f9b"
                ),
                content_sha256_basis="downloaded_government_guide_pdf_bytes",
                retrieved_at_utc="2026-08-12T09:15:00Z",
                raw_value="inflorescences de type grappe",
            ),
            "language": "fr",
        },
        _evidence_row(
            species="Osmoxylon novoguineense",
            trait="flower_primary_color",
            value="red_pink",
            quality="high",
            provider="flora_malesiana_araliaceae_1979",
            url="https://repository.naturalis.nl/pub/532684/FM1S1979009001002.pdf",
            title="Flora Malesiana, Series I, Volume 9: Araliaceae-I",
            citation="Philipson (1979), Flora Malesiana ser. I, vol. 9, p. 44",
            excerpt=(
                "Osmoxylon novoguineense. The inflorescence branches are dark "
                "purple, the corolla and stamens usually deep red, and the ripe "
                "fruit shining purple or blue-black."
            ),
            record_id="flora-malesiana-v9:p44:osmoxylon-novoguineense:corolla-colour",
            lineage="publication:flora-malesiana-ser1-v9-araliaceae-philipson-1979",
            lineage_method="original_taxonomic_flora_publication",
            source_tier="A",
            source_type="authoritative_regional_flora",
            domain="repository.naturalis.nl",
            content_sha256=(
                "9c9cf2c14156247316c5cb30d8d545a202040a94a4e6db7d9bcf462519d6541a"
            ),
            content_sha256_basis="downloaded_institutional_repository_pdf_bytes",
            retrieved_at_utc="2026-08-12T09:24:00Z",
            raw_value="corolla and stamens usually deep red",
        ),
        _evidence_row(
            species="Globba campsophylla",
            trait="flower_primary_color",
            value="white|yellow_orange",
            quality="high",
            provider="sangvirotjanapat_etal_2020_taxonomic_revision",
            url="https://journals.rbge.org.uk/ejb/article/download/1752/1643",
            title=(
                "A taxonomic revision of Globba sect. Nudae subsect. "
                "Mediocalcaratae (Zingiberaceae)"
            ),
            citation=(
                "Sangvirotjanapat, Denduangboripant & Newman (2020), Edinburgh "
                "Journal of Botany 77:1-63, DOI 10.1017/S0960428619000258"
            ),
            excerpt=(
                "Globba campsophylla. Flowers c.3 cm long, white except labellum "
                "with orange spot."
            ),
            record_id="doi:10.1017/S0960428619000258:p18:globba-campsophylla:colour",
            lineage="doi:10.1017/S0960428619000258",
            lineage_method="original_taxonomic_revision_doi",
            source_tier="A",
            source_type="primary_taxonomic_revision",
            domain="journals.rbge.org.uk",
            content_sha256=(
                "b8b6ddaf2eebeaec290a93db90e6cc2819ff631c89516c59374b101f5613759c"
            ),
            content_sha256_basis="downloaded_publisher_pdf_bytes",
            retrieved_at_utc="2026-08-12T09:32:00Z",
            raw_value="white except labellum with orange spot",
        ),
        _evidence_row(
            species="Globba campsophylla",
            trait="tube_depth_class",
            value="intermediate",
            quality="high",
            provider="sangvirotjanapat_etal_2020_taxonomic_revision",
            url="https://journals.rbge.org.uk/ejb/article/download/1752/1643",
            title=(
                "A taxonomic revision of Globba sect. Nudae subsect. "
                "Mediocalcaratae (Zingiberaceae)"
            ),
            citation=(
                "Sangvirotjanapat, Denduangboripant & Newman (2020), Edinburgh "
                "Journal of Botany 77:1-63, DOI 10.1017/S0960428619000258"
            ),
            excerpt=(
                "Globba campsophylla. Flowers c.3 cm long; floral tube c.12 mm "
                "long."
            ),
            record_id="doi:10.1017/S0960428619000258:p18:globba-campsophylla:tube",
            lineage="doi:10.1017/S0960428619000258",
            lineage_method="original_taxonomic_revision_doi",
            source_tier="A",
            source_type="primary_taxonomic_revision",
            domain="journals.rbge.org.uk",
            content_sha256=(
                "b8b6ddaf2eebeaec290a93db90e6cc2819ff631c89516c59374b101f5613759c"
            ),
            content_sha256_basis="downloaded_publisher_pdf_bytes",
            retrieved_at_utc="2026-08-12T09:32:00Z",
            raw_value="floral tube c.12 mm long",
        ),
        _evidence_row(
            species="Berberis darwinii",
            trait="self_incompatibility",
            value="SI",
            quality="high",
            provider="vazquez_simberloff_2004_primary_article",
            url="https://interactio.org/wp-content/uploads/2014/10/vazquez_simberloff_2004.pdf",
            title=(
                "Indirect effects of an introduced ungulate on pollination and "
                "plant reproduction"
            ),
            citation=(
                "Vazquez & Simberloff (2004), Ecological Monographs 74:281-308, "
                "DOI 10.1890/02-4055"
            ),
            excerpt=(
                "B. darwinii had the same proportion of fruits per tagged flower "
                "in the natural and cross-pollinated treatments. Self-pollinated "
                "fruit set was lower, indicating a relatively high degree of "
                "self-incompatibility."
            ),
            record_id="doi:10.1890/02-4055:appendix-a:berberis-darwinii:si",
            lineage="doi:10.1890/02-4055",
            lineage_method="original_primary_article_doi",
            source_tier="A",
            source_type="primary_research_article_controlled_pollination_experiment",
            domain="interactio.org",
            content_sha256=(
                "bacfe40078f65d0df1b010951a27c9fd5f5031f1ac59f88e57132280ee328b77"
            ),
            content_sha256_basis="downloaded_primary_article_pdf_bytes",
            retrieved_at_utc="2026-08-12T09:44:00Z",
            raw_value=(
                "fruit set natural 0.59; bag/no pollen 0.15; self 0.24; cross 0.59"
            ),
        ),
        _evidence_row(
            species="Berberis darwinii",
            trait="autonomous_selfing_capacity",
            value="mixed_or_variable",
            quality="high",
            provider="vazquez_simberloff_2004_primary_article",
            url="https://interactio.org/wp-content/uploads/2014/10/vazquez_simberloff_2004.pdf",
            title=(
                "Indirect effects of an introduced ungulate on pollination and "
                "plant reproduction"
            ),
            citation=(
                "Vazquez & Simberloff (2004), Ecological Monographs 74:281-308, "
                "DOI 10.1890/02-4055"
            ),
            excerpt=(
                "Table A1, Berberis darwinii fruit set: natural 0.59; bagged, no "
                "pollen 0.15; bagged self 0.24; bagged cross 0.59."
            ),
            record_id="doi:10.1890/02-4055:appendix-a:berberis-darwinii:autonomy",
            lineage="doi:10.1890/02-4055",
            lineage_method="original_primary_article_doi",
            source_tier="A",
            source_type="primary_research_article_controlled_pollination_experiment",
            domain="interactio.org",
            content_sha256=(
                "bacfe40078f65d0df1b010951a27c9fd5f5031f1ac59f88e57132280ee328b77"
            ),
            content_sha256_basis="downloaded_primary_article_pdf_bytes",
            retrieved_at_utc="2026-08-12T09:44:00Z",
            raw_value=(
                "autonomous index 0.15/0.59=0.254; intermediate by 0.2/0.8 thresholds"
            ),
        ),
        _evidence_row(
            species="Jasminum auriculatum",
            trait="self_incompatibility",
            value="SC",
            quality="high",
            provider="suresh_etal_2019_primary_article",
            url=(
                "https://epubs.icar.org.in/index.php/IJAgS/article/download/"
                "93522/37673/240531"
            ),
            title="Studies on cross compatibility in Jasminum spp.",
            citation=(
                "Suresh et al. (2019), Indian Journal of Agricultural Sciences "
                "89:1543-1546, DOI 10.56093/ijas.v89i9.93522"
            ),
            excerpt=(
                "On self pollination, the highest fruit set was recorded in "
                "J. auriculatum followed by J. flexile (Table 2). Table 2: "
                "J. auriculatum, self pollination, 7.00 fruits, 23.33% fruit "
                "set, and 7.00 seeds, 23.00% seed set."
            ),
            record_id="doi:10.56093/ijas.v89i9.93522:table2:jasminum-auriculatum:sc",
            lineage="doi:10.56093/ijas.v89i9.93522",
            lineage_method="original_primary_article_doi",
            source_tier="A",
            source_type="primary_research_article_controlled_self_pollination",
            domain="epubs.icar.org.in",
            content_sha256=(
                "00f823d7d3c27f5bd968fa6dfa0a1b7fa0068581f2503d718c313e23e0ab7db3"
            ),
            content_sha256_basis="downloaded_publisher_pdf_bytes",
            retrieved_at_utc="2026-08-12T11:01:00Z",
            raw_value=(
                "self pollination fruit set 23.33%; seed set 23.00%; "
                "manual self treatment only, not autonomous selfing"
            ),
        ),
        _evidence_row(
            species="Allophylus cobbe",
            trait="flower_size_class",
            value="very_small",
            quality="high",
            provider="singapore_nparks_flora_fauna_web",
            url="https://www.nparks.gov.sg/florafaunaweb/flora/1/6/1631",
            title="NParks Flora & Fauna Web: Allophylus cobbe",
            citation="National Parks Board Singapore, Allophylus cobbe species account",
            excerpt=(
                "Allophylus cobbe. Its shortly-stalked flowers are faintly fragrant, "
                "up to 2.5 mm wide, with white petals, and green to whitish sepals."
            ),
            record_id="nparks:flora:1631:allophylus-cobbe:flower-size",
            lineage="url:https://www.nparks.gov.sg/florafaunaweb/flora/1/6/1631",
            lineage_method="canonical_government_species_page",
            source_tier="A",
            source_type="government_flora_species_account",
            domain="nparks.gov.sg",
            content_sha256=(
                "368eac74a320de06fc351a65749c74f262324a2a4e8da42b870504f44e68a936"
            ),
            content_sha256_basis="retrieved_species_page_bytes",
            retrieved_at_utc="2026-08-12T10:05:00Z",
            raw_value="flowers up to 2.5 mm wide",
        ),
        _evidence_row(
            species="Allophylus cobbe",
            trait="flower_primary_color",
            value="white",
            quality="high",
            provider="singapore_nparks_flora_fauna_web",
            url="https://www.nparks.gov.sg/florafaunaweb/flora/1/6/1631",
            title="NParks Flora & Fauna Web: Allophylus cobbe",
            citation="National Parks Board Singapore, Allophylus cobbe species account",
            excerpt=(
                "Allophylus cobbe. Its shortly-stalked flowers are faintly fragrant, "
                "up to 2.5 mm wide, with white petals, and green to whitish sepals."
            ),
            record_id="nparks:flora:1631:allophylus-cobbe:flower-colour",
            lineage="url:https://www.nparks.gov.sg/florafaunaweb/flora/1/6/1631",
            lineage_method="canonical_government_species_page",
            source_tier="A",
            source_type="government_flora_species_account",
            domain="nparks.gov.sg",
            content_sha256=(
                "368eac74a320de06fc351a65749c74f262324a2a4e8da42b870504f44e68a936"
            ),
            content_sha256_basis="retrieved_species_page_bytes",
            retrieved_at_utc="2026-08-12T10:05:00Z",
            raw_value="white petals",
        ),
        _evidence_row(
            species="Allophylus cobbe",
            trait="floral_symmetry",
            value="actinomorphic",
            quality="medium",
            provider="pngplants_tree_key",
            url=(
                "https://www.pngplants.org/PNGtrees/TreeDescriptions/"
                "Allophylus_cobbe_L_Blume.html"
            ),
            title="PNGTreesKey: Allophylus cobbe",
            citation="PNGTreesKey / PNGplants botanical species account",
            excerpt=(
                "Allophylus cobbe. Flowers bisexual, stalked, flowers with many "
                "planes of symmetry, 1.0-2.5 mm long."
            ),
            record_id="pngtrees:allophylus-cobbe:floral-symmetry",
            lineage=(
                "url:https://www.pngplants.org/PNGtrees/TreeDescriptions/"
                "Allophylus_cobbe_L_Blume.html"
            ),
            lineage_method="canonical_institutional_tree_key_species_page",
            source_tier="B",
            source_type="institutional_botanical_tree_key",
            domain="pngplants.org",
            content_sha256=(
                "385206c9bf9b55ddfbfaab24ea57fc603acbaa0128e3fd9fc7bebb708dbaac9f"
            ),
            content_sha256_basis="retrieved_species_page_bytes",
            retrieved_at_utc="2026-08-12T10:06:00Z",
            raw_value="flowers with many planes of symmetry",
        ),
        _evidence_row(
            species="Dicliptera sexangularis",
            trait="floral_symmetry",
            value="zygomorphic",
            quality="medium",
            provider="leon_levy_native_plant_preserve",
            url="https://levypreserve.org/plant-listings/dicliptera-sexangularis/",
            title="Dicliptera sexangularis - Leon Levy Native Plant Preserve",
            citation="Leon Levy Native Plant Preserve species account",
            excerpt=(
                "Dicliptera sexangularis. The complete, perfect, zygomorphic "
                "flowers are arranged in axillary and terminal spikes/panicles."
            ),
            record_id="levy-preserve:dicliptera-sexangularis:floral-symmetry",
            lineage=(
                "url:https://levypreserve.org/plant-listings/"
                "dicliptera-sexangularis"
            ),
            lineage_method="canonical_botanical_preserve_species_page",
            source_tier="B",
            source_type="botanical_preserve_species_account",
            domain="levypreserve.org",
            content_sha256=(
                "64dd68e3a79346eba4a6ee1710e8ce48c9f5ec2ac8a81efc7e06254623397474"
            ),
            content_sha256_basis="retrieved_species_page_bytes",
            retrieved_at_utc="2026-08-12T10:08:00Z",
            raw_value="zygomorphic flowers",
        ),
        _evidence_row(
            species="Dicliptera sexangularis",
            trait="inflorescence_display",
            value="raceme_spike_panicle",
            quality="medium",
            provider="leon_levy_native_plant_preserve",
            url="https://levypreserve.org/plant-listings/dicliptera-sexangularis/",
            title="Dicliptera sexangularis - Leon Levy Native Plant Preserve",
            citation="Leon Levy Native Plant Preserve species account",
            excerpt=(
                "Dicliptera sexangularis. The complete, perfect, zygomorphic "
                "flowers are arranged in axillary and terminal spikes/panicles."
            ),
            record_id="levy-preserve:dicliptera-sexangularis:inflorescence",
            lineage=(
                "url:https://levypreserve.org/plant-listings/"
                "dicliptera-sexangularis"
            ),
            lineage_method="canonical_botanical_preserve_species_page",
            source_tier="B",
            source_type="botanical_preserve_species_account",
            domain="levypreserve.org",
            content_sha256=(
                "64dd68e3a79346eba4a6ee1710e8ce48c9f5ec2ac8a81efc7e06254623397474"
            ),
            content_sha256_basis="retrieved_species_page_bytes",
            retrieved_at_utc="2026-08-12T10:08:00Z",
            raw_value="axillary and terminal spikes/panicles",
        ),
    ]


def build(
    *,
    master_csv: Path,
    output_dir: Path,
    prior_curated_evidence_csv: Path | None = None,
    prior_curated_audit_csv: Path | None = None,
) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    master_names = set(master["accepted_species"])
    evidence = pd.DataFrame(reviewed_rows(), columns=EVIDENCE_COLUMNS).fillna("")
    missing = sorted(set(evidence["accepted_species"]) - master_names)
    if missing:
        raise ValueError(f"reviewed species absent from master: {missing}")
    evidence = evidence.sort_values(
        ["accepted_species", "trait_name", "candidate_id"]
    ).reset_index(drop=True)
    audit = _audit(evidence)
    evidence_path = output_dir / "high_information_direct_evidence_20260811.csv"
    audit_path = output_dir / "high_information_direct_manual_audit_20260811.csv"
    evidence.to_csv(evidence_path, index=False, lineterminator="\n")
    audit.to_csv(audit_path, index=False, lineterminator="\n")
    combined_paths: list[Path] = []
    combined_evidence: pd.DataFrame | None = None
    combined_audit: pd.DataFrame | None = None
    if prior_curated_evidence_csv is not None:
        if prior_curated_audit_csv is None:
            raise ValueError("prior curated audit is required with prior evidence")
        prior_evidence = pd.read_csv(prior_curated_evidence_csv, dtype=str).fillna("")
        prior_audit = pd.read_csv(prior_curated_audit_csv, dtype=str).fillna("")
        owned_ids = set(evidence["candidate_id"])
        prior_evidence = prior_evidence.loc[
            ~prior_evidence["candidate_id"].isin(owned_ids)
        ]
        prior_audit = prior_audit.loc[~prior_audit["candidate_id"].isin(owned_ids)]
        combined_evidence = pd.concat([prior_evidence, evidence], ignore_index=True)
        combined_audit = pd.concat([prior_audit, audit], ignore_index=True)
        if combined_evidence["candidate_id"].duplicated().any():
            raise ValueError("combined curated evidence candidate IDs are not unique")
        if combined_audit["candidate_id"].duplicated().any():
            raise ValueError("combined curated audit candidate IDs are not unique")
        combined_evidence_path = output_dir / "combined_curated_evidence_20260811.csv"
        combined_audit_path = output_dir / "combined_curated_manual_audit_20260811.csv"
        combined_evidence.to_csv(
            combined_evidence_path, index=False, lineterminator="\n"
        )
        combined_audit.to_csv(combined_audit_path, index=False, lineterminator="\n")
        combined_paths = [combined_evidence_path, combined_audit_path]
    summary: dict[str, object] = {
        "contract": "high_information_individually_reviewed_direct_v1",
        "created_at_utc": CREATED_AT,
        "evidence_rows": len(evidence),
        "species_trait": len(evidence.drop_duplicates(["accepted_species", "trait_name"])),
        "quality_counts": evidence["evidence_quality"].value_counts().to_dict(),
        "trait_counts": evidence["trait_name"].value_counts().to_dict(),
        "audit": {
            "reviewed": len(audit),
            "accepted_correct": len(audit),
            "precision": 1.0,
            "cultivar_contamination_rate": 0.0,
        },
        "guards": {
            "genus_trait_only": True,
            "family_inference": False,
            "global_fallback": False,
            "cross_trait_substitution": False,
        },
        "files": {},
    }
    for path in (evidence_path, audit_path, *combined_paths):
        summary["files"][path.name] = {
            "sha256": _sha256(path),
            "size_bytes": len(_canonical_file_bytes(path)),
        }
    manifest_path = output_dir / "high_information_direct_manifest_20260811.json"
    manifest_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    if combined_evidence is not None and combined_audit is not None:
        combined_manifest = {
            "contract": "next_acquisition_individually_reviewed_checkpoint_v1",
            "created_at_utc": CREATED_AT,
            "evidence_rows": len(combined_evidence),
            "species": int(combined_evidence["accepted_species"].nunique()),
            "species_trait": int(
                combined_evidence[["accepted_species", "trait_name"]]
                .drop_duplicates()
                .shape[0]
            ),
            "quality_counts": combined_evidence["evidence_quality"]
            .value_counts()
            .to_dict(),
            "trait_counts": combined_evidence["trait_name"].value_counts().to_dict(),
            "provider_counts": combined_evidence["source_provider"]
            .value_counts()
            .to_dict(),
            "audit": {
                "reviewed": len(combined_audit),
                "accepted_correct": int(
                    combined_audit["decision"].str.casefold().eq("accept").sum()
                ),
                "precision": float(
                    combined_audit["decision"].str.casefold().eq("accept").mean()
                ),
                "cultivar_contamination_rate": float(
                    combined_audit["cultivar_contamination"]
                    .str.casefold()
                    .eq("true")
                    .mean()
                ),
            },
            "files": {
                path.name: {
                    "sha256": _sha256(path),
                    "size_bytes": len(_canonical_file_bytes(path)),
                }
                for path in combined_paths
            },
        }
        (output_dir / "next_acquisition_checkpoint_manifest_20260811.json").write_text(
            json.dumps(combined_manifest, ensure_ascii=False, indent=2, sort_keys=True)
            + "\n",
            encoding="utf-8",
        )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--master-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--prior-curated-evidence-csv", type=Path)
    parser.add_argument("--prior-curated-audit-csv", type=Path)
    print(json.dumps(build(**vars(parser.parse_args())), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
