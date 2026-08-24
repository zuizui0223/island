"""Freeze reviewed visible-flower evidence targeted at support-2 genus rules.

Every row is species-direct and backed by the retrieved original page or PDF.
This checkpoint never emits genus inference; the shared all-evidence rebuild
decides whether a genus x trait_name rule passes dominance and masked tests.
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

CREATED_AT = "2026-08-13T15:49:21Z"
SOURCE_GROUP = "visible_morphology_rule_checkpoint_20260813"

ROWS = [
    {
        "species": "Glyceria declinata",
        "trait": "inflorescence_display",
        "value": "raceme_spike_panicle",
        "quality": "high",
        "provider": "VicFlora, Royal Botanic Gardens Victoria",
        "url": "https://vicflora.rbg.vic.gov.au/flora/taxon/27aea429-98f2-448a-ac06-d82159e0a6ae",
        "title": "Glyceria declinata",
        "citation": "Walsh (1994), Poaceae, Flora of Victoria 2: 356-627",
        "excerpt": "Inflorescence a loose, slender panicle, often appearing one-sided, 4-40 cm long",
        "record_id": "vicflora:27aea429-98f2-448a-ac06-d82159e0a6ae:inflorescence",
        "lineage": "flora_treatment:flora_of_victoria:glyceria_declinata",
        "source_tier": "A",
        "source_type": "official_taxonomic_flora",
        "domain": "vicflora.rbg.vic.gov.au",
        "content_sha256": "c45a3be22b9fa84cde4a5a27301bfc615fe68d7ba0d66caed2a9c53a71c18798",
        "raw_value": "loose, slender panicle",
    },
    {
        "species": "Crambe maritima",
        "trait": "inflorescence_display",
        "value": "raceme_spike_panicle",
        "quality": "high",
        "provider": "Cambridge University Botanic Garden",
        "url": "https://www.botanic.cam.ac.uk/the-garden/plant-list/crambe-maritima/",
        "title": "Crambe maritima",
        "citation": "Cambridge University Botanic Garden species profile",
        "excerpt": "Stout, woody, edible stems can reach 75cm in height and bear large racemes of small, white, honey-scented, cruciform flowers.",
        "record_id": "cambridge-botanic-garden:crambe-maritima:racemes",
        "lineage": "botanic_garden_profile:cambridge:crambe_maritima",
        "source_tier": "A",
        "source_type": "university_botanic_garden_species_profile",
        "domain": "botanic.cam.ac.uk",
        "content_sha256": "d071809a10160a37cc8068d17a8aa38aa8d802d0cbbc55cdf9692badc0af8af5",
        "raw_value": "large racemes",
    },
    {
        "species": "Crambe maritima",
        "trait": "floral_symmetry",
        "value": "actinomorphic",
        "quality": "high",
        "provider": "Pladias, Database of the Czech Flora and Vegetation",
        "url": "https://pladias.cz/taxon/data/Crambe%20maritima",
        "title": "Crambe maritima - katran primorsky",
        "citation": "Pladias trait sourced to Kvetena Ceske republiky",
        "excerpt": "Symetrie kvetu: aktinomorfni (dve a vice rovin soumernosti)",
        "record_id": "pladias:crambe-maritima:flower-symmetry",
        "lineage": "provider_treatment:pladias:crambe_maritima",
        "source_tier": "A",
        "source_type": "national_flora_trait_database",
        "domain": "pladias.cz",
        "content_sha256": "710120556f9b77c7f18e52ed97d4c8d1f364f0623e5479e78b54a6c89824dcbf",
        "raw_value": "aktinomorfni",
    },
    {
        "species": "Herniaria glabra",
        "trait": "inflorescence_display",
        "value": "umbel_corymb",
        "quality": "high",
        "provider": "World Flora Online / Northeastern Flora",
        "url": "https://www.worldfloraonline.org/taxon/wfo-0000720677",
        "title": "Herniaria glabra L.",
        "citation": "Northeastern Flora local description distributed by World Flora Online",
        "excerpt": "fls congested in small cymose axillary (or seemingly lf-opposed) clusters",
        "record_id": "wfo:wfo-0000720677:northeastern-flora:cymose-clusters",
        "lineage": "provider_treatment:northeastern_flora:wfo-0000720677",
        "source_tier": "A",
        "source_type": "official_flora_species_treatment",
        "domain": "worldfloraonline.org",
        "content_sha256": "67ec059b4eeece8ecc40a3ad1497b38b3a80a90798c26037905b074203e1a1c1",
        "raw_value": "small cymose axillary clusters",
    },
    {
        "species": "Plumeria obtusa",
        "trait": "flower_size_class",
        "value": "large",
        "quality": "high",
        "provider": "World Flora Online / Flora of China",
        "url": "https://www.worldfloraonline.org/taxon/wfo-0000279164",
        "title": "Plumeria obtusa L.",
        "citation": "Flora of China species treatment distributed by World Flora Online",
        "excerpt": "Corolla white, ca. 4 cm in diam., throat yellow; lobes spreading, slightly recurved.",
        "record_id": "wfo:wfo-0000279164:flora-of-china:corolla-diameter",
        "lineage": "provider_treatment:flora_of_china:wfo-0000279164",
        "source_tier": "A",
        "source_type": "official_flora_species_treatment",
        "domain": "worldfloraonline.org",
        "content_sha256": "dfc5984706fe75dc00a8ac97f8d848004ebb6640032048afe5b56755e9d318f0",
        "raw_value": "ca. 4 cm in diameter",
    },
    {
        "species": "Plumeria obtusa",
        "trait": "tube_depth_class",
        "value": "intermediate",
        "quality": "high",
        "provider": "World Flora Online / e-Flora of Thailand",
        "url": "https://www.worldfloraonline.org/taxon/wfo-0000279164",
        "title": "Plumeria obtusa L.",
        "citation": "e-Flora of Thailand species treatment distributed by World Flora Online",
        "excerpt": "Corolla white with a yellow mouth; tube 1.7-2.3 cm long; lobes 2.3-6 cm long",
        "record_id": "wfo:wfo-0000279164:eflora-thailand:corolla-tube",
        "lineage": "provider_treatment:eflora_thailand:wfo-0000279164",
        "source_tier": "A",
        "source_type": "official_flora_species_treatment",
        "domain": "worldfloraonline.org",
        "content_sha256": "dfc5984706fe75dc00a8ac97f8d848004ebb6640032048afe5b56755e9d318f0",
        "raw_value": "1.7-2.3 cm",
    },
    {
        "species": "Plumeria obtusa",
        "trait": "inflorescence_display",
        "value": "umbel_corymb",
        "quality": "high",
        "provider": "University of Kerala Digital Garden",
        "url": "https://www.keralauniversity.ac.in/digital_garden/170.html",
        "title": "Plumeria obtusa L.",
        "citation": "University of Kerala Digital Garden species treatment",
        "excerpt": "Flowers bisexual, white, in terminal corymbose stout cymes",
        "record_id": "kerala-university-digital-garden:170:inflorescence",
        "lineage": "university_botanical_garden:kerala:plumeria_obtusa",
        "source_tier": "A",
        "source_type": "university_botanical_garden_species_profile",
        "domain": "keralauniversity.ac.in",
        "content_sha256": "4fe05c73a9bd780a4ab883536f29303647a15b0667d8696cdc67df3825dd6668",
        "raw_value": "terminal corymbose stout cymes",
    },
    {
        "species": "Odontochilus yakushimensis",
        "trait": "inflorescence_display",
        "value": "few_flowered",
        "quality": "medium",
        "provider": "Yamakei botanical species guide",
        "url": "https://ymkr-okurimono.sakura.ne.jp/plants/angiosperm/f70202_ran/yksm-hm-aridoshiran.html",
        "title": "Yakushimahimearidoshiran",
        "citation": "Photographic Japanese specialist botanical species page",
        "excerpt": "茎の先に1 - 3個の白色の花を穂状花序につける。",
        "record_id": "ymkr:odontochilus-yakushimensis:one-to-three-flowers",
        "lineage": "specialist_species_profile:ymkr:odontochilus_yakushimensis",
        "source_tier": "B",
        "source_type": "specialist_personal_botanical_species_guide",
        "domain": "ymkr-okurimono.sakura.ne.jp",
        "content_sha256": "8a3589d86e375fd5d2096222f0d02cbe70bd516052e9416751b599ea1a345412",
        "raw_value": "1-3 flowers",
    },
    {
        "species": "Sericolea arfakensis",
        "trait": "inflorescence_display",
        "value": "raceme_spike_panicle",
        "quality": "high",
        "provider": "Naturalis / Blumea",
        "url": "https://repository.naturalis.nl/pub/524701/BLUM1982028001006.pdf",
        "title": "A revision of Sericolea Schlechter (Elaeocarpaceae)",
        "citation": "van Balgooy (1982), Blumea 28(1): 103-141",
        "excerpt": "Inflorescence an axillary raceme up to 3 cm long, bearing 2-4 flowers",
        "record_id": "blumea:1982:28:1:sericolea-arfakensis:inflorescence",
        "lineage": "publication:blumea:1982:28:1:sericolea_revision",
        "source_tier": "A",
        "source_type": "primary_taxonomic_revision",
        "domain": "repository.naturalis.nl",
        "content_sha256": "51918f6482176542a36e1395fbee5d1eda5508ef998ac0325446b0eecb9fb9af",
        "raw_value": "axillary raceme bearing 2-4 flowers",
    },
    {
        "species": "Sericolea arfakensis",
        "trait": "flower_primary_color",
        "value": "red_pink",
        "quality": "high",
        "provider": "Naturalis / Blumea",
        "url": "https://repository.naturalis.nl/pub/524701/BLUM1982028001006.pdf",
        "title": "A revision of Sericolea Schlechter (Elaeocarpaceae)",
        "citation": "van Balgooy (1982), Blumea 28(1): 103-141",
        "excerpt": "petals oblong obovate, 2.5-4 X 1-2 mm, base cuneate, apex truncate, weakly 3-lobed, glabrous except for a small tuft of hairs inside at base, pink or red in vivo",
        "record_id": "blumea:1982:28:1:sericolea-arfakensis:petal-colour",
        "lineage": "publication:blumea:1982:28:1:sericolea_revision",
        "source_tier": "A",
        "source_type": "primary_taxonomic_revision",
        "domain": "repository.naturalis.nl",
        "content_sha256": "51918f6482176542a36e1395fbee5d1eda5508ef998ac0325446b0eecb9fb9af",
        "raw_value": "petals pink or red in vivo",
    },
    {
        "species": "Aporosa octandra",
        "trait": "floral_symmetry",
        "value": "actinomorphic",
        "quality": "medium",
        "provider": "Shiu Ying Hu Herbarium, Chinese University of Hong Kong",
        "url": "https://syhuherbarium.sls.cuhk.edu.hk/collections/factsheet-pro/aporusa-dioica/",
        "title": "Aporosa octandra var. octandra",
        "citation": "CUHK Shiu Ying Hu Herbarium Pro-Factsheet; references Flora of Hong Kong, Flora of China and FRPS",
        "excerpt": "Flower symmetry: Radial (Actinomorphic)",
        "record_id": "cuhk-herbarium:aporosa-octandra:flower-symmetry",
        "lineage": "provider_treatment:cuhk_shiu_ying_hu_herbarium:aporosa_octandra",
        "source_tier": "B",
        "source_type": "university_herbarium_species_profile",
        "domain": "syhuherbarium.sls.cuhk.edu.hk",
        "content_sha256": "38c8a9fb63fb35a78ad7934e3784d84e9fbb88e7ca6ac233d60d6e3d35ecd740",
        "raw_value": "Radial (Actinomorphic)",
    },
    {
        "species": "Aporosa octandra",
        "trait": "inflorescence_display",
        "value": "raceme_spike_panicle",
        "quality": "medium",
        "provider": "Shiu Ying Hu Herbarium, Chinese University of Hong Kong",
        "url": "https://syhuherbarium.sls.cuhk.edu.hk/collections/factsheet-pro/aporusa-dioica/",
        "title": "Aporosa octandra var. octandra",
        "citation": "CUHK Shiu Ying Hu Herbarium Pro-Factsheet; references Flora of Hong Kong, Flora of China and FRPS",
        "excerpt": "Inflorescences - Type: Spike",
        "record_id": "cuhk-herbarium:aporosa-octandra:inflorescence",
        "lineage": "provider_treatment:cuhk_shiu_ying_hu_herbarium:aporosa_octandra",
        "source_tier": "B",
        "source_type": "university_herbarium_species_profile",
        "domain": "syhuherbarium.sls.cuhk.edu.hk",
        "content_sha256": "38c8a9fb63fb35a78ad7934e3784d84e9fbb88e7ca6ac233d60d6e3d35ecd740",
        "raw_value": "Spike",
    },
    {
        "species": "Aporosa octandra",
        "trait": "flower_primary_color",
        "value": "yellow_orange",
        "quality": "medium",
        "provider": "Shiu Ying Hu Herbarium, Chinese University of Hong Kong",
        "url": "https://syhuherbarium.sls.cuhk.edu.hk/collections/factsheet-pro/aporusa-dioica/",
        "title": "Aporosa octandra var. octandra",
        "citation": "CUHK Shiu Ying Hu Herbarium Pro-Factsheet; references Flora of Hong Kong, Flora of China and FRPS",
        "excerpt": "Flowers (general) - Major Color: Yellow",
        "record_id": "cuhk-herbarium:aporosa-octandra:flower-colour",
        "lineage": "provider_treatment:cuhk_shiu_ying_hu_herbarium:aporosa_octandra",
        "source_tier": "B",
        "source_type": "university_herbarium_species_profile",
        "domain": "syhuherbarium.sls.cuhk.edu.hk",
        "content_sha256": "38c8a9fb63fb35a78ad7934e3784d84e9fbb88e7ca6ac233d60d6e3d35ecd740",
        "raw_value": "Yellow",
    },
    {
        "species": "Angelica gigas",
        "trait": "floral_form",
        "value": "open_radial",
        "quality": "medium",
        "provider": "Kwantlen Polytechnic University Plant Database",
        "url": "https://plantdatabase.kpu.ca/Plant/angi",
        "title": "Angelica gigas - giant angelica",
        "citation": "KPU School of Horticulture Plant Database species profile",
        "excerpt": "Corolla Shape: Rotate/stellate",
        "record_id": "kpu-plant-database:angi:corolla-shape",
        "lineage": "university_horticulture_profile:kpu:angelica_gigas",
        "source_tier": "B",
        "source_type": "university_horticulture_species_profile",
        "domain": "plantdatabase.kpu.ca",
        "content_sha256": "29c023a7be06d150faacfc529d4daa6666b89fac9977eac86fbac7506139ae9d",
        "raw_value": "Rotate/stellate",
    },
    {
        "species": "Angelica gigas",
        "trait": "inflorescence_display",
        "value": "umbel_corymb",
        "quality": "medium",
        "provider": "Kwantlen Polytechnic University Plant Database",
        "url": "https://plantdatabase.kpu.ca/Plant/angi",
        "title": "Angelica gigas - giant angelica",
        "citation": "KPU School of Horticulture Plant Database species profile",
        "excerpt": "Inflorescence Type: Umbel",
        "record_id": "kpu-plant-database:angi:inflorescence",
        "lineage": "university_horticulture_profile:kpu:angelica_gigas",
        "source_tier": "B",
        "source_type": "university_horticulture_species_profile",
        "domain": "plantdatabase.kpu.ca",
        "content_sha256": "29c023a7be06d150faacfc529d4daa6666b89fac9977eac86fbac7506139ae9d",
        "raw_value": "Umbel",
    },
    {
        "species": "Eumachia straminea",
        "trait": "tube_depth_class",
        "value": "shallow",
        "quality": "high",
        "provider": "PhytoKeys / PubMed Central",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12048816/",
        "title": "Eumachia brevipedunculata, a new species from Yunnan, China",
        "citation": "Quan & Tan (2025), PhytoKeys 255: 235-246, doi:10.3897/phytokeys.255.143380",
        "excerpt": "Morphological comparison of Eumachia brevipedunculata and E. straminea: Corolla tube length ... 1.5-2 mm",
        "record_id": "doi:10.3897/phytokeys.255.143380:eumachia-straminea:corolla-tube",
        "lineage": "doi:10.3897/phytokeys.255.143380",
        "source_tier": "A",
        "source_type": "primary_taxonomic_publication",
        "domain": "pmc.ncbi.nlm.nih.gov",
        "content_sha256": "a9fee5b29d48252de678cbd30ebdc1ededed0277c2ee0e4f25ea8be0c2b4a24c",
        "raw_value": "corolla tube 1.5-2 mm",
    },
    {
        "species": "Micromeria imbricata",
        "trait": "tube_depth_class",
        "value": "shallow|intermediate",
        "quality": "high",
        "provider": "South African National Biodiversity Institute",
        "url": "https://biodiversityadvisor.sanbi.org/search/detail/404f6c17-d3e9-49f6-9a4e-933b5f5b9369",
        "title": "Micromeria imbricata (Forssk.) C.Chr.",
        "citation": "Ryding (2007), Botanical Journal of the Linnean Society 155: 427-446, doi:10.1111/j.1095-8339.2007.00702.x",
        "excerpt": "Corolla strongly two-lipped, 3-12 mm long; tube tubular in the proximal part and funnel-shaped in the distal part, 2-7.5 mm long",
        "record_id": "sanbi:micromeria-imbricata:corolla-tube",
        "lineage": "doi:10.1111/j.1095-8339.2007.00702.x",
        "source_tier": "A",
        "source_type": "primary_taxonomic_revision_republished_by_government",
        "domain": "biodiversityadvisor.sanbi.org",
        "content_sha256": "c6375cc9c62ad1b7a6e1ebc9c2b7ba230dc05d68e6ba75b0165be81e93b4896e",
        "raw_value": "corolla tube 2-7.5 mm",
    },
    {
        "species": "Micromeria imbricata",
        "trait": "cleistogamy",
        "value": "facultative",
        "quality": "high",
        "provider": "South African National Biodiversity Institute",
        "url": "https://biodiversityadvisor.sanbi.org/search/detail/404f6c17-d3e9-49f6-9a4e-933b5f5b9369",
        "title": "Micromeria imbricata (Forssk.) C.Chr.",
        "citation": "Ryding (2007), Botanical Journal of the Linnean Society 155: 427-446, doi:10.1111/j.1095-8339.2007.00702.x",
        "excerpt": "Flowers chasmogamous and/or cleistogamous, mostly with both types in the same inflorescence.",
        "record_id": "sanbi:micromeria-imbricata:cleistogamy",
        "lineage": "doi:10.1111/j.1095-8339.2007.00702.x",
        "source_tier": "A",
        "source_type": "primary_taxonomic_revision_republished_by_government",
        "domain": "biodiversityadvisor.sanbi.org",
        "content_sha256": "c6375cc9c62ad1b7a6e1ebc9c2b7ba230dc05d68e6ba75b0165be81e93b4896e",
        "raw_value": "chasmogamous and/or cleistogamous; mostly both types",
    },
    {
        "species": "Capparis arborea",
        "trait": "flower_size_class",
        "value": "small",
        "quality": "medium",
        "provider": "Townsville State of the Environment species guide",
        "url": "https://www.soe-townsville.org/schoolsshadetrees/trees.html",
        "title": "School shade trees - Capparis arborea",
        "citation": "Townsville State of the Environment local species account",
        "excerpt": "Capparis arborea ... Flowers: Small, cream-green, 0.2-0.3 cm diameter; borne in small clusters on old wood.",
        "record_id": "townsville-soe:capparis-arborea:flower-diameter",
        "lineage": "local_government_species_profile:townsville:capparis_arborea",
        "source_tier": "B",
        "source_type": "local_government_environment_species_profile",
        "domain": "soe-townsville.org",
        "content_sha256": "a398288b8d0bbe99e55f99241c1369144eb2491aff3a9094b901c192421dcafc",
        "raw_value": "0.2-0.3 cm diameter",
    },
    {
        "species": "Colea sytsmae",
        "trait": "flower_primary_color",
        "value": "white",
        "quality": "high",
        "provider": "Annales Botanici Fennici / Finnish Zoological and Botanical Publishing Board",
        "url": "https://web.archive.org/web/20160305021040id_/http://www.sekj.org/PDF/anb43-free/anb43-225.pdf",
        "title": "New taxa of Coleeae (Bignoniaceae) from Madagascar. I. A collection from Masoala Peninsula",
        "citation": "Zjhra (2006), Annales Botanici Fennici 43: 225-239",
        "excerpt": "Corolla campanulate, medium-sized (6–11 mm length, 4–8 mm width), white, throat without nectar guides; corolla lobes mostly patent, pubescent.",
        "record_id": "zjhra-2006:colea-sytsmae:corolla-colour",
        "lineage": "publication:zjhra2006_coleeae_masoala",
        "source_tier": "A",
        "source_type": "primary_taxonomic_publication",
        "domain": "web.archive.org",
        "content_sha256": "ac7e7119a615cf3c2876b57ecb5691671792cf924dc8871fb653da18b50d2c6a",
        "raw_value": "corolla white",
    },
    {
        "species": "Ziziphus spina-christi",
        "trait": "floral_symmetry",
        "value": "actinomorphic",
        "quality": "high",
        "provider": "Frontiers in Plant Science",
        "url": "https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2019.01315/full",
        "title": "Pollinator Behavior Drives Sexual Specializations in the Hermaphrodite Flowers of a Heterodichogamous Tree",
        "citation": "Wajnberg et al. (2019), Frontiers in Plant Science 10:1315, doi:10.3389/fpls.2019.01315",
        "excerpt": "Flowers are small, radially symmetrical, and lack a corolla tube.",
        "record_id": "doi:10.3389/fpls.2019.01315:ziziphus-spina-christi:symmetry",
        "lineage": "doi:10.3389/fpls.2019.01315",
        "source_tier": "A",
        "source_type": "primary_field_study",
        "domain": "frontiersin.org",
        "content_sha256": "5b83b7756c53745f8c92677e69327b4e271222605df7fd8babdc415cd3d65a3a",
        "raw_value": "radially symmetrical",
    },
    {
        "species": "Emilia sonchifolia",
        "trait": "floral_symmetry",
        "value": "actinomorphic",
        "quality": "high",
        "provider": "Journal of BioScience and Biotechnology / University of Plovdiv",
        "url": "https://doi.uni-plovdiv.bg/bitstreams/80be634c-dfb1-40ec-84df-329645fdc227/download",
        "title": "Pump mechanism, secondary pollen presentation, psychophily and anemochory in Emilia sonchifolia",
        "citation": "Medabalimi, Aluri & Kunuku (2017), Journal of BioScience and Biotechnology 6(2):129-137, doi:10.69085/jbb20172129",
        "excerpt": "The disc florets are small (9.1 ± 1.2 mm long, 1.6 ± 0.4 mm wide), pink or purplish, odourless, actinomorphic, bisexual and nectariferous.",
        "record_id": "doi:10.69085/jbb20172129:emilia-sonchifolia:symmetry",
        "lineage": "doi:10.69085/jbb20172129",
        "source_tier": "A",
        "source_type": "primary_field_study",
        "domain": "doi.uni-plovdiv.bg",
        "content_sha256": "d1b6cb1d0bf09f6909640450637357d92e8cdb62352f4eca81477fcb0c9c7837",
        "raw_value": "disc florets actinomorphic",
    },
    {
        "species": "Emilia coccinea",
        "trait": "floral_symmetry",
        "value": "actinomorphic",
        "quality": "high",
        "provider": "PROTA species monograph / Pl@ntUse",
        "url": "https://plantuse.plantnet.org/en/Emilia_coccinea_%28PROTA%29",
        "title": "Emilia coccinea (PROTA)",
        "citation": "PROTA4U species account for Emilia coccinea",
        "excerpt": "Flowers bisexual, regular, 5-merous; corolla tubular, 5–9.5 mm long, yellow to orange or red (sometimes scarlet-red in ornamental cultivars).",
        "record_id": "prota:emilia-coccinea:flower-symmetry",
        "lineage": "provider_treatment:prota:Emilia coccinea",
        "source_tier": "A",
        "source_type": "expert_species_monograph",
        "domain": "plantuse.plantnet.org",
        "content_sha256": "c6bd430250f6f8d3abdbec5aeea376a810c8181df4dab38f901ec0171df509d0",
        "raw_value": "flowers regular",
    },
    {
        "species": "Emilia coccinea",
        "trait": "flower_size_class",
        "value": "small",
        "quality": "high",
        "provider": "PROTA species monograph / Pl@ntUse",
        "url": "https://plantuse.plantnet.org/en/Emilia_coccinea_%28PROTA%29",
        "title": "Emilia coccinea (PROTA)",
        "citation": "PROTA4U species account for Emilia coccinea",
        "excerpt": "Flowers bisexual, regular, 5-merous; corolla tubular, 5–9.5 mm long, yellow to orange or red (sometimes scarlet-red in ornamental cultivars).",
        "record_id": "prota:emilia-coccinea:corolla-length",
        "lineage": "provider_treatment:prota:Emilia coccinea",
        "source_tier": "A",
        "source_type": "expert_species_monograph",
        "domain": "plantuse.plantnet.org",
        "content_sha256": "c6bd430250f6f8d3abdbec5aeea376a810c8181df4dab38f901ec0171df509d0",
        "raw_value": "corolla 5-9.5 mm long",
    },
    {
        "species": "Paris quadrifolia",
        "trait": "inflorescence_display",
        "value": "solitary",
        "quality": "high",
        "provider": "World Flora Online / Flora Helvetica",
        "url": "https://www.worldfloraonline.org/taxon/wfo-0000715208",
        "title": "Paris quadrifolia L.",
        "citation": "Flora Helvetica species treatment distributed by World Flora Online",
        "excerpt": "Fleur verdâtre, solitaire, terminale ; pédoncule dressé, long de 3–6 cm, au milieu du verticille.",
        "record_id": "wfo:wfo-0000715208:flora-helvetica:solitary-flower",
        "lineage": "provider_treatment:flora_helvetica:wfo-0000715208",
        "source_tier": "A",
        "source_type": "official_flora_species_treatment",
        "domain": "worldfloraonline.org",
        "content_sha256": "c6aaccf532b64aace025ec023e59c3da8b388b6b7f62bab924e741fd687148f8",
        "raw_value": "fleur solitaire, terminale",
    },
    {
        "species": "Chrysanthemum indicum",
        "trait": "floral_form",
        "value": "composite_head",
        "quality": "high",
        "provider": "World Flora Online / Flora of China",
        "url": "https://www.worldfloraonline.org/taxon/wfo-0000041005",
        "title": "Chrysanthemum indicum L.",
        "citation": "Flora of China species treatment distributed by World Flora Online",
        "excerpt": "Synflorescence a lax terminal flat-topped cyme. Capitula many or few.",
        "record_id": "wfo:wfo-0000041005:flora-of-china:capitula",
        "lineage": "provider_treatment:flora_of_china:wfo-0000041005",
        "source_tier": "A",
        "source_type": "official_flora_species_treatment",
        "domain": "worldfloraonline.org",
        "content_sha256": "4814b52c768bb1b91b6bb53ad3e727838606748d250d7362d84a6d90f2a99bfa",
        "raw_value": "capitula many or few",
    },
    {
        "species": "Nymphoides indica",
        "trait": "floral_symmetry",
        "value": "actinomorphic",
        "quality": "high",
        "provider": "World Flora Online / Manual de Plantas de Costa Rica",
        "url": "https://www.worldfloraonline.org/taxon/wfo-0000382493",
        "title": "Nymphoides indica (L.) Kuntze",
        "citation": "Manual de Plantas de Costa Rica species treatment distributed by World Flora Online",
        "excerpt": "Fls. bisexuales, actinomorfas, con pedicelo 20–110 mm.",
        "record_id": "wfo:wfo-0000382493:manual-costa-rica:actinomorphic",
        "lineage": "provider_treatment:manual_de_plantas_de_costa_rica:wfo-0000382493",
        "source_tier": "A",
        "source_type": "official_flora_species_treatment",
        "domain": "worldfloraonline.org",
        "content_sha256": "a3ef4f75a398c4f24383401c313a833f49f4a467a3b9b90f36fd60d38e8fc9da",
        "raw_value": "flores actinomorfas",
    },
    {
        "species": "Tournefortia gnaphalodes",
        "trait": "inflorescence_display",
        "value": "umbel_corymb",
        "quality": "high",
        "provider": "World Flora Online / Flora de Nicaragua",
        "url": "https://www.worldfloraonline.org/taxon/wfo-0000410543",
        "title": "Tournefortia gnaphalodes (L.) R. Br. ex Roem. & Schult.",
        "citation": "Flora de Nicaragua species treatment distributed by World Flora Online",
        "excerpt": "Inflorescencias cimas densas y esparcidamente ramificadas, más o menos terminales.",
        "record_id": "wfo:wfo-0000410543:flora-de-nicaragua:cymes",
        "lineage": "provider_treatment:flora_de_nicaragua:wfo-0000410543",
        "source_tier": "A",
        "source_type": "official_flora_species_treatment",
        "domain": "worldfloraonline.org",
        "content_sha256": "85a21a8dcf21ad6bc63eb2934cf1ccb891ae720a5f909e6b98aa5cae723a8721",
        "raw_value": "cimas densas, terminales",
    },
    {
        "species": "Chisocheton cumingianus",
        "trait": "floral_form",
        "value": "tubular",
        "quality": "high",
        "provider": "PROSEA species monograph / Pl@ntUse",
        "url": "https://plantuse.plantnet.org/en/Chisocheton_cumingianus_%28PROSEA%29",
        "title": "Chisocheton cumingianus (PROSEA)",
        "citation": "PROSEA species account for Chisocheton cumingianus",
        "excerpt": "Flowers unisexual or bisexual, tubular, nearly 2 cm long; calyx campanulate, 1-3 mm long; petals (3-)4(-5), spatulate, 1-2.5 cm × 2.5 mm, pale yellow to white.",
        "record_id": "prosea:chisocheton-cumingianus:flower-form",
        "lineage": "provider_treatment:prosea:Chisocheton cumingianus",
        "source_tier": "A",
        "source_type": "expert_species_monograph",
        "domain": "plantuse.plantnet.org",
        "content_sha256": "70602b5bf88e15300a259170d45fdffdd0f0c8b0abe032e681af8cd60dd90a82",
        "raw_value": "flowers tubular",
    },
    {
        "species": "Onosma visianii",
        "trait": "floral_form",
        "value": "tubular",
        "quality": "high",
        "provider": "Flora Europaea treatment reproduced in Molecules / PubMed Central",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6152756/",
        "title": "Identification of Onosma visianii Roots Extract and Purified Shikonin Derivatives as Potential Acaricidal Agents against Tetranychus urticae",
        "citation": "Tutin (1964), Flora Europaea; species description reproduced in Molecules (2017) 22:1002, doi:10.3390/molecules22061002",
        "excerpt": "Flowers are sympetalous and heterochlamydeous, with short pedicels and pale-yellow tubular corolla.",
        "record_id": "doi:10.3390/molecules22061002:onosma-visianii:corolla-form",
        "lineage": "flora_treatment:flora_europaea:onosma_visianii",
        "source_tier": "A",
        "source_type": "official_flora_species_treatment_reproduced_in_peer_reviewed_article",
        "domain": "pmc.ncbi.nlm.nih.gov",
        "content_sha256": "a567b722d489ef10888f195e5574bfc1be9d5c081ac84af00f4f890e93e04b84",
        "raw_value": "pale-yellow tubular corolla",
    },
]


def reviewed_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for item in ROWS:
        evidence = _evidence_row(
            species=item["species"],
            trait=item["trait"],
            value=item["value"],
            quality=item["quality"],
            provider=item["provider"],
            url=item["url"],
            title=item["title"],
            citation=item["citation"],
            excerpt=item["excerpt"],
            record_id=item["record_id"],
            lineage=item["lineage"],
            lineage_method="original_species_treatment_or_publication_lineage",
            source_tier=item["source_tier"],
            source_type=item["source_type"],
            domain=item["domain"],
            content_sha256=item["content_sha256"],
            content_sha256_basis="retrieved_original_page_or_pdf_bytes",
            retrieved_at_utc=CREATED_AT,
            raw_value=item["raw_value"],
        )
        evidence["source_group"] = SOURCE_GROUP
        evidence["query"] = "support2_visible_morphology_original_species_source"
        rows.append(evidence)
    return rows


def _sha256(path: Path) -> str:
    return hashlib.sha256(_canonical_file_bytes(path)).hexdigest()


def build(
    *,
    master_csv: Path,
    prior_curated_evidence_csv: Path,
    prior_curated_audit_csv: Path,
    output_dir: Path,
) -> dict[str, object]:
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    evidence = pd.DataFrame(reviewed_rows(), columns=EVIDENCE_COLUMNS).fillna("")
    missing = sorted(set(evidence["accepted_species"]) - set(master["accepted_species"]))
    if missing:
        raise ValueError(f"checkpoint species missing from master: {missing}")
    if evidence[["accepted_species", "trait_name"]].duplicated().any():
        raise ValueError("checkpoint species-trait pairs must be unique")
    if not evidence["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all():
        raise ValueError("every checkpoint row requires a 64-character content hash")

    audit = _audit(evidence)
    audit["reviewer"] = "Codex source-backed visible morphology audit"
    audit["reviewed_at_utc"] = CREATED_AT
    audit["decision_reason"] = (
        "Accepted after target-master identity, exact organ and trait statement, "
        "value ontology, source lineage and cultivar status were rechecked. "
        "No genus, family or cross-trait inference is emitted here."
    )

    prior_evidence = pd.read_csv(prior_curated_evidence_csv, dtype=str).fillna("")
    prior_audit = pd.read_csv(prior_curated_audit_csv, dtype=str).fillna("")
    owned = prior_evidence["source_group"].eq(SOURCE_GROUP)
    prior_owned_ids = set(prior_evidence.loc[owned, "candidate_id"])
    current_ids = set(evidence["candidate_id"])
    combined_evidence = pd.concat([prior_evidence.loc[~owned], evidence], ignore_index=True)
    combined_audit = pd.concat(
        [
            prior_audit.loc[
                ~prior_audit["candidate_id"].isin(prior_owned_ids | current_ids)
            ],
            audit,
        ],
        ignore_index=True,
    )
    for label, frame in (("evidence", combined_evidence), ("audit", combined_audit)):
        if frame["candidate_id"].duplicated().any():
            raise ValueError(f"combined {label} candidate IDs must be unique")

    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "evidence": output_dir / "visible_morphology_rule_evidence_20260813.csv",
        "audit": output_dir / "visible_morphology_rule_manual_audit_20260813.csv",
        "combined_evidence": output_dir / "combined_curated_evidence_20260813.csv",
        "combined_audit": output_dir / "combined_curated_manual_audit_20260813.csv",
    }
    evidence.to_csv(paths["evidence"], index=False, lineterminator="\n")
    audit.to_csv(paths["audit"], index=False, lineterminator="\n")
    combined_evidence.to_csv(paths["combined_evidence"], index=False, lineterminator="\n")
    combined_audit.to_csv(paths["combined_audit"], index=False, lineterminator="\n")

    summary: dict[str, object] = {
        "contract": "visible_morphology_support2_rule_checkpoint_v3",
        "created_at_utc": CREATED_AT,
        "accepted_rows": len(evidence),
        "species_trait": int(
            evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
        ),
        "quality_counts": evidence["evidence_quality"].value_counts().to_dict(),
        "trait_counts": evidence["trait_name"].value_counts().to_dict(),
        "source_lineages": int(evidence["source_lineage"].nunique()),
        "theoretical_rule_cells_before_formal_validation": 443,
        "theoretical_only_not_formal_coverage": True,
        "explicit_counterevidence": [
            "Sericolea arfakensis|flower_primary_color|red_pink",
            "Micromeria imbricata|tube_depth_class|shallow|intermediate",
            "Capparis arborea|flower_size_class|small",
        ],
        "combined": {
            "evidence_rows": len(combined_evidence),
            "audit_rows": len(combined_audit),
        },
        "guardrails": {
            "trait_specific_records": True,
            "genus_inference_emitted_here": False,
            "family_inference": False,
            "global_fallback": False,
            "cross_trait_substitution": False,
            "search_snippet_evidence": False,
        },
        "inputs": {
            "master_csv": {"path": str(master_csv), "sha256": _sha256(master_csv)},
            "prior_curated_evidence_csv": {
                "path": str(prior_curated_evidence_csv),
                "sha256": _sha256(prior_curated_evidence_csv),
            },
            "prior_curated_audit_csv": {
                "path": str(prior_curated_audit_csv),
                "sha256": _sha256(prior_curated_audit_csv),
            },
        },
        "files": {},
    }
    for path in paths.values():
        summary["files"][path.name] = {
            "sha256": _sha256(path),
            "size_bytes": len(_canonical_file_bytes(path)),
        }
    manifest_path = output_dir / "visible_morphology_rule_manifest_20260813.json"
    manifest_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--master-csv", type=Path, required=True)
    parser.add_argument("--prior-curated-evidence-csv", type=Path, required=True)
    parser.add_argument("--prior-curated-audit-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    print(json.dumps(build(**vars(parser.parse_args())), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
