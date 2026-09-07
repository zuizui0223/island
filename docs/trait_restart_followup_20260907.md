# Reproductive source follow-up

Baseline remains recovered Run 34093932899, artifact 10007872352. Reviewer: Codex, automated source screening on 2026-09-07. No completed new acquisition artifact, no accepted gain claimed.

## Callicarpa: useful evidence hidden behind a synonym

Source: Shiuchi and Maruno, 2017, Ex situ conservation of Callicarpa longissima, Bull. Bot. Gard. Toyama 23. Original PDF: https://www.bgtym.org/_wp/wp-content/uploads/2022/02/research23.pdf . Relevant printed page 42 (PDF page index 43), section on growth in cultivation. The authors report seed and seedling production from separately isolated, wild-origin plants, not a named horticultural cultivar. Short exact excerpt: 「個体別に隔離して栽培した高隈山地周辺のタカクマムラサキでも、同様に自殖により種子が得られ、実生を得ている。」

Candidate lineage: citation:Shiuchi-Maruno-2017-Bull-Bot-Gard-Toyama-23. Source name Callicarpa longissima is absent as a key in the baseline, but Callicarpa dolichophylla is present: structure High, colour unresolved, reproduction unresolved. GBIF strict match of Callicarpa longissima returns EXACT, species rank, SYNONYM, acceptedUsageKey 5609799; that accepted record is Callicarpa dolichophylla Merr., family Lamiaceae. These two API endpoints were read directly:

- https://api.gbif.org/v1/species/match?name=Callicarpa%20longissima&strict=true
- https://api.gbif.org/v1/species/5609799

WFO indexed result wfo-0000768673 lists the same synonym, but direct page fetch timed out. POWO original synonym record 861386-1 returned HTTP403. Therefore the strict second-backbone original-response receipt is NOT yet complete; search results are not substituted for that gate.

Candidate trait: autonomous_selfing_capacity, candidate value autonomous; review status pending_second_backbone_and_admission_review. Keep original observational conditions and distinguish this from controlled bagging or self-compatibility experiments. Promotion=false; genus_rule_training=false. This source is also potential counterevidence against extrapolating autonomous-selfing absence throughout Callicarpa, not a third agreeing vote for the existing two-species absent pattern. Do not silently invalidate existing Low before a trait-specific review.

## Other inspected leads and limits

- Euonymus chloranthoides: Chinese-language original research page https://html.rhhz.net/linyekexue/html/20070507.htm cites Zhang et al. 2006 (DOI 10.3969/j.issn.1000-3142.2006.03.017) for compatibility and pollinator requirement. The exact species key is not in the baseline; primary experiment and identity are unresolved. No direct gain or mating-system substitution.
- Sideroxylon obtusifolium: a 2015 study at https://www.scielo.br/j/rarv/a/Tms99bXBb8qbFWwVwLSNTmL/?lang=pt reports low hand-self fruiting, higher cross fruiting, and no spontaneous selfing. Discovery remains unadmitted; do not convert partial fruiting into a conflict-resolved SC value or replace the existing SI record from a different study without review.
- Callicarpa americana: the located 2014 source already supplies its SC record. Cultivar colour segregation cannot be transferred to wild species colours or used as evidence for autonomous selfing.
- Lithocarpus pollen-tube/arrest papers do not by themselves establish self-incompatibility. No values were inferred from the mechanism or from article-title matches.

## Next action

Complete the Callicarpa synonym receipt through an available official backbone endpoint, preserve the species-level paper quote, and run admission against the baseline and full source lineage. In parallel with subsequent source work, prioritize the existing reviewed Wave56-58 collision audit over repeatedly searching genera with only abstracts or duplicate records. Two consecutive zero-gain completed batches must stop a lane; the current screen is unfinished, not an excuse to relaunch a completed global search.
