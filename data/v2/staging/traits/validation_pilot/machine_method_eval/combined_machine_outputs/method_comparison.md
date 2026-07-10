# Machine Method Evaluation

Human review was not required for this run. Outputs are machine layers, not human-adjudicated truth.

## Recommendation

- Use `web_reported_species_direct` as the current main no-review trait lane: 23 rows across 12 species.
- Use `globi_interaction_claims` as the current no-review interaction lane: 21 records across 5 taxa.
- Use `machine_pollinator_guild_index` as a derived no-review visitor guild layer: 4 species rows.
- Keep `rule_based_trait_proxy` and indirect/descriptive text as sensitivity layers.
- Build the next LLM extractor for reproductive mode, selfing, and visitor statements; do not use images for those.
- Use images later only for visible floral traits, because image acquisition/inference is heavier.

## Output Counts

- `machine_trait_selected.csv`: 12 selected species-trait rows
- `machine_trait_sensitivity.csv`: 25 sensitivity/proxy rows
- `machine_interaction_claims.csv`: 21 interaction-claim rows
- `machine_pollinator_guild_index.csv`: 4 pollinator-guild species rows

## Method Table

| method | n_rows | n_species | species_coverage_rate | n_valid_machine_values | n_conflict_species_trait_groups | recommended_role | reason |
| --- | --- | --- | --- | --- | --- | --- | --- |
| web_reported_species_direct | 23 | 12 | 0.24 | 23 | 6 | best_current_no_review_trait_lane | Explicit species-page text with source URLs/excerpts; highest precision among current trait methods. |
| web_reported_species_indirect | 4 | 3 | 0.06 | 4 | 1 | sensitivity_only | Useful lead, but genus/indirect pages are weaker without manual checking. |
| gbif_descriptive_scout | 1 | 1 | 0.02 | 1 | 0 | not_recommended_as_main | Very sparse in the validation pilot and prone to description-context noise. |
| rule_based_trait_proxy | 17 | 17 | 0.34 | 17 | 0 | proxy_sensitivity_only | Expands coverage but is phenotype/family inference, not source-stated trait evidence. |
| llm_reported_ecology_excerpt | 3 | 1 | 0.02 | 3 | 0 | next_main_text_lane_for_reproduction_and_visitors | Source-excerpted reproductive, selfing, pollen-vector, or visitor statements; usable as main no-review text evidence when species-direct. |
| visual_image_candidate | 0 | 0 | 0.0 | 0 | 0 | later_auxiliary_visible_traits | Useful for visible M0 traits only; heavier than text and not suitable for reproductive or visitor states. |
| globi_interaction_claims | 21 | 5 | 0.1 | 21 | 0 | best_current_no_review_interaction_lane | Source-backed pollination/flower-visit claims; sparse and not effectiveness or guild by itself. |
| machine_pollinator_guild_index | 4 | 4 | 0.08 | 4 | 0 | derived_no_review_pollinator_guild_layer | Derived from source-backed GloBI partner taxa; use as a machine guild index, not human-verified ecology. |
