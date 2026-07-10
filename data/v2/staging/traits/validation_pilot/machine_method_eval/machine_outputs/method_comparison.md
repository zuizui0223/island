# Machine Method Evaluation

Human review was not required for this run. Outputs are machine layers, not human-adjudicated truth.

## Recommendation

- Use `web_reported_species_direct` as the current main no-review trait lane: 23 rows across 12 species.
- Use `globi_interaction_claims` as the current no-review interaction lane: 21 records across 5 taxa.
- Keep `rule_based_trait_proxy` and indirect/descriptive text as sensitivity layers.
- Build the next LLM extractor for reproductive mode, selfing, and visitor statements; do not use images for those.
- Use images later only for visible floral traits, because image acquisition/inference is heavier.

## Output Counts

- `machine_trait_selected.csv`: 12 selected species-trait rows
- `machine_trait_sensitivity.csv`: 22 sensitivity/proxy rows
- `machine_interaction_claims.csv`: 21 interaction-claim rows

## Method Table

| method | n_rows | n_species | species_coverage_rate | n_valid_machine_values | n_conflict_species_trait_groups | recommended_role | reason |
| --- | --- | --- | --- | --- | --- | --- | --- |
| web_reported_species_direct | 23 | 12 | 0.24 | 23 | 6 | best_current_no_review_trait_lane | Explicit species-page text with source URLs/excerpts; highest precision among current trait methods. |
| web_reported_species_indirect | 4 | 3 | 0.06 | 4 | 1 | sensitivity_only | Useful lead, but genus/indirect pages are weaker without manual checking. |
| gbif_descriptive_scout | 1 | 1 | 0.02 | 1 | 0 | not_recommended_as_main | Very sparse in the validation pilot and prone to description-context noise. |
| rule_based_trait_proxy | 17 | 17 | 0.34 | 17 | 0 | proxy_sensitivity_only | Expands coverage but is phenotype/family inference, not source-stated trait evidence. |
| globi_interaction_claims | 21 | 5 | 0.1 | 21 | 0 | best_current_no_review_interaction_lane | Source-backed pollination/flower-visit claims; sparse and not effectiveness or guild by itself. |
| visual_image_candidate | 0 | 0 | 0.0 | 0 | 0 | later_auxiliary_visible_traits | Potentially useful for visible M0 traits, but heavier and no local image manifest was available in this run. |
| llm_reported_ecology_excerpt | 0 | 0 | 0.0 | 0 | 0 | next_method_to_build_for_reproduction_and_visitors | Best fit for selfing/reproductive-mode/visitor statements if the output preserves source excerpts. |
