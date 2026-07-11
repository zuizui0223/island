# v2 trait recovery benchmark

The 25-species validation pilot is treated as a fixed benchmark for acquisition methods. New recovery methods should be compared on the same species before global scaling.

Current benchmark sequence:

1. fixed-site public-source baseline;
2. evidence-locked LLM extraction over the frozen fixed-site corpus;
3. query-driven web discovery with GBIF synonym expansion;
4. evidence-locked LLM extraction over the dynamic corpus;
5. visual evidence lane for flower colour and gross floral form.

The objective is not to maximize raw row count. A method wins only when it increases species-level trait coverage while preserving species identity, source URLs, frozen evidence, ontology validation, and an audit trail.

For the first fixed-site pilot, the baseline reached 7/25 species for any reported trait and the evidence-locked LLM lane increased this to 9/25. This is insufficient for global scaling and motivates dynamic discovery.
