# Result ownership by actual biological input — 25 November 2026

[Machine-readable result registry](../research/result_input_provenance_20261125.json) maps 51 rows from the pinned #347 evidence ledger plus 2 continuation/readiness entries. Counts are sample support, never source identifiers. Each result retains its original status, claim limit and URL. Family-level audit entries must be checked against exact runner and input hashes before promotion. Unknown dependencies cannot qualify for external-only status.

| Data that determine the result, including cohort selection and training | Owner |
|---|---|
| Project GBIF incidence and/or integrated floral traits (including derived public aggregates) | island |
| Project biological input + official GIFT checklists, mainland donor pools or animal/experimental inputs | island, mixed |
| Official GIFT checklists + independent publisher trait values, no own-derived input | fd |
| External network/experiment only, no own-derived selection | fd, external source |
| Unknown ancestry and no established own dependency | hold |

Shared generic code and outcome-blind geography/climate do not by themselves create biological dependence. A geographic sample selected on own flora/trait coverage DOES. External-only does not prove independence of original species, literature or observations.

**Important correction:** the old EuPPollNet guild-removal simulation has external network observations but is restricted to eleven islands from our 928 triple-Direct cohort; the entire selected result is mixed and belongs to island. A fresh network-only cohort could be an external fd study. Existing frozen GEB H3 (GloPL isolation gradient) is external-source supporting evidence, not a second new fd claim; H4 joins own trait values and is mixed. Do not move or modify the frozen paper.

The new GIFT publisher pilot (fd #72, [source-audit run](https://github.com/zuizui0223/fd/actions/runs/45560668861)) is external-only but coverage only. Official strict colour × maximum length has 404 shared work IDs; in the bounded 24-list island pilot only two lists have any overlap and the maximum is four carriers. Do not infer a null ecological effect from low support. First check official checklist geography and trait-definition/reference overlap, then predeclare a comparable estimand. Use separate data and output manifests for the four flora × trait source combinations; mixed cross-checks stay island.

Registry classifications are declared ancestry checks, not automatic code instrumentation. Test with `python scripts/check_result_input_provenance.py` and `python -m pytest tests/test_result_input_provenance.py -q`. No model was rerun; no raw protected data copied; no GEB lock changed.
