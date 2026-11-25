# Island floral beta implementation plan

Goal: species turnover versus floral composition similarity, not within-species divergence or adaptation. The user approved this question after explicitly noting the lack of island/mainland population trait measurements.

Main population: all 4505 original GIFT floras. Ten independent trait masks and three genuinely observed two-trait masks; all-recorded and source-native scopes remain separate. No third trait, shared-genus, photo or guild-state prerequisite. Each endpoint retains its original >=5 observed species rule, checked against archived prior counts and raw Q.

Methods and contrasts fixed before new original-data results:

1. Categorical species-proportion affinity B=sum sqrt(p*q), range 0..1. Exact mean across ALL unordered eligible island pairs using sums of square roots. Profile Jaccard/Sorensen with turnover/nestedness partition, not historical loss and not hypervolume.
2. Every eligible island is a focal in a deterministic sample of 16 uniform comparison islands selected without category values. Other pair distributions are estimated from this graph and focal-island equally weighted. Three geographic-nearest comparators are a separate descriptive graph. Never treat pair counts as independent ecological replicates.
3. In independently typed geology, enumerate all within-broad-region OO/CC/OC pairs and summarize symmetrically with island-equal weighting. These are descriptive contrasts, not climate-adjusted geological causal effects.
4. Taxonomic Jaccard is calculated for both the complete recorded flora and the exact trait-known subset. In the preselected global-uniform graph, isolate pairs whose complete recorded floras share NO species. This is a secondary contrast, not a restriction on the main population.
5. Three primary endpoints (colour, form, colour-form) each receive 199 global and 199 within-genus whole-profile permutations, preserving original species occurrence/missingness and global category frequencies. Repeated species get the SAME permuted label at every island; pair profiles are shuffled intact, never fabricated by cross-products. Test excess disjoint-flora affinity, fixed six-slot exploratory BH family. Lower-tail diagnostic separately labelled. Null pool is global observed species, not ancestors and not a missingness correction.
6. Test identities and cohort retention locally, then pinned source-protected original CI. Reconcile all 26 baseline masks/raw Q. Export only anonymous aggregates, never taxon/island/pair keys or coordinates. Add Japanese results with mechanism limits.

Progress: local seven tests first failed with missing implementation, then all seven passed. Original data results pending when this plan was committed. Species-incidence and profile null tests do not establish pollination effectiveness or intraspecific evolution.
