# Visible morphology support=2 checkpoint (2026-08-14)

This checkpoint adds 11 individually reviewed species/synonym-direct trait
records for 10 species. Sources are original or curated species treatments from
PROSEA, WFO-hosted floras, Kew-hosted floras, a university flora, a primary
taxonomic revision, and a peer-reviewed botanical-garden morphology study.

The records remain at trait grain. They do not emit genus inference and do not
count acquisition-queue potential as coverage. Cassipourea elliptica retains the
explicit `solitary|umbel_corymb` variation; it is counterevidence, not a forced
rule unlock. Scolopia rhamniphylla contributes only its stated inflorescence;
petal dimensions were not treated as whole-flower size. Pollen vector and reward
traits remain separate from the strict flower-structure axis.

Run `python -m island_v2.visible_morphology_wave6_checkpoint` with the previous
combined evidence and audit checkpoint to reproduce the four CSV files and
manifest in this directory. Formal High/Medium/Validated-Low coverage is valid
only after the shared all-evidence workflow rebuilds its artifact.
