# Wave56 reproductive-assurance recovery

## Goal

Wave56 is dedicated to the Chapter 1 bottleneck: `reproductive_assurance`.

Canonical working checkpoint before Wave56 promotion:

- fixed analysis universe: 106,295 species;
- reproductive-assurance resolved cells: **48,527 / 106,295 = 45.65%**;
- unresolved reproductive cells: **57,768**;
- first milestone: **50% = 53,148 resolved species**;
- net new cells required for the 50% milestone: **4,621**.

The overall Chapter 1 checkpoint remains 222,759 / 318,885 = 69.8556%. Wave56 does not change that value until a record passes the strict evidence contract and resolves a previously unresolved cell in the current private-plus-public canonical ledger.

## Acquisition rule

This packet is intentionally conservative.

1. Species-level evidence and exact accepted-name reconciliation come first.
2. `self-compatible` is not automatically recoded as autonomous selfing.
3. Mixed mating, self-compatibility, autonomous selfing and self-incompatibility remain distinct source traits even though they contribute to the same reproductive-assurance axis.
4. External congeners can train a trait-specific genus rule only after underlying source lineage, accepted identity, dominance, species LOO and source-lineage LOO all pass.
5. Family inference and global fallback are prohibited.
6. Contradictory evidence is retained as a repair target, not averaged away or used to inflate coverage.

## Current Wave56 packet

`reviewed_direct_evidence.csv` contains source-grounded species-level records that are worth collision auditing against the current canonical ledger:

- `Schoenoplectiella juncoides`: successful self-pollination seed production supports `self_incompatibility=SC`; this ports the already reviewed Wave54 record into the current recovery lane.
- `Mesosphaerum suaveolens` (source name `Hyptis suaveolens`): the original breeding-system study reports both autogamy and allogamy, supporting `mating_system=mixed_mating`.
- `Ligustrum lucidum`: a species-level reproductive study describes a mixed mating system with self-compatibility; retained at Medium because the breeding-system statement cites earlier work rather than being the focal hand-cross experiment of that paper.
- `Illicium anisatum`: controlled self- and cross-pollination showed no significant difference in fruit set, supporting `self_incompatibility=SC` without inferring autonomous selfing.

All rows remain `promotion_allowed=false` until collision/conflict audit against the current canonical ledger.

`external_congener_candidates.csv` contains useful but not yet rule-trainable evidence. In particular, the 2018 `Cyanotis longifolia` paper states that the species is self-compatible but attributes that statement to earlier literature; the underlying lineage must be recovered before treating it as independent genus-rule support.

`conflict_repair_queue.csv` records evidence that should reduce false confidence rather than raise coverage. Two important examples are `Sideroxylon obtusifolium`, where controlled field reproduction supports facultative autogamy despite an existing SI-coded lineage, and `Illicium floridanum`, where a later controlled study found no prezygotic SI and attributed self-sterility primarily to early inbreeding depression while leaving late post-zygotic SI unresolved.

## High-yield queue

`support2_priority_queue.csv` is a recovery map derived from the Wave52 support-two diagnostic. Its potential-cell counts are **stale diagnostic upper bounds**, not current guaranteed gains. It is used only to choose which genera deserve a third independent species-level source search.

The main targets are currently `Spermacoce`, `Dicliptera`, `Portulaca`, `Myriophyllum`, `Schoenoplectiella`, `Lindernia` and `Cyanotis`. `Sideroxylon` and `Illicium` are explicitly excluded from genus-wide promotion while their conflicts are unresolved.

## Source anchors

- Schoenoplectiella juncoides: Sada et al. 2013, DOI `10.1016/j.pestbp.2013.05.013`.
- Mesosphaerum suaveolens / Hyptis suaveolens: Aluri 1990, DOI `10.1111/j.1442-1984.1990.tb00183.x`.
- Ligustrum lucidum: Aguirre-Acosta et al. 2014, DOI `10.1007/s10530-013-0577-x`.
- Illicium anisatum: Takahashi & Yamauchi, *Pollination biology of four ANITA-group species*, Gifu University repository, handle `20.500.12099/2822`.
- Cyanotis longifolia candidate: Delhaye et al. 2018, DOI `10.1016/j.envexpbot.2018.09.001`.
- Sideroxylon obtusifolium repair: Kiill et al. 2014, *Revista Arvore* 38(6):1015-1025.
- Illicium floridanum repair: Koehl et al. 2004, DOI `10.1093/aob/mch113`.
