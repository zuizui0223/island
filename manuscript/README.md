# Chapter 1 submission package

## Current state

The scientific analysis is frozen and locally verified against the exact artifact named in
`config/chapter1_submission_freeze_lock.json`. This directory is a journal-neutral writing
package, not yet a journal-formatted submission.

Files:

- `chapter1_manuscript.md`: main-text draft with the frozen claim ceiling;
- `supporting_information.md`: analysis hierarchy, robustness matrix, and attrition rules;
- `figure_table_plan.md`: every planned display mapped to locked source files.

## Submission gates

- [x] one immutable analysis lock identifies the PR #138 -> PR #139 -> area chain;
- [x] automated row-count, result-state, digest, and claim-ceiling validation exists;
- [x] main-text and SI narrative drafts exist;
- [ ] copy the expiring Actions artifact to a durable release/repository and record its DOI
  or permanent URL;
- [ ] generate and visually inspect final figures and tables from the locked artifact;
- [ ] complete reference-manager import and verify every in-text citation;
- [ ] choose a journal and apply its word-count, section, figure, data-policy, and reporting
  requirements;
- [ ] merge or otherwise archive the currently stacked draft PR history without changing
  the locked estimands.

PR #140 and subsequent trait acquisition are post-freeze and are not submission inputs.
