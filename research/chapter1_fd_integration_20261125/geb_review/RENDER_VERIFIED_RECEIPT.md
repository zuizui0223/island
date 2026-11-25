# Integrated GEB review candidate — verified render receipt

**Author review only; not ready for external submission.**

[Render run 45593987569](https://github.com/zuizui0223/island/actions/runs/45593987569) succeeded at ee985b06c0ccb733df8267af96830ec2f512ea0b. [Text QA run 45593987376](https://github.com/zuizui0223/island/actions/runs/45593987376) also succeeded.

| Download | GitHub artifact | Retention |
|---|---|---|
| [DOCX/PDF review bundle](https://github.com/zuizui0223/island/actions/runs/45593987569/artifacts/15012142877) | 15012142877 | 2026-12-25 10:21 UTC |
| [All-page QA images](https://github.com/zuizui0223/island/actions/runs/45593987569/artifacts/15013135181) | 15013135181 | 2026-12-09 10:21 UTC |

Both outer ZIP archives were actually retrieved into the working session. Their byte SHA256 digests matched the GitHub metadata. The inner author-review bundle was extracted and all 42 entries listed in its review manifest were independently rehashed: 42/42 match. The anonymous manuscript and SI PDFs have author metadata "Anonymous". Manuscript 17 pages, SI 19, cover and title one page each. All 38 pages were rasterized by CI; no automated bounds or possible-blank-page flags. Contact sheets of all 36 manuscript/SI pages were inspected; no gross layout defects were observed. S8b and Figure S8 were also inspected enlarged; the new figure and caption share a page. This is not a completed full-resolution human proof of every page, and does not guarantee perfect typesetting or anonymity.

No further exploratory fit, source-ledger redistribution, final DOI, journal upload or email has occurred. Existing science and GEB submission selectors remain unchanged.

## Next blocking work: preserve the original source lawfully

The original protected artifact (run 34232450884, artifact 10058653212, 105,740,636 bytes, SHA256 b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999) is still listed as available and expires **2026-12-07 13:31 UTC (22:31 JST)**. The source contains analysis-used third-party material; do not place its bytes in a public GitHub tree, release, Zenodo deposit or newly shared artifact.

For the authorized source custodian, before expiry:
```sh
gh api repos/zuizui0223/island/actions/artifacts/10058653212/zip > chapter1-original-10058653212.zip
echo 'b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999  chapter1-original-10058653212.zip' | sha256sum --check
```
Run inside existing access-controlled, rights-appropriate storage only; do not infer authorization to copy the source to other providers. Retain the original archive, provenance and permission inventory, and verify the exact internal source hashes before using it again. The raw source ZIP has **not** been transferred or newly preserved in this turn. Its expiry date and digest were rechecked from GitHub metadata. The original WCVP sensitivity artifacts remain unavailable, and this preservation alone does not recover them.

The next approval gates are the source-level rights determination (176,414 cells restricted or still under review), any editor-approved lawful reviewer-access route, final author information and all-author approval, anonymous access and archive authorization. Existing review outputs and scripts are ready for internal author review. Do not treat this receipt as a submission authorization.
