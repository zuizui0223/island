# Supplemental private backup and offline restoration

All four current-lineage source archives remain preserved under the existing [Drive backup receipt](DRIVE_BACKUP_VERIFIED.json). On 26 November 2026 JST, three additional original artifacts were downloaded, stored in separate owner-only subfolders of the authorized private Drive, then downloaded again and verified against their exact GitHub ZIP SHA256, sizes and CRCs:

| Artifact | Purpose | Bytes |
|---|---|---:|
| 10288721725 | DB 2.0 alpha1 — developmental, not the active Chapter 1 database | 83,300,729 |
| 15012142877 | Integrated GEB author-review bundle — not submission-approved | 3,708,112 |
| 15013135181 | All-page QA images | 9,023,416 |

The additional 96,032,257 bytes are private preservation copies, not public redistribution. No changes to source licences, review access or the frozen selector.

## Restore without GitHub

Download the two original parts and their `.parts.json` manifest to an authorized local directory. Create a separate authorized local destination. Run:

```sh
python restore_island_archive.py --source-dir /authorized/downloaded-parts --output-dir /authorized/restored
```

The script verifies both part pins and their order against the supplied manifest and pinned original identity, concatenates into a temporary file, verifies the complete 105,740,636-byte archive SHA256 and CRC, and never replaces an existing file. A .bin suffix added by Drive download is accepted only if its bytes match the expected part. It uses no network. Tested on actual Drive-readback parts, including the second part's .bin suffix, in a clean temporary output directory; exact restoration succeeded, repeat restoration was idempotent, and reversing the manifest part order was rejected. The script and test receipt are stored privately alongside the original parts. No source records are published.

Exact-name Drive searches for the missing WCVP native-range and TDWG mapping files, and the two original artifact names, returned no matches. This does not prove that other copies are absent; the historical WCVP replay remains unavailable. No substitution with alpha1 or a new WCVP snapshot. Scientific review and external approvals are unchanged.

See [verification receipt](SUPPLEMENTAL_DRIVE_BACKUP_VERIFIED.json).
